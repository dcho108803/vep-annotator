/**
 * dbNSFP Annotation Source
 *
 * Provides pathogenicity predictions and conservation scores from dbNSFP database.
 * dbNSFP is a database of all potential non-synonymous single nucleotide variants (nsSNVs).
 */

#include "annotation_source.hpp"
#include "file_parsers.hpp"
#include "dbnsfp_fields.hpp"
#include "vep_annotator.hpp"
#include <sstream>
#include <algorithm>
#include <cmath>

namespace vep {

/**
 * dbNSFP Annotation Source
 *
 * Queries dbNSFP database for pathogenicity predictions and conservation scores.
 * Requires a tabix-indexed dbNSFP file (.txt.gz + .tbi).
 */
class DbNSFPSource : public VariantAnnotationSource {
public:
    /**
     * Construct dbNSFP source
     * @param path Path to tabix-indexed dbNSFP file
     * @param fields Fields to extract (empty = essential preset)
     */
    DbNSFPSource(const std::string& path, const std::string& field_spec = "essential")
        : path_(path), field_spec_(field_spec) {
        requested_fields_ = parse_dbnsfp_fields(field_spec);
    }

    std::string name() const override { return "dbnsfp"; }
    std::string type() const override { return "pathogenicity"; }
    std::string description() const override {
        return "dbNSFP pathogenicity predictions and conservation scores";
    }

    bool is_ready() const override { return reader_ != nullptr && reader_->is_valid(); }
    bool is_thread_safe() const override { return false; }  // Tabix file handles not thread-safe

    void initialize() override {
        std::lock_guard<std::recursive_mutex> lock(mutex_);
        if (reader_) return;

        log(LogLevel::INFO, "Loading dbNSFP from: " + path_);

        // dbNSFP uses 1-based chromosome column (0) and position column (1)
        // Columns: #chr, pos(1-based), ref, alt, ...
        reader_ = std::make_unique<TabixTSVReader>(path_, 0, 1);

        if (!reader_->is_valid()) {
            log(LogLevel::ERROR, "Failed to open dbNSFP file: " + path_);
            reader_.reset();
            return;
        }

        // Get available columns and build index
        auto columns = reader_->get_columns();
        for (size_t i = 0; i < columns.size(); ++i) {
            column_index_[columns[i]] = i;
        }

        // Map requested fields to column indices
        for (const auto& field : requested_fields_) {
            auto it = column_index_.find(field.column);
            if (it != column_index_.end()) {
                field_to_column_[field.name] = it->second;
            } else {
                log(LogLevel::WARNING, "dbNSFP column not found: " + field.column);
            }
        }

        log(LogLevel::INFO, "dbNSFP loaded with " + std::to_string(field_to_column_.size()) +
                           " fields from " + std::to_string(columns.size()) + " columns");
    }

    void annotate(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript* transcript,
        std::unordered_map<std::string, std::string>& annotations
    ) override {
        ensure_initialized();
        if (!reader_) return;

        // Query dbNSFP at this position
        auto records = reader_->query(chrom, pos);

        // Find matching record by ref/alt
        for (const auto& record : records) {
            // Check ref allele match
            auto ref_it = record.find("ref");
            if (ref_it != record.end() && ref_it->second != ref) {
                continue;
            }

            // Check alt allele match
            auto alt_it = record.find("alt");
            if (alt_it != record.end()) {
                // dbNSFP may have multiple alts separated by ';'
                std::string alts = alt_it->second;
                bool found_alt = false;
                size_t alt_index = 0;

                std::istringstream iss(alts);
                std::string single_alt;
                size_t idx = 0;
                while (std::getline(iss, single_alt, ';')) {
                    if (single_alt == alt) {
                        found_alt = true;
                        alt_index = idx;
                        break;
                    }
                    ++idx;
                }

                if (!found_alt) continue;

                // Determine per-transcript index from Ensembl_transcriptid column
                // Many dbNSFP score fields (SIFT, PolyPhen, etc.) are semicolon-separated
                // per-transcript, NOT per-alt. The Ensembl_transcriptid column maps indices.
                int transcript_index = -1;
                if (transcript) {
                    auto tid_it = record.find("Ensembl_transcriptid");
                    if (tid_it != record.end() && !tid_it->second.empty()) {
                        std::string query_id = strip_version(transcript->id);
                        std::istringstream tss(tid_it->second);
                        std::string tid;
                        int tidx = 0;
                        while (std::getline(tss, tid, ';')) {
                            if (strip_version(tid) == query_id) {
                                transcript_index = tidx;
                                break;
                            }
                            ++tidx;
                        }
                    }
                }

                // Extract requested fields
                for (const auto& field : requested_fields_) {
                    auto col_it = record.find(field.column);
                    if (col_it != record.end()) {
                        std::string value = col_it->second;

                        // Handle multi-value fields (separated by ';')
                        if (value.find(';') != std::string::npos) {
                            std::vector<std::string> values;
                            std::istringstream vss(value);
                            std::string v;
                            while (std::getline(vss, v, ';')) {
                                values.push_back(v);
                            }

                            value = select_value_for_transcript(
                                values, transcript_index, alt_index);
                        }

                        // Skip missing values
                        if (value != "." && !value.empty()) {
                            annotations["dbnsfp:" + field.name] = value;
                        }
                    }
                }

                // Found matching record, stop searching
                return;
            }
        }
    }

    std::vector<std::string> get_fields() const override {
        std::vector<std::string> fields;
        for (const auto& f : requested_fields_) {
            fields.push_back("dbnsfp:" + f.name);
        }
        return fields;
    }

    bool requires_allele_match() const override { return true; }

    std::unordered_map<std::string, std::string> query(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt
    ) const override {
        std::unordered_map<std::string, std::string> result;

        if (!reader_) return result;

        auto records = reader_->query(chrom, pos);

        for (const auto& record : records) {
            auto ref_it = record.find("ref");
            if (ref_it != record.end() && ref_it->second != ref) continue;

            auto alt_it = record.find("alt");
            if (alt_it != record.end()) {
                // Properly match alt allele using semicolon-split (like annotate())
                std::string alts = alt_it->second;
                bool found_alt = false;
                size_t alt_index = 0;

                std::istringstream aiss(alts);
                std::string single_alt;
                size_t aidx = 0;
                while (std::getline(aiss, single_alt, ';')) {
                    if (single_alt == alt) {
                        found_alt = true;
                        alt_index = aidx;
                        break;
                    }
                    ++aidx;
                }

                if (!found_alt) continue;

                // Found match - copy all requested fields
                // query() has no transcript context, so fall back to first valid value
                for (const auto& field : requested_fields_) {
                    auto col_it = record.find(field.column);
                    if (col_it != record.end() && col_it->second != "." && !col_it->second.empty()) {
                        std::string value = col_it->second;
                        // Handle multi-value fields (separated by ';')
                        if (value.find(';') != std::string::npos) {
                            std::vector<std::string> values;
                            std::istringstream vss(value);
                            std::string v;
                            while (std::getline(vss, v, ';')) {
                                values.push_back(v);
                            }
                            // No transcript context available, use first valid value
                            value = select_value_for_transcript(
                                values, -1, alt_index);
                        }
                        if (value != "." && !value.empty()) {
                            result["dbnsfp:" + field.name] = value;
                        }
                    }
                }
                break;
            }
        }

        return result;
    }

    std::string get_data_path() const override { return path_; }

    size_t memory_usage() const override {
        // Approximate memory usage
        return column_index_.size() * 64 + field_to_column_.size() * 32;
    }

    /**
     * Get list of requested fields
     */
    const std::vector<DbNSFPField>& get_requested_fields() const {
        return requested_fields_;
    }

    /**
     * Get interpretation for a pathogenicity score
     */
    static std::string interpret_score(const DbNSFPField& field, double score) {
        if (field.damaging_threshold < 0) return "";

        if (field.higher_is_damaging) {
            return score >= field.damaging_threshold ? "D" : "T";
        } else {
            return score <= field.damaging_threshold ? "D" : "T";
        }
    }

private:
    /**
     * Strip version suffix from transcript ID (e.g., "ENST00000269305.8" -> "ENST00000269305")
     */
    static std::string strip_version(const std::string& id) {
        auto dot = id.rfind('.');
        if (dot != std::string::npos) {
            // Only strip if what follows the dot is numeric (a version number)
            bool all_digits = true;
            for (size_t i = dot + 1; i < id.size(); ++i) {
                if (!std::isdigit(static_cast<unsigned char>(id[i]))) {
                    all_digits = false;
                    break;
                }
            }
            if (all_digits && dot + 1 < id.size()) {
                return id.substr(0, dot);
            }
        }
        return id;
    }

    /**
     * Select the best value from a semicolon-split multi-value field.
     *
     * Many dbNSFP fields (SIFT_score, Polyphen2_HDIV_score, etc.) are per-transcript:
     * the semicolons separate values for different transcripts (matching Ensembl_transcriptid).
     *
     * Strategy:
     * 1. If a matching transcript index was found, use that value (if valid).
     * 2. Otherwise, fall back to the first non-empty/non-"." value.
     */
    static std::string select_value_for_transcript(
        const std::vector<std::string>& values,
        int transcript_index,
        size_t /*alt_index*/)
    {
        // Try transcript-matched index first
        if (transcript_index >= 0 &&
            static_cast<size_t>(transcript_index) < values.size()) {
            const std::string& v = values[static_cast<size_t>(transcript_index)];
            if (!v.empty() && v != ".") {
                return v;
            }
        }

        // Fall back to first valid value
        for (const auto& v : values) {
            if (!v.empty() && v != ".") {
                return v;
            }
        }

        // All values are missing
        return ".";
    }

    std::string path_;
    std::string field_spec_;
    std::vector<DbNSFPField> requested_fields_;
    std::unique_ptr<TabixTSVReader> reader_;
    std::map<std::string, size_t> column_index_;
    std::map<std::string, size_t> field_to_column_;
};

/**
 * Factory function to create dbNSFP source
 */
std::shared_ptr<AnnotationSource> create_dbnsfp_source(
    const std::string& path,
    const std::string& fields
) {
    return std::make_shared<DbNSFPSource>(path, fields);
}

} // namespace vep
