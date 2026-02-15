/**
 * Output Writer - Multiple Output Format Support
 *
 * Supports TSV (default), JSON, and VCF output formats.
 */

#ifndef OUTPUT_WRITER_HPP
#define OUTPUT_WRITER_HPP

#include "vep_annotator.hpp"
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <cctype>
#include <zlib.h>

namespace vep {

/**
 * Output format types
 */
enum class OutputFormat {
    TSV,    // Tab-separated values (default)
    JSON,   // JSON format
    VCF     // VCF with CSQ INFO field
};

/**
 * Parse output format from string
 */
inline OutputFormat parse_output_format(const std::string& format) {
    std::string lower = format;
    for (size_t i = 0; i < lower.size(); ++i) {
        lower[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(lower[i])));
    }

    if (lower == "json") return OutputFormat::JSON;
    if (lower == "vcf") return OutputFormat::VCF;
    return OutputFormat::TSV;
}

/**
 * Statistics collector for annotation summary
 */
struct AnnotationStats {
    int total_variants = 0;
    int annotated_variants = 0;
    std::map<std::string, int> consequence_counts;
    std::map<std::string, int> impact_counts;
    std::map<std::string, int> biotype_counts;

    void add(const VariantAnnotation& ann) {
        total_variants++;
        if (!ann.gene_symbol.empty()) {
            annotated_variants++;
        }

        for (const auto& csq : ann.consequences) {
            consequence_counts[consequence_to_string(csq)]++;
        }
        impact_counts[impact_to_string(ann.impact)]++;
        if (!ann.biotype.empty()) {
            biotype_counts[ann.biotype]++;
        }
    }

    std::string to_string() const {
        std::ostringstream oss;
        oss << "=== Annotation Statistics ===\n";
        oss << "Total variants: " << total_variants << "\n";
        oss << "Annotated variants: " << annotated_variants << "\n";
        oss << "\nConsequence counts:\n";
        for (const auto& pair : consequence_counts) {
            oss << "  " << pair.first << ": " << pair.second << "\n";
        }
        oss << "\nImpact counts:\n";
        for (const auto& pair : impact_counts) {
            oss << "  " << pair.first << ": " << pair.second << "\n";
        }
        return oss.str();
    }

    std::string to_json() const {
        std::ostringstream oss;
        oss << "{\n";
        oss << "  \"total_variants\": " << total_variants << ",\n";
        oss << "  \"annotated_variants\": " << annotated_variants << ",\n";
        oss << "  \"consequence_counts\": {";
        bool first = true;
        for (const auto& pair : consequence_counts) {
            if (!first) oss << ",";
            oss << "\n    \"" << pair.first << "\": " << pair.second;
            first = false;
        }
        oss << "\n  },\n";
        oss << "  \"impact_counts\": {";
        first = true;
        for (const auto& pair : impact_counts) {
            if (!first) oss << ",";
            oss << "\n    \"" << pair.first << "\": " << pair.second;
            first = false;
        }
        oss << "\n  }\n";
        oss << "}";
        return oss.str();
    }
};

/**
 * Abstract base class for output writers
 */
class OutputWriter {
public:
    virtual ~OutputWriter() = default;

    virtual void write_header(const std::vector<std::string>& custom_columns) = 0;
    virtual void write_annotation(const VariantAnnotation& ann,
                                  const std::vector<std::string>& custom_columns) = 0;
    virtual void write_annotations(const std::vector<VariantAnnotation>& anns,
                                   const std::vector<std::string>& custom_columns) = 0;
    virtual void write_footer() = 0;
    virtual void close() = 0;

    // Convenience overloads without custom columns
    void write_header() { write_header(std::vector<std::string>()); }
    void write_annotation(const VariantAnnotation& ann) {
        write_annotation(ann, std::vector<std::string>());
    }
    void write_annotations(const std::vector<VariantAnnotation>& anns) {
        write_annotations(anns, std::vector<std::string>());
    }

    const AnnotationStats& get_stats() const { return stats_; }

    void set_skip_header(bool v) { skip_header_ = v; }
    void set_term_style(const std::string& style) { term_style_ = style; }

protected:
    AnnotationStats stats_;
    bool skip_header_ = false;
    std::string term_style_ = "SO";

    // Get consequence string respecting term style
    std::string format_consequence(const VariantAnnotation& ann) const {
        if (term_style_ == "display") {
            std::string result;
            for (size_t i = 0; i < ann.consequences.size(); ++i) {
                if (i > 0) result += "&";
                result += consequence_to_display_term(ann.consequences[i]);
            }
            return result;
        }
        return ann.get_consequence_string();
    }
};

// Helper to check if path ends with .gz
inline bool ends_with_gz(const std::string& path) {
    return path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
}

// Format position with optional range and /total: "start-end/total", "start/total", or "start"
// empty_val is returned when position <= 0 (e.g., "-" for TSV, "" for VCF)
inline std::string format_position_with_total(
    int start, int end, const std::string& total_key,
    const VariantAnnotation& ann, const std::string& empty_val) {
    if (start <= 0) return empty_val;
    std::string s = std::to_string(start);
    if (end > 0 && end != start) s += "-" + std::to_string(end);
    auto it = ann.custom_annotations.find(total_key);
    if (it != ann.custom_annotations.end() && !it->second.empty()) s += "/" + it->second;
    return s;
}

// Format protein position with optional range and /total (total = CDS_LENGTH/3)
inline std::string format_protein_position(const VariantAnnotation& ann, const std::string& empty_val) {
    if (ann.protein_position <= 0) return empty_val;
    std::string s = std::to_string(ann.protein_position);
    if (ann.protein_end > 0 && ann.protein_end != ann.protein_position)
        s += "-" + std::to_string(ann.protein_end);
    auto cl = ann.custom_annotations.find("CDS_LENGTH");
    if (cl != ann.custom_annotations.end() && !cl->second.empty()) {
        try {
            int cds_len = std::stoi(cl->second);
            if (cds_len > 0) s += "/" + std::to_string(cds_len / 3);
        } catch (...) {}
    }
    return s;
}

/**
 * TSV output writer (default format)
 */
class TSVWriter : public OutputWriter {
public:
    explicit TSVWriter(const std::string& output_path, bool compress = false)
        : output_path_(output_path), compress_(compress), gz_file_(nullptr),
          use_stdout_(output_path.empty() || output_path == "-" || output_path == "STDOUT") {

        if (use_stdout_) {
            compress_ = false;  // Cannot compress stdout
        } else if (compress_ || ends_with_gz(output_path_)) {
            compress_ = true;
            gz_file_ = gzopen(output_path_.c_str(), "wb");
            if (!gz_file_) {
                throw std::runtime_error("Cannot open output file: " + output_path_);
            }
        } else {
            output_.open(output_path_);
            if (!output_.is_open()) {
                throw std::runtime_error("Cannot open output file: " + output_path_);
            }
        }
    }

    ~TSVWriter() override {
        close();
    }

    void write_header(const std::vector<std::string>& custom_columns) override {
        if (skip_header_) return;

        std::ostringstream header;
        // Perl VEP header comment line
        header << "## ENSEMBL VARIANT EFFECT PREDICTOR\n";
        // Perl VEP column names
        header << "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\t"
               << "Feature_type\tConsequence\tcDNA_position\tCDS_position\t"
               << "Protein_position\tAmino_acids\tCodons\tExisting_variation\tExtra";

        header << "\n";

        write_string(header.str());
        custom_columns_ = custom_columns;
    }

    void write_annotation(const VariantAnnotation& ann,
                          const std::vector<std::string>& custom_columns) override {
        stats_.add(ann);

        std::ostringstream line;

        // #Uploaded_variation: Use VCF ID (rs#) if available, else CHR_POS_ALLELES
        if (!ann.vcf_id.empty() && ann.vcf_id != ".") {
            line << ann.vcf_id << "\t";
        } else {
            // Ensembl format: CHR_POS_REF/ALT using display alleles
            std::string d_ref = ann.display_ref.empty() ? ann.ref_allele : ann.display_ref;
            std::string d_alt = ann.display_alt.empty() ? ann.alt_allele : ann.display_alt;
            int d_start = ann.display_start > 0 ? ann.display_start : ann.position;
            line << ann.chromosome << "_" << d_start << "_" << d_ref << "/" << d_alt << "\t";
        }

        // Location: CHROM:POS or CHROM:POS-END using display coords
        {
            int d_start = ann.display_start > 0 ? ann.display_start : ann.position;
            int d_end = ann.display_end > 0 ? ann.display_end :
                        (ann.position + static_cast<int>(ann.ref_allele.size()) - 1);
            line << ann.chromosome << ":" << d_start;
            if (d_end > d_start) {
                line << "-" << d_end;
            }
        }
        line << "\t";

        // Allele: the display alt allele (Ensembl format: "-" for deletions)
        line << (ann.display_alt.empty() ? ann.alt_allele : ann.display_alt) << "\t";

        // Gene
        line << ann.gene_id << "\t";

        // Feature (transcript ID)
        line << ann.transcript_id << "\t";

        // Feature_type
        line << ann.feature_type << "\t";

        // Consequence
        line << format_consequence(ann) << "\t";

        // cDNA_position: "start-end/total" for multi-base, "start/total" for single, "start" without --total-length
        line << format_position_with_total(ann.cdna_position, ann.cdna_end, "TRANSCRIPT_LENGTH", ann, "-") << "\t";

        // CDS_position
        line << format_position_with_total(ann.cds_position, ann.cds_end, "CDS_LENGTH", ann, "-") << "\t";

        // Protein_position (total = CDS_LENGTH/3)
        line << format_protein_position(ann, "-") << "\t";

        // Amino_acids
        line << (ann.amino_acids.empty() ? "-" : ann.amino_acids) << "\t";

        // Codons
        line << (ann.codons.empty() ? "-" : ann.codons) << "\t";

        // Existing_variation
        line << (ann.existing_variation.empty() ? "-" : ann.existing_variation) << "\t";

        // Extra: key=value pairs separated by semicolons
        std::ostringstream extra;
        bool first_extra = true;

        auto append_extra = [&](const std::string& key, const std::string& value) {
            if (!value.empty()) {
                if (!first_extra) extra << ";";
                extra << key << "=" << value;
                first_extra = false;
            }
        };

        // Perl VEP Extra field ordering: IMPACT, DISTANCE, STRAND, FLAGS first (from 'user' flag),
        // then SYMBOL group, then BIOTYPE, CANONICAL, then others
        append_extra("IMPACT", impact_to_string(ann.impact));
        if (ann.distance > 0) {
            append_extra("DISTANCE", std::to_string(ann.distance));
        }
        if (ann.strand != '\0') {
            append_extra("STRAND", ann.strand == '+' ? "1" : "-1");
        }
        // FLAGS (cds_start_NF, cds_end_NF) from transcript metadata
        {
            auto flags_it = ann.custom_annotations.find("FLAGS");
            if (flags_it != ann.custom_annotations.end() && !flags_it->second.empty()) {
                append_extra("FLAGS", flags_it->second);
            }
        }
        if (!ann.gene_symbol.empty()) {
            append_extra("SYMBOL", ann.gene_symbol);
            auto ss_it = ann.custom_annotations.find("SYMBOL_SOURCE");
            if (ss_it != ann.custom_annotations.end() && !ss_it->second.empty()) {
                append_extra("SYMBOL_SOURCE", ss_it->second);
            }
            auto hgnc_it = ann.custom_annotations.find("HGNC_ID");
            if (hgnc_it != ann.custom_annotations.end() && !hgnc_it->second.empty()) {
                append_extra("HGNC_ID", hgnc_it->second);
            }
        }
        if (!ann.biotype.empty()) {
            append_extra("BIOTYPE", ann.biotype);
        }
        if (ann.is_canonical) {
            append_extra("CANONICAL", "YES");
        }
        if (ann.exon_number > 0) {
            append_extra("EXON", std::to_string(ann.exon_number) + "/" + std::to_string(ann.total_exons));
        }
        if (ann.intron_number > 0) {
            append_extra("INTRON", std::to_string(ann.intron_number) + "/" + std::to_string(ann.total_introns));
        }
        if (!ann.hgvsc.empty()) {
            append_extra("HGVSc", ann.hgvsc);
        }
        if (!ann.hgvsp.empty()) {
            append_extra("HGVSp", ann.hgvsp);
        }
        // HGVS_OFFSET (shift distance when HGVS 3' normalization applied)
        {
            auto offset_it = ann.custom_annotations.find("HGVS_OFFSET");
            if (offset_it != ann.custom_annotations.end() && !offset_it->second.empty()) {
                append_extra("HGVS_OFFSET", offset_it->second);
            }
        }
        if (!ann.hgvsg.empty()) {
            append_extra("HGVSg", ann.hgvsg);
        }
        if (!ann.source.empty()) {
            append_extra("SOURCE", ann.source);
        }

        // Append custom annotation columns as Extra key=value pairs.
        // Skip fields already output explicitly above to prevent duplication.
        static const std::set<std::string> tsv_phase1_fields = {
            "IMPACT", "DISTANCE", "STRAND", "FLAGS", "SYMBOL", "SYMBOL_SOURCE",
            "HGNC_ID", "BIOTYPE", "CANONICAL", "EXON", "INTRON",
            "HGVSc", "HGVSp", "HGVS_OFFSET", "HGVSg", "SOURCE"
        };
        for (const auto& col : custom_columns) {
            if (tsv_phase1_fields.count(col)) continue;  // Already output above
            auto it = ann.custom_annotations.find(col);
            if (it != ann.custom_annotations.end() && !it->second.empty()) {
                append_extra(col, it->second);
            }
        }

        std::string extra_str = extra.str();
        line << (extra_str.empty() ? "-" : extra_str);

        line << "\n";

        write_string(line.str());
    }

    void write_annotations(const std::vector<VariantAnnotation>& anns,
                           const std::vector<std::string>& custom_columns) override {
        for (const auto& ann : anns) {
            write_annotation(ann, custom_columns);
        }
    }

    void write_footer() override {
        // TSV has no footer
    }

    void close() override {
        if (gz_file_) {
            gzclose(gz_file_);
            gz_file_ = nullptr;
        }
        if (output_.is_open()) {
            output_.close();
        }
    }

private:
    std::string output_path_;
    bool compress_;
    std::ofstream output_;
    gzFile gz_file_;
    bool use_stdout_;
    std::vector<std::string> custom_columns_;

    void write_string(const std::string& s) {
        if (use_stdout_) {
            std::cout << s;
        } else if (compress_ && gz_file_) {
            gzwrite(gz_file_, s.c_str(), static_cast<unsigned int>(s.size()));
        } else {
            output_ << s;
        }
    }
};

/**
 * JSON output writer - Perl VEP compatible nested structure.
 *
 * Groups annotations by variant (chrom+pos+ref+alt) into per-variant objects
 * with a transcript_consequences array, matching Perl VEP JSON output.
 */
class JSONWriter : public OutputWriter {
public:
    explicit JSONWriter(const std::string& output_path, bool compress = false)
        : output_path_(output_path), compress_(compress), gz_file_(nullptr), first_variant_(true),
          use_stdout_(output_path.empty() || output_path == "-" || output_path == "STDOUT") {

        if (use_stdout_) {
            compress_ = false;
        } else if (compress_ || ends_with_gz(output_path_)) {
            compress_ = true;
            gz_file_ = gzopen(output_path_.c_str(), "wb");
            if (!gz_file_) {
                throw std::runtime_error("Cannot open output file: " + output_path_);
            }
        } else {
            output_.open(output_path_);
            if (!output_.is_open()) {
                throw std::runtime_error("Cannot open output file: " + output_path_);
            }
        }
    }

    ~JSONWriter() override {
        close();
    }

    void set_assembly_name(const std::string& name) { assembly_name_ = name; }

    void write_header(const std::vector<std::string>& /*custom_columns*/) override {
        if (!skip_header_) write_string("[\n");
    }

    void write_annotation(const VariantAnnotation& ann,
                          const std::vector<std::string>& /*custom_columns*/) override {
        stats_.add(ann);

        // Build variant key for grouping: chrom:pos:ref:alt
        std::string variant_key = ann.chromosome + ":" + std::to_string(ann.position)
                                  + ":" + ann.ref_allele + ":" + ann.alt_allele;

        if (variant_key != current_variant_key_) {
            // Flush the previous variant group
            flush_current_variant();
            current_variant_key_ = variant_key;
        }

        buffered_annotations_.push_back(ann);
    }

    void write_annotations(const std::vector<VariantAnnotation>& anns,
                           const std::vector<std::string>& custom_columns) override {
        for (const auto& ann : anns) {
            write_annotation(ann, custom_columns);
        }
    }

    void write_footer() override {
        flush_current_variant();
        if (!skip_header_) write_string("\n]\n");
    }

    void close() override {
        if (gz_file_) {
            gzclose(gz_file_);
            gz_file_ = nullptr;
        }
        if (output_.is_open()) {
            output_.close();
        }
    }

private:
    std::string output_path_;
    bool compress_;
    std::ofstream output_;
    gzFile gz_file_;
    bool first_variant_;
    bool use_stdout_;

    std::string assembly_name_ = "GRCh38";

    // Buffering for per-variant grouping
    std::string current_variant_key_;
    std::vector<VariantAnnotation> buffered_annotations_;

    void flush_current_variant() {
        if (buffered_annotations_.empty()) return;

        if (!first_variant_) {
            write_string(",\n");
        }
        first_variant_ = false;

        const auto& first = buffered_annotations_[0];

        // Determine most_severe_consequence across all transcript annotations
        Impact best_impact = Impact::MODIFIER;
        ConsequenceType most_severe_csq = ConsequenceType::UNKNOWN;
        for (const auto& ann : buffered_annotations_) {
            for (const auto& csq : ann.consequences) {
                Impact ci = get_impact(csq);
                if (static_cast<int>(ci) < static_cast<int>(best_impact) ||
                    (static_cast<int>(ci) == static_cast<int>(best_impact) &&
                     static_cast<int>(csq) < static_cast<int>(most_severe_csq))) {
                    best_impact = ci;
                    most_severe_csq = csq;
                }
            }
        }

        // Compute variant-level end position
        int end_pos = first.position + static_cast<int>(first.ref_allele.size()) - 1;
        if (end_pos < first.position) end_pos = first.position;

        std::ostringstream json;
        json << "  {\n";
        json << "    \"input\": \"" << escape_json(first.input_variant) << "\",\n";
        // Use VCF ID (e.g., rs#) if available, otherwise synthetic CHR_POS_REF/ALT
        if (!first.vcf_id.empty() && first.vcf_id != ".") {
            json << "    \"id\": \"" << escape_json(first.vcf_id) << "\",\n";
        } else {
            std::string id_ref = first.ref_allele.empty() ? "-" : first.ref_allele;
            std::string id_alt = first.alt_allele.empty() ? "-" : first.alt_allele;
            json << "    \"id\": \"" << escape_json(first.chromosome) << "_"
                 << first.position << "_"
                 << escape_json(id_ref) << "/" << escape_json(id_alt) << "\",\n";
        }
        json << "    \"assembly_name\": \"" << escape_json(assembly_name_) << "\",\n";
        json << "    \"seq_region_name\": \"" << escape_json(first.chromosome) << "\",\n";
        json << "    \"start\": " << first.position << ",\n";
        json << "    \"end\": " << end_pos << ",\n";
        {
            std::string as_ref = first.ref_allele.empty() ? "-" : first.ref_allele;
            std::string as_alt = first.alt_allele.empty() ? "-" : first.alt_allele;
            json << "    \"allele_string\": \"" << escape_json(as_ref) << "/"
                 << escape_json(as_alt) << "\",\n";
        }
        json << "    \"strand\": 1,\n";
        json << "    \"most_severe_consequence\": \""
             << (term_style_ == "display" ? consequence_to_display_term(most_severe_csq) : consequence_to_string(most_severe_csq))
             << "\",\n";

        // Variant class (when --variant-class enabled, stored as VARIANT_CLASS in custom_annotations)
        {
            auto vc_it = first.custom_annotations.find("VARIANT_CLASS");
            if (vc_it != first.custom_annotations.end() && !vc_it->second.empty()) {
                json << "    \"variant_class\": \"" << escape_json(vc_it->second) << "\",\n";
            }
        }

        // Co-located variants (Perl VEP structure: array of objects, one per co-located variant)
        if (!first.existing_variation.empty() && first.existing_variation != "-") {
            // Split comma-separated IDs into individual co-located variant objects
            std::vector<std::string> coloc_ids;
            {
                std::istringstream id_ss(first.existing_variation);
                std::string single_id;
                while (std::getline(id_ss, single_id, ',')) {
                    if (!single_id.empty()) coloc_ids.push_back(single_id);
                }
            }
            if (coloc_ids.empty()) coloc_ids.push_back(first.existing_variation);

            json << "    \"colocated_variants\": [";
            for (size_t ci = 0; ci < coloc_ids.size(); ++ci) {
                if (ci > 0) json << ", ";
                json << "{\n";
                json << "      \"id\": \"" << escape_json(coloc_ids[ci]) << "\",\n";
                json << "      \"seq_region_name\": \"" << escape_json(first.chromosome) << "\",\n";
                json << "      \"start\": " << first.position << ",\n";
                json << "      \"end\": " << end_pos << ",\n";
                json << "      \"allele_string\": \"" << escape_json(first.ref_allele) << "/"
                     << escape_json(first.alt_allele) << "\",\n";
                json << "      \"strand\": 1";
                // CLIN_SIG on first variant if available
                if (ci == 0) {
                    auto clin_it = first.custom_annotations.find("CLIN_SIG");
                    if (clin_it != first.custom_annotations.end() && !clin_it->second.empty()) {
                        json << ",\n      \"clin_sig\": [";
                        std::string cs = clin_it->second;
                        bool first_cs = true;
                        size_t start_cs = 0;
                        for (size_t i = 0; i <= cs.size(); ++i) {
                            if (i == cs.size() || cs[i] == ',' || cs[i] == '&' || cs[i] == '/') {
                                if (i > start_cs) {
                                    if (!first_cs) json << ", ";
                                    json << "\"" << escape_json(cs.substr(start_cs, i - start_cs)) << "\"";
                                    first_cs = false;
                                }
                                start_cs = i + 1;
                            }
                        }
                        json << "]";
                    }
                }
                json << "\n    }";
            }
            json << "],\n";
        }

        // Separate annotations by type: transcript, regulatory, intergenic
        std::vector<const VariantAnnotation*> transcript_anns;
        std::vector<const VariantAnnotation*> regulatory_anns;
        std::vector<const VariantAnnotation*> intergenic_anns;
        for (const auto& ann : buffered_annotations_) {
            bool is_intergenic = false;
            bool is_regulatory = false;
            for (const auto& csq : ann.consequences) {
                if (csq == ConsequenceType::INTERGENIC_VARIANT) is_intergenic = true;
                if (csq == ConsequenceType::REGULATORY_REGION_VARIANT ||
                    csq == ConsequenceType::REGULATORY_REGION_ABLATION ||
                    csq == ConsequenceType::REGULATORY_REGION_AMPLIFICATION ||
                    csq == ConsequenceType::TF_BINDING_SITE_VARIANT ||
                    csq == ConsequenceType::TFBS_ABLATION ||
                    csq == ConsequenceType::TFBS_AMPLIFICATION) is_regulatory = true;
            }
            if (is_intergenic) intergenic_anns.push_back(&ann);
            else if (is_regulatory && ann.transcript_id.empty()) regulatory_anns.push_back(&ann);
            else transcript_anns.push_back(&ann);
        }

        // transcript_consequences array
        if (!transcript_anns.empty()) {
            json << "    \"transcript_consequences\": [\n";
            bool first_tc = true;
            for (const auto* ann_ptr : transcript_anns) {
                const auto& ann = *ann_ptr;
                if (!first_tc) json << ",\n";
                first_tc = false;
                write_transcript_consequence(json, ann);
            }
            json << "\n    ]";
            if (!regulatory_anns.empty() || !intergenic_anns.empty()) json << ",";
            json << "\n";
        }

        // regulatory_feature_consequences array
        if (!regulatory_anns.empty()) {
            json << "    \"regulatory_feature_consequences\": [\n";
            bool first_rc = true;
            for (const auto* ann_ptr : regulatory_anns) {
                const auto& ann = *ann_ptr;
                if (!first_rc) json << ",\n";
                first_rc = false;
                json << "      {\n";
                json << "        \"consequence_terms\": [";
                bool first_csq = true;
                for (const auto& csq : ann.consequences) {
                    if (!first_csq) json << ", ";
                    json << "\"" << (term_style_ == "display" ? consequence_to_display_term(csq) : consequence_to_string(csq)) << "\"";
                    first_csq = false;
                }
                json << "],\n";
                json << "        \"impact\": \"" << impact_to_string(ann.impact) << "\",\n";
                auto reg_it = ann.custom_annotations.find("regulatory:feature_type");
                if (reg_it != ann.custom_annotations.end()) {
                    json << "        \"biotype\": \"" << escape_json(reg_it->second) << "\",\n";
                }
                json << "        \"variant_allele\": \"" << escape_json(ann.display_alt.empty() ? ann.alt_allele : ann.display_alt) << "\",\n";
                json << "        \"regulatory_feature_id\": \"" << escape_json(ann.gene_id) << "\"\n";
                json << "      }";
            }
            json << "\n    ]";
            if (!intergenic_anns.empty()) json << ",";
            json << "\n";
        }

        // intergenic_consequences array
        if (!intergenic_anns.empty()) {
            json << "    \"intergenic_consequences\": [\n";
            bool first_ic = true;
            for (const auto* ann_ptr : intergenic_anns) {
                const auto& ann = *ann_ptr;
                if (!first_ic) json << ",\n";
                first_ic = false;
                json << "      {\n";
                json << "        \"consequence_terms\": [";
                {
                    bool first_csq = true;
                    for (const auto& csq : ann.consequences) {
                        if (!first_csq) json << ", ";
                        json << "\"" << (term_style_ == "display" ? consequence_to_display_term(csq) : consequence_to_string(csq)) << "\"";
                        first_csq = false;
                    }
                }
                json << "],\n";
                json << "        \"impact\": \"" << impact_to_string(ann.impact) << "\",\n";
                json << "        \"variant_allele\": \"" << escape_json(ann.display_alt.empty() ? ann.alt_allele : ann.display_alt) << "\"\n";
                json << "      }";
            }
            json << "\n    ]\n";
        }

        json << "  }";

        write_string(json.str());
        buffered_annotations_.clear();
        current_variant_key_.clear();
    }

    void write_transcript_consequence(std::ostringstream& json, const VariantAnnotation& ann) {
        json << "      {\n";
        json << "        \"gene_id\": \"" << escape_json(ann.gene_id) << "\",\n";
        json << "        \"gene_symbol\": \"" << escape_json(ann.gene_symbol) << "\",\n";
        {
            auto ss_it = ann.custom_annotations.find("SYMBOL_SOURCE");
            if (ss_it != ann.custom_annotations.end() && !ss_it->second.empty()) {
                json << "        \"gene_symbol_source\": \"" << escape_json(ss_it->second) << "\",\n";
            }
            auto hgnc_it = ann.custom_annotations.find("HGNC_ID");
            if (hgnc_it != ann.custom_annotations.end() && !hgnc_it->second.empty()) {
                json << "        \"hgnc_id\": \"" << escape_json(hgnc_it->second) << "\",\n";
            }
        }
        json << "        \"transcript_id\": \"" << escape_json(ann.transcript_id) << "\",\n";
        if (!ann.source.empty()) {
            json << "        \"source\": \"" << escape_json(ann.source) << "\",\n";
        }
        json << "        \"biotype\": \"" << escape_json(ann.biotype) << "\",\n";

        // Perl VEP uses "canonical": 1 (integer) not boolean
        if (ann.is_canonical) {
            json << "        \"canonical\": 1,\n";
        }

        // Variant allele at transcript level (Perl VEP uses Ensembl display format)
        json << "        \"variant_allele\": \"" << escape_json(ann.display_alt.empty() ? ann.alt_allele : ann.display_alt) << "\",\n";

        // consequence_terms as array of strings
        json << "        \"consequence_terms\": [";
        bool first_csq = true;
        for (const auto& csq : ann.consequences) {
            if (!first_csq) json << ", ";
            json << "\"" << (term_style_ == "display" ? consequence_to_display_term(csq) : consequence_to_string(csq)) << "\"";
            first_csq = false;
        }
        json << "],\n";

        json << "        \"impact\": \"" << impact_to_string(ann.impact) << "\",\n";

        if (ann.strand != '\0') {
            json << "        \"strand\": " << (ann.strand == '+' ? "1" : "-1") << ",\n";
        }
        if (ann.distance > 0) {
            json << "        \"distance\": " << ann.distance << ",\n";
        }

        // Position details - Perl VEP uses start/end pairs
        if (ann.cdna_position > 0) {
            int cdna_e = (ann.cdna_end > 0) ? ann.cdna_end : ann.cdna_position;
            json << "        \"cdna_start\": " << ann.cdna_position << ",\n";
            json << "        \"cdna_end\": " << cdna_e << ",\n";
        }
        if (ann.cds_position > 0) {
            int cds_e = (ann.cds_end > 0) ? ann.cds_end : ann.cds_position;
            json << "        \"cds_start\": " << ann.cds_position << ",\n";
            json << "        \"cds_end\": " << cds_e << ",\n";
        }
        if (ann.protein_position > 0) {
            int prot_e = (ann.protein_end > 0) ? ann.protein_end : ann.protein_position;
            json << "        \"protein_start\": " << ann.protein_position << ",\n";
            json << "        \"protein_end\": " << prot_e << ",\n";
        }

        // Exon/intron
        if (ann.exon_number > 0) {
            json << "        \"exon\": \"" << ann.exon_number << "/" << ann.total_exons << "\",\n";
        }
        if (ann.intron_number > 0) {
            json << "        \"intron\": \"" << ann.intron_number << "/" << ann.total_introns << "\",\n";
        }

        // Sequence changes
        if (!ann.amino_acids.empty()) {
            json << "        \"amino_acids\": \"" << escape_json(ann.amino_acids) << "\",\n";
        }
        if (!ann.codons.empty()) {
            json << "        \"codons\": \"" << escape_json(ann.codons) << "\",\n";
        }
        if (!ann.hgvsc.empty()) {
            json << "        \"hgvsc\": \"" << escape_json(ann.hgvsc) << "\",\n";
        }
        if (!ann.hgvsp.empty()) {
            json << "        \"hgvsp\": \"" << escape_json(ann.hgvsp) << "\",\n";
        }
        if (!ann.hgvsg.empty()) {
            json << "        \"hgvsg\": \"" << escape_json(ann.hgvsg) << "\",\n";
        }

        // FLAGS as array (Perl VEP outputs array of flag strings)
        {
            auto flags_it = ann.custom_annotations.find("FLAGS");
            if (flags_it != ann.custom_annotations.end() && !flags_it->second.empty()) {
                json << "        \"flags\": [";
                // Split by comma
                std::string fl = flags_it->second;
                bool first_f = true;
                size_t start_f = 0;
                for (size_t i = 0; i <= fl.size(); ++i) {
                    if (i == fl.size() || fl[i] == ',') {
                        if (i > start_f) {
                            if (!first_f) json << ", ";
                            json << "\"" << escape_json(fl.substr(start_f, i - start_f)) << "\"";
                            first_f = false;
                        }
                        start_f = i + 1;
                    }
                }
                json << "],\n";
            }
        }

        // SIFT as separate prediction + score (Perl VEP format)
        {
            auto sift_it = ann.custom_annotations.find("SIFT");
            if (sift_it != ann.custom_annotations.end() && !sift_it->second.empty()) {
                // Parse "prediction(score)" format
                std::string val = sift_it->second;
                size_t paren = val.find('(');
                if (paren != std::string::npos && val.back() == ')') {
                    std::string pred = val.substr(0, paren);
                    std::string score = val.substr(paren + 1, val.size() - paren - 2);
                    json << "        \"sift_prediction\": \"" << escape_json(pred) << "\",\n";
                    json << "        \"sift_score\": " << score << ",\n";
                } else {
                    json << "        \"sift_prediction\": \"" << escape_json(val) << "\",\n";
                }
            }
        }

        // PolyPhen as separate prediction + score (Perl VEP format)
        {
            auto pp_it = ann.custom_annotations.find("PolyPhen");
            if (pp_it != ann.custom_annotations.end() && !pp_it->second.empty()) {
                std::string val = pp_it->second;
                size_t paren = val.find('(');
                if (paren != std::string::npos && val.back() == ')') {
                    std::string pred = val.substr(0, paren);
                    std::string score = val.substr(paren + 1, val.size() - paren - 2);
                    json << "        \"polyphen_prediction\": \"" << escape_json(pred) << "\",\n";
                    json << "        \"polyphen_score\": " << score << ",\n";
                } else {
                    json << "        \"polyphen_prediction\": \"" << escape_json(val) << "\",\n";
                }
            }
        }

        // Domains as array of objects (Perl VEP format: [{db: "pfam", name: "PF00001"}])
        {
            bool has_domains = false;
            std::ostringstream domains;
            domains << "        \"domains\": [";
            bool first_d = true;
            // Check for pfam domains
            auto pfam_id_it = ann.custom_annotations.find("pfam:domain_id");
            if (pfam_id_it != ann.custom_annotations.end() && !pfam_id_it->second.empty()) {
                has_domains = true;
                domains << "{\"db\": \"Pfam\", \"name\": \"" << escape_json(pfam_id_it->second) << "\"}";

                first_d = false;
            }
            // Check for interpro domains
            auto ip_id_it = ann.custom_annotations.find("interpro:domain_id");
            if (ip_id_it != ann.custom_annotations.end() && !ip_id_it->second.empty()) {
                has_domains = true;
                if (!first_d) domains << ", ";
                domains << "{\"db\": \"Interpro\", \"name\": \"" << escape_json(ip_id_it->second) << "\"}";

            }
            domains << "],\n";
            if (has_domains) {
                json << domains.str();
            }
        }

        // MANE, TSL, APPRIS, ENSP as explicit fields
        {
            auto mane_it = ann.custom_annotations.find("MANE_SELECT");
            if (mane_it != ann.custom_annotations.end() && !mane_it->second.empty()) {
                json << "        \"mane_select\": \"" << escape_json(mane_it->second) << "\",\n";
            }
            auto mane_plus_it = ann.custom_annotations.find("MANE_PLUS_CLINICAL");
            if (mane_plus_it != ann.custom_annotations.end() && !mane_plus_it->second.empty()) {
                json << "        \"mane_plus_clinical\": \"" << escape_json(mane_plus_it->second) << "\",\n";
            }
            auto tsl_it = ann.custom_annotations.find("TSL");
            if (tsl_it != ann.custom_annotations.end() && !tsl_it->second.empty()) {
                json << "        \"tsl\": " << tsl_it->second << ",\n";
            }
            auto appris_it = ann.custom_annotations.find("APPRIS");
            if (appris_it != ann.custom_annotations.end() && !appris_it->second.empty()) {
                json << "        \"appris\": \"" << escape_json(appris_it->second) << "\",\n";
            }
            auto ensp_it = ann.custom_annotations.find("ENSP");
            if (ensp_it != ann.custom_annotations.end() && !ensp_it->second.empty()) {
                json << "        \"protein_id\": \"" << escape_json(ensp_it->second) << "\",\n";
            }
        }

        // Remaining custom annotations (skip internal/already-output fields)
        static const std::set<std::string> skip_custom = {
            // Already output as explicit JSON fields
            "SYMBOL_SOURCE", "HGNC_ID", "FLAGS", "SIFT", "PolyPhen",
            "SIFT_score", "SIFT_pred", "Polyphen2_HDIV_score", "Polyphen2_HDIV_pred",
            "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred",
            "pfam:domain_id", "pfam:domain_name", "interpro:domain_id", "interpro:domain_name",
            "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", "APPRIS", "ENSP",
            "HGVS_OFFSET", "regulatory_type", "regulatory:feature_type", "CLIN_SIG",
            "_colocated:ID", "_colocated:CLNSIG", "EXISTING_CSQ",
            // Variant-level fields (output at variant level, not per-transcript)
            "VARIANT_CLASS", "MINIMISED",
            // Fields already on the VariantAnnotation struct (would be duplicated)
            "CANONICAL", "BIOTYPE", "STRAND", "_consequences"
        };
        for (const auto& pair : ann.custom_annotations) {
            if (skip_custom.count(pair.first) || pair.first.substr(0, 11) == "_colocated:") continue;
            if (!pair.second.empty()) {
                // Perl VEP lowercases all JSON keys
                std::string key = pair.first;
                for (auto& c : key) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
                // Try to output numeric values as numbers
                bool is_numeric = !pair.second.empty();
                bool has_dot = false;
                for (size_t i = 0; i < pair.second.size(); ++i) {
                    char c = pair.second[i];
                    if (c == '-' && i == 0) continue;
                    if (c == '.' && !has_dot) { has_dot = true; continue; }
                    if (!std::isdigit(static_cast<unsigned char>(c))) { is_numeric = false; break; }
                }
                if (is_numeric && !pair.second.empty() && pair.second != "-" && pair.second != ".") {
                    json << "        \"" << escape_json(key) << "\": " << pair.second << ",\n";
                } else {
                    json << "        \"" << escape_json(key) << "\": \"" << escape_json(pair.second) << "\",\n";
                }
            }
        }

        // Remove trailing comma from last field and close object
        {
            std::string s = json.str();
            if (s.size() >= 2 && s.back() == '\n' && s[s.size()-2] == ',') {
                s.resize(s.size() - 2);
                s += '\n';
            }
            json.str(std::move(s));
            json.seekp(0, std::ios_base::end);
        }
        json << "      }";
    }

    void write_string(const std::string& s) {
        if (use_stdout_) {
            std::cout << s;
        } else if (compress_ && gz_file_) {
            gzwrite(gz_file_, s.c_str(), static_cast<unsigned int>(s.size()));
        } else {
            output_ << s;
        }
    }

    static std::string escape_json(const std::string& s) {
        std::string result;
        result.reserve(s.size());
        for (char c : s) {
            switch (c) {
                case '"': result += "\\\""; break;
                case '\\': result += "\\\\"; break;
                case '\b': result += "\\b"; break;
                case '\f': result += "\\f"; break;
                case '\n': result += "\\n"; break;
                case '\r': result += "\\r"; break;
                case '\t': result += "\\t"; break;
                default: result += c; break;
            }
        }
        return result;
    }
};

/**
 * VCF output writer - writes VCF with CSQ INFO field
 */
class VCFWriter : public OutputWriter {
public:
    explicit VCFWriter(const std::string& output_path, bool compress = false)
        : output_path_(output_path), compress_(compress), gz_file_(nullptr),
          info_field_name_("CSQ"),
          use_stdout_(output_path.empty() || output_path == "-" || output_path == "STDOUT") {

        if (use_stdout_) {
            compress_ = false;
        } else if (compress_ || ends_with_gz(output_path_)) {
            compress_ = true;
            gz_file_ = gzopen(output_path_.c_str(), "wb");
            if (!gz_file_) {
                throw std::runtime_error("Cannot open output file: " + output_path_);
            }
        } else {
            output_.open(output_path_);
            if (!output_.is_open()) {
                throw std::runtime_error("Cannot open output file: " + output_path_);
            }
        }
    }

    ~VCFWriter() override {
        close();
    }

    void set_info_field_name(const std::string& name) { info_field_name_ = name; }
    void set_no_escape(bool v) { no_escape_ = v; }
    void set_keep_csq(bool v) { keep_csq_ = v; }

    /** Add a passthrough VCF header line (## meta-information lines from input VCF) */
    void add_passthrough_header(const std::string& line) { passthrough_headers_.push_back(line); }

    /** Set the original column header line from input VCF (used for sample column names) */
    void set_column_header(const std::string& line) { original_column_header_ = line; }

    /**
     * Set custom field order for CSQ output.
     * Fields: Allele, Consequence, IMPACT, SYMBOL, Gene, Feature_type, Feature,
     * BIOTYPE, EXON, INTRON, HGVSc, HGVSp, cDNA_position, CDS_position,
     * Protein_position, Amino_acids, Codons, CANONICAL, plus any custom columns.
     */
    void set_field_order(const std::vector<std::string>& fields) { field_order_ = fields; }

    void write_header(const std::vector<std::string>& custom_columns) override {
        custom_columns_ = custom_columns;
        if (skip_header_) return;

        std::ostringstream header;

        // Write passthrough headers from input VCF (## meta-lines)
        if (!passthrough_headers_.empty()) {
            for (const auto& h : passthrough_headers_) {
                // Skip ##INFO CSQ if we're adding our own
                if (h.find("##INFO=<ID=" + info_field_name_ + ",") == std::string::npos) {
                    header << h << "\n";
                }
            }
        } else {
            header << "##fileformat=VCFv4.2\n";
        }

        // Add CSQ INFO header
        header << "##INFO=<ID=" << info_field_name_ << ",Number=.,Type=String,Description=\"Consequence annotations from VEP. "
               << "Format: ";

        if (!field_order_.empty()) {
            for (size_t i = 0; i < field_order_.size(); ++i) {
                if (i > 0) header << "|";
                header << field_order_[i];
            }
        } else {
            header << "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|"
                   << "EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|"
                   << "Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS";
            for (const auto& col : custom_columns_) {
                header << "|" << col;
            }
        }
        header << "\">\n";

        // Add filter header if not already present from passthrough
        if (passthrough_headers_.empty()) {
            header << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
        }

        // Column header - use original if available (preserves sample columns)
        if (!original_column_header_.empty()) {
            header << original_column_header_ << "\n";
        } else {
            header << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        }

        write_string(header.str());
    }

    void write_annotation(const VariantAnnotation& ann,
                          const std::vector<std::string>& /*custom_columns*/) override {
        stats_.add(ann);
        current_annotations_.push_back(ann);
    }

    void write_annotations(const std::vector<VariantAnnotation>& anns,
                           const std::vector<std::string>& custom_columns) override {
        for (const auto& ann : anns) {
            write_annotation(ann, custom_columns);
        }
    }

    // Call this after adding all annotations for a variant position
    void flush_variant() {
        if (current_annotations_.empty()) return;

        const auto& first = current_annotations_[0];

        // Use preserved VCF fields if available, otherwise fall back to defaults
        std::string out_id = first.vcf_id.empty() ? "." : first.vcf_id;
        std::string out_qual = first.vcf_qual.empty() ? "." : first.vcf_qual;
        std::string out_filter = first.vcf_filter.empty() ? "." : first.vcf_filter;

        // Build INFO field: preserve original INFO fields and append/replace CSQ
        std::string info_str;
        if (!first.vcf_info.empty() && first.vcf_info != ".") {
            // Rebuild INFO by removing any existing CSQ field and appending new one
            std::istringstream info_ss(first.vcf_info);
            std::string info_token;
            bool first_info = true;
            while (std::getline(info_ss, info_token, ';')) {
                // Skip the existing CSQ/info_field_name_ entry
                size_t eq = info_token.find('=');
                std::string key = (eq != std::string::npos) ? info_token.substr(0, eq) : info_token;
                if (key == info_field_name_) continue;
                if (!first_info) info_str += ";";
                info_str += info_token;
                first_info = false;
            }
            if (!info_str.empty()) info_str += ";";
        }
        info_str += info_field_name_ + "=";

        // Use original VCF REF/ALT if available (preserves multi-allele ALT)
        std::string vcf_ref = !first.vcf_ref.empty() ? first.vcf_ref :
                              (first.ref_allele.empty() ? "-" : first.ref_allele);
        std::string vcf_alt = !first.vcf_alt.empty() ? first.vcf_alt :
                              (first.alt_allele.empty() ? "-" : first.alt_allele);

        std::ostringstream line;
        line << first.chromosome << "\t"
             << first.position << "\t"
             << out_id << "\t"
             << vcf_ref << "\t"
             << vcf_alt << "\t"
             << out_qual << "\t"
             << out_filter << "\t"
             << info_str;

        // Prepend existing CSQ entries when --keep-csq is active
        bool first_csq = true;
        if (keep_csq_) {
            auto ec_it = first.custom_annotations.find("EXISTING_CSQ");
            if (ec_it != first.custom_annotations.end() && !ec_it->second.empty()) {
                line << ec_it->second;
                first_csq = false;
            }
        }
        for (const auto& ann : current_annotations_) {
            if (!first_csq) line << ",";
            line << format_csq(ann);
            first_csq = false;
        }

        // Append sample columns (FORMAT + samples) if present
        if (!first.vcf_sample_columns.empty()) {
            line << "\t" << first.vcf_sample_columns;
        }

        line << "\n";
        write_string(line.str());

        current_annotations_.clear();
    }

    void write_footer() override {
        flush_variant();  // Ensure last variant is written
    }

    void close() override {
        flush_variant();
        if (gz_file_) {
            gzclose(gz_file_);
            gz_file_ = nullptr;
        }
        if (output_.is_open()) {
            output_.close();
        }
    }

private:
    std::string output_path_;
    bool compress_;
    std::ofstream output_;
    gzFile gz_file_;
    std::string info_field_name_;
    bool no_escape_ = false;
    bool keep_csq_ = false;
    bool use_stdout_;
    std::vector<std::string> field_order_;
    std::vector<std::string> custom_columns_;
    std::vector<VariantAnnotation> current_annotations_;
    std::vector<std::string> passthrough_headers_;
    std::string original_column_header_;

    void write_string(const std::string& s) {
        if (use_stdout_) {
            std::cout << s;
        } else if (compress_ && gz_file_) {
            gzwrite(gz_file_, s.c_str(), static_cast<unsigned int>(s.size()));
        } else {
            output_ << s;
        }
    }

    // Get a named field value from annotation
    std::string get_csq_field(const VariantAnnotation& ann, const std::string& field) const {
        if (field == "Allele") return escape_vcf(ann.display_alt.empty() ? ann.alt_allele : ann.display_alt);
        if (field == "Consequence") return escape_vcf(format_consequence(ann));
        if (field == "IMPACT") return escape_vcf(impact_to_string(ann.impact));
        if (field == "SYMBOL") return escape_vcf(ann.gene_symbol);
        if (field == "Gene") return escape_vcf(ann.gene_id);
        if (field == "Feature_type") return escape_vcf(ann.feature_type);
        if (field == "Feature") return escape_vcf(ann.transcript_id);
        if (field == "BIOTYPE") return escape_vcf(ann.biotype);
        if (field == "EXON") return ann.exon_number > 0 ? std::to_string(ann.exon_number) + "/" + std::to_string(ann.total_exons) : "";
        if (field == "INTRON") return ann.intron_number > 0 ? std::to_string(ann.intron_number) + "/" + std::to_string(ann.total_introns) : "";
        if (field == "HGVSc") return escape_vcf(ann.hgvsc);
        if (field == "HGVSp") return escape_vcf(ann.hgvsp);
        if (field == "HGVSg") return escape_vcf(ann.hgvsg);
        if (field == "cDNA_position") return format_position_with_total(ann.cdna_position, ann.cdna_end, "TRANSCRIPT_LENGTH", ann, "");
        if (field == "CDS_position") return format_position_with_total(ann.cds_position, ann.cds_end, "CDS_LENGTH", ann, "");
        if (field == "Protein_position") return format_protein_position(ann, "");
        if (field == "Amino_acids") return escape_vcf(ann.amino_acids);
        if (field == "Codons") return escape_vcf(ann.codons);
        if (field == "Existing_variation") return escape_vcf(ann.existing_variation);
        if (field == "CANONICAL") return ann.is_canonical ? "YES" : "";
        if (field == "STRAND") return ann.strand != '\0' ? (ann.strand == '+' ? "1" : "-1") : "";
        if (field == "DISTANCE") return ann.distance > 0 ? std::to_string(ann.distance) : "";
        if (field == "SOURCE") return escape_vcf(ann.source);
        // Fall through to custom annotations
        auto it = ann.custom_annotations.find(field);
        if (it != ann.custom_annotations.end()) return escape_vcf(it->second);
        return "";
    }

    std::string format_csq(const VariantAnnotation& ann) {
        // Use custom field order if set via --fields
        if (!field_order_.empty()) {
            std::ostringstream csq;
            for (size_t i = 0; i < field_order_.size(); ++i) {
                if (i > 0) csq << "|";
                csq << get_csq_field(ann, field_order_[i]);
            }
            return csq.str();
        }

        // Default format
        std::ostringstream csq;

        csq << escape_vcf(ann.display_alt.empty() ? ann.alt_allele : ann.display_alt) << "|"
            << escape_vcf(format_consequence(ann)) << "|"
            << escape_vcf(impact_to_string(ann.impact)) << "|"
            << escape_vcf(ann.gene_symbol) << "|"
            << escape_vcf(ann.gene_id) << "|"
            << escape_vcf(ann.feature_type) << "|"
            << escape_vcf(ann.transcript_id) << "|"
            << escape_vcf(ann.biotype) << "|"
            << (ann.exon_number > 0 ? std::to_string(ann.exon_number) + "/" + std::to_string(ann.total_exons) : "") << "|"
            << (ann.intron_number > 0 ? std::to_string(ann.intron_number) + "/" + std::to_string(ann.total_introns) : "") << "|"
            << escape_vcf(ann.hgvsc) << "|"
            << escape_vcf(ann.hgvsp) << "|"
            << format_position_with_total(ann.cdna_position, ann.cdna_end, "TRANSCRIPT_LENGTH", ann, "") << "|"
            << format_position_with_total(ann.cds_position, ann.cds_end, "CDS_LENGTH", ann, "") << "|"
            << format_protein_position(ann, "") << "|"
            << escape_vcf(ann.amino_acids) << "|"
            << escape_vcf(ann.codons) << "|"
            << escape_vcf(ann.existing_variation) << "|"
            << (ann.distance > 0 ? std::to_string(ann.distance) : "") << "|"
            << (ann.strand != '\0' ? (ann.strand == '+' ? "1" : "-1") : "") << "|";
        // FLAGS from custom_annotations
        {
            auto flags_it = ann.custom_annotations.find("FLAGS");
            csq << (flags_it != ann.custom_annotations.end() ? escape_vcf(flags_it->second) : "");
        }

        for (const auto& col : custom_columns_) {
            csq << "|";
            auto it = ann.custom_annotations.find(col);
            if (it != ann.custom_annotations.end()) {
                csq << escape_vcf(it->second);
            }
        }

        return csq.str();
    }

    std::string escape_vcf(const std::string& s) const {
        if (no_escape_) return s;
        std::string result;
        result.reserve(s.size());
        for (char c : s) {
            switch (c) {
                case '%': result += "%25"; break;
                case '|': result += "&"; break;
                case ',': result += "&"; break;
                case ';': result += "%3B"; break;
                case '=': result += "%3D"; break;
                case ' ': result += "_"; break;
                case '\t': result += "_"; break;
                default: result += c; break;
            }
        }
        return result;
    }
};

/**
 * Factory function to create appropriate writer
 */
inline std::unique_ptr<OutputWriter> create_output_writer(
    const std::string& output_path,
    OutputFormat format,
    bool compress = false) {

    if (format == OutputFormat::JSON) {
        return std::make_unique<JSONWriter>(output_path, compress);
    } else if (format == OutputFormat::VCF) {
        return std::make_unique<VCFWriter>(output_path, compress);
    } else {
        return std::make_unique<TSVWriter>(output_path, compress);
    }
}

} // namespace vep

#endif // OUTPUT_WRITER_HPP
