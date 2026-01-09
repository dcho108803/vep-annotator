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
#include <map>
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

protected:
    AnnotationStats stats_;
};

// Helper to check if path ends with .gz
inline bool ends_with_gz(const std::string& path) {
    return path.size() > 3 && path.substr(path.size() - 3) == ".gz";
}

/**
 * TSV output writer (default format)
 */
class TSVWriter : public OutputWriter {
public:
    explicit TSVWriter(const std::string& output_path, bool compress = false)
        : output_path_(output_path), compress_(compress), gz_file_(nullptr) {

        if (compress_ || ends_with_gz(output_path_)) {
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
        std::ostringstream header;
        header << "CHROM\tPOS\tREF\tALT\tGENE\tTRANSCRIPT\tCONSEQUENCE\tIMPACT\t"
               << "CDS_POS\tPROTEIN_POS\tAMINO_ACIDS\tCODONS\tHGVSc\tHGVSp\t"
               << "EXON\tINTRON\tBIOTYPE\tCANONICAL";

        for (const auto& col : custom_columns) {
            header << "\t" << col;
        }
        header << "\n";

        write_string(header.str());
    }

    void write_annotation(const VariantAnnotation& ann,
                          const std::vector<std::string>& custom_columns) override {
        stats_.add(ann);

        std::ostringstream line;
        line << ann.chromosome << "\t"
             << ann.position << "\t"
             << ann.ref_allele << "\t"
             << ann.alt_allele << "\t"
             << ann.gene_symbol << "\t"
             << ann.transcript_id << "\t"
             << ann.get_consequence_string() << "\t"
             << impact_to_string(ann.impact) << "\t"
             << (ann.cds_position > 0 ? std::to_string(ann.cds_position) : "") << "\t"
             << (ann.protein_position > 0 ? std::to_string(ann.protein_position) : "") << "\t"
             << ann.amino_acids << "\t"
             << ann.codons << "\t"
             << ann.hgvsc << "\t"
             << ann.hgvsp << "\t"
             << (ann.exon_number > 0 ? std::to_string(ann.exon_number) : "") << "\t"
             << (ann.intron_number > 0 ? std::to_string(ann.intron_number) : "") << "\t"
             << ann.biotype << "\t"
             << (ann.is_canonical ? "YES" : "");

        for (const auto& col : custom_columns) {
            line << "\t";
            auto it = ann.custom_annotations.find(col);
            if (it != ann.custom_annotations.end()) {
                line << it->second;
            }
        }
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

    void write_string(const std::string& s) {
        if (compress_ && gz_file_) {
            gzwrite(gz_file_, s.c_str(), static_cast<unsigned int>(s.size()));
        } else {
            output_ << s;
        }
    }
};

/**
 * JSON output writer
 */
class JSONWriter : public OutputWriter {
public:
    explicit JSONWriter(const std::string& output_path, bool compress = false)
        : output_path_(output_path), compress_(compress), gz_file_(nullptr), first_record_(true) {

        if (compress_ || ends_with_gz(output_path_)) {
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

    void write_header(const std::vector<std::string>& /*custom_columns*/) override {
        write_string("[\n");
    }

    void write_annotation(const VariantAnnotation& ann,
                          const std::vector<std::string>& /*custom_columns*/) override {
        stats_.add(ann);

        if (!first_record_) {
            write_string(",\n");
        }
        first_record_ = false;

        std::ostringstream json;
        json << "  {\n";
        json << "    \"input\": \"" << escape_json(ann.input_variant) << "\",\n";
        json << "    \"chromosome\": \"" << escape_json(ann.chromosome) << "\",\n";
        json << "    \"position\": " << ann.position << ",\n";
        json << "    \"ref_allele\": \"" << escape_json(ann.ref_allele) << "\",\n";
        json << "    \"alt_allele\": \"" << escape_json(ann.alt_allele) << "\",\n";
        json << "    \"gene_symbol\": \"" << escape_json(ann.gene_symbol) << "\",\n";
        json << "    \"gene_id\": \"" << escape_json(ann.gene_id) << "\",\n";
        json << "    \"transcript_id\": \"" << escape_json(ann.transcript_id) << "\",\n";
        json << "    \"biotype\": \"" << escape_json(ann.biotype) << "\",\n";
        json << "    \"is_canonical\": " << (ann.is_canonical ? "true" : "false") << ",\n";

        // Consequences array
        json << "    \"consequences\": [";
        bool first_csq = true;
        for (const auto& csq : ann.consequences) {
            if (!first_csq) json << ", ";
            json << "\"" << consequence_to_string(csq) << "\"";
            first_csq = false;
        }
        json << "],\n";

        json << "    \"impact\": \"" << impact_to_string(ann.impact) << "\",\n";
        json << "    \"exon_number\": " << (ann.exon_number > 0 ? std::to_string(ann.exon_number) : "null") << ",\n";
        json << "    \"intron_number\": " << (ann.intron_number > 0 ? std::to_string(ann.intron_number) : "null") << ",\n";
        json << "    \"cds_position\": " << (ann.cds_position > 0 ? std::to_string(ann.cds_position) : "null") << ",\n";
        json << "    \"protein_position\": " << (ann.protein_position > 0 ? std::to_string(ann.protein_position) : "null") << ",\n";
        json << "    \"codons\": \"" << escape_json(ann.codons) << "\",\n";
        json << "    \"amino_acids\": \"" << escape_json(ann.amino_acids) << "\",\n";
        json << "    \"hgvsc\": \"" << escape_json(ann.hgvsc) << "\",\n";
        json << "    \"hgvsp\": \"" << escape_json(ann.hgvsp) << "\"";

        // Custom annotations
        if (!ann.custom_annotations.empty()) {
            json << ",\n    \"annotations\": {\n";
            bool first_ann = true;
            for (const auto& pair : ann.custom_annotations) {
                if (!first_ann) json << ",\n";
                json << "      \"" << escape_json(pair.first) << "\": \"" << escape_json(pair.second) << "\"";
                first_ann = false;
            }
            json << "\n    }";
        }

        json << "\n  }";

        write_string(json.str());
    }

    void write_annotations(const std::vector<VariantAnnotation>& anns,
                           const std::vector<std::string>& custom_columns) override {
        for (const auto& ann : anns) {
            write_annotation(ann, custom_columns);
        }
    }

    void write_footer() override {
        write_string("\n]\n");
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
    bool first_record_;

    void write_string(const std::string& s) {
        if (compress_ && gz_file_) {
            gzwrite(gz_file_, s.c_str(), static_cast<unsigned int>(s.size()));
        } else {
            output_ << s;
        }
    }

    static std::string escape_json(const std::string& s) {
        std::ostringstream oss;
        for (char c : s) {
            switch (c) {
                case '"': oss << "\\\""; break;
                case '\\': oss << "\\\\"; break;
                case '\b': oss << "\\b"; break;
                case '\f': oss << "\\f"; break;
                case '\n': oss << "\\n"; break;
                case '\r': oss << "\\r"; break;
                case '\t': oss << "\\t"; break;
                default: oss << c; break;
            }
        }
        return oss.str();
    }
};

/**
 * VCF output writer - writes VCF with CSQ INFO field
 */
class VCFWriter : public OutputWriter {
public:
    explicit VCFWriter(const std::string& output_path, bool compress = false)
        : output_path_(output_path), compress_(compress), gz_file_(nullptr) {

        if (compress_ || ends_with_gz(output_path_)) {
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

    void write_header(const std::vector<std::string>& custom_columns) override {
        custom_columns_ = custom_columns;

        std::ostringstream header;
        header << "##fileformat=VCFv4.2\n";
        header << "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from VEP. "
               << "Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|"
               << "EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|"
               << "Amino_acids|Codons|CANONICAL";

        for (const auto& col : custom_columns_) {
            header << "|" << col;
        }
        header << "\">\n";

        // Add filter header
        header << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";

        // Column header
        header << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

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

        std::ostringstream line;
        line << first.chromosome << "\t"
             << first.position << "\t"
             << "." << "\t"
             << first.ref_allele << "\t"
             << first.alt_allele << "\t"
             << "." << "\t"
             << "PASS" << "\t"
             << "CSQ=";

        bool first_csq = true;
        for (const auto& ann : current_annotations_) {
            if (!first_csq) line << ",";
            line << format_csq(ann);
            first_csq = false;
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
    std::vector<std::string> custom_columns_;
    std::vector<VariantAnnotation> current_annotations_;

    void write_string(const std::string& s) {
        if (compress_ && gz_file_) {
            gzwrite(gz_file_, s.c_str(), static_cast<unsigned int>(s.size()));
        } else {
            output_ << s;
        }
    }

    std::string format_csq(const VariantAnnotation& ann) {
        std::ostringstream csq;

        csq << escape_vcf(ann.alt_allele) << "|"
            << escape_vcf(ann.get_consequence_string()) << "|"
            << escape_vcf(impact_to_string(ann.impact)) << "|"
            << escape_vcf(ann.gene_symbol) << "|"
            << escape_vcf(ann.gene_id) << "|"
            << "Transcript" << "|"
            << escape_vcf(ann.transcript_id) << "|"
            << escape_vcf(ann.biotype) << "|"
            << (ann.exon_number > 0 ? std::to_string(ann.exon_number) : "") << "|"
            << (ann.intron_number > 0 ? std::to_string(ann.intron_number) : "") << "|"
            << escape_vcf(ann.hgvsc) << "|"
            << escape_vcf(ann.hgvsp) << "|"
            << "|"  // cDNA_position (not implemented yet)
            << (ann.cds_position > 0 ? std::to_string(ann.cds_position) : "") << "|"
            << (ann.protein_position > 0 ? std::to_string(ann.protein_position) : "") << "|"
            << escape_vcf(ann.amino_acids) << "|"
            << escape_vcf(ann.codons) << "|"
            << (ann.is_canonical ? "YES" : "");

        for (const auto& col : custom_columns_) {
            csq << "|";
            auto it = ann.custom_annotations.find(col);
            if (it != ann.custom_annotations.end()) {
                csq << escape_vcf(it->second);
            }
        }

        return csq.str();
    }

    static std::string escape_vcf(const std::string& s) {
        std::ostringstream oss;
        for (char c : s) {
            switch (c) {
                case '|': oss << "%7C"; break;
                case ';': oss << "%3B"; break;
                case '=': oss << "%3D"; break;
                case ',': oss << "%2C"; break;
                case ' ': oss << "%20"; break;
                case '\t': oss << "%09"; break;
                default: oss << c; break;
            }
        }
        return oss.str();
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
        return std::unique_ptr<OutputWriter>(new JSONWriter(output_path, compress));
    } else if (format == OutputFormat::VCF) {
        return std::unique_ptr<OutputWriter>(new VCFWriter(output_path, compress));
    } else {
        return std::unique_ptr<OutputWriter>(new TSVWriter(output_path, compress));
    }
}

} // namespace vep

#endif // OUTPUT_WRITER_HPP
