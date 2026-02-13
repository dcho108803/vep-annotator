/**
 * dbscSNV Annotation Source
 *
 * Provides splice predictions from dbscSNV database.
 * dbscSNV contains all potential human SNVs within splicing consensus regions.
 */

#include "annotation_source.hpp"
#include "file_parsers.hpp"
#include "vep_annotator.hpp"
#include <sstream>
#include <algorithm>

namespace vep {

/**
 * dbscSNV Annotation Source
 *
 * Queries dbscSNV for ada_score and rf_score at each splice site.
 */
class DbscSNVSource : public VariantAnnotationSource {
public:
    explicit DbscSNVSource(const std::string& path)
        : path_(path) {}

    std::string name() const override { return "dbscsnv"; }
    std::string type() const override { return "splice"; }
    std::string description() const override {
        return "dbscSNV splice site predictions (ada_score, rf_score)";
    }

    bool is_ready() const override { return reader_ != nullptr && reader_->is_valid(); }

    void initialize() override {
        std::lock_guard<std::recursive_mutex> lock(mutex_);
        if (reader_) return;

        log(LogLevel::INFO, "Loading dbscSNV from: " + path_);

        reader_ = std::make_unique<TabixTSVReader>(path_, 0, 1);

        if (!reader_->is_valid()) {
            log(LogLevel::ERROR, "Failed to open dbscSNV file: " + path_);
            reader_.reset();
            return;
        }

        log(LogLevel::INFO, "dbscSNV loaded successfully");
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

        (void)transcript;

        auto records = reader_->query(chrom, pos);

        for (const auto& record : records) {
            // Check for ref/alt match
            auto ref_it = record.find("ref");
            auto alt_it = record.find("alt");

            if (ref_it != record.end() && ref_it->second != ref) continue;
            if (alt_it != record.end() && alt_it->second != alt) continue;

            // Extract scores
            auto ada_it = record.find("ada_score");
            auto rf_it = record.find("rf_score");

            if (ada_it != record.end() && ada_it->second != ".") {
                annotations["dbscsnv:ada_score"] = ada_it->second;
            }
            if (rf_it != record.end() && rf_it->second != ".") {
                annotations["dbscsnv:rf_score"] = rf_it->second;
            }

            // Compute max score for interpretation
            double ada = 0, rf = 0;
            try {
                if (ada_it != record.end() && ada_it->second != ".") {
                    ada = std::stod(ada_it->second);
                }
                if (rf_it != record.end() && rf_it->second != ".") {
                    rf = std::stod(rf_it->second);
                }
            } catch (...) {}

            double max_score = std::max(ada, rf);
            if (max_score > 0) {
                annotations["dbscsnv:max_score"] = std::to_string(max_score);
                std::string pred = interpret_score(max_score);
                if (!pred.empty()) {
                    annotations["dbscsnv:prediction"] = pred;
                }
            }

            return;
        }
    }

    std::vector<std::string> get_fields() const override {
        return {
            "dbscsnv:ada_score",   // Adaptive boosting score
            "dbscsnv:rf_score",    // Random forest score
            "dbscsnv:max_score",   // Maximum of both scores
            "dbscsnv:prediction"   // Interpretation
        };
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
            auto alt_it = record.find("alt");

            if (ref_it != record.end() && ref_it->second != ref) continue;
            if (alt_it != record.end() && alt_it->second != alt) continue;

            auto ada_it = record.find("ada_score");
            auto rf_it = record.find("rf_score");

            if (ada_it != record.end() && ada_it->second != ".") {
                result["dbscsnv:ada_score"] = ada_it->second;
            }
            if (rf_it != record.end() && rf_it->second != ".") {
                result["dbscsnv:rf_score"] = rf_it->second;
            }
            break;
        }

        return result;
    }

    std::string get_data_path() const override { return path_; }

    /**
     * Interpret dbscSNV score
     * @return "splice_altering" (>= 0.6), "possible" (>= 0.4), or ""
     */
    static std::string interpret_score(double score) {
        if (score >= 0.6) return "splice_altering";
        if (score >= 0.4) return "possible";
        return "";
    }

private:
    std::string path_;
    std::unique_ptr<TabixTSVReader> reader_;
};

/**
 * Factory function to create dbscSNV source
 */
std::shared_ptr<AnnotationSource> create_dbscsnv_source(const std::string& path) {
    return std::make_shared<DbscSNVSource>(path);
}

} // namespace vep
