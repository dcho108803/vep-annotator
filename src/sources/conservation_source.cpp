/**
 * Conservation Score Annotation Sources
 *
 * Provides conservation scores from bigWig files:
 * - PhyloP: Phylogenetic p-values (measures evolutionary conservation)
 * - PhastCons: Probability of negative selection
 * - GERP++: Genomic Evolutionary Rate Profiling
 */

#include "annotation_source.hpp"
#include "file_parsers.hpp"
#include "vep_annotator.hpp"
#include <cmath>
#include <sstream>

namespace vep {

/**
 * Generic bigWig-based conservation score source
 */
class BigWigScoreSource : public ScoreAnnotationSource {
public:
    BigWigScoreSource(const std::string& path,
                      const std::string& source_name,
                      const std::string& field_name,
                      const std::string& desc)
        : path_(path),
          source_name_(source_name),
          field_name_(field_name),
          description_(desc) {}

    std::string name() const override { return source_name_; }
    std::string type() const override { return "conservation"; }
    std::string description() const override { return description_; }

    bool is_ready() const override {
        return reader_ != nullptr && reader_->is_valid();
    }

    void initialize() override {
        std::lock_guard<std::recursive_mutex> lock(mutex_);
        if (reader_) return;

        log(LogLevel::INFO, "Loading " + source_name_ + " from: " + path_);

        reader_ = std::make_unique<BigWigReader>(path_);

        if (!reader_->is_valid()) {
            log(LogLevel::ERROR, "Failed to open bigWig file: " + path_);
            reader_.reset();
            return;
        }

        log(LogLevel::INFO, source_name_ + " loaded with " +
                           std::to_string(reader_->get_chromosomes().size()) + " chromosomes");
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

        // For SNPs, get single position score
        if (ref.length() == 1 && alt.length() == 1) {
            auto score = reader_->get_value(chrom, pos);
            if (score.has_value()) {
                annotations[source_name_ + ":" + field_name_] =
                    format_score(score.value());
            }
        } else {
            // For indels, get mean score across affected region
            int start = pos;
            int end = pos + static_cast<int>(ref.length()) - 1;

            auto mean = reader_->get_mean(chrom, start, end);
            if (mean.has_value()) {
                annotations[source_name_ + ":" + field_name_] =
                    format_score(mean.value());
            }
        }
    }

    std::vector<std::string> get_fields() const override {
        return { source_name_ + ":" + field_name_ };
    }

    std::optional<double> get_score(const std::string& chrom, int pos) const override {
        if (!reader_) return std::nullopt;
        return reader_->get_value(chrom, pos);
    }

    std::vector<double> get_scores(
        const std::string& chrom,
        int start,
        int end
    ) const override {
        if (!reader_) return {};
        return reader_->get_values(chrom, start, end);
    }

    std::string get_data_path() const override { return path_; }

    bool is_thread_safe() const override { return true; }

protected:
    std::string path_;
    std::string source_name_;
    std::string field_name_;
    std::string description_;
    std::unique_ptr<BigWigReader> reader_;

    static std::string format_score(double score) {
        if (std::isnan(score)) return ".";
        char buf[32];
        snprintf(buf, sizeof(buf), "%.4f", score);
        return buf;
    }
};

/**
 * PhyloP Conservation Score Source
 *
 * PhyloP measures evolutionary conservation using phylogenetic p-values.
 * Positive scores indicate conservation (slower evolution).
 * Negative scores indicate acceleration (faster evolution).
 */
class PhyloPSource : public BigWigScoreSource {
public:
    explicit PhyloPSource(const std::string& path)
        : BigWigScoreSource(path, "phylop", "score",
                           "PhyloP evolutionary conservation score") {}

    void annotate(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript* transcript,
        std::unordered_map<std::string, std::string>& annotations
    ) override {
        // Call parent implementation
        BigWigScoreSource::annotate(chrom, pos, ref, alt, transcript, annotations);

        // Add interpretation
        auto it = annotations.find("phylop:score");
        if (it != annotations.end() && it->second != ".") {
            try {
                double score = std::stod(it->second);
                std::string interp = interpret_phylop(score);
                if (!interp.empty()) {
                    annotations["phylop:interpretation"] = interp;
                }
            } catch (...) {}
        }
    }

    std::vector<std::string> get_fields() const override {
        return { "phylop:score", "phylop:interpretation" };
    }

    /**
     * Interpret PhyloP score
     * @return "highly_conserved" (>= 2), "conserved" (>= 0.5), "neutral", or "accelerated" (< -1)
     */
    static std::string interpret_phylop(double score) {
        if (score >= 2.0) return "highly_conserved";
        if (score >= 0.5) return "conserved";
        if (score < -1.0) return "accelerated";
        return "neutral";
    }
};

/**
 * PhastCons Conservation Score Source
 *
 * PhastCons measures probability of negative selection.
 * Scores range from 0 to 1, where 1 indicates highest conservation.
 */
class PhastConsSource : public BigWigScoreSource {
public:
    explicit PhastConsSource(const std::string& path)
        : BigWigScoreSource(path, "phastcons", "score",
                           "PhastCons conservation probability") {}

    void annotate(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript* transcript,
        std::unordered_map<std::string, std::string>& annotations
    ) override {
        BigWigScoreSource::annotate(chrom, pos, ref, alt, transcript, annotations);

        auto it = annotations.find("phastcons:score");
        if (it != annotations.end() && it->second != ".") {
            try {
                double score = std::stod(it->second);
                std::string interp = interpret_phastcons(score);
                if (!interp.empty()) {
                    annotations["phastcons:interpretation"] = interp;
                }
            } catch (...) {}
        }
    }

    std::vector<std::string> get_fields() const override {
        return { "phastcons:score", "phastcons:interpretation" };
    }

    /**
     * Interpret PhastCons score
     * @return "highly_conserved" (>= 0.9), "conserved" (>= 0.5), "moderately_conserved" (>= 0.2), or "not_conserved"
     */
    static std::string interpret_phastcons(double score) {
        if (score >= 0.9) return "highly_conserved";
        if (score >= 0.5) return "conserved";
        if (score >= 0.2) return "moderately_conserved";
        return "not_conserved";
    }
};

/**
 * GERP++ Conservation Score Source
 *
 * GERP++ (Genomic Evolutionary Rate Profiling) measures conservation
 * as rejected substitutions - the difference between expected and observed substitutions.
 * Higher scores indicate more conservation.
 */
class GERPSource : public BigWigScoreSource {
public:
    explicit GERPSource(const std::string& path)
        : BigWigScoreSource(path, "gerp", "RS",
                           "GERP++ rejected substitutions score") {}

    void annotate(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript* transcript,
        std::unordered_map<std::string, std::string>& annotations
    ) override {
        BigWigScoreSource::annotate(chrom, pos, ref, alt, transcript, annotations);

        auto it = annotations.find("gerp:RS");
        if (it != annotations.end() && it->second != ".") {
            try {
                double score = std::stod(it->second);
                std::string interp = interpret_gerp(score);
                if (!interp.empty()) {
                    annotations["gerp:interpretation"] = interp;
                }
            } catch (...) {}
        }
    }

    std::vector<std::string> get_fields() const override {
        return { "gerp:RS", "gerp:interpretation" };
    }

    /**
     * Interpret GERP++ RS score
     * @return "highly_conserved" (>= 4), "conserved" (>= 2), "moderately_conserved" (>= 0), or "not_conserved"
     */
    static std::string interpret_gerp(double score) {
        if (score >= 4.0) return "highly_conserved";
        if (score >= 2.0) return "conserved";
        if (score >= 0.0) return "moderately_conserved";
        return "not_conserved";
    }
};

/**
 * Factory functions
 */
std::shared_ptr<AnnotationSource> create_phylop_source(const std::string& path) {
    return std::make_shared<PhyloPSource>(path);
}

std::shared_ptr<AnnotationSource> create_phastcons_source(const std::string& path) {
    return std::make_shared<PhastConsSource>(path);
}

std::shared_ptr<AnnotationSource> create_gerp_source(const std::string& path) {
    return std::make_shared<GERPSource>(path);
}

} // namespace vep
