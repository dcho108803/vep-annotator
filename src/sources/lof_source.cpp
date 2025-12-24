/**
 * Loss-of-Function (LoF) Annotation Sources
 *
 * Provides LoF classifications and predictions:
 * - LOFTEE: Loss-of-function Transcript Effect Estimator
 * - NMD: Nonsense-Mediated Decay prediction
 * - LoFtool: Gene-level constraint scores
 */

#include "annotation_source.hpp"
#include "file_parsers.hpp"
#include "vep_annotator.hpp"
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>

namespace vep {

/**
 * LOFTEE-style LoF Annotation Source
 *
 * Classifies loss-of-function variants as High Confidence (HC) or Low Confidence (LC)
 * based on various features like position within transcript, NMD susceptibility, etc.
 */
class LOFTEESource : public AnnotationSource {
public:
    LOFTEESource() = default;

    std::string name() const override { return "loftee"; }
    std::string type() const override { return "lof"; }
    std::string description() const override {
        return "LOFTEE-style LoF classification (HC/LC)";
    }

    bool is_ready() const override { return true; }  // Algorithmic - always ready

    void initialize() override {
        log(LogLevel::INFO, "LOFTEE LoF classifier initialized");
    }

    void annotate(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript* transcript,
        std::map<std::string, std::string>& annotations
    ) override {
        if (!transcript) return;

        (void)chrom;
        (void)ref;
        (void)alt;

        // Check if this is a potential LoF variant based on consequences
        // (This would typically be done based on the consequence determination
        // from the main annotator, but we can check position-based features here)

        // Check for truncating location in last exon
        bool is_last_exon = false;
        bool is_single_exon = (transcript->exons.size() == 1);
        int cds_length = 0;
        int pos_in_cds = 0;

        // Calculate position within CDS
        for (const auto& exon : transcript->exons) {
            if (pos >= exon.start && pos <= exon.end) {
                // Found the exon
                if (&exon == &transcript->exons.back()) {
                    is_last_exon = true;
                }
            }
            // Approximate CDS length (would need actual CDS coordinates for accuracy)
            cds_length += exon.end - exon.start + 1;
        }

        // Determine relative position in transcript
        double relative_pos = 0;
        if (transcript->end > transcript->start) {
            relative_pos = static_cast<double>(pos - transcript->start) /
                          (transcript->end - transcript->start);
        }

        // LOFTEE flags and classification
        std::vector<std::string> flags;
        bool is_hc = true;  // Start as HC, downgrade based on flags

        // Flag: Last exon
        if (is_last_exon && !is_single_exon) {
            flags.push_back("END_TRUNC");
            // Variants in last 50bp of last exon may escape NMD
        }

        // Flag: Single exon gene
        if (is_single_exon) {
            flags.push_back("SINGLE_EXON");
        }

        // Flag: Incomplete CDS
        if (transcript->cds_start <= 0 || transcript->cds_end <= 0) {
            flags.push_back("INCOMPLETE_CDS");
            is_hc = false;
        }

        // Flag: Non-canonical splice (would need actual consequence)
        // This is typically done based on splice consequence type

        // Flag: NAGNAG splice site (would need sequence context)

        // Set annotations
        if (is_hc) {
            annotations["loftee:classification"] = "HC";
        } else {
            annotations["loftee:classification"] = "LC";
        }

        if (!flags.empty()) {
            std::ostringstream oss;
            for (size_t i = 0; i < flags.size(); ++i) {
                if (i > 0) oss << ",";
                oss << flags[i];
            }
            annotations["loftee:flags"] = oss.str();
        }

        // Confidence score (simplified)
        double confidence = is_hc ? 0.9 : 0.5;
        if (!flags.empty()) {
            confidence -= flags.size() * 0.1;
        }
        annotations["loftee:confidence"] = std::to_string(std::max(0.1, confidence));
    }

    std::vector<std::string> get_fields() const override {
        return {
            "loftee:classification",
            "loftee:flags",
            "loftee:confidence"
        };
    }

    bool is_thread_safe() const override { return true; }
};

/**
 * NMD Prediction Source
 *
 * Predicts susceptibility to Nonsense-Mediated Decay based on
 * the 50-55bp rule (PTC must be >50-55bp upstream of last exon-exon junction).
 */
class NMDSource : public AnnotationSource {
public:
    NMDSource() = default;

    std::string name() const override { return "nmd"; }
    std::string type() const override { return "lof"; }
    std::string description() const override {
        return "Nonsense-Mediated Decay prediction";
    }

    bool is_ready() const override { return true; }

    void initialize() override {
        log(LogLevel::INFO, "NMD predictor initialized");
    }

    void annotate(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript* transcript,
        std::map<std::string, std::string>& annotations
    ) override {
        if (!transcript) return;

        (void)chrom;
        (void)ref;
        (void)alt;

        // NMD prediction based on 50bp rule
        // A premature termination codon (PTC) triggers NMD if:
        // 1. It is >50-55bp upstream of the last exon-exon junction
        // 2. The gene has more than one exon

        // Single exon genes don't undergo NMD
        if (transcript->exons.size() <= 1) {
            annotations["nmd:susceptible"] = "false";
            annotations["nmd:reason"] = "single_exon_gene";
            return;
        }

        // Find the last exon-exon junction
        int last_junction = 0;
        if (transcript->strand == '+') {
            // For + strand, junction is at start of last exon
            if (transcript->exons.size() >= 2) {
                last_junction = transcript->exons.back().start;
            }
        } else {
            // For - strand, junction is at end of first exon (on - strand)
            if (transcript->exons.size() >= 2) {
                last_junction = transcript->exons.front().end;
            }
        }

        // Calculate distance from variant to last junction
        int distance = 0;
        if (transcript->strand == '+') {
            distance = last_junction - pos;
        } else {
            distance = pos - last_junction;
        }

        // Apply 50bp rule
        bool susceptible = (distance > 50);

        annotations["nmd:susceptible"] = susceptible ? "true" : "false";
        annotations["nmd:distance_to_junction"] = std::to_string(distance);

        if (susceptible) {
            annotations["nmd:reason"] = "ptc_upstream_of_last_junction";
        } else if (distance <= 0) {
            annotations["nmd:reason"] = "in_last_exon";
        } else {
            annotations["nmd:reason"] = "within_50bp_of_junction";
        }
    }

    std::vector<std::string> get_fields() const override {
        return {
            "nmd:susceptible",
            "nmd:distance_to_junction",
            "nmd:reason"
        };
    }

    bool is_thread_safe() const override { return true; }
};

/**
 * LoFtool Gene Tolerance Source
 *
 * Provides gene-level constraint/tolerance scores from LoFtool.
 * Lower scores indicate less tolerance to LoF variants (more constrained).
 */
class LoFtoolSource : public AnnotationSource {
public:
    explicit LoFtoolSource(const std::string& path)
        : path_(path) {}

    std::string name() const override { return "loftool"; }
    std::string type() const override { return "lof"; }
    std::string description() const override {
        return "LoFtool gene constraint scores";
    }

    bool is_ready() const override { return loaded_; }

    void initialize() override {
        std::lock_guard<std::mutex> lock(mutex_);
        if (loaded_) return;

        log(LogLevel::INFO, "Loading LoFtool scores from: " + path_);

        std::ifstream file(path_);
        if (!file.is_open()) {
            log(LogLevel::ERROR, "Failed to open LoFtool file: " + path_);
            return;
        }

        std::string line;
        size_t count = 0;

        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;

            auto fields = split_line(line, '\t');
            if (fields.size() < 2) continue;

            std::string gene = fields[0];
            try {
                double score = std::stod(fields[1]);
                scores_[gene] = score;
                count++;
            } catch (...) {
                continue;
            }
        }

        loaded_ = true;
        log(LogLevel::INFO, "LoFtool loaded " + std::to_string(count) + " gene scores");
    }

    void annotate(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript* transcript,
        std::map<std::string, std::string>& annotations
    ) override {
        ensure_initialized();
        if (!transcript) return;

        (void)chrom;
        (void)pos;
        (void)ref;
        (void)alt;

        // Look up score by gene name
        auto it = scores_.find(transcript->gene_name);
        if (it == scores_.end()) {
            // Try gene ID
            it = scores_.find(transcript->gene_id);
        }

        if (it != scores_.end()) {
            annotations["loftool:score"] = std::to_string(it->second);

            // Interpret score (lower = less tolerant = more constrained)
            std::string interpretation = interpret_score(it->second);
            if (!interpretation.empty()) {
                annotations["loftool:interpretation"] = interpretation;
            }
        }
    }

    std::vector<std::string> get_fields() const override {
        return {
            "loftool:score",
            "loftool:interpretation"
        };
    }

    std::string get_data_path() const override { return path_; }

    /**
     * Interpret LoFtool score
     * Lower score = more intolerant to LoF
     */
    static std::string interpret_score(double score) {
        // LoFtool scores are percentiles (0-1)
        if (score <= 0.1) return "highly_intolerant";
        if (score <= 0.35) return "intolerant";
        if (score <= 0.65) return "neutral";
        return "tolerant";
    }

private:
    std::string path_;
    bool loaded_ = false;
    std::map<std::string, double> scores_;
};

/**
 * Factory functions
 */
std::shared_ptr<AnnotationSource> create_loftee_source() {
    return std::make_shared<LOFTEESource>();
}

std::shared_ptr<AnnotationSource> create_nmd_source() {
    return std::make_shared<NMDSource>();
}

std::shared_ptr<AnnotationSource> create_loftool_source(const std::string& path) {
    return std::make_shared<LoFtoolSource>(path);
}

} // namespace vep
