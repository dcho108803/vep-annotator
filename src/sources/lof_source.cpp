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
#include <cstdio>
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
        std::unordered_map<std::string, std::string>& annotations
    ) override {
        if (!transcript) return;
        if (!transcript->is_coding()) return;

        (void)chrom;
        (void)ref;
        (void)alt;

        // Only classify LoF for truncating consequence types.
        // Check if the variant is at a position that would cause a LoF consequence:
        // - In splice donor (2bp after exon end) or acceptor (2bp before exon start)
        // - In CDS region (potential stop_gained or frameshift)
        bool is_splice_site = false;
        bool is_in_cds = false;
        for (size_t i = 0; i < transcript->exons.size(); ++i) {
            const auto& exon = transcript->exons[i];
            // Splice acceptor: 2bp before exon start
            if (i > 0 && pos >= exon.start - 2 && pos <= exon.start - 1) {
                is_splice_site = true;
            }
            // Splice donor: 2bp after exon end
            if (i < transcript->exons.size() - 1 && pos >= exon.end + 1 && pos <= exon.end + 2) {
                is_splice_site = true;
            }
        }
        for (const auto& cds : transcript->cds_regions) {
            if (pos >= cds.start && pos <= cds.end) {
                is_in_cds = true;
                break;
            }
        }

        // Only proceed for splice-site or CDS variants (potential LoF)
        if (!is_splice_site && !is_in_cds) return;

        // Check for truncating location in last exon
        bool is_last_exon = false;
        bool is_single_exon = (transcript->exons.size() == 1);

        for (const auto& exon : transcript->exons) {
            if (pos >= exon.start && pos <= exon.end) {
                // Exons are sorted by genomic position (ascending).
                // For plus strand: last transcript-order exon = exons.back() (highest genomic pos)
                // For minus strand: last transcript-order exon = exons.front() (lowest genomic pos)
                const auto& last_exon = (transcript->strand == '-') ?
                    transcript->exons.front() : transcript->exons.back();
                if (&exon == &last_exon) {
                    is_last_exon = true;
                }
            }
        }

        // LOFTEE flags and classification
        std::vector<std::string> flags;
        bool is_hc = true;  // Start as HC, downgrade based on flags

        // Flag: Last exon (stop_gained in last exon may escape NMD)
        if (is_last_exon && !is_single_exon && is_in_cds) {
            flags.push_back("END_TRUNC");
            is_hc = false;  // NMD escape - downgrade to LC
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

        // Flag: Non-protein-coding
        if (transcript->biotype != "protein_coding") {
            flags.push_back("NON_CODING");
            is_hc = false;
        }

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
        char conf_buf[32]; std::snprintf(conf_buf, sizeof(conf_buf), "%.4f", std::max(0.1, confidence));
        annotations["loftee:confidence"] = conf_buf;
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
        std::unordered_map<std::string, std::string>& annotations
    ) override {
        if (!transcript) return;
        if (!transcript->is_coding()) return;

        (void)chrom;
        (void)ref;
        (void)alt;

        // NMD prediction based on 50bp rule (in CDS coordinates)
        // A premature termination codon (PTC) triggers NMD if:
        // 1. It is >50-55bp upstream of the last exon-exon junction IN CDS COORDINATES
        // 2. The gene has more than one exon

        // Single exon genes don't undergo NMD
        if (transcript->exons.size() <= 1) {
            annotations["nmd:susceptible"] = "false";
            annotations["nmd:reason"] = "single_exon_gene";
            return;
        }

        // Calculate CDS position of the variant
        int variant_cds_pos = calculate_cds_pos(pos, *transcript);
        if (variant_cds_pos <= 0) {
            annotations["nmd:susceptible"] = "false";
            annotations["nmd:reason"] = "not_in_cds";
            return;
        }

        // Calculate CDS position of the last exon-exon junction
        // The last junction is at the start of the last exon (+ strand)
        // or end of the first exon (- strand).
        // In CDS terms, this is the CDS position where the last exon begins.
        int junction_genomic = 0;
        if (transcript->strand == '+') {
            junction_genomic = transcript->exons.back().start;
        } else {
            junction_genomic = transcript->exons.front().end;
        }

        int junction_cds_pos = calculate_cds_pos(junction_genomic, *transcript);

        // If the last exon-exon junction is not in CDS (last exon is entirely UTR),
        // then all CDS is upstream of the junction - variant is effectively in last coding exon
        if (junction_cds_pos <= 0) {
            annotations["nmd:susceptible"] = "false";
            annotations["nmd:reason"] = "last_junction_in_utr";
            annotations["nmd:distance_to_junction"] = "0";
            return;
        }

        // Distance in CDS coordinates (how far PTC is upstream of last junction)
        int cds_distance = junction_cds_pos - variant_cds_pos;

        // Apply 50bp rule in CDS coordinates
        bool susceptible = (cds_distance > 50);

        annotations["nmd:susceptible"] = susceptible ? "true" : "false";
        annotations["nmd:distance_to_junction"] = std::to_string(cds_distance);

        if (susceptible) {
            annotations["nmd:reason"] = "ptc_upstream_of_last_junction";
        } else if (cds_distance <= 0) {
            annotations["nmd:reason"] = "in_last_exon";
        } else {
            annotations["nmd:reason"] = "within_50bp_of_junction";
        }
    }

    // Calculate CDS position from genomic position using transcript CDS regions
    static int calculate_cds_pos(int genomic_pos, const Transcript& transcript) {
        int cds_pos = 0;
        if (transcript.strand == '+') {
            for (const auto& cds : transcript.cds_regions) {
                if (genomic_pos >= cds.start && genomic_pos <= cds.end) {
                    return cds_pos + (genomic_pos - cds.start + 1);
                }
                cds_pos += cds.end - cds.start + 1;
            }
        } else {
            // For minus strand, iterate CDS regions in reverse
            std::vector<CDS> sorted_cds(transcript.cds_regions.rbegin(),
                                        transcript.cds_regions.rend());
            for (const auto& cds : sorted_cds) {
                if (genomic_pos >= cds.start && genomic_pos <= cds.end) {
                    return cds_pos + (cds.end - genomic_pos + 1);
                }
                cds_pos += cds.end - cds.start + 1;
            }
        }
        return 0; // Not in CDS
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
        std::lock_guard<std::recursive_mutex> lock(mutex_);
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
        std::unordered_map<std::string, std::string>& annotations
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
            char score_buf[32]; std::snprintf(score_buf, sizeof(score_buf), "%.4f", it->second);
            annotations["loftool:score"] = score_buf;

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
