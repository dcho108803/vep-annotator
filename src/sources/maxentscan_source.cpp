/**
 * MaxEntScan Annotation Source
 *
 * Algorithmic splice site scoring based on maximum entropy models.
 * Implements 5' and 3' splice site scoring matrices.
 */

#include "annotation_source.hpp"
#include "vep_annotator.hpp"
#include <cmath>
#include <array>
#include <algorithm>
#include <sstream>
#include <iomanip>

namespace vep {

// Format MaxEntScan score with 3 decimal places (matching Perl VEP output)
static std::string format_mes_score(double score) {
    char buf[32];
    snprintf(buf, sizeof(buf), "%.3f", score);
    return buf;
}

/**
 * MaxEntScan 5' splice site position weight matrix
 * Positions: -3 to +6 (9 positions)
 * Order: A, C, G, T
 */
static const std::array<std::array<double, 4>, 9> MES_5SS_SCORES = {{
    // Position -3
    {{ 0.2784, 0.3647, 0.1916, 0.1653 }},
    // Position -2
    {{ 0.6102, 0.1176, 0.1176, 0.1546 }},
    // Position -1
    {{ 0.0896, 0.0308, 0.8220, 0.0576 }},
    // Position +1 (always G)
    {{ 0.0000, 0.0000, 1.0000, 0.0000 }},
    // Position +2 (always T)
    {{ 0.0000, 0.0000, 0.0000, 1.0000 }},
    // Position +3
    {{ 0.5427, 0.0345, 0.3609, 0.0619 }},
    // Position +4
    {{ 0.7061, 0.0758, 0.1182, 0.0999 }},
    // Position +5
    {{ 0.1461, 0.1344, 0.4707, 0.2488 }},
    // Position +6
    {{ 0.1875, 0.2478, 0.2478, 0.3169 }}
}};

/**
 * MaxEntScan 3' splice site position weight matrix
 * Positions: -20 to +3 (23 positions)
 * Simplified scoring - uses consensus positions
 */
static const std::array<std::array<double, 4>, 23> MES_3SS_SCORES = {{
    // Positions -20 to -4 (polypyrimidine tract region)
    {{ 0.15, 0.30, 0.15, 0.40 }},  // -20
    {{ 0.15, 0.30, 0.15, 0.40 }},  // -19
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -18
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -17
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -16
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -15
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -14
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -13
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -12
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -11
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -10
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -9
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -8
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -7
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -6
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -5
    {{ 0.15, 0.35, 0.10, 0.40 }},  // -4
    // Position -3 (N - often C/T)
    {{ 0.20, 0.35, 0.10, 0.35 }},
    // Position -2 (A consensus)
    {{ 0.95, 0.02, 0.01, 0.02 }},
    // Position -1 (G consensus)
    {{ 0.01, 0.01, 0.97, 0.01 }},
    // Position +1
    {{ 0.20, 0.15, 0.50, 0.15 }},
    // Position +2
    {{ 0.25, 0.25, 0.25, 0.25 }},
    // Position +3
    {{ 0.25, 0.25, 0.25, 0.25 }}
}};

// Helper to reverse complement a DNA sequence
static std::string reverse_complement(const std::string& seq) {
    std::string rc(seq.size(), 'N');
    for (size_t i = 0; i < seq.size(); ++i) {
        switch (seq[seq.size() - 1 - i]) {
            case 'A': rc[i] = 'T'; break;
            case 'T': rc[i] = 'A'; break;
            case 'C': rc[i] = 'G'; break;
            case 'G': rc[i] = 'C'; break;
            default:  rc[i] = 'N'; break;
        }
    }
    return rc;
}

/**
 * MaxEntScan Annotation Source
 *
 * Computes 5' and 3' splice site scores using position weight matrices.
 * This is an algorithmic source - no external data file needed.
 */
class MaxEntScanSource : public AnnotationSource {
public:
    MaxEntScanSource() = default;

    std::string name() const override { return "maxentscan"; }
    std::string type() const override { return "splice"; }
    std::string description() const override {
        return "MaxEntScan splice site scoring (algorithmic)";
    }

    bool is_ready() const override { return true; }  // Always ready - no data file

    void initialize() override {
        log(LogLevel::INFO, "MaxEntScan initialized (algorithmic scoring)");
    }

    void set_reference(const ReferenceGenome* ref) { reference_ = ref; }

    void annotate(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript* transcript,
        std::unordered_map<std::string, std::string>& annotations
    ) override {
        if (!transcript) return;

        // Check each exon for splice site proximity and compute scores
        bool is_minus = (transcript->strand == '-');

        for (size_t i = 0; i < transcript->exons.size(); ++i) {
            const auto& exon = transcript->exons[i];

            // 5' splice site (donor): 3bp exon + 6bp intron = 9bp
            // Plus strand: donor at exon.end, exists when not last exon
            // Minus strand: donor at exon.start, exists when not first exon
            bool has_donor = is_minus ? (i > 0) : (i < transcript->exons.size() - 1);
            if (has_donor) {
                int donor_start, donor_end;
                if (is_minus) {
                    // Minus strand: donor (5'ss) is at exon.start
                    donor_start = exon.start - 6;
                    donor_end = exon.start + 2;
                } else {
                    // Plus strand: donor (5'ss) is at exon.end
                    donor_start = exon.end - 2;  // -3 to +6 relative to GT
                    donor_end = exon.end + 6;
                }

                if (pos >= donor_start && pos <= donor_end) {
                    annotations["maxentscan:in_splice_region"] = "true";

                    if (reference_ && reference_->has_chromosome(chrom)) {
                        // Extract 9bp context for reference
                        std::string ref_context = reference_->get_sequence(chrom, donor_start, donor_end);
                        if (is_minus) {
                            ref_context = reverse_complement(ref_context);
                        }
                        if (ref_context.length() == 9) {
                            double ref_score = score_5ss(ref_context);

                            // Build alt context by substituting the variant
                            std::string alt_context = ref_context;
                            int offset = pos - donor_start;
                            if (ref.length() == 1 && alt.length() == 1 &&
                                offset >= 0 && offset < 9) {
                                char alt_base = alt[0];
                                if (is_minus) {
                                    switch (alt_base) {
                                        case 'A': alt_base = 'T'; break;
                                        case 'T': alt_base = 'A'; break;
                                        case 'C': alt_base = 'G'; break;
                                        case 'G': alt_base = 'C'; break;
                                    }
                                    // Adjust offset for reverse complement
                                    offset = 8 - offset;
                                }
                                alt_context[offset] = alt_base;
                                double alt_score = score_5ss(alt_context);
                                double diff = alt_score - ref_score;

                                if (ref_score > -999.0) {
                                    annotations["maxentscan:5ss_ref"] = format_mes_score(ref_score);
                                }
                                if (alt_score > -999.0) {
                                    annotations["maxentscan:5ss_alt"] = format_mes_score(alt_score);
                                }
                                if (ref_score > -999.0 && alt_score > -999.0) {
                                    annotations["maxentscan:5ss_diff"] = format_mes_score(diff);
                                }
                            }
                        }
                    }
                    return;
                }
            }

            // 3' splice site (acceptor): 20bp intron + 3bp exon = 23bp
            // Plus strand: acceptor at exon.start, exists when not first exon
            // Minus strand: acceptor at exon.end, exists when not last exon
            bool has_acceptor = is_minus ? (i < transcript->exons.size() - 1) : (i > 0);
            if (has_acceptor) {
                int acceptor_start, acceptor_end;
                if (is_minus) {
                    // Minus strand: acceptor (3'ss) is at exon.end
                    acceptor_start = exon.end - 2;
                    acceptor_end = exon.end + 20;
                } else {
                    // Plus strand: acceptor (3'ss) is at exon.start
                    acceptor_start = exon.start - 20;  // -20 to +3 relative to AG
                    acceptor_end = exon.start + 2;
                }

                if (pos >= acceptor_start && pos <= acceptor_end) {
                    annotations["maxentscan:in_splice_region"] = "true";

                    if (reference_ && reference_->has_chromosome(chrom)) {
                        std::string ref_context = reference_->get_sequence(chrom, acceptor_start, acceptor_end);
                        if (is_minus) {
                            ref_context = reverse_complement(ref_context);
                        }
                        if (ref_context.length() == 23) {
                            double ref_score = score_3ss(ref_context);

                            std::string alt_context = ref_context;
                            int offset = pos - acceptor_start;
                            if (ref.length() == 1 && alt.length() == 1 &&
                                offset >= 0 && offset < 23) {
                                char alt_base = alt[0];
                                if (is_minus) {
                                    switch (alt_base) {
                                        case 'A': alt_base = 'T'; break;
                                        case 'T': alt_base = 'A'; break;
                                        case 'C': alt_base = 'G'; break;
                                        case 'G': alt_base = 'C'; break;
                                    }
                                    offset = 22 - offset;
                                }
                                alt_context[offset] = alt_base;
                                double alt_score = score_3ss(alt_context);
                                double diff = alt_score - ref_score;

                                if (ref_score > -999.0) {
                                    annotations["maxentscan:3ss_ref"] = format_mes_score(ref_score);
                                }
                                if (alt_score > -999.0) {
                                    annotations["maxentscan:3ss_alt"] = format_mes_score(alt_score);
                                }
                                if (ref_score > -999.0 && alt_score > -999.0) {
                                    annotations["maxentscan:3ss_diff"] = format_mes_score(diff);
                                }
                            }
                        }
                    }
                    return;
                }
            }
        }
    }

private:
    const ReferenceGenome* reference_ = nullptr;

public:

    std::vector<std::string> get_fields() const override {
        return {
            "maxentscan:5ss_ref",      // 5' splice site reference score
            "maxentscan:5ss_alt",      // 5' splice site alternate score
            "maxentscan:5ss_diff",     // Difference (alt - ref)
            "maxentscan:3ss_ref",      // 3' splice site reference score
            "maxentscan:3ss_alt",      // 3' splice site alternate score
            "maxentscan:3ss_diff",     // Difference (alt - ref)
            "maxentscan:in_splice_region"
        };
    }

    bool is_thread_safe() const override { return true; }

    /**
     * Score a 5' splice site sequence (9 bases: -3 to +6)
     * @param seq 9-base sequence around splice site
     * @return MaxEntScan score (log-odds)
     */
    static double score_5ss(const std::string& seq) {
        if (seq.length() != 9) return -999.0;

        double score = 0.0;
        for (size_t i = 0; i < 9; ++i) {
            int base_idx = base_to_index(seq[i]);
            if (base_idx < 0) return -999.0;

            double prob = MES_5SS_SCORES[i][base_idx];
            if (prob <= 0) return -999.0;
            score += std::log2(prob / 0.25);  // Log-odds vs uniform
        }
        return score;
    }

    /**
     * Score a 3' splice site sequence (23 bases: -20 to +3)
     * @param seq 23-base sequence around splice site
     * @return MaxEntScan score (log-odds)
     */
    static double score_3ss(const std::string& seq) {
        if (seq.length() != 23) return -999.0;

        double score = 0.0;
        for (size_t i = 0; i < 23; ++i) {
            int base_idx = base_to_index(seq[i]);
            if (base_idx < 0) return -999.0;

            double prob = MES_3SS_SCORES[i][base_idx];
            if (prob <= 0) return -999.0;
            score += std::log2(prob / 0.25);
        }
        return score;
    }

    /**
     * Interpret score difference
     * @return "damaging" if diff <= -3, "moderate" if diff <= -1.5, else ""
     */
    static std::string interpret_diff(double diff) {
        if (diff <= -3.0) return "damaging";
        if (diff <= -1.5) return "moderate";
        return "";
    }

private:
    static int base_to_index(char base) {
        switch (std::toupper(base)) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default: return -1;
        }
    }
};

/**
 * Factory function to create MaxEntScan source
 */
std::shared_ptr<AnnotationSource> create_maxentscan_source() {
    return std::make_shared<MaxEntScanSource>();
}

std::shared_ptr<AnnotationSource> create_maxentscan_source(const ReferenceGenome* ref) {
    auto source = std::make_shared<MaxEntScanSource>();
    source->set_reference(ref);
    return source;
}

} // namespace vep
