/**
 * Transcript Filter - Transcript Selection and Filtering Options
 *
 * Implements --pick, --pick_allele, --most_severe, --per_gene, --canonical,
 * --mane, --biotype filtering, and other transcript selection options.
 */

#ifndef TRANSCRIPT_FILTER_HPP
#define TRANSCRIPT_FILTER_HPP

#include "vep_annotator.hpp"
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <functional>
#include <sstream>
#include <cctype>

namespace vep {

/**
 * Transcript pick criteria for ranking
 */
enum class PickCriteria {
    CANONICAL,      // Is canonical transcript
    MANE_SELECT,    // Is MANE Select transcript
    MANE_PLUS,      // Is MANE Plus Clinical transcript
    APPRIS,         // APPRIS principal isoform
    TSL,            // Transcript Support Level (lower is better)
    BIOTYPE,        // protein_coding preferred
    CCDS,           // Has CCDS annotation
    RANK,           // Consequence severity rank
    LENGTH          // Transcript length (longer preferred)
};

/**
 * Configuration for transcript filtering
 */
struct TranscriptFilterConfig {
    // Selection mode (only one should be true)
    bool pick = false;              // One annotation per variant
    bool pick_allele = false;       // One annotation per allele
    bool pick_allele_gene = false;  // One annotation per allele per gene
    bool per_gene = false;          // One annotation per gene
    bool most_severe = false;       // Only most severe consequence
    bool flag_pick = false;         // Flag picked without filtering
    bool flag_pick_allele = false;  // Flag picked per allele without filtering
    bool flag_pick_allele_gene = false; // Flag picked per allele+gene without filtering

    // Filtering options
    bool canonical_only = false;    // Only canonical transcripts
    bool mane_only = false;         // Only MANE Select transcripts
    bool coding_only = false;       // Only protein_coding transcripts
    bool gencode_basic = false;     // Only GENCODE basic transcripts
    bool all_refseq = false;        // Include all RefSeq transcripts
    std::set<std::string> biotypes; // Allowed biotypes (empty = all)

    // Consequence filtering
    std::set<std::string> include_consequences;  // Only these consequences
    std::set<std::string> exclude_consequences;  // Exclude these consequences
    std::set<Impact> include_impacts;            // Only these impact levels

    // Pick order
    std::vector<PickCriteria> pick_order;

    // Output options
    bool show_canonical = false;    // Add CANONICAL column
    bool show_mane = false;         // Add MANE column
    bool show_tsl = false;          // Add TSL column
    bool show_appris = false;       // Add APPRIS column

    // Frequency filtering
    bool check_frequency = false;   // Enable frequency filtering
    std::string freq_pop;           // Population to check
    double freq_threshold = 0.01;   // Maximum frequency
    bool freq_gt = false;           // true = greater than, false = less than

    // No intergenic
    bool no_intergenic = false;     // Skip intergenic variants

    TranscriptFilterConfig() {
        // Default pick order
        pick_order.push_back(PickCriteria::CANONICAL);
        pick_order.push_back(PickCriteria::MANE_SELECT);
        pick_order.push_back(PickCriteria::APPRIS);
        pick_order.push_back(PickCriteria::TSL);
        pick_order.push_back(PickCriteria::BIOTYPE);
        pick_order.push_back(PickCriteria::CCDS);
        pick_order.push_back(PickCriteria::RANK);
        pick_order.push_back(PickCriteria::LENGTH);
        freq_pop = "gnomAD_AF";
    }
};

/**
 * Get consequence rank (lower = more severe)
 */
inline int get_consequence_rank(ConsequenceType csq) {
    if (csq == ConsequenceType::TRANSCRIPT_ABLATION) return 1;
    if (csq == ConsequenceType::SPLICE_ACCEPTOR_VARIANT) return 2;
    if (csq == ConsequenceType::SPLICE_DONOR_VARIANT) return 3;
    if (csq == ConsequenceType::STOP_GAINED) return 4;
    if (csq == ConsequenceType::FRAMESHIFT_VARIANT) return 5;
    if (csq == ConsequenceType::STOP_LOST) return 6;
    if (csq == ConsequenceType::START_LOST) return 7;
    if (csq == ConsequenceType::INFRAME_INSERTION) return 8;
    if (csq == ConsequenceType::INFRAME_DELETION) return 9;
    if (csq == ConsequenceType::MISSENSE_VARIANT) return 10;
    if (csq == ConsequenceType::PROTEIN_ALTERING_VARIANT) return 11;
    if (csq == ConsequenceType::SPLICE_REGION_VARIANT) return 12;
    if (csq == ConsequenceType::INCOMPLETE_TERMINAL_CODON_VARIANT) return 13;
    if (csq == ConsequenceType::START_RETAINED_VARIANT) return 14;
    if (csq == ConsequenceType::STOP_RETAINED_VARIANT) return 15;
    if (csq == ConsequenceType::SYNONYMOUS_VARIANT) return 16;
    if (csq == ConsequenceType::CODING_SEQUENCE_VARIANT) return 17;
    if (csq == ConsequenceType::MATURE_MIRNA_VARIANT) return 18;
    if (csq == ConsequenceType::FIVE_PRIME_UTR_VARIANT) return 19;
    if (csq == ConsequenceType::THREE_PRIME_UTR_VARIANT) return 20;
    if (csq == ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT) return 21;
    if (csq == ConsequenceType::INTRON_VARIANT) return 22;
    if (csq == ConsequenceType::NMD_TRANSCRIPT_VARIANT) return 23;
    if (csq == ConsequenceType::NON_CODING_TRANSCRIPT_VARIANT) return 24;
    if (csq == ConsequenceType::UPSTREAM_GENE_VARIANT) return 25;
    if (csq == ConsequenceType::DOWNSTREAM_GENE_VARIANT) return 26;
    if (csq == ConsequenceType::INTERGENIC_VARIANT) return 27;
    return 50;
}

/**
 * Extended annotation with additional metadata
 */
struct AnnotationWithMeta {
    VariantAnnotation annotation;
    bool is_picked = false;
    bool is_mane_select = false;
    bool is_mane_plus = false;
    int tsl = 0;                    // Transcript Support Level (1-5, lower is better)
    std::string appris;             // APPRIS annotation
    bool has_ccds = false;
    int transcript_length = 0;

    // Consequence rank (lower = more severe)
    int get_rank() const {
        if (annotation.consequences.empty()) return 100;

        int min_rank = 100;
        for (const auto& csq : annotation.consequences) {
            int rank = get_consequence_rank(csq);
            if (rank < min_rank) min_rank = rank;
        }
        return min_rank;
    }
};

/**
 * Compare APPRIS annotations
 */
inline int compare_appris(const std::string& a, const std::string& b) {
    auto score = [](const std::string& s) -> int {
        if (s.empty()) return 100;
        if (s.find("principal1") != std::string::npos) return 1;
        if (s.find("principal2") != std::string::npos) return 2;
        if (s.find("principal3") != std::string::npos) return 3;
        if (s.find("principal4") != std::string::npos) return 4;
        if (s.find("principal5") != std::string::npos) return 5;
        if (s.find("alternative1") != std::string::npos) return 6;
        if (s.find("alternative2") != std::string::npos) return 7;
        return 50;
    };
    return score(a) - score(b);
}

/**
 * Transcript filter class
 */
class TranscriptFilter {
public:
    explicit TranscriptFilter(const TranscriptFilterConfig& config = TranscriptFilterConfig())
        : config_(config) {}

    /**
     * Filter and rank annotations based on configuration
     * @param annotations Input annotations
     * @return Filtered and potentially picked annotations
     */
    std::vector<VariantAnnotation> filter(
        const std::vector<VariantAnnotation>& annotations) const {

        if (annotations.empty()) return annotations;

        // Convert to extended annotations
        std::vector<AnnotationWithMeta> extended;
        for (const auto& ann : annotations) {
            extended.push_back(extend_annotation(ann));
        }

        // Apply basic filters
        std::vector<AnnotationWithMeta> filtered = apply_basic_filters(extended);

        // Apply frequency filter if enabled
        if (config_.check_frequency) {
            filtered = apply_frequency_filter(filtered);
        }

        // Apply pick/selection logic
        if (config_.pick) {
            filtered = pick_one(filtered);
        } else if (config_.pick_allele) {
            filtered = pick_per_allele(filtered);
        } else if (config_.pick_allele_gene) {
            filtered = pick_per_allele_gene(filtered);
        } else if (config_.per_gene) {
            filtered = pick_per_gene(filtered);
        } else if (config_.most_severe) {
            filtered = pick_most_severe(filtered);
        } else if (config_.flag_pick) {
            flag_picked(filtered);
        } else if (config_.flag_pick_allele) {
            flag_picked_per_allele(filtered);
        } else if (config_.flag_pick_allele_gene) {
            flag_picked_per_allele_gene(filtered);
        }

        // Convert back to VariantAnnotation
        std::vector<VariantAnnotation> result;
        for (auto& ext : filtered) {
            if (config_.flag_pick && ext.is_picked) {
                ext.annotation.custom_annotations["PICK"] = "1";
            }
            result.push_back(ext.annotation);
        }

        return result;
    }

    /**
     * Check if annotation passes filters
     */
    bool passes_filter(const VariantAnnotation& ann) const {
        AnnotationWithMeta ext = extend_annotation(ann);
        return passes_basic_filter(ext);
    }

    const TranscriptFilterConfig& config() const { return config_; }
    void set_config(const TranscriptFilterConfig& config) { config_ = config; }

private:
    TranscriptFilterConfig config_;

    AnnotationWithMeta extend_annotation(const VariantAnnotation& ann) const {
        AnnotationWithMeta ext;
        ext.annotation = ann;
        ext.is_mane_select = false;
        ext.is_mane_plus = false;
        ext.tsl = 0;
        ext.has_ccds = false;

        // Extract from custom annotations if present
        auto mane_it = ann.custom_annotations.find("MANE_SELECT");
        if (mane_it != ann.custom_annotations.end() && !mane_it->second.empty()) {
            ext.is_mane_select = true;
        }

        auto tsl_it = ann.custom_annotations.find("TSL");
        if (tsl_it != ann.custom_annotations.end()) {
            try {
                ext.tsl = std::stoi(tsl_it->second);
            } catch (...) {}
        }

        auto appris_it = ann.custom_annotations.find("APPRIS");
        if (appris_it != ann.custom_annotations.end()) {
            ext.appris = appris_it->second;
        }

        return ext;
    }

    bool passes_basic_filter(const AnnotationWithMeta& ext) const {
        const VariantAnnotation& ann = ext.annotation;

        // Canonical only
        if (config_.canonical_only && !ann.is_canonical) {
            return false;
        }

        // MANE only
        if (config_.mane_only && !ext.is_mane_select && !ext.is_mane_plus) {
            return false;
        }

        // Coding only
        if (config_.coding_only && ann.biotype != "protein_coding") {
            return false;
        }

        // Biotype filter
        if (!config_.biotypes.empty() &&
            config_.biotypes.find(ann.biotype) == config_.biotypes.end()) {
            return false;
        }

        // No intergenic
        if (config_.no_intergenic) {
            bool is_intergenic = false;
            for (const auto& csq : ann.consequences) {
                if (csq == ConsequenceType::INTERGENIC_VARIANT) {
                    is_intergenic = true;
                    break;
                }
            }
            if (is_intergenic) return false;
        }

        // Consequence filter
        if (!config_.include_consequences.empty()) {
            bool has_included = false;
            for (const auto& csq : ann.consequences) {
                if (config_.include_consequences.find(consequence_to_string(csq))
                    != config_.include_consequences.end()) {
                    has_included = true;
                    break;
                }
            }
            if (!has_included) return false;
        }

        if (!config_.exclude_consequences.empty()) {
            for (const auto& csq : ann.consequences) {
                if (config_.exclude_consequences.find(consequence_to_string(csq))
                    != config_.exclude_consequences.end()) {
                    return false;
                }
            }
        }

        // Impact filter
        if (!config_.include_impacts.empty() &&
            config_.include_impacts.find(ann.impact) == config_.include_impacts.end()) {
            return false;
        }

        return true;
    }

    std::vector<AnnotationWithMeta> apply_basic_filters(
        const std::vector<AnnotationWithMeta>& annotations) const {

        std::vector<AnnotationWithMeta> result;
        for (const auto& ext : annotations) {
            if (passes_basic_filter(ext)) {
                result.push_back(ext);
            }
        }
        return result;
    }

    std::vector<AnnotationWithMeta> apply_frequency_filter(
        const std::vector<AnnotationWithMeta>& annotations) const {

        std::vector<AnnotationWithMeta> result;
        for (const auto& ext : annotations) {
            auto it = ext.annotation.custom_annotations.find(config_.freq_pop);
            if (it != ext.annotation.custom_annotations.end()) {
                try {
                    double freq = std::stod(it->second);
                    bool passes = config_.freq_gt ?
                        (freq > config_.freq_threshold) :
                        (freq < config_.freq_threshold);
                    if (!passes) continue;
                } catch (...) {
                    // Keep if can't parse frequency
                }
            }
            result.push_back(ext);
        }
        return result;
    }

    // Compare two annotations based on pick order
    bool is_better(const AnnotationWithMeta& a, const AnnotationWithMeta& b) const {
        for (const auto& criteria : config_.pick_order) {
            int cmp = compare_by_criteria(a, b, criteria);
            if (cmp != 0) return cmp < 0;
        }
        return false;
    }

    int compare_by_criteria(const AnnotationWithMeta& a, const AnnotationWithMeta& b,
                            PickCriteria criteria) const {
        if (criteria == PickCriteria::CANONICAL) {
            return (b.annotation.is_canonical ? 1 : 0) - (a.annotation.is_canonical ? 1 : 0);
        }
        if (criteria == PickCriteria::MANE_SELECT) {
            return (b.is_mane_select ? 1 : 0) - (a.is_mane_select ? 1 : 0);
        }
        if (criteria == PickCriteria::MANE_PLUS) {
            return (b.is_mane_plus ? 1 : 0) - (a.is_mane_plus ? 1 : 0);
        }
        if (criteria == PickCriteria::TSL) {
            // Lower TSL is better, 0 means unknown (worst)
            int tsl_a = a.tsl > 0 ? a.tsl : 10;
            int tsl_b = b.tsl > 0 ? b.tsl : 10;
            return tsl_a - tsl_b;
        }
        if (criteria == PickCriteria::BIOTYPE) {
            // protein_coding is best
            int score_a = (a.annotation.biotype == "protein_coding") ? 0 : 1;
            int score_b = (b.annotation.biotype == "protein_coding") ? 0 : 1;
            return score_a - score_b;
        }
        if (criteria == PickCriteria::CCDS) {
            return (b.has_ccds ? 1 : 0) - (a.has_ccds ? 1 : 0);
        }
        if (criteria == PickCriteria::RANK) {
            return a.get_rank() - b.get_rank();
        }
        if (criteria == PickCriteria::LENGTH) {
            return b.transcript_length - a.transcript_length;
        }
        if (criteria == PickCriteria::APPRIS) {
            return compare_appris(a.appris, b.appris);
        }
        return 0;
    }

    std::vector<AnnotationWithMeta> pick_one(
        std::vector<AnnotationWithMeta> annotations) const {

        if (annotations.empty()) return annotations;

        // Sort by pick criteria
        std::sort(annotations.begin(), annotations.end(),
            [this](const AnnotationWithMeta& a, const AnnotationWithMeta& b) {
                return is_better(a, b);
            });

        std::vector<AnnotationWithMeta> result;
        result.push_back(annotations[0]);
        return result;
    }

    std::vector<AnnotationWithMeta> pick_per_allele(
        std::vector<AnnotationWithMeta> annotations) const {

        if (annotations.empty()) return annotations;

        // Group by allele
        std::map<std::string, std::vector<AnnotationWithMeta>> by_allele;
        for (const auto& ann : annotations) {
            std::string key = ann.annotation.ref_allele + ">" + ann.annotation.alt_allele;
            by_allele[key].push_back(ann);
        }

        std::vector<AnnotationWithMeta> result;
        for (auto& pair : by_allele) {
            std::vector<AnnotationWithMeta> picked = pick_one(pair.second);
            if (!picked.empty()) {
                result.push_back(picked[0]);
            }
        }

        return result;
    }

    std::vector<AnnotationWithMeta> pick_per_gene(
        std::vector<AnnotationWithMeta> annotations) const {

        if (annotations.empty()) return annotations;

        // Group by gene
        std::map<std::string, std::vector<AnnotationWithMeta>> by_gene;
        for (const auto& ann : annotations) {
            by_gene[ann.annotation.gene_symbol].push_back(ann);
        }

        std::vector<AnnotationWithMeta> result;
        for (auto& pair : by_gene) {
            std::vector<AnnotationWithMeta> picked = pick_one(pair.second);
            if (!picked.empty()) {
                result.push_back(picked[0]);
            }
        }

        return result;
    }

    std::vector<AnnotationWithMeta> pick_most_severe(
        std::vector<AnnotationWithMeta> annotations) const {

        if (annotations.empty()) return annotations;

        // Find minimum rank
        int min_rank = 100;
        for (const auto& ann : annotations) {
            int rank = ann.get_rank();
            if (rank < min_rank) min_rank = rank;
        }

        // Return all with minimum rank
        std::vector<AnnotationWithMeta> result;
        for (const auto& ann : annotations) {
            if (ann.get_rank() == min_rank) {
                result.push_back(ann);
            }
        }

        // If still multiple, pick the best one
        if (result.size() > 1) {
            return pick_one(result);
        }

        return result;
    }

    void flag_picked(std::vector<AnnotationWithMeta>& annotations) const {
        if (annotations.empty()) return;

        // Sort to find the best one
        std::sort(annotations.begin(), annotations.end(),
            [this](const AnnotationWithMeta& a, const AnnotationWithMeta& b) {
                return is_better(a, b);
            });

        annotations[0].is_picked = true;
    }

    std::vector<AnnotationWithMeta> pick_per_allele_gene(
        std::vector<AnnotationWithMeta> annotations) const {

        if (annotations.empty()) return annotations;

        // Group by allele and gene
        std::map<std::string, std::vector<AnnotationWithMeta> > by_allele_gene;
        for (const auto& ann : annotations) {
            std::string key = ann.annotation.ref_allele + ">" + ann.annotation.alt_allele +
                              ":" + ann.annotation.gene_symbol;
            by_allele_gene[key].push_back(ann);
        }

        std::vector<AnnotationWithMeta> result;
        for (auto& pair : by_allele_gene) {
            std::vector<AnnotationWithMeta> picked = pick_one(pair.second);
            if (!picked.empty()) {
                result.push_back(picked[0]);
            }
        }

        return result;
    }

    void flag_picked_per_allele(std::vector<AnnotationWithMeta>& annotations) const {
        if (annotations.empty()) return;

        // Group by allele
        std::map<std::string, std::vector<AnnotationWithMeta*> > by_allele;
        for (size_t i = 0; i < annotations.size(); ++i) {
            std::string key = annotations[i].annotation.ref_allele + ">" +
                              annotations[i].annotation.alt_allele;
            by_allele[key].push_back(&annotations[i]);
        }

        // Pick best for each allele
        for (auto& pair : by_allele) {
            if (pair.second.empty()) continue;

            std::vector<AnnotationWithMeta> temp;
            for (auto* ptr : pair.second) {
                temp.push_back(*ptr);
            }

            std::sort(temp.begin(), temp.end(),
                [this](const AnnotationWithMeta& a, const AnnotationWithMeta& b) {
                    return is_better(a, b);
                });

            // Find and flag the best one in the original vector
            for (size_t i = 0; i < annotations.size(); ++i) {
                if (annotations[i].annotation.ref_allele == temp[0].annotation.ref_allele &&
                    annotations[i].annotation.alt_allele == temp[0].annotation.alt_allele &&
                    annotations[i].annotation.transcript_id == temp[0].annotation.transcript_id) {
                    annotations[i].is_picked = true;
                    break;
                }
            }
        }
    }

    void flag_picked_per_allele_gene(std::vector<AnnotationWithMeta>& annotations) const {
        if (annotations.empty()) return;

        // Group by allele and gene
        std::map<std::string, std::vector<AnnotationWithMeta*> > by_allele_gene;
        for (size_t i = 0; i < annotations.size(); ++i) {
            std::string key = annotations[i].annotation.ref_allele + ">" +
                              annotations[i].annotation.alt_allele + ":" +
                              annotations[i].annotation.gene_symbol;
            by_allele_gene[key].push_back(&annotations[i]);
        }

        // Pick best for each allele+gene
        for (auto& pair : by_allele_gene) {
            if (pair.second.empty()) continue;

            std::vector<AnnotationWithMeta> temp;
            for (auto* ptr : pair.second) {
                temp.push_back(*ptr);
            }

            std::sort(temp.begin(), temp.end(),
                [this](const AnnotationWithMeta& a, const AnnotationWithMeta& b) {
                    return is_better(a, b);
                });

            // Find and flag the best one in the original vector
            for (size_t i = 0; i < annotations.size(); ++i) {
                if (annotations[i].annotation.ref_allele == temp[0].annotation.ref_allele &&
                    annotations[i].annotation.alt_allele == temp[0].annotation.alt_allele &&
                    annotations[i].annotation.gene_symbol == temp[0].annotation.gene_symbol &&
                    annotations[i].annotation.transcript_id == temp[0].annotation.transcript_id) {
                    annotations[i].is_picked = true;
                    break;
                }
            }
        }
    }
};

/**
 * Parse pick order from string
 */
inline std::vector<PickCriteria> parse_pick_order(const std::string& order_str) {
    std::vector<PickCriteria> result;

    std::istringstream iss(order_str);
    std::string token;

    while (std::getline(iss, token, ',')) {
        // Trim whitespace
        size_t start = token.find_first_not_of(" \t");
        size_t end = token.find_last_not_of(" \t");
        if (start != std::string::npos && end != std::string::npos) {
            token = token.substr(start, end - start + 1);
        }

        // Convert to lowercase
        for (size_t i = 0; i < token.size(); ++i) {
            token[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(token[i])));
        }

        if (token == "canonical") result.push_back(PickCriteria::CANONICAL);
        else if (token == "mane" || token == "mane_select") result.push_back(PickCriteria::MANE_SELECT);
        else if (token == "mane_plus") result.push_back(PickCriteria::MANE_PLUS);
        else if (token == "appris") result.push_back(PickCriteria::APPRIS);
        else if (token == "tsl") result.push_back(PickCriteria::TSL);
        else if (token == "biotype") result.push_back(PickCriteria::BIOTYPE);
        else if (token == "ccds") result.push_back(PickCriteria::CCDS);
        else if (token == "rank") result.push_back(PickCriteria::RANK);
        else if (token == "length") result.push_back(PickCriteria::LENGTH);
    }

    return result;
}

/**
 * Parse biotypes from comma-separated string
 */
inline std::set<std::string> parse_biotypes(const std::string& biotypes_str) {
    std::set<std::string> result;

    std::istringstream iss(biotypes_str);
    std::string token;

    while (std::getline(iss, token, ',')) {
        size_t start = token.find_first_not_of(" \t");
        size_t end = token.find_last_not_of(" \t");
        if (start != std::string::npos && end != std::string::npos) {
            token = token.substr(start, end - start + 1);
            if (!token.empty()) {
                result.insert(token);
            }
        }
    }

    return result;
}

} // namespace vep

#endif // TRANSCRIPT_FILTER_HPP
