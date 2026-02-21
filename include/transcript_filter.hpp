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

// Case-insensitive consequence set lookup
inline bool consequence_set_contains(const std::set<std::string>& s, const std::string& term) {
    if (s.find(term) != s.end()) return true;
    // Try lowercase comparison
    std::string lower = term;
    for (auto& c : lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s.find(lower) != s.end();
}

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
    LENGTH,         // Transcript length (longer preferred)
    ENSEMBL,        // Prefer Ensembl transcripts
    REFSEQ          // Prefer RefSeq transcripts
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
    bool exclude_predicted = false; // Exclude XM_/XR_ predicted transcripts
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
        // Default pick order (matches Perl VEP: MANE before CANONICAL)
        pick_order.push_back(PickCriteria::MANE_SELECT);
        pick_order.push_back(PickCriteria::MANE_PLUS);
        pick_order.push_back(PickCriteria::CANONICAL);
        pick_order.push_back(PickCriteria::APPRIS);
        pick_order.push_back(PickCriteria::TSL);
        pick_order.push_back(PickCriteria::BIOTYPE);
        pick_order.push_back(PickCriteria::CCDS);
        pick_order.push_back(PickCriteria::RANK);
        pick_order.push_back(PickCriteria::LENGTH);
        freq_pop = "";  // Set by --filter-common (MAX_AF) or --freq-pop
    }
};

/**
 * Get consequence rank (lower = more severe)
 */
inline int get_consequence_rank(ConsequenceType csq) {
    // Lookup array indexed by ConsequenceType enum value for O(1) access.
    // Enum order: TRANSCRIPT_ABLATION(0) through UNKNOWN(38)
    static const int ranks[] = {
        1,   // TRANSCRIPT_ABLATION
        2,   // SPLICE_ACCEPTOR_VARIANT
        3,   // SPLICE_DONOR_VARIANT
        4,   // STOP_GAINED
        5,   // FRAMESHIFT_VARIANT
        6,   // STOP_LOST
        7,   // START_LOST
        8,   // TRANSCRIPT_AMPLIFICATION
        9,   // FEATURE_ELONGATION
        10,  // FEATURE_TRUNCATION
        11,  // INFRAME_INSERTION
        12,  // INFRAME_DELETION
        13,  // MISSENSE_VARIANT
        14,  // PROTEIN_ALTERING_VARIANT
        15,  // SPLICE_DONOR_5TH_BASE_VARIANT
        17,  // SPLICE_DONOR_REGION_VARIANT
        18,  // SPLICE_POLYPYRIMIDINE_TRACT_VARIANT
        16,  // SPLICE_REGION_VARIANT
        19,  // INCOMPLETE_TERMINAL_CODON_VARIANT
        20,  // START_RETAINED_VARIANT
        21,  // STOP_RETAINED_VARIANT
        22,  // SYNONYMOUS_VARIANT
        23,  // CODING_SEQUENCE_VARIANT
        31,  // CODING_TRANSCRIPT_VARIANT
        24,  // MATURE_MIRNA_VARIANT
        25,  // FIVE_PRIME_UTR_VARIANT
        26,  // THREE_PRIME_UTR_VARIANT
        27,  // NON_CODING_TRANSCRIPT_EXON_VARIANT
        28,  // INTRON_VARIANT
        29,  // NMD_TRANSCRIPT_VARIANT
        30,  // NON_CODING_TRANSCRIPT_VARIANT
        32,  // UPSTREAM_GENE_VARIANT
        33,  // DOWNSTREAM_GENE_VARIANT
        34,  // TFBS_ABLATION
        35,  // TFBS_AMPLIFICATION
        36,  // TF_BINDING_SITE_VARIANT
        37,  // REGULATORY_REGION_ABLATION
        38,  // REGULATORY_REGION_AMPLIFICATION
        39,  // REGULATORY_REGION_VARIANT
        40,  // INTERGENIC_VARIANT
        41,  // SEQUENCE_VARIANT
        50,  // UNKNOWN
    };
    int idx = static_cast<int>(csq);
    if (idx >= 0 && idx < static_cast<int>(sizeof(ranks) / sizeof(ranks[0])))
        return ranks[idx];
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
        if (s.find("principal_1") != std::string::npos) return 1;
        if (s.find("principal_2") != std::string::npos) return 2;
        if (s.find("principal_3") != std::string::npos) return 3;
        if (s.find("principal_4") != std::string::npos) return 4;
        if (s.find("principal_5") != std::string::npos) return 5;
        if (s.find("alternative_1") != std::string::npos) return 6;
        if (s.find("alternative_2") != std::string::npos) return 7;
        if (s.find("alternative_3") != std::string::npos) return 8;
        if (s.find("alternative_4") != std::string::npos) return 9;
        if (s.find("alternative_5") != std::string::npos) return 10;
        // Any other APPRIS value
        if (s.find("principal") != std::string::npos) return 5;
        if (s.find("alternative") != std::string::npos) return 10;
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

        if (annotations.empty()) return {};

        // Convert to extended annotations (move into metadata wrappers)
        std::vector<AnnotationWithMeta> extended;
        extended.reserve(annotations.size());
        for (const auto& ann : annotations) {
            extended.push_back(extend_annotation(ann));
        }

        // Apply basic filters (in-place via erase-remove)
        extended.erase(
            std::remove_if(extended.begin(), extended.end(),
                [this](const AnnotationWithMeta& ext) { return !passes_basic_filter(ext); }),
            extended.end());

        // Apply frequency filter if enabled (in-place)
        if (config_.check_frequency) {
            apply_frequency_filter_inplace(extended);
        }

        // Apply pick/selection logic
        if (config_.pick) {
            extended = pick_one(std::move(extended));
        } else if (config_.pick_allele) {
            extended = pick_per_allele(std::move(extended));
        } else if (config_.pick_allele_gene) {
            extended = pick_per_allele_gene(std::move(extended));
        } else if (config_.per_gene) {
            extended = pick_per_gene(std::move(extended));
        } else if (config_.most_severe) {
            extended = pick_most_severe(std::move(extended));
        } else if (config_.flag_pick) {
            flag_picked(extended);
        } else if (config_.flag_pick_allele) {
            flag_picked_per_allele(extended);
        } else if (config_.flag_pick_allele_gene) {
            flag_picked_per_allele_gene(extended);
        }

        // Convert back to VariantAnnotation (move out of wrappers)
        std::vector<VariantAnnotation> result;
        result.reserve(extended.size());
        for (auto& ext : extended) {
            if ((config_.flag_pick || config_.flag_pick_allele || config_.flag_pick_allele_gene) && ext.is_picked) {
                ext.annotation.custom_annotations["PICK"] = "1";
            }
            result.push_back(std::move(ext.annotation));
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

        auto mane_plus_it = ann.custom_annotations.find("MANE_PLUS_CLINICAL");
        if (mane_plus_it != ann.custom_annotations.end() && !mane_plus_it->second.empty()) {
            ext.is_mane_plus = true;
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

        auto ccds_it = ann.custom_annotations.find("CCDS");
        if (ccds_it != ann.custom_annotations.end() && !ccds_it->second.empty()) {
            ext.has_ccds = true;
        }

        auto tl_it = ann.custom_annotations.find("TRANSCRIPT_LENGTH");
        if (tl_it != ann.custom_annotations.end() && !tl_it->second.empty()) {
            try {
                ext.transcript_length = std::stoi(tl_it->second);
            } catch (...) {}
        }

        return ext;
    }

    bool passes_basic_filter(const AnnotationWithMeta& ext) const {
        const VariantAnnotation& ann = ext.annotation;

        // Canonical only
        if (config_.canonical_only && !ann.is_canonical) {
            return false;
        }

        // MANE only (Perl VEP: only MANE Select, not MANE Plus Clinical)
        if (config_.mane_only && !ext.is_mane_select) {
            return false;
        }

        // Exclude predicted (XM_/XR_) RefSeq transcripts
        if (config_.exclude_predicted && !ann.transcript_id.empty()) {
            if (ann.transcript_id.substr(0, 3) == "XM_" ||
                ann.transcript_id.substr(0, 3) == "XR_") {
                return false;
            }
        }

        // Coding only
        if (config_.coding_only && ann.biotype != "protein_coding") {
            return false;
        }

        // GENCODE basic filter
        if (config_.gencode_basic) {
            auto gb_it = ann.custom_annotations.find("GENCODE_BASIC");
            if (gb_it == ann.custom_annotations.end() || gb_it->second != "YES") {
                return false;
            }
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
                if (consequence_set_contains(config_.include_consequences,
                                             consequence_to_string(csq))) {
                    has_included = true;
                    break;
                }
            }
            if (!has_included) return false;
        }

        if (!config_.exclude_consequences.empty()) {
            // Only exclude if the most severe consequence is in the exclude set
            int best_rank = 999;
            std::string most_severe;
            for (const auto& csq : ann.consequences) {
                int r = get_consequence_rank(csq);
                if (r < best_rank) {
                    best_rank = r;
                    most_severe = consequence_to_string(csq);
                }
            }
            if (consequence_set_contains(config_.exclude_consequences, most_severe)) {
                return false;
            }
        }

        // Impact filter
        if (!config_.include_impacts.empty() &&
            config_.include_impacts.find(ann.impact) == config_.include_impacts.end()) {
            return false;
        }

        return true;
    }

    void apply_frequency_filter_inplace(
        std::vector<AnnotationWithMeta>& annotations) const {

        annotations.erase(
            std::remove_if(annotations.begin(), annotations.end(),
                [this](const AnnotationWithMeta& ext) {
                    auto it = ext.annotation.custom_annotations.find(config_.freq_pop);
                    if (it != ext.annotation.custom_annotations.end()) {
                        try {
                            double freq = std::stod(it->second);
                            bool passes = config_.freq_gt ?
                                (freq > config_.freq_threshold) :
                                (freq < config_.freq_threshold);
                            return !passes;
                        } catch (...) {}
                    }
                    return false;
                }),
            annotations.end());
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
            // Perl VEP biotype ranking hierarchy
            auto biotype_score = [](const std::string& bt) -> int {
                if (bt == "protein_coding") return 0;
                if (bt == "lncRNA" || bt == "lincRNA") return 1;
                if (bt == "miRNA") return 2;
                if (bt == "snRNA") return 3;
                if (bt == "snoRNA") return 4;
                if (bt == "rRNA") return 5;
                if (bt == "misc_RNA") return 6;
                if (bt == "nonsense_mediated_decay") return 7;
                if (bt.find("pseudogene") != std::string::npos) return 8;
                return 9;
            };
            return biotype_score(a.annotation.biotype) - biotype_score(b.annotation.biotype);
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
        if (criteria == PickCriteria::ENSEMBL) {
            // Prefer Ensembl source (source contains "Ensembl" or "ensembl")
            auto is_ensembl = [](const std::string& s) {
                return s.find("Ensembl") != std::string::npos || s.find("ensembl") != std::string::npos;
            };
            return (is_ensembl(b.annotation.source) ? 1 : 0) - (is_ensembl(a.annotation.source) ? 1 : 0);
        }
        if (criteria == PickCriteria::REFSEQ) {
            // Prefer RefSeq source
            auto is_refseq = [](const std::string& s) {
                return s.find("RefSeq") != std::string::npos || s.find("refseq") != std::string::npos;
            };
            return (is_refseq(b.annotation.source) ? 1 : 0) - (is_refseq(a.annotation.source) ? 1 : 0);
        }
        return 0;
    }

    std::vector<AnnotationWithMeta> pick_one(
        std::vector<AnnotationWithMeta> annotations) const {

        if (annotations.empty()) return annotations;

        std::sort(annotations.begin(), annotations.end(),
            [this](const AnnotationWithMeta& a, const AnnotationWithMeta& b) {
                return is_better(a, b);
            });

        std::vector<AnnotationWithMeta> result;
        result.push_back(std::move(annotations[0]));
        return result;
    }

    std::vector<AnnotationWithMeta> pick_per_allele(
        std::vector<AnnotationWithMeta> annotations) const {

        if (annotations.empty()) return annotations;

        std::map<std::string, std::vector<AnnotationWithMeta>> by_allele;
        for (auto& ann : annotations) {
            std::string key = ann.annotation.ref_allele + ">" + ann.annotation.alt_allele;
            by_allele[key].push_back(std::move(ann));
        }

        std::vector<AnnotationWithMeta> result;
        for (auto& pair : by_allele) {
            auto picked = pick_one(std::move(pair.second));
            if (!picked.empty()) {
                result.push_back(std::move(picked[0]));
            }
        }

        return result;
    }

    std::vector<AnnotationWithMeta> pick_per_gene(
        std::vector<AnnotationWithMeta> annotations) const {

        if (annotations.empty()) return annotations;

        std::map<std::string, std::vector<AnnotationWithMeta>> by_gene;
        for (auto& ann : annotations) {
            std::string gene_key = ann.annotation.gene_id.empty() ? ann.annotation.gene_symbol : ann.annotation.gene_id;
            by_gene[gene_key].push_back(std::move(ann));
        }

        std::vector<AnnotationWithMeta> result;
        for (auto& pair : by_gene) {
            auto picked = pick_one(std::move(pair.second));
            if (!picked.empty()) {
                result.push_back(std::move(picked[0]));
            }
        }

        return result;
    }

    std::vector<AnnotationWithMeta> pick_most_severe(
        std::vector<AnnotationWithMeta> annotations) const {

        if (annotations.empty()) return annotations;

        int min_rank = 100;
        for (const auto& ann : annotations) {
            int rank = ann.get_rank();
            if (rank < min_rank) min_rank = rank;
        }

        std::vector<AnnotationWithMeta> result;
        for (auto& ann : annotations) {
            if (ann.get_rank() == min_rank) {
                result.push_back(std::move(ann));
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

        // Group by allele and gene (prefer gene_id for uniqueness)
        std::map<std::string, std::vector<AnnotationWithMeta> > by_allele_gene;
        for (auto& ann : annotations) {
            std::string gene_key = ann.annotation.gene_id.empty() ? ann.annotation.gene_symbol : ann.annotation.gene_id;
            std::string key = ann.annotation.ref_allele + ">" + ann.annotation.alt_allele +
                              ":" + gene_key;
            by_allele_gene[key].push_back(std::move(ann));
        }

        std::vector<AnnotationWithMeta> result;
        for (auto& pair : by_allele_gene) {
            auto picked = pick_one(std::move(pair.second));
            if (!picked.empty()) {
                result.push_back(std::move(picked[0]));
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

        // Group by allele and gene (prefer gene_id for uniqueness)
        std::map<std::string, std::vector<AnnotationWithMeta*> > by_allele_gene;
        for (size_t i = 0; i < annotations.size(); ++i) {
            std::string gene_key = annotations[i].annotation.gene_id.empty() ?
                                   annotations[i].annotation.gene_symbol : annotations[i].annotation.gene_id;
            std::string key = annotations[i].annotation.ref_allele + ">" +
                              annotations[i].annotation.alt_allele + ":" +
                              gene_key;
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
                    annotations[i].annotation.gene_id == temp[0].annotation.gene_id &&
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
        else if (token == "mane_plus" || token == "mane_plus_clinical") result.push_back(PickCriteria::MANE_PLUS);
        else if (token == "appris") result.push_back(PickCriteria::APPRIS);
        else if (token == "tsl") result.push_back(PickCriteria::TSL);
        else if (token == "biotype") result.push_back(PickCriteria::BIOTYPE);
        else if (token == "ccds") result.push_back(PickCriteria::CCDS);
        else if (token == "rank") result.push_back(PickCriteria::RANK);
        else if (token == "length") result.push_back(PickCriteria::LENGTH);
        else if (token == "ensembl") result.push_back(PickCriteria::ENSEMBL);
        else if (token == "refseq") result.push_back(PickCriteria::REFSEQ);
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
