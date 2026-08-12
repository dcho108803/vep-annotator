/**
 * ACMG/ClinGen CNV classification (--acmg-cnv)
 *
 * Implements the computable subset of the ACMG/ClinGen 2020 technical
 * standard for CNV interpretation (Riggs et al., Genet Med 2020): the
 * evidence categories derivable from gene overlap, ClinGen dosage
 * sensitivity curations (--clingen-dosage) and population SV frequency
 * (--sv-frequency). Family/inheritance and literature evidence (sections
 * 4A-4N and 5) require case-level data and are not scored; classifications
 * here are a triage aid, not a substitute for manual curation.
 *
 * Scoring (loss = copy-number loss / deletion; gain = duplication):
 *   Section 1  content:      1A genes overlapped (0) / 1B none (-0.60)
 *   Section 2  dosage:       2A full overlap of an established HI/TS gene
 *                            (score 3) (+1.00); 2C-2E partial overlap of an
 *                            HI gene disrupting coding sequence (+0.90,
 *                            losses only); 2F CNV fully contained within a
 *                            dosage-insensitive gene (score 40) (-1.00)
 *   Section 3  gene count:   losses: 25-34 genes +0.45, >=35 +0.90
 *                            gains:  35-49 genes +0.45, >=50 +0.90
 *   Section 4  population:   4O co-located SV with AF >= 1% (-1.00)
 *
 * Classification thresholds (per the standard):
 *   >= 0.99 pathogenic; 0.90..0.98 likely_pathogenic;
 *   <= -0.99 benign; -0.98..-0.90 likely_benign; else uncertain_significance
 */

#ifndef ACMG_CNV_HPP
#define ACMG_CNV_HPP

#include <string>
#include <vector>
#include <map>
#include <cstdio>

#include "structural_variant.hpp"
#include "gene_constraint.hpp"

namespace vep {

struct AcmgCnvResult {
    bool applicable = false;        // false for INV/BND/INS or unknown-CN CNVs
    double score = 0.0;
    std::string classification;     // pathogenic .. benign
    std::vector<std::string> evidence;  // e.g. "2A:BRCA1(+1.00)"
};

inline std::string acmg_score_to_class(double s) {
    if (s >= 0.99) return "pathogenic";
    if (s >= 0.90) return "likely_pathogenic";
    if (s <= -0.99) return "benign";
    if (s <= -0.90) return "likely_benign";
    return "uncertain_significance";
}

inline AcmgCnvResult classify_acmg_cnv(
    const StructuralVariant& sv,
    const std::vector<const Transcript*>& transcripts,
    const GeneConstraintDB& constraint_db,
    double population_af = -1.0) {

    AcmgCnvResult res;

    // The rubric applies to copy-number losses and gains only
    bool loss = (sv.sv_type == SVType::DEL);
    bool gain = (sv.sv_type == SVType::DUP || sv.sv_type == SVType::TDUP);
    if (sv.sv_type == SVType::CNV && sv.copy_number >= 0) {
        if (sv.copy_number < 1.5) loss = true;
        else if (sv.copy_number > 2.0) gain = true;
    }
    if (!loss && !gain) return res;
    res.applicable = true;

    auto add = [&](const std::string& label, double points) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%s(%+.2f)", label.c_str(), points);
        res.evidence.push_back(buf);
        res.score += points;
    };

    // Aggregate transcripts into genes (gene span approximated by the union
    // of its transcript spans)
    struct GeneSpan {
        int start = 0;
        int end = 0;
        std::string symbol;
        bool cds_overlap = false;
    };
    std::map<std::string, GeneSpan> genes;
    for (const Transcript* tx : transcripts) {
        if (!tx) continue;
        const std::string& key = tx->gene_id.empty() ? tx->gene_name : tx->gene_id;
        if (key.empty()) continue;
        auto it = genes.find(key);
        if (it == genes.end()) {
            GeneSpan g;
            g.start = tx->start;
            g.end = tx->end;
            g.symbol = tx->gene_name.empty() ? key : tx->gene_name;
            it = genes.emplace(key, g).first;
        } else {
            it->second.start = std::min(it->second.start, tx->start);
            it->second.end = std::max(it->second.end, tx->end);
        }
        const bool has_cds = tx->is_coding() || (tx->cds_start > 0 && tx->cds_end > 0);
        if (has_cds && sv.overlaps(tx->cds_start, tx->cds_end)) {
            it->second.cds_overlap = true;
        }
    }

    // Section 1: genomic content
    if (genes.empty()) {
        add("1B", -0.60);
        res.classification = acmg_score_to_class(res.score);
        return res;
    }
    add("1A", 0.0);

    // Section 2: established dosage sensitivity. The strongest applicable
    // pathogenic item is scored once; the benign item (2F) applies only when
    // no dosage-sensitive gene is involved.
    double best_dosage = 0.0;
    std::string best_label;
    std::string benign_container;
    for (const auto& kv : genes) {
        const GeneSpan& g = kv.second;
        GeneConstraint gc = constraint_db.get_by_symbol(g.symbol);
        if (!gc.has_data()) gc = constraint_db.get_by_gene_id(kv.first);
        double dosage = loss ? gc.hi_score : gc.ts_score;
        if (dosage == 3.0) {
            if (sv.contains(g.start, g.end)) {
                if (best_dosage < 1.0) {
                    best_dosage = 1.0;
                    best_label = "2A:" + g.symbol;
                }
            } else if (loss && g.cds_overlap && best_dosage < 0.90) {
                // Partial overlap of an established HI gene with coding
                // sequence involved (2C-2E collapsed)
                best_dosage = 0.90;
                best_label = "2C-2E:" + g.symbol;
            }
        }
        // Dosage-insensitive gene (ClinGen code 40) containing the whole CNV
        if (dosage == 40.0 && g.start <= sv.start && g.end >= sv.end) {
            benign_container = g.symbol;
        }
    }
    if (best_dosage > 0.0) {
        add(best_label, best_dosage);
    } else if (!benign_container.empty()) {
        add("2F:" + benign_container, -1.0);
    }

    // Section 3: number of protein-coding genes (approximated by all genes
    // with any overlapping transcript)
    size_t n = genes.size();
    if (loss) {
        if (n >= 35) add("3C", 0.90);
        else if (n >= 25) add("3B", 0.45);
        else add("3A", 0.0);
    } else {
        if (n >= 50) add("3C", 0.90);
        else if (n >= 35) add("3B", 0.45);
        else add("3A", 0.0);
    }

    // Section 4: population evidence — a co-located population SV at >= 1%
    // allele frequency is common variation (4O)
    if (population_af >= 0.01) {
        add("4O", -1.0);
    }

    res.classification = acmg_score_to_class(res.score);
    return res;
}

/**
 * Join evidence labels for a single output column value.
 */
inline std::string format_acmg_evidence(const std::vector<std::string>& evidence) {
    std::string out;
    for (size_t i = 0; i < evidence.size(); ++i) {
        if (i) out += "&";
        out += evidence[i];
    }
    return out;
}

} // namespace vep

#endif // ACMG_CNV_HPP
