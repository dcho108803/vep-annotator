/**
 * Structural Variant Support
 *
 * Handles structural variants: INS, DEL, DUP, INV, CNV, BND (breakends)
 */

#ifndef STRUCTURAL_VARIANT_HPP
#define STRUCTURAL_VARIANT_HPP

#include "vep_annotator.hpp"
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>

namespace vep {

/**
 * Structural variant types
 */
enum class SVType {
    SNV,        // Not an SV (regular variant)
    INS,        // Insertion
    DEL,        // Deletion
    DUP,        // Duplication
    TDUP,       // Tandem duplication
    INV,        // Inversion
    CNV,        // Copy number variant
    BND,        // Breakend (translocation)
    UNKNOWN
};

/**
 * Get string representation of SV type
 */
inline std::string sv_type_to_string(SVType type) {
    if (type == SVType::SNV) return "SNV";
    if (type == SVType::INS) return "INS";
    if (type == SVType::DEL) return "DEL";
    if (type == SVType::DUP) return "DUP";
    if (type == SVType::TDUP) return "TDUP";
    if (type == SVType::INV) return "INV";
    if (type == SVType::CNV) return "CNV";
    if (type == SVType::BND) return "BND";
    return "UNKNOWN";
}

/**
 * Parse SV type from string
 */
inline SVType parse_sv_type(const std::string& type_str) {
    std::string upper = type_str;
    for (size_t i = 0; i < upper.size(); ++i) {
        upper[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(upper[i])));
    }

    // Handle VCF symbolic alleles like <DEL>, <DUP>, etc.
    std::string clean = upper;
    if (clean.size() > 2 && clean[0] == '<' && clean[clean.size()-1] == '>') {
        clean = clean.substr(1, clean.size() - 2);
    }

    if (clean == "INS") return SVType::INS;
    if (clean == "DEL") return SVType::DEL;
    if (clean == "DUP") return SVType::DUP;
    if (clean == "TDUP" || clean == "DUP:TANDEM") return SVType::TDUP;
    if (clean == "INV") return SVType::INV;
    if (clean == "CNV") return SVType::CNV;
    if (clean == "BND") return SVType::BND;

    // Handle copy number values: CN0, CN1, CN2, etc.
    if (clean.size() >= 2 && clean.substr(0, 2) == "CN") {
        return SVType::CNV;
    }

    return SVType::UNKNOWN;
}

/**
 * Structural variant representation
 */
struct StructuralVariant {
    std::string chromosome;
    int start = 0;              // 1-based start position
    int end = 0;                // 1-based end position (inclusive)
    std::string ref_allele;
    std::string alt_allele;
    SVType sv_type = SVType::UNKNOWN;
    int sv_len = 0;             // Length of the SV (positive for INS, negative for DEL)
    int copy_number = -1;       // Copy number for CNV (-1 = unknown)

    // Breakend information (for BND)
    std::string bnd_mate_chrom;
    int bnd_mate_pos = 0;
    bool bnd_mate_forward = true;

    // Computed properties
    int length() const {
        if (sv_len != 0) return std::abs(sv_len);
        return end - start + 1;
    }

    bool is_sv() const {
        return sv_type != SVType::SNV && sv_type != SVType::UNKNOWN;
    }

    bool overlaps(int s, int e) const {
        return start <= e && end >= s;
    }

    // Check if this SV completely contains a region
    bool contains(int s, int e) const {
        return start <= s && end >= e;
    }
};

/**
 * SV-specific consequence types
 */
inline std::vector<ConsequenceType> get_sv_consequences(
    const StructuralVariant& sv,
    const Transcript& transcript) {

    std::vector<ConsequenceType> consequences;

    // Check if SV overlaps the transcript at all
    if (!sv.overlaps(transcript.start, transcript.end)) {
        // Check upstream/downstream
        if (sv.end < transcript.start) {
            int distance = transcript.start - sv.end;
            if (distance <= 5000) {
                if (transcript.strand == '+') {
                    consequences.push_back(ConsequenceType::UPSTREAM_GENE_VARIANT);
                } else {
                    consequences.push_back(ConsequenceType::DOWNSTREAM_GENE_VARIANT);
                }
            }
        } else if (sv.start > transcript.end) {
            int distance = sv.start - transcript.end;
            if (distance <= 5000) {
                if (transcript.strand == '+') {
                    consequences.push_back(ConsequenceType::DOWNSTREAM_GENE_VARIANT);
                } else {
                    consequences.push_back(ConsequenceType::UPSTREAM_GENE_VARIANT);
                }
            }
        }

        if (consequences.empty()) {
            consequences.push_back(ConsequenceType::INTERGENIC_VARIANT);
        }
        return consequences;
    }

    // Check for transcript ablation (entire transcript deleted)
    if (sv.sv_type == SVType::DEL && sv.contains(transcript.start, transcript.end)) {
        consequences.push_back(ConsequenceType::TRANSCRIPT_ABLATION);
        return consequences;
    }

    // Check for transcript amplification (entire transcript duplicated)
    if ((sv.sv_type == SVType::DUP || sv.sv_type == SVType::TDUP) &&
        sv.contains(transcript.start, transcript.end)) {
        consequences.push_back(ConsequenceType::TRANSCRIPT_AMPLIFICATION);
        return consequences;
    }

    // Check for feature elongation (insertion extends the feature)
    if (sv.sv_type == SVType::INS && sv.overlaps(transcript.start, transcript.end)) {
        // INS within a transcript extends it
        consequences.push_back(ConsequenceType::FEATURE_ELONGATION);
    }

    // Check for feature truncation (partial deletion)
    if (sv.sv_type == SVType::DEL && sv.overlaps(transcript.start, transcript.end) &&
        !sv.contains(transcript.start, transcript.end)) {
        consequences.push_back(ConsequenceType::FEATURE_TRUNCATION);
    }

    // Check for coding sequence effects
    if (transcript.is_coding()) {
        bool affects_cds = false;
        bool affects_utr5 = false;
        bool affects_utr3 = false;
        bool affects_splice = false;

        // Check each CDS region for coding impact
        for (const auto& cds : transcript.cds_regions) {
            if (sv.overlaps(cds.start, cds.end)) {
                affects_cds = true;
            }
        }

        // Check splice sites at exon boundaries (not CDS boundaries)
        for (size_t i = 0; i < transcript.exons.size(); ++i) {
            const auto& exon = transcript.exons[i];
            if (i > 0 && sv.overlaps(exon.start - 2, exon.start + 1)) {
                affects_splice = true;
            }
            if (i + 1 < transcript.exons.size() && sv.overlaps(exon.end - 1, exon.end + 2)) {
                affects_splice = true;
            }
        }

        // Check UTRs
        if (transcript.strand == '+') {
            if (sv.overlaps(transcript.start, transcript.cds_start - 1)) {
                affects_utr5 = true;
            }
            if (sv.overlaps(transcript.cds_end + 1, transcript.end)) {
                affects_utr3 = true;
            }
        } else {
            if (sv.overlaps(transcript.cds_end + 1, transcript.end)) {
                affects_utr5 = true;
            }
            if (sv.overlaps(transcript.start, transcript.cds_start - 1)) {
                affects_utr3 = true;
            }
        }

        // Determine consequences based on SV type and affected regions
        if (sv.sv_type == SVType::DEL) {
            if (affects_splice) {
                consequences.push_back(ConsequenceType::SPLICE_DONOR_VARIANT);
            }
            if (affects_cds) {
                // Check if frameshift
                int deleted_bp = sv.length();
                if (deleted_bp % 3 != 0) {
                    consequences.push_back(ConsequenceType::FRAMESHIFT_VARIANT);
                } else {
                    consequences.push_back(ConsequenceType::INFRAME_DELETION);
                }
            }
        } else if (sv.sv_type == SVType::INS || sv.sv_type == SVType::DUP || sv.sv_type == SVType::TDUP) {
            if (affects_splice) {
                consequences.push_back(ConsequenceType::SPLICE_REGION_VARIANT);
            }
            if (affects_cds) {
                if (sv.sv_type == SVType::INS && sv.sv_len == 0) {
                    // Unknown insertion length - can't determine frameshift/inframe
                    consequences.push_back(ConsequenceType::CODING_SEQUENCE_VARIANT);
                } else {
                    int inserted_bp = sv.length();
                    if (inserted_bp % 3 != 0) {
                        consequences.push_back(ConsequenceType::FRAMESHIFT_VARIANT);
                    } else {
                        consequences.push_back(ConsequenceType::INFRAME_INSERTION);
                    }
                }
            }
        } else if (sv.sv_type == SVType::INV) {
            if (affects_cds) {
                consequences.push_back(ConsequenceType::CODING_SEQUENCE_VARIANT);
            }
        } else if (sv.sv_type == SVType::CNV) {
            if (sv.copy_number == 0) {
                // Complete deletion
                if (sv.contains(transcript.start, transcript.end)) {
                    consequences.push_back(ConsequenceType::TRANSCRIPT_ABLATION);
                } else if (affects_cds) {
                    consequences.push_back(ConsequenceType::FRAMESHIFT_VARIANT);
                }
            } else if (sv.copy_number > 2) {
                // Duplication/amplification
                if (affects_cds) {
                    consequences.push_back(ConsequenceType::CODING_SEQUENCE_VARIANT);
                }
            }
        }

        // Add UTR variants
        if (affects_utr5) {
            consequences.push_back(ConsequenceType::FIVE_PRIME_UTR_VARIANT);
        }
        if (affects_utr3) {
            consequences.push_back(ConsequenceType::THREE_PRIME_UTR_VARIANT);
        }

    } else {
        // Non-coding transcript
        consequences.push_back(ConsequenceType::NON_CODING_TRANSCRIPT_VARIANT);
    }

    // Check for intron variants if no exonic effects
    if (consequences.empty()) {
        bool in_exon = false;
        for (const auto& exon : transcript.exons) {
            if (sv.overlaps(exon.start, exon.end)) {
                in_exon = true;
                break;
            }
        }

        if (in_exon) {
            if (transcript.is_coding()) {
                consequences.push_back(ConsequenceType::CODING_SEQUENCE_VARIANT);
            } else {
                consequences.push_back(ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT);
            }
        } else {
            consequences.push_back(ConsequenceType::INTRON_VARIANT);
        }
    }

    return consequences;
}

/**
 * Parse structural variant from VCF record
 */
inline StructuralVariant parse_sv_from_vcf(
    const std::string& chrom,
    int pos,
    const std::string& ref,
    const std::string& alt,
    const std::map<std::string, std::string>& info) {

    StructuralVariant sv;
    sv.chromosome = chrom;
    sv.start = pos;
    sv.ref_allele = ref;
    sv.alt_allele = alt;

    // Try to get SV type from SVTYPE INFO field
    auto svtype_it = info.find("SVTYPE");
    if (svtype_it != info.end()) {
        sv.sv_type = parse_sv_type(svtype_it->second);
    } else {
        // Try to infer from ALT allele
        sv.sv_type = parse_sv_type(alt);
    }

    // Get END position
    auto end_it = info.find("END");
    if (end_it != info.end()) {
        try {
            sv.end = std::stoi(end_it->second);
        } catch (...) {
            sv.end = pos;
        }
    } else {
        // Estimate end from ref/alt lengths
        if (sv.sv_type == SVType::DEL) {
            sv.end = pos + static_cast<int>(ref.size()) - 1;
        } else {
            sv.end = pos;
        }
    }

    // Get SVLEN
    auto svlen_it = info.find("SVLEN");
    if (svlen_it != info.end()) {
        try {
            sv.sv_len = std::stoi(svlen_it->second);
        } catch (...) {}
        // If END was not found, derive it from SVLEN
        if (end_it == info.end() && sv.sv_len != 0) {
            if (sv.sv_type == SVType::DEL) {
                sv.end = pos + std::abs(sv.sv_len) - 1;
            } else if (sv.sv_type == SVType::DUP || sv.sv_type == SVType::INV) {
                sv.end = pos + std::abs(sv.sv_len) - 1;
            }
        }
    } else {
        // Calculate from positions
        if (sv.sv_type == SVType::DEL) {
            sv.sv_len = -(sv.end - sv.start + 1);
        } else if (sv.sv_type == SVType::DUP) {
            sv.sv_len = sv.end - sv.start + 1;
        } else if (sv.sv_type == SVType::INS && sv.end != sv.start) {
            // Only calculate INS length when END was explicitly provided
            sv.sv_len = sv.end - sv.start + 1;
        }
        // For symbolic INS without SVLEN and without END, sv_len stays 0 (unknown)
    }

    // Get copy number for CNV
    auto cn_it = info.find("CN");
    if (cn_it != info.end()) {
        try {
            sv.copy_number = std::stoi(cn_it->second);
        } catch (...) {}
    } else {
        // Try to parse from alt allele (e.g., <CN0>, <CN3>)
        if (alt.size() > 3 && alt.substr(0, 3) == "<CN") {
            try {
                sv.copy_number = std::stoi(alt.substr(3, alt.size() - 4));
            } catch (...) {}
        }
    }

    // Handle BND (breakend) notation
    if (sv.sv_type == SVType::BND || alt.find('[') != std::string::npos || alt.find(']') != std::string::npos) {
        sv.sv_type = SVType::BND;
        // Parse BND notation: N[chr:pos[ or ]chr:pos]N
        size_t bracket_pos = alt.find('[');
        if (bracket_pos == std::string::npos) {
            bracket_pos = alt.find(']');
        }

        if (bracket_pos != std::string::npos) {
            size_t colon_pos = alt.find(':', bracket_pos);
            size_t end_bracket = alt.find_first_of("[]", bracket_pos + 1);

            if (colon_pos != std::string::npos && end_bracket != std::string::npos) {
                sv.bnd_mate_chrom = alt.substr(bracket_pos + 1, colon_pos - bracket_pos - 1);
                try {
                    sv.bnd_mate_pos = std::stoi(alt.substr(colon_pos + 1, end_bracket - colon_pos - 1));
                } catch (...) {}
            }

            sv.bnd_mate_forward = (alt.find('[') != std::string::npos);
        }
    }

    // If still no type, check if it's a regular indel
    if (sv.sv_type == SVType::UNKNOWN) {
        if (ref.size() > alt.size()) {
            sv.sv_type = SVType::DEL;
            sv.sv_len = static_cast<int>(ref.size()) - static_cast<int>(alt.size());
            sv.end = pos + static_cast<int>(ref.size()) - 1;
        } else if (alt.size() > ref.size()) {
            sv.sv_type = SVType::INS;
            sv.sv_len = static_cast<int>(alt.size()) - static_cast<int>(ref.size());
            sv.end = pos;
        } else {
            sv.sv_type = SVType::SNV;
            sv.end = pos + static_cast<int>(ref.size()) - 1;
        }
    }

    return sv;
}

/**
 * Check if a variant is a structural variant (vs small indel)
 * Default threshold: 50bp
 */
inline bool is_structural_variant(const std::string& ref, const std::string& alt, int threshold = 50) {
    // Check for symbolic alleles
    if (!alt.empty() && alt[0] == '<') {
        return true;
    }

    // Check for BND notation
    if (alt.find('[') != std::string::npos || alt.find(']') != std::string::npos) {
        return true;
    }

    // Check length
    int len_diff = std::abs(static_cast<int>(ref.size()) - static_cast<int>(alt.size()));
    return len_diff >= threshold;
}

/**
 * Calculate overlap percentage between SV and a region
 */
inline double calculate_overlap_percentage(
    int sv_start, int sv_end,
    int region_start, int region_end,
    bool reciprocal = false) {

    int overlap_start = std::max(sv_start, region_start);
    int overlap_end = std::min(sv_end, region_end);

    if (overlap_start > overlap_end) {
        return 0.0;
    }

    int overlap_len = overlap_end - overlap_start + 1;
    int sv_len = sv_end - sv_start + 1;
    int region_len = region_end - region_start + 1;

    double sv_overlap = static_cast<double>(overlap_len) / sv_len;
    double region_overlap = static_cast<double>(overlap_len) / region_len;

    if (reciprocal) {
        return std::min(sv_overlap, region_overlap);
    }

    return sv_overlap;
}

/**
 * Overlap types for custom annotations
 */
enum class OverlapType {
    ANY,        // Any overlap
    WITHIN,     // SV completely within region
    SURROUNDING, // SV completely surrounds region
    EXACT,      // Exact position match
    RECIPROCAL  // Reciprocal overlap percentage
};

inline OverlapType parse_overlap_type(const std::string& type_str) {
    std::string lower = type_str;
    for (size_t i = 0; i < lower.size(); ++i) {
        lower[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(lower[i])));
    }

    if (lower == "within") return OverlapType::WITHIN;
    if (lower == "surrounding") return OverlapType::SURROUNDING;
    if (lower == "exact") return OverlapType::EXACT;
    if (lower == "reciprocal") return OverlapType::RECIPROCAL;
    return OverlapType::ANY;
}

/**
 * Get SV regulatory consequences based on overlap with regulatory features
 * @param sv Structural variant
 * @param feature_type Type of regulatory feature (e.g., "promoter", "enhancer", "TF_binding_site")
 * @param feature_start Start of the regulatory feature
 * @param feature_end End of the regulatory feature
 * @return Vector of regulatory consequence types
 */
inline std::vector<ConsequenceType> get_sv_regulatory_consequences(
    const StructuralVariant& sv,
    const std::string& feature_type,
    int feature_start, int feature_end) {

    std::vector<ConsequenceType> consequences;

    if (!sv.overlaps(feature_start, feature_end)) {
        return consequences;
    }

    bool is_tfbs = (feature_type == "TF_binding_site" || feature_type == "TFBS" ||
                    feature_type.find("TF_binding") != std::string::npos);

    if (sv.sv_type == SVType::DEL) {
        if (sv.contains(feature_start, feature_end)) {
            // Complete deletion of the feature
            if (is_tfbs) {
                consequences.push_back(ConsequenceType::TFBS_ABLATION);
            } else {
                consequences.push_back(ConsequenceType::REGULATORY_REGION_ABLATION);
            }
        } else {
            // Partial overlap
            if (is_tfbs) {
                consequences.push_back(ConsequenceType::TF_BINDING_SITE_VARIANT);
            } else {
                consequences.push_back(ConsequenceType::REGULATORY_REGION_VARIANT);
            }
        }
    } else if (sv.sv_type == SVType::DUP || sv.sv_type == SVType::TDUP) {
        if (sv.contains(feature_start, feature_end)) {
            // Complete duplication of the feature
            if (is_tfbs) {
                consequences.push_back(ConsequenceType::TFBS_AMPLIFICATION);
            } else {
                consequences.push_back(ConsequenceType::REGULATORY_REGION_AMPLIFICATION);
            }
        } else {
            if (is_tfbs) {
                consequences.push_back(ConsequenceType::TF_BINDING_SITE_VARIANT);
            } else {
                consequences.push_back(ConsequenceType::REGULATORY_REGION_VARIANT);
            }
        }
    } else {
        // Other SV types overlapping regulatory features
        if (is_tfbs) {
            consequences.push_back(ConsequenceType::TF_BINDING_SITE_VARIANT);
        } else {
            consequences.push_back(ConsequenceType::REGULATORY_REGION_VARIANT);
        }
    }

    return consequences;
}

} // namespace vep

#endif // STRUCTURAL_VARIANT_HPP
