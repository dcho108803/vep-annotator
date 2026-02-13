/**
 * Exon/Intron Numbering
 *
 * Provides exon and intron numbering for variant annotations.
 * Supports the --numbers flag functionality.
 */

#ifndef EXON_INTRON_NUMBERS_HPP
#define EXON_INTRON_NUMBERS_HPP

#include <string>
#include <vector>
#include <algorithm>

namespace vep {

// Forward declaration - actual Transcript struct should be in vep_annotator.hpp
struct TranscriptRegion {
    int start;
    int end;
    int number;  // 1-based exon/intron number
};

/**
 * Result of exon/intron lookup
 */
struct ExonIntronInfo {
    bool found = false;
    bool is_exon = false;
    int number = 0;           // 1-based number
    int total_exons = 0;      // Total number of exons
    int total_introns = 0;    // Total number of introns

    // Position within exon/intron (1-based)
    int position_in_feature = 0;
    int feature_length = 0;

    // Phase information
    int phase = -1;           // CDS phase at exon start (0, 1, 2, or -1)

    std::string to_string() const {
        if (!found) return "";

        std::string result;
        if (is_exon) {
            result = std::to_string(number) + "/" + std::to_string(total_exons);
        } else {
            result = std::to_string(number) + "/" + std::to_string(total_introns);
        }
        return result;
    }

    std::string feature_type() const {
        if (!found) return "";
        return is_exon ? "exon" : "intron";
    }
};

/**
 * Calculate exon number for a position in a transcript
 *
 * @param position Genomic position (1-based)
 * @param exon_starts Vector of exon start positions (1-based)
 * @param exon_ends Vector of exon end positions (1-based)
 * @param strand '+' or '-'
 * @return ExonIntronInfo with numbering information
 */
inline ExonIntronInfo get_exon_intron_number(
    int position,
    const std::vector<int>& exon_starts,
    const std::vector<int>& exon_ends,
    char strand) {

    ExonIntronInfo info;

    if (exon_starts.empty() || exon_starts.size() != exon_ends.size()) {
        return info;
    }

    int num_exons = static_cast<int>(exon_starts.size());
    info.total_exons = num_exons;
    info.total_introns = num_exons > 1 ? num_exons - 1 : 0;

    // Exons are already sorted by start position from GTF loading;
    // use the input arrays directly without copying or sorting.

    // Check each exon
    for (int i = 0; i < num_exons; ++i) {
        int start = exon_starts[i];
        int end = exon_ends[i];

        if (position >= start && position <= end) {
            // Position is in this exon
            info.found = true;
            info.is_exon = true;

            // Exon number depends on strand
            if (strand == '-') {
                info.number = num_exons - i;
            } else {
                info.number = i + 1;
            }

            info.position_in_feature = position - start + 1;
            info.feature_length = end - start + 1;
            return info;
        }

        // Check if in intron (between this exon and next)
        if (i < num_exons - 1) {
            int intron_start = end + 1;
            int intron_end = exon_starts[i + 1] - 1;

            if (position >= intron_start && position <= intron_end) {
                // Position is in this intron
                info.found = true;
                info.is_exon = false;

                // Intron number depends on strand
                if (strand == '-') {
                    info.number = num_exons - i - 1;
                } else {
                    info.number = i + 1;
                }

                info.position_in_feature = position - intron_start + 1;
                info.feature_length = intron_end - intron_start + 1;
                return info;
            }
        }
    }

    return info;
}

/**
 * Format exon number string (e.g., "5/12")
 */
inline std::string format_exon_number(const ExonIntronInfo& info) {
    if (!info.found || !info.is_exon) {
        return "";
    }
    return std::to_string(info.number) + "/" + std::to_string(info.total_exons);
}

/**
 * Format intron number string (e.g., "4/11")
 */
inline std::string format_intron_number(const ExonIntronInfo& info) {
    if (!info.found || info.is_exon) {
        return "";
    }
    return std::to_string(info.number) + "/" + std::to_string(info.total_introns);
}

/**
 * Get CDS exon number (exons that overlap CDS)
 *
 * @param position Genomic position (1-based)
 * @param exon_starts Vector of exon start positions
 * @param exon_ends Vector of exon end positions
 * @param cds_start CDS start position
 * @param cds_end CDS end position
 * @param strand '+' or '-'
 */
inline ExonIntronInfo get_cds_exon_number(
    int position,
    const std::vector<int>& exon_starts,
    const std::vector<int>& exon_ends,
    int cds_start,
    int cds_end,
    char strand) {

    ExonIntronInfo info;

    if (exon_starts.empty() || exon_starts.size() != exon_ends.size()) {
        return info;
    }

    // Count coding exons and find position
    std::vector<std::pair<int, int> > coding_exons;
    int num_exons = static_cast<int>(exon_starts.size());

    for (int i = 0; i < num_exons; ++i) {
        int ex_start = exon_starts[i];
        int ex_end = exon_ends[i];

        // Check if exon overlaps CDS
        if (ex_end >= cds_start && ex_start <= cds_end) {
            // Calculate coding portion of exon
            int coding_start = std::max(ex_start, cds_start);
            int coding_end = std::min(ex_end, cds_end);
            coding_exons.push_back(std::make_pair(coding_start, coding_end));
        }
    }

    if (coding_exons.empty()) {
        return info;
    }

    // coding_exons are already in sorted order because exon_starts/exon_ends
    // are sorted by start position from GTF loading.

    info.total_exons = static_cast<int>(coding_exons.size());
    info.total_introns = info.total_exons > 1 ? info.total_exons - 1 : 0;

    // Find position in coding exons
    for (size_t i = 0; i < coding_exons.size(); ++i) {
        int start = coding_exons[i].first;
        int end = coding_exons[i].second;

        if (position >= start && position <= end) {
            info.found = true;
            info.is_exon = true;

            if (strand == '-') {
                info.number = static_cast<int>(coding_exons.size()) - static_cast<int>(i);
            } else {
                info.number = static_cast<int>(i) + 1;
            }

            info.position_in_feature = position - start + 1;
            info.feature_length = end - start + 1;
            return info;
        }

        // Check intron between coding exons
        if (i + 1 < coding_exons.size()) {
            int intron_start = end + 1;
            int intron_end = coding_exons[i + 1].first - 1;

            if (position >= intron_start && position <= intron_end) {
                info.found = true;
                info.is_exon = false;

                if (strand == '-') {
                    info.number = static_cast<int>(coding_exons.size()) - static_cast<int>(i) - 1;
                } else {
                    info.number = static_cast<int>(i) + 1;
                }

                info.position_in_feature = position - intron_start + 1;
                info.feature_length = intron_end - intron_start + 1;
                return info;
            }
        }
    }

    return info;
}

/**
 * Calculate total CDS length
 */
inline int calculate_cds_length(
    const std::vector<int>& exon_starts,
    const std::vector<int>& exon_ends,
    int cds_start,
    int cds_end) {

    int total_length = 0;
    int num_exons = static_cast<int>(exon_starts.size());

    for (int i = 0; i < num_exons; ++i) {
        int ex_start = exon_starts[i];
        int ex_end = exon_ends[i];

        // Check if exon overlaps CDS
        if (ex_end >= cds_start && ex_start <= cds_end) {
            int coding_start = std::max(ex_start, cds_start);
            int coding_end = std::min(ex_end, cds_end);
            total_length += coding_end - coding_start + 1;
        }
    }

    return total_length;
}

/**
 * Calculate position within CDS
 */
inline int calculate_cds_position(
    int genomic_position,
    const std::vector<int>& exon_starts,
    const std::vector<int>& exon_ends,
    int cds_start,
    int cds_end,
    char strand) {

    // Create list of coding segments
    std::vector<std::pair<int, int> > coding_segments;
    int num_exons = static_cast<int>(exon_starts.size());

    for (int i = 0; i < num_exons; ++i) {
        int ex_start = exon_starts[i];
        int ex_end = exon_ends[i];

        if (ex_end >= cds_start && ex_start <= cds_end) {
            int coding_start = std::max(ex_start, cds_start);
            int coding_end = std::min(ex_end, cds_end);
            coding_segments.push_back(std::make_pair(coding_start, coding_end));
        }
    }

    if (coding_segments.empty()) {
        return -1;
    }

    // coding_segments are already in sorted order because exon_starts/exon_ends
    // are sorted by start position from GTF loading.

    if (strand == '-') {
        // Reverse for minus strand
        std::reverse(coding_segments.begin(), coding_segments.end());
    }

    int cds_pos = 0;
    for (size_t i = 0; i < coding_segments.size(); ++i) {
        int start = coding_segments[i].first;
        int end = coding_segments[i].second;

        if (strand == '-') {
            if (genomic_position >= start && genomic_position <= end) {
                return cds_pos + (end - genomic_position + 1);
            }
            cds_pos += end - start + 1;
        } else {
            if (genomic_position >= start && genomic_position <= end) {
                return cds_pos + (genomic_position - start + 1);
            }
            cds_pos += end - start + 1;
        }
    }

    return -1;
}

/**
 * Format distance to nearest splice site
 * Returns string like "+5" for 5bp into intron from donor,
 * or "-10" for 10bp before acceptor
 */
inline std::string format_splice_distance(
    int position,
    int exon_end,
    int next_exon_start) {

    // Distance from donor site (end of exon)
    int donor_dist = position - exon_end;

    // Distance to acceptor site (start of next exon)
    int acceptor_dist = next_exon_start - position;

    if (donor_dist > 0 && donor_dist <= acceptor_dist) {
        return "+" + std::to_string(donor_dist);
    } else if (acceptor_dist > 0) {
        return "-" + std::to_string(acceptor_dist);
    }

    return "";
}

} // namespace vep

#endif // EXON_INTRON_NUMBERS_HPP
