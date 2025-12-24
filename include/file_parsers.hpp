/**
 * File Format Parsers
 *
 * Utilities for reading various annotation file formats:
 * - TabixTSVReader: Tab-delimited files with tabix index
 * - BigWigReader: bigWig format for conservation scores
 * - GFF3Parser: GFF3 annotation files
 * - IntervalTree: Efficient range queries
 */

#ifndef FILE_PARSERS_HPP
#define FILE_PARSERS_HPP

#include <string>
#include <vector>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <functional>

namespace vep {

/**
 * Normalize chromosome name (remove "chr" prefix for consistency)
 */
inline std::string normalize_chrom(const std::string& chrom) {
    if (chrom.length() > 3 && chrom.substr(0, 3) == "chr") {
        return chrom.substr(3);
    }
    return chrom;
}

// ============================================================================
// Tabix TSV Reader
// ============================================================================

/**
 * Tab-delimited file reader with tabix index support
 * Used for dbNSFP, LoFtool, and other TSV annotation files
 */
class TabixTSVReader {
public:
    /**
     * Open a tabix-indexed TSV file
     * @param path Path to .tsv.gz or .txt.gz file (must have .tbi index)
     * @param chrom_col 0-based column index for chromosome
     * @param pos_col 0-based column index for position
     * @param columns Column names to extract (empty = use header)
     */
    TabixTSVReader(
        const std::string& path,
        int chrom_col = 0,
        int pos_col = 1,
        const std::vector<std::string>& columns = {}
    );

    ~TabixTSVReader();

    // Prevent copying
    TabixTSVReader(const TabixTSVReader&) = delete;
    TabixTSVReader& operator=(const TabixTSVReader&) = delete;

    /**
     * Query records at a specific position
     * @return Vector of row maps (column_name -> value)
     */
    std::vector<std::map<std::string, std::string>> query(
        const std::string& chrom,
        int pos
    );

    /**
     * Query records in a range
     */
    std::vector<std::map<std::string, std::string>> query_range(
        const std::string& chrom,
        int start,
        int end
    );

    /**
     * Get all column names
     */
    std::vector<std::string> get_columns() const;

    /**
     * Check if file is open and indexed
     */
    bool is_valid() const;

    /**
     * Get the file path
     */
    std::string get_path() const { return path_; }

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
    std::string path_;
};

// ============================================================================
// BigWig Reader
// ============================================================================

/**
 * bigWig file reader for conservation scores (PhyloP, PhastCons, GERP)
 * Requires libBigWig library
 */
class BigWigReader {
public:
    /**
     * Open a bigWig file
     * @param path Path to .bw or .bigWig file
     */
    explicit BigWigReader(const std::string& path);

    ~BigWigReader();

    // Prevent copying
    BigWigReader(const BigWigReader&) = delete;
    BigWigReader& operator=(const BigWigReader&) = delete;

    /**
     * Get value at a specific position
     * @return Score value or nullopt if not available
     */
    std::optional<double> get_value(const std::string& chrom, int pos) const;

    /**
     * Get values for a range
     * @return Vector of values, NaN for missing positions
     */
    std::vector<double> get_values(
        const std::string& chrom,
        int start,
        int end
    ) const;

    /**
     * Get mean value for a range
     */
    std::optional<double> get_mean(
        const std::string& chrom,
        int start,
        int end
    ) const;

    /**
     * Get max value for a range
     */
    std::optional<double> get_max(
        const std::string& chrom,
        int start,
        int end
    ) const;

    /**
     * Check if chromosome exists in file
     */
    bool has_chromosome(const std::string& chrom) const;

    /**
     * Get list of chromosomes in file
     */
    std::vector<std::string> get_chromosomes() const;

    /**
     * Check if file is valid
     */
    bool is_valid() const;

    /**
     * Get file path
     */
    std::string get_path() const { return path_; }

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
    std::string path_;
};

// ============================================================================
// GFF3 Parser
// ============================================================================

/**
 * GFF3 feature structure
 */
struct GFF3Feature {
    std::string seqid;      // Chromosome/contig
    std::string source;     // Data source
    std::string type;       // Feature type
    int start;              // 1-based start
    int end;                // 1-based end (inclusive)
    double score = -1;      // Score or -1 if not present
    char strand = '.';      // +, -, or .
    int phase = -1;         // 0, 1, 2 for CDS, -1 otherwise

    // Parsed attributes
    std::string id;
    std::string name;
    std::string parent;
    std::map<std::string, std::string> attributes;

    bool overlaps(int pos) const { return pos >= start && pos <= end; }
    bool overlaps(int s, int e) const { return s <= end && e >= start; }
};

/**
 * GFF3 file parser with efficient interval queries
 */
class GFF3Database {
public:
    /**
     * Load GFF3 file
     * @param path Path to .gff3 or .gff3.gz file
     * @param feature_types Feature types to load (empty = all)
     */
    explicit GFF3Database(
        const std::string& path,
        const std::set<std::string>& feature_types = {}
    );

    ~GFF3Database();

    /**
     * Query features overlapping a position
     */
    std::vector<const GFF3Feature*> query(
        const std::string& chrom,
        int pos
    ) const;

    /**
     * Query features overlapping a range
     */
    std::vector<const GFF3Feature*> query(
        const std::string& chrom,
        int start,
        int end
    ) const;

    /**
     * Query features by type
     */
    std::vector<const GFF3Feature*> query_by_type(
        const std::string& chrom,
        int start,
        int end,
        const std::string& type
    ) const;

    /**
     * Get all feature types loaded
     */
    std::set<std::string> get_feature_types() const;

    /**
     * Get total number of features
     */
    size_t feature_count() const;

    /**
     * Get file path
     */
    std::string get_path() const { return path_; }

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
    std::string path_;
};

// ============================================================================
// Interval Tree
// ============================================================================

/**
 * Interval tree for efficient range queries
 * Template parameter T is the data type to store with each interval
 */
template<typename T>
class IntervalTree {
public:
    IntervalTree() = default;

    /**
     * Insert an interval with associated data
     * @param start Interval start (inclusive)
     * @param end Interval end (inclusive)
     * @param data Data to associate with interval
     */
    void insert(int start, int end, T data);

    /**
     * Build the tree (call after all insertions)
     * Must be called before querying
     */
    void build();

    /**
     * Query intervals overlapping a point
     */
    std::vector<T> query(int point) const;

    /**
     * Query intervals overlapping a range
     */
    std::vector<T> query(int start, int end) const;

    /**
     * Check if tree is built
     */
    bool is_built() const { return built_; }

    /**
     * Get number of intervals
     */
    size_t size() const { return intervals_.size(); }

    /**
     * Clear all intervals
     */
    void clear();

private:
    struct Interval {
        int start;
        int end;
        T data;
    };

    std::vector<Interval> intervals_;
    bool built_ = false;

    // Tree structure for efficient queries
    struct Node {
        int center;
        std::vector<size_t> overlapping;  // Indices sorted by start
        std::unique_ptr<Node> left;
        std::unique_ptr<Node> right;
    };

    std::unique_ptr<Node> root_;

    std::unique_ptr<Node> build_tree(std::vector<size_t>& indices, int depth = 0);
    void query_node(const Node* node, int start, int end, std::vector<T>& results) const;
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * Parse a line into fields by delimiter
 */
std::vector<std::string> split_line(const std::string& line, char delim = '\t');

/**
 * Parse GFF3 attributes string
 */
std::map<std::string, std::string> parse_gff3_attributes(const std::string& attrs);

/**
 * URL decode a string (for GFF3 attribute values)
 */
std::string url_decode(const std::string& str);

/**
 * Check if a file exists
 */
bool file_exists(const std::string& path);

/**
 * Get file extension (handles .gz)
 */
std::string get_extension(const std::string& path);

} // namespace vep

#endif // FILE_PARSERS_HPP
