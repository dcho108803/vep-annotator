/**
 * Annotation Source Interface
 *
 * Base classes for all annotation sources (pathogenicity, conservation,
 * regulatory, etc.). Provides a unified interface for querying annotations.
 */

#ifndef ANNOTATION_SOURCE_HPP
#define ANNOTATION_SOURCE_HPP

#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <optional>
#include <mutex>
#include <shared_mutex>

namespace vep {

// Forward declarations
struct Transcript;
struct VariantAnnotation;

/**
 * Abstract base class for all annotation sources
 */
class AnnotationSource {
public:
    virtual ~AnnotationSource() = default;

    /**
     * Get the source name (e.g., "dbnsfp", "spliceai", "phylop")
     */
    virtual std::string name() const = 0;

    /**
     * Get the type of source for CLI help
     * (e.g., "pathogenicity", "conservation", "splice", "regulatory")
     */
    virtual std::string type() const = 0;

    /**
     * Get a description of this source
     */
    virtual std::string description() const = 0;

    /**
     * Check if the source is initialized and ready
     */
    virtual bool is_ready() const = 0;

    /**
     * Initialize the source (lazy loading)
     * Called automatically on first use if not manually initialized
     */
    virtual void initialize() = 0;

    /**
     * Annotate a variant, adding results to the annotations map
     * @param chrom Chromosome
     * @param pos Position (1-based)
     * @param ref Reference allele
     * @param alt Alternate allele
     * @param transcript Optional transcript context (for transcript-specific annotations)
     * @param annotations Output map to populate with "source:field" -> "value" pairs
     */
    virtual void annotate(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript* transcript,
        std::map<std::string, std::string>& annotations
    ) = 0;

    /**
     * Get list of fields this source provides
     */
    virtual std::vector<std::string> get_fields() const = 0;

    /**
     * Get memory usage in bytes (for stats/debugging)
     */
    virtual size_t memory_usage() const { return 0; }

    /**
     * Thread-safe check - if false, source will be locked during annotation
     */
    virtual bool is_thread_safe() const { return true; }

    /**
     * Get the data file path (for debugging/info)
     */
    virtual std::string get_data_path() const { return ""; }

protected:
    mutable std::mutex mutex_;  // For thread safety if needed

    /**
     * Ensure the source is initialized
     */
    void ensure_initialized() {
        if (!is_ready()) {
            std::lock_guard<std::mutex> lock(mutex_);
            if (!is_ready()) {
                initialize();
            }
        }
    }
};

/**
 * Score-based annotation source (for continuous scores like PhyloP)
 */
class ScoreAnnotationSource : public AnnotationSource {
public:
    /**
     * Get score at a specific position
     * @return Score value or nullopt if not available
     */
    virtual std::optional<double> get_score(
        const std::string& chrom,
        int pos
    ) const = 0;

    /**
     * Get scores for a range (for indels)
     * @return Vector of scores, may contain NaN for missing values
     */
    virtual std::vector<double> get_scores(
        const std::string& chrom,
        int start,
        int end
    ) const = 0;

    /**
     * Get aggregated score for a range (mean, max, min)
     */
    enum class Aggregation { MEAN, MAX, MIN, FIRST, LAST };

    virtual std::optional<double> get_aggregated_score(
        const std::string& chrom,
        int start,
        int end,
        Aggregation method = Aggregation::MEAN
    ) const;
};

/**
 * Interval-based annotation source (for BED, GFF3 features)
 */
class IntervalAnnotationSource : public AnnotationSource {
public:
    /**
     * Feature structure for interval-based annotations
     */
    struct Feature {
        std::string chrom;
        int start;
        int end;
        std::string type;           // Feature type (e.g., "promoter", "enhancer")
        std::string id;             // Feature ID
        char strand = '.';
        std::map<std::string, std::string> attributes;
    };

    /**
     * Query overlapping features
     */
    virtual std::vector<Feature> query(
        const std::string& chrom,
        int start,
        int end
    ) const = 0;

    /**
     * Query overlapping features of a specific type
     */
    virtual std::vector<Feature> query_by_type(
        const std::string& chrom,
        int start,
        int end,
        const std::string& feature_type
    ) const = 0;
};

/**
 * Variant-based annotation source (for VCF-style annotations)
 * Matches by position and optionally by allele
 */
class VariantAnnotationSource : public AnnotationSource {
public:
    /**
     * Whether to require exact allele match (REF/ALT)
     */
    virtual bool requires_allele_match() const { return true; }

    /**
     * Query annotations for a specific variant
     */
    virtual std::map<std::string, std::string> query(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt
    ) const = 0;
};

/**
 * Annotation source manager - handles multiple sources
 */
class AnnotationSourceManager {
public:
    /**
     * Register an annotation source
     */
    void add_source(std::shared_ptr<AnnotationSource> source);

    /**
     * Get all registered sources
     */
    std::vector<std::shared_ptr<AnnotationSource>> get_sources() const;

    /**
     * Get source by name
     */
    std::shared_ptr<AnnotationSource> get_source(const std::string& name) const;

    /**
     * Enable/disable a source by name
     */
    void set_enabled(const std::string& name, bool enabled);

    /**
     * Check if a source is enabled
     */
    bool is_enabled(const std::string& name) const;

    /**
     * Initialize all sources
     */
    void initialize_all();

    /**
     * Annotate a variant with all enabled sources
     */
    void annotate_all(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript* transcript,
        std::map<std::string, std::string>& annotations
    );

    /**
     * Get all available fields from all sources
     */
    std::vector<std::pair<std::string, std::string>> get_all_fields() const;

    /**
     * Get memory usage of all sources
     */
    size_t total_memory_usage() const;

    /**
     * Get stats string
     */
    std::string get_stats() const;

private:
    std::vector<std::shared_ptr<AnnotationSource>> sources_;
    std::set<std::string> disabled_;
    mutable std::shared_mutex mutex_;
};

} // namespace vep

#endif // ANNOTATION_SOURCE_HPP
