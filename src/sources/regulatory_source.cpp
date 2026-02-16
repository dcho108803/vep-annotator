/**
 * Regulatory Annotation Source
 *
 * Provides regulatory annotations from Ensembl Regulatory Build GFF3 files.
 * Includes: promoters, enhancers, TFBS, open chromatin regions, CTCF binding sites.
 */

#include "annotation_source.hpp"
#include "file_parsers.hpp"
#include "vep_annotator.hpp"
#include <sstream>
#include <algorithm>

namespace vep {

/**
 * Regulatory Annotation Source
 *
 * Uses GFF3Database to query regulatory features overlapping variants.
 */
class RegulatorySource : public IntervalAnnotationSource {
public:
    RegulatorySource(const std::string& path,
                     const std::set<std::string>& cell_types = {})
        : path_(path), cell_types_(cell_types) {}

    std::string name() const override { return "regulatory"; }
    std::string type() const override { return "regulatory"; }
    std::string description() const override {
        return "Ensembl Regulatory Build annotations";
    }

    bool is_ready() const override { return db_ != nullptr; }

    void initialize() override {
        std::lock_guard<std::recursive_mutex> lock(mutex_);
        if (db_) return;

        log(LogLevel::INFO, "Loading regulatory annotations from: " + path_);

        // Load ALL regulatory features (no type filtering)
        // Cell type filtering is done at query time via attributes
        db_ = std::make_unique<GFF3Database>(path_);

        log(LogLevel::INFO, "Regulatory features loaded: " +
                           std::to_string(db_->feature_count()));
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
        if (!db_) return;

        (void)transcript;

        // Query region based on variant extent
        int start = pos;
        int end = pos + static_cast<int>(ref.length()) - 1;

        auto features = db_->query(chrom, start, end);

        if (features.empty()) return;

        // Collect feature types and IDs, optionally filtering by cell type
        std::set<std::string> types;
        std::vector<std::string> ids;
        std::vector<std::string> names;

        for (const auto* feature : features) {
            // If cell types are specified, filter by cell_type attribute
            if (!cell_types_.empty()) {
                auto ct_it = feature->attributes.find("cell_type");
                if (ct_it != feature->attributes.end()) {
                    bool match = false;
                    for (const auto& ct : cell_types_) {
                        if (ct_it->second.find(ct) != std::string::npos) {
                            match = true;
                            break;
                        }
                    }
                    if (!match) continue;
                } else {
                    // Feature has no cell_type attribute - skip when filtering by cell type
                    continue;
                }
            }

            types.insert(feature->type);
            if (!feature->id.empty()) {
                ids.push_back(feature->id);
            }
            if (!feature->name.empty()) {
                names.push_back(feature->name);
            }
        }

        // Format annotations
        if (!types.empty()) {
            std::ostringstream oss;
            bool first = true;
            for (const auto& t : types) {
                if (!first) oss << ",";
                oss << t;
                first = false;
            }
            annotations["regulatory:feature_type"] = oss.str();
        }

        if (!ids.empty()) {
            std::ostringstream oss;
            for (size_t i = 0; i < ids.size() && i < 5; ++i) {
                if (i > 0) oss << ",";
                oss << ids[i];
            }
            if (ids.size() > 5) oss << ",...";
            annotations["regulatory:feature_id"] = oss.str();
        }

        annotations["regulatory:count"] = std::to_string(ids.size());

        // Add specific feature type flags
        if (types.count("promoter") || types.count("Promoter")) {
            annotations["regulatory:in_promoter"] = "true";
        }
        if (types.count("enhancer") || types.count("Enhancer")) {
            annotations["regulatory:in_enhancer"] = "true";
        }
        if (types.count("TFBS") || types.count("TF_binding_site")) {
            annotations["regulatory:in_tfbs"] = "true";
        }
        if (types.count("open_chromatin") || types.count("open_chromatin_region")) {
            annotations["regulatory:in_open_chromatin"] = "true";
        }
        if (types.count("CTCF_binding_site")) {
            annotations["regulatory:in_ctcf"] = "true";
        }
    }

    std::vector<std::string> get_fields() const override {
        return {
            "regulatory:feature_type",
            "regulatory:feature_id",
            "regulatory:count",
            "regulatory:in_promoter",
            "regulatory:in_enhancer",
            "regulatory:in_tfbs",
            "regulatory:in_open_chromatin",
            "regulatory:in_ctcf"
        };
    }

    std::vector<Feature> query(
        const std::string& chrom,
        int start,
        int end
    ) const override {
        std::vector<Feature> result;
        if (!db_) return result;

        auto gff_features = db_->query(chrom, start, end);
        for (const auto* gf : gff_features) {
            Feature f;
            f.chrom = gf->seqid;
            f.start = gf->start;
            f.end = gf->end;
            f.type = gf->type;
            f.id = gf->id;
            f.strand = gf->strand;
            f.attributes = gf->attributes;
            result.push_back(f);
        }
        return result;
    }

    std::vector<Feature> query_by_type(
        const std::string& chrom,
        int start,
        int end,
        const std::string& feature_type
    ) const override {
        std::vector<Feature> result;
        if (!db_) return result;

        auto gff_features = db_->query_by_type(chrom, start, end, feature_type);
        for (const auto* gf : gff_features) {
            Feature f;
            f.chrom = gf->seqid;
            f.start = gf->start;
            f.end = gf->end;
            f.type = gf->type;
            f.id = gf->id;
            f.strand = gf->strand;
            f.attributes = gf->attributes;
            result.push_back(f);
        }
        return result;
    }

    std::string get_data_path() const override { return path_; }

    size_t memory_usage() const override {
        if (db_) {
            return db_->feature_count() * sizeof(GFF3Feature) * 2;  // Approximate
        }
        return 0;
    }

private:
    std::string path_;
    std::set<std::string> cell_types_;
    std::unique_ptr<GFF3Database> db_;
};

/**
 * Factory function
 */
std::shared_ptr<AnnotationSource> create_regulatory_source(
    const std::string& path,
    const std::set<std::string>& cell_types
) {
    return std::make_shared<RegulatorySource>(path, cell_types);
}

} // namespace vep
