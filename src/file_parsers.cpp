/**
 * File Format Parsers - Implementation
 */

#include "file_parsers.hpp"
#include "annotation_source.hpp"
#include "vep_annotator.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <sys/stat.h>

#ifdef HAVE_HTSLIB
#include <htslib/tbx.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#endif

#ifdef HAVE_BIGWIG
#include <bigWig.h>
#endif

#include <zlib.h>

namespace vep {

// ============================================================================
// Utility Functions
// ============================================================================

std::vector<std::string> split_line(const std::string& line, char delim) {
    std::vector<std::string> result;
    size_t start = 0;
    size_t pos = line.find(delim);
    while (pos != std::string::npos) {
        result.emplace_back(line, start, pos - start);
        start = pos + 1;
        pos = line.find(delim, start);
    }
    result.emplace_back(line, start);
    return result;
}

std::unordered_map<std::string, std::string> parse_gff3_attributes(const std::string& attrs) {
    std::unordered_map<std::string, std::string> result;

    auto fields = split_line(attrs, ';');
    for (const auto& field : fields) {
        size_t eq = field.find('=');
        if (eq != std::string::npos) {
            std::string key = field.substr(0, eq);
            std::string value = url_decode(field.substr(eq + 1));
            result[key] = value;
        }
    }

    return result;
}

std::string url_decode(const std::string& str) {
    std::string result;
    result.reserve(str.size());

    auto hex_val = [](char c) -> int {
        if (c >= '0' && c <= '9') return c - '0';
        if (c >= 'a' && c <= 'f') return c - 'a' + 10;
        if (c >= 'A' && c <= 'F') return c - 'A' + 10;
        return -1;
    };

    for (size_t i = 0; i < str.size(); ++i) {
        if (str[i] == '%' && i + 2 < str.size()) {
            int hi = hex_val(str[i + 1]);
            int lo = hex_val(str[i + 2]);
            if (hi >= 0 && lo >= 0) {
                result += static_cast<char>(hi * 16 + lo);
                i += 2;
                continue;
            }
        }
        result += str[i];
    }

    return result;
}

bool file_exists(const std::string& path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

std::string get_extension(const std::string& path) {
    std::string ext;
    size_t dot = path.rfind('.');

    if (dot != std::string::npos) {
        ext = path.substr(dot);

        // Handle .gz compression
        if (ext == ".gz" && dot > 0) {
            size_t dot2 = path.rfind('.', dot - 1);
            if (dot2 != std::string::npos) {
                ext = path.substr(dot2);
            }
        }
    }

    return ext;
}

// ============================================================================
// TabixTSVReader Implementation
// ============================================================================

#ifdef HAVE_HTSLIB

struct TabixTSVReader::Impl {
    htsFile* fp = nullptr;
    tbx_t* tbx = nullptr;
    std::vector<std::string> columns;
    std::map<std::string, int> column_index;
    int chrom_col;
    int pos_col;
    bool valid = false;

    ~Impl() {
        if (tbx) tbx_destroy(tbx);
        if (fp) hts_close(fp);
    }
};

TabixTSVReader::TabixTSVReader(
    const std::string& path,
    int chrom_col,
    int pos_col,
    const std::vector<std::string>& columns
) : pimpl_(std::make_unique<Impl>()), path_(path) {

    pimpl_->chrom_col = chrom_col;
    pimpl_->pos_col = pos_col;

    // Open file
    pimpl_->fp = hts_open(path.c_str(), "r");
    if (!pimpl_->fp) {
        log(LogLevel::ERROR, "Cannot open TSV file: " + path);
        return;
    }

    // Load tabix index
    pimpl_->tbx = tbx_index_load(path.c_str());
    if (!pimpl_->tbx) {
        log(LogLevel::ERROR, "Cannot load tabix index for: " + path);
        return;
    }

    // Read header to get column names
    kstring_t str = {0, 0, nullptr};
    while (hts_getline(pimpl_->fp, KS_SEP_LINE, &str) >= 0) {
        if (str.l == 0) continue;
        if (str.s[0] == '#') {
            // Parse header line
            std::string header(str.s);
            if (header[0] == '#') header = header.substr(1);

            pimpl_->columns = split_line(header, '\t');
            for (size_t i = 0; i < pimpl_->columns.size(); ++i) {
                pimpl_->column_index[pimpl_->columns[i]] = i;
            }
        } else {
            break;  // End of header
        }
    }
    free(str.s);

    // Use provided columns if specified
    if (!columns.empty()) {
        pimpl_->columns = columns;
        // Re-index
        pimpl_->column_index.clear();
        for (size_t i = 0; i < columns.size(); ++i) {
            pimpl_->column_index[columns[i]] = i;
        }
    }

    pimpl_->valid = true;
    log(LogLevel::INFO, "Opened tabix TSV: " + path + " (" +
        std::to_string(pimpl_->columns.size()) + " columns)");
}

TabixTSVReader::~TabixTSVReader() = default;

std::vector<std::map<std::string, std::string>> TabixTSVReader::query(
    const std::string& chrom,
    int pos
) {
    return query_range(chrom, pos, pos);
}

std::vector<std::map<std::string, std::string>> TabixTSVReader::query_range(
    const std::string& chrom,
    int start,
    int end
) {
    std::vector<std::map<std::string, std::string>> results;

    if (!pimpl_->valid) return results;

    // Try different chromosome formats
    std::vector<std::string> chrom_variants;
    if (chrom.substr(0, 3) == "chr") {
        chrom_variants.push_back(chrom);
        chrom_variants.push_back(chrom.substr(3));
    } else {
        chrom_variants.push_back("chr" + chrom);
        chrom_variants.push_back(chrom);
    }

    for (const auto& try_chrom : chrom_variants) {
        std::string region = try_chrom + ":" + std::to_string(start) + "-" + std::to_string(end);

        hts_itr_t* itr = tbx_itr_querys(pimpl_->tbx, region.c_str());
        if (!itr) continue;

        kstring_t str = {0, 0, nullptr};

        while (tbx_itr_next(pimpl_->fp, pimpl_->tbx, itr, &str) >= 0) {
            auto fields = split_line(std::string(str.s, str.l), '\t');

            std::map<std::string, std::string> row;
            for (size_t i = 0; i < fields.size() && i < pimpl_->columns.size(); ++i) {
                row[pimpl_->columns[i]] = fields[i];
            }

            results.push_back(std::move(row));
        }

        free(str.s);
        tbx_itr_destroy(itr);

        if (!results.empty()) break;
    }

    return results;
}

std::vector<std::string> TabixTSVReader::get_columns() const {
    return pimpl_->columns;
}

bool TabixTSVReader::is_valid() const {
    return pimpl_->valid;
}

#else  // No HTSLIB

struct TabixTSVReader::Impl {
    bool valid = false;
};

TabixTSVReader::TabixTSVReader(
    const std::string& path,
    int chrom_col,
    int pos_col,
    const std::vector<std::string>& columns
) : pimpl_(std::make_unique<Impl>()), path_(path) {
    log(LogLevel::WARNING, "TabixTSVReader requires htslib. Build with -DHAVE_HTSLIB");
}

TabixTSVReader::~TabixTSVReader() = default;

std::vector<std::map<std::string, std::string>> TabixTSVReader::query(
    const std::string&, int) {
    return {};
}

std::vector<std::map<std::string, std::string>> TabixTSVReader::query_range(
    const std::string&, int, int) {
    return {};
}

std::vector<std::string> TabixTSVReader::get_columns() const {
    return {};
}

bool TabixTSVReader::is_valid() const {
    return false;
}

#endif  // HAVE_HTSLIB

// ============================================================================
// BigWigReader Implementation
// ============================================================================

#ifdef HAVE_BIGWIG

struct BigWigReader::Impl {
    bigWigFile_t* bw = nullptr;
    std::vector<std::string> chroms;
    std::set<std::string> chrom_set;
    bool valid = false;

    ~Impl() {
        if (bw) bwClose(bw);
    }
};

BigWigReader::BigWigReader(const std::string& path)
    : pimpl_(std::make_unique<Impl>()), path_(path) {

    // Initialize libBigWig
    if (bwInit(1 << 17) != 0) {
        log(LogLevel::ERROR, "Failed to initialize libBigWig");
        return;
    }

    // Open file
    pimpl_->bw = bwOpen(const_cast<char*>(path.c_str()), nullptr, "r");
    if (!pimpl_->bw) {
        log(LogLevel::ERROR, "Cannot open bigWig file: " + path);
        return;
    }

    // Get chromosome list
    if (pimpl_->bw->cl) {
        for (int64_t i = 0; i < pimpl_->bw->cl->nKeys; ++i) {
            std::string chrom(pimpl_->bw->cl->chrom[i]);
            pimpl_->chroms.push_back(chrom);
            pimpl_->chrom_set.insert(chrom);
            pimpl_->chrom_set.insert(normalize_chrom(chrom));
        }
    }

    pimpl_->valid = true;
    log(LogLevel::INFO, "Opened bigWig: " + path + " (" +
        std::to_string(pimpl_->chroms.size()) + " chromosomes)");
}

BigWigReader::~BigWigReader() {
    bwCleanup();
}

std::optional<double> BigWigReader::get_value(
    const std::string& chrom,
    int pos
) const {
    if (!pimpl_->valid) return std::nullopt;

    // Try different chromosome formats
    std::vector<std::string> try_chroms = {chrom, "chr" + chrom};
    if (chrom.substr(0, 3) == "chr") {
        try_chroms.push_back(chrom.substr(3));
    }

    for (const auto& try_chrom : try_chroms) {
        if (pimpl_->chrom_set.count(try_chrom) == 0) continue;

        double* vals = bwStats(pimpl_->bw, const_cast<char*>(try_chrom.c_str()),
                               pos - 1, pos, 1, mean);
        if (vals) {
            double result = vals[0];
            free(vals);
            if (!std::isnan(result)) {
                return result;
            }
        }
    }

    return std::nullopt;
}

std::vector<double> BigWigReader::get_values(
    const std::string& chrom,
    int start,
    int end
) const {
    std::vector<double> results;

    if (!pimpl_->valid) return results;

    std::vector<std::string> try_chroms = {chrom, "chr" + chrom};
    if (chrom.substr(0, 3) == "chr") {
        try_chroms.push_back(chrom.substr(3));
    }

    for (const auto& try_chrom : try_chroms) {
        if (pimpl_->chrom_set.count(try_chrom) == 0) continue;

        bwOverlappingIntervals_t* intervals = bwGetValues(
            pimpl_->bw, const_cast<char*>(try_chrom.c_str()),
            start - 1, end, 0
        );

        if (intervals) {
            for (uint32_t i = 0; i < intervals->l; ++i) {
                results.push_back(intervals->value[i]);
            }
            bwDestroyOverlappingIntervals(intervals);
            break;
        }
    }

    return results;
}

std::optional<double> BigWigReader::get_mean(
    const std::string& chrom,
    int start,
    int end
) const {
    if (!pimpl_->valid) return std::nullopt;

    std::vector<std::string> try_chroms = {chrom, "chr" + chrom};
    if (chrom.substr(0, 3) == "chr") {
        try_chroms.push_back(chrom.substr(3));
    }

    for (const auto& try_chrom : try_chroms) {
        if (pimpl_->chrom_set.count(try_chrom) == 0) continue;

        double* vals = bwStats(pimpl_->bw, const_cast<char*>(try_chrom.c_str()),
                               start - 1, end, 1, mean);
        if (vals) {
            double result = vals[0];
            free(vals);
            if (!std::isnan(result)) {
                return result;
            }
        }
    }

    return std::nullopt;
}

std::optional<double> BigWigReader::get_max(
    const std::string& chrom,
    int start,
    int end
) const {
    if (!pimpl_->valid) return std::nullopt;

    std::vector<std::string> try_chroms = {chrom, "chr" + chrom};
    if (chrom.substr(0, 3) == "chr") {
        try_chroms.push_back(chrom.substr(3));
    }

    for (const auto& try_chrom : try_chroms) {
        if (pimpl_->chrom_set.count(try_chrom) == 0) continue;

        double* vals = bwStats(pimpl_->bw, const_cast<char*>(try_chrom.c_str()),
                               start - 1, end, 1, max);
        if (vals) {
            double result = vals[0];
            free(vals);
            if (!std::isnan(result)) {
                return result;
            }
        }
    }

    return std::nullopt;
}

bool BigWigReader::has_chromosome(const std::string& chrom) const {
    return pimpl_->chrom_set.count(chrom) > 0 ||
           pimpl_->chrom_set.count("chr" + chrom) > 0 ||
           pimpl_->chrom_set.count(normalize_chrom(chrom)) > 0;
}

std::vector<std::string> BigWigReader::get_chromosomes() const {
    return pimpl_->chroms;
}

bool BigWigReader::is_valid() const {
    return pimpl_->valid;
}

#else  // No BIGWIG

struct BigWigReader::Impl {
    bool valid = false;
};

BigWigReader::BigWigReader(const std::string& path)
    : pimpl_(std::make_unique<Impl>()), path_(path) {
    log(LogLevel::WARNING, "BigWigReader requires libBigWig. Build with -DHAVE_BIGWIG");
}

BigWigReader::~BigWigReader() = default;

std::optional<double> BigWigReader::get_value(const std::string&, int) const {
    return std::nullopt;
}

std::vector<double> BigWigReader::get_values(const std::string&, int, int) const {
    return {};
}

std::optional<double> BigWigReader::get_mean(const std::string&, int, int) const {
    return std::nullopt;
}

std::optional<double> BigWigReader::get_max(const std::string&, int, int) const {
    return std::nullopt;
}

bool BigWigReader::has_chromosome(const std::string&) const {
    return false;
}

std::vector<std::string> BigWigReader::get_chromosomes() const {
    return {};
}

bool BigWigReader::is_valid() const {
    return false;
}

#endif  // HAVE_BIGWIG

// ============================================================================
// GFF3Database Implementation
// ============================================================================

struct GFF3Database::Impl {
    std::vector<GFF3Feature> features;
    std::map<std::string, IntervalTree<size_t>> chrom_trees;
    std::set<std::string> feature_types;
    bool built = false;
};

GFF3Database::GFF3Database(
    const std::string& path,
    const std::set<std::string>& feature_types
) : pimpl_(std::make_unique<Impl>()), path_(path) {

    log(LogLevel::INFO, "Loading GFF3 file: " + path);

    bool is_gzipped = (path.size() >= 3 && path.substr(path.size() - 3) == ".gz");

    std::function<bool(std::string&)> read_line;
    gzFile gz = nullptr;
    std::ifstream file;

    if (is_gzipped) {
        gz = gzopen(path.c_str(), "rb");
        if (!gz) {
            log(LogLevel::ERROR, "Cannot open GFF3 file: " + path);
            return;
        }

        char buffer[65536];
        read_line = [&gz, &buffer](std::string& line) -> bool {
            if (gzgets(gz, buffer, sizeof(buffer)) == nullptr) return false;
            line = buffer;
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
                line.pop_back();
            }
            return true;
        };
    } else {
        file.open(path);
        if (!file.is_open()) {
            log(LogLevel::ERROR, "Cannot open GFF3 file: " + path);
            return;
        }

        read_line = [&file](std::string& line) -> bool {
            if (!std::getline(file, line)) return false;
            while (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            return true;
        };
    }

    std::string line;
    size_t count = 0;

    while (read_line(line)) {
        if (line.empty() || line[0] == '#') continue;

        auto fields = split_line(line, '\t');
        if (fields.size() < 9) continue;

        std::string type = fields[2];

        // Filter by feature types if specified
        if (!feature_types.empty() && feature_types.count(type) == 0) {
            continue;
        }

        GFF3Feature feat;
        feat.seqid = normalize_chrom(fields[0]);
        feat.source = fields[1];
        feat.type = type;
        try { feat.start = std::stoi(fields[3]); }
        catch (...) { continue; }
        try { feat.end = std::stoi(fields[4]); }
        catch (...) { continue; }

        if (fields[5] != ".") {
            try { feat.score = std::stod(fields[5]); }
            catch (...) {}
        }

        feat.strand = fields[6].empty() ? '.' : fields[6][0];

        if (fields[7] != ".") {
            try { feat.phase = std::stoi(fields[7]); }
            catch (...) {}
        }

        // Parse attributes
        feat.attributes = parse_gff3_attributes(fields[8]);
        feat.id = feat.attributes.count("ID") ? feat.attributes["ID"] : "";
        feat.name = feat.attributes.count("Name") ? feat.attributes["Name"] : "";
        feat.parent = feat.attributes.count("Parent") ? feat.attributes["Parent"] : "";

        pimpl_->feature_types.insert(type);

        size_t idx = pimpl_->features.size();
        pimpl_->features.push_back(std::move(feat));

        // Add to interval tree
        pimpl_->chrom_trees[pimpl_->features[idx].seqid].insert(
            pimpl_->features[idx].start,
            pimpl_->features[idx].end,
            idx
        );

        count++;
    }

    if (gz) gzclose(gz);

    // Build all interval trees
    for (auto& [chrom, tree] : pimpl_->chrom_trees) {
        tree.build();
    }

    pimpl_->built = true;

    log(LogLevel::INFO, "Loaded " + std::to_string(count) + " GFF3 features from " + path);
}

GFF3Database::~GFF3Database() = default;

std::vector<const GFF3Feature*> GFF3Database::query(
    const std::string& chrom,
    int pos
) const {
    return query(chrom, pos, pos);
}

std::vector<const GFF3Feature*> GFF3Database::query(
    const std::string& chrom,
    int start,
    int end
) const {
    std::vector<const GFF3Feature*> results;
    std::string norm_chrom = normalize_chrom(chrom);

    auto it = pimpl_->chrom_trees.find(norm_chrom);
    if (it == pimpl_->chrom_trees.end()) {
        return results;
    }

    auto indices = it->second.query(start, end);
    for (size_t idx : indices) {
        results.push_back(&pimpl_->features[idx]);
    }

    return results;
}

std::vector<const GFF3Feature*> GFF3Database::query_by_type(
    const std::string& chrom,
    int start,
    int end,
    const std::string& type
) const {
    auto all = query(chrom, start, end);
    std::vector<const GFF3Feature*> results;

    for (const auto* feat : all) {
        if (feat->type == type) {
            results.push_back(feat);
        }
    }

    return results;
}

std::set<std::string> GFF3Database::get_feature_types() const {
    return pimpl_->feature_types;
}

size_t GFF3Database::feature_count() const {
    return pimpl_->features.size();
}

// ============================================================================
// IntervalTree Implementation
// ============================================================================

template<typename T>
void IntervalTree<T>::insert(int start, int end, T data) {
    intervals_.push_back({start, end, std::move(data)});
    built_ = false;
}

template<typename T>
void IntervalTree<T>::build() {
    if (intervals_.empty()) {
        built_ = true;
        return;
    }

    std::vector<size_t> indices(intervals_.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }

    root_ = build_tree(indices);
    built_ = true;
}

template<typename T>
std::unique_ptr<typename IntervalTree<T>::Node> IntervalTree<T>::build_tree(
    std::vector<size_t>& indices,
    int depth
) {
    if (indices.empty()) return nullptr;

    // Find center point
    int min_start = INT_MAX, max_end = INT_MIN;
    for (size_t idx : indices) {
        min_start = std::min(min_start, intervals_[idx].start);
        max_end = std::max(max_end, intervals_[idx].end);
    }
    int center = (min_start + max_end) / 2;

    auto node = std::make_unique<Node>();
    node->center = center;

    std::vector<size_t> left_indices, right_indices;

    for (size_t idx : indices) {
        const auto& interval = intervals_[idx];

        if (interval.end < center) {
            left_indices.push_back(idx);
        } else if (interval.start > center) {
            right_indices.push_back(idx);
        } else {
            node->overlapping.push_back(idx);
        }
    }

    // Sort overlapping by start
    std::sort(node->overlapping.begin(), node->overlapping.end(),
              [this](size_t a, size_t b) {
                  return intervals_[a].start < intervals_[b].start;
              });

    if (depth < 20) {  // Prevent deep recursion
        node->left = build_tree(left_indices, depth + 1);
        node->right = build_tree(right_indices, depth + 1);
    } else {
        // Flatten remaining intervals into current node to avoid silently dropping them
        for (size_t idx : left_indices) node->overlapping.push_back(idx);
        for (size_t idx : right_indices) node->overlapping.push_back(idx);
    }

    return node;
}

template<typename T>
std::vector<T> IntervalTree<T>::query(int point) const {
    return query(point, point);
}

template<typename T>
std::vector<T> IntervalTree<T>::query(int start, int end) const {
    std::vector<T> results;

    if (!built_ || !root_) return results;

    query_node(root_.get(), start, end, results);
    return results;
}

template<typename T>
void IntervalTree<T>::query_node(
    const Node* node,
    int start,
    int end,
    std::vector<T>& results
) const {
    if (!node) return;

    // Check overlapping intervals at this node
    for (size_t idx : node->overlapping) {
        const auto& interval = intervals_[idx];
        if (interval.start <= end && interval.end >= start) {
            results.push_back(interval.data);
        }
    }

    // Recurse to children
    if (start < node->center && node->left) {
        query_node(node->left.get(), start, end, results);
    }
    if (end > node->center && node->right) {
        query_node(node->right.get(), start, end, results);
    }
}

template<typename T>
void IntervalTree<T>::clear() {
    intervals_.clear();
    root_.reset();
    built_ = false;
}

// Explicit template instantiations
template class IntervalTree<size_t>;
template class IntervalTree<int>;
template class IntervalTree<std::string>;

// ============================================================================
// AnnotationSourceManager Implementation
// ============================================================================

void AnnotationSourceManager::add_source(std::shared_ptr<AnnotationSource> source) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    sources_.push_back(std::move(source));
}

std::vector<std::shared_ptr<AnnotationSource>> AnnotationSourceManager::get_sources() const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    return sources_;
}

std::shared_ptr<AnnotationSource> AnnotationSourceManager::get_source(
    const std::string& name
) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    for (const auto& source : sources_) {
        if (source->name() == name) {
            return source;
        }
    }
    return nullptr;
}

void AnnotationSourceManager::set_enabled(const std::string& name, bool enabled) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    if (enabled) {
        disabled_.erase(name);
    } else {
        disabled_.insert(name);
    }
}

bool AnnotationSourceManager::is_enabled(const std::string& name) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    return disabled_.count(name) == 0;
}

void AnnotationSourceManager::initialize_all() {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    for (auto& source : sources_) {
        if (!source->is_ready() && disabled_.count(source->name()) == 0) {
            source->initialize();
        }
    }
}

void AnnotationSourceManager::annotate_all(
    const std::string& chrom,
    int pos,
    const std::string& ref,
    const std::string& alt,
    const Transcript* transcript,
    std::unordered_map<std::string, std::string>& annotations
) {
    std::shared_lock<std::shared_mutex> lock(mutex_);

    for (auto& source : sources_) {
        if (disabled_.count(source->name()) != 0) continue;

        if (!source->is_ready()) {
            // Initialize on first use
            lock.unlock();
            source->initialize();
            lock.lock();
        }

        try {
            source->annotate(chrom, pos, ref, alt, transcript, annotations);
        } catch (const std::exception& e) {
            log(LogLevel::WARNING, "Annotation source '" + source->name() +
                "' failed: " + e.what());
        }
    }
}

std::vector<std::pair<std::string, std::string>> AnnotationSourceManager::get_all_fields() const {
    std::vector<std::pair<std::string, std::string>> result;
    std::shared_lock<std::shared_mutex> lock(mutex_);

    for (const auto& source : sources_) {
        for (const auto& field : source->get_fields()) {
            result.emplace_back(field, source->type());
        }
    }

    return result;
}

size_t AnnotationSourceManager::total_memory_usage() const {
    size_t total = 0;
    std::shared_lock<std::shared_mutex> lock(mutex_);

    for (const auto& source : sources_) {
        total += source->memory_usage();
    }

    return total;
}

std::string AnnotationSourceManager::get_stats() const {
    std::ostringstream oss;
    std::shared_lock<std::shared_mutex> lock(mutex_);

    oss << "Annotation Sources (" << sources_.size() << "):\n";

    for (const auto& source : sources_) {
        oss << "  " << source->name() << " [" << source->type() << "]";

        if (disabled_.count(source->name())) {
            oss << " (disabled)";
        } else if (source->is_ready()) {
            oss << " (ready)";
        } else {
            oss << " (not loaded)";
        }

        size_t mem = source->memory_usage();
        if (mem > 0) {
            oss << " - " << (mem / 1024 / 1024) << " MB";
        }

        oss << "\n";
    }

    return oss.str();
}

// ============================================================================
// ScoreAnnotationSource Default Implementation
// ============================================================================

std::optional<double> ScoreAnnotationSource::get_aggregated_score(
    const std::string& chrom,
    int start,
    int end,
    Aggregation method
) const {
    auto scores = get_scores(chrom, start, end);

    if (scores.empty()) return std::nullopt;

    // Remove NaN values
    scores.erase(
        std::remove_if(scores.begin(), scores.end(),
                       [](double d) { return std::isnan(d); }),
        scores.end()
    );

    if (scores.empty()) return std::nullopt;

    switch (method) {
        case Aggregation::MEAN: {
            double sum = 0;
            for (double s : scores) sum += s;
            return sum / scores.size();
        }
        case Aggregation::MAX:
            return *std::max_element(scores.begin(), scores.end());
        case Aggregation::MIN:
            return *std::min_element(scores.begin(), scores.end());
        case Aggregation::FIRST:
            return scores.front();
        case Aggregation::LAST:
            return scores.back();
    }

    return std::nullopt;
}

} // namespace vep
