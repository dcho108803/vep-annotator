/**
 * Protein Domain Annotation Sources
 *
 * Provides protein domain annotations from Pfam and InterPro databases.
 */

#include "annotation_source.hpp"
#include "file_parsers.hpp"
#include "vep_annotator.hpp"
#include <sstream>
#include <fstream>
#include <algorithm>

namespace vep {

/**
 * Domain entry structure
 */
struct DomainEntry {
    std::string transcript_id;
    std::string protein_id;
    int aa_start;
    int aa_end;
    std::string domain_id;
    std::string domain_name;
    std::string domain_desc;
    double evalue = 0;
};

/**
 * Protein Domain Annotation Source
 *
 * Matches variants to overlapping protein domains based on amino acid position.
 */
class DomainSource : public AnnotationSource {
public:
    DomainSource(const std::string& path,
                 const std::string& source_name,
                 const std::string& desc)
        : path_(path), source_name_(source_name), description_(desc) {}

    std::string name() const override { return source_name_; }
    std::string type() const override { return "domain"; }
    std::string description() const override { return description_; }

    bool is_ready() const override { return loaded_; }

    void initialize() override {
        std::lock_guard<std::mutex> lock(mutex_);
        if (loaded_) return;

        log(LogLevel::INFO, "Loading " + source_name_ + " domains from: " + path_);

        std::ifstream file(path_);
        if (!file.is_open()) {
            log(LogLevel::ERROR, "Failed to open domain file: " + path_);
            return;
        }

        std::string line;
        size_t count = 0;

        // Skip header
        if (std::getline(file, line)) {
            // Check if header line
            if (line[0] != '#' && line.find("transcript") == std::string::npos) {
                // Not a header, parse it
                parse_domain_line(line);
                count++;
            }
        }

        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            parse_domain_line(line);
            count++;
        }

        loaded_ = true;
        log(LogLevel::INFO, source_name_ + " loaded " + std::to_string(count) +
                           " domain entries for " + std::to_string(domains_.size()) + " transcripts");
    }

    void annotate(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript* transcript,
        std::map<std::string, std::string>& annotations
    ) override {
        ensure_initialized();
        if (!transcript) return;

        (void)chrom;
        (void)pos;
        (void)ref;
        (void)alt;

        // Look up domains for this transcript
        auto it = domains_.find(transcript->id);
        if (it == domains_.end()) {
            // Try gene name as fallback
            it = domains_.find(transcript->gene_name);
        }

        if (it == domains_.end()) return;

        // Find overlapping domains based on amino acid position
        // Note: This requires the caller to have computed protein_position
        // For now, we return all domains for the transcript

        const auto& domain_list = it->second;
        if (domain_list.empty()) return;

        // Collect domain info
        std::vector<std::string> domain_ids;
        std::vector<std::string> domain_names;

        for (const auto& d : domain_list) {
            domain_ids.push_back(d.domain_id);
            if (!d.domain_name.empty()) {
                domain_names.push_back(d.domain_name);
            }
        }

        // Format output
        if (!domain_ids.empty()) {
            std::ostringstream oss;
            for (size_t i = 0; i < domain_ids.size() && i < 5; ++i) {
                if (i > 0) oss << ",";
                oss << domain_ids[i];
            }
            if (domain_ids.size() > 5) oss << ",...";
            annotations[source_name_ + ":domain_id"] = oss.str();
        }

        if (!domain_names.empty()) {
            std::ostringstream oss;
            for (size_t i = 0; i < domain_names.size() && i < 3; ++i) {
                if (i > 0) oss << ",";
                oss << domain_names[i];
            }
            annotations[source_name_ + ":domain_name"] = oss.str();
        }

        annotations[source_name_ + ":domain_count"] = std::to_string(domain_list.size());
    }

    std::vector<std::string> get_fields() const override {
        return {
            source_name_ + ":domain_id",
            source_name_ + ":domain_name",
            source_name_ + ":domain_count"
        };
    }

    std::string get_data_path() const override { return path_; }

    size_t memory_usage() const override {
        size_t size = 0;
        for (const auto& [key, domains] : domains_) {
            size += key.size() + domains.size() * sizeof(DomainEntry);
        }
        return size;
    }

protected:
    std::string path_;
    std::string source_name_;
    std::string description_;
    bool loaded_ = false;
    std::map<std::string, std::vector<DomainEntry>> domains_;

    virtual void parse_domain_line(const std::string& line) {
        // Default TSV format: transcript_id, aa_start, aa_end, domain_id, domain_name, evalue
        auto fields = split_line(line, '\t');
        if (fields.size() < 4) return;

        DomainEntry entry;
        entry.transcript_id = fields[0];

        try {
            entry.aa_start = std::stoi(fields[1]);
            entry.aa_end = std::stoi(fields[2]);
        } catch (...) {
            return;
        }

        entry.domain_id = fields[3];
        if (fields.size() > 4) entry.domain_name = fields[4];
        if (fields.size() > 5) {
            try {
                entry.evalue = std::stod(fields[5]);
            } catch (...) {}
        }

        domains_[entry.transcript_id].push_back(entry);
    }
};

/**
 * Pfam Domain Source
 */
class PfamSource : public DomainSource {
public:
    explicit PfamSource(const std::string& path)
        : DomainSource(path, "pfam", "Pfam protein domain annotations") {}

protected:
    void parse_domain_line(const std::string& line) override {
        // Pfam format varies, try to parse flexibly
        auto fields = split_line(line, '\t');
        if (fields.size() < 4) return;

        DomainEntry entry;

        // Common Pfam TSV format: seq_id, alignment_start, alignment_end, pfam_acc, pfam_name, ...
        entry.transcript_id = fields[0];

        try {
            entry.aa_start = std::stoi(fields[1]);
            entry.aa_end = std::stoi(fields[2]);
        } catch (...) {
            return;
        }

        entry.domain_id = fields[3];
        if (fields.size() > 4) entry.domain_name = fields[4];
        if (fields.size() > 5) entry.domain_desc = fields[5];

        // Look for E-value
        for (size_t i = 6; i < fields.size(); ++i) {
            if (fields[i].find('e') != std::string::npos ||
                fields[i].find('E') != std::string::npos) {
                try {
                    entry.evalue = std::stod(fields[i]);
                    break;
                } catch (...) {}
            }
        }

        domains_[entry.transcript_id].push_back(entry);
    }
};

/**
 * InterPro Domain Source
 */
class InterProSource : public DomainSource {
public:
    explicit InterProSource(const std::string& path)
        : DomainSource(path, "interpro", "InterPro protein domain annotations") {}

protected:
    void parse_domain_line(const std::string& line) override {
        // InterPro TSV format: protein_id, MD5, length, analysis, signature_acc, signature_desc, start, end, evalue, status, date, interpro_acc, interpro_desc, GO_terms, pathways
        auto fields = split_line(line, '\t');
        if (fields.size() < 8) return;

        DomainEntry entry;

        entry.protein_id = fields[0];
        entry.transcript_id = fields[0];  // Use protein_id as key

        try {
            entry.aa_start = std::stoi(fields[6]);
            entry.aa_end = std::stoi(fields[7]);
        } catch (...) {
            return;
        }

        // Use InterPro accession if available, otherwise signature
        if (fields.size() > 11 && !fields[11].empty() && fields[11] != "-") {
            entry.domain_id = fields[11];  // InterPro accession
            if (fields.size() > 12) entry.domain_name = fields[12];
        } else {
            entry.domain_id = fields[4];   // Signature accession
            if (fields.size() > 5) entry.domain_name = fields[5];
        }

        // E-value
        if (fields.size() > 8 && fields[8] != "-") {
            try {
                entry.evalue = std::stod(fields[8]);
            } catch (...) {}
        }

        domains_[entry.transcript_id].push_back(entry);
    }
};

/**
 * Factory functions
 */
std::shared_ptr<AnnotationSource> create_pfam_source(const std::string& path) {
    return std::make_shared<PfamSource>(path);
}

std::shared_ptr<AnnotationSource> create_interpro_source(const std::string& path) {
    return std::make_shared<InterProSource>(path);
}

} // namespace vep
