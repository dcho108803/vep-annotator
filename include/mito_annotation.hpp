/**
 * Mitochondrial annotation helpers
 *
 * Heteroplasmy: mitochondrial variants are usually present in only a
 * fraction of the cell's mtDNA copies. Callers report that fraction as the
 * per-sample FORMAT AF field (DRAGEN mito mode, GATK Mutect2 --mitochondria);
 * when AF is absent it can be derived from allelic depths (AD).
 */

#ifndef MITO_ANNOTATION_HPP
#define MITO_ANNOTATION_HPP

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include <cstdio>
#include <cctype>

#include "file_parsers.hpp"

namespace vep {

/**
 * True if the chromosome names the mitochondrial contig (chrM/M/chrMT/MT).
 */
inline bool is_mito_chrom(const std::string& chrom) {
    return normalize_chrom(chrom) == "MT";
}

/**
 * Compute the heteroplasmy fraction for one alt allele of a VCF record.
 *
 * @param sample_columns  FORMAT column followed by tab-separated sample
 *                        columns, exactly as read from the VCF line
 *                        (e.g. "GT:AF:AD\t0/1:0.85:15,85"). Values are taken
 *                        from the first sample; mitochondrial calling is
 *                        effectively always single-sample.
 * @param allele_index    1-based alt allele number within the record
 *                        (1 for the first ALT), matching ALLELE_NUM.
 * @return The fraction as a string ("0.85"), or "" when unavailable.
 *
 * FORMAT AF (Number=A, one value per alt) is preferred and passed through
 * verbatim; AD (Number=R: ref,alt1,...) is the fallback, reported as
 * AD[allele_index] / sum(AD) with 4 significant digits.
 */
inline std::string compute_heteroplasmy(const std::string& sample_columns,
                                        int allele_index) {
    if (sample_columns.empty() || allele_index < 1) return "";

    size_t tab = sample_columns.find('\t');
    if (tab == std::string::npos) return "";
    std::string format = sample_columns.substr(0, tab);
    std::string first_sample = sample_columns.substr(tab + 1);
    size_t next_tab = first_sample.find('\t');
    if (next_tab != std::string::npos) {
        first_sample = first_sample.substr(0, next_tab);
    }

    auto split_colon = [](const std::string& s) {
        std::vector<std::string> out;
        std::istringstream ss(s);
        std::string tok;
        while (std::getline(ss, tok, ':')) out.push_back(tok);
        return out;
    };
    auto split_comma = [](const std::string& s) {
        std::vector<std::string> out;
        std::istringstream ss(s);
        std::string tok;
        while (std::getline(ss, tok, ',')) out.push_back(tok);
        return out;
    };

    std::vector<std::string> keys = split_colon(format);
    std::vector<std::string> values = split_colon(first_sample);

    auto value_of = [&](const std::string& key) -> std::string {
        for (size_t i = 0; i < keys.size() && i < values.size(); ++i) {
            if (keys[i] == key) return values[i];
        }
        return "";
    };

    // Preferred: caller-reported allele fraction (Number=A)
    std::string af = value_of("AF");
    if (!af.empty() && af != ".") {
        std::vector<std::string> entries = split_comma(af);
        size_t idx = static_cast<size_t>(allele_index - 1);
        if (idx < entries.size() && !entries[idx].empty() && entries[idx] != ".") {
            // Validate numeric before passing through verbatim
            try {
                (void)std::stod(entries[idx]);
                return entries[idx];
            } catch (...) {
                return "";
            }
        }
        return "";
    }

    // Fallback: derive from allelic depths (Number=R: ref,alt1,alt2,...)
    std::string ad = value_of("AD");
    if (!ad.empty() && ad != ".") {
        std::vector<std::string> entries = split_comma(ad);
        size_t idx = static_cast<size_t>(allele_index);  // ref is entry 0
        if (idx < entries.size()) {
            double total = 0.0, alt_depth = 0.0;
            try {
                for (size_t i = 0; i < entries.size(); ++i) {
                    if (entries[i].empty() || entries[i] == ".") continue;
                    double d = std::stod(entries[i]);
                    if (d < 0) return "";
                    total += d;
                    if (i == idx) alt_depth = d;
                }
            } catch (...) {
                return "";
            }
            if (total > 0.0) {
                char buf[32];
                std::snprintf(buf, sizeof(buf), "%.4g", alt_depth / total);
                return buf;
            }
        }
    }

    return "";
}

/**
 * Mitochondrial disease variant database (--mitomap)
 *
 * Position+allele keyed lookup for MITOMAP / HmtVar style variant lists,
 * emitting the associated disease and curation status for matching MT
 * variants. Accepts:
 * - TSV with a header row naming pos/ref/alt (or a combined allele column
 *   like "A3243G" or "m.3243A>G") plus disease and status columns
 * - headerless TSV in the order POS REF ALT DISEASE [STATUS]
 * - VCF keyed on POS/REF/ALT with Disease/CLNDN and Status/CLNSIG INFO keys
 */
class MitoDiseaseDB {
public:
    struct Entry {
        std::string disease;
        std::string status;
    };

    bool load(const std::string& path) {
        std::ifstream file(path);
        if (!file.is_open()) return false;

        bool is_vcf = path.find(".vcf") != std::string::npos;
        // Column indices for TSV mode (headerless defaults)
        int pos_col = 0, ref_col = 1, alt_col = 2, disease_col = 3, status_col = 4;
        int allele_col = -1;
        bool first_line = true;
        bool have_header = false;

        std::string line;
        while (std::getline(file, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            if (first_line && line.compare(0, 16, "##fileformat=VCF") == 0) is_vcf = true;
            first_line = false;
            if (line[0] == '#' && line.compare(0, 2, "##") == 0) continue;

            std::vector<std::string> f = split_line(line, '\t');

            if (is_vcf) {
                if (line[0] == '#') continue;
                load_vcf_fields(f);
                continue;
            }

            // Header row: map named columns (case-insensitive contains)
            std::string maybe_header = line[0] == '#' ? line.substr(1) : line;
            if (!have_header && !is_numeric(f.empty() ? "" : strip_hash(f[0]))) {
                std::vector<std::string> cols = split_line(maybe_header, '\t');
                pos_col = ref_col = alt_col = disease_col = status_col = -1;
                for (size_t i = 0; i < cols.size(); ++i) {
                    std::string c = lower(cols[i]);
                    int idx = static_cast<int>(i);
                    if (pos_col < 0 && (c == "pos" || c == "position")) pos_col = idx;
                    else if (ref_col < 0 && c.compare(0, 3, "ref") == 0) ref_col = idx;
                    else if (alt_col < 0 && (c.compare(0, 3, "alt") == 0 || c == "variant")) alt_col = idx;
                    else if (allele_col < 0 && c == "allele") allele_col = idx;
                    else if (disease_col < 0 && (c.find("disease") != std::string::npos ||
                                                 c.find("phenotype") != std::string::npos ||
                                                 c == "clndn")) disease_col = idx;
                    else if (status_col < 0 && (c.find("status") != std::string::npos ||
                                                c.find("significance") != std::string::npos ||
                                                c == "clnsig")) status_col = idx;
                }
                have_header = true;
                continue;
            }
            if (line[0] == '#') continue;

            int pos = 0;
            std::string ref, alt;
            if (allele_col >= 0 && allele_col < static_cast<int>(f.size()) &&
                parse_allele_notation(f[allele_col], pos, ref, alt)) {
                // combined allele column handled
            } else if (pos_col >= 0 && ref_col >= 0 && alt_col >= 0 &&
                       std::max(pos_col, std::max(ref_col, alt_col)) < static_cast<int>(f.size())) {
                try {
                    pos = std::stoi(f[pos_col]);
                } catch (const std::exception&) {
                    continue;
                }
                ref = f[ref_col];
                alt = f[alt_col];
            } else {
                continue;
            }
            if (pos <= 0 || ref.empty() || alt.empty()) continue;

            Entry e;
            if (disease_col >= 0 && disease_col < static_cast<int>(f.size())) {
                e.disease = sanitize(f[disease_col]);
            }
            if (status_col >= 0 && status_col < static_cast<int>(f.size())) {
                e.status = sanitize(f[status_col]);
            }
            if (e.disease.empty() && e.status.empty()) continue;
            entries_[key(pos, ref, alt)] = e;
        }
        loaded_ = true;
        return true;
    }

    bool is_loaded() const { return loaded_; }
    size_t size() const { return entries_.size(); }

    const Entry* lookup(int pos, const std::string& ref, const std::string& alt) const {
        auto it = entries_.find(key(pos, ref, alt));
        return it == entries_.end() ? nullptr : &it->second;
    }

private:
    static std::string lower(const std::string& s) {
        std::string out = s;
        for (char& c : out) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return out;
    }
    static std::string upper(const std::string& s) {
        std::string out = s;
        for (char& c : out) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        return out;
    }
    static std::string strip_hash(const std::string& s) {
        return (!s.empty() && s[0] == '#') ? s.substr(1) : s;
    }
    static bool is_numeric(const std::string& s) {
        if (s.empty()) return false;
        for (char c : s) {
            if (!std::isdigit(static_cast<unsigned char>(c))) return false;
        }
        return true;
    }
    // Output fields are ;/= delimited downstream; keep values single-token
    static std::string sanitize(const std::string& s) {
        std::string out;
        for (char c : s) {
            if (c == ' ' || c == '\t') out += '_';
            else if (c == ';') out += ',';
            else if (c == '=') out += ':';
            else out += c;
        }
        return out;
    }
    static std::string key(int pos, const std::string& ref, const std::string& alt) {
        return std::to_string(pos) + ":" + upper(ref) + ":" + upper(alt);
    }

    /**
     * Parse combined allele spellings: "A3243G", "m.3243A>G", "3243A>G".
     */
    static bool parse_allele_notation(const std::string& raw, int& pos,
                                      std::string& ref, std::string& alt) {
        std::string s = upper(raw);
        if (s.compare(0, 2, "M.") == 0) s = s.substr(2);

        size_t gt = s.find('>');
        if (gt != std::string::npos) {
            // 3243A>G: digits, ref base(s), '>', alt base(s)
            size_t i = 0;
            while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) ++i;
            if (i == 0 || i >= gt) return false;
            try {
                pos = std::stoi(s.substr(0, i));
            } catch (const std::exception&) {
                return false;
            }
            ref = s.substr(i, gt - i);
            alt = s.substr(gt + 1);
            return !ref.empty() && !alt.empty();
        }

        // A3243G: ref base, digits, alt base
        if (s.size() >= 3 && std::isalpha(static_cast<unsigned char>(s.front())) &&
            std::isalpha(static_cast<unsigned char>(s.back()))) {
            std::string digits = s.substr(1, s.size() - 2);
            if (!is_numeric(digits)) return false;
            try {
                pos = std::stoi(digits);
            } catch (const std::exception&) {
                return false;
            }
            ref = s.substr(0, 1);
            alt = s.substr(s.size() - 1);
            return true;
        }
        return false;
    }

    void load_vcf_fields(const std::vector<std::string>& f) {
        if (f.size() < 8) return;
        int pos = 0;
        try {
            pos = std::stoi(f[1]);
        } catch (const std::exception&) {
            return;
        }
        std::map<std::string, std::string> info;
        std::istringstream info_ss(f[7]);
        std::string tok;
        while (std::getline(info_ss, tok, ';')) {
            size_t eq = tok.find('=');
            if (eq != std::string::npos) info[lower(tok.substr(0, eq))] = tok.substr(eq + 1);
        }
        Entry e;
        for (const char* k : {"disease", "clndn", "phenotype"}) {
            auto it = info.find(k);
            if (it != info.end()) { e.disease = sanitize(it->second); break; }
        }
        for (const char* k : {"status", "clnsig", "clinicalsignificance"}) {
            auto it = info.find(k);
            if (it != info.end()) { e.status = sanitize(it->second); break; }
        }
        if (e.disease.empty() && e.status.empty()) return;
        // One entry per alt allele
        std::istringstream alts(f[4]);
        std::string alt;
        while (std::getline(alts, alt, ',')) {
            if (!alt.empty()) entries_[key(pos, f[3], alt)] = e;
        }
    }

    std::map<std::string, Entry> entries_;
    bool loaded_ = false;
};

} // namespace vep

#endif // MITO_ANNOTATION_HPP
