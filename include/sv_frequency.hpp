/**
 * Population SV frequency database (--sv-frequency)
 *
 * Matches input structural variants against a population SV callset —
 * a gnomAD-SV style VCF or a DGV-style TSV/BED — by SV class and overlap,
 * so output can carry the population allele frequency and identifier of
 * the co-located known SV. Matching uses the --overlaps mode (default
 * reciprocal) and --overlap-cutoff fraction.
 */

#ifndef SV_FREQUENCY_HPP
#define SV_FREQUENCY_HPP

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "structural_variant.hpp"
#include "file_parsers.hpp"

namespace vep {

class SVFrequencyDB {
public:
    struct Record {
        int start = 0;              // 1-based inclusive
        int end = 0;
        SVType type = SVType::UNKNOWN;
        double af = -1.0;           // population allele frequency (-1 unknown)
        std::string id;
    };

    struct Match {
        Record record;
        double overlap = 0.0;       // reciprocal overlap fraction [0,1]
    };

    /**
     * Load a population SV file. Auto-detects format:
     * - VCF (".vcf" extension or ##fileformat=VCF header): SVTYPE/END from
     *   INFO (or symbolic ALT), AF from INFO AF (first value).
     * - TSV (DGV-style, 1-based): CHROM\tSTART\tEND\tTYPE[\tAF[\tID]]
     * - BED (".bed" extension, 0-based start): same columns as TSV.
     * Plain text only (decompress .gz files first).
     */
    bool load(const std::string& path) {
        std::ifstream file(path);
        if (!file.is_open()) return false;

        bool is_vcf = path.find(".vcf") != std::string::npos;
        bool is_bed = !is_vcf && path.size() > 4 &&
                      path.compare(path.size() - 4, 4, ".bed") == 0;

        std::string line;
        bool first_line = true;
        while (std::getline(file, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            if (first_line && line.compare(0, 16, "##fileformat=VCF") == 0) {
                is_vcf = true;
            }
            first_line = false;
            if (line[0] == '#') continue;

            if (is_vcf) {
                load_vcf_line(line);
            } else {
                load_tsv_line(line, is_bed);
            }
        }

        for (auto& kv : by_chrom_) {
            std::sort(kv.second.begin(), kv.second.end(),
                      [](const Record& a, const Record& b) { return a.start < b.start; });
        }
        loaded_ = true;
        return true;
    }

    bool is_loaded() const { return loaded_; }

    size_t size() const {
        size_t n = 0;
        for (const auto& kv : by_chrom_) n += kv.second.size();
        return n;
    }

    /**
     * Best match for an SV: among same-class records satisfying the overlap
     * mode and cutoff, the one with the highest reciprocal overlap.
     * @param min_overlap fraction [0,1]; applies to RECIPROCAL mode
     */
    bool best_match(const StructuralVariant& sv, Match& out,
                    double min_overlap = 0.7,
                    OverlapType mode = OverlapType::RECIPROCAL) const {
        auto it = by_chrom_.find(normalize_chrom(sv.chromosome));
        if (it == by_chrom_.end()) return false;

        const std::vector<Record>& records = it->second;
        // Records are sorted by start; everything starting after sv.end is out
        auto upper = std::upper_bound(
            records.begin(), records.end(), sv.end,
            [](int end, const Record& r) { return end < r.start; });

        bool found = false;
        for (auto rit = records.begin(); rit != upper; ++rit) {
            if (rit->end < sv.start) continue;
            if (!class_compatible(sv, rit->type)) continue;

            bool mode_ok = false;
            switch (mode) {
                case OverlapType::ANY:
                    mode_ok = true;
                    break;
                case OverlapType::WITHIN:
                    mode_ok = rit->start <= sv.start && rit->end >= sv.end;
                    break;
                case OverlapType::SURROUNDING:
                    mode_ok = sv.start <= rit->start && sv.end >= rit->end;
                    break;
                case OverlapType::EXACT:
                    mode_ok = sv.start == rit->start && sv.end == rit->end;
                    break;
                case OverlapType::RECIPROCAL:
                    mode_ok = calculate_overlap_percentage(
                                  sv.start, sv.end, rit->start, rit->end,
                                  /*reciprocal=*/true) >= min_overlap;
                    break;
            }
            if (!mode_ok) continue;

            double recip = calculate_overlap_percentage(
                sv.start, sv.end, rit->start, rit->end, /*reciprocal=*/true);
            if (!found || recip > out.overlap) {
                out.record = *rit;
                out.overlap = recip;
                found = true;
            }
        }
        return found;
    }

private:
    /**
     * SV class compatibility for matching: TDUP matches the DUP class, a CNV
     * with known copy number matches its loss/gain class, and CNV records in
     * the database (e.g. DGV gain+loss entries) match either class.
     */
    static bool class_compatible(const StructuralVariant& sv, SVType db_type) {
        auto cls = [](SVType t) { return t == SVType::TDUP ? SVType::DUP : t; };
        SVType q = cls(sv.sv_type);
        SVType d = cls(db_type);
        if (q == SVType::CNV && sv.copy_number >= 0) {
            if (sv.copy_number < 1.5) q = SVType::DEL;
            else if (sv.copy_number > 2.0) q = SVType::DUP;
        }
        if (q == d) return true;
        bool q_cnv_like = (q == SVType::DEL || q == SVType::DUP || q == SVType::CNV);
        bool d_cnv_like = (d == SVType::DEL || d == SVType::DUP || d == SVType::CNV);
        return (d == SVType::CNV && q_cnv_like) || (q == SVType::CNV && d_cnv_like);
    }

    void load_vcf_line(const std::string& line) {
        std::vector<std::string> f = split_line(line, '\t');
        if (f.size() < 8) return;

        std::map<std::string, std::string> info;
        std::istringstream info_ss(f[7]);
        std::string tok;
        while (std::getline(info_ss, tok, ';')) {
            size_t eq = tok.find('=');
            if (eq != std::string::npos) info[tok.substr(0, eq)] = tok.substr(eq + 1);
            else info[tok] = "";
        }

        int pos = 0;
        try {
            pos = std::stoi(f[1]);
        } catch (const std::exception&) {
            return;
        }

        // Reuse the SV record parser for SVTYPE/END/symbolic-ALT handling
        StructuralVariant sv = parse_sv_from_vcf(f[0], pos, f[3], f[4], info);
        if (!sv.is_sv()) return;

        Record rec;
        rec.start = sv.start;
        rec.end = sv.end;
        rec.type = sv.sv_type;
        rec.id = (f[2] == ".") ? "" : f[2];
        auto af_it = info.find("AF");
        if (af_it != info.end()) {
            // Multi-allelic AF is comma-separated; take the first value
            std::string af = af_it->second.substr(0, af_it->second.find(','));
            try {
                rec.af = std::stod(af);
            } catch (const std::exception&) {}
        }
        by_chrom_[normalize_chrom(f[0])].push_back(rec);
    }

    void load_tsv_line(const std::string& line, bool zero_based) {
        std::vector<std::string> f = split_line(line, '\t');
        if (f.size() < 4) return;
        // Tolerate a headerless header row ("chr\tstart\t...")
        Record rec;
        try {
            rec.start = std::stoi(f[1]) + (zero_based ? 1 : 0);
            rec.end = std::stoi(f[2]);
        } catch (const std::exception&) {
            return;
        }
        rec.type = parse_sv_type(f[3]);
        if (rec.type == SVType::UNKNOWN || rec.type == SVType::SNV) return;
        if (f.size() >= 5 && !f[4].empty() && f[4] != ".") {
            try {
                rec.af = std::stod(f[4]);
            } catch (const std::exception&) {}
        }
        if (f.size() >= 6) rec.id = f[5];
        by_chrom_[normalize_chrom(f[0])].push_back(rec);
    }

    std::map<std::string, std::vector<Record>> by_chrom_;
    bool loaded_ = false;
};

} // namespace vep

#endif // SV_FREQUENCY_HPP
