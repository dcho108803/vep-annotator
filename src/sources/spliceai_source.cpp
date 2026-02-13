/**
 * SpliceAI Annotation Source
 *
 * Provides splice predictions from SpliceAI VCF files.
 * Supports separate SNV/indel files (matching Perl VEP plugin)
 * and a unified file for backward compatibility.
 *
 * Output fields match Perl VEP SpliceAI plugin:
 *   SpliceAI_pred, SpliceAI_pred_SYMBOL,
 *   SpliceAI_pred_DS_AG/AL/DG/DL, SpliceAI_pred_DP_AG/AL/DG/DL
 */

#include "annotation_source.hpp"
#include "file_parsers.hpp"
#include "vep_annotator.hpp"
#include <sstream>
#include <algorithm>
#include <cmath>

namespace vep {

class SpliceAISource : public VariantAnnotationSource {
public:
    SpliceAISource(const std::string& snv_path,
                   const std::string& indel_path,
                   const std::string& unified_path,
                   double cutoff = -1.0)
        : snv_path_(snv_path),
          indel_path_(indel_path),
          unified_path_(unified_path),
          cutoff_(cutoff) {}

    std::string name() const override { return "spliceai"; }
    std::string type() const override { return "splice"; }
    std::string description() const override {
        return "SpliceAI deep learning splice predictions";
    }

    bool is_ready() const override {
        return (snv_reader_ && snv_reader_->is_valid()) ||
               (indel_reader_ && indel_reader_->is_valid()) ||
               (unified_reader_ && unified_reader_->is_valid());
    }

    void initialize() override {
        std::lock_guard<std::recursive_mutex> lock(mutex_);
        if (snv_reader_ || indel_reader_ || unified_reader_) return;

        if (!snv_path_.empty()) {
            log(LogLevel::INFO, "Loading SpliceAI SNV from: " + snv_path_);
            snv_reader_ = std::make_unique<TabixTSVReader>(snv_path_, 0, 1);
            if (!snv_reader_->is_valid()) {
                log(LogLevel::ERROR, "Failed to open SpliceAI SNV VCF: " + snv_path_);
                snv_reader_.reset();
            }
        }

        if (!indel_path_.empty()) {
            log(LogLevel::INFO, "Loading SpliceAI indel from: " + indel_path_);
            indel_reader_ = std::make_unique<TabixTSVReader>(indel_path_, 0, 1);
            if (!indel_reader_->is_valid()) {
                log(LogLevel::ERROR, "Failed to open SpliceAI indel VCF: " + indel_path_);
                indel_reader_.reset();
            }
        }

        if (!unified_path_.empty()) {
            log(LogLevel::INFO, "Loading SpliceAI from: " + unified_path_);
            unified_reader_ = std::make_unique<TabixTSVReader>(unified_path_, 0, 1);
            if (!unified_reader_->is_valid()) {
                log(LogLevel::ERROR, "Failed to open SpliceAI VCF: " + unified_path_);
                unified_reader_.reset();
            }
        }

        log(LogLevel::INFO, "SpliceAI initialized");
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
        (void)transcript;

        // Query all active readers, return first match
        std::vector<TabixTSVReader*> readers;
        if (snv_reader_ && snv_reader_->is_valid()) readers.push_back(snv_reader_.get());
        if (indel_reader_ && indel_reader_->is_valid()) readers.push_back(indel_reader_.get());
        if (unified_reader_ && unified_reader_->is_valid()) readers.push_back(unified_reader_.get());

        for (auto* reader : readers) {
            if (query_reader(*reader, chrom, pos, ref, alt, annotations)) {
                return;
            }
        }
    }

    std::vector<std::string> get_fields() const override {
        std::vector<std::string> fields = {
            "SpliceAI_pred",           // Combined: SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
            "SpliceAI_pred_SYMBOL",    // Gene symbol
            "SpliceAI_pred_DS_AG",     // Delta score - acceptor gain
            "SpliceAI_pred_DS_AL",     // Delta score - acceptor loss
            "SpliceAI_pred_DS_DG",     // Delta score - donor gain
            "SpliceAI_pred_DS_DL",     // Delta score - donor loss
            "SpliceAI_pred_DP_AG",     // Delta position - acceptor gain
            "SpliceAI_pred_DP_AL",     // Delta position - acceptor loss
            "SpliceAI_pred_DP_DG",     // Delta position - donor gain
            "SpliceAI_pred_DP_DL",     // Delta position - donor loss
            "SpliceAI_pred_DS_max"     // Maximum delta score
        };
        if (cutoff_ >= 0) {
            fields.push_back("SpliceAI_cutoff");
        }
        return fields;
    }

    bool requires_allele_match() const override { return true; }

    std::unordered_map<std::string, std::string> query(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt
    ) const override {
        std::unordered_map<std::string, std::string> result;

        std::vector<TabixTSVReader*> readers;
        if (snv_reader_ && snv_reader_->is_valid()) readers.push_back(snv_reader_.get());
        if (indel_reader_ && indel_reader_->is_valid()) readers.push_back(indel_reader_.get());
        if (unified_reader_ && unified_reader_->is_valid()) readers.push_back(unified_reader_.get());

        for (auto* reader : readers) {
            auto records = reader->query(chrom, pos);
            for (const auto& record : records) {
                auto ref_it = record.find("REF");
                auto alt_it = record.find("ALT");
                if (ref_it == record.end() || alt_it == record.end()) continue;
                if (ref_it->second != ref) continue;

                int alt_index = find_alt_index(alt_it->second, alt);
                if (alt_index < 0) continue;

                auto info_it = record.find("INFO");
                if (info_it != record.end()) {
                    parse_spliceai_info(info_it->second, alt_index, result);
                    return result;
                }
            }
        }
        return result;
    }

    std::string get_data_path() const override {
        std::string paths;
        if (!snv_path_.empty()) paths += "snv=" + snv_path_;
        if (!indel_path_.empty()) {
            if (!paths.empty()) paths += ",";
            paths += "indel=" + indel_path_;
        }
        if (!unified_path_.empty()) {
            if (!paths.empty()) paths += ",";
            paths += unified_path_;
        }
        return paths;
    }

    bool is_thread_safe() const override { return true; }

    double get_cutoff() const { return cutoff_; }

private:
    std::string snv_path_;
    std::string indel_path_;
    std::string unified_path_;
    double cutoff_;

    mutable std::unique_ptr<TabixTSVReader> snv_reader_;
    mutable std::unique_ptr<TabixTSVReader> indel_reader_;
    mutable std::unique_ptr<TabixTSVReader> unified_reader_;

    static int find_alt_index(const std::string& alts, const std::string& target) {
        std::istringstream iss(alts);
        std::string a;
        int idx = 0;
        while (std::getline(iss, a, ',')) {
            if (a == target) return idx;
            idx++;
        }
        return -1;
    }

    bool query_reader(TabixTSVReader& reader,
                      const std::string& chrom, int pos,
                      const std::string& ref, const std::string& alt,
                      std::unordered_map<std::string, std::string>& annotations) const {
        auto records = reader.query(chrom, pos);
        for (const auto& record : records) {
            auto ref_it = record.find("REF");
            auto alt_it = record.find("ALT");
            if (ref_it == record.end() || alt_it == record.end()) continue;
            if (ref_it->second != ref) continue;

            int alt_index = find_alt_index(alt_it->second, alt);
            if (alt_index < 0) continue;

            auto info_it = record.find("INFO");
            if (info_it == record.end()) continue;

            parse_spliceai_info(info_it->second, alt_index, annotations);
            return true;
        }
        return false;
    }

    void parse_spliceai_info(const std::string& info, int alt_index,
                             std::unordered_map<std::string, std::string>& annotations) const {
        // SpliceAI INFO format: SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
        // Multiple alleles separated by ','

        size_t start = info.find("SpliceAI=");
        if (start == std::string::npos) return;

        std::string spliceai_value = info.substr(start + 9);
        size_t end = spliceai_value.find(';');
        if (end != std::string::npos) {
            spliceai_value = spliceai_value.substr(0, end);
        }

        // Split by comma for multiple alleles
        std::vector<std::string> allele_values;
        {
            std::istringstream iss(spliceai_value);
            std::string v;
            while (std::getline(iss, v, ',')) {
                allele_values.push_back(v);
            }
        }

        if (alt_index >= static_cast<int>(allele_values.size())) return;

        // Parse pipe-delimited fields for this allele
        std::string value = allele_values[alt_index];
        std::vector<std::string> fields;
        {
            std::istringstream iss(value);
            std::string f;
            while (std::getline(iss, f, '|')) {
                fields.push_back(f);
            }
        }

        if (fields.size() < 10) return;

        // fields: [0]=ALLELE, [1]=SYMBOL, [2]=DS_AG, [3]=DS_AL, [4]=DS_DG, [5]=DS_DL,
        //         [6]=DP_AG, [7]=DP_AL, [8]=DP_DG, [9]=DP_DL

        // Combined field (matches Perl VEP SpliceAI_pred format): SYMBOL|DS_AG|...|DP_DL
        annotations["SpliceAI_pred"] = fields[1] + "|" + fields[2] + "|" + fields[3] + "|" +
                                       fields[4] + "|" + fields[5] + "|" + fields[6] + "|" +
                                       fields[7] + "|" + fields[8] + "|" + fields[9];

        // Individual fields
        if (fields[1] != ".") annotations["SpliceAI_pred_SYMBOL"] = fields[1];

        double ds_ag = 0, ds_al = 0, ds_dg = 0, ds_dl = 0;

        if (fields[2] != ".") {
            annotations["SpliceAI_pred_DS_AG"] = fields[2];
            try { ds_ag = std::stod(fields[2]); } catch (...) {}
        }
        if (fields[3] != ".") {
            annotations["SpliceAI_pred_DS_AL"] = fields[3];
            try { ds_al = std::stod(fields[3]); } catch (...) {}
        }
        if (fields[4] != ".") {
            annotations["SpliceAI_pred_DS_DG"] = fields[4];
            try { ds_dg = std::stod(fields[4]); } catch (...) {}
        }
        if (fields[5] != ".") {
            annotations["SpliceAI_pred_DS_DL"] = fields[5];
            try { ds_dl = std::stod(fields[5]); } catch (...) {}
        }
        if (fields[6] != ".") annotations["SpliceAI_pred_DP_AG"] = fields[6];
        if (fields[7] != ".") annotations["SpliceAI_pred_DP_AL"] = fields[7];
        if (fields[8] != ".") annotations["SpliceAI_pred_DP_DG"] = fields[8];
        if (fields[9] != ".") annotations["SpliceAI_pred_DP_DL"] = fields[9];

        // Max delta score - preserve original string formatting
        double max_ds = std::max({ds_ag, ds_al, ds_dg, ds_dl});
        if (max_ds > 0) {
            std::string max_str;
            if (max_ds == ds_ag && fields[2] != ".") max_str = fields[2];
            else if (max_ds == ds_al && fields[3] != ".") max_str = fields[3];
            else if (max_ds == ds_dg && fields[4] != ".") max_str = fields[4];
            else if (max_ds == ds_dl && fields[5] != ".") max_str = fields[5];
            if (!max_str.empty()) {
                annotations["SpliceAI_pred_DS_max"] = max_str;
            }
        }

        // Cutoff PASS/FAIL
        if (cutoff_ >= 0) {
            annotations["SpliceAI_cutoff"] = (max_ds >= cutoff_) ? "PASS" : "FAIL";
        }
    }
};

// Legacy single-file factory (backward compatibility)
std::shared_ptr<AnnotationSource> create_spliceai_source(const std::string& path) {
    return std::make_shared<SpliceAISource>("", "", path);
}

// New factory with separate SNV/indel files and cutoff
std::shared_ptr<AnnotationSource> create_spliceai_source(
    const std::string& snv_path,
    const std::string& indel_path,
    const std::string& unified_path,
    double cutoff) {
    return std::make_shared<SpliceAISource>(snv_path, indel_path, unified_path, cutoff);
}

} // namespace vep
