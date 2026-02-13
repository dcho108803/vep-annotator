/**
 * SpliceAI Annotation Source
 *
 * Provides splice predictions from SpliceAI VCF file.
 * SpliceAI predicts splice-altering variants using deep learning.
 */

#include "annotation_source.hpp"
#include "file_parsers.hpp"
#include "vep_annotator.hpp"
#include <sstream>
#include <algorithm>
#include <cmath>

namespace vep {

/**
 * SpliceAI Annotation Source
 *
 * Queries SpliceAI VCF for delta scores at each position.
 * Scores include: DS_AG (acceptor gain), DS_AL (acceptor loss),
 *                 DS_DG (donor gain), DS_DL (donor loss)
 */
class SpliceAISource : public VariantAnnotationSource {
public:
    explicit SpliceAISource(const std::string& path)
        : path_(path) {}

    std::string name() const override { return "spliceai"; }
    std::string type() const override { return "splice"; }
    std::string description() const override {
        return "SpliceAI deep learning splice predictions";
    }

    bool is_ready() const override { return reader_ != nullptr && reader_->is_valid(); }

    void initialize() override {
        std::lock_guard<std::recursive_mutex> lock(mutex_);
        if (reader_) return;

        log(LogLevel::INFO, "Loading SpliceAI from: " + path_);

        // SpliceAI VCF format: CHROM, POS, ID, REF, ALT, ...
        // INFO field contains SpliceAI scores
        reader_ = std::make_unique<TabixTSVReader>(path_, 0, 1);

        if (!reader_->is_valid()) {
            log(LogLevel::ERROR, "Failed to open SpliceAI VCF: " + path_);
            reader_.reset();
            return;
        }

        log(LogLevel::INFO, "SpliceAI loaded successfully");
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
        if (!reader_) return;

        (void)transcript;  // Not used

        auto records = reader_->query(chrom, pos);

        for (const auto& record : records) {
            // Check REF/ALT match
            auto ref_it = record.find("REF");
            auto alt_it = record.find("ALT");

            if (ref_it == record.end() || alt_it == record.end()) continue;
            if (ref_it->second != ref) continue;

            // Check if ALT matches (may have multiple alts)
            std::string alts = alt_it->second;
            int alt_index = -1;
            {
                std::istringstream iss(alts);
                std::string a;
                int idx = 0;
                while (std::getline(iss, a, ',')) {
                    if (a == alt) {
                        alt_index = idx;
                        break;
                    }
                    idx++;
                }
            }

            if (alt_index < 0) continue;

            // Parse SpliceAI INFO field
            auto info_it = record.find("INFO");
            if (info_it == record.end()) continue;

            parse_spliceai_info(info_it->second, alt_index, annotations);
            return;
        }
    }

    std::vector<std::string> get_fields() const override {
        return {
            "spliceai:DS_AG",   // Delta score - acceptor gain
            "spliceai:DS_AL",   // Delta score - acceptor loss
            "spliceai:DS_DG",   // Delta score - donor gain
            "spliceai:DS_DL",   // Delta score - donor loss
            "spliceai:DP_AG",   // Delta position - acceptor gain
            "spliceai:DP_AL",   // Delta position - acceptor loss
            "spliceai:DP_DG",   // Delta position - donor gain
            "spliceai:DP_DL",   // Delta position - donor loss
            "spliceai:max_DS",  // Maximum delta score
            "spliceai:SYMBOL"   // Gene symbol
        };
    }

    bool requires_allele_match() const override { return true; }

    std::unordered_map<std::string, std::string> query(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt
    ) const override {
        std::unordered_map<std::string, std::string> result;

        if (!reader_) return result;

        auto records = reader_->query(chrom, pos);

        for (const auto& record : records) {
            auto ref_it = record.find("REF");
            auto alt_it = record.find("ALT");

            if (ref_it == record.end() || alt_it == record.end()) continue;
            if (ref_it->second != ref) continue;

            // Find exact ALT match by splitting comma-delimited alts
            int alt_index = -1;
            {
                std::istringstream iss(alt_it->second);
                std::string a;
                int idx = 0;
                while (std::getline(iss, a, ',')) {
                    if (a == alt) {
                        alt_index = idx;
                        break;
                    }
                    idx++;
                }
            }
            if (alt_index < 0) continue;

            auto info_it = record.find("INFO");
            if (info_it != record.end()) {
                std::unordered_map<std::string, std::string> tmp;
                parse_spliceai_info(info_it->second, alt_index, tmp);
                for (auto& kv : tmp) {
                    result[kv.first] = std::move(kv.second);
                }
                break;
            }
        }

        return result;
    }

    std::string get_data_path() const override { return path_; }

    /**
     * Interpret SpliceAI max delta score
     * @return "high_impact" (>=0.8), "moderate_impact" (>=0.5), "low_impact" (>=0.2), or empty
     */
    static std::string interpret_max_score(double max_ds) {
        if (max_ds >= 0.8) return "high_impact";
        if (max_ds >= 0.5) return "moderate_impact";
        if (max_ds >= 0.2) return "low_impact";
        return "";
    }

private:
    std::string path_;
    std::unique_ptr<TabixTSVReader> reader_;

    void parse_spliceai_info(const std::string& info, int alt_index,
                             std::unordered_map<std::string, std::string>& annotations) const {
        // SpliceAI INFO format: SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
        // May have multiple alleles separated by ','

        size_t pos = info.find("SpliceAI=");
        if (pos == std::string::npos) return;

        std::string spliceai_value = info.substr(pos + 9);

        // Find end of value (next ';' or end of string)
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

        if (alt_index >= static_cast<int>(allele_values.size())) {
            return;  // No matching allele data available
        }

        // Parse the value for this allele
        std::string value = allele_values[alt_index];

        // Split by '|'
        std::vector<std::string> fields;
        {
            std::istringstream iss(value);
            std::string f;
            while (std::getline(iss, f, '|')) {
                fields.push_back(f);
            }
        }

        if (fields.size() >= 10) {
            // fields: ALLELE, SYMBOL, DS_AG, DS_AL, DS_DG, DS_DL, DP_AG, DP_AL, DP_DG, DP_DL
            if (fields[1] != ".") annotations["spliceai:SYMBOL"] = fields[1];

            double ds_ag = 0, ds_al = 0, ds_dg = 0, ds_dl = 0;

            if (fields[2] != ".") {
                annotations["spliceai:DS_AG"] = fields[2];
                try { ds_ag = std::stod(fields[2]); } catch (...) {}
            }
            if (fields[3] != ".") {
                annotations["spliceai:DS_AL"] = fields[3];
                try { ds_al = std::stod(fields[3]); } catch (...) {}
            }
            if (fields[4] != ".") {
                annotations["spliceai:DS_DG"] = fields[4];
                try { ds_dg = std::stod(fields[4]); } catch (...) {}
            }
            if (fields[5] != ".") {
                annotations["spliceai:DS_DL"] = fields[5];
                try { ds_dl = std::stod(fields[5]); } catch (...) {}
            }
            if (fields[6] != ".") annotations["spliceai:DP_AG"] = fields[6];
            if (fields[7] != ".") annotations["spliceai:DP_AL"] = fields[7];
            if (fields[8] != ".") annotations["spliceai:DP_DG"] = fields[8];
            if (fields[9] != ".") annotations["spliceai:DP_DL"] = fields[9];

            // Calculate max delta score - use original string to preserve formatting
            double max_ds = std::max({ds_ag, ds_al, ds_dg, ds_dl});
            if (max_ds > 0) {
                // Find which field string corresponds to max value
                std::string max_str;
                if (max_ds == ds_ag && fields[2] != ".") max_str = fields[2];
                else if (max_ds == ds_al && fields[3] != ".") max_str = fields[3];
                else if (max_ds == ds_dg && fields[4] != ".") max_str = fields[4];
                else if (max_ds == ds_dl && fields[5] != ".") max_str = fields[5];
                if (!max_str.empty()) {
                    annotations["spliceai:max_DS"] = max_str;
                }
            }
        }
    }
};

/**
 * Factory function to create SpliceAI source
 */
std::shared_ptr<AnnotationSource> create_spliceai_source(const std::string& path) {
    return std::make_shared<SpliceAISource>(path);
}

} // namespace vep
