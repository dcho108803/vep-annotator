/**
 * Gene Constraint Scores
 *
 * Supports pLI (probability of loss-of-function intolerance),
 * LOEUF (loss-of-function observed/expected upper bound fraction),
 * and other gene-level constraint metrics.
 */

#ifndef GENE_CONSTRAINT_HPP
#define GENE_CONSTRAINT_HPP

#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>

namespace vep {

/**
 * Gene constraint data
 */
struct GeneConstraint {
    std::string gene_symbol;
    std::string gene_id;        // Ensembl gene ID

    // gnomAD constraint metrics
    double pLI = -1.0;          // Probability of LoF intolerance
    double pRec = -1.0;         // Probability of recessive
    double pNull = -1.0;        // Probability of tolerant

    // LOEUF (LoF observed/expected upper confidence fraction)
    double oe_lof = -1.0;       // Observed/expected ratio for LoF
    double oe_lof_upper = -1.0; // Upper bound (LOEUF)
    double oe_lof_lower = -1.0; // Lower bound

    // Missense constraint
    double oe_mis = -1.0;       // O/E for missense
    double oe_mis_upper = -1.0;
    double oe_mis_lower = -1.0;
    double mis_z = -1.0;        // Missense Z-score

    // Synonymous (control)
    double oe_syn = -1.0;
    double oe_syn_upper = -1.0;
    double oe_syn_lower = -1.0;
    double syn_z = -1.0;

    // Expected/observed counts
    int exp_lof = 0;
    int obs_lof = 0;
    int exp_mis = 0;
    int obs_mis = 0;
    int exp_syn = 0;
    int obs_syn = 0;

    // CDS metrics
    int cds_length = 0;
    int num_coding_exons = 0;

    // Haploinsufficiency/triplosensitivity
    double hi_score = -1.0;     // Haploinsufficiency score
    double ts_score = -1.0;     // Triplosensitivity score

    // Other metrics
    double s_het = -1.0;        // Selection coefficient for heterozygotes

    bool has_data() const {
        return pLI >= 0 || oe_lof_upper >= 0;
    }

    bool is_constrained() const {
        // Common thresholds: pLI > 0.9 or LOEUF < 0.35
        if (pLI >= 0.9) return true;
        if (oe_lof_upper >= 0 && oe_lof_upper < 0.35) return true;
        return false;
    }

    std::string get_constraint_level() const {
        if (pLI >= 0.9 || (oe_lof_upper >= 0 && oe_lof_upper < 0.35)) {
            return "highly_constrained";
        }
        if (pLI >= 0.5 || (oe_lof_upper >= 0 && oe_lof_upper < 0.6)) {
            return "moderately_constrained";
        }
        if (pLI >= 0 || oe_lof_upper >= 0) {
            return "tolerant";
        }
        return "unknown";
    }
};

/**
 * Gene Constraint Database
 *
 * Loads constraint data from gnomAD constraint files or similar sources.
 */
class GeneConstraintDB {
public:
    GeneConstraintDB() : loaded_(false) {}

    /**
     * Load gnomAD constraint file
     * Expected format (TSV with header):
     * gene    transcript    pLI    oe_lof_upper    ...
     */
    bool load_gnomad_constraint(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            return false;
        }

        std::string line;
        std::map<std::string, int> col_map;

        // Read header
        if (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string col;
            int idx = 0;
            while (std::getline(iss, col, '\t')) {
                col_map[col] = idx++;
            }
        }

        // Validate required columns
        bool has_gene = col_map.count("gene") > 0 || col_map.count("gene_symbol") > 0;
        if (!has_gene) {
            return false;
        }

        // Read data rows
        while (std::getline(file, line)) {
            if (line.empty()) continue;

            std::vector<std::string> fields;
            std::istringstream iss(line);
            std::string field;
            while (std::getline(iss, field, '\t')) {
                fields.push_back(field);
            }

            GeneConstraint constraint;

            // Get gene symbol
            int gene_col = -1;
            if (col_map.count("gene") > 0) {
                gene_col = col_map["gene"];
            } else if (col_map.count("gene_symbol") > 0) {
                gene_col = col_map["gene_symbol"];
            }

            if (gene_col >= 0 && gene_col < static_cast<int>(fields.size())) {
                constraint.gene_symbol = fields[gene_col];
            } else {
                continue;
            }

            // Get gene ID
            if (col_map.count("gene_id") > 0) {
                int idx = col_map["gene_id"];
                if (idx < static_cast<int>(fields.size())) {
                    constraint.gene_id = fields[idx];
                }
            }

            // Parse constraint metrics
            auto parse_double = [&](const std::string& col_name) -> double {
                if (col_map.count(col_name) == 0) return -1.0;
                int idx = col_map[col_name];
                if (idx >= static_cast<int>(fields.size())) return -1.0;
                const std::string& val = fields[idx];
                if (val.empty() || val == "NA" || val == "." || val == "NaN") return -1.0;
                try {
                    return std::stod(val);
                } catch (...) {
                    return -1.0;
                }
            };

            auto parse_int = [&](const std::string& col_name) -> int {
                if (col_map.count(col_name) == 0) return 0;
                int idx = col_map[col_name];
                if (idx >= static_cast<int>(fields.size())) return 0;
                const std::string& val = fields[idx];
                if (val.empty() || val == "NA" || val == ".") return 0;
                try {
                    return std::stoi(val);
                } catch (...) {
                    return 0;
                }
            };

            // pLI and related
            constraint.pLI = parse_double("pLI");
            constraint.pRec = parse_double("pRec");
            constraint.pNull = parse_double("pNull");

            // LOEUF
            constraint.oe_lof = parse_double("oe_lof");
            constraint.oe_lof_upper = parse_double("oe_lof_upper");
            constraint.oe_lof_lower = parse_double("oe_lof_lower");

            // Missense
            constraint.oe_mis = parse_double("oe_mis");
            constraint.oe_mis_upper = parse_double("oe_mis_upper");
            constraint.oe_mis_lower = parse_double("oe_mis_lower");
            constraint.mis_z = parse_double("mis_z");

            // Synonymous
            constraint.oe_syn = parse_double("oe_syn");
            constraint.oe_syn_upper = parse_double("oe_syn_upper");
            constraint.oe_syn_lower = parse_double("oe_syn_lower");
            constraint.syn_z = parse_double("syn_z");

            // Counts
            constraint.exp_lof = parse_int("exp_lof");
            constraint.obs_lof = parse_int("obs_lof");
            constraint.exp_mis = parse_int("exp_mis");
            constraint.obs_mis = parse_int("obs_mis");
            constraint.exp_syn = parse_int("exp_syn");
            constraint.obs_syn = parse_int("obs_syn");

            // CDS metrics
            constraint.cds_length = parse_int("cds_length");
            constraint.num_coding_exons = parse_int("num_coding_exons");

            // Haploinsufficiency
            constraint.hi_score = parse_double("hi_score");
            constraint.ts_score = parse_double("ts_score");

            // Selection coefficient
            constraint.s_het = parse_double("s_het");

            // Store by gene symbol
            if (!constraint.gene_symbol.empty()) {
                gene_data_[constraint.gene_symbol] = constraint;
            }

            // Also store by gene ID
            if (!constraint.gene_id.empty()) {
                gene_id_data_[constraint.gene_id] = constraint;
            }
        }

        loaded_ = true;
        return true;
    }

    /**
     * Load pLI scores from simple format
     * Format: GENE\tpLI\n
     */
    bool load_pli_scores(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            return false;
        }

        std::string line;
        bool has_header = false;

        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;

            // Skip header
            if (!has_header && (line.find("gene") != std::string::npos ||
                                line.find("pLI") != std::string::npos)) {
                has_header = true;
                continue;
            }

            std::istringstream iss(line);
            std::string gene;
            double pli;

            if (!(iss >> gene >> pli)) continue;

            GeneConstraint constraint;
            constraint.gene_symbol = gene;
            constraint.pLI = pli;

            gene_data_[gene] = constraint;
        }

        loaded_ = true;
        return true;
    }

    /**
     * Load LOEUF scores
     * Format: GENE\tLOEUF\n
     */
    bool load_loeuf_scores(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file.is_open()) {
            return false;
        }

        std::string line;
        bool has_header = false;

        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;

            if (!has_header && (line.find("gene") != std::string::npos ||
                                line.find("LOEUF") != std::string::npos ||
                                line.find("oe_lof") != std::string::npos)) {
                has_header = true;
                continue;
            }

            std::istringstream iss(line);
            std::string gene;
            double loeuf;

            if (!(iss >> gene >> loeuf)) continue;

            // Check if we already have data for this gene
            if (gene_data_.count(gene) > 0) {
                gene_data_[gene].oe_lof_upper = loeuf;
            } else {
                GeneConstraint constraint;
                constraint.gene_symbol = gene;
                constraint.oe_lof_upper = loeuf;
                gene_data_[gene] = constraint;
            }
        }

        loaded_ = true;
        return true;
    }

    /**
     * Get constraint data by gene symbol
     */
    GeneConstraint get_by_symbol(const std::string& gene_symbol) const {
        auto it = gene_data_.find(gene_symbol);
        if (it != gene_data_.end()) {
            return it->second;
        }
        return GeneConstraint();
    }

    /**
     * Get constraint data by Ensembl gene ID
     */
    GeneConstraint get_by_gene_id(const std::string& gene_id) const {
        auto it = gene_id_data_.find(gene_id);
        if (it != gene_id_data_.end()) {
            return it->second;
        }
        return GeneConstraint();
    }

    /**
     * Check if database is loaded
     */
    bool is_loaded() const { return loaded_; }

    /**
     * Get number of genes
     */
    size_t size() const { return gene_data_.size(); }

    /**
     * Get all highly constrained genes (pLI > 0.9 or LOEUF < 0.35)
     */
    std::vector<std::string> get_constrained_genes(double pli_threshold = 0.9,
                                                    double loeuf_threshold = 0.35) const {
        std::vector<std::string> result;
        for (auto it = gene_data_.begin(); it != gene_data_.end(); ++it) {
            const GeneConstraint& c = it->second;
            if (c.pLI >= pli_threshold ||
                (c.oe_lof_upper >= 0 && c.oe_lof_upper < loeuf_threshold)) {
                result.push_back(it->first);
            }
        }
        return result;
    }

private:
    std::map<std::string, GeneConstraint> gene_data_;
    std::map<std::string, GeneConstraint> gene_id_data_;
    bool loaded_;
};

/**
 * Singleton accessor for gene constraint database
 */
inline GeneConstraintDB& get_gene_constraint_db() {
    static GeneConstraintDB db;
    return db;
}

/**
 * Format constraint score for output
 */
inline std::string format_constraint_score(double value, int precision = 4) {
    if (value < 0 || std::isnan(value)) {
        return ".";
    }

    std::ostringstream oss;
    oss.precision(precision);
    oss << std::fixed << value;

    std::string result = oss.str();

    // Remove trailing zeros
    size_t dot_pos = result.find('.');
    if (dot_pos != std::string::npos) {
        size_t last_nonzero = result.find_last_not_of('0');
        if (last_nonzero != std::string::npos && last_nonzero > dot_pos) {
            result = result.substr(0, last_nonzero + 1);
        }
        if (result[result.size() - 1] == '.') {
            result = result.substr(0, result.size() - 1);
        }
    }

    return result;
}

} // namespace vep

#endif // GENE_CONSTRAINT_HPP
