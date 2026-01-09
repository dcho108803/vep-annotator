/**
 * filter_vep Equivalent
 *
 * Post-processing filter utility for VEP annotations.
 * Supports filtering by consequence, impact, gene list, frequency, etc.
 */

#ifndef FILTER_VEP_HPP
#define FILTER_VEP_HPP

#include <string>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cctype>
#include <cmath>

namespace vep {

/**
 * Filter operators
 */
enum class FilterOperator {
    EQUALS,         // eq, =, is
    NOT_EQUALS,     // ne, !=
    GREATER,        // gt, >
    GREATER_EQ,     // ge, >=
    LESS,           // lt, <
    LESS_EQ,        // le, <=
    CONTAINS,       // contains, match
    NOT_CONTAINS,   // not contains
    IN,             // in (list)
    NOT_IN,         // not in
    EXISTS,         // exists, defined
    NOT_EXISTS,     // not exists
    REGEX           // regex match
};

/**
 * Parse filter operator from string
 */
inline FilterOperator parse_filter_operator(const std::string& op_str) {
    std::string lower = op_str;
    for (size_t i = 0; i < lower.size(); ++i) {
        lower[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(lower[i])));
    }

    if (lower == "eq" || lower == "=" || lower == "is") return FilterOperator::EQUALS;
    if (lower == "ne" || lower == "!=") return FilterOperator::NOT_EQUALS;
    if (lower == "gt" || lower == ">") return FilterOperator::GREATER;
    if (lower == "ge" || lower == ">=") return FilterOperator::GREATER_EQ;
    if (lower == "lt" || lower == "<") return FilterOperator::LESS;
    if (lower == "le" || lower == "<=") return FilterOperator::LESS_EQ;
    if (lower == "contains" || lower == "match") return FilterOperator::CONTAINS;
    if (lower == "in") return FilterOperator::IN;
    if (lower == "exists" || lower == "defined") return FilterOperator::EXISTS;
    if (lower == "regex" || lower == "re") return FilterOperator::REGEX;

    return FilterOperator::EQUALS;
}

/**
 * Single filter condition
 */
struct FilterCondition {
    std::string field;
    FilterOperator op;
    std::string value;
    std::vector<std::string> value_list;  // For IN operator
    bool negated = false;

    FilterCondition() : op(FilterOperator::EQUALS), negated(false) {}

    FilterCondition(const std::string& f, FilterOperator o, const std::string& v)
        : field(f), op(o), value(v), negated(false) {}
};

/**
 * Filter configuration
 */
struct FilterConfig {
    std::vector<FilterCondition> conditions;
    bool match_all = true;  // AND vs OR

    // Quick filters
    std::set<std::string> consequence_filter;  // Filter by consequence types
    std::set<std::string> impact_filter;       // Filter by impact (HIGH, MODERATE, etc.)
    std::set<std::string> gene_filter;         // Filter by gene list
    std::set<std::string> biotype_filter;      // Filter by biotype

    // Numeric filters
    double min_af = -1.0;
    double max_af = -1.0;
    double min_cadd = -1.0;
    double min_revel = -1.0;

    // Boolean flags
    bool coding_only = false;
    bool exclude_intergenic = false;
    bool exclude_intronic = false;
    bool canonical_only = false;
    bool mane_only = false;
    bool pick_one = false;

    bool has_any_filter() const {
        return !conditions.empty() ||
               !consequence_filter.empty() ||
               !impact_filter.empty() ||
               !gene_filter.empty() ||
               !biotype_filter.empty() ||
               min_af >= 0 || max_af >= 0 ||
               min_cadd >= 0 || min_revel >= 0 ||
               coding_only || exclude_intergenic ||
               exclude_intronic || canonical_only ||
               mane_only || pick_one;
    }
};

/**
 * Represents a single annotation record for filtering
 */
struct FilterableRecord {
    std::map<std::string, std::string> fields;
    std::string original_line;

    std::string get(const std::string& field) const {
        auto it = fields.find(field);
        if (it != fields.end()) {
            return it->second;
        }
        return "";
    }

    bool has(const std::string& field) const {
        return fields.count(field) > 0 && !fields.at(field).empty();
    }

    double get_numeric(const std::string& field) const {
        std::string val = get(field);
        if (val.empty() || val == "." || val == "NA" || val == "NaN") {
            return std::numeric_limits<double>::quiet_NaN();
        }
        try {
            return std::stod(val);
        } catch (...) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
};

/**
 * Apply a single filter condition
 */
inline bool apply_condition(const FilterableRecord& record, const FilterCondition& cond) {
    std::string value = record.get(cond.field);
    bool result = false;

    // Handle EXISTS operator specially
    if (cond.op == FilterOperator::EXISTS) {
        result = record.has(cond.field);
        return cond.negated ? !result : result;
    }

    if (cond.op == FilterOperator::NOT_EXISTS) {
        result = !record.has(cond.field);
        return cond.negated ? !result : result;
    }

    // Try numeric comparison first
    bool is_numeric = true;
    double num_value = 0, num_target = 0;

    if (value.empty() || value == "." || value == "NA") {
        is_numeric = false;
    } else {
        try {
            num_value = std::stod(value);
            num_target = std::stod(cond.value);
        } catch (...) {
            is_numeric = false;
        }
    }

    if (cond.op == FilterOperator::EQUALS) {
        if (is_numeric) {
            result = (std::abs(num_value - num_target) < 1e-9);
        } else {
            result = (value == cond.value);
        }
    } else if (cond.op == FilterOperator::NOT_EQUALS) {
        if (is_numeric) {
            result = (std::abs(num_value - num_target) >= 1e-9);
        } else {
            result = (value != cond.value);
        }
    } else if (cond.op == FilterOperator::GREATER) {
        if (is_numeric) {
            result = (num_value > num_target);
        }
    } else if (cond.op == FilterOperator::GREATER_EQ) {
        if (is_numeric) {
            result = (num_value >= num_target);
        }
    } else if (cond.op == FilterOperator::LESS) {
        if (is_numeric) {
            result = (num_value < num_target);
        }
    } else if (cond.op == FilterOperator::LESS_EQ) {
        if (is_numeric) {
            result = (num_value <= num_target);
        }
    } else if (cond.op == FilterOperator::CONTAINS) {
        result = (value.find(cond.value) != std::string::npos);
    } else if (cond.op == FilterOperator::NOT_CONTAINS) {
        result = (value.find(cond.value) == std::string::npos);
    } else if (cond.op == FilterOperator::IN) {
        for (size_t i = 0; i < cond.value_list.size(); ++i) {
            if (value == cond.value_list[i]) {
                result = true;
                break;
            }
        }
    } else if (cond.op == FilterOperator::NOT_IN) {
        result = true;
        for (size_t i = 0; i < cond.value_list.size(); ++i) {
            if (value == cond.value_list[i]) {
                result = false;
                break;
            }
        }
    } else if (cond.op == FilterOperator::REGEX) {
        // Simple regex support (just contains for now)
        result = (value.find(cond.value) != std::string::npos);
    }

    return cond.negated ? !result : result;
}

/**
 * Apply all filter conditions
 */
inline bool apply_filter(const FilterableRecord& record, const FilterConfig& config) {
    // Quick filters first
    if (!config.consequence_filter.empty()) {
        std::string consequence = record.get("CONSEQUENCE");
        if (consequence.empty()) consequence = record.get("Consequence");

        bool found = false;
        for (auto it = config.consequence_filter.begin(); it != config.consequence_filter.end(); ++it) {
            if (consequence.find(*it) != std::string::npos) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }

    if (!config.impact_filter.empty()) {
        std::string impact = record.get("IMPACT");
        if (impact.empty()) impact = record.get("Impact");

        if (config.impact_filter.count(impact) == 0) {
            return false;
        }
    }

    if (!config.gene_filter.empty()) {
        std::string gene = record.get("GENE");
        if (gene.empty()) gene = record.get("Gene");
        if (gene.empty()) gene = record.get("SYMBOL");

        if (config.gene_filter.count(gene) == 0) {
            return false;
        }
    }

    if (!config.biotype_filter.empty()) {
        std::string biotype = record.get("BIOTYPE");
        if (biotype.empty()) biotype = record.get("Biotype");

        if (config.biotype_filter.count(biotype) == 0) {
            return false;
        }
    }

    // Numeric filters
    if (config.min_af >= 0 || config.max_af >= 0) {
        double af = record.get_numeric("AF");
        if (std::isnan(af)) {
            af = record.get_numeric("gnomAD_AF");
        }
        if (std::isnan(af)) {
            af = record.get_numeric("gnomAD:AF");
        }

        if (!std::isnan(af)) {
            if (config.min_af >= 0 && af < config.min_af) return false;
            if (config.max_af >= 0 && af > config.max_af) return false;
        }
    }

    if (config.min_cadd >= 0) {
        double cadd = record.get_numeric("CADD_phred");
        if (std::isnan(cadd)) {
            cadd = record.get_numeric("CADD");
        }
        if (!std::isnan(cadd) && cadd < config.min_cadd) {
            return false;
        }
    }

    if (config.min_revel >= 0) {
        double revel = record.get_numeric("REVEL_score");
        if (std::isnan(revel)) {
            revel = record.get_numeric("REVEL");
        }
        if (!std::isnan(revel) && revel < config.min_revel) {
            return false;
        }
    }

    // Boolean filters
    if (config.coding_only) {
        std::string biotype = record.get("BIOTYPE");
        if (biotype.empty()) biotype = record.get("Biotype");
        if (biotype.find("protein_coding") == std::string::npos) {
            return false;
        }
    }

    if (config.exclude_intergenic) {
        std::string consequence = record.get("CONSEQUENCE");
        if (consequence.empty()) consequence = record.get("Consequence");
        if (consequence.find("intergenic") != std::string::npos) {
            return false;
        }
    }

    if (config.exclude_intronic) {
        std::string consequence = record.get("CONSEQUENCE");
        if (consequence.empty()) consequence = record.get("Consequence");
        if (consequence.find("intron_variant") != std::string::npos) {
            // Only filter if it's ONLY an intron variant
            if (consequence == "intron_variant") {
                return false;
            }
        }
    }

    if (config.canonical_only) {
        std::string canonical = record.get("CANONICAL");
        if (canonical != "YES" && canonical != "1" && canonical != "true") {
            return false;
        }
    }

    if (config.mane_only) {
        std::string mane = record.get("MANE_SELECT");
        if (mane.empty()) mane = record.get("MANE");
        if (mane.empty() || mane == "." || mane == "NA") {
            return false;
        }
    }

    // Custom conditions
    if (!config.conditions.empty()) {
        if (config.match_all) {
            // AND logic
            for (size_t i = 0; i < config.conditions.size(); ++i) {
                if (!apply_condition(record, config.conditions[i])) {
                    return false;
                }
            }
        } else {
            // OR logic
            bool any_match = false;
            for (size_t i = 0; i < config.conditions.size(); ++i) {
                if (apply_condition(record, config.conditions[i])) {
                    any_match = true;
                    break;
                }
            }
            if (!any_match) return false;
        }
    }

    return true;
}

/**
 * Parse filter expression
 * Format: FIELD OP VALUE
 * Examples:
 *   "IMPACT is HIGH"
 *   "AF < 0.01"
 *   "Consequence contains missense"
 */
inline FilterCondition parse_filter_expression(const std::string& expr) {
    FilterCondition cond;

    // Find operator
    std::vector<std::string> operators;
    operators.push_back(" is ");
    operators.push_back(" eq ");
    operators.push_back(" ne ");
    operators.push_back(" gt ");
    operators.push_back(" ge ");
    operators.push_back(" lt ");
    operators.push_back(" le ");
    operators.push_back(" contains ");
    operators.push_back(" in ");
    operators.push_back(" match ");
    operators.push_back(" exists");
    operators.push_back(">=");
    operators.push_back("<=");
    operators.push_back("!=");
    operators.push_back(">");
    operators.push_back("<");
    operators.push_back("=");

    size_t op_pos = std::string::npos;
    std::string found_op;

    for (size_t i = 0; i < operators.size(); ++i) {
        size_t pos = expr.find(operators[i]);
        if (pos != std::string::npos && (op_pos == std::string::npos || pos < op_pos)) {
            op_pos = pos;
            found_op = operators[i];
        }
    }

    if (op_pos == std::string::npos) {
        // No operator found, treat as field exists check
        cond.field = expr;
        cond.op = FilterOperator::EXISTS;
        return cond;
    }

    // Extract field
    cond.field = expr.substr(0, op_pos);
    // Trim whitespace
    while (!cond.field.empty() && std::isspace(cond.field[cond.field.size() - 1])) {
        cond.field.erase(cond.field.size() - 1);
    }
    while (!cond.field.empty() && std::isspace(cond.field[0])) {
        cond.field.erase(0, 1);
    }

    // Handle "not" prefix
    if (cond.field.size() > 4 && cond.field.substr(0, 4) == "not ") {
        cond.field = cond.field.substr(4);
        cond.negated = true;
    }

    // Parse operator
    std::string op_str = found_op;
    // Trim whitespace
    while (!op_str.empty() && std::isspace(op_str[op_str.size() - 1])) {
        op_str.erase(op_str.size() - 1);
    }
    while (!op_str.empty() && std::isspace(op_str[0])) {
        op_str.erase(0, 1);
    }

    cond.op = parse_filter_operator(op_str);

    // Extract value
    if (op_pos + found_op.size() < expr.size()) {
        cond.value = expr.substr(op_pos + found_op.size());
        // Trim whitespace
        while (!cond.value.empty() && std::isspace(cond.value[cond.value.size() - 1])) {
            cond.value.erase(cond.value.size() - 1);
        }
        while (!cond.value.empty() && std::isspace(cond.value[0])) {
            cond.value.erase(0, 1);
        }

        // Handle IN operator - parse comma-separated list
        if (cond.op == FilterOperator::IN || cond.op == FilterOperator::NOT_IN) {
            std::istringstream iss(cond.value);
            std::string item;
            while (std::getline(iss, item, ',')) {
                // Trim whitespace
                while (!item.empty() && std::isspace(item[item.size() - 1])) {
                    item.erase(item.size() - 1);
                }
                while (!item.empty() && std::isspace(item[0])) {
                    item.erase(0, 1);
                }
                if (!item.empty()) {
                    cond.value_list.push_back(item);
                }
            }
        }
    }

    return cond;
}

/**
 * Load gene list from file
 */
inline std::set<std::string> load_gene_list(const std::string& filepath) {
    std::set<std::string> genes;
    std::ifstream file(filepath);

    if (!file.is_open()) {
        return genes;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        // Trim whitespace
        while (!line.empty() && std::isspace(line[line.size() - 1])) {
            line.erase(line.size() - 1);
        }
        while (!line.empty() && std::isspace(line[0])) {
            line.erase(0, 1);
        }

        if (!line.empty()) {
            // Handle TSV format (take first column)
            size_t tab_pos = line.find('\t');
            if (tab_pos != std::string::npos) {
                line = line.substr(0, tab_pos);
            }
            genes.insert(line);
        }
    }

    return genes;
}

/**
 * Parse TSV header to get column indices
 */
inline std::map<std::string, int> parse_tsv_header(const std::string& header_line) {
    std::map<std::string, int> col_map;
    std::istringstream iss(header_line);
    std::string col;
    int idx = 0;

    while (std::getline(iss, col, '\t')) {
        col_map[col] = idx++;
    }

    return col_map;
}

/**
 * Parse TSV line into FilterableRecord
 */
inline FilterableRecord parse_tsv_record(const std::string& line,
                                          const std::map<std::string, int>& col_map) {
    FilterableRecord record;
    record.original_line = line;

    std::vector<std::string> fields;
    std::istringstream iss(line);
    std::string field;
    while (std::getline(iss, field, '\t')) {
        fields.push_back(field);
    }

    for (auto it = col_map.begin(); it != col_map.end(); ++it) {
        if (it->second < static_cast<int>(fields.size())) {
            record.fields[it->first] = fields[it->second];
        }
    }

    return record;
}

/**
 * Filter TSV file
 */
inline int filter_tsv_file(const std::string& input_path,
                           const std::string& output_path,
                           const FilterConfig& config) {
    std::ifstream input(input_path);
    if (!input.is_open()) {
        return -1;
    }

    std::ofstream output(output_path);
    if (!output.is_open()) {
        return -1;
    }

    std::string line;
    std::map<std::string, int> col_map;
    int lines_passed = 0;
    int lines_total = 0;

    while (std::getline(input, line)) {
        // Handle header
        if (line.empty()) continue;
        if (line[0] == '#') {
            output << line << "\n";
            continue;
        }

        // First non-comment line is header
        if (col_map.empty()) {
            col_map = parse_tsv_header(line);
            output << line << "\n";
            continue;
        }

        lines_total++;

        FilterableRecord record = parse_tsv_record(line, col_map);

        if (apply_filter(record, config)) {
            output << line << "\n";
            lines_passed++;
        }
    }

    return lines_passed;
}

} // namespace vep

#endif // FILTER_VEP_HPP
