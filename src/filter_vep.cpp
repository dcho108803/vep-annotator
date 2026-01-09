/**
 * filter_vep - Standalone Filter Tool for VEP Annotations
 *
 * Filters VEP annotation output files by consequence, impact, gene list,
 * frequency, and custom filter expressions.
 */

#include "filter_vep.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <set>
#include <string>

void print_filter_usage(const char* program_name) {
    std::cout << "filter_vep - Filter VEP Annotation Output\n"
              << "==========================================\n\n"
              << "Usage: " << program_name << " [OPTIONS] -i INPUT -o OUTPUT\n\n"
              << "Required:\n"
              << "  -i, --input FILE         Input TSV file from vep_annotator\n"
              << "  -o, --output FILE        Output filtered TSV file\n\n"
              << "Filter Options:\n"
              << "  -f, --filter EXPR        Filter expression (can be used multiple times)\n"
              << "                           Format: FIELD OPERATOR VALUE\n"
              << "                           Operators: is, eq, ne, gt, ge, lt, le, contains, in\n"
              << "                           Examples:\n"
              << "                             'IMPACT is HIGH'\n"
              << "                             'AF < 0.01'\n"
              << "                             'Consequence contains missense'\n"
              << "                             'GENE in BRCA1,BRCA2,TP53'\n\n"
              << "Quick Filters:\n"
              << "  --consequence LIST       Filter by consequence types (comma-separated)\n"
              << "  --impact LIST            Filter by impact levels (HIGH,MODERATE,LOW,MODIFIER)\n"
              << "  --gene-list FILE         Filter by genes from file (one per line)\n"
              << "  --biotype LIST           Filter by biotype (comma-separated)\n\n"
              << "Numeric Filters:\n"
              << "  --min-af VALUE           Minimum allele frequency\n"
              << "  --max-af VALUE           Maximum allele frequency\n"
              << "  --min-cadd VALUE         Minimum CADD score\n"
              << "  --min-revel VALUE        Minimum REVEL score\n\n"
              << "Boolean Filters:\n"
              << "  --coding-only            Only protein_coding transcripts\n"
              << "  --exclude-intergenic     Exclude intergenic variants\n"
              << "  --exclude-intronic       Exclude intron variants\n"
              << "  --canonical-only         Only canonical transcripts\n"
              << "  --mane-only              Only MANE Select transcripts\n"
              << "  --pick                   Output only one annotation per variant\n\n"
              << "Logic Options:\n"
              << "  --and                    All filter conditions must match (default)\n"
              << "  --or                     Any filter condition can match\n\n"
              << "Other Options:\n"
              << "  -h, --help               Show this help message\n"
              << "  --count                  Only print count of matching records\n"
              << "  --list-columns           List available columns in input file\n\n"
              << "Examples:\n"
              << "  # Filter for HIGH impact variants\n"
              << "  " << program_name << " -i vep_output.tsv -o filtered.tsv --impact HIGH\n\n"
              << "  # Filter for rare missense variants\n"
              << "  " << program_name << " -i vep_output.tsv -o filtered.tsv \\\n"
              << "      --consequence missense_variant --max-af 0.01\n\n"
              << "  # Filter for specific genes with CADD > 20\n"
              << "  " << program_name << " -i vep_output.tsv -o filtered.tsv \\\n"
              << "      --gene-list genes.txt --min-cadd 20\n\n"
              << "  # Complex filter expression\n"
              << "  " << program_name << " -i vep_output.tsv -o filtered.tsv \\\n"
              << "      -f 'IMPACT in HIGH,MODERATE' -f 'AF < 0.01' --canonical-only\n"
              << std::endl;
}

std::set<std::string> parse_comma_list(const std::string& list_str) {
    std::set<std::string> result;
    std::istringstream iss(list_str);
    std::string item;

    while (std::getline(iss, item, ',')) {
        // Trim whitespace
        size_t start = item.find_first_not_of(" \t");
        size_t end = item.find_last_not_of(" \t");
        if (start != std::string::npos && end != std::string::npos) {
            result.insert(item.substr(start, end - start + 1));
        }
    }

    return result;
}

void list_columns(const std::string& input_path) {
    std::ifstream input(input_path);
    if (!input.is_open()) {
        std::cerr << "Error: Cannot open file: " << input_path << std::endl;
        return;
    }

    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;

        // First non-comment line is header
        std::cout << "Available columns:\n";
        std::istringstream iss(line);
        std::string col;
        int idx = 1;
        while (std::getline(iss, col, '\t')) {
            std::cout << "  " << idx++ << ". " << col << "\n";
        }
        break;
    }
}

int main(int argc, char* argv[]) {
    std::string input_path;
    std::string output_path;
    vep::FilterConfig config;
    bool count_only = false;
    bool list_cols = false;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_filter_usage(argv[0]);
            return 0;
        } else if ((arg == "-i" || arg == "--input") && i + 1 < argc) {
            input_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            output_path = argv[++i];
        } else if ((arg == "-f" || arg == "--filter") && i + 1 < argc) {
            std::string expr = argv[++i];
            vep::FilterCondition cond = vep::parse_filter_expression(expr);
            config.conditions.push_back(cond);
        } else if (arg == "--consequence" && i + 1 < argc) {
            config.consequence_filter = parse_comma_list(argv[++i]);
        } else if (arg == "--impact" && i + 1 < argc) {
            config.impact_filter = parse_comma_list(argv[++i]);
        } else if (arg == "--gene-list" && i + 1 < argc) {
            config.gene_filter = vep::load_gene_list(argv[++i]);
            if (config.gene_filter.empty()) {
                std::cerr << "Warning: No genes loaded from file" << std::endl;
            }
        } else if (arg == "--biotype" && i + 1 < argc) {
            config.biotype_filter = parse_comma_list(argv[++i]);
        } else if (arg == "--min-af" && i + 1 < argc) {
            config.min_af = std::stod(argv[++i]);
        } else if (arg == "--max-af" && i + 1 < argc) {
            config.max_af = std::stod(argv[++i]);
        } else if (arg == "--min-cadd" && i + 1 < argc) {
            config.min_cadd = std::stod(argv[++i]);
        } else if (arg == "--min-revel" && i + 1 < argc) {
            config.min_revel = std::stod(argv[++i]);
        } else if (arg == "--coding-only") {
            config.coding_only = true;
        } else if (arg == "--exclude-intergenic") {
            config.exclude_intergenic = true;
        } else if (arg == "--exclude-intronic") {
            config.exclude_intronic = true;
        } else if (arg == "--canonical-only") {
            config.canonical_only = true;
        } else if (arg == "--mane-only") {
            config.mane_only = true;
        } else if (arg == "--pick") {
            config.pick_one = true;
        } else if (arg == "--and") {
            config.match_all = true;
        } else if (arg == "--or") {
            config.match_all = false;
        } else if (arg == "--count") {
            count_only = true;
        } else if (arg == "--list-columns") {
            list_cols = true;
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            print_filter_usage(argv[0]);
            return 1;
        }
    }

    // Validate required arguments
    if (input_path.empty()) {
        std::cerr << "Error: Input file (-i) is required.\n" << std::endl;
        print_filter_usage(argv[0]);
        return 1;
    }

    // Handle --list-columns
    if (list_cols) {
        list_columns(input_path);
        return 0;
    }

    if (output_path.empty() && !count_only) {
        std::cerr << "Error: Output file (-o) is required.\n" << std::endl;
        print_filter_usage(argv[0]);
        return 1;
    }

    // Check if any filter is specified
    if (!config.has_any_filter()) {
        std::cerr << "Warning: No filters specified. All records will be output.\n" << std::endl;
    }

    // Open input file
    std::ifstream input(input_path);
    if (!input.is_open()) {
        std::cerr << "Error: Cannot open input file: " << input_path << std::endl;
        return 1;
    }

    // Open output file (if not count-only mode)
    std::ofstream output;
    if (!count_only) {
        output.open(output_path);
        if (!output.is_open()) {
            std::cerr << "Error: Cannot open output file: " << output_path << std::endl;
            return 1;
        }
    }

    std::string line;
    std::map<std::string, int> col_map;
    int lines_passed = 0;
    int lines_total = 0;

    // Track variants for --pick
    std::set<std::string> seen_variants;

    while (std::getline(input, line)) {
        // Handle empty lines and comments
        if (line.empty()) continue;
        if (line[0] == '#') {
            if (!count_only) output << line << "\n";
            continue;
        }

        // First non-comment line is header
        if (col_map.empty()) {
            col_map = vep::parse_tsv_header(line);
            if (!count_only) output << line << "\n";
            continue;
        }

        lines_total++;

        vep::FilterableRecord record = vep::parse_tsv_record(line, col_map);

        // Apply --pick filter
        if (config.pick_one) {
            std::string variant_key = record.get("CHROM") + ":" +
                                      record.get("POS") + ":" +
                                      record.get("REF") + ":" +
                                      record.get("ALT");
            if (variant_key.empty()) {
                // Try alternative column names
                variant_key = record.get("Chromosome") + ":" +
                              record.get("Position") + ":" +
                              record.get("Ref_allele") + ":" +
                              record.get("Alt_allele");
            }

            if (seen_variants.count(variant_key) > 0) {
                continue;  // Skip - already output one for this variant
            }

            if (vep::apply_filter(record, config)) {
                seen_variants.insert(variant_key);
            }
        }

        if (vep::apply_filter(record, config)) {
            if (!count_only) output << line << "\n";
            lines_passed++;
        }
    }

    // Print summary
    std::cout << "Filtered " << lines_passed << " of " << lines_total << " records";
    if (lines_total > 0) {
        double pct = 100.0 * lines_passed / lines_total;
        std::cout << " (" << std::fixed << std::setprecision(1) << pct << "%)";
    }
    std::cout << std::endl;

    if (!count_only) {
        std::cout << "Output written to: " << output_path << std::endl;
    }

    return 0;
}
