/**
 * VEP Variant Annotator - Pure C++ Local Implementation
 */

#include "vep_annotator.hpp"
#include "annotation_source.hpp"
#include "file_parsers.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <ctime>
#include <iomanip>
#include <regex>
#include <zlib.h>

// Tabix support via htslib (optional - compile with -DHAVE_HTSLIB)
#ifdef HAVE_HTSLIB
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>
#endif

namespace vep {

// ============================================================================
// Logging
// ============================================================================

static LogLevel g_log_level = LogLevel::INFO;

void set_log_level(LogLevel level) {
    g_log_level = level;
}

void log(LogLevel level, const std::string& message) {
    if (level < g_log_level) return;

    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);

    const char* level_str;
    switch (level) {
        case LogLevel::DEBUG:   level_str = "DEBUG"; break;
        case LogLevel::INFO:    level_str = "INFO"; break;
        case LogLevel::WARNING: level_str = "WARNING"; break;
        case LogLevel::ERROR:   level_str = "ERROR"; break;
        default:                level_str = "UNKNOWN"; break;
    }

    std::cerr << std::put_time(std::localtime(&time_t_now), "%Y-%m-%d %H:%M:%S")
              << " - " << level_str << " - " << message << std::endl;
}

// ============================================================================
// Consequence type utilities
// ============================================================================

std::string consequence_to_string(ConsequenceType type) {
    switch (type) {
        case ConsequenceType::TRANSCRIPT_ABLATION: return "transcript_ablation";
        case ConsequenceType::SPLICE_ACCEPTOR_VARIANT: return "splice_acceptor_variant";
        case ConsequenceType::SPLICE_DONOR_VARIANT: return "splice_donor_variant";
        case ConsequenceType::STOP_GAINED: return "stop_gained";
        case ConsequenceType::FRAMESHIFT_VARIANT: return "frameshift_variant";
        case ConsequenceType::STOP_LOST: return "stop_lost";
        case ConsequenceType::START_LOST: return "start_lost";
        case ConsequenceType::INFRAME_INSERTION: return "inframe_insertion";
        case ConsequenceType::INFRAME_DELETION: return "inframe_deletion";
        case ConsequenceType::MISSENSE_VARIANT: return "missense_variant";
        case ConsequenceType::PROTEIN_ALTERING_VARIANT: return "protein_altering_variant";
        case ConsequenceType::SPLICE_REGION_VARIANT: return "splice_region_variant";
        case ConsequenceType::INCOMPLETE_TERMINAL_CODON_VARIANT: return "incomplete_terminal_codon_variant";
        case ConsequenceType::START_RETAINED_VARIANT: return "start_retained_variant";
        case ConsequenceType::STOP_RETAINED_VARIANT: return "stop_retained_variant";
        case ConsequenceType::SYNONYMOUS_VARIANT: return "synonymous_variant";
        case ConsequenceType::CODING_SEQUENCE_VARIANT: return "coding_sequence_variant";
        case ConsequenceType::MATURE_MIRNA_VARIANT: return "mature_miRNA_variant";
        case ConsequenceType::FIVE_PRIME_UTR_VARIANT: return "5_prime_UTR_variant";
        case ConsequenceType::THREE_PRIME_UTR_VARIANT: return "3_prime_UTR_variant";
        case ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT: return "non_coding_transcript_exon_variant";
        case ConsequenceType::INTRON_VARIANT: return "intron_variant";
        case ConsequenceType::NMD_TRANSCRIPT_VARIANT: return "NMD_transcript_variant";
        case ConsequenceType::NON_CODING_TRANSCRIPT_VARIANT: return "non_coding_transcript_variant";
        case ConsequenceType::UPSTREAM_GENE_VARIANT: return "upstream_gene_variant";
        case ConsequenceType::DOWNSTREAM_GENE_VARIANT: return "downstream_gene_variant";
        case ConsequenceType::INTERGENIC_VARIANT: return "intergenic_variant";
        default: return "unknown";
    }
}

Impact get_impact(ConsequenceType type) {
    switch (type) {
        case ConsequenceType::TRANSCRIPT_ABLATION:
        case ConsequenceType::SPLICE_ACCEPTOR_VARIANT:
        case ConsequenceType::SPLICE_DONOR_VARIANT:
        case ConsequenceType::STOP_GAINED:
        case ConsequenceType::FRAMESHIFT_VARIANT:
        case ConsequenceType::STOP_LOST:
        case ConsequenceType::START_LOST:
            return Impact::HIGH;

        case ConsequenceType::INFRAME_INSERTION:
        case ConsequenceType::INFRAME_DELETION:
        case ConsequenceType::MISSENSE_VARIANT:
        case ConsequenceType::PROTEIN_ALTERING_VARIANT:
            return Impact::MODERATE;

        case ConsequenceType::SPLICE_REGION_VARIANT:
        case ConsequenceType::INCOMPLETE_TERMINAL_CODON_VARIANT:
        case ConsequenceType::START_RETAINED_VARIANT:
        case ConsequenceType::STOP_RETAINED_VARIANT:
        case ConsequenceType::SYNONYMOUS_VARIANT:
            return Impact::LOW;

        default:
            return Impact::MODIFIER;
    }
}

std::string impact_to_string(Impact impact) {
    switch (impact) {
        case Impact::HIGH: return "HIGH";
        case Impact::MODERATE: return "MODERATE";
        case Impact::LOW: return "LOW";
        case Impact::MODIFIER: return "MODIFIER";
        default: return "UNKNOWN";
    }
}

// ============================================================================
// Codon Table
// ============================================================================

char CodonTable::translate(const std::string& codon) {
    static const std::unordered_map<std::string, char> table = {
        {"TTT", 'F'}, {"TTC", 'F'}, {"TTA", 'L'}, {"TTG", 'L'},
        {"CTT", 'L'}, {"CTC", 'L'}, {"CTA", 'L'}, {"CTG", 'L'},
        {"ATT", 'I'}, {"ATC", 'I'}, {"ATA", 'I'}, {"ATG", 'M'},
        {"GTT", 'V'}, {"GTC", 'V'}, {"GTA", 'V'}, {"GTG", 'V'},
        {"TCT", 'S'}, {"TCC", 'S'}, {"TCA", 'S'}, {"TCG", 'S'},
        {"CCT", 'P'}, {"CCC", 'P'}, {"CCA", 'P'}, {"CCG", 'P'},
        {"ACT", 'T'}, {"ACC", 'T'}, {"ACA", 'T'}, {"ACG", 'T'},
        {"GCT", 'A'}, {"GCC", 'A'}, {"GCA", 'A'}, {"GCG", 'A'},
        {"TAT", 'Y'}, {"TAC", 'Y'}, {"TAA", '*'}, {"TAG", '*'},
        {"CAT", 'H'}, {"CAC", 'H'}, {"CAA", 'Q'}, {"CAG", 'Q'},
        {"AAT", 'N'}, {"AAC", 'N'}, {"AAA", 'K'}, {"AAG", 'K'},
        {"GAT", 'D'}, {"GAC", 'D'}, {"GAA", 'E'}, {"GAG", 'E'},
        {"TGT", 'C'}, {"TGC", 'C'}, {"TGA", '*'}, {"TGG", 'W'},
        {"CGT", 'R'}, {"CGC", 'R'}, {"CGA", 'R'}, {"CGG", 'R'},
        {"AGT", 'S'}, {"AGC", 'S'}, {"AGA", 'R'}, {"AGG", 'R'},
        {"GGT", 'G'}, {"GGC", 'G'}, {"GGA", 'G'}, {"GGG", 'G'}
    };

    if (codon.length() != 3) return 'X';

    std::string upper_codon = codon;
    std::transform(upper_codon.begin(), upper_codon.end(), upper_codon.begin(), ::toupper);

    auto it = table.find(upper_codon);
    return (it != table.end()) ? it->second : 'X';
}

std::string CodonTable::get_three_letter(char aa) {
    static const std::unordered_map<char, std::string> table = {
        {'A', "Ala"}, {'R', "Arg"}, {'N', "Asn"}, {'D', "Asp"},
        {'C', "Cys"}, {'E', "Glu"}, {'Q', "Gln"}, {'G', "Gly"},
        {'H', "His"}, {'I', "Ile"}, {'L', "Leu"}, {'K', "Lys"},
        {'M', "Met"}, {'F', "Phe"}, {'P', "Pro"}, {'S', "Ser"},
        {'T', "Thr"}, {'W', "Trp"}, {'Y', "Tyr"}, {'V', "Val"},
        {'*', "Ter"}, {'X', "Xaa"}
    };

    auto it = table.find(aa);
    return (it != table.end()) ? it->second : "Xaa";
}

bool CodonTable::is_start_codon(const std::string& codon) {
    std::string upper = codon;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
    return upper == "ATG";
}

bool CodonTable::is_stop_codon(const std::string& codon) {
    std::string upper = codon;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
    return upper == "TAA" || upper == "TAG" || upper == "TGA";
}

// ============================================================================
// Transcript methods
// ============================================================================

int Transcript::get_cds_length() const {
    int length = 0;
    for (const auto& cds : cds_regions) {
        length += cds.end - cds.start + 1;
    }
    return length;
}

// ============================================================================
// VariantAnnotation methods
// ============================================================================

std::string VariantAnnotation::get_consequence_string() const {
    std::ostringstream oss;
    bool first = true;
    for (const auto& c : consequences) {
        if (!first) oss << ",";
        oss << consequence_to_string(c);
        first = false;
    }
    return oss.str();
}

std::map<std::string, std::string> VariantAnnotation::get_summary() const {
    std::map<std::string, std::string> summary;

    summary["input"] = input_variant;
    summary["location"] = chromosome + ":" + std::to_string(position);
    summary["alleles"] = ref_allele + "/" + alt_allele;
    summary["gene"] = gene_symbol;
    summary["transcript"] = transcript_id;
    summary["consequence"] = get_consequence_string();
    summary["impact"] = impact_to_string(impact);
    summary["biotype"] = biotype;

    if (cds_position > 0) {
        summary["cds_position"] = std::to_string(cds_position);
    }
    if (protein_position > 0) {
        summary["protein_position"] = std::to_string(protein_position);
    }
    if (!amino_acids.empty()) {
        summary["amino_acids"] = amino_acids;
    }
    if (!codons.empty()) {
        summary["codons"] = codons;
    }
    if (!hgvsc.empty()) {
        summary["hgvsc"] = hgvsc;
    }
    if (!hgvsp.empty()) {
        summary["hgvsp"] = hgvsp;
    }

    // Add custom annotations to summary
    for (const auto& [key, value] : custom_annotations) {
        summary[key] = value;
    }

    return summary;
}

std::optional<std::string> VariantAnnotation::get_annotation(
    const std::string& source, const std::string& field) const {
    std::string key = source + ":" + field;
    auto it = custom_annotations.find(key);
    if (it != custom_annotations.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<double> VariantAnnotation::get_annotation_double(
    const std::string& source, const std::string& field) const {
    auto value = get_annotation(source, field);
    if (value.has_value()) {
        try {
            return std::stod(value.value());
        } catch (...) {
            return std::nullopt;
        }
    }
    return std::nullopt;
}

// ============================================================================
// TabixVCFReader - On-disk queries for large VCF files
// ============================================================================

#ifdef HAVE_HTSLIB

/**
 * TabixVCFReader - Query tabix-indexed VCF files without loading into memory
 */
class TabixVCFReader {
public:
    TabixVCFReader(const std::string& vcf_path, const std::vector<std::string>& info_fields)
        : vcf_path_(vcf_path), info_fields_(info_fields) {

        // Open the VCF file
        hts_file_ = hts_open(vcf_path.c_str(), "r");
        if (!hts_file_) {
            throw std::runtime_error("Cannot open VCF file: " + vcf_path);
        }

        // Load the tabix index
        tbx_ = tbx_index_load(vcf_path.c_str());
        if (!tbx_) {
            hts_close(hts_file_);
            throw std::runtime_error("Cannot load tabix index for: " + vcf_path +
                                    ". Make sure the file is bgzip compressed and indexed with tabix.");
        }

        extract_all_ = info_fields_.empty();
        for (const auto& f : info_fields_) {
            fields_to_extract_.insert(f);
        }

        log(LogLevel::INFO, "Opened tabix-indexed VCF: " + vcf_path);
    }

    ~TabixVCFReader() {
        if (tbx_) tbx_destroy(tbx_);
        if (hts_file_) hts_close(hts_file_);
    }

    // Prevent copying
    TabixVCFReader(const TabixVCFReader&) = delete;
    TabixVCFReader& operator=(const TabixVCFReader&) = delete;

    /**
     * Query for records at a specific position
     */
    std::vector<VCFRecord> query(const std::string& chrom, int pos) {
        std::vector<VCFRecord> results;

        // Normalize chromosome name
        std::string query_chrom = chrom;

        // Try with "chr" prefix first, then without
        std::vector<std::string> chrom_variants;
        if (chrom.substr(0, 3) == "chr") {
            chrom_variants.push_back(chrom);
            chrom_variants.push_back(chrom.substr(3));
        } else {
            chrom_variants.push_back("chr" + chrom);
            chrom_variants.push_back(chrom);
        }

        for (const auto& try_chrom : chrom_variants) {
            // Build region string: chr:pos-pos
            std::string region = try_chrom + ":" + std::to_string(pos) + "-" + std::to_string(pos);

            hts_itr_t* itr = tbx_itr_querys(tbx_, region.c_str());
            if (!itr) continue;

            kstring_t str = {0, 0, nullptr};

            while (tbx_itr_next(hts_file_, tbx_, itr, &str) >= 0) {
                VCFRecord record = parse_vcf_line(str.s);
                if (record.pos == pos) {
                    results.push_back(std::move(record));
                }
            }

            free(str.s);
            tbx_itr_destroy(itr);

            if (!results.empty()) break;
        }

        return results;
    }

private:
    std::string vcf_path_;
    std::vector<std::string> info_fields_;
    std::set<std::string> fields_to_extract_;
    bool extract_all_ = false;

    htsFile* hts_file_ = nullptr;
    tbx_t* tbx_ = nullptr;

    VCFRecord parse_vcf_line(const char* line) {
        VCFRecord record;

        std::istringstream iss(line);
        std::string qual, filter, info_str;

        if (!(iss >> record.chrom >> record.pos >> record.id >> record.ref)) {
            return record;
        }

        std::string alt;
        if (!(iss >> alt >> qual >> filter >> info_str)) {
            return record;
        }

        // Normalize chromosome
        if (record.chrom.length() > 3 && record.chrom.substr(0, 3) == "chr") {
            record.chrom = record.chrom.substr(3);
        }

        // Parse alt alleles
        std::istringstream alt_iss(alt);
        std::string single_alt;
        while (std::getline(alt_iss, single_alt, ',')) {
            record.alts.push_back(single_alt);
        }

        // Parse INFO field
        if (info_str != "." && !info_str.empty()) {
            std::istringstream info_iss(info_str);
            std::string item;

            while (std::getline(info_iss, item, ';')) {
                size_t eq_pos = item.find('=');
                if (eq_pos != std::string::npos) {
                    std::string key = item.substr(0, eq_pos);
                    std::string value = item.substr(eq_pos + 1);

                    if (extract_all_ || fields_to_extract_.count(key)) {
                        record.info[key] = value;
                    }
                } else {
                    // Flag field
                    if (extract_all_ || fields_to_extract_.count(item)) {
                        record.info[item] = "1";
                    }
                }
            }
        }

        return record;
    }
};

#endif // HAVE_HTSLIB

// ============================================================================
// VCFAnnotationDatabase implementation
// ============================================================================

struct VCFAnnotationDatabase::Impl {
    // Source configurations
    std::vector<VCFAnnotationConfig> configs;

    // Index: source -> chrom -> pos -> vector of records (for in-memory sources)
    std::unordered_map<std::string,
        std::unordered_map<std::string,
            std::map<int, std::vector<VCFRecord>>>> index;

    // Track which fields are available per source
    std::unordered_map<std::string, std::set<std::string>> available_fields;

    // Statistics
    std::unordered_map<std::string, size_t> record_counts;

#ifdef HAVE_HTSLIB
    // Tabix readers for on-disk sources (large files)
    std::unordered_map<std::string, std::unique_ptr<TabixVCFReader>> tabix_readers;
#endif

    // Track which sources use tabix mode
    std::set<std::string> tabix_sources;

    static std::string normalize_chrom(const std::string& chrom) {
        if (chrom.length() > 3 && chrom.substr(0, 3) == "chr") {
            return chrom.substr(3);
        }
        return chrom;
    }

    static std::map<std::string, std::string> parse_info(const std::string& info_str) {
        std::map<std::string, std::string> info;
        if (info_str == "." || info_str.empty()) return info;

        std::istringstream iss(info_str);
        std::string item;

        while (std::getline(iss, item, ';')) {
            size_t eq_pos = item.find('=');
            if (eq_pos != std::string::npos) {
                std::string key = item.substr(0, eq_pos);
                std::string value = item.substr(eq_pos + 1);
                info[key] = value;
            } else {
                // Flag field (no value)
                info[item] = "1";
            }
        }

        return info;
    }

    static std::vector<std::string> split(const std::string& s, char delim) {
        std::vector<std::string> result;
        std::istringstream iss(s);
        std::string item;
        while (std::getline(iss, item, delim)) {
            result.push_back(item);
        }
        return result;
    }
};

VCFAnnotationDatabase::VCFAnnotationDatabase()
    : pimpl_(std::make_unique<Impl>()) {}

VCFAnnotationDatabase::~VCFAnnotationDatabase() = default;

void VCFAnnotationDatabase::add_source(const std::string& name,
                                        const std::string& vcf_path,
                                        const std::string& info_fields) {
    VCFAnnotationConfig config;
    config.name = name;
    config.vcf_path = vcf_path;

    if (!info_fields.empty()) {
        std::istringstream iss(info_fields);
        std::string field;
        while (std::getline(iss, field, ',')) {
            // Trim whitespace
            field.erase(0, field.find_first_not_of(" \t"));
            field.erase(field.find_last_not_of(" \t") + 1);
            if (!field.empty()) {
                config.info_fields.push_back(field);
            }
        }
    }

    add_source(config);
}

void VCFAnnotationDatabase::add_source(const VCFAnnotationConfig& config) {
    pimpl_->configs.push_back(config);

#ifdef HAVE_HTSLIB
    // Use tabix mode for large files
    if (config.use_tabix) {
        log(LogLevel::INFO, "Opening tabix-indexed VCF source '" + config.name + "': " + config.vcf_path);

        try {
            auto reader = std::make_unique<TabixVCFReader>(config.vcf_path, config.info_fields);
            pimpl_->tabix_readers[config.name] = std::move(reader);
            pimpl_->tabix_sources.insert(config.name);
            pimpl_->record_counts[config.name] = 0;  // Unknown count for tabix

            // Store requested fields as available (we don't know all fields without scanning)
            for (const auto& field : config.info_fields) {
                pimpl_->available_fields[config.name].insert(field);
            }

            log(LogLevel::INFO, "Tabix source '" + config.name + "' ready for on-disk queries");
            return;
        } catch (const std::exception& e) {
            log(LogLevel::WARNING, "Failed to open tabix index: " + std::string(e.what()) +
                                  ". Falling back to in-memory loading.");
        }
    }
#else
    if (config.use_tabix) {
        log(LogLevel::WARNING, "Tabix support not compiled in. Build with -DHAVE_HTSLIB and link htslib. "
                              "Falling back to in-memory loading for: " + config.vcf_path);
    }
#endif

    // In-memory loading (original behavior)
    log(LogLevel::INFO, "Loading VCF annotation source '" + config.name + "' from: " + config.vcf_path);

    bool is_gzipped = (config.vcf_path.size() >= 3 &&
                       config.vcf_path.substr(config.vcf_path.size() - 3) == ".gz");

    std::function<bool(std::string&)> read_line;
    gzFile gz = nullptr;
    std::ifstream file;

    if (is_gzipped) {
        gz = gzopen(config.vcf_path.c_str(), "rb");
        if (!gz) {
            throw std::runtime_error("Cannot open VCF file: " + config.vcf_path);
        }

        static char buffer[65536];
        read_line = [&gz](std::string& line) -> bool {
            if (gzgets(gz, buffer, sizeof(buffer)) == nullptr) return false;
            line = buffer;
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
                line.pop_back();
            }
            return true;
        };
    } else {
        file.open(config.vcf_path);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open VCF file: " + config.vcf_path);
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
    size_t record_count = 0;
    std::set<std::string> fields_to_extract(config.info_fields.begin(), config.info_fields.end());
    bool extract_all = config.info_fields.empty();

    while (read_line(line)) {
        if (line.empty() || line[0] == '#') continue;

        // Parse VCF line
        std::istringstream iss(line);
        std::string chrom, id, ref, alt, qual, filter, info_str;
        int pos;

        if (!(iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info_str)) {
            continue;
        }

        chrom = Impl::normalize_chrom(chrom);

        // Parse INFO field
        auto all_info = Impl::parse_info(info_str);

        // Filter to requested fields
        std::map<std::string, std::string> filtered_info;
        for (const auto& [key, value] : all_info) {
            if (extract_all || fields_to_extract.count(key)) {
                filtered_info[key] = value;
                pimpl_->available_fields[config.name].insert(key);
            }
        }

        // Parse alt alleles
        std::vector<std::string> alts = Impl::split(alt, ',');

        // Create record
        VCFRecord record;
        record.chrom = chrom;
        record.pos = pos;
        record.id = id;
        record.ref = ref;
        record.alts = alts;
        record.info = filtered_info;

        // Add to index
        pimpl_->index[config.name][chrom][pos].push_back(std::move(record));
        record_count++;

        if (record_count % 1000000 == 0) {
            log(LogLevel::DEBUG, "Loaded " + std::to_string(record_count) + " records from " + config.name);
        }
    }

    if (gz) gzclose(gz);

    pimpl_->record_counts[config.name] = record_count;

    log(LogLevel::INFO, "Loaded " + std::to_string(record_count) + " records from '" + config.name + "'");
}

std::map<std::string, std::string> VCFAnnotationDatabase::get_annotations(
    const std::string& chrom, int pos,
    const std::string& ref, const std::string& alt) const {

    std::map<std::string, std::string> result;
    std::string norm_chrom = Impl::normalize_chrom(chrom);

    // Helper lambda to process records
    auto process_records = [&](const std::string& source_name, const VCFAnnotationConfig& config,
                               const std::vector<VCFRecord>& records) {
        for (const auto& record : records) {
            // Check allele match if required
            if (config.match_allele) {
                if (record.ref != ref) continue;

                bool alt_match = false;
                int alt_index = -1;
                for (size_t i = 0; i < record.alts.size(); ++i) {
                    if (record.alts[i] == alt) {
                        alt_match = true;
                        alt_index = static_cast<int>(i);
                        break;
                    }
                }
                if (!alt_match) continue;

                // Extract annotations, handling multi-allelic sites
                for (const auto& [field, value] : record.info) {
                    std::string key = source_name + ":" + field;

                    // Check if value is comma-separated (per-allele)
                    if (alt_index >= 0 && value.find(',') != std::string::npos) {
                        auto values = Impl::split(value, ',');
                        if (alt_index < static_cast<int>(values.size())) {
                            result[key] = values[alt_index];
                        } else if (!values.empty()) {
                            result[key] = values[0];
                        }
                    } else {
                        result[key] = value;
                    }
                }
            } else {
                // Position-only match
                for (const auto& [field, value] : record.info) {
                    std::string key = source_name + ":" + field;
                    result[key] = value;
                }
            }
        }
    };

    for (const auto& config : pimpl_->configs) {
#ifdef HAVE_HTSLIB
        // Check if this source uses tabix
        if (pimpl_->tabix_sources.count(config.name)) {
            auto reader_it = pimpl_->tabix_readers.find(config.name);
            if (reader_it != pimpl_->tabix_readers.end()) {
                // Query tabix index on-the-fly
                auto records = reader_it->second->query(chrom, pos);
                process_records(config.name, config, records);
            }
            continue;
        }
#endif

        // In-memory lookup
        auto source_it = pimpl_->index.find(config.name);
        if (source_it == pimpl_->index.end()) continue;

        auto chrom_it = source_it->second.find(norm_chrom);
        if (chrom_it == source_it->second.end()) continue;

        auto pos_it = chrom_it->second.find(pos);
        if (pos_it == chrom_it->second.end()) continue;

        process_records(config.name, config, pos_it->second);
    }

    return result;
}

std::vector<const VCFRecord*> VCFAnnotationDatabase::get_records_at(
    const std::string& chrom, int pos) const {

    std::vector<const VCFRecord*> result;
    std::string norm_chrom = Impl::normalize_chrom(chrom);

    for (const auto& [source, chrom_map] : pimpl_->index) {
        auto chrom_it = chrom_map.find(norm_chrom);
        if (chrom_it == chrom_map.end()) continue;

        auto pos_it = chrom_it->second.find(pos);
        if (pos_it == chrom_it->second.end()) continue;

        for (const auto& record : pos_it->second) {
            result.push_back(&record);
        }
    }

    return result;
}

std::vector<std::string> VCFAnnotationDatabase::get_sources() const {
    std::vector<std::string> sources;
    for (const auto& config : pimpl_->configs) {
        sources.push_back(config.name);
    }
    return sources;
}

std::vector<std::string> VCFAnnotationDatabase::get_fields(const std::string& source) const {
    std::vector<std::string> fields;
    auto it = pimpl_->available_fields.find(source);
    if (it != pimpl_->available_fields.end()) {
        fields.assign(it->second.begin(), it->second.end());
    }
    return fields;
}

std::string VCFAnnotationDatabase::get_stats() const {
    std::ostringstream oss;
    oss << "VCF Annotation Sources:\n";
    for (const auto& config : pimpl_->configs) {
        oss << "  " << config.name << ": ";

        // Check if using tabix mode
        if (pimpl_->tabix_sources.count(config.name)) {
            oss << "[tabix indexed, on-disk queries]";
        } else {
            oss << pimpl_->record_counts.at(config.name) << " records";
        }

        auto fields = get_fields(config.name);
        if (!fields.empty()) {
            oss << " (fields: ";
            for (size_t i = 0; i < fields.size() && i < 5; ++i) {
                if (i > 0) oss << ", ";
                oss << fields[i];
            }
            if (fields.size() > 5) {
                oss << ", ... +" << (fields.size() - 5) << " more";
            }
            oss << ")";
        }
        oss << "\n";
    }
    return oss.str();
}

// ============================================================================
// ReferenceGenome implementation
// ============================================================================

struct ReferenceGenome::Impl {
    std::unordered_map<std::string, std::string> sequences;
    std::unordered_map<std::string, int> lengths;
    std::string fasta_path;
    bool loaded_all = false;

    // For indexed access (not loading all into memory)
    std::unordered_map<std::string, std::streampos> index;
};

ReferenceGenome::ReferenceGenome(const std::string& fasta_path, bool load_all)
    : pimpl_(std::make_unique<Impl>()) {

    pimpl_->fasta_path = fasta_path;
    pimpl_->loaded_all = load_all;

    log(LogLevel::INFO, "Loading reference genome from: " + fasta_path);

    // Check if file is gzipped
    bool is_gzipped = (fasta_path.size() >= 3 &&
                       fasta_path.substr(fasta_path.size() - 3) == ".gz");

    if (is_gzipped) {
        gzFile gz = gzopen(fasta_path.c_str(), "rb");
        if (!gz) {
            throw std::runtime_error("Cannot open gzipped FASTA file: " + fasta_path);
        }

        std::string current_chrom;
        std::string current_seq;
        char buffer[8192];

        while (gzgets(gz, buffer, sizeof(buffer)) != nullptr) {
            std::string line(buffer);
            // Remove trailing newline
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
                line.pop_back();
            }

            if (line.empty()) continue;

            if (line[0] == '>') {
                // Save previous chromosome
                if (!current_chrom.empty() && !current_seq.empty()) {
                    pimpl_->sequences[current_chrom] = std::move(current_seq);
                    pimpl_->lengths[current_chrom] = pimpl_->sequences[current_chrom].length();
                }

                // Parse new chromosome name
                size_t space_pos = line.find(' ');
                current_chrom = line.substr(1, space_pos - 1);

                // Normalize chromosome name (remove "chr" prefix for consistency)
                if (current_chrom.substr(0, 3) == "chr") {
                    current_chrom = current_chrom.substr(3);
                }

                current_seq.clear();
                current_seq.reserve(300000000);  // Reserve ~300MB for human chromosomes
            } else {
                // Convert to uppercase
                std::transform(line.begin(), line.end(), line.begin(), ::toupper);
                current_seq += line;
            }
        }

        // Save last chromosome
        if (!current_chrom.empty() && !current_seq.empty()) {
            pimpl_->sequences[current_chrom] = std::move(current_seq);
            pimpl_->lengths[current_chrom] = pimpl_->sequences[current_chrom].length();
        }

        gzclose(gz);
    } else {
        std::ifstream file(fasta_path);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open FASTA file: " + fasta_path);
        }

        std::string current_chrom;
        std::string current_seq;
        std::string line;

        while (std::getline(file, line)) {
            // Remove trailing carriage return if present
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }

            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!current_chrom.empty() && !current_seq.empty()) {
                    pimpl_->sequences[current_chrom] = std::move(current_seq);
                    pimpl_->lengths[current_chrom] = pimpl_->sequences[current_chrom].length();
                }

                size_t space_pos = line.find(' ');
                current_chrom = line.substr(1, space_pos - 1);

                if (current_chrom.substr(0, 3) == "chr") {
                    current_chrom = current_chrom.substr(3);
                }

                current_seq.clear();
                current_seq.reserve(300000000);
            } else {
                std::transform(line.begin(), line.end(), line.begin(), ::toupper);
                current_seq += line;
            }
        }

        if (!current_chrom.empty() && !current_seq.empty()) {
            pimpl_->sequences[current_chrom] = std::move(current_seq);
            pimpl_->lengths[current_chrom] = pimpl_->sequences[current_chrom].length();
        }
    }

    log(LogLevel::INFO, "Loaded " + std::to_string(pimpl_->sequences.size()) + " chromosomes");
}

ReferenceGenome::~ReferenceGenome() = default;

std::string ReferenceGenome::get_sequence(const std::string& chrom, int start, int end) const {
    std::string normalized_chrom = chrom;
    if (normalized_chrom.substr(0, 3) == "chr") {
        normalized_chrom = normalized_chrom.substr(3);
    }

    auto it = pimpl_->sequences.find(normalized_chrom);
    if (it == pimpl_->sequences.end()) {
        return "";
    }

    const std::string& seq = it->second;

    // Convert to 0-based indexing
    int start_idx = start - 1;
    int end_idx = end;

    if (start_idx < 0) start_idx = 0;
    if (end_idx > static_cast<int>(seq.length())) end_idx = seq.length();

    if (start_idx >= end_idx || start_idx >= static_cast<int>(seq.length())) {
        return "";
    }

    return seq.substr(start_idx, end_idx - start_idx);
}

char ReferenceGenome::get_base(const std::string& chrom, int pos) const {
    std::string seq = get_sequence(chrom, pos, pos);
    return seq.empty() ? 'N' : seq[0];
}

bool ReferenceGenome::has_chromosome(const std::string& chrom) const {
    std::string normalized = chrom;
    if (normalized.substr(0, 3) == "chr") {
        normalized = normalized.substr(3);
    }
    return pimpl_->sequences.find(normalized) != pimpl_->sequences.end();
}

int ReferenceGenome::get_chromosome_length(const std::string& chrom) const {
    std::string normalized = chrom;
    if (normalized.substr(0, 3) == "chr") {
        normalized = normalized.substr(3);
    }
    auto it = pimpl_->lengths.find(normalized);
    return (it != pimpl_->lengths.end()) ? it->second : 0;
}

// ============================================================================
// TranscriptDatabase implementation
// ============================================================================

struct TranscriptDatabase::Impl {
    std::unordered_map<std::string, Transcript> transcripts;
    std::unordered_map<std::string, Gene> genes;

    // Spatial index: chromosome -> sorted list of (start, end, transcript_id)
    std::unordered_map<std::string, std::vector<std::tuple<int, int, std::string>>> chrom_index;

    void build_index() {
        for (const auto& [tid, transcript] : transcripts) {
            chrom_index[transcript.chromosome].emplace_back(
                transcript.start, transcript.end, tid
            );
        }

        // Sort by start position
        for (auto& [chrom, entries] : chrom_index) {
            std::sort(entries.begin(), entries.end());
        }
    }
};

// Helper to parse GTF attribute string
static std::unordered_map<std::string, std::string> parse_gtf_attributes(const std::string& attr_str) {
    std::unordered_map<std::string, std::string> attrs;

    std::regex attr_regex(R"((\w+)\s+\"([^\"]*)\")");
    std::sregex_iterator it(attr_str.begin(), attr_str.end(), attr_regex);
    std::sregex_iterator end;

    while (it != end) {
        attrs[(*it)[1].str()] = (*it)[2].str();
        ++it;
    }

    return attrs;
}

TranscriptDatabase::TranscriptDatabase(const std::string& gtf_path)
    : pimpl_(std::make_unique<Impl>()) {

    log(LogLevel::INFO, "Loading annotations from: " + gtf_path);

    bool is_gzipped = (gtf_path.size() >= 3 &&
                       gtf_path.substr(gtf_path.size() - 3) == ".gz");

    std::function<bool(std::string&)> read_line;
    gzFile gz = nullptr;
    std::ifstream file;

    if (is_gzipped) {
        gz = gzopen(gtf_path.c_str(), "rb");
        if (!gz) {
            throw std::runtime_error("Cannot open GTF file: " + gtf_path);
        }

        static char buffer[65536];
        read_line = [&gz](std::string& line) -> bool {
            if (gzgets(gz, buffer, sizeof(buffer)) == nullptr) return false;
            line = buffer;
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
                line.pop_back();
            }
            return true;
        };
    } else {
        file.open(gtf_path);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open GTF file: " + gtf_path);
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
    int line_count = 0;

    while (read_line(line)) {
        if (line.empty() || line[0] == '#') continue;

        line_count++;
        if (line_count % 500000 == 0) {
            log(LogLevel::DEBUG, "Processed " + std::to_string(line_count) + " lines...");
        }

        // Parse GTF line
        std::istringstream iss(line);
        std::string chrom, source, feature;
        int start, end;
        std::string score, strand_str, frame, attributes;

        if (!(iss >> chrom >> source >> feature >> start >> end >> score >> strand_str >> frame)) {
            continue;
        }

        std::getline(iss, attributes);

        // Normalize chromosome
        if (chrom.substr(0, 3) == "chr") {
            chrom = chrom.substr(3);
        }

        char strand = strand_str[0];
        auto attrs = parse_gtf_attributes(attributes);

        std::string gene_id = attrs.count("gene_id") ? attrs["gene_id"] : "";
        std::string gene_name = attrs.count("gene_name") ? attrs["gene_name"] : "";
        std::string transcript_id = attrs.count("transcript_id") ? attrs["transcript_id"] : "";
        std::string biotype = attrs.count("gene_biotype") ? attrs["gene_biotype"] :
                             (attrs.count("transcript_biotype") ? attrs["transcript_biotype"] : "");

        if (feature == "gene" && !gene_id.empty()) {
            Gene& gene = pimpl_->genes[gene_id];
            gene.id = gene_id;
            gene.name = gene_name;
            gene.chromosome = chrom;
            gene.start = start;
            gene.end = end;
            gene.strand = strand;
            gene.biotype = biotype;
        }
        else if (feature == "transcript" && !transcript_id.empty()) {
            Transcript& tr = pimpl_->transcripts[transcript_id];
            tr.id = transcript_id;
            tr.gene_id = gene_id;
            tr.gene_name = gene_name;
            tr.chromosome = chrom;
            tr.start = start;
            tr.end = end;
            tr.strand = strand;
            tr.biotype = biotype;

            // Check for canonical tag
            if (attrs.count("tag")) {
                tr.is_canonical = (attrs["tag"].find("canonical") != std::string::npos ||
                                  attrs["tag"].find("MANE") != std::string::npos);
            }

            // Add transcript to gene
            if (pimpl_->genes.count(gene_id)) {
                pimpl_->genes[gene_id].transcript_ids.push_back(transcript_id);
            }
        }
        else if (feature == "exon" && !transcript_id.empty()) {
            if (pimpl_->transcripts.count(transcript_id)) {
                Exon exon;
                exon.start = start;
                exon.end = end;
                exon.exon_number = attrs.count("exon_number") ?
                    std::stoi(attrs["exon_number"]) : 0;

                pimpl_->transcripts[transcript_id].exons.push_back(exon);
            }
        }
        else if (feature == "CDS" && !transcript_id.empty()) {
            if (pimpl_->transcripts.count(transcript_id)) {
                CDS cds;
                cds.start = start;
                cds.end = end;
                cds.phase = (frame != ".") ? std::stoi(frame) : 0;

                Transcript& tr = pimpl_->transcripts[transcript_id];
                tr.cds_regions.push_back(cds);

                // Update CDS boundaries
                if (tr.cds_start == 0 || start < tr.cds_start) {
                    tr.cds_start = start;
                }
                if (end > tr.cds_end) {
                    tr.cds_end = end;
                }
            }
        }
        else if ((feature == "start_codon" || feature == "stop_codon") && !transcript_id.empty()) {
            if (pimpl_->transcripts.count(transcript_id)) {
                Transcript& tr = pimpl_->transcripts[transcript_id];
                if (feature == "start_codon") {
                    if (tr.strand == '+') {
                        tr.cds_start = start;
                    } else {
                        tr.cds_end = end;
                    }
                }
            }
        }
    }

    if (gz) gzclose(gz);

    // Sort exons and CDS regions by position
    for (auto& [tid, tr] : pimpl_->transcripts) {
        std::sort(tr.exons.begin(), tr.exons.end(),
                  [](const Exon& a, const Exon& b) { return a.start < b.start; });
        std::sort(tr.cds_regions.begin(), tr.cds_regions.end(),
                  [](const CDS& a, const CDS& b) { return a.start < b.start; });

        // Assign exon numbers if not present
        for (size_t i = 0; i < tr.exons.size(); ++i) {
            if (tr.exons[i].exon_number == 0) {
                tr.exons[i].exon_number = (tr.strand == '+') ? (i + 1) : (tr.exons.size() - i);
            }
        }
    }

    // Build spatial index
    pimpl_->build_index();

    log(LogLevel::INFO, "Loaded " + std::to_string(pimpl_->genes.size()) + " genes, "
                      + std::to_string(pimpl_->transcripts.size()) + " transcripts");
}

TranscriptDatabase::~TranscriptDatabase() = default;

std::vector<const Transcript*> TranscriptDatabase::get_transcripts_at(
    const std::string& chrom, int pos) const {
    return get_transcripts_in_region(chrom, pos, pos);
}

std::vector<const Transcript*> TranscriptDatabase::get_transcripts_in_region(
    const std::string& chrom, int start, int end) const {

    std::string normalized = chrom;
    if (normalized.substr(0, 3) == "chr") {
        normalized = normalized.substr(3);
    }

    std::vector<const Transcript*> results;

    auto it = pimpl_->chrom_index.find(normalized);
    if (it == pimpl_->chrom_index.end()) {
        return results;
    }

    for (const auto& [tr_start, tr_end, tid] : it->second) {
        // Check overlap
        if (tr_start <= end && tr_end >= start) {
            auto tr_it = pimpl_->transcripts.find(tid);
            if (tr_it != pimpl_->transcripts.end()) {
                results.push_back(&tr_it->second);
            }
        }
        // Since sorted by start, can break early
        if (tr_start > end) break;
    }

    return results;
}

const Transcript* TranscriptDatabase::get_transcript(const std::string& transcript_id) const {
    auto it = pimpl_->transcripts.find(transcript_id);
    return (it != pimpl_->transcripts.end()) ? &it->second : nullptr;
}

const Gene* TranscriptDatabase::get_gene(const std::string& gene_id) const {
    auto it = pimpl_->genes.find(gene_id);
    return (it != pimpl_->genes.end()) ? &it->second : nullptr;
}

std::vector<const Transcript*> TranscriptDatabase::get_gene_transcripts(
    const std::string& gene_id) const {

    std::vector<const Transcript*> results;
    auto gene_it = pimpl_->genes.find(gene_id);
    if (gene_it == pimpl_->genes.end()) {
        return results;
    }

    for (const auto& tid : gene_it->second.transcript_ids) {
        auto tr_it = pimpl_->transcripts.find(tid);
        if (tr_it != pimpl_->transcripts.end()) {
            results.push_back(&tr_it->second);
        }
    }

    return results;
}

std::vector<const Gene*> TranscriptDatabase::get_nearby_genes(
    const std::string& chrom, int pos, int distance) const {

    std::vector<const Gene*> results;
    std::string normalized = chrom;
    if (normalized.substr(0, 3) == "chr") {
        normalized = normalized.substr(3);
    }

    for (const auto& [gene_id, gene] : pimpl_->genes) {
        if (gene.chromosome != normalized) continue;

        int dist = 0;
        if (pos < gene.start) {
            dist = gene.start - pos;
        } else if (pos > gene.end) {
            dist = pos - gene.end;
        }

        if (dist <= distance) {
            results.push_back(&gene);
        }
    }

    return results;
}

size_t TranscriptDatabase::transcript_count() const {
    return pimpl_->transcripts.size();
}

size_t TranscriptDatabase::gene_count() const {
    return pimpl_->genes.size();
}

// ============================================================================
// VEPAnnotator implementation
// ============================================================================

VEPAnnotator::VEPAnnotator(const std::string& gtf_path, const std::string& fasta_path)
    : transcript_db_(std::make_unique<TranscriptDatabase>(gtf_path)),
      reference_(std::make_unique<ReferenceGenome>(fasta_path)),
      vcf_annotations_(std::make_unique<VCFAnnotationDatabase>()),
      source_manager_(std::make_unique<AnnotationSourceManager>()) {

    log(LogLevel::INFO, "VEP Annotator initialized");
}

VEPAnnotator::~VEPAnnotator() = default;

void VEPAnnotator::add_annotation_source(const std::string& name,
                                          const std::string& vcf_path,
                                          const std::string& info_fields) {
    vcf_annotations_->add_source(name, vcf_path, info_fields);
}

void VEPAnnotator::add_annotation_source(const VCFAnnotationConfig& config) {
    vcf_annotations_->add_source(config);
}

std::vector<std::string> VEPAnnotator::get_annotation_sources() const {
    return vcf_annotations_->get_sources();
}

std::vector<VariantAnnotation> VEPAnnotator::annotate(
    const std::string& chrom,
    int pos,
    const std::string& ref,
    const std::string& alt) {

    std::vector<VariantAnnotation> results;

    // Get overlapping transcripts
    int var_end = pos + static_cast<int>(ref.length()) - 1;
    auto transcripts = transcript_db_->get_transcripts_in_region(chrom, pos, var_end);

    if (transcripts.empty()) {
        // Intergenic variant - check for nearby genes
        VariantAnnotation ann;
        ann.input_variant = chrom + "-" + std::to_string(pos) + "-" + ref + "-" + alt;
        ann.chromosome = chrom;
        ann.position = pos;
        ann.ref_allele = ref;
        ann.alt_allele = alt;
        ann.consequences.push_back(ConsequenceType::INTERGENIC_VARIANT);
        ann.impact = Impact::MODIFIER;

        // Find nearby genes for upstream/downstream annotation
        auto nearby = transcript_db_->get_nearby_genes(chrom, pos, 5000);
        if (!nearby.empty()) {
            const Gene* closest = nearby[0];
            ann.gene_symbol = closest->name;
            ann.gene_id = closest->id;

            if (pos < closest->start) {
                ann.consequences[0] = (closest->strand == '+') ?
                    ConsequenceType::UPSTREAM_GENE_VARIANT :
                    ConsequenceType::DOWNSTREAM_GENE_VARIANT;
            } else {
                ann.consequences[0] = (closest->strand == '+') ?
                    ConsequenceType::DOWNSTREAM_GENE_VARIANT :
                    ConsequenceType::UPSTREAM_GENE_VARIANT;
            }
        }

        results.push_back(std::move(ann));
    } else {
        for (const auto* transcript : transcripts) {
            auto ann = annotate_transcript(chrom, pos, ref, alt, *transcript);
            results.push_back(std::move(ann));
        }
    }

    // Apply custom VCF annotations to all results
    auto vcf_annots = vcf_annotations_->get_annotations(chrom, pos, ref, alt);
    if (!vcf_annots.empty()) {
        for (auto& ann : results) {
            for (const auto& [key, value] : vcf_annots) {
                ann.custom_annotations[key] = value;
            }
        }
    }

    // Apply annotation sources (dbNSFP, conservation, etc.)
    if (source_manager_) {
        for (auto& ann : results) {
            // Get transcript pointer for transcript-specific annotations
            const Transcript* tx_ptr = nullptr;
            if (!ann.transcript_id.empty()) {
                tx_ptr = transcript_db_->get_transcript(ann.transcript_id);
            }

            // Call all enabled annotation sources
            source_manager_->annotate_all(chrom, pos, ref, alt, tx_ptr, ann.custom_annotations);
        }
    }

    return results;
}

VariantAnnotation VEPAnnotator::annotate_most_severe(
    const std::string& chrom,
    int pos,
    const std::string& ref,
    const std::string& alt) {

    auto all_annotations = annotate(chrom, pos, ref, alt);

    if (all_annotations.empty()) {
        VariantAnnotation empty;
        empty.input_variant = chrom + "-" + std::to_string(pos) + "-" + ref + "-" + alt;
        return empty;
    }

    // Sort by impact severity
    std::sort(all_annotations.begin(), all_annotations.end(),
              [](const VariantAnnotation& a, const VariantAnnotation& b) {
                  return static_cast<int>(a.impact) < static_cast<int>(b.impact);
              });

    return all_annotations[0];
}

VariantAnnotation VEPAnnotator::annotate_transcript(
    const std::string& chrom,
    int pos,
    const std::string& ref,
    const std::string& alt,
    const Transcript& transcript) {

    VariantAnnotation ann;
    ann.input_variant = chrom + "-" + std::to_string(pos) + "-" + ref + "-" + alt;
    ann.chromosome = chrom;
    ann.position = pos;
    ann.ref_allele = ref;
    ann.alt_allele = alt;
    ann.gene_symbol = transcript.gene_name;
    ann.gene_id = transcript.gene_id;
    ann.transcript_id = transcript.id;
    ann.biotype = transcript.biotype;
    ann.is_canonical = transcript.is_canonical;

    // Determine consequences
    ann.consequences = determine_consequences(pos, ref, alt, transcript);

    // Get most severe impact
    ann.impact = Impact::MODIFIER;
    for (const auto& c : ann.consequences) {
        Impact i = get_impact(c);
        if (static_cast<int>(i) < static_cast<int>(ann.impact)) {
            ann.impact = i;
        }
    }

    // Calculate CDS position and codon changes for coding variants
    if (transcript.is_coding()) {
        int cds_pos = calculate_cds_position(pos, transcript);
        if (cds_pos > 0) {
            ann.cds_position = cds_pos;
            ann.protein_position = (cds_pos - 1) / 3 + 1;

            // Get affected codons
            auto [ref_codon, alt_codon] = get_affected_codons(cds_pos, ref, alt, transcript);
            if (!ref_codon.empty() && !alt_codon.empty()) {
                ann.codons = ref_codon + "/" + alt_codon;

                char ref_aa = CodonTable::translate(ref_codon);
                char alt_aa = CodonTable::translate(alt_codon);
                ann.amino_acids = std::string(1, ref_aa) + "/" + std::string(1, alt_aa);

                // Generate HGVS notations
                ann.hgvsc = generate_hgvsc(cds_pos, ref, alt, transcript);
                ann.hgvsp = generate_hgvsp(
                    CodonTable::get_three_letter(ref_aa),
                    CodonTable::get_three_letter(alt_aa),
                    ann.protein_position,
                    transcript
                );
            }
        }
    }

    // Find exon/intron number
    for (size_t i = 0; i < transcript.exons.size(); ++i) {
        const auto& exon = transcript.exons[i];
        if (pos >= exon.start && pos <= exon.end) {
            ann.exon_number = exon.exon_number;
            break;
        }
        if (i + 1 < transcript.exons.size()) {
            const auto& next_exon = transcript.exons[i + 1];
            if (pos > exon.end && pos < next_exon.start) {
                ann.intron_number = (transcript.strand == '+') ? (i + 1) : (transcript.exons.size() - i - 1);
                break;
            }
        }
    }

    return ann;
}

std::vector<ConsequenceType> VEPAnnotator::determine_consequences(
    int pos,
    const std::string& ref,
    const std::string& alt,
    const Transcript& transcript) {

    std::vector<ConsequenceType> consequences;
    int var_end = pos + static_cast<int>(ref.length()) - 1;

    // Check if variant is in exon or intron
    bool in_exon = false;
    bool in_intron = false;
    bool near_splice = false;

    for (size_t i = 0; i < transcript.exons.size(); ++i) {
        const auto& exon = transcript.exons[i];

        if (pos >= exon.start && pos <= exon.end) {
            in_exon = true;
        }

        // Check splice sites (2bp at each end)
        if (i > 0) {
            // Acceptor site (3' end of intron, 5' end of exon)
            if (pos >= exon.start - 2 && pos <= exon.start - 1) {
                consequences.push_back(ConsequenceType::SPLICE_ACCEPTOR_VARIANT);
                near_splice = true;
            }
        }
        if (i < transcript.exons.size() - 1) {
            // Donor site (3' end of exon, 5' end of intron)
            if (pos >= exon.end + 1 && pos <= exon.end + 2) {
                consequences.push_back(ConsequenceType::SPLICE_DONOR_VARIANT);
                near_splice = true;
            }
        }

        // Check splice region (3-8 bp into intron, 1-3 bp into exon)
        if (!near_splice) {
            if ((pos >= exon.start - 8 && pos <= exon.start - 3) ||
                (pos >= exon.start && pos <= exon.start + 2) ||
                (pos >= exon.end - 2 && pos <= exon.end) ||
                (pos >= exon.end + 3 && pos <= exon.end + 8)) {
                consequences.push_back(ConsequenceType::SPLICE_REGION_VARIANT);
            }
        }

        // Check introns
        if (i + 1 < transcript.exons.size()) {
            const auto& next_exon = transcript.exons[i + 1];
            if (pos > exon.end && pos < next_exon.start) {
                in_intron = true;
            }
        }
    }

    // Non-coding transcript
    if (!transcript.is_coding()) {
        if (in_exon) {
            consequences.push_back(ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT);
        } else if (in_intron) {
            consequences.push_back(ConsequenceType::INTRON_VARIANT);
        } else {
            consequences.push_back(ConsequenceType::NON_CODING_TRANSCRIPT_VARIANT);
        }
        return consequences;
    }

    // For coding transcripts
    if (in_intron && consequences.empty()) {
        consequences.push_back(ConsequenceType::INTRON_VARIANT);
        return consequences;
    }

    // Check UTR regions
    if (in_exon) {
        bool in_5utr = false;
        bool in_3utr = false;

        if (transcript.strand == '+') {
            if (pos < transcript.cds_start) in_5utr = true;
            if (pos > transcript.cds_end) in_3utr = true;
        } else {
            if (pos > transcript.cds_end) in_5utr = true;
            if (pos < transcript.cds_start) in_3utr = true;
        }

        if (in_5utr) {
            consequences.push_back(ConsequenceType::FIVE_PRIME_UTR_VARIANT);
            return consequences;
        }
        if (in_3utr) {
            consequences.push_back(ConsequenceType::THREE_PRIME_UTR_VARIANT);
            return consequences;
        }
    }

    // Coding region variant - determine specific consequence
    if (in_exon && pos >= transcript.cds_start && pos <= transcript.cds_end) {
        int cds_pos = calculate_cds_position(pos, transcript);

        if (cds_pos > 0) {
            // Check for frameshift
            int ref_len = ref.length();
            int alt_len = alt.length();
            int length_diff = alt_len - ref_len;

            if (length_diff != 0 && length_diff % 3 != 0) {
                consequences.push_back(ConsequenceType::FRAMESHIFT_VARIANT);
                return consequences;
            }

            // Inframe indel
            if (length_diff > 0 && length_diff % 3 == 0) {
                consequences.push_back(ConsequenceType::INFRAME_INSERTION);
            } else if (length_diff < 0 && length_diff % 3 == 0) {
                consequences.push_back(ConsequenceType::INFRAME_DELETION);
            }

            // SNV - check for missense/synonymous/stop
            if (ref.length() == 1 && alt.length() == 1) {
                auto [ref_codon, alt_codon] = get_affected_codons(cds_pos, ref, alt, transcript);

                if (!ref_codon.empty() && !alt_codon.empty()) {
                    char ref_aa = CodonTable::translate(ref_codon);
                    char alt_aa = CodonTable::translate(alt_codon);

                    // Check for start codon
                    if (cds_pos <= 3 && CodonTable::is_start_codon(ref_codon)) {
                        if (!CodonTable::is_start_codon(alt_codon)) {
                            consequences.push_back(ConsequenceType::START_LOST);
                            return consequences;
                        }
                    }

                    // Check for stop gained
                    if (alt_aa == '*' && ref_aa != '*') {
                        consequences.push_back(ConsequenceType::STOP_GAINED);
                        return consequences;
                    }

                    // Check for stop lost
                    if (ref_aa == '*' && alt_aa != '*') {
                        consequences.push_back(ConsequenceType::STOP_LOST);
                        return consequences;
                    }

                    // Synonymous or missense
                    if (ref_aa == alt_aa) {
                        consequences.push_back(ConsequenceType::SYNONYMOUS_VARIANT);
                    } else {
                        consequences.push_back(ConsequenceType::MISSENSE_VARIANT);
                    }
                } else {
                    consequences.push_back(ConsequenceType::CODING_SEQUENCE_VARIANT);
                }
            }
        }
    }

    if (consequences.empty()) {
        consequences.push_back(ConsequenceType::CODING_SEQUENCE_VARIANT);
    }

    return consequences;
}

int VEPAnnotator::calculate_cds_position(int genomic_pos, const Transcript& transcript) const {
    if (!transcript.is_coding()) return 0;
    if (genomic_pos < transcript.cds_start || genomic_pos > transcript.cds_end) return 0;

    int cds_pos = 0;

    if (transcript.strand == '+') {
        for (const auto& cds : transcript.cds_regions) {
            if (genomic_pos > cds.end) {
                cds_pos += cds.end - cds.start + 1;
            } else if (genomic_pos >= cds.start) {
                cds_pos += genomic_pos - cds.start + 1;
                break;
            }
        }
    } else {
        // Reverse strand - iterate backwards
        for (auto it = transcript.cds_regions.rbegin(); it != transcript.cds_regions.rend(); ++it) {
            if (genomic_pos < it->start) {
                cds_pos += it->end - it->start + 1;
            } else if (genomic_pos <= it->end) {
                cds_pos += it->end - genomic_pos + 1;
                break;
            }
        }
    }

    return cds_pos;
}

std::pair<std::string, std::string> VEPAnnotator::get_affected_codons(
    int cds_pos,
    const std::string& ref,
    const std::string& alt,
    const Transcript& transcript) const {

    if (cds_pos <= 0) return {"", ""};

    // Calculate codon position (0-based codon number)
    int codon_num = (cds_pos - 1) / 3;
    int pos_in_codon = (cds_pos - 1) % 3;

    // Get the CDS sequence
    std::string cds_seq;
    if (transcript.strand == '+') {
        for (const auto& cds : transcript.cds_regions) {
            cds_seq += reference_->get_sequence(transcript.chromosome, cds.start, cds.end);
        }
    } else {
        // Reverse strand - need to reverse complement
        for (auto it = transcript.cds_regions.rbegin(); it != transcript.cds_regions.rend(); ++it) {
            std::string seq = reference_->get_sequence(transcript.chromosome, it->start, it->end);
            // Reverse complement
            std::reverse(seq.begin(), seq.end());
            for (char& c : seq) {
                switch (c) {
                    case 'A': c = 'T'; break;
                    case 'T': c = 'A'; break;
                    case 'G': c = 'C'; break;
                    case 'C': c = 'G'; break;
                }
            }
            cds_seq += seq;
        }
    }

    if (cds_seq.empty()) return {"", ""};

    // Get the reference codon
    int codon_start = codon_num * 3;
    if (codon_start + 3 > static_cast<int>(cds_seq.length())) {
        return {"", ""};
    }

    std::string ref_codon = cds_seq.substr(codon_start, 3);

    // Create the alternate codon
    std::string alt_codon = ref_codon;
    if (ref.length() == 1 && alt.length() == 1) {
        char alt_base = alt[0];
        // For reverse strand, complement the alternate base
        if (transcript.strand == '-') {
            switch (alt_base) {
                case 'A': alt_base = 'T'; break;
                case 'T': alt_base = 'A'; break;
                case 'G': alt_base = 'C'; break;
                case 'C': alt_base = 'G'; break;
            }
        }
        alt_codon[pos_in_codon] = alt_base;
    }

    return {ref_codon, alt_codon};
}

std::string VEPAnnotator::generate_hgvsc(
    int cds_pos,
    const std::string& ref,
    const std::string& alt,
    const Transcript& transcript) const {

    std::ostringstream oss;
    oss << transcript.id << ":c.";

    if (ref.length() == 1 && alt.length() == 1) {
        // SNV
        oss << cds_pos << ref << ">" << alt;
    } else if (ref.length() > alt.length()) {
        // Deletion
        int del_len = ref.length() - alt.length();
        if (del_len == 1) {
            oss << (cds_pos + alt.length()) << "del";
        } else {
            oss << (cds_pos + alt.length()) << "_" << (cds_pos + ref.length() - 1) << "del";
        }
    } else {
        // Insertion
        oss << cds_pos << "_" << (cds_pos + 1) << "ins" << alt.substr(ref.length());
    }

    return oss.str();
}

std::string VEPAnnotator::generate_hgvsp(
    const std::string& ref_aa,
    const std::string& alt_aa,
    int protein_pos,
    const Transcript& transcript) const {

    if (ref_aa.empty() || alt_aa.empty()) return "";

    std::ostringstream oss;
    // Use transcript ID for now (ideally would use protein ID)
    oss << transcript.id << ":p.";

    if (ref_aa == alt_aa) {
        oss << ref_aa << protein_pos << "=";
    } else if (alt_aa == "Ter") {
        oss << ref_aa << protein_pos << "Ter";
    } else {
        oss << ref_aa << protein_pos << alt_aa;
    }

    return oss.str();
}

std::string VEPAnnotator::get_stats() const {
    std::ostringstream oss;
    oss << "Transcripts: " << transcript_db_->transcript_count() << "\n"
        << "Genes: " << transcript_db_->gene_count();

    auto sources = vcf_annotations_->get_sources();
    if (!sources.empty()) {
        oss << "\n" << vcf_annotations_->get_stats();
    }

    // Add annotation sources stats
    if (source_manager_) {
        auto ann_sources = source_manager_->get_sources();
        if (!ann_sources.empty()) {
            oss << "\n" << source_manager_->get_stats();
        }
    }

    return oss.str();
}

// ============================================================================
// Extended Annotation Source System
// ============================================================================

void VEPAnnotator::add_source(std::shared_ptr<AnnotationSource> source) {
    if (source_manager_) {
        source_manager_->add_source(std::move(source));
    }
}

std::vector<std::shared_ptr<AnnotationSource>> VEPAnnotator::get_sources() const {
    if (source_manager_) {
        return source_manager_->get_sources();
    }
    return {};
}

void VEPAnnotator::set_source_enabled(const std::string& name, bool enabled) {
    if (source_manager_) {
        source_manager_->set_enabled(name, enabled);
    }
}

std::vector<std::pair<std::string, std::string>> VEPAnnotator::get_available_fields() const {
    if (source_manager_) {
        return source_manager_->get_all_fields();
    }
    return {};
}

void VEPAnnotator::initialize_sources() {
    if (source_manager_) {
        source_manager_->initialize_all();
    }
}

// ============================================================================
// VCF file annotation
// ============================================================================

void annotate_vcf_file(
    const std::string& vcf_input,
    const std::string& output_path,
    const std::string& gtf_path,
    const std::string& fasta_path,
    const std::vector<std::tuple<std::string, std::string, std::string, bool>>& annotation_vcfs) {

    log(LogLevel::INFO, "Annotating VCF file: " + vcf_input);

    VEPAnnotator annotator(gtf_path, fasta_path);

    // Add custom annotation sources
    std::vector<std::string> custom_columns;
    for (const auto& [name, vcf_path, info_fields, use_tabix] : annotation_vcfs) {
        VCFAnnotationConfig config;
        config.name = name;
        config.vcf_path = vcf_path;
        config.use_tabix = use_tabix;

        // Parse info fields
        if (!info_fields.empty()) {
            std::istringstream iss(info_fields);
            std::string field;
            while (std::getline(iss, field, ',')) {
                field.erase(0, field.find_first_not_of(" \t"));
                field.erase(field.find_last_not_of(" \t") + 1);
                if (!field.empty()) {
                    config.info_fields.push_back(field);
                    custom_columns.push_back(name + ":" + field);
                }
            }
        } else {
            custom_columns.push_back(name + ":*");  // All fields
        }

        annotator.add_annotation_source(config);
    }

    // Check if input is gzipped
    bool is_gzipped = (vcf_input.length() > 3 &&
                       vcf_input.substr(vcf_input.length() - 3) == ".gz");

    gzFile gz_file = nullptr;
    std::ifstream plain_file;

    if (is_gzipped) {
        gz_file = gzopen(vcf_input.c_str(), "rb");
        if (!gz_file) {
            throw std::runtime_error("Cannot open gzipped VCF file: " + vcf_input);
        }
        log(LogLevel::INFO, "Reading gzipped VCF file");
    } else {
        plain_file.open(vcf_input);
        if (!plain_file.is_open()) {
            throw std::runtime_error("Cannot open VCF file: " + vcf_input);
        }
    }

    std::ofstream output(output_path);
    if (!output.is_open()) {
        if (gz_file) gzclose(gz_file);
        throw std::runtime_error("Cannot open output file: " + output_path);
    }

    // Write header
    output << "CHROM\tPOS\tREF\tALT\tGENE\tTRANSCRIPT\tCONSEQUENCE\tIMPACT\t"
           << "CDS_POS\tPROTEIN_POS\tAMINO_ACIDS\tCODONS\tHGVSc\tHGVSp";

    // Add custom annotation columns to header
    for (const auto& col : custom_columns) {
        output << "\t" << col;
    }
    output << "\n";

    std::string line;
    int variant_count = 0;
    char gz_buffer[65536];

    // Helper lambda to read next line
    auto read_line = [&]() -> bool {
        if (is_gzipped) {
            if (gzgets(gz_file, gz_buffer, sizeof(gz_buffer)) == nullptr) {
                return false;
            }
            line = gz_buffer;
            // Remove trailing newline
            if (!line.empty() && line.back() == '\n') {
                line.pop_back();
            }
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            return true;
        } else {
            return static_cast<bool>(std::getline(plain_file, line));
        }
    };

    while (read_line()) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string chrom, id, ref, alt, qual, filter, info;
        int pos;

        if (!(iss >> chrom >> pos >> id >> ref >> alt)) continue;

        // Handle multiple alt alleles
        std::istringstream alt_iss(alt);
        std::string single_alt;

        while (std::getline(alt_iss, single_alt, ',')) {
            auto annotation = annotator.annotate_most_severe(chrom, pos, ref, single_alt);

            output << chrom << "\t" << pos << "\t" << ref << "\t" << single_alt << "\t"
                   << annotation.gene_symbol << "\t" << annotation.transcript_id << "\t"
                   << annotation.get_consequence_string() << "\t"
                   << impact_to_string(annotation.impact) << "\t"
                   << (annotation.cds_position > 0 ? std::to_string(annotation.cds_position) : "") << "\t"
                   << (annotation.protein_position > 0 ? std::to_string(annotation.protein_position) : "") << "\t"
                   << annotation.amino_acids << "\t"
                   << annotation.codons << "\t"
                   << annotation.hgvsc << "\t"
                   << annotation.hgvsp;

            // Add custom annotation values
            for (const auto& col : custom_columns) {
                auto it = annotation.custom_annotations.find(col);
                output << "\t";
                if (it != annotation.custom_annotations.end()) {
                    output << it->second;
                }
            }
            output << "\n";

            variant_count++;
        }

        if (variant_count % 1000 == 0) {
            log(LogLevel::INFO, "Processed " + std::to_string(variant_count) + " variants...");
        }
    }

    // Cleanup
    if (gz_file) {
        gzclose(gz_file);
    }

    log(LogLevel::INFO, "Annotation complete. " + std::to_string(variant_count) + " variants written to " + output_path);
}

} // namespace vep
