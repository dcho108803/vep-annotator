/**
 * VEP Variant Annotator - Pure C++ Local Implementation
 */

#include "vep_annotator.hpp"
#include "hgvs_parser.hpp"
#include "annotation_source.hpp"
#include "file_parsers.hpp"
#include "transcript_filter.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <ctime>
#include <iomanip>
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
// Complement helper
// ============================================================================

static char complement_base(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'G': return 'C';
        case 'C': return 'G';
        default:  return base;
    }
}

static std::string complement_sequence(const std::string& seq) {
    std::string result;
    result.reserve(seq.size());
    for (char c : seq) {
        result += complement_base(c);
    }
    return result;
}

static std::string reverse_complement_sequence(const std::string& seq) {
    std::string result;
    result.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        result += complement_base(*it);
    }
    return result;
}

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
        case ConsequenceType::TRANSCRIPT_AMPLIFICATION: return "transcript_amplification";
        case ConsequenceType::FEATURE_ELONGATION: return "feature_elongation";
        case ConsequenceType::FEATURE_TRUNCATION: return "feature_truncation";
        case ConsequenceType::INFRAME_INSERTION: return "inframe_insertion";
        case ConsequenceType::INFRAME_DELETION: return "inframe_deletion";
        case ConsequenceType::MISSENSE_VARIANT: return "missense_variant";
        case ConsequenceType::PROTEIN_ALTERING_VARIANT: return "protein_altering_variant";
        case ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT: return "splice_donor_5th_base_variant";
        case ConsequenceType::SPLICE_DONOR_REGION_VARIANT: return "splice_donor_region_variant";
        case ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT: return "splice_polypyrimidine_tract_variant";
        case ConsequenceType::SPLICE_REGION_VARIANT: return "splice_region_variant";
        case ConsequenceType::INCOMPLETE_TERMINAL_CODON_VARIANT: return "incomplete_terminal_codon_variant";
        case ConsequenceType::START_RETAINED_VARIANT: return "start_retained_variant";
        case ConsequenceType::STOP_RETAINED_VARIANT: return "stop_retained_variant";
        case ConsequenceType::SYNONYMOUS_VARIANT: return "synonymous_variant";
        case ConsequenceType::CODING_SEQUENCE_VARIANT: return "coding_sequence_variant";
        case ConsequenceType::CODING_TRANSCRIPT_VARIANT: return "coding_transcript_variant";
        case ConsequenceType::MATURE_MIRNA_VARIANT: return "mature_miRNA_variant";
        case ConsequenceType::FIVE_PRIME_UTR_VARIANT: return "5_prime_UTR_variant";
        case ConsequenceType::THREE_PRIME_UTR_VARIANT: return "3_prime_UTR_variant";
        case ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT: return "non_coding_transcript_exon_variant";
        case ConsequenceType::INTRON_VARIANT: return "intron_variant";
        case ConsequenceType::NMD_TRANSCRIPT_VARIANT: return "NMD_transcript_variant";
        case ConsequenceType::NON_CODING_TRANSCRIPT_VARIANT: return "non_coding_transcript_variant";
        case ConsequenceType::UPSTREAM_GENE_VARIANT: return "upstream_gene_variant";
        case ConsequenceType::DOWNSTREAM_GENE_VARIANT: return "downstream_gene_variant";
        case ConsequenceType::TFBS_ABLATION: return "TFBS_ablation";
        case ConsequenceType::TFBS_AMPLIFICATION: return "TFBS_amplification";
        case ConsequenceType::TF_BINDING_SITE_VARIANT: return "TF_binding_site_variant";
        case ConsequenceType::REGULATORY_REGION_ABLATION: return "regulatory_region_ablation";
        case ConsequenceType::REGULATORY_REGION_AMPLIFICATION: return "regulatory_region_amplification";
        case ConsequenceType::REGULATORY_REGION_VARIANT: return "regulatory_region_variant";
        case ConsequenceType::INTERGENIC_VARIANT: return "intergenic_variant";
        case ConsequenceType::SEQUENCE_VARIANT: return "sequence_variant";
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
        case ConsequenceType::TRANSCRIPT_AMPLIFICATION:
        case ConsequenceType::FEATURE_ELONGATION:
        case ConsequenceType::FEATURE_TRUNCATION:
            return Impact::HIGH;

        case ConsequenceType::INFRAME_INSERTION:
        case ConsequenceType::INFRAME_DELETION:
        case ConsequenceType::MISSENSE_VARIANT:
        case ConsequenceType::PROTEIN_ALTERING_VARIANT:
            return Impact::MODERATE;

        case ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT:
        case ConsequenceType::SPLICE_DONOR_REGION_VARIANT:
        case ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT:
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
        {'U', "Sec"}, {'*', "Ter"}, {'X', "Xaa"}
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

bool CodonTable::is_stop_codon(const std::string& codon, const std::string& chromosome) {
    std::string upper = codon;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

    // Normalize chromosome name for MT detection
    std::string chr = chromosome;
    if (chr.substr(0, 3) == "chr") chr = chr.substr(3);

    if (chr == "M" || chr == "MT") {
        // Vertebrate mitochondrial code: AGA and AGG are stops; TGA is NOT a stop (codes for Trp)
        return upper == "TAA" || upper == "TAG" || upper == "AGA" || upper == "AGG";
    }
    return upper == "TAA" || upper == "TAG" || upper == "TGA";
}

char CodonTable::translate_mt(const std::string& codon) {
    // Vertebrate mitochondrial code differences from standard:
    // AGA -> Ter (not Arg), AGG -> Ter (not Arg)
    // ATA -> Met (not Ile), TGA -> Trp (not Ter)
    static const std::unordered_map<std::string, char> mt_overrides = {
        {"AGA", '*'}, {"AGG", '*'}, {"ATA", 'M'}, {"TGA", 'W'}
    };

    if (codon.length() != 3) return 'X';
    std::string upper = codon;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

    auto it = mt_overrides.find(upper);
    if (it != mt_overrides.end()) return it->second;
    return translate(codon);
}

char CodonTable::translate(const std::string& codon, const std::string& chromosome) {
    std::string norm = chromosome;
    if (norm.length() > 3 && norm.substr(0, 3) == "chr") {
        norm = norm.substr(3);
    }
    if (norm == "MT" || norm == "M" || norm == "chrM") {
        return translate_mt(codon);
    }
    return translate(codon);
}

// Known human selenoprotein genes where TGA encodes selenocysteine (U)
static const std::set<std::string> SELENOPROTEIN_GENES = {
    "DIO1", "DIO2", "DIO3",            // Iodothyronine deiodinases
    "GPX1", "GPX2", "GPX3", "GPX4", "GPX6", // Glutathione peroxidases
    "MSRB1",                            // Methionine sulfoxide reductase
    "SELENOF", "SELENOH", "SELENOI", "SELENOK", "SELENOM",
    "SELENON", "SELENOO", "SELENOP", "SELENOS", "SELENOT",
    "SELENOV", "SELENOW",
    "SEPHS2",                           // Selenophosphate synthetase
    "TXNRD1", "TXNRD2", "TXNRD3"       // Thioredoxin reductases
};

/**
 * Check if a gene is a known selenoprotein gene
 */
static bool is_selenoprotein_gene(const std::string& gene_name) {
    return SELENOPROTEIN_GENES.count(gene_name) > 0;
}

/**
 * Translate codon with selenocysteine awareness
 * If gene is a selenoprotein, TGA -> U (selenocysteine) instead of * (stop)
 */
static char translate_with_sec(const std::string& codon, const std::string& chromosome,
                                const std::string& gene_name) {
    char aa = CodonTable::translate(codon, chromosome);
    if (aa == '*' && !gene_name.empty()) {
        std::string upper = codon;
        std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
        if (upper == "TGA" && is_selenoprotein_gene(gene_name)) {
            return 'U';  // Selenocysteine
        }
    }
    return aa;
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
    std::string result;
    for (const auto& c : consequences) {
        if (!result.empty()) result += '&';
        result += consequence_to_string(c);
    }
    return result;
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
    if (!hgvsg.empty()) {
        summary["hgvsg"] = hgvsg;
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
    // Also extracts VCF ID column and builds allele_string for co-located variant output
    auto process_records = [&](const std::string& source_name, const VCFAnnotationConfig& config,
                               const std::vector<VCFRecord>& records) {
        std::string collected_ids;  // Comma-separated IDs from all matched records
        for (const auto& record : records) {
            // Check allele match if required
            bool matched = false;
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
                matched = true;

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
                matched = true;
                for (const auto& [field, value] : record.info) {
                    std::string key = source_name + ":" + field;
                    result[key] = value;
                }
            }

            // Collect VCF ID column from matched records
            if (matched && !record.id.empty() && record.id != ".") {
                if (!collected_ids.empty()) collected_ids += ",";
                collected_ids += record.id;
            }
        }
        // Store collected IDs (comma-separated if multiple co-located variants)
        if (!collected_ids.empty()) {
            result[source_name + ":ID"] = collected_ids;
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
            auto rc_it = pimpl_->record_counts.find(config.name);
            oss << (rc_it != pimpl_->record_counts.end() ? rc_it->second : 0) << " records";
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
    std::string normalized_chrom = chrom;
    if (normalized_chrom.size() >= 4 && normalized_chrom[0] == 'c' && normalized_chrom[1] == 'h' && normalized_chrom[2] == 'r') {
        normalized_chrom = normalized_chrom.substr(3);
    }
    auto it = pimpl_->sequences.find(normalized_chrom);
    if (it == pimpl_->sequences.end()) return 'N';
    int idx = pos - 1; // Convert to 0-based
    if (idx < 0 || idx >= static_cast<int>(it->second.size())) return 'N';
    return it->second[idx];
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

    // Gene spatial index: chromosome -> sorted vector of (gene.start, gene_ptr)
    std::unordered_map<std::string, std::vector<std::pair<int, const Gene*>>> gene_chrom_index;

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

        // Build gene spatial index (genes map is fully populated at this point)
        for (auto& [gene_id, gene] : genes) {
            gene_chrom_index[gene.chromosome].emplace_back(gene.start, &gene);
        }

        // Sort gene index by start position
        for (auto& [chrom, entries] : gene_chrom_index) {
            std::sort(entries.begin(), entries.end());
        }
    }
};

// Helper to parse GTF attribute string
static std::unordered_map<std::string, std::string> parse_gtf_attributes(const std::string& attr_str) {
    std::unordered_map<std::string, std::string> attrs;
    size_t pos = 0;
    const size_t len = attr_str.size();

    while (pos < len) {
        // Skip whitespace and semicolons
        while (pos < len && (attr_str[pos] == ' ' || attr_str[pos] == ';' || attr_str[pos] == '\t'))
            ++pos;
        if (pos >= len) break;

        // Read key (word characters)
        size_t key_start = pos;
        while (pos < len && attr_str[pos] != ' ' && attr_str[pos] != '\t' && attr_str[pos] != '"')
            ++pos;
        std::string key = attr_str.substr(key_start, pos - key_start);

        // Skip to opening quote
        while (pos < len && attr_str[pos] != '"') ++pos;
        if (pos >= len) break;
        ++pos; // skip opening quote

        // Read value until closing quote
        size_t val_start = pos;
        while (pos < len && attr_str[pos] != '"') ++pos;
        attrs[std::move(key)] = attr_str.substr(val_start, pos - val_start);
        if (pos < len) ++pos; // skip closing quote
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
            // Parse gene source and HGNC ID from GTF attributes
            if (attrs.count("gene_source")) gene.source = attrs["gene_source"];
            if (attrs.count("hgnc_id")) gene.hgnc_id = attrs["hgnc_id"];
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

            // Parse tags (GTF "tag" attribute may contain multiple semicolon- or
            // comma-separated values: "basic", "Ensembl_canonical", "MANE_Select", etc.)
            if (attrs.count("tag")) {
                const std::string& tags = attrs["tag"];
                // Helper: check if a tag value appears as a complete token
                // (bounded by start/end of string, comma, semicolon, or space)
                auto has_tag = [&tags](const std::string& tag) -> bool {
                    size_t pos = 0;
                    while ((pos = tags.find(tag, pos)) != std::string::npos) {
                        // Check left boundary
                        bool left_ok = (pos == 0 || tags[pos - 1] == ',' ||
                                       tags[pos - 1] == ';' || tags[pos - 1] == ' ');
                        // Check right boundary
                        size_t end = pos + tag.length();
                        bool right_ok = (end >= tags.length() || tags[end] == ',' ||
                                        tags[end] == ';' || tags[end] == ' ');
                        if (left_ok && right_ok) return true;
                        pos = end;
                    }
                    return false;
                };
                tr.is_canonical = (has_tag("Ensembl_canonical") || has_tag("canonical"));
                if (has_tag("MANE_Select")) {
                    tr.mane_select = "YES";
                    tr.is_canonical = true; // MANE implies canonical
                }
                if (has_tag("MANE_Plus_Clinical")) {
                    tr.mane_plus_clinical = "YES";
                }
                if (has_tag("basic")) {
                    tr.gencode_basic = true;
                }
                if (has_tag("cds_start_NF")) {
                    tr.cds_start_NF = true;
                }
                if (has_tag("cds_end_NF")) {
                    tr.cds_end_NF = true;
                }
            }

            // Parse APPRIS annotation
            if (attrs.count("tag")) {
                const std::string& tags = attrs["tag"];
                // APPRIS tags: appris_principal_1 through appris_alternative_2
                for (const auto& prefix : {"appris_principal", "appris_alternative"}) {
                    size_t apos = tags.find(prefix);
                    if (apos != std::string::npos) {
                        // Extract the full appris tag (e.g., "appris_principal_1")
                        size_t end = tags.find_first_of(",; ", apos);
                        tr.appris = tags.substr(apos, end == std::string::npos ?
                            std::string::npos : end - apos);
                        break;
                    }
                }
            }

            // Parse transcript support level
            if (attrs.count("transcript_support_level")) {
                try {
                    // Value may be "1 (assigned to previous version)"
                    std::string tsl_str = attrs["transcript_support_level"];
                    tr.tsl = std::stoi(tsl_str);
                } catch (...) {}
            }

            // Extract CCDS and protein IDs from GTF attributes
            if (attrs.count("ccds_id")) {
                tr.ccds_id = attrs["ccds_id"];
            } else if (attrs.count("ccdsid")) {
                tr.ccds_id = attrs["ccdsid"];
            }
            if (attrs.count("protein_id")) {
                tr.protein_id = attrs["protein_id"];
            }
            if (attrs.count("transcript_version")) {
                tr.version = attrs["transcript_version"];
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
                exon.exon_number = 0;
                if (attrs.count("exon_number")) {
                    try { exon.exon_number = std::stoi(attrs["exon_number"]); }
                    catch (...) {}
                }

                pimpl_->transcripts[transcript_id].exons.push_back(exon);
            }
        }
        else if (feature == "CDS" && !transcript_id.empty()) {
            if (pimpl_->transcripts.count(transcript_id)) {
                CDS cds;
                cds.start = start;
                cds.end = end;
                cds.phase = 0;
                if (frame != ".") {
                    try { cds.phase = std::stoi(frame); }
                    catch (...) {}
                }

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

    // Look up the gene spatial index for this chromosome
    auto chrom_it = pimpl_->gene_chrom_index.find(normalized);
    if (chrom_it == pimpl_->gene_chrom_index.end()) {
        return results;
    }

    const auto& gene_entries = chrom_it->second;

    int left_bound = pos - distance;
    int right_bound = pos + distance;

    // Find first gene with start >= left_bound via binary search
    auto it = std::lower_bound(
        gene_entries.begin(), gene_entries.end(),
        std::make_pair(left_bound, static_cast<const Gene*>(nullptr)),
        [](const std::pair<int, const Gene*>& entry, const std::pair<int, const Gene*>& val) {
            return entry.first < val.first;
        }
    );

    // Scan backwards from `it` to find genes that start before left_bound
    // but whose end extends into range (i.e., gene.end >= pos - distance)
    if (it != gene_entries.begin()) {
        auto back_it = it;
        do {
            --back_it;
            const Gene* gene = back_it->second;
            // Gene starts before left_bound. Check if it's within distance of pos.
            int dist = 0;
            if (pos < gene->start) {
                dist = gene->start - pos;
            } else if (pos > gene->end) {
                dist = pos - gene->end;
            }
            if (dist <= distance) {
                results.push_back(gene);
            } else {
                // This gene's end is too far left. Since genes are sorted by start,
                // all further-left genes will also be too far, so stop.
                break;
            }
        } while (back_it != gene_entries.begin());
    }

    // Scan forward from `it` while gene.start <= right_bound
    for (auto fwd_it = it; fwd_it != gene_entries.end(); ++fwd_it) {
        if (fwd_it->first > right_bound) {
            break;
        }
        const Gene* gene = fwd_it->second;
        // Gene starts within [left_bound, right_bound]. Check actual distance.
        int dist = 0;
        if (pos < gene->start) {
            dist = gene->start - pos;
        } else if (pos > gene->end) {
            dist = pos - gene->end;
        }
        if (dist <= distance) {
            results.push_back(gene);
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
    const std::string& ref_in,
    const std::string& alt_in) {

    // Normalize dash-encoded empty alleles (from --minimal mode)
    std::string ref = (ref_in == "-") ? "" : ref_in;
    std::string alt = (alt_in == "-") ? "" : alt_in;

    std::vector<VariantAnnotation> results;

    // Get overlapping transcripts
    int var_end = pos + std::max(0, static_cast<int>(ref.length()) - 1);
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
        ann.feature_type = "Intergenic";

        // Find nearby genes for upstream/downstream annotation
        int max_dist = std::max(upstream_distance_, downstream_distance_);
        auto nearby = transcript_db_->get_nearby_genes(chrom, pos, max_dist);
        if (!nearby.empty()) {
            // Find the actual closest gene by distance
            const Gene* closest = nearby[0];
            int closest_dist = (var_end < closest->start) ? (closest->start - var_end) : (pos > closest->end ? pos - closest->end : 0);
            for (size_t gi = 1; gi < nearby.size(); ++gi) {
                int d = (var_end < nearby[gi]->start) ? (nearby[gi]->start - var_end) : (pos > nearby[gi]->end ? pos - nearby[gi]->end : 0);
                if (d < closest_dist) {
                    closest = nearby[gi];
                    closest_dist = d;
                }
            }
            ann.gene_symbol = closest->name;
            ann.gene_id = closest->id;

            // Determine direction and apply per-direction distance thresholds
            bool is_upstream_of_gene = (pos < closest->start);
            ConsequenceType dir_csq;

            if (is_upstream_of_gene) {
                dir_csq = (closest->strand == '+') ?
                    ConsequenceType::UPSTREAM_GENE_VARIANT :
                    ConsequenceType::DOWNSTREAM_GENE_VARIANT;
            } else {
                dir_csq = (closest->strand == '+') ?
                    ConsequenceType::DOWNSTREAM_GENE_VARIANT :
                    ConsequenceType::UPSTREAM_GENE_VARIANT;
            }

            // Check per-direction distance limit
            int actual_dist = is_upstream_of_gene ?
                (closest->start - var_end) : (pos - closest->end);
            int limit = (dir_csq == ConsequenceType::UPSTREAM_GENE_VARIANT) ?
                upstream_distance_ : downstream_distance_;

            if (actual_dist <= limit) {
                ann.consequences[0] = dir_csq;
                ann.distance = actual_dist;
                ann.custom_annotations["DISTANCE"] = std::to_string(actual_dist);
            }
            // else remains INTERGENIC_VARIANT

            // Store nearest gene for --nearest output
            ann.custom_annotations["NEAREST"] = closest->name;
        }

        results.push_back(std::move(ann));
    } else {
        for (const auto* transcript : transcripts) {
            auto ann = annotate_transcript(chrom, pos, ref, alt, *transcript);
            results.push_back(std::move(ann));
        }

        // Also check for upstream/downstream of non-overlapping nearby transcripts
        int max_dist = std::max(upstream_distance_, downstream_distance_);
        auto wider_transcripts = transcript_db_->get_transcripts_in_region(
            chrom, pos - max_dist, var_end + max_dist);

        for (const auto* transcript : wider_transcripts) {
            // Skip transcripts that already overlap the variant (already annotated above)
            if (transcript->start <= var_end && transcript->end >= pos) continue;

            // This transcript is nearby but non-overlapping - generate upstream/downstream
            bool is_upstream_of_gene = (pos < transcript->start);
            ConsequenceType dir_csq;

            if (is_upstream_of_gene) {
                dir_csq = (transcript->strand == '+') ?
                    ConsequenceType::UPSTREAM_GENE_VARIANT :
                    ConsequenceType::DOWNSTREAM_GENE_VARIANT;
            } else {
                dir_csq = (transcript->strand == '+') ?
                    ConsequenceType::DOWNSTREAM_GENE_VARIANT :
                    ConsequenceType::UPSTREAM_GENE_VARIANT;
            }

            int actual_dist = is_upstream_of_gene ?
                (transcript->start - var_end) : (pos - transcript->end);
            int limit = (dir_csq == ConsequenceType::UPSTREAM_GENE_VARIANT) ?
                upstream_distance_ : downstream_distance_;

            if (actual_dist > 0 && actual_dist <= limit) {
                VariantAnnotation ann;
                ann.input_variant = chrom + "-" + std::to_string(pos) + "-" + ref + "-" + alt;
                ann.chromosome = chrom;
                ann.position = pos;
                ann.ref_allele = ref;
                ann.alt_allele = alt;
                ann.gene_symbol = transcript->gene_name;
                ann.gene_id = transcript->gene_id;
                ann.transcript_id = transcript->id;
                ann.biotype = transcript->biotype;
                ann.is_canonical = transcript->is_canonical;
                ann.consequences.push_back(dir_csq);
                ann.impact = Impact::MODIFIER;
                ann.strand = transcript->strand;
                ann.distance = actual_dist;
                ann.custom_annotations["DISTANCE"] = std::to_string(actual_dist);

                // Set transcript source
                if (transcript->id.substr(0, 4) == "ENST") {
                    ann.source = "Ensembl";
                } else if (transcript->id.substr(0, 2) == "NM" || transcript->id.substr(0, 2) == "NR" ||
                           transcript->id.substr(0, 2) == "XM" || transcript->id.substr(0, 2) == "XR") {
                    ann.source = "RefSeq";
                }
                results.push_back(std::move(ann));
            }
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

            // Inject consequence string so annotation sources can check it
            {
                std::ostringstream csq_oss;
                bool first_csq = true;
                for (const auto& c : ann.consequences) {
                    if (!first_csq) csq_oss << "&";
                    csq_oss << consequence_to_string(c);
                    first_csq = false;
                }
                ann.custom_annotations["_consequences"] = csq_oss.str();
            }

            // Call all enabled annotation sources
            source_manager_->annotate_all(chrom, pos, ref, alt, tx_ptr, ann.custom_annotations);

            // Inject regulatory consequences based on source annotations
            append_regulatory_consequences(ann);
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

    // Sort by consequence severity using Perl VEP severity ranking
    std::sort(all_annotations.begin(), all_annotations.end(),
              [](const VariantAnnotation& a, const VariantAnnotation& b) {
                  auto get_min_rank = [](const VariantAnnotation& ann) -> int {
                      int min_val = 999;
                      for (const auto& c : ann.consequences) {
                          int v = get_consequence_rank(c);
                          if (v < min_val) min_val = v;
                      }
                      return min_val;
                  };
                  return get_min_rank(a) < get_min_rank(b);
              });

    return all_annotations[0];
}

// Forward declarations for static helper functions used by annotate_transcript
static std::string calculate_intronic_hgvsc_position(
    int genomic_pos,
    const Transcript& transcript,
    const std::function<int(int, const Transcript&)>& calc_cds_pos);

static std::string generate_hgvsc_intronic(
    const std::string& intronic_pos,
    const std::string& ref,
    const std::string& alt,
    const Transcript& transcript,
    bool include_version = false);

static std::string generate_hgvsc_utr(
    const std::string& utr_pos_str,
    const std::string& ref,
    const std::string& alt,
    const Transcript& transcript,
    bool include_version = false);

VariantAnnotation VEPAnnotator::annotate_transcript(
    const std::string& chrom,
    int pos,
    const std::string& ref,
    const std::string& alt,
    const Transcript& transcript) {

    // Compute 3' shifted coordinates for genomic output and HGVS
    int shifted_pos = pos;
    std::string shifted_ref = ref;
    std::string shifted_alt = alt;

    if ((shift_3prime_ || shift_genomic_) && ref.length() != alt.length()) {
        auto [sp, sr, sa] = right_normalize(chrom, pos, ref, alt, transcript);
        shifted_pos = sp;
        shifted_ref = sr;
        shifted_alt = sa;
    }

    int hgvs_offset = shifted_pos - pos; // Track shift distance for HGVS_OFFSET

    VariantAnnotation ann;
    ann.input_variant = chrom + "-" + std::to_string(pos) + "-" + ref + "-" + alt;
    ann.chromosome = chrom;
    ann.position = shift_genomic_ ? shifted_pos : pos;
    ann.ref_allele = shift_genomic_ ? shifted_ref : ref;
    ann.alt_allele = shift_genomic_ ? shifted_alt : alt;
    ann.gene_symbol = transcript.gene_name;
    ann.gene_id = transcript.gene_id;
    ann.transcript_id = transcript.id;
    ann.biotype = transcript.biotype;
    ann.is_canonical = transcript.is_canonical;
    ann.strand = transcript.strand;
    ann.custom_annotations["STRAND"] = (transcript.strand == '+') ? "1" : "-1";

    // Set transcript source based on ID prefix
    if (transcript.id.substr(0, 4) == "ENST") {
        ann.source = "Ensembl";
    } else if (transcript.id.substr(0, 2) == "NM" || transcript.id.substr(0, 2) == "NR" ||
               transcript.id.substr(0, 2) == "XM" || transcript.id.substr(0, 2) == "XR") {
        ann.source = "RefSeq";
    }

    // Determine consequences (always use original coordinates per Perl VEP behavior)
    ann.consequences = determine_consequences(pos, ref, alt, transcript);

    // Add NMD_TRANSCRIPT_VARIANT for NMD-targeted transcripts
    if (transcript.biotype == "nonsense_mediated_decay") {
        ann.consequences.push_back(ConsequenceType::NMD_TRANSCRIPT_VARIANT);
    }

    // Get most severe impact
    ann.impact = Impact::MODIFIER;
    for (const auto& c : ann.consequences) {
        Impact i = get_impact(c);
        if (static_cast<int>(i) < static_cast<int>(ann.impact)) {
            ann.impact = i;
        }
    }

    // Calculate cDNA position
    ann.cdna_position = calculate_cdna_position(pos, transcript);

    // Calculate CDS position and codon changes for coding variants
    if (transcript.is_coding()) {
        const std::string cached_cds = build_cds_sequence(chrom, transcript);
        int cds_pos = calculate_cds_position(pos, transcript);
        // For minus strand with multi-base ref, calculate_cds_position(pos) returns the
        // HIGHEST affected CDS position. Adjust to the FIRST (lowest) affected position.
        if (transcript.strand == '-' && ref.length() > 1) {
            cds_pos = cds_pos - static_cast<int>(ref.length()) + 1;
            if (cds_pos < 1) cds_pos = 1;
        }
        if (cds_pos > 0) {
            annotate_coding_region(chrom, pos, ref, alt, transcript, cached_cds, cds_pos, ann);
        } else {
            annotate_noncds_hgvsc(pos, ref, alt, transcript, ann);
        }
    }

    // Populate metadata and return
    populate_transcript_metadata(chrom, pos, ref, alt, transcript, hgvs_offset, ann);

    return ann;
}

// ---------------------------------------------------------------------------
// annotate_coding_region: CDS position, codons, HGVSc, HGVSp
// ---------------------------------------------------------------------------
void VEPAnnotator::annotate_coding_region(
    const std::string& chrom, int pos,
    const std::string& ref, const std::string& alt,
    const Transcript& transcript,
    const std::string& cached_cds,
    int cds_pos,
    VariantAnnotation& ann) {
    (void)pos; // cds_pos passed directly, pos retained for potential future use

    ann.cds_position = cds_pos;
    ann.protein_position = (cds_pos - 1) / 3 + 1;

    // Get affected codons (using cached CDS)
    auto [ref_codon, alt_codon] = get_affected_codons(cds_pos, ref, alt, transcript, cached_cds);
    if (ref_codon.empty() || alt_codon.empty()) return;

    ann.codons = ref_codon + "/" + alt_codon;

    char ref_aa = translate_with_sec(ref_codon, transcript.chromosome, transcript.gene_name);
    char alt_aa = translate_with_sec(alt_codon, transcript.chromosome, transcript.gene_name);
    ann.amino_acids = std::string(1, ref_aa) + "/" + std::string(1, alt_aa);

    // Apply 3' shifting in CDS space for HGVS generation (indels only)
    int hgvs_cds_pos = cds_pos;
    std::string hgvs_ref = ref;
    std::string hgvs_alt = alt;
    if (ref.length() != alt.length()) {
        std::string shorter_allele = (ref.length() < alt.length()) ? ref : alt;
        std::string longer_allele = (ref.length() >= alt.length()) ? ref : alt;
        if (!shorter_allele.empty() && shorter_allele[0] == longer_allele[0]) {
            std::string indel_bases = longer_allele.substr(shorter_allele.length());
            if (transcript.strand == '-') {
                indel_bases = reverse_complement_sequence(indel_bases);
            }

            int shift_pos = cds_pos;
            int cds_len = static_cast<int>(cached_cds.size());
            int indel_len = static_cast<int>(indel_bases.size());
            bool is_cds_insertion = (ref.length() < alt.length());
            for (int i = 0; i < 1000; ++i) {
                int check_idx = is_cds_insertion ? shift_pos : (shift_pos + indel_len - 1);
                if (check_idx >= cds_len) break;
                if (std::toupper(cached_cds[check_idx]) == std::toupper(indel_bases[0])) {
                    std::rotate(indel_bases.begin(), indel_bases.begin() + 1, indel_bases.end());
                    shift_pos++;
                } else {
                    break;
                }
            }
            hgvs_cds_pos = shift_pos;
        }
    }

    ann.hgvsc = generate_hgvsc(hgvs_cds_pos, hgvs_ref, hgvs_alt, transcript, cached_cds);

    // Compute HGVSp extras for indels
    int hgvs_protein_end = 0;
    int hgvs_fs_ter_dist = 0;
    std::string hgvs_end_ref_aa;

    if (ref.length() != alt.length()) {
        // protein_end and end AA for multi-codon changes (deletions and delins)
        if (ref.length() > 1) {
            int cds_end_pos = cds_pos + static_cast<int>(ref.length()) - 1;
            hgvs_protein_end = (cds_end_pos - 1) / 3 + 1;
            if (hgvs_protein_end > ann.protein_position && !cached_cds.empty()) {
                int end_codon_start = (hgvs_protein_end - 1) * 3;
                if (end_codon_start + 3 <= static_cast<int>(cached_cds.length())) {
                    std::string end_codon = cached_cds.substr(end_codon_start, 3);
                    for (auto& c : end_codon) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                    char end_aa_ch = CodonTable::translate(end_codon, transcript.chromosome);
                    hgvs_end_ref_aa = CodonTable::get_three_letter(end_aa_ch);
                }
            }
        }

        // For inframe insertions: get flanking AA at protein_pos+1
        if (alt.length() > ref.length() && (alt.length() - ref.length()) % 3 == 0 && !cached_cds.empty()) {
            int next_pos = ann.protein_position;
            int next_codon_start = next_pos * 3;
            if (next_codon_start + 3 <= static_cast<int>(cached_cds.length())) {
                std::string next_codon = cached_cds.substr(next_codon_start, 3);
                for (auto& c : next_codon) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                char next_aa_ch = CodonTable::translate(next_codon, transcript.chromosome);
                hgvs_end_ref_aa = CodonTable::get_three_letter(next_aa_ch);
            }
        }

        // Frameshift fsTer distance: build mutated CDS, scan for stop
        bool is_fs = false;
        for (const auto& c : ann.consequences) {
            if (c == ConsequenceType::FRAMESHIFT_VARIANT) { is_fs = true; break; }
        }
        if (is_fs && !cached_cds.empty()) {
            std::string mut_cds;
            int ref_len = static_cast<int>(ref.length());
            if (transcript.strand == '+') {
                int cds_offset = cds_pos - 1;
                if (cds_offset + ref_len <= static_cast<int>(cached_cds.length()))
                    mut_cds = cached_cds.substr(0, cds_offset) + alt + cached_cds.substr(cds_offset + ref_len);
            } else {
                std::string rc_alt = reverse_complement_sequence(alt);
                int cds_offset = cds_pos - 1;
                if (cds_offset >= 0 && cds_offset + ref_len <= static_cast<int>(cached_cds.length()))
                    mut_cds = cached_cds.substr(0, cds_offset) + rc_alt + cached_cds.substr(cds_offset + ref_len);
            }

            if (!mut_cds.empty()) {
                int scan_start = (ann.protein_position - 1) * 3;
                int dist = 0;
                for (int i = scan_start; i + 2 < static_cast<int>(mut_cds.length()); i += 3) {
                    dist++;
                    std::string codon = mut_cds.substr(i, 3);
                    for (auto& c : codon) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                    if (CodonTable::is_stop_codon(codon, transcript.chromosome)) {
                        hgvs_fs_ter_dist = dist;
                        break;
                    }
                }
            }
        }
    }

    ann.hgvsp = generate_hgvsp(
        CodonTable::get_three_letter(ref_aa),
        CodonTable::get_three_letter(alt_aa),
        ann.protein_position,
        transcript,
        ann.consequences,
        hgvs_protein_end,
        hgvs_fs_ter_dist,
        hgvs_end_ref_aa
    );
}

// ---------------------------------------------------------------------------
// annotate_noncds_hgvsc: intronic and UTR HGVSc notation
// ---------------------------------------------------------------------------
void VEPAnnotator::annotate_noncds_hgvsc(
    int pos,
    const std::string& ref, const std::string& alt,
    const Transcript& transcript,
    VariantAnnotation& ann) {

    // Check if variant is intronic (between exons)
    bool is_intronic = false;
    std::vector<Exon> sorted_exons = transcript.exons;
    std::sort(sorted_exons.begin(), sorted_exons.end(),
              [](const Exon& a, const Exon& b) { return a.start < b.start; });
    for (size_t i = 0; i + 1 < sorted_exons.size(); ++i) {
        if (pos > sorted_exons[i].end && pos < sorted_exons[i + 1].start) {
            is_intronic = true;
            break;
        }
    }

    if (is_intronic) {
        auto calc_cds_fn = [this](int gpos, const Transcript& t) -> int {
            return this->calculate_cds_position(gpos, t);
        };
        std::string intronic_pos = calculate_intronic_hgvsc_position(
            pos, transcript, calc_cds_fn);
        if (!intronic_pos.empty()) {
            ann.hgvsc = generate_hgvsc_intronic(intronic_pos, ref, alt, transcript, transcript_version_);
        }
        // Splice donor/acceptor variants in coding transcripts get p.?
        for (const auto& c : ann.consequences) {
            if (c == ConsequenceType::SPLICE_DONOR_VARIANT ||
                c == ConsequenceType::SPLICE_ACCEPTOR_VARIANT) {
                std::string protein_ref;
                if (!transcript.protein_id.empty()) {
                    protein_ref = transcript.protein_id;
                } else {
                    protein_ref = transcript.id;
                    if (transcript_version_ && !transcript.version.empty()) {
                        protein_ref += "." + transcript.version;
                    }
                }
                ann.hgvsp = protein_ref + ":p.?";
                break;
            }
        }
        return;
    }

    // Check for UTR variants
    bool in_5utr = false;
    bool in_3utr = false;

    if (transcript.strand == '+') {
        if (pos < transcript.cds_start && pos >= transcript.start) {
            in_5utr = true;
        } else if (pos > transcript.cds_end && pos <= transcript.end) {
            in_3utr = true;
        }
    } else {
        if (pos > transcript.cds_end && pos <= transcript.end) {
            in_5utr = true;
        } else if (pos < transcript.cds_start && pos >= transcript.start) {
            in_3utr = true;
        }
    }

    if (in_5utr) {
        int distance = 0;
        if (transcript.strand == '+') {
            for (const auto& exon : sorted_exons) {
                if (exon.end < pos) continue;
                if (exon.start >= transcript.cds_start) break;
                int ex_start = std::max(exon.start, pos);
                int ex_end = std::min(exon.end, transcript.cds_start - 1);
                if (ex_start <= ex_end) {
                    distance += ex_end - ex_start + 1;
                }
            }
        } else {
            for (auto it = sorted_exons.rbegin(); it != sorted_exons.rend(); ++it) {
                if (it->start > pos) continue;
                if (it->end <= transcript.cds_end) break;
                int ex_start = std::max(it->start, transcript.cds_end + 1);
                int ex_end = std::min(it->end, pos);
                if (ex_start <= ex_end) {
                    distance += ex_end - ex_start + 1;
                }
            }
        }
        if (distance > 0) {
            std::string utr_pos_str = "-" + std::to_string(distance);
            ann.hgvsc = generate_hgvsc_utr(utr_pos_str, ref, alt, transcript, transcript_version_);
        }
    } else if (in_3utr) {
        int distance = 0;
        if (transcript.strand == '+') {
            for (const auto& exon : sorted_exons) {
                if (exon.end <= transcript.cds_end) continue;
                if (exon.start > pos) break;
                int ex_start = std::max(exon.start, transcript.cds_end + 1);
                int ex_end = std::min(exon.end, pos);
                if (ex_start <= ex_end) {
                    distance += ex_end - ex_start + 1;
                }
            }
        } else {
            for (auto it = sorted_exons.rbegin(); it != sorted_exons.rend(); ++it) {
                if (it->start >= transcript.cds_start) continue;
                if (it->end < pos) break;  // past the variant going downward
                int ex_start = std::max(it->start, pos);
                int ex_end = std::min(it->end, transcript.cds_start - 1);
                if (ex_start <= ex_end) {
                    distance += ex_end - ex_start + 1;
                }
            }
        }
        if (distance > 0) {
            std::string utr_pos_str = "*" + std::to_string(distance);
            ann.hgvsc = generate_hgvsc_utr(utr_pos_str, ref, alt, transcript, transcript_version_);
        }
    }
}

// ---------------------------------------------------------------------------
// populate_transcript_metadata: HGVSg, CCDS, ENSP, CANONICAL, MANE, etc.
// ---------------------------------------------------------------------------
void VEPAnnotator::populate_transcript_metadata(
    const std::string& chrom, int pos,
    const std::string& ref, const std::string& alt,
    const Transcript& transcript,
    int hgvs_offset,
    VariantAnnotation& ann) {

    // Generate HGVSg (genomic-level, same for all transcripts)
    ann.hgvsg = generate_hgvsg(chrom, pos, ref, alt);

    // Add CCDS and protein IDs if available
    if (!transcript.ccds_id.empty()) {
        ann.custom_annotations["CCDS"] = transcript.ccds_id;
    }
    if (!transcript.protein_id.empty()) {
        ann.custom_annotations["ENSP"] = transcript.protein_id;
    }

    // Populate display annotations from transcript metadata
    if (ann.is_canonical) {
        ann.custom_annotations["CANONICAL"] = "YES";
    }
    if (!transcript.mane_select.empty()) {
        ann.custom_annotations["MANE_SELECT"] = transcript.mane_select;
    }
    if (!transcript.mane_plus_clinical.empty()) {
        ann.custom_annotations["MANE_PLUS_CLINICAL"] = transcript.mane_plus_clinical;
    }
    if (transcript.tsl > 0) {
        ann.custom_annotations["TSL"] = std::to_string(transcript.tsl);
    }
    if (!transcript.appris.empty()) {
        ann.custom_annotations["APPRIS"] = transcript.appris;
    }
    ann.custom_annotations["BIOTYPE"] = transcript.biotype;
    if (transcript.gencode_basic) {
        ann.custom_annotations["GENCODE_BASIC"] = "YES";
    }

    // Compute transcript and CDS lengths if requested
    if (include_total_length_) {
        int transcript_length = 0;
        for (const auto& exon : transcript.exons) {
            transcript_length += exon.end - exon.start + 1;
        }
        ann.custom_annotations["TRANSCRIPT_LENGTH"] = std::to_string(transcript_length);

        if (transcript.is_coding()) {
            int cds_length = transcript.get_cds_length();
            ann.custom_annotations["CDS_LENGTH"] = std::to_string(cds_length);
        }
    }

    // Set total exon/intron counts
    ann.total_exons = static_cast<int>(transcript.exons.size());
    ann.total_introns = ann.total_exons > 1 ? ann.total_exons - 1 : 0;

    // Append transcript version if enabled
    if (transcript_version_ && !transcript.version.empty()) {
        ann.transcript_id = transcript.id + "." + transcript.version;
    }

    // Find exon/intron number (only populate when --numbers is set)
    if (include_numbers_) {
        for (size_t i = 0; i < transcript.exons.size(); ++i) {
            const auto& exon = transcript.exons[i];
            if (pos >= exon.start && pos <= exon.end) {
                ann.exon_number = exon.exon_number;
                break;
            }
            if (i + 1 < transcript.exons.size()) {
                const auto& next_exon = transcript.exons[i + 1];
                if (pos > exon.end && pos < next_exon.start) {
                    ann.intron_number = (transcript.strand == '+') ? static_cast<int>(i + 1) : static_cast<int>(transcript.exons.size() - i - 1);
                    break;
                }
            }
        }
    }

    // Add HGVS_OFFSET if variant was shifted
    if (hgvs_offset != 0) {
        ann.custom_annotations["HGVS_OFFSET"] = std::to_string(hgvs_offset);
    }

    // Add SYMBOL_SOURCE and HGNC_ID from gene record
    if (!transcript.gene_id.empty() && transcript_db_) {
        const Gene* gene = transcript_db_->get_gene(transcript.gene_id);
        if (gene) {
            if (!gene->source.empty()) {
                ann.custom_annotations["SYMBOL_SOURCE"] = gene->source;
            } else if (gene->id.substr(0, 4) == "ENSG") {
                ann.custom_annotations["SYMBOL_SOURCE"] = "HGNC";
            } else {
                ann.custom_annotations["SYMBOL_SOURCE"] = "EntrezGene";
            }
            if (!gene->hgnc_id.empty()) {
                ann.custom_annotations["HGNC_ID"] = gene->hgnc_id;
            }
        }
    }

    // Add FLAGS for incomplete CDS transcripts
    if (transcript.cds_start_NF || transcript.cds_end_NF) {
        std::string flags;
        if (transcript.cds_start_NF) flags += "cds_start_NF";
        if (transcript.cds_end_NF) {
            if (!flags.empty()) flags += "&";
            flags += "cds_end_NF";
        }
        ann.custom_annotations["FLAGS"] = flags;
    }
}

std::vector<ConsequenceType> VEPAnnotator::determine_consequences(
    int pos,
    const std::string& ref,
    const std::string& alt,
    const Transcript& transcript) {

    std::vector<ConsequenceType> consequences;
    int var_end = pos + std::max(0, static_cast<int>(ref.length()) - 1);

    // Check if variant is in exon or intron
    bool in_exon = false;
    bool in_intron = false;
    bool near_splice = false;

    for (size_t i = 0; i < transcript.exons.size(); ++i) {
        const auto& exon = transcript.exons[i];
        bool near_splice_this_exon = false;

        // Check full variant range for exon overlap
        if (var_end >= exon.start && pos <= exon.end) {
            in_exon = true;
        }

        // Check splice sites (2bp at each end) with strand awareness
        bool is_minus = (transcript.strand == '-');

        if (i > 0) {
            // Intron to the left of this exon
            if (!is_minus) {
                // Plus: left side of intron is end of prev exon (we check acceptor at exon.start)
                // Acceptor site at exon.start-2 to exon.start-1
                if (var_end >= exon.start - 2 && pos <= exon.start - 1) {
                    consequences.push_back(ConsequenceType::SPLICE_ACCEPTOR_VARIANT);
                    near_splice = true;
                    near_splice_this_exon = true;
                }
                // Polypyrimidine tract
                if (var_end >= exon.start - 17 && pos <= exon.start - 3) {
                    consequences.push_back(ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT);
                }
            } else {
                // Minus: exon.start boundary is the donor side
                if (var_end >= exon.start - 2 && pos <= exon.start - 1) {
                    consequences.push_back(ConsequenceType::SPLICE_DONOR_VARIANT);
                    near_splice = true;
                    near_splice_this_exon = true;
                }
                // Donor 5th base: exon.start - 5
                if (pos <= exon.start - 5 && var_end >= exon.start - 5) {
                    consequences.push_back(ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT);
                }
                // Donor region: exon.start - 6 to exon.start - 3
                if (var_end >= exon.start - 6 && pos <= exon.start - 3) {
                    consequences.push_back(ConsequenceType::SPLICE_DONOR_REGION_VARIANT);
                }
            }
        }
        if (i < transcript.exons.size() - 1) {
            // Intron to the right of this exon
            if (!is_minus) {
                // Plus: donor at exon.end+1 to exon.end+2
                if (var_end >= exon.end + 1 && pos <= exon.end + 2) {
                    consequences.push_back(ConsequenceType::SPLICE_DONOR_VARIANT);
                    near_splice = true;
                    near_splice_this_exon = true;
                }
                // Donor 5th base
                if (pos <= exon.end + 5 && var_end >= exon.end + 5) {
                    consequences.push_back(ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT);
                }
                // Donor region
                if (var_end >= exon.end + 3 && pos <= exon.end + 6) {
                    consequences.push_back(ConsequenceType::SPLICE_DONOR_REGION_VARIANT);
                }
            } else {
                // Minus: exon.end boundary is the acceptor side
                if (var_end >= exon.end + 1 && pos <= exon.end + 2) {
                    consequences.push_back(ConsequenceType::SPLICE_ACCEPTOR_VARIANT);
                    near_splice = true;
                    near_splice_this_exon = true;
                }
                // Polypyrimidine tract: 3-17bp into intron from acceptor
                if (var_end >= exon.end + 3 && pos <= exon.end + 17) {
                    consequences.push_back(ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT);
                }
            }
        }

        // Check splice region (3-8 bp into intron, 1-3 bp into exon)
        if (!near_splice_this_exon) {
            bool is_splice_region = false;
            // Intronic splice region at exon start (3-8bp before) - skip first exon
            if (i > 0 && var_end >= exon.start - 8 && pos <= exon.start - 3) {
                is_splice_region = true;
            }
            // Exonic splice region at exon start (1-3bp into exon) - skip first exon
            if (i > 0 && var_end >= exon.start && pos <= exon.start + 2) {
                is_splice_region = true;
            }
            // Exonic splice region at exon end (1-3bp into exon end) - skip last exon
            if (i < transcript.exons.size() - 1 && var_end >= exon.end - 2 && pos <= exon.end) {
                is_splice_region = true;
            }
            // Intronic splice region at exon end (3-8bp after) - skip last exon
            if (i < transcript.exons.size() - 1 && var_end >= exon.end + 3 && pos <= exon.end + 8) {
                is_splice_region = true;
            }
            if (is_splice_region) {
                consequences.push_back(ConsequenceType::SPLICE_REGION_VARIANT);
            }
        }

        // Check introns
        if (i + 1 < transcript.exons.size()) {
            const auto& next_exon = transcript.exons[i + 1];
            if (pos > exon.end && var_end < next_exon.start) {
                in_intron = true;
            }
        }
    }

    // Deduplicate splice consequences for variants spanning multiple exon boundaries
    if (consequences.size() > 1) {
        std::sort(consequences.begin(), consequences.end());
        consequences.erase(std::unique(consequences.begin(), consequences.end()), consequences.end());
    }

    // Non-coding transcript
    if (!transcript.is_coding()) {
        if (in_exon) {
            // Detect mature miRNA variant
            if (transcript.biotype == "miRNA") {
                consequences.push_back(ConsequenceType::MATURE_MIRNA_VARIANT);
            } else {
                consequences.push_back(ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT);
            }
        } else if (in_intron) {
            if (consequences.empty()) {
                consequences.push_back(ConsequenceType::INTRON_VARIANT);
            }
        } else {
            consequences.push_back(ConsequenceType::NON_CODING_TRANSCRIPT_VARIANT);
        }
        return consequences;
    }

    // For coding transcripts - intronic
    if (in_intron && !near_splice && consequences.empty()) {
        consequences.push_back(ConsequenceType::INTRON_VARIANT);
        return consequences;
    }
    // If only splice-related consequences were added for intronic variant, add intron_variant too
    if (in_intron && !consequences.empty()) {
        bool has_non_splice = false;
        for (const auto& c : consequences) {
            if (c != ConsequenceType::SPLICE_ACCEPTOR_VARIANT &&
                c != ConsequenceType::SPLICE_DONOR_VARIANT &&
                c != ConsequenceType::SPLICE_REGION_VARIANT &&
                c != ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT &&
                c != ConsequenceType::SPLICE_DONOR_REGION_VARIANT &&
                c != ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT) {
                has_non_splice = true;
                break;
            }
        }
        if (!has_non_splice) {
            return consequences; // Return splice consequences only
        }
    }

    // Check UTR regions
    if (in_exon) {
        bool in_5utr = false;
        bool in_3utr = false;

        if (transcript.strand == '+') {
            if (pos < transcript.cds_start) in_5utr = true;
            if (var_end > transcript.cds_end) in_3utr = true;
        } else {
            if (var_end > transcript.cds_end) in_5utr = true;
            if (pos < transcript.cds_start) in_3utr = true;
        }

        if (in_5utr) {
            consequences.push_back(ConsequenceType::FIVE_PRIME_UTR_VARIANT);
        }
        if (in_3utr) {
            consequences.push_back(ConsequenceType::THREE_PRIME_UTR_VARIANT);
        }
    }

    // Coding region variant - determine specific consequence
    if (in_exon && pos <= transcript.cds_end && var_end >= transcript.cds_start) {
        int cds_pos = calculate_cds_position(pos, transcript);
        // For minus strand with multi-base ref, adjust from last to first affected CDS position
        if (transcript.strand == '-' && ref.length() > 1) {
            cds_pos = cds_pos - static_cast<int>(ref.length()) + 1;
            if (cds_pos < 1) cds_pos = 1;
        }

        if (cds_pos > 0) {
            // Check for incomplete terminal codon
            int cds_length = transcript.get_cds_length();
            if (cds_length % 3 != 0) {
                int incomplete_start = cds_length - (cds_length % 3) + 1;
                if (cds_pos >= incomplete_start) {
                    consequences.push_back(ConsequenceType::INCOMPLETE_TERMINAL_CODON_VARIANT);
                    return consequences;
                }
            }

            int ref_len = static_cast<int>(ref.length());
            int alt_len = static_cast<int>(alt.length());
            int length_diff = alt_len - ref_len;

            // Start codon check for indels - start_lost if variant deletes/disrupts
            // bases within the start codon (CDS pos 1-3)
            // Only for deletions/frameshifts that actually overlap start codon bases
            if (length_diff != 0 && cds_pos >= 1 && cds_pos <= 3) {
                // For deletions: variant removes bases from start codon
                // For insertions at cds_pos 1 or 2: disrupts start codon
                // But insertions after pos 3 (i.e., cds_pos > 3) don't affect start codon
                bool affects_start = false;
                if (length_diff < 0) {
                    // Deletion overlapping start codon bases
                    affects_start = true;
                } else if (cds_pos <= 2) {
                    // Insertion within start codon (between base 1-2 or 2-3)
                    affects_start = true;
                }
                // Insertions at cds_pos == 3 (after last base of start codon) don't affect it
                if (affects_start) {
                    consequences.push_back(ConsequenceType::START_LOST);
                    return consequences;
                }
            }

            // Stop codon check for inframe deletions overlapping the stop codon
            if (length_diff < 0 && length_diff % 3 == 0) {
                int stop_codon_start = cds_length - 2;
                int del_end = cds_pos + ref_len - 1;
                if (stop_codon_start > 0 && del_end >= stop_codon_start) {
                    consequences.push_back(ConsequenceType::STOP_LOST);
                    return consequences;
                }
            }

            // Frameshift
            if (length_diff != 0 && length_diff % 3 != 0) {
                consequences.push_back(ConsequenceType::FRAMESHIFT_VARIANT);
                return consequences;
            }

            // Inframe indel
            if (length_diff != 0 && length_diff % 3 == 0) {
                if (ref_len > 1 && alt_len > 1 && ref_len != alt_len) {
                    // Complex inframe change (delins)
                    consequences.push_back(ConsequenceType::PROTEIN_ALTERING_VARIANT);
                } else if (length_diff > 0) {
                    consequences.push_back(ConsequenceType::INFRAME_INSERTION);
                } else {
                    consequences.push_back(ConsequenceType::INFRAME_DELETION);
                }
                return consequences;
            }

            // SNV or MNV (same length substitution)
            if (ref.length() == alt.length()) {
                // Build CDS sequence for codon extraction
                std::string cds_seq;
                if (transcript.strand == '+') {
                    for (const auto& cds : transcript.cds_regions) {
                        cds_seq += reference_->get_sequence(transcript.chromosome, cds.start, cds.end);
                    }
                } else {
                    for (auto it = transcript.cds_regions.rbegin(); it != transcript.cds_regions.rend(); ++it) {
                        std::string seq = reference_->get_sequence(transcript.chromosome, it->start, it->end);
                        std::reverse(seq.begin(), seq.end());
                        for (char& c : seq) c = complement_base(c);
                        cds_seq += seq;
                    }
                }

                if (cds_seq.empty()) {
                    consequences.push_back(ConsequenceType::CODING_SEQUENCE_VARIANT);
                    return consequences;
                }

                // Determine affected codon range
                int cds_end_pos = cds_pos + static_cast<int>(ref.length()) - 1;
                int first_codon = (cds_pos - 1) / 3;
                int last_codon = (cds_end_pos - 1) / 3;

                int codon_region_start = first_codon * 3;
                int codon_region_end = (last_codon + 1) * 3;
                if (codon_region_end > static_cast<int>(cds_seq.length())) {
                    codon_region_end = static_cast<int>(cds_seq.length());
                }

                if (codon_region_start >= static_cast<int>(cds_seq.length())) {
                    consequences.push_back(ConsequenceType::CODING_SEQUENCE_VARIANT);
                    return consequences;
                }

                // Create mutated CDS segment
                std::string ref_segment = cds_seq.substr(codon_region_start, codon_region_end - codon_region_start);
                std::string alt_segment = ref_segment;

                // Apply the substitution
                int offset_in_segment = (cds_pos - 1) - codon_region_start;
                for (size_t k = 0; k < ref.length() && (offset_in_segment + static_cast<int>(k)) < static_cast<int>(alt_segment.length()); ++k) {
                    // For minus strand, reverse the index: genomic order is reversed in CDS
                    size_t alt_idx = (transcript.strand == '-') ? (alt.length() - 1 - k) : k;
                    char alt_base = alt[alt_idx];
                    if (transcript.strand == '-') {
                        alt_base = complement_base(alt_base);
                    }
                    alt_segment[offset_in_segment + k] = alt_base;
                }

                // Translate both segments
                std::string ref_protein, alt_protein;
                for (int c = 0; c + 2 < static_cast<int>(ref_segment.length()); c += 3) {
                    ref_protein += CodonTable::translate(ref_segment.substr(c, 3), transcript.chromosome);
                }
                for (int c = 0; c + 2 < static_cast<int>(alt_segment.length()); c += 3) {
                    alt_protein += CodonTable::translate(alt_segment.substr(c, 3), transcript.chromosome);
                }

                // Check start codon (first codon of CDS)
                if (first_codon == 0 && !ref_segment.empty() && !alt_segment.empty()) {
                    std::string ref_start = ref_segment.substr(0, 3);
                    std::string alt_start = alt_segment.substr(0, 3);
                    if (CodonTable::is_start_codon(ref_start)) {
                        if (!CodonTable::is_start_codon(alt_start)) {
                            consequences.push_back(ConsequenceType::START_LOST);
                            return consequences;
                        } else if (ref_start != alt_start) {
                            consequences.push_back(ConsequenceType::START_RETAINED_VARIANT);
                            return consequences;
                        }
                    }
                }

                // Check for stop gained/lost/retained
                bool ref_has_stop = (ref_protein.find('*') != std::string::npos);
                bool alt_has_stop = (alt_protein.find('*') != std::string::npos);

                if (alt_has_stop && !ref_has_stop) {
                    consequences.push_back(ConsequenceType::STOP_GAINED);
                    return consequences;
                }
                if (ref_has_stop && !alt_has_stop) {
                    consequences.push_back(ConsequenceType::STOP_LOST);
                    return consequences;
                }
                if (ref_has_stop && alt_has_stop && ref_segment != alt_segment) {
                    consequences.push_back(ConsequenceType::STOP_RETAINED_VARIANT);
                    return consequences;
                }

                // Synonymous or missense
                if (ref_protein == alt_protein) {
                    consequences.push_back(ConsequenceType::SYNONYMOUS_VARIANT);
                } else {
                    consequences.push_back(ConsequenceType::MISSENSE_VARIANT);
                }
            }
        }
    }

    if (consequences.empty()) {
        if (transcript.is_coding()) {
            consequences.push_back(ConsequenceType::CODING_TRANSCRIPT_VARIANT);
        } else {
            consequences.push_back(ConsequenceType::SEQUENCE_VARIANT);
        }
    }

    return consequences;
}

int VEPAnnotator::calculate_cds_position(int genomic_pos, const Transcript& transcript) const {
    if (!transcript.is_coding()) return 0;
    if (genomic_pos < transcript.cds_start || genomic_pos > transcript.cds_end) return 0;

    int cds_pos = 0;
    bool in_cds = false;

    if (transcript.strand == '+') {
        for (const auto& cds : transcript.cds_regions) {
            if (genomic_pos > cds.end) {
                cds_pos += cds.end - cds.start + 1;
            } else if (genomic_pos >= cds.start) {
                cds_pos += genomic_pos - cds.start + 1;
                in_cds = true;
                break;
            } else {
                // Position is before this CDS region (in an intron)
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
                in_cds = true;
                break;
            } else {
                // Position is past this CDS region (in an intron)
                break;
            }
        }
    }

    // Adjust for CDS phase in cds_start_NF transcripts
    if (in_cds && transcript.cds_start_NF && !transcript.cds_regions.empty()) {
        int first_phase = (transcript.strand == '+')
            ? transcript.cds_regions.front().phase
            : transcript.cds_regions.back().phase;
        if (first_phase > 0 && first_phase <= 2) {
            cds_pos -= first_phase;
            if (cds_pos <= 0) return 0; // Position is in the partial codon
        }
    }

    return in_cds ? cds_pos : 0;
}

std::string VEPAnnotator::build_cds_sequence(
    const std::string& chrom, const Transcript& transcript) const {

    std::string cds_seq;
    if (transcript.strand == '+') {
        for (const auto& cds : transcript.cds_regions) {
            cds_seq += reference_->get_sequence(chrom, cds.start, cds.end);
        }
    } else {
        for (auto it = transcript.cds_regions.rbegin(); it != transcript.cds_regions.rend(); ++it) {
            std::string seg = reference_->get_sequence(chrom, it->start, it->end);
            std::reverse(seg.begin(), seg.end());
            for (char& c : seg) c = complement_base(c);
            cds_seq += seg;
        }
    }

    // For cds_start_NF transcripts, trim leading bases per first CDS phase
    if (transcript.cds_start_NF && !transcript.cds_regions.empty()) {
        int first_phase = (transcript.strand == '+')
            ? transcript.cds_regions.front().phase
            : transcript.cds_regions.back().phase;
        if (first_phase > 0 && first_phase <= 2 &&
            static_cast<int>(cds_seq.size()) > first_phase) {
            cds_seq = cds_seq.substr(first_phase);
        }
    }

    return cds_seq;
}

std::pair<std::string, std::string> VEPAnnotator::get_affected_codons(
    int cds_pos,
    const std::string& ref,
    const std::string& alt,
    const Transcript& transcript) const {

    return get_affected_codons(cds_pos, ref, alt, transcript,
                               build_cds_sequence(transcript.chromosome, transcript));
}

std::pair<std::string, std::string> VEPAnnotator::get_affected_codons(
    int cds_pos,
    const std::string& ref,
    const std::string& alt,
    const Transcript& transcript,
    const std::string& cds_seq) const {

    if (cds_pos <= 0) return {"", ""};

    // Calculate codon position (0-based codon number)
    int codon_num = (cds_pos - 1) / 3;
    int pos_in_codon = (cds_pos - 1) % 3;

    if (cds_seq.empty()) return {"", ""};

    // Get the reference codon
    int codon_start = codon_num * 3;
    if (codon_start + 3 > static_cast<int>(cds_seq.length())) {
        return {"", ""};
    }

    std::string ref_codon = cds_seq.substr(codon_start, 3);

    // Create the alternate codon
    std::string alt_codon = ref_codon;
    if (ref.length() == alt.length()) {
        // Handle SNVs and MNVs (same-length substitutions)
        for (size_t i = 0; i < alt.length(); ++i) {
            int codon_offset = pos_in_codon + static_cast<int>(i);
            if (codon_offset >= 0 && codon_offset < 3) {
                // For minus strand, reverse the index: genomic order is reversed in CDS
                size_t alt_idx = (transcript.strand == '-') ? (alt.length() - 1 - i) : i;
                char alt_base = alt[alt_idx];
                if (transcript.strand == '-') {
                    alt_base = complement_base(alt_base);
                }
                alt_codon[codon_offset] = alt_base;
            }
        }
    } else {
        // Handle indels (insertions and deletions)
        // Build mutated CDS by applying the indel, then extract the alt codon
        std::string mut_cds;
        int ref_len = static_cast<int>(ref.length());

        if (transcript.strand == '+') {
            int cds_offset = cds_pos - 1; // 0-based position in CDS
            if (cds_offset < 0 || cds_offset + ref_len > static_cast<int>(cds_seq.length())) {
                return {ref_codon, ref_codon};
            }
            mut_cds = cds_seq.substr(0, cds_offset) + alt + cds_seq.substr(cds_offset + ref_len);
        } else {
            // For minus strand: cds_pos is the 1-based CDS position of the first
            // affected base. The 0-based offset is cds_pos - 1, same as plus strand.
            // Replace ref_len bases starting at cds_offset with rc(alt).
            std::string rc_alt = reverse_complement_sequence(alt);
            int cds_offset = cds_pos - 1; // 0-based start of affected region
            if (cds_offset < 0 || cds_offset + ref_len > static_cast<int>(cds_seq.length())) {
                return {ref_codon, ref_codon};
            }
            mut_cds = cds_seq.substr(0, cds_offset) + rc_alt + cds_seq.substr(cds_offset + ref_len);
        }

        // Extract the alt codon at the same codon position from the mutated CDS
        if (codon_start + 3 <= static_cast<int>(mut_cds.length())) {
            alt_codon = mut_cds.substr(codon_start, 3);
        }
        // If mutated CDS is too short (e.g., large deletion), alt_codon stays as ref_codon
    }

    return {ref_codon, alt_codon};
}

std::string VEPAnnotator::generate_hgvsc(
    int cds_pos,
    const std::string& ref,
    const std::string& alt,
    const Transcript& transcript) const {

    return generate_hgvsc(cds_pos, ref, alt, transcript,
                          build_cds_sequence(transcript.chromosome, transcript));
}

std::string VEPAnnotator::generate_hgvsc(
    int cds_pos,
    const std::string& ref,
    const std::string& alt,
    const Transcript& transcript,
    const std::string& cds_seq) const {

    std::ostringstream oss;
    std::string tid = transcript.id;
    if (transcript_version_ && !transcript.version.empty()) {
        tid += "." + transcript.version;
    }
    oss << tid << ":c.";

    // For minus-strand transcripts, complement the alleles to transcript orientation
    std::string hgvs_ref = ref;
    std::string hgvs_alt = alt;
    if (transcript.strand == '-') {
        hgvs_ref = reverse_complement_sequence(ref);
        hgvs_alt = reverse_complement_sequence(alt);
    }

    if (hgvs_ref.length() == 1 && hgvs_alt.length() == 1) {
        // SNV
        oss << cds_pos << hgvs_ref << ">" << hgvs_alt;
    } else if (hgvs_ref.length() > hgvs_alt.length()) {
        // Deletion
        int del_len = static_cast<int>(hgvs_ref.length() - hgvs_alt.length());
        if (del_len == 1) {
            oss << (cds_pos + static_cast<int>(hgvs_alt.length())) << "del";
        } else {
            oss << (cds_pos + static_cast<int>(hgvs_alt.length())) << "_"
                << (cds_pos + static_cast<int>(hgvs_ref.length()) - 1) << "del";
        }
    } else if (hgvs_alt.length() > hgvs_ref.length()) {
        // Insertion - check for duplication
        std::string inserted = hgvs_alt.substr(hgvs_ref.length());
        int ins_len = static_cast<int>(inserted.length());

        bool is_dup = false;

        // Check if inserted bases match the preceding reference context (0-based)
        int cds_idx = cds_pos - 1; // Convert to 0-based
        if (cds_idx >= ins_len && cds_idx <= static_cast<int>(cds_seq.length())) {
            std::string preceding = cds_seq.substr(cds_idx - ins_len, ins_len);
            std::string ins_upper = inserted;
            std::string prec_upper = preceding;
            for (auto& c : ins_upper) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            for (auto& c : prec_upper) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            if (ins_upper == prec_upper) {
                is_dup = true;
            }
        }

        if (is_dup) {
            int dup_start = cds_pos - ins_len + 1;
            int dup_end = cds_pos;
            if (ins_len == 1) {
                oss << dup_end << "dup";
            } else {
                oss << dup_start << "_" << dup_end << "dup";
            }
        } else {
            oss << cds_pos << "_" << (cds_pos + 1) << "ins" << inserted;
        }
    } else if (hgvs_ref.length() > 1 && hgvs_alt.length() > 1 &&
               hgvs_ref.length() != hgvs_alt.length()) {
        std::string delins_ref = hgvs_ref;
        std::string delins_alt = hgvs_alt;
        int delins_pos = cds_pos;
        if (delins_ref.length() > 1 && delins_alt.length() > 1 &&
            delins_ref[0] == delins_alt[0]) {
            delins_ref = delins_ref.substr(1);
            delins_alt = delins_alt.substr(1);
            delins_pos += 1;
        }
        int del_end = delins_pos + static_cast<int>(delins_ref.length()) - 1;
        if (delins_pos == del_end) {
            oss << delins_pos << "delins" << delins_alt;
        } else {
            oss << delins_pos << "_" << del_end << "delins" << delins_alt;
        }
    } else {
        int end_pos = cds_pos + static_cast<int>(hgvs_ref.length()) - 1;
        if (cds_pos == end_pos) {
            oss << cds_pos << hgvs_ref << ">" << hgvs_alt;
        } else {
            oss << cds_pos << "_" << end_pos << "delins" << hgvs_alt;
        }
    }

    return oss.str();
}

/**
 * Calculate intronic HGVSc position string for a variant in an intron.
 * Returns a string like "123+5" or "124-3" representing the intronic offset
 * from the nearest exon boundary in CDS coordinates, or empty if not intronic.
 */
static std::string calculate_intronic_hgvsc_position(
    int genomic_pos,
    const Transcript& transcript,
    const std::function<int(int, const Transcript&)>& calc_cds_pos) {

    // Need sorted exons in genomic order
    std::vector<Exon> sorted_exons = transcript.exons;
    std::sort(sorted_exons.begin(), sorted_exons.end(),
              [](const Exon& a, const Exon& b) { return a.start < b.start; });

    // Find which intron the variant falls in (between sorted_exons[i] and sorted_exons[i+1])
    for (size_t i = 0; i + 1 < sorted_exons.size(); ++i) {
        int intron_start = sorted_exons[i].end + 1;
        int intron_end = sorted_exons[i + 1].start - 1;

        if (genomic_pos >= intron_start && genomic_pos <= intron_end) {
            // Found the intron
            int dist_to_upstream_exon_end = genomic_pos - sorted_exons[i].end;
            int dist_to_downstream_exon_start = sorted_exons[i + 1].start - genomic_pos;

            if (transcript.strand == '+') {
                if (dist_to_upstream_exon_end <= dist_to_downstream_exon_start) {
                    // Closer to upstream exon end: c.CDS_POS+offset
                    int cds_of_exon_end = calc_cds_pos(sorted_exons[i].end, transcript);
                    if (cds_of_exon_end > 0) {
                        return std::to_string(cds_of_exon_end) + "+" + std::to_string(dist_to_upstream_exon_end);
                    }
                } else {
                    // Closer to downstream exon start: c.CDS_POS-offset
                    int cds_of_exon_start = calc_cds_pos(sorted_exons[i + 1].start, transcript);
                    if (cds_of_exon_start > 0) {
                        return std::to_string(cds_of_exon_start) + "-" + std::to_string(dist_to_downstream_exon_start);
                    }
                }
            } else {
                // Minus strand: genomic upstream exon end is transcript downstream, and vice versa
                if (dist_to_downstream_exon_start <= dist_to_upstream_exon_end) {
                    // Closer to the genomic downstream exon start (which is transcript upstream)
                    int cds_of_exon_start = calc_cds_pos(sorted_exons[i + 1].start, transcript);
                    if (cds_of_exon_start > 0) {
                        return std::to_string(cds_of_exon_start) + "+" + std::to_string(dist_to_downstream_exon_start);
                    }
                } else {
                    // Closer to the genomic upstream exon end (which is transcript downstream)
                    int cds_of_exon_end = calc_cds_pos(sorted_exons[i].end, transcript);
                    if (cds_of_exon_end > 0) {
                        return std::to_string(cds_of_exon_end) + "-" + std::to_string(dist_to_upstream_exon_end);
                    }
                }
            }
            break;
        }
    }
    return "";
}

/**
 * Generate HGVSc notation for intronic variants using pre-calculated intronic position.
 */
static std::string generate_hgvsc_intronic(
    const std::string& intronic_pos,
    const std::string& ref,
    const std::string& alt,
    const Transcript& transcript,
    bool include_version) {

    std::ostringstream oss;
    std::string tid = transcript.id;
    if (include_version && !transcript.version.empty()) {
        tid += "." + transcript.version;
    }
    oss << tid << ":c.";

    // For minus-strand transcripts, complement the alleles to transcript orientation
    std::string hgvs_ref = ref;
    std::string hgvs_alt = alt;
    if (transcript.strand == '-') {
        hgvs_ref = reverse_complement_sequence(ref);
        hgvs_alt = reverse_complement_sequence(alt);
    }

    if (hgvs_ref.length() == 1 && hgvs_alt.length() == 1) {
        // SNV
        oss << intronic_pos << hgvs_ref << ">" << hgvs_alt;
    } else if (hgvs_ref.length() > hgvs_alt.length()) {
        // Deletion
        oss << intronic_pos << "del";
    } else if (hgvs_alt.length() > hgvs_ref.length()) {
        // Insertion
        std::string inserted = hgvs_alt.substr(hgvs_ref.length());
        oss << intronic_pos << "ins" << inserted;
    } else {
        // Substitution (MNV)
        oss << intronic_pos << "delins" << hgvs_alt;
    }

    return oss.str();
}

/**
 * Generate HGVSc notation for UTR variants.
 * 5'UTR uses c.-N notation, 3'UTR uses c.*N notation.
 */
static std::string generate_hgvsc_utr(
    const std::string& utr_pos_str,
    const std::string& ref,
    const std::string& alt,
    const Transcript& transcript,
    bool include_version) {

    std::ostringstream oss;
    std::string tid = transcript.id;
    if (include_version && !transcript.version.empty()) {
        tid += "." + transcript.version;
    }
    oss << tid << ":c.";

    std::string hgvs_ref = ref;
    std::string hgvs_alt = alt;
    if (transcript.strand == '-') {
        hgvs_ref = reverse_complement_sequence(ref);
        hgvs_alt = reverse_complement_sequence(alt);
    }

    if (hgvs_ref.length() == 1 && hgvs_alt.length() == 1) {
        oss << utr_pos_str << hgvs_ref << ">" << hgvs_alt;
    } else if (hgvs_ref.length() > hgvs_alt.length()) {
        oss << utr_pos_str << "del";
    } else if (hgvs_alt.length() > hgvs_ref.length()) {
        std::string inserted = hgvs_alt.substr(hgvs_ref.length());
        oss << utr_pos_str << "ins" << inserted;
    } else {
        oss << utr_pos_str << "delins" << hgvs_alt;
    }

    return oss.str();
}

std::string VEPAnnotator::generate_hgvsp(
    const std::string& ref_aa,
    const std::string& alt_aa,
    int protein_pos,
    const Transcript& transcript,
    const std::vector<ConsequenceType>& consequences,
    int protein_end,
    int fs_ter_distance,
    const std::string& end_ref_aa) const {

    if (ref_aa.empty() || alt_aa.empty()) return "";

    // Check consequence types for special notation
    bool is_frameshift = false;
    bool is_inframe_del = false;
    bool is_inframe_ins = false;
    bool is_stop_lost = false;
    bool is_start_lost = false;
    bool is_protein_altering = false;
    for (const auto& c : consequences) {
        if (c == ConsequenceType::FRAMESHIFT_VARIANT) is_frameshift = true;
        if (c == ConsequenceType::INFRAME_DELETION) is_inframe_del = true;
        if (c == ConsequenceType::INFRAME_INSERTION) is_inframe_ins = true;
        if (c == ConsequenceType::STOP_LOST) is_stop_lost = true;
        if (c == ConsequenceType::START_LOST) is_start_lost = true;
        if (c == ConsequenceType::PROTEIN_ALTERING_VARIANT) is_protein_altering = true;
    }

    std::ostringstream oss;
    // Use protein ID when available, fall back to transcript ID (with version if enabled)
    std::string protein_ref;
    if (!transcript.protein_id.empty()) {
        protein_ref = transcript.protein_id;
    } else {
        protein_ref = transcript.id;
        if (transcript_version_ && !transcript.version.empty()) {
            protein_ref += "." + transcript.version;
        }
    }
    oss << protein_ref << ":p.";

    if (is_frameshift) {
        // Frameshift: p.Ala123GlyfsTer5 (Perl VEP format)
        oss << ref_aa << protein_pos << alt_aa << "fs";
        if (fs_ter_distance > 0) {
            oss << "Ter" << fs_ter_distance;
        } else {
            oss << "Ter?";
        }
    } else if (is_stop_lost && ref_aa == "Ter") {
        // Stop lost / extension: p.Ter123AlaextTer?
        oss << ref_aa << protein_pos << alt_aa << "extTer?";
    } else if (is_start_lost) {
        // Start lost: p.Met1?
        oss << ref_aa << protein_pos << "?";
    } else if (is_inframe_del) {
        // Inframe deletion: p.Ala123del or p.Ala123_Gly125del (range)
        oss << ref_aa << protein_pos;
        if (protein_end > 0 && protein_end > protein_pos) {
            // Multi-residue deletion: show range with end AA name
            std::string end_aa_name = end_ref_aa.empty() ? ref_aa : end_ref_aa;
            oss << "_" << end_aa_name << protein_end;
        }
        oss << "del";
    } else if (is_inframe_ins) {
        // Inframe insertion: p.Ala123_Gly124insLeu
        std::string next_aa = end_ref_aa.empty() ? ref_aa : end_ref_aa;
        oss << ref_aa << protein_pos << "_" << next_aa << (protein_pos + 1) << "ins" << alt_aa;
    } else if (is_protein_altering) {
        // Complex inframe (delins): p.Ala123_Gly125delinsVal
        oss << ref_aa << protein_pos;
        if (protein_end > 0 && protein_end > protein_pos) {
            std::string end_aa_name = end_ref_aa.empty() ? ref_aa : end_ref_aa;
            oss << "_" << end_aa_name << protein_end;
        }
        oss << "delins" << alt_aa;
    } else if (ref_aa == alt_aa) {
        // Synonymous: p.Ala123=
        oss << ref_aa << protein_pos << "=";
    } else if (alt_aa == "Ter") {
        // Stop gained: p.Gln123Ter (or p.Gln123*)
        oss << ref_aa << protein_pos << "Ter";
    } else {
        // Missense: p.Ala123Gly
        oss << ref_aa << protein_pos << alt_aa;
    }

    return oss.str();
}

const Transcript* VEPAnnotator::get_transcript(const std::string& transcript_id) const {
    return transcript_db_->get_transcript(transcript_id);
}

int VEPAnnotator::map_cds_to_genomic(int cds_pos, const Transcript& transcript) const {
    if (!transcript.is_coding() || cds_pos <= 0) return 0;

    int remaining = cds_pos;
    if (transcript.strand == '+') {
        for (const auto& cds : transcript.cds_regions) {
            int cds_len = cds.end - cds.start + 1;
            if (remaining <= cds_len) {
                return cds.start + remaining - 1;
            }
            remaining -= cds_len;
        }
    } else {
        for (auto it = transcript.cds_regions.rbegin(); it != transcript.cds_regions.rend(); ++it) {
            int cds_len = it->end - it->start + 1;
            if (remaining <= cds_len) {
                return it->end - remaining + 1;
            }
            remaining -= cds_len;
        }
    }
    return 0;
}

int VEPAnnotator::calculate_cdna_position(int genomic_pos, const Transcript& transcript) const {
    int cdna_pos = 0;
    if (transcript.strand == '+') {
        for (const auto& exon : transcript.exons) {
            if (genomic_pos > exon.end) {
                cdna_pos += exon.end - exon.start + 1;
            } else if (genomic_pos >= exon.start) {
                cdna_pos += genomic_pos - exon.start + 1;
                return cdna_pos;
            }
        }
    } else {
        for (auto it = transcript.exons.rbegin(); it != transcript.exons.rend(); ++it) {
            if (genomic_pos < it->start) {
                cdna_pos += it->end - it->start + 1;
            } else if (genomic_pos <= it->end) {
                cdna_pos += it->end - genomic_pos + 1;
                return cdna_pos;
            }
        }
    }
    return 0; // Not in an exon
}

void VEPAnnotator::append_regulatory_consequences(VariantAnnotation& ann) {
    auto tfbs_it = ann.custom_annotations.find("regulatory:in_tfbs");
    if (tfbs_it != ann.custom_annotations.end() && tfbs_it->second == "true") {
        ann.consequences.push_back(ConsequenceType::TF_BINDING_SITE_VARIANT);
    }

    auto reg_type_it = ann.custom_annotations.find("regulatory:feature_type");
    if (reg_type_it != ann.custom_annotations.end()) {
        // Case-insensitive comparison for GFF3 feature type variations
        std::string ft = reg_type_it->second;
        for (auto& c : ft) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (ft.find("promoter") != std::string::npos ||
            ft.find("enhancer") != std::string::npos ||
            ft.find("open_chromatin") != std::string::npos ||
            ft.find("ctcf") != std::string::npos) {
            ann.consequences.push_back(ConsequenceType::REGULATORY_REGION_VARIANT);
        }
    }

    // Recalculate impact if regulatory consequences added
    if (!ann.consequences.empty()) {
        ann.impact = Impact::MODIFIER;
        for (const auto& c : ann.consequences) {
            Impact i = get_impact(c);
            if (static_cast<int>(i) < static_cast<int>(ann.impact)) {
                ann.impact = i;
            }
        }
    }
}

std::tuple<int, std::string, std::string> VEPAnnotator::right_normalize(
    const std::string& chrom, int pos,
    const std::string& ref, const std::string& alt,
    const Transcript& transcript) const {

    if (ref.length() == alt.length()) {
        return {pos, ref, alt};
    }

    // Determine if insertion or deletion
    std::string shorter = (ref.length() < alt.length()) ? ref : alt;
    std::string longer = (ref.length() >= alt.length()) ? ref : alt;

    // Simple VCF-style indels: shared leading base
    if (shorter.length() < 1 || shorter[0] != longer[0]) {
        return {pos, ref, alt};
    }

    std::string indel_seq = longer.substr(shorter.length());
    if (indel_seq.empty()) {
        return {pos, ref, alt};
    }

    bool is_insertion = (alt.length() > ref.length());

    // 3' shift: slide rightward while bases match
    int shifted_pos = pos;
    std::string shifted_indel = indel_seq;
    int chrom_len = reference_->get_chromosome_length(chrom);
    int anchor_offset = static_cast<int>(shorter.length());

    for (int i = 0; i < 1000; ++i) {
        // For deletions: check base after deleted region
        // For insertions: check base right after anchor (don't add indel length)
        int check_pos = shifted_pos + anchor_offset + (is_insertion ? 0 : static_cast<int>(shifted_indel.length()));
        if (check_pos > chrom_len || check_pos <= 0) break;

        char next_base = reference_->get_base(chrom, check_pos);
        if (static_cast<unsigned char>(std::toupper(next_base)) ==
            static_cast<unsigned char>(std::toupper(shifted_indel[0]))) {
            std::rotate(shifted_indel.begin(), shifted_indel.begin() + 1, shifted_indel.end());
            shifted_pos++;
        } else {
            break;
        }
    }

    // Reconstruct
    char anchor = reference_->get_base(chrom, shifted_pos);
    std::string new_shorter(1, anchor);
    std::string new_longer = std::string(1, anchor) + shifted_indel;

    if (is_insertion) {
        return {shifted_pos, new_shorter, new_longer};
    } else {
        return {shifted_pos, new_longer, new_shorter};
    }
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
