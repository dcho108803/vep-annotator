/**
 * HGVS Notation Parser
 *
 * Parses HGVS notation into variant coordinates:
 * - c. notation (coding DNA): ENST00000366667:c.803C>T
 * - g. notation (genomic): NC_000007.14:g.140753336A>T
 * - p. notation (protein): BRAF:p.Val600Glu
 * - n. notation (non-coding): NR_024540.1:n.1234A>G
 */

#ifndef HGVS_PARSER_HPP
#define HGVS_PARSER_HPP

#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <algorithm>

namespace vep {

/**
 * HGVS coordinate types
 */
enum class HGVSType {
    GENOMIC,    // g.
    CODING,     // c.
    NONCODING,  // n.
    PROTEIN,    // p.
    MITOCHONDRIAL, // m.
    RNA,        // r.
    UNKNOWN
};

/**
 * Variant type from HGVS
 */
enum class HGVSVariantType {
    SUBSTITUTION,   // >
    DELETION,       // del
    INSERTION,      // ins
    DUPLICATION,    // dup
    INVERSION,      // inv
    DELINS,         // delins
    REPEAT,         // []
    UNKNOWN
};

/**
 * Parsed HGVS result
 */
struct HGVSParseResult {
    bool valid = false;
    std::string error_message;

    // Reference sequence
    std::string reference_id;    // ENST00000366667, NC_000007.14, NM_000546.5, etc.
    std::string gene_symbol;     // Optional gene symbol

    // HGVS type
    HGVSType hgvs_type = HGVSType::UNKNOWN;
    HGVSVariantType variant_type = HGVSVariantType::UNKNOWN;

    // Coordinates (may need transcript mapping for c./p. notation)
    int start_pos = 0;
    int end_pos = 0;
    int intron_offset = 0;       // For intronic positions like c.123+5

    // Alleles
    std::string ref_allele;
    std::string alt_allele;

    // For protein notation
    std::string ref_aa;
    std::string alt_aa;
    int protein_pos = 0;

    // Converted genomic coordinates (after transcript mapping)
    std::string chromosome;
    int genomic_pos = 0;
    std::string genomic_ref;
    std::string genomic_alt;
};

/**
 * Three-letter to one-letter amino acid conversion
 */
inline std::string aa_three_to_one(const std::string& three_letter) {
    static const std::unordered_map<std::string, std::string> aa_map = [] {
        std::unordered_map<std::string, std::string> m;
        m["Ala"] = "A"; m["ALA"] = "A";
        m["Arg"] = "R"; m["ARG"] = "R";
        m["Asn"] = "N"; m["ASN"] = "N";
        m["Asp"] = "D"; m["ASP"] = "D";
        m["Cys"] = "C"; m["CYS"] = "C";
        m["Gln"] = "Q"; m["GLN"] = "Q";
        m["Glu"] = "E"; m["GLU"] = "E";
        m["Gly"] = "G"; m["GLY"] = "G";
        m["His"] = "H"; m["HIS"] = "H";
        m["Ile"] = "I"; m["ILE"] = "I";
        m["Leu"] = "L"; m["LEU"] = "L";
        m["Lys"] = "K"; m["LYS"] = "K";
        m["Met"] = "M"; m["MET"] = "M";
        m["Phe"] = "F"; m["PHE"] = "F";
        m["Pro"] = "P"; m["PRO"] = "P";
        m["Ser"] = "S"; m["SER"] = "S";
        m["Thr"] = "T"; m["THR"] = "T";
        m["Trp"] = "W"; m["TRP"] = "W";
        m["Tyr"] = "Y"; m["TYR"] = "Y";
        m["Val"] = "V"; m["VAL"] = "V";
        m["Ter"] = "*"; m["TER"] = "*"; m["Stop"] = "*";
        m["Sec"] = "U"; m["SEC"] = "U";  // Selenocysteine
        m["Pyl"] = "O"; m["PYL"] = "O";  // Pyrrolysine
        return m;
    }();

    auto it = aa_map.find(three_letter);
    if (it != aa_map.end()) {
        return it->second;
    }

    // Already one letter or unknown
    if (three_letter.size() == 1) {
        return three_letter;
    }

    return "X";  // Unknown
}

/**
 * One-letter to three-letter amino acid conversion
 */
inline std::string aa_one_to_three(char one_letter) {
    static const std::unordered_map<char, std::string> aa_map = [] {
        std::unordered_map<char, std::string> m;
        m['A'] = "Ala"; m['R'] = "Arg"; m['N'] = "Asn";
        m['D'] = "Asp"; m['C'] = "Cys"; m['Q'] = "Gln";
        m['E'] = "Glu"; m['G'] = "Gly"; m['H'] = "His";
        m['I'] = "Ile"; m['L'] = "Leu"; m['K'] = "Lys";
        m['M'] = "Met"; m['F'] = "Phe"; m['P'] = "Pro";
        m['S'] = "Ser"; m['T'] = "Thr"; m['W'] = "Trp";
        m['Y'] = "Tyr"; m['V'] = "Val"; m['*'] = "Ter";
        m['U'] = "Sec"; m['O'] = "Pyl"; m['X'] = "Xaa";
        return m;
    }();

    auto it = aa_map.find(one_letter);
    if (it != aa_map.end()) {
        return it->second;
    }
    return "Xaa";
}

/**
 * Parse HGVS type from notation
 */
inline HGVSType parse_hgvs_type(const std::string& type_str) {
    if (type_str == "g") return HGVSType::GENOMIC;
    if (type_str == "c") return HGVSType::CODING;
    if (type_str == "n") return HGVSType::NONCODING;
    if (type_str == "p") return HGVSType::PROTEIN;
    if (type_str == "m") return HGVSType::MITOCHONDRIAL;
    if (type_str == "r") return HGVSType::RNA;
    return HGVSType::UNKNOWN;
}

/**
 * Get HGVS type string
 */
inline std::string hgvs_type_to_string(HGVSType type) {
    if (type == HGVSType::GENOMIC) return "g";
    if (type == HGVSType::CODING) return "c";
    if (type == HGVSType::NONCODING) return "n";
    if (type == HGVSType::PROTEIN) return "p";
    if (type == HGVSType::MITOCHONDRIAL) return "m";
    if (type == HGVSType::RNA) return "r";
    return "?";
}

/**
 * Parse genomic HGVS notation (g.)
 * Examples:
 *   NC_000007.14:g.140753336A>T
 *   chr7:g.140753336A>T
 *   7:g.140753336A>T
 */
inline HGVSParseResult parse_genomic_hgvs(const std::string& reference, const std::string& change) {
    HGVSParseResult result;
    result.reference_id = reference;
    result.hgvs_type = HGVSType::GENOMIC;

    // Parse substitution: 140753336A>T
    static const std::regex sub_regex("(\\d+)([ACGTacgt])>([ACGTacgt])");
    std::smatch match;

    if (std::regex_match(change, match, sub_regex)) {
        result.variant_type = HGVSVariantType::SUBSTITUTION;
        result.start_pos = std::stoi(match[1].str());
        result.end_pos = result.start_pos;
        result.ref_allele = match[2].str();
        result.alt_allele = match[3].str();

        // Convert to uppercase
        for (size_t i = 0; i < result.ref_allele.size(); ++i) {
            result.ref_allele[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(result.ref_allele[i])));
        }
        for (size_t i = 0; i < result.alt_allele.size(); ++i) {
            result.alt_allele[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(result.alt_allele[i])));
        }

        result.valid = true;
        return result;
    }

    // Parse deletion: 140753336_140753340del or 140753336del
    static const std::regex del_regex("(\\d+)(?:_(\\d+))?del([ACGTacgt]*)");
    if (std::regex_match(change, match, del_regex)) {
        result.variant_type = HGVSVariantType::DELETION;
        result.start_pos = std::stoi(match[1].str());
        if (match[2].matched) {
            result.end_pos = std::stoi(match[2].str());
        } else {
            result.end_pos = result.start_pos;
        }
        result.ref_allele = match[3].str();
        result.alt_allele = "";
        result.valid = true;
        return result;
    }

    // Parse insertion: 140753336_140753337insACGT
    static const std::regex ins_regex("(\\d+)_(\\d+)ins([ACGTacgt]+)");
    if (std::regex_match(change, match, ins_regex)) {
        result.variant_type = HGVSVariantType::INSERTION;
        result.start_pos = std::stoi(match[1].str());
        result.end_pos = std::stoi(match[2].str());
        result.ref_allele = "";
        result.alt_allele = match[3].str();
        result.valid = true;
        return result;
    }

    // Parse duplication: 140753336_140753340dup
    static const std::regex dup_regex("(\\d+)(?:_(\\d+))?dup([ACGTacgt]*)");
    if (std::regex_match(change, match, dup_regex)) {
        result.variant_type = HGVSVariantType::DUPLICATION;
        result.start_pos = std::stoi(match[1].str());
        if (match[2].matched) {
            result.end_pos = std::stoi(match[2].str());
        } else {
            result.end_pos = result.start_pos;
        }
        result.ref_allele = match[3].str();
        result.valid = true;
        return result;
    }

    // Parse delins: 140753336_140753340delinsACGT
    static const std::regex delins_regex("(\\d+)(?:_(\\d+))?delins([ACGTacgt]+)");
    if (std::regex_match(change, match, delins_regex)) {
        result.variant_type = HGVSVariantType::DELINS;
        result.start_pos = std::stoi(match[1].str());
        if (match[2].matched) {
            result.end_pos = std::stoi(match[2].str());
        } else {
            result.end_pos = result.start_pos;
        }
        result.alt_allele = match[3].str();
        result.valid = true;
        return result;
    }

    result.error_message = "Unrecognized genomic HGVS pattern: " + change;
    return result;
}

/**
 * Parse coding HGVS notation (c.)
 * Examples:
 *   ENST00000366667:c.803C>T
 *   NM_000546.5:c.215C>G
 *   c.123+5G>A (intronic)
 *   c.123-10del
 */
inline HGVSParseResult parse_coding_hgvs(const std::string& reference, const std::string& change) {
    HGVSParseResult result;
    result.reference_id = reference;
    result.hgvs_type = HGVSType::CODING;

    // Parse substitution with optional intronic offset: 803C>T or 123+5G>A
    static const std::regex sub_regex("(-?\\*?\\d+)([+-]\\d+)?([ACGTacgt])>([ACGTacgt])");
    std::smatch match;

    if (std::regex_match(change, match, sub_regex)) {
        result.variant_type = HGVSVariantType::SUBSTITUTION;

        std::string pos_str = match[1].str();
        // Handle UTR positions (*123 for 3'UTR, -123 for 5'UTR)
        if (!pos_str.empty() && pos_str[0] == '*') {
            result.start_pos = std::stoi(pos_str.substr(1));
            // Mark as 3'UTR position (positive value with flag)
        } else {
            result.start_pos = std::stoi(pos_str);
        }
        result.end_pos = result.start_pos;

        if (match[2].matched) {
            result.intron_offset = std::stoi(match[2].str());
        }

        result.ref_allele = match[3].str();
        result.alt_allele = match[4].str();

        // Convert to uppercase
        for (size_t i = 0; i < result.ref_allele.size(); ++i) {
            result.ref_allele[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(result.ref_allele[i])));
        }
        for (size_t i = 0; i < result.alt_allele.size(); ++i) {
            result.alt_allele[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(result.alt_allele[i])));
        }

        result.valid = true;
        return result;
    }

    // Parse deletion with optional intronic offset
    static const std::regex del_regex("(-?\\*?\\d+)([+-]\\d+)?(?:_(-?\\*?\\d+)([+-]\\d+)?)?del([ACGTacgt]*)");
    if (std::regex_match(change, match, del_regex)) {
        result.variant_type = HGVSVariantType::DELETION;

        std::string pos_str = match[1].str();
        if (!pos_str.empty() && pos_str[0] == '*') {
            result.start_pos = std::stoi(pos_str.substr(1));
        } else {
            result.start_pos = std::stoi(pos_str);
        }

        if (match[2].matched) {
            result.intron_offset = std::stoi(match[2].str());
        }

        if (match[3].matched) {
            std::string end_str = match[3].str();
            if (!end_str.empty() && end_str[0] == '*') {
                result.end_pos = std::stoi(end_str.substr(1));
            } else {
                result.end_pos = std::stoi(end_str);
            }
        } else {
            result.end_pos = result.start_pos;
        }

        result.ref_allele = match[5].str();
        result.alt_allele = "";
        result.valid = true;
        return result;
    }

    // Parse insertion
    static const std::regex ins_regex("(-?\\*?\\d+)([+-]\\d+)?_(-?\\*?\\d+)([+-]\\d+)?ins([ACGTacgt]+)");
    if (std::regex_match(change, match, ins_regex)) {
        result.variant_type = HGVSVariantType::INSERTION;

        std::string pos_str = match[1].str();
        if (!pos_str.empty() && pos_str[0] == '*') {
            result.start_pos = std::stoi(pos_str.substr(1));
        } else {
            result.start_pos = std::stoi(pos_str);
        }

        if (match[2].matched) {
            result.intron_offset = std::stoi(match[2].str());
        }

        std::string end_str = match[3].str();
        if (!end_str.empty() && end_str[0] == '*') {
            result.end_pos = std::stoi(end_str.substr(1));
        } else {
            result.end_pos = std::stoi(end_str);
        }

        result.ref_allele = "";
        result.alt_allele = match[5].str();
        result.valid = true;
        return result;
    }

    // Parse duplication
    static const std::regex dup_regex("(-?\\*?\\d+)([+-]\\d+)?(?:_(-?\\*?\\d+)([+-]\\d+)?)?dup([ACGTacgt]*)");
    if (std::regex_match(change, match, dup_regex)) {
        result.variant_type = HGVSVariantType::DUPLICATION;

        std::string pos_str = match[1].str();
        if (!pos_str.empty() && pos_str[0] == '*') {
            result.start_pos = std::stoi(pos_str.substr(1));
        } else {
            result.start_pos = std::stoi(pos_str);
        }

        if (match[2].matched) {
            result.intron_offset = std::stoi(match[2].str());
        }

        if (match[3].matched) {
            std::string end_str = match[3].str();
            if (!end_str.empty() && end_str[0] == '*') {
                result.end_pos = std::stoi(end_str.substr(1));
            } else {
                result.end_pos = std::stoi(end_str);
            }
        } else {
            result.end_pos = result.start_pos;
        }

        result.ref_allele = match[5].str();
        result.valid = true;
        return result;
    }

    result.error_message = "Unrecognized coding HGVS pattern: " + change;
    return result;
}

/**
 * Parse protein HGVS notation (p.)
 * Examples:
 *   p.Val600Glu
 *   p.V600E
 *   p.Arg234Ter (nonsense)
 *   p.Gly12_Gly13insVal (insertion)
 *   p.Met1? (start lost)
 */
inline HGVSParseResult parse_protein_hgvs(const std::string& reference, const std::string& change) {
    HGVSParseResult result;
    result.reference_id = reference;
    result.hgvs_type = HGVSType::PROTEIN;

    // Parse missense: Val600Glu or V600E
    // Three-letter: ([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})
    static const std::regex missense_3letter("([A-Z][a-z]{2})(\\d+)([A-Z][a-z]{2}|\\?)");
    std::smatch match;

    if (std::regex_match(change, match, missense_3letter)) {
        result.variant_type = HGVSVariantType::SUBSTITUTION;
        result.ref_aa = match[1].str();
        result.protein_pos = std::stoi(match[2].str());
        result.alt_aa = match[3].str();
        result.valid = true;
        return result;
    }

    // One-letter: ([A-Z*])(\d+)([A-Z*?])
    static const std::regex missense_1letter("([A-Z*])(\\d+)([A-Z*?])");
    if (std::regex_match(change, match, missense_1letter)) {
        result.variant_type = HGVSVariantType::SUBSTITUTION;
        result.ref_aa = aa_one_to_three(match[1].str()[0]);
        result.protein_pos = std::stoi(match[2].str());
        result.alt_aa = match[3].str() == "?" ? "?" : aa_one_to_three(match[3].str()[0]);
        result.valid = true;
        return result;
    }

    // Parse frameshift: Val600fs or Val600GlufsTer12
    static const std::regex fs_regex("([A-Z][a-z]{2}|[A-Z*])(\\d+)(?:([A-Z][a-z]{2}|[A-Z*]))?fs(?:Ter|\\*)?(\\d+)?");
    if (std::regex_match(change, match, fs_regex)) {
        result.variant_type = HGVSVariantType::DELINS;  // Frameshift is complex
        if (match[1].str().size() == 1) {
            result.ref_aa = aa_one_to_three(match[1].str()[0]);
        } else {
            result.ref_aa = match[1].str();
        }
        result.protein_pos = std::stoi(match[2].str());
        result.alt_aa = "fs";
        result.valid = true;
        return result;
    }

    // Parse deletion: Gly12del or Gly12_Val14del
    static const std::regex del_regex("([A-Z][a-z]{2}|[A-Z*])(\\d+)(?:_([A-Z][a-z]{2}|[A-Z*])(\\d+))?del");
    if (std::regex_match(change, match, del_regex)) {
        result.variant_type = HGVSVariantType::DELETION;
        if (match[1].str().size() == 1) {
            result.ref_aa = aa_one_to_three(match[1].str()[0]);
        } else {
            result.ref_aa = match[1].str();
        }
        result.protein_pos = std::stoi(match[2].str());
        result.start_pos = result.protein_pos;
        if (match[4].matched) {
            result.end_pos = std::stoi(match[4].str());
        } else {
            result.end_pos = result.start_pos;
        }
        result.valid = true;
        return result;
    }

    // Parse insertion: Gly12_Ala13insVal
    static const std::regex ins_regex("([A-Z][a-z]{2}|[A-Z*])(\\d+)_([A-Z][a-z]{2}|[A-Z*])(\\d+)ins([A-Z][a-z]{2}|[A-Z*])+");
    if (std::regex_search(change, match, ins_regex)) {
        result.variant_type = HGVSVariantType::INSERTION;
        if (match[1].str().size() == 1) {
            result.ref_aa = aa_one_to_three(match[1].str()[0]);
        } else {
            result.ref_aa = match[1].str();
        }
        result.protein_pos = std::stoi(match[2].str());
        result.start_pos = result.protein_pos;
        result.end_pos = std::stoi(match[4].str());
        // Extract inserted AAs
        size_t ins_pos = change.find("ins");
        if (ins_pos != std::string::npos) {
            result.alt_aa = change.substr(ins_pos + 3);
        }
        result.valid = true;
        return result;
    }

    // Parse synonymous: Val600= or V600=
    static const std::regex syn_regex("([A-Z][a-z]{2}|[A-Z*])(\\d+)=");
    if (std::regex_match(change, match, syn_regex)) {
        result.variant_type = HGVSVariantType::SUBSTITUTION;
        if (match[1].str().size() == 1) {
            result.ref_aa = aa_one_to_three(match[1].str()[0]);
        } else {
            result.ref_aa = match[1].str();
        }
        result.protein_pos = std::stoi(match[2].str());
        result.alt_aa = result.ref_aa;  // Synonymous
        result.valid = true;
        return result;
    }

    result.error_message = "Unrecognized protein HGVS pattern: " + change;
    return result;
}

/**
 * Main HGVS parser function
 * Parses full HGVS notation: REFERENCE:TYPE.CHANGE
 * Examples:
 *   ENST00000366667:c.803C>T
 *   NC_000007.14:g.140753336A>T
 *   BRAF:p.Val600Glu
 */
inline HGVSParseResult parse_hgvs(const std::string& hgvs) {
    HGVSParseResult result;

    // Find the colon separator
    size_t colon_pos = hgvs.find(':');
    if (colon_pos == std::string::npos) {
        result.error_message = "Invalid HGVS format: missing colon separator";
        return result;
    }

    std::string reference = hgvs.substr(0, colon_pos);
    std::string notation = hgvs.substr(colon_pos + 1);

    if (notation.size() < 3) {
        result.error_message = "Invalid HGVS format: notation too short";
        return result;
    }

    // Parse type (c., g., p., n., m., r.)
    char type_char = notation[0];
    if (notation[1] != '.') {
        result.error_message = "Invalid HGVS format: expected '.' after type";
        return result;
    }

    std::string change = notation.substr(2);
    HGVSType hgvs_type = parse_hgvs_type(std::string(1, type_char));

    try {
        if (hgvs_type == HGVSType::GENOMIC || hgvs_type == HGVSType::MITOCHONDRIAL) {
            result = parse_genomic_hgvs(reference, change);
        } else if (hgvs_type == HGVSType::CODING || hgvs_type == HGVSType::NONCODING) {
            result = parse_coding_hgvs(reference, change);
        } else if (hgvs_type == HGVSType::PROTEIN) {
            result = parse_protein_hgvs(reference, change);
        } else if (hgvs_type == HGVSType::RNA) {
            // RNA notation is similar to coding
            result = parse_coding_hgvs(reference, change);
            result.hgvs_type = HGVSType::RNA;
        } else {
            result.error_message = "Unsupported HGVS type: " + std::string(1, type_char);
            return result;
        }
    } catch (const std::out_of_range&) {
        result.error_message = "HGVS position out of range: " + hgvs;
        result.valid = false;
        return result;
    } catch (const std::invalid_argument&) {
        result.error_message = "Invalid HGVS position format: " + hgvs;
        result.valid = false;
        return result;
    }

    result.reference_id = reference;
    return result;
}

/**
 * Check if a string looks like HGVS notation
 */
inline bool is_hgvs_notation(const std::string& input) {
    // Must contain a colon
    size_t colon_pos = input.find(':');
    if (colon_pos == std::string::npos || colon_pos == 0) {
        return false;
    }

    // Check for type indicator after colon
    if (input.size() > colon_pos + 2) {
        char type_char = input[colon_pos + 1];
        char dot_char = input[colon_pos + 2];

        if (dot_char == '.') {
            return (type_char == 'c' || type_char == 'g' || type_char == 'p' ||
                    type_char == 'n' || type_char == 'm' || type_char == 'r');
        }
    }

    return false;
}

/**
 * RefSeq chromosome accession to chromosome name mapping
 */
inline std::string refseq_to_chromosome(const std::string& refseq) {
    static const std::unordered_map<std::string, std::string> refseq_map = [] {
        std::unordered_map<std::string, std::string> m;
        // GRCh38 chromosome accessions
        m["NC_000001.11"] = "1";
        m["NC_000002.12"] = "2";
        m["NC_000003.12"] = "3";
        m["NC_000004.12"] = "4";
        m["NC_000005.10"] = "5";
        m["NC_000006.12"] = "6";
        m["NC_000007.14"] = "7";
        m["NC_000008.11"] = "8";
        m["NC_000009.12"] = "9";
        m["NC_000010.11"] = "10";
        m["NC_000011.10"] = "11";
        m["NC_000012.12"] = "12";
        m["NC_000013.11"] = "13";
        m["NC_000014.9"] = "14";
        m["NC_000015.10"] = "15";
        m["NC_000016.10"] = "16";
        m["NC_000017.11"] = "17";
        m["NC_000018.10"] = "18";
        m["NC_000019.10"] = "19";
        m["NC_000020.11"] = "20";
        m["NC_000021.9"] = "21";
        m["NC_000022.11"] = "22";
        m["NC_000023.11"] = "X";
        m["NC_000024.10"] = "Y";
        m["NC_012920.1"] = "MT";

        // GRCh37 chromosome accessions
        m["NC_000001.10"] = "1";
        m["NC_000002.11"] = "2";
        m["NC_000003.11"] = "3";
        m["NC_000004.11"] = "4";
        m["NC_000005.9"] = "5";
        m["NC_000006.11"] = "6";
        m["NC_000007.13"] = "7";
        m["NC_000008.10"] = "8";
        m["NC_000009.11"] = "9";
        m["NC_000010.10"] = "10";
        m["NC_000011.9"] = "11";
        m["NC_000012.11"] = "12";
        m["NC_000013.10"] = "13";
        m["NC_000014.8"] = "14";
        m["NC_000015.9"] = "15";
        m["NC_000016.9"] = "16";
        m["NC_000017.10"] = "17";
        m["NC_000018.9"] = "18";
        m["NC_000019.9"] = "19";
        m["NC_000020.10"] = "20";
        m["NC_000021.8"] = "21";
        m["NC_000022.10"] = "22";
        m["NC_000023.10"] = "X";
        m["NC_000024.9"] = "Y";
        return m;
    }();

    // Check for exact match
    auto it = refseq_map.find(refseq);
    if (it != refseq_map.end()) {
        return it->second;
    }

    // Check for versioned accession (match without version)
    size_t dot_pos = refseq.find('.');
    if (dot_pos != std::string::npos) {
        std::string base = refseq.substr(0, dot_pos);
        for (auto map_it = refseq_map.begin(); map_it != refseq_map.end(); ++map_it) {
            size_t map_dot = map_it->first.find('.');
            if (map_dot != std::string::npos && map_it->first.substr(0, map_dot) == base) {
                return map_it->second;
            }
        }
    }

    // Check if already a chromosome name
    if (refseq.size() <= 2 ||
        refseq == "MT" || refseq == "X" || refseq == "Y" ||
        (refseq.size() > 3 && refseq.substr(0, 3) == "chr")) {
        return refseq;
    }

    return "";  // Unknown
}

/**
 * Shared chromosome-to-RefSeq accession mapping (GRCh38)
 */
inline const std::unordered_map<std::string, std::string>& get_chrom_to_refseq() {
    static const std::unordered_map<std::string, std::string> chrom_to_refseq = [] {
        std::unordered_map<std::string, std::string> m;
        m["1"] = "NC_000001.11";
        m["2"] = "NC_000002.12";
        m["3"] = "NC_000003.12";
        m["4"] = "NC_000004.12";
        m["5"] = "NC_000005.10";
        m["6"] = "NC_000006.12";
        m["7"] = "NC_000007.14";
        m["8"] = "NC_000008.11";
        m["9"] = "NC_000009.12";
        m["10"] = "NC_000010.11";
        m["11"] = "NC_000011.10";
        m["12"] = "NC_000012.12";
        m["13"] = "NC_000013.11";
        m["14"] = "NC_000014.9";
        m["15"] = "NC_000015.10";
        m["16"] = "NC_000016.10";
        m["17"] = "NC_000017.11";
        m["18"] = "NC_000018.10";
        m["19"] = "NC_000019.10";
        m["20"] = "NC_000020.11";
        m["21"] = "NC_000021.9";
        m["22"] = "NC_000022.11";
        m["X"] = "NC_000023.11";
        m["Y"] = "NC_000024.10";
        m["MT"] = "NC_012920.1";
        m["M"] = "NC_012920.1";  // chrM -> M after prefix strip
        return m;
    }();
    return chrom_to_refseq;
}

/**
 * Look up RefSeq accession for a chromosome name, normalizing chr prefix
 */
inline std::string chrom_to_refseq_lookup(const std::string& chrom) {
    std::string norm_chrom = chrom;
    if (norm_chrom.size() > 3 && norm_chrom.substr(0, 3) == "chr") {
        norm_chrom = norm_chrom.substr(3);
    }

    const auto& chrom_to_refseq = get_chrom_to_refseq();
    auto it = chrom_to_refseq.find(norm_chrom);
    if (it != chrom_to_refseq.end()) {
        return it->second;
    }
    return chrom;  // Use as-is
}

/**
 * Generate genomic HGVS notation (g.)
 */
inline std::string generate_hgvsg(const std::string& chrom, int pos,
                                   const std::string& ref, const std::string& alt) {
    std::string result;

    // Get RefSeq accession via shared helper
    std::string refseq = chrom_to_refseq_lookup(chrom);

    result = refseq + ":g.";

    if (ref.size() == 1 && alt.size() == 1) {
        // Substitution
        result += std::to_string(pos) + ref + ">" + alt;
    } else if (alt.empty() || (ref.size() > alt.size() && alt.size() <= 1)) {
        // Deletion
        int del_start = pos;
        int del_end = pos + static_cast<int>(ref.size()) - 1;
        // Strip VCF anchor base if present
        if (!alt.empty() && alt.size() == 1 && ref.size() > 1 && ref[0] == alt[0]) {
            del_start = pos + 1;
        }
        if (del_start == del_end) {
            result += std::to_string(del_start) + "del";
        } else {
            result += std::to_string(del_start) + "_" + std::to_string(del_end) + "del";
        }
    } else if (ref.empty() || (alt.size() > ref.size() && ref.size() <= 1)) {
        // Insertion - check for duplication heuristic
        std::string inserted = alt;
        if (!ref.empty() && ref.size() == 1) {
            inserted = alt.substr(1);
        }
        // Heuristic duplication detection: if the inserted sequence is a single base
        // and matches the anchor/ref base, it's likely a duplication. For multi-base
        // insertions, we cannot determine duplication without reference context, so
        // we use standard ins notation.
        bool is_dup = false;
        if (!ref.empty() && ref.size() == 1 && inserted == ref) {
            // Single-base duplication (e.g., ref=A, alt=AA -> dup of A at pos)
            is_dup = true;
        }
        if (is_dup) {
            int ins_len = static_cast<int>(inserted.size());
            if (ins_len == 1) {
                result += std::to_string(pos) + "dup";
            } else {
                result += std::to_string(pos) + "_" + std::to_string(pos + ins_len - 1) + "dup";
            }
        } else {
            result += std::to_string(pos) + "_" + std::to_string(pos + 1) + "ins" + inserted;
        }
    } else {
        // Complex (delins)
        int delins_start = pos;
        int delins_end = pos + static_cast<int>(ref.size()) - 1;
        std::string delins_alt = alt;
        // Strip VCF anchor base if present
        if (ref.size() > 1 && alt.size() > 1 && ref[0] == alt[0]) {
            delins_start = pos + 1;
            delins_alt = alt.substr(1);
        }
        if (delins_start == delins_end) {
            result += std::to_string(delins_start) + "delins" + delins_alt;
        } else {
            result += std::to_string(delins_start) + "_" + std::to_string(delins_end) + "delins" + delins_alt;
        }
    }

    return result;
}

/**
 * Generate SPDI format
 * Format: SEQ_ID:POSITION:DELETED_SEQUENCE:INSERTED_SEQUENCE
 */
inline std::string generate_spdi(const std::string& chrom, int pos,
                                  const std::string& ref, const std::string& alt) {
    // SPDI uses 0-based coordinates
    int spdi_pos = pos - 1;

    // Get RefSeq accession via shared helper
    std::string refseq = chrom_to_refseq_lookup(chrom);

    return refseq + ":" + std::to_string(spdi_pos) + ":" + ref + ":" + alt;
}

/**
 * Parse SPDI notation
 * Format: SEQ_ID:POSITION:DELETED_SEQUENCE:INSERTED_SEQUENCE
 * Note: SPDI uses 0-based positions
 * Examples:
 *   NC_000007.14:140753335:A:T
 *   7:140753335:A:T
 */
struct SPDIParseResult {
    bool valid = false;
    std::string error_message;
    std::string chromosome;
    int position = 0;      // 1-based (converted from 0-based SPDI)
    std::string ref_allele;
    std::string alt_allele;
};

inline SPDIParseResult parse_spdi(const std::string& spdi) {
    SPDIParseResult result;

    // Split on colons - expect exactly 4 parts
    std::vector<std::string> parts;
    size_t start = 0;
    size_t pos = spdi.find(':');
    while (pos != std::string::npos) {
        parts.push_back(spdi.substr(start, pos - start));
        start = pos + 1;
        pos = spdi.find(':', start);
    }
    parts.push_back(spdi.substr(start));

    if (parts.size() != 4) {
        result.error_message = "Invalid SPDI format: expected 4 colon-separated fields";
        return result;
    }

    // Map RefSeq accession to chromosome name
    std::string chrom = refseq_to_chromosome(parts[0]);
    if (chrom.empty()) {
        chrom = parts[0];
    }

    // Normalize chr prefix
    if (chrom.size() > 3 && chrom.substr(0, 3) == "chr") {
        chrom = chrom.substr(3);
    }

    result.chromosome = chrom;

    // Convert 0-based SPDI position to 1-based
    try {
        result.position = std::stoi(parts[1]) + 1;
    } catch (...) {
        result.error_message = "Invalid SPDI position: " + parts[1];
        return result;
    }

    result.ref_allele = parts[2];
    result.alt_allele = parts[3];

    // Convert to uppercase
    for (size_t i = 0; i < result.ref_allele.size(); ++i) {
        result.ref_allele[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(result.ref_allele[i])));
    }
    for (size_t i = 0; i < result.alt_allele.size(); ++i) {
        result.alt_allele[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(result.alt_allele[i])));
    }

    result.valid = true;
    return result;
}

/**
 * Check if a string looks like SPDI notation
 * SPDI has exactly 3 colons and the second field is numeric
 */
inline bool is_spdi_notation(const std::string& input) {
    int colon_count = 0;
    size_t first_colon = std::string::npos;
    size_t second_colon = std::string::npos;
    for (size_t i = 0; i < input.size(); ++i) {
        if (input[i] == ':') {
            colon_count++;
            if (colon_count == 1) first_colon = i;
            if (colon_count == 2) second_colon = i;
        }
    }
    if (colon_count != 3) return false;

    // Check if the second field is numeric (the position)
    if (first_colon == std::string::npos || second_colon == std::string::npos) return false;
    std::string pos_str = input.substr(first_colon + 1, second_colon - first_colon - 1);
    if (pos_str.empty()) return false;
    for (size_t i = 0; i < pos_str.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(pos_str[i]))) return false;
    }

    // Distinguish from CHR:POS:REF:ALT by checking if first field looks like a RefSeq accession
    // or if the position field is very long (SPDI uses 0-based)
    std::string first_field = input.substr(0, first_colon);
    // If it starts with NC_ or NT_ or NW_, it's SPDI
    if (first_field.size() >= 3 &&
        (first_field.substr(0, 3) == "NC_" || first_field.substr(0, 3) == "NT_" ||
         first_field.substr(0, 3) == "NW_")) {
        return true;
    }
    return false;
}

/**
 * Parse Ensembl default format line
 * Format: CHR START END ALLELE STRAND
 * Examples:
 *   1 100 100 A/T +
 *   7 140753336 140753336 A/T 1
 *   X 100 102 ATG/- +
 */
struct EnsemblFormatResult {
    bool valid = false;
    std::string error_message;
    std::string chromosome;
    int position = 0;      // 1-based
    std::string ref_allele;
    std::string alt_allele;
};

inline EnsemblFormatResult parse_ensembl_format(const std::string& line) {
    EnsemblFormatResult result;

    std::istringstream iss(line);
    std::string chrom, allele, strand;
    int start_pos, end_pos;

    if (!(iss >> chrom >> start_pos >> end_pos >> allele)) {
        result.error_message = "Invalid Ensembl format: need at least 4 fields";
        return result;
    }

    // Strand is optional
    if (!(iss >> strand)) {
        strand = "+";
    }

    // Parse allele: REF/ALT
    size_t slash_pos = allele.find('/');
    if (slash_pos == std::string::npos) {
        result.error_message = "Invalid allele format: expected REF/ALT";
        return result;
    }

    result.chromosome = chrom;
    result.position = start_pos;
    result.ref_allele = allele.substr(0, slash_pos);
    result.alt_allele = allele.substr(slash_pos + 1);

    // Handle - for deletions/insertions
    if (result.ref_allele == "-") result.ref_allele = "";
    if (result.alt_allele == "-") result.alt_allele = "";

    // Convert to uppercase
    for (size_t i = 0; i < result.ref_allele.size(); ++i) {
        result.ref_allele[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(result.ref_allele[i])));
    }
    for (size_t i = 0; i < result.alt_allele.size(); ++i) {
        result.alt_allele[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(result.alt_allele[i])));
    }

    // Handle strand
    bool negative_strand = (strand == "-" || strand == "-1");
    if (negative_strand) {
        // Complement and reverse the alleles
        for (int si = 0; si < 2; ++si) {
            std::string& seq = (si == 0) ? result.ref_allele : result.alt_allele;
            for (size_t i = 0; i < seq.size(); ++i) {
                switch (seq[i]) {
                    case 'A': seq[i] = 'T'; break;
                    case 'T': seq[i] = 'A'; break;
                    case 'C': seq[i] = 'G'; break;
                    case 'G': seq[i] = 'C'; break;
                    default: break;
                }
            }
            std::reverse(seq.begin(), seq.end());
        }
    }

    result.valid = true;
    return result;
}

/**
 * Check if a line looks like Ensembl default format
 * 5 whitespace-separated fields where fields 2 and 3 are numeric
 * and field 4 contains a /
 */
inline bool is_ensembl_format(const std::string& line) {
    std::istringstream iss(line);
    std::string f1, f2, f3, f4;
    if (!(iss >> f1 >> f2 >> f3 >> f4)) return false;

    // f2 and f3 must be numeric
    for (size_t i = 0; i < f2.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(f2[i]))) return false;
    }
    for (size_t i = 0; i < f3.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(f3[i]))) return false;
    }

    // f4 must contain a /
    if (f4.find('/') == std::string::npos) return false;

    return true;
}

/**
 * REST-style region format: CHR:START-END:STRAND/ALLELE
 * e.g., "7:140753336-140753336:1/T" or "1:12345-12345:-1/A"
 */
struct RESTRegionResult {
    bool valid = false;
    std::string error_message;
    std::string chromosome;
    int position = 0;
    int end_position = 0;
    int strand = 1;         // 1 or -1
    std::string ref_allele;
    std::string alt_allele;
};

inline bool is_rest_region_format(const std::string& input) {
    // Pattern: CHR:START-END:STRAND/ALLELE
    size_t first_colon = input.find(':');
    if (first_colon == std::string::npos) return false;

    size_t dash = input.find('-', first_colon + 1);
    if (dash == std::string::npos) return false;

    size_t second_colon = input.find(':', dash + 1);
    if (second_colon == std::string::npos) return false;

    size_t slash = input.find('/', second_colon + 1);
    if (slash == std::string::npos) return false;

    return true;
}

inline RESTRegionResult parse_rest_region(const std::string& input) {
    RESTRegionResult result;

    size_t first_colon = input.find(':');
    size_t dash = input.find('-', first_colon + 1);
    size_t second_colon = input.find(':', dash + 1);
    size_t slash = input.find('/', second_colon + 1);

    if (first_colon == std::string::npos || dash == std::string::npos ||
        second_colon == std::string::npos || slash == std::string::npos) {
        result.error_message = "Invalid REST region format";
        return result;
    }

    try {
        result.chromosome = input.substr(0, first_colon);
        result.position = std::stoi(input.substr(first_colon + 1, dash - first_colon - 1));
        result.end_position = std::stoi(input.substr(dash + 1, second_colon - dash - 1));
        result.strand = std::stoi(input.substr(second_colon + 1, slash - second_colon - 1));
        result.alt_allele = input.substr(slash + 1);

        // For single-base substitutions, ref is inferred from position
        // For insertions (start > end), ref is "-"
        if (result.position > result.end_position) {
            result.ref_allele = "-";
        } else {
            result.ref_allele = "-";  // Will be filled by reference lookup
        }

        result.valid = true;
    } catch (const std::exception& e) {
        result.error_message = std::string("Failed to parse REST region: ") + e.what();
    }

    return result;
}

} // namespace vep

#endif // HGVS_PARSER_HPP
