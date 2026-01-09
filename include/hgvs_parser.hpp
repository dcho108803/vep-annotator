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
#include <vector>
#include <regex>
#include <stdexcept>
#include <cctype>

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
    static std::map<std::string, std::string> aa_map;
    if (aa_map.empty()) {
        aa_map["Ala"] = "A"; aa_map["ALA"] = "A";
        aa_map["Arg"] = "R"; aa_map["ARG"] = "R";
        aa_map["Asn"] = "N"; aa_map["ASN"] = "N";
        aa_map["Asp"] = "D"; aa_map["ASP"] = "D";
        aa_map["Cys"] = "C"; aa_map["CYS"] = "C";
        aa_map["Gln"] = "Q"; aa_map["GLN"] = "Q";
        aa_map["Glu"] = "E"; aa_map["GLU"] = "E";
        aa_map["Gly"] = "G"; aa_map["GLY"] = "G";
        aa_map["His"] = "H"; aa_map["HIS"] = "H";
        aa_map["Ile"] = "I"; aa_map["ILE"] = "I";
        aa_map["Leu"] = "L"; aa_map["LEU"] = "L";
        aa_map["Lys"] = "K"; aa_map["LYS"] = "K";
        aa_map["Met"] = "M"; aa_map["MET"] = "M";
        aa_map["Phe"] = "F"; aa_map["PHE"] = "F";
        aa_map["Pro"] = "P"; aa_map["PRO"] = "P";
        aa_map["Ser"] = "S"; aa_map["SER"] = "S";
        aa_map["Thr"] = "T"; aa_map["THR"] = "T";
        aa_map["Trp"] = "W"; aa_map["TRP"] = "W";
        aa_map["Tyr"] = "Y"; aa_map["TYR"] = "Y";
        aa_map["Val"] = "V"; aa_map["VAL"] = "V";
        aa_map["Ter"] = "*"; aa_map["TER"] = "*"; aa_map["Stop"] = "*";
        aa_map["Sec"] = "U"; aa_map["SEC"] = "U";  // Selenocysteine
        aa_map["Pyl"] = "O"; aa_map["PYL"] = "O";  // Pyrrolysine
    }

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
    static std::map<char, std::string> aa_map;
    if (aa_map.empty()) {
        aa_map['A'] = "Ala"; aa_map['R'] = "Arg"; aa_map['N'] = "Asn";
        aa_map['D'] = "Asp"; aa_map['C'] = "Cys"; aa_map['Q'] = "Gln";
        aa_map['E'] = "Glu"; aa_map['G'] = "Gly"; aa_map['H'] = "His";
        aa_map['I'] = "Ile"; aa_map['L'] = "Leu"; aa_map['K'] = "Lys";
        aa_map['M'] = "Met"; aa_map['F'] = "Phe"; aa_map['P'] = "Pro";
        aa_map['S'] = "Ser"; aa_map['T'] = "Thr"; aa_map['W'] = "Trp";
        aa_map['Y'] = "Tyr"; aa_map['V'] = "Val"; aa_map['*'] = "Ter";
        aa_map['U'] = "Sec"; aa_map['O'] = "Pyl"; aa_map['X'] = "Xaa";
    }

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
    std::regex sub_regex("(\\d+)([ACGTacgt])>([ACGTacgt])");
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
    std::regex del_regex("(\\d+)(?:_(\\d+))?del([ACGTacgt]*)");
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
    std::regex ins_regex("(\\d+)_(\\d+)ins([ACGTacgt]+)");
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
    std::regex dup_regex("(\\d+)(?:_(\\d+))?dup([ACGTacgt]*)");
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
    std::regex delins_regex("(\\d+)(?:_(\\d+))?delins([ACGTacgt]+)");
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
    std::regex sub_regex("(-?\\*?\\d+)([+-]\\d+)?([ACGTacgt])>([ACGTacgt])");
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
    std::regex del_regex("(-?\\*?\\d+)([+-]\\d+)?(?:_(-?\\*?\\d+)([+-]\\d+)?)?del([ACGTacgt]*)");
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
    std::regex ins_regex("(-?\\*?\\d+)([+-]\\d+)?_(-?\\*?\\d+)([+-]\\d+)?ins([ACGTacgt]+)");
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
    std::regex dup_regex("(-?\\*?\\d+)([+-]\\d+)?(?:_(-?\\*?\\d+)([+-]\\d+)?)?dup([ACGTacgt]*)");
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
    std::regex missense_3letter("([A-Z][a-z]{2})(\\d+)([A-Z][a-z]{2}|\\?)");
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
    std::regex missense_1letter("([A-Z*])(\\d+)([A-Z*?])");
    if (std::regex_match(change, match, missense_1letter)) {
        result.variant_type = HGVSVariantType::SUBSTITUTION;
        result.ref_aa = aa_one_to_three(match[1].str()[0]);
        result.protein_pos = std::stoi(match[2].str());
        result.alt_aa = match[3].str() == "?" ? "?" : aa_one_to_three(match[3].str()[0]);
        result.valid = true;
        return result;
    }

    // Parse frameshift: Val600fs or Val600GlufsTer12
    std::regex fs_regex("([A-Z][a-z]{2}|[A-Z*])(\\d+)(?:([A-Z][a-z]{2}|[A-Z*]))?fs(?:Ter|\\*)?(\\d+)?");
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
    std::regex del_regex("([A-Z][a-z]{2}|[A-Z*])(\\d+)(?:_([A-Z][a-z]{2}|[A-Z*])(\\d+))?del");
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
    std::regex ins_regex("([A-Z][a-z]{2}|[A-Z*])(\\d+)_([A-Z][a-z]{2}|[A-Z*])(\\d+)ins([A-Z][a-z]{2}|[A-Z*])+");
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
    std::regex syn_regex("([A-Z][a-z]{2}|[A-Z*])(\\d+)=");
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
    static std::map<std::string, std::string> refseq_map;
    if (refseq_map.empty()) {
        // GRCh38 chromosome accessions
        refseq_map["NC_000001.11"] = "1";
        refseq_map["NC_000002.12"] = "2";
        refseq_map["NC_000003.12"] = "3";
        refseq_map["NC_000004.12"] = "4";
        refseq_map["NC_000005.10"] = "5";
        refseq_map["NC_000006.12"] = "6";
        refseq_map["NC_000007.14"] = "7";
        refseq_map["NC_000008.11"] = "8";
        refseq_map["NC_000009.12"] = "9";
        refseq_map["NC_000010.11"] = "10";
        refseq_map["NC_000011.10"] = "11";
        refseq_map["NC_000012.12"] = "12";
        refseq_map["NC_000013.11"] = "13";
        refseq_map["NC_000014.9"] = "14";
        refseq_map["NC_000015.10"] = "15";
        refseq_map["NC_000016.10"] = "16";
        refseq_map["NC_000017.11"] = "17";
        refseq_map["NC_000018.10"] = "18";
        refseq_map["NC_000019.10"] = "19";
        refseq_map["NC_000020.11"] = "20";
        refseq_map["NC_000021.9"] = "21";
        refseq_map["NC_000022.11"] = "22";
        refseq_map["NC_000023.11"] = "X";
        refseq_map["NC_000024.10"] = "Y";
        refseq_map["NC_012920.1"] = "MT";

        // GRCh37 chromosome accessions
        refseq_map["NC_000001.10"] = "1";
        refseq_map["NC_000002.11"] = "2";
        refseq_map["NC_000003.11"] = "3";
        refseq_map["NC_000004.11"] = "4";
        refseq_map["NC_000005.9"] = "5";
        refseq_map["NC_000006.11"] = "6";
        refseq_map["NC_000007.13"] = "7";
        refseq_map["NC_000008.10"] = "8";
        refseq_map["NC_000009.11"] = "9";
        refseq_map["NC_000010.10"] = "10";
        refseq_map["NC_000011.9"] = "11";
        refseq_map["NC_000012.11"] = "12";
        refseq_map["NC_000013.10"] = "13";
        refseq_map["NC_000014.8"] = "14";
        refseq_map["NC_000015.9"] = "15";
        refseq_map["NC_000016.9"] = "16";
        refseq_map["NC_000017.10"] = "17";
        refseq_map["NC_000018.9"] = "18";
        refseq_map["NC_000019.9"] = "19";
        refseq_map["NC_000020.10"] = "20";
        refseq_map["NC_000021.8"] = "21";
        refseq_map["NC_000022.10"] = "22";
        refseq_map["NC_000023.10"] = "X";
        refseq_map["NC_000024.9"] = "Y";
    }

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
 * Generate genomic HGVS notation (g.)
 */
inline std::string generate_hgvsg(const std::string& chrom, int pos,
                                   const std::string& ref, const std::string& alt) {
    std::string result;

    // Map chromosome to RefSeq accession
    static std::map<std::string, std::string> chrom_to_refseq;
    if (chrom_to_refseq.empty()) {
        chrom_to_refseq["1"] = "NC_000001.11";
        chrom_to_refseq["2"] = "NC_000002.12";
        chrom_to_refseq["3"] = "NC_000003.12";
        chrom_to_refseq["4"] = "NC_000004.12";
        chrom_to_refseq["5"] = "NC_000005.10";
        chrom_to_refseq["6"] = "NC_000006.12";
        chrom_to_refseq["7"] = "NC_000007.14";
        chrom_to_refseq["8"] = "NC_000008.11";
        chrom_to_refseq["9"] = "NC_000009.12";
        chrom_to_refseq["10"] = "NC_000010.11";
        chrom_to_refseq["11"] = "NC_000011.10";
        chrom_to_refseq["12"] = "NC_000012.12";
        chrom_to_refseq["13"] = "NC_000013.11";
        chrom_to_refseq["14"] = "NC_000014.9";
        chrom_to_refseq["15"] = "NC_000015.10";
        chrom_to_refseq["16"] = "NC_000016.10";
        chrom_to_refseq["17"] = "NC_000017.11";
        chrom_to_refseq["18"] = "NC_000018.10";
        chrom_to_refseq["19"] = "NC_000019.10";
        chrom_to_refseq["20"] = "NC_000020.11";
        chrom_to_refseq["21"] = "NC_000021.9";
        chrom_to_refseq["22"] = "NC_000022.11";
        chrom_to_refseq["X"] = "NC_000023.11";
        chrom_to_refseq["Y"] = "NC_000024.10";
        chrom_to_refseq["MT"] = "NC_012920.1";
    }

    // Normalize chromosome name
    std::string norm_chrom = chrom;
    if (norm_chrom.size() > 3 && norm_chrom.substr(0, 3) == "chr") {
        norm_chrom = norm_chrom.substr(3);
    }

    // Get RefSeq accession
    std::string refseq;
    auto it = chrom_to_refseq.find(norm_chrom);
    if (it != chrom_to_refseq.end()) {
        refseq = it->second;
    } else {
        refseq = chrom;  // Use as-is
    }

    result = refseq + ":g.";

    if (ref.size() == 1 && alt.size() == 1) {
        // Substitution
        result += std::to_string(pos) + ref + ">" + alt;
    } else if (alt.empty() || (ref.size() > alt.size() && alt.size() <= 1)) {
        // Deletion
        if (ref.size() == 1) {
            result += std::to_string(pos) + "del";
        } else {
            result += std::to_string(pos) + "_" + std::to_string(pos + static_cast<int>(ref.size()) - 1) + "del";
        }
    } else if (ref.empty() || (alt.size() > ref.size() && ref.size() <= 1)) {
        // Insertion
        std::string inserted = alt;
        if (!ref.empty() && ref.size() == 1) {
            inserted = alt.substr(1);
        }
        result += std::to_string(pos) + "_" + std::to_string(pos + 1) + "ins" + inserted;
    } else {
        // Complex (delins)
        if (ref.size() == 1) {
            result += std::to_string(pos) + "delins" + alt;
        } else {
            result += std::to_string(pos) + "_" + std::to_string(pos + static_cast<int>(ref.size()) - 1) + "delins" + alt;
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

    // Map chromosome to RefSeq accession
    static std::map<std::string, std::string> chrom_to_refseq;
    if (chrom_to_refseq.empty()) {
        chrom_to_refseq["1"] = "NC_000001.11";
        chrom_to_refseq["2"] = "NC_000002.12";
        chrom_to_refseq["3"] = "NC_000003.12";
        chrom_to_refseq["4"] = "NC_000004.12";
        chrom_to_refseq["5"] = "NC_000005.10";
        chrom_to_refseq["6"] = "NC_000006.12";
        chrom_to_refseq["7"] = "NC_000007.14";
        chrom_to_refseq["8"] = "NC_000008.11";
        chrom_to_refseq["9"] = "NC_000009.12";
        chrom_to_refseq["10"] = "NC_000010.11";
        chrom_to_refseq["11"] = "NC_000011.10";
        chrom_to_refseq["12"] = "NC_000012.12";
        chrom_to_refseq["13"] = "NC_000013.11";
        chrom_to_refseq["14"] = "NC_000014.9";
        chrom_to_refseq["15"] = "NC_000015.10";
        chrom_to_refseq["16"] = "NC_000016.10";
        chrom_to_refseq["17"] = "NC_000017.11";
        chrom_to_refseq["18"] = "NC_000018.10";
        chrom_to_refseq["19"] = "NC_000019.10";
        chrom_to_refseq["20"] = "NC_000020.11";
        chrom_to_refseq["21"] = "NC_000021.9";
        chrom_to_refseq["22"] = "NC_000022.11";
        chrom_to_refseq["X"] = "NC_000023.11";
        chrom_to_refseq["Y"] = "NC_000024.10";
        chrom_to_refseq["MT"] = "NC_012920.1";
    }

    // Normalize chromosome name
    std::string norm_chrom = chrom;
    if (norm_chrom.size() > 3 && norm_chrom.substr(0, 3) == "chr") {
        norm_chrom = norm_chrom.substr(3);
    }

    // Get RefSeq accession
    std::string refseq;
    auto it = chrom_to_refseq.find(norm_chrom);
    if (it != chrom_to_refseq.end()) {
        refseq = it->second;
    } else {
        refseq = chrom;  // Use as-is
    }

    return refseq + ":" + std::to_string(spdi_pos) + ":" + ref + ":" + alt;
}

} // namespace vep

#endif // HGVS_PARSER_HPP
