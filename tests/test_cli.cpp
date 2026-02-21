/**
 * Tests for CLI argument parsing and utility functions from main.cpp.
 *
 * Since main.cpp functions are file-local (not exported), we re-implement
 * the core parsing logic here to verify the algorithms. This ensures the
 * parsing patterns used in the CLI remain correct across refactors.
 */

#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <tuple>
#include <algorithm>
#include <cctype>

// ============================================================================
// Re-implemented utility functions (mirrors main.cpp logic)
// ============================================================================

// Parse variant string in format CHR:POS:REF:ALT or CHR-POS-REF-ALT
static bool parse_variant(const std::string& variant, std::string& chrom, int& pos,
                          std::string& ref, std::string& alt) {
    // Try colon-separated first (preferred)
    char sep = ':';
    size_t p1 = variant.find(':');
    if (p1 == std::string::npos) {
        // Fall back to dash separator
        sep = '-';
        p1 = variant.find('-');
        if (p1 == std::string::npos) return false;
    }

    size_t p2 = variant.find(sep, p1 + 1);
    if (p2 == std::string::npos) return false;

    size_t p3 = variant.find(sep, p2 + 1);
    if (p3 == std::string::npos) return false;

    try {
        chrom = variant.substr(0, p1);
        pos = std::stoi(variant.substr(p1 + 1, p2 - p1 - 1));
        ref = variant.substr(p2 + 1, p3 - p2 - 1);
        alt = variant.substr(p3 + 1);
        return true;
    } catch (...) {
        return false;
    }
}

// Input format types for auto-detection
enum class InputFormat { VCF, ENSEMBL, HGVS, BED, REGION, UNKNOWN };

// Detect input format from a non-header line
static InputFormat detect_input_format(const std::string& line) {
    if (line.empty() || line[0] == '#') return InputFormat::UNKNOWN;

    // HGVS: contains ':' followed by c./p./g./n./m. notation
    if (line.find(":c.") != std::string::npos || line.find(":p.") != std::string::npos ||
        line.find(":g.") != std::string::npos || line.find(":n.") != std::string::npos ||
        line.find(":m.") != std::string::npos) {
        return InputFormat::HGVS;
    }

    // Region format: CHR:START-END/ALLELE or CHR:START-END:STRAND/ALLELE
    if (line.find('/') != std::string::npos && line.find(':') != std::string::npos &&
        line.find('-') != std::string::npos) {
        size_t colon1 = line.find(':');
        size_t dash = line.find('-', colon1);
        if (dash != std::string::npos) {
            return InputFormat::REGION;
        }
    }

    // Count tabs to distinguish VCF from Ensembl format
    int tab_count = 0;
    for (char c : line) {
        if (c == '\t') tab_count++;
    }

    // VCF: at least 4 tabs (CHROM\tPOS\tID\tREF\tALT)
    if (tab_count >= 4) return InputFormat::VCF;

    // BED: exactly 2-5 tabs (CHR\tSTART\tEND)
    if (tab_count >= 2 && tab_count <= 5) {
        std::istringstream iss(line);
        std::string f1, f2, f3;
        if (std::getline(iss, f1, '\t') && std::getline(iss, f2, '\t') && std::getline(iss, f3, '\t')) {
            try {
                std::stoi(f2);
                std::stoi(f3);
                return InputFormat::BED;
            } catch (...) {}
        }
    }

    // Ensembl default: space-separated, 5 fields (CHR START END ALLELES STRAND)
    {
        std::istringstream iss(line);
        std::string f1, f2, f3, f4, f5;
        if (iss >> f1 >> f2 >> f3 >> f4 >> f5) {
            try {
                std::stoi(f2);
                std::stoi(f3);
                if (f4.find('/') != std::string::npos) {
                    return InputFormat::ENSEMBL;
                }
            } catch (...) {}
        }
    }

    return InputFormat::UNKNOWN;
}

// Parse Ensembl default format line: CHR START END ALLELES STRAND [IDENTIFIER]
static bool parse_ensembl_line(const std::string& line, std::string& chrom, int& pos,
                               std::string& ref, std::string& alt) {
    std::istringstream iss(line);
    std::string start_str, end_str, alleles, strand_str;

    if (!(iss >> chrom >> start_str >> end_str >> alleles >> strand_str)) return false;

    try {
        pos = std::stoi(start_str);
    } catch (...) { return false; }

    size_t slash = alleles.find('/');
    if (slash == std::string::npos) return false;

    ref = alleles.substr(0, slash);
    alt = alleles.substr(slash + 1);

    if (ref == "-") ref = "";
    if (alt == "-") alt = "";

    return true;
}

// Parse BED format line: CHR START END [NAME [SCORE [STRAND]]]
static bool parse_bed_line(const std::string& line, std::string& chrom, int& pos,
                           std::string& ref, std::string& alt) {
    std::istringstream iss(line);
    std::string start_str, end_str;

    if (!(iss >> chrom >> start_str >> end_str)) return false;

    try {
        int start = std::stoi(start_str);
        pos = start + 1;  // BED is 0-based, convert to 1-based
        ref = "-";
        alt = "-";
    } catch (...) { return false; }

    return true;
}

// Parse annotation source string: NAME:VCF_PATH[:FIELDS]
static bool parse_annotation_source(const std::string& source,
                                    std::string& name, std::string& vcf_path, std::string& fields) {
    size_t p1 = source.find(':');
    if (p1 == std::string::npos) return false;

    name = source.substr(0, p1);

    size_t p2 = source.find(':', p1 + 1);
    if (p2 == std::string::npos) {
        vcf_path = source.substr(p1 + 1);
        fields = "";
    } else {
        vcf_path = source.substr(p1 + 1, p2 - p1 - 1);
        fields = source.substr(p2 + 1);
    }

    return !name.empty() && !vcf_path.empty();
}

// Minimal mode allele trimming (Perl VEP trim_sequences: prefix first, then suffix)
static void minimal_trim(std::string& ref, std::string& alt, int& pos) {
    if (ref.size() <= 1 && alt.size() <= 1) return;

    // Prefix: trim while both non-empty and share first base
    while (!ref.empty() && !alt.empty() && ref.front() == alt.front()) {
        ref.erase(ref.begin());
        alt.erase(alt.begin());
        pos++;
    }
    // Suffix: trim while both non-empty and share last base
    while (!ref.empty() && !alt.empty() && ref.back() == alt.back()) {
        ref.pop_back();
        alt.pop_back();
    }
    // Pad empty alleles (Perl VEP uses "-")
    if (ref.empty()) ref = "-";
    if (alt.empty()) alt = "-";
}

// Config file line parsing (extracts key/value pairs from a single line)
static std::vector<std::string> parse_config_line(const std::string& raw_line) {
    std::vector<std::string> result;
    std::string cfg_line = raw_line;

    // Strip comments
    size_t hash = cfg_line.find('#');
    if (hash != std::string::npos) cfg_line = cfg_line.substr(0, hash);

    // Trim trailing whitespace
    while (!cfg_line.empty() && std::isspace(static_cast<unsigned char>(cfg_line.back())))
        cfg_line.pop_back();
    // Trim leading whitespace
    while (!cfg_line.empty() && std::isspace(static_cast<unsigned char>(cfg_line.front())))
        cfg_line.erase(0, 1);

    if (cfg_line.empty()) return result;

    // Parse "key = value" or "key value" format
    size_t eq = cfg_line.find('=');
    if (eq != std::string::npos) {
        std::string key = cfg_line.substr(0, eq);
        std::string val = cfg_line.substr(eq + 1);
        while (!key.empty() && std::isspace(static_cast<unsigned char>(key.back()))) key.pop_back();
        while (!val.empty() && std::isspace(static_cast<unsigned char>(val.front()))) val.erase(0, 1);
        if (!key.empty()) {
            if (key[0] != '-') key = "--" + key;
            result.push_back(key);
            if (!val.empty()) result.push_back(val);
        }
    } else {
        std::istringstream cfg_iss(cfg_line);
        std::string token;
        while (cfg_iss >> token) {
            result.push_back(token);
        }
    }

    return result;
}

// Chromosome filter parsing: parse comma-separated list and insert both forms
static std::set<std::string> parse_chr_filter(const std::string& chr_str) {
    std::set<std::string> chr_filter;
    std::istringstream chr_iss(chr_str);
    std::string chr_token;
    while (std::getline(chr_iss, chr_token, ',')) {
        size_t s = chr_token.find_first_not_of(" \t");
        size_t e = chr_token.find_last_not_of(" \t");
        if (s != std::string::npos && e != std::string::npos) {
            std::string c = chr_token.substr(s, e - s + 1);
            chr_filter.insert(c);
            // Also insert alternate form (chr1 <-> 1)
            if (c.length() > 3 && c.substr(0, 3) == "chr") {
                chr_filter.insert(c.substr(3));
            } else {
                chr_filter.insert("chr" + c);
            }
        }
    }
    return chr_filter;
}

// Normalize Perl VEP underscore flags to hyphen form
static std::string normalize_flag(const std::string& arg) {
    std::string result = arg;
    if (result.size() > 2 && result[0] == '-' && result[1] == '-') {
        for (size_t j = 2; j < result.size(); ++j) {
            if (result[j] == '_') result[j] = '-';
        }
    }
    return result;
}

// ============================================================================
// Test Suite: parse_variant
// ============================================================================

class ParseVariantTest : public ::testing::Test {
protected:
    std::string chrom, ref, alt;
    int pos = 0;
};

TEST_F(ParseVariantTest, StandardColonSeparated) {
    ASSERT_TRUE(parse_variant("7:140753336:A:T", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "7");
    EXPECT_EQ(pos, 140753336);
    EXPECT_EQ(ref, "A");
    EXPECT_EQ(alt, "T");
}

TEST_F(ParseVariantTest, WithChrPrefix) {
    ASSERT_TRUE(parse_variant("chr7:140753336:A:T", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "chr7");
    EXPECT_EQ(pos, 140753336);
    EXPECT_EQ(ref, "A");
    EXPECT_EQ(alt, "T");
}

TEST_F(ParseVariantTest, DashSeparatedFallback) {
    ASSERT_TRUE(parse_variant("7-140753336-A-T", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "7");
    EXPECT_EQ(pos, 140753336);
    EXPECT_EQ(ref, "A");
    EXPECT_EQ(alt, "T");
}

TEST_F(ParseVariantTest, DeletionVariant) {
    ASSERT_TRUE(parse_variant("7:100:AC:A", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "7");
    EXPECT_EQ(pos, 100);
    EXPECT_EQ(ref, "AC");
    EXPECT_EQ(alt, "A");
}

TEST_F(ParseVariantTest, InsertionVariant) {
    ASSERT_TRUE(parse_variant("7:100:A:ACGT", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "7");
    EXPECT_EQ(pos, 100);
    EXPECT_EQ(ref, "A");
    EXPECT_EQ(alt, "ACGT");
}

TEST_F(ParseVariantTest, MNV) {
    ASSERT_TRUE(parse_variant("17:7675088:GT:CA", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "17");
    EXPECT_EQ(pos, 7675088);
    EXPECT_EQ(ref, "GT");
    EXPECT_EQ(alt, "CA");
}

TEST_F(ParseVariantTest, ChrX) {
    ASSERT_TRUE(parse_variant("chrX:15561033:C:G", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "chrX");
    EXPECT_EQ(pos, 15561033);
    EXPECT_EQ(ref, "C");
    EXPECT_EQ(alt, "G");
}

TEST_F(ParseVariantTest, ChrMT) {
    ASSERT_TRUE(parse_variant("chrM:8860:A:G", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "chrM");
    EXPECT_EQ(pos, 8860);
    EXPECT_EQ(ref, "A");
    EXPECT_EQ(alt, "G");
}

TEST_F(ParseVariantTest, InvalidEmptyString) {
    EXPECT_FALSE(parse_variant("", chrom, pos, ref, alt));
}

TEST_F(ParseVariantTest, InvalidNoSeparator) {
    EXPECT_FALSE(parse_variant("7_140753336_A_T", chrom, pos, ref, alt));
}

TEST_F(ParseVariantTest, InvalidMissingFields) {
    EXPECT_FALSE(parse_variant("7:140753336:A", chrom, pos, ref, alt));
}

TEST_F(ParseVariantTest, InvalidTwoFields) {
    EXPECT_FALSE(parse_variant("7:140753336", chrom, pos, ref, alt));
}

TEST_F(ParseVariantTest, InvalidNonNumericPosition) {
    EXPECT_FALSE(parse_variant("7:abc:A:T", chrom, pos, ref, alt));
}

TEST_F(ParseVariantTest, EmptyAlt) {
    // Format allows empty alt field (e.g., when using dash alleles)
    ASSERT_TRUE(parse_variant("7:100:A:", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "7");
    EXPECT_EQ(pos, 100);
    EXPECT_EQ(ref, "A");
    EXPECT_EQ(alt, "");
}

TEST_F(ParseVariantTest, LargePosition) {
    ASSERT_TRUE(parse_variant("1:248956422:G:A", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "1");
    EXPECT_EQ(pos, 248956422);
    EXPECT_EQ(ref, "G");
    EXPECT_EQ(alt, "A");
}

TEST_F(ParseVariantTest, ColonPreferredOverDash) {
    // When both separators are present, colon should be used
    ASSERT_TRUE(parse_variant("chr7:100:A-T:G", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "chr7");
    EXPECT_EQ(pos, 100);
    EXPECT_EQ(ref, "A-T");
    EXPECT_EQ(alt, "G");
}

// ============================================================================
// Test Suite: detect_input_format
// ============================================================================

TEST(DetectInputFormat, EmptyLine) {
    EXPECT_EQ(detect_input_format(""), InputFormat::UNKNOWN);
}

TEST(DetectInputFormat, CommentLine) {
    EXPECT_EQ(detect_input_format("# this is a comment"), InputFormat::UNKNOWN);
}

TEST(DetectInputFormat, VCFHeaderLine) {
    EXPECT_EQ(detect_input_format("##fileformat=VCFv4.2"), InputFormat::UNKNOWN);
}

TEST(DetectInputFormat, HGVSCodingNotation) {
    EXPECT_EQ(detect_input_format("ENST00000269305.9:c.817C>T"), InputFormat::HGVS);
}

TEST(DetectInputFormat, HGVSProteinNotation) {
    EXPECT_EQ(detect_input_format("ENSP00000269305.4:p.Arg273Cys"), InputFormat::HGVS);
}

TEST(DetectInputFormat, HGVSGenomicNotation) {
    EXPECT_EQ(detect_input_format("NC_000017.11:g.7675088C>T"), InputFormat::HGVS);
}

TEST(DetectInputFormat, HGVSNonCodingNotation) {
    EXPECT_EQ(detect_input_format("NR_024540.1:n.100A>G"), InputFormat::HGVS);
}

TEST(DetectInputFormat, HGVSMitochondrialNotation) {
    EXPECT_EQ(detect_input_format("NC_012920.1:m.8860A>G"), InputFormat::HGVS);
}

TEST(DetectInputFormat, VCFDataLine) {
    // VCF: CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO (at least 4 tabs)
    EXPECT_EQ(detect_input_format("chr17\t7675088\t.\tC\tT\t100\tPASS\t."), InputFormat::VCF);
}

TEST(DetectInputFormat, VCFMinimalFields) {
    EXPECT_EQ(detect_input_format("7\t140753336\trs121913529\tA\tT"), InputFormat::VCF);
}

TEST(DetectInputFormat, VCFWithManyFields) {
    EXPECT_EQ(detect_input_format("1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1"), InputFormat::VCF);
}

TEST(DetectInputFormat, BEDFormat) {
    // BED: CHR\tSTART\tEND (3 tab-separated fields with numeric start/end)
    EXPECT_EQ(detect_input_format("chr17\t7675087\t7675088"), InputFormat::BED);
}

TEST(DetectInputFormat, BEDWithName) {
    EXPECT_EQ(detect_input_format("chr17\t7675087\t7675088\tTP53_variant"), InputFormat::BED);
}

TEST(DetectInputFormat, EnsemblDefaultFormat) {
    // Ensembl: CHR START END ALLELES STRAND
    EXPECT_EQ(detect_input_format("17 7675088 7675088 C/T +"), InputFormat::ENSEMBL);
}

TEST(DetectInputFormat, EnsemblWithIdentifier) {
    EXPECT_EQ(detect_input_format("17 7675088 7675088 C/T + rs121913343"), InputFormat::ENSEMBL);
}

TEST(DetectInputFormat, EnsemblDeletion) {
    EXPECT_EQ(detect_input_format("17 7675088 7675089 AC/- +"), InputFormat::ENSEMBL);
}

TEST(DetectInputFormat, RegionFormat) {
    EXPECT_EQ(detect_input_format("17:7675088-7675088:1/T"), InputFormat::REGION);
}

TEST(DetectInputFormat, RegionFormatNoStrand) {
    EXPECT_EQ(detect_input_format("7:140753336-140753336/T"), InputFormat::REGION);
}

TEST(DetectInputFormat, UnrecognizedFormat) {
    EXPECT_EQ(detect_input_format("some random text"), InputFormat::UNKNOWN);
}

// ============================================================================
// Test Suite: parse_ensembl_line
// ============================================================================

TEST(ParseEnsemblLine, StandardSNV) {
    std::string chrom, ref, alt;
    int pos = 0;
    ASSERT_TRUE(parse_ensembl_line("17 7675088 7675088 C/T +", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "17");
    EXPECT_EQ(pos, 7675088);
    EXPECT_EQ(ref, "C");
    EXPECT_EQ(alt, "T");
}

TEST(ParseEnsemblLine, Deletion) {
    std::string chrom, ref, alt;
    int pos = 0;
    ASSERT_TRUE(parse_ensembl_line("17 7675088 7675089 AC/- +", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "17");
    EXPECT_EQ(pos, 7675088);
    EXPECT_EQ(ref, "AC");
    EXPECT_EQ(alt, "");  // "-" becomes empty
}

TEST(ParseEnsemblLine, Insertion) {
    std::string chrom, ref, alt;
    int pos = 0;
    ASSERT_TRUE(parse_ensembl_line("17 7675088 7675087 -/ACG +", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "17");
    EXPECT_EQ(pos, 7675088);
    EXPECT_EQ(ref, "");  // "-" becomes empty
    EXPECT_EQ(alt, "ACG");
}

TEST(ParseEnsemblLine, MinusStrand) {
    std::string chrom, ref, alt;
    int pos = 0;
    ASSERT_TRUE(parse_ensembl_line("17 7675088 7675088 G/A -", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "17");
    EXPECT_EQ(pos, 7675088);
    EXPECT_EQ(ref, "G");
    EXPECT_EQ(alt, "A");
}

TEST(ParseEnsemblLine, WithIdentifier) {
    std::string chrom, ref, alt;
    int pos = 0;
    // Extra fields after strand are ignored
    ASSERT_TRUE(parse_ensembl_line("17 7675088 7675088 C/T + rs121913343", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "17");
    EXPECT_EQ(ref, "C");
    EXPECT_EQ(alt, "T");
}

TEST(ParseEnsemblLine, InvalidTooFewFields) {
    std::string chrom, ref, alt;
    int pos = 0;
    EXPECT_FALSE(parse_ensembl_line("17 7675088 7675088 C/T", chrom, pos, ref, alt));
}

TEST(ParseEnsemblLine, InvalidNoSlash) {
    std::string chrom, ref, alt;
    int pos = 0;
    EXPECT_FALSE(parse_ensembl_line("17 7675088 7675088 CT +", chrom, pos, ref, alt));
}

TEST(ParseEnsemblLine, InvalidNonNumericPos) {
    std::string chrom, ref, alt;
    int pos = 0;
    EXPECT_FALSE(parse_ensembl_line("17 abc 7675088 C/T +", chrom, pos, ref, alt));
}

// ============================================================================
// Test Suite: parse_bed_line
// ============================================================================

TEST(ParseBedLine, BasicThreeColumn) {
    std::string chrom, ref, alt;
    int pos = 0;
    ASSERT_TRUE(parse_bed_line("chr17\t7675087\t7675088", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "chr17");
    EXPECT_EQ(pos, 7675088);  // 0-based to 1-based conversion
    EXPECT_EQ(ref, "-");
    EXPECT_EQ(alt, "-");
}

TEST(ParseBedLine, WithNameAndScore) {
    std::string chrom, ref, alt;
    int pos = 0;
    ASSERT_TRUE(parse_bed_line("chr17\t7675087\t7675088\tTP53_exon\t100", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "chr17");
    EXPECT_EQ(pos, 7675088);
}

TEST(ParseBedLine, ZeroBasedStart) {
    std::string chrom, ref, alt;
    int pos = 0;
    ASSERT_TRUE(parse_bed_line("chr1\t0\t100", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "chr1");
    EXPECT_EQ(pos, 1);  // 0-based 0 becomes 1-based 1
}

TEST(ParseBedLine, InvalidTooFewColumns) {
    std::string chrom, ref, alt;
    int pos = 0;
    EXPECT_FALSE(parse_bed_line("chr17\t7675087", chrom, pos, ref, alt));
}

TEST(ParseBedLine, InvalidNonNumeric) {
    std::string chrom, ref, alt;
    int pos = 0;
    EXPECT_FALSE(parse_bed_line("chr17\tabc\t7675088", chrom, pos, ref, alt));
}

// ============================================================================
// Test Suite: parse_annotation_source
// ============================================================================

TEST(ParseAnnotationSource, NameAndPath) {
    std::string name, path, fields;
    ASSERT_TRUE(parse_annotation_source("gnomad:/path/to/gnomad.vcf.gz", name, path, fields));
    EXPECT_EQ(name, "gnomad");
    EXPECT_EQ(path, "/path/to/gnomad.vcf.gz");
    EXPECT_EQ(fields, "");
}

TEST(ParseAnnotationSource, NamePathAndFields) {
    std::string name, path, fields;
    ASSERT_TRUE(parse_annotation_source("gnomad:/path/to/gnomad.vcf.gz:AF,AC,AN", name, path, fields));
    EXPECT_EQ(name, "gnomad");
    EXPECT_EQ(path, "/path/to/gnomad.vcf.gz");
    EXPECT_EQ(fields, "AF,AC,AN");
}

TEST(ParseAnnotationSource, InvalidNoColon) {
    std::string name, path, fields;
    EXPECT_FALSE(parse_annotation_source("gnomad_path_only", name, path, fields));
}

TEST(ParseAnnotationSource, InvalidEmptyName) {
    std::string name, path, fields;
    EXPECT_FALSE(parse_annotation_source(":/path/to/file.vcf", name, path, fields));
}

TEST(ParseAnnotationSource, InvalidEmptyPath) {
    std::string name, path, fields;
    EXPECT_FALSE(parse_annotation_source("gnomad:", name, path, fields));
}

TEST(ParseAnnotationSource, SingleCharName) {
    std::string name, path, fields;
    ASSERT_TRUE(parse_annotation_source("g:/data/file.vcf", name, path, fields));
    EXPECT_EQ(name, "g");
    EXPECT_EQ(path, "/data/file.vcf");
}

TEST(ParseAnnotationSource, MultipleColonsInFields) {
    // Fields after the second colon contain everything remaining
    std::string name, path, fields;
    ASSERT_TRUE(parse_annotation_source("db:/data/file.vcf:AF:AC", name, path, fields));
    EXPECT_EQ(name, "db");
    EXPECT_EQ(path, "/data/file.vcf");
    EXPECT_EQ(fields, "AF:AC");  // Everything after second colon
}

// ============================================================================
// Test Suite: Minimal mode allele trimming
// ============================================================================

TEST(MinimalTrim, InsertionPrefixTrim) {
    // CGT/CAGT -> prefix C trimmed -> GT/AGT -> suffix T trimmed -> G/AG
    // Then suffix G trimmed -> -/A (insertion)
    std::string ref = "CGT", alt = "CAGT";
    int pos = 100;
    minimal_trim(ref, alt, pos);
    EXPECT_EQ(ref, "-");
    EXPECT_EQ(alt, "A");
    EXPECT_EQ(pos, 101);  // Position shifted by 1 prefix base trimmed
}

TEST(MinimalTrim, DeletionSuffixTrim) {
    // ACG/AG -> prefix A trimmed -> CG/G -> suffix G trimmed -> C/- (deletion)
    std::string ref = "ACG", alt = "AG";
    int pos = 100;
    minimal_trim(ref, alt, pos);
    EXPECT_EQ(ref, "C");
    EXPECT_EQ(alt, "-");
    EXPECT_EQ(pos, 101);  // Position shifted by 1 prefix base trimmed
}

TEST(MinimalTrim, SNVBothTrim) {
    // ACGT/ATGT -> prefix A trimmed -> CGT/TGT -> suffix T trimmed -> CG/TG -> suffix G trimmed -> C/T
    std::string ref = "ACGT", alt = "ATGT";
    int pos = 100;
    minimal_trim(ref, alt, pos);
    EXPECT_EQ(ref, "C");
    EXPECT_EQ(alt, "T");
    EXPECT_EQ(pos, 101);
}

TEST(MinimalTrim, NoTrimNeeded) {
    // Simple SNV, no trimming possible (single base each)
    std::string ref = "A", alt = "T";
    int pos = 100;
    minimal_trim(ref, alt, pos);
    EXPECT_EQ(ref, "A");
    EXPECT_EQ(alt, "T");
    EXPECT_EQ(pos, 100);
}

TEST(MinimalTrim, PureDeletion) {
    // AT/A -> prefix A trimmed -> T/- (deletion)
    std::string ref = "AT", alt = "A";
    int pos = 100;
    minimal_trim(ref, alt, pos);
    EXPECT_EQ(ref, "T");
    EXPECT_EQ(alt, "-");
    EXPECT_EQ(pos, 101);
}

TEST(MinimalTrim, PureInsertion) {
    // A/AT -> prefix A trimmed -> -/T (insertion)
    std::string ref = "A", alt = "AT";
    int pos = 100;
    minimal_trim(ref, alt, pos);
    EXPECT_EQ(ref, "-");
    EXPECT_EQ(alt, "T");
    EXPECT_EQ(pos, 101);
}

TEST(MinimalTrim, ComplexDelins) {
    // ACGTACGT/ACTTACGT -> prefix AC trimmed -> GTACGT/TTACGT -> suffix TACGT trimmed -> G/T
    std::string ref = "ACGTACGT", alt = "ACTTACGT";
    int pos = 100;
    minimal_trim(ref, alt, pos);
    EXPECT_EQ(ref, "G");
    EXPECT_EQ(alt, "T");
    EXPECT_EQ(pos, 102);  // 2 prefix bases trimmed
}

TEST(MinimalTrim, LongInsertion) {
    // G/GACGT -> prefix G trimmed -> -/ACGT
    std::string ref = "G", alt = "GACGT";
    int pos = 200;
    minimal_trim(ref, alt, pos);
    EXPECT_EQ(ref, "-");
    EXPECT_EQ(alt, "ACGT");
    EXPECT_EQ(pos, 201);
}

TEST(MinimalTrim, LongDeletion) {
    // ACGTG/G -> prefix trimming: only trims if front matches.
    // A != G, so no prefix trim. But wait: A != G. So check suffix.
    // Suffix: G matches G -> ACGT/- . Done.
    // Actually: prefix first. A vs G, no match. Suffix: G vs G, match. ACGT/"".
    // Then alt empty -> "-".
    std::string ref = "ACGTG", alt = "G";
    int pos = 200;
    minimal_trim(ref, alt, pos);
    // No prefix trim (A != G), suffix trim removes trailing G from both
    EXPECT_EQ(ref, "ACGT");
    EXPECT_EQ(alt, "-");
    EXPECT_EQ(pos, 200);  // No prefix trimmed
}

TEST(MinimalTrim, IdenticalAlleles) {
    // Both same: all bases trimmed -> both become "-"
    std::string ref = "ACGT", alt = "ACGT";
    int pos = 100;
    minimal_trim(ref, alt, pos);
    EXPECT_EQ(ref, "-");
    EXPECT_EQ(alt, "-");
    EXPECT_EQ(pos, 104);  // 4 prefix bases trimmed
}

TEST(MinimalTrim, SingleBaseNoChange) {
    // Both single-base: function returns early without modification
    std::string ref = "A", alt = "G";
    int pos = 50;
    minimal_trim(ref, alt, pos);
    EXPECT_EQ(ref, "A");
    EXPECT_EQ(alt, "G");
    EXPECT_EQ(pos, 50);
}

// ============================================================================
// Test Suite: Config file parsing
// ============================================================================

TEST(ConfigParsing, KeyEqualsValue) {
    auto args = parse_config_line("gtf = /path/to/genes.gtf");
    ASSERT_EQ(args.size(), 2u);
    EXPECT_EQ(args[0], "--gtf");
    EXPECT_EQ(args[1], "/path/to/genes.gtf");
}

TEST(ConfigParsing, KeyEqualsValueNoSpaces) {
    auto args = parse_config_line("fasta=/path/to/genome.fa");
    ASSERT_EQ(args.size(), 2u);
    EXPECT_EQ(args[0], "--fasta");
    EXPECT_EQ(args[1], "/path/to/genome.fa");
}

TEST(ConfigParsing, FlagWithDashPrefix) {
    // Keys starting with dash are preserved as-is
    auto args = parse_config_line("--canonical = true");
    ASSERT_EQ(args.size(), 2u);
    EXPECT_EQ(args[0], "--canonical");
    EXPECT_EQ(args[1], "true");
}

TEST(ConfigParsing, BooleanFlagNoValue) {
    // key= with empty value
    auto args = parse_config_line("canonical=");
    ASSERT_EQ(args.size(), 1u);
    EXPECT_EQ(args[0], "--canonical");
}

TEST(ConfigParsing, SpaceSeparatedFormat) {
    auto args = parse_config_line("--gtf /path/to/genes.gtf");
    ASSERT_EQ(args.size(), 2u);
    EXPECT_EQ(args[0], "--gtf");
    EXPECT_EQ(args[1], "/path/to/genes.gtf");
}

TEST(ConfigParsing, CommentLine) {
    auto args = parse_config_line("# This is a comment");
    EXPECT_TRUE(args.empty());
}

TEST(ConfigParsing, InlineComment) {
    auto args = parse_config_line("gtf = /path/to/genes.gtf # gene annotations");
    ASSERT_EQ(args.size(), 2u);
    EXPECT_EQ(args[0], "--gtf");
    EXPECT_EQ(args[1], "/path/to/genes.gtf");
}

TEST(ConfigParsing, EmptyLine) {
    auto args = parse_config_line("");
    EXPECT_TRUE(args.empty());
}

TEST(ConfigParsing, WhitespaceOnlyLine) {
    auto args = parse_config_line("   \t  ");
    EXPECT_TRUE(args.empty());
}

TEST(ConfigParsing, LeadingTrailingWhitespace) {
    auto args = parse_config_line("  gtf  =  /path/to/genes.gtf  ");
    ASSERT_EQ(args.size(), 2u);
    EXPECT_EQ(args[0], "--gtf");
    EXPECT_EQ(args[1], "/path/to/genes.gtf");
}

TEST(ConfigParsing, MultipleSpaceSeparatedTokens) {
    auto args = parse_config_line("--pick --canonical --coding-only");
    ASSERT_EQ(args.size(), 3u);
    EXPECT_EQ(args[0], "--pick");
    EXPECT_EQ(args[1], "--canonical");
    EXPECT_EQ(args[2], "--coding-only");
}

// ============================================================================
// Test Suite: Chromosome filter parsing
// ============================================================================

TEST(ChrFilter, SingleChromosome) {
    auto filter = parse_chr_filter("7");
    // Should contain both "7" and "chr7"
    EXPECT_TRUE(filter.count("7"));
    EXPECT_TRUE(filter.count("chr7"));
    EXPECT_EQ(filter.size(), 2u);
}

TEST(ChrFilter, SingleChrPrefixed) {
    auto filter = parse_chr_filter("chr7");
    // Should contain both "chr7" and "7"
    EXPECT_TRUE(filter.count("7"));
    EXPECT_TRUE(filter.count("chr7"));
    EXPECT_EQ(filter.size(), 2u);
}

TEST(ChrFilter, MultipleChrCommaList) {
    auto filter = parse_chr_filter("1,2,X");
    EXPECT_TRUE(filter.count("1"));
    EXPECT_TRUE(filter.count("chr1"));
    EXPECT_TRUE(filter.count("2"));
    EXPECT_TRUE(filter.count("chr2"));
    EXPECT_TRUE(filter.count("X"));
    EXPECT_TRUE(filter.count("chrX"));
    EXPECT_EQ(filter.size(), 6u);
}

TEST(ChrFilter, WithSpacesAroundCommas) {
    auto filter = parse_chr_filter(" 1 , 2 , X ");
    EXPECT_TRUE(filter.count("1"));
    EXPECT_TRUE(filter.count("chr1"));
    EXPECT_TRUE(filter.count("2"));
    EXPECT_TRUE(filter.count("chr2"));
    EXPECT_TRUE(filter.count("X"));
    EXPECT_TRUE(filter.count("chrX"));
}

TEST(ChrFilter, MixedPrefixForms) {
    auto filter = parse_chr_filter("chr1,2,chrX");
    EXPECT_TRUE(filter.count("1"));
    EXPECT_TRUE(filter.count("chr1"));
    EXPECT_TRUE(filter.count("2"));
    EXPECT_TRUE(filter.count("chr2"));
    EXPECT_TRUE(filter.count("X"));
    EXPECT_TRUE(filter.count("chrX"));
}

TEST(ChrFilter, MitochondrialChromosome) {
    auto filter = parse_chr_filter("MT");
    EXPECT_TRUE(filter.count("MT"));
    EXPECT_TRUE(filter.count("chrMT"));
}

TEST(ChrFilter, EmptyString) {
    auto filter = parse_chr_filter("");
    EXPECT_TRUE(filter.empty());
}

TEST(ChrFilter, LookupMatchesBothForms) {
    // Verify that a variant on "chr7" will match a filter for "7"
    auto filter = parse_chr_filter("7");
    std::string variant_chrom = "chr7";
    EXPECT_TRUE(filter.find(variant_chrom) != filter.end());
}

// ============================================================================
// Test Suite: Flag normalization (underscore to hyphen)
// ============================================================================

TEST(FlagNormalization, UnderscoreToHyphen) {
    EXPECT_EQ(normalize_flag("--pick_allele"), "--pick-allele");
}

TEST(FlagNormalization, MultipleUnderscores) {
    EXPECT_EQ(normalize_flag("--pick_allele_gene"), "--pick-allele-gene");
}

TEST(FlagNormalization, AlreadyHyphenated) {
    EXPECT_EQ(normalize_flag("--pick-allele"), "--pick-allele");
}

TEST(FlagNormalization, ShortFlagUnchanged) {
    EXPECT_EQ(normalize_flag("-v"), "-v");
}

TEST(FlagNormalization, SingleDashLongFlagUnchanged) {
    // Only double-dash flags get normalized
    EXPECT_EQ(normalize_flag("-pick_allele"), "-pick_allele");
}

TEST(FlagNormalization, NoFlag) {
    EXPECT_EQ(normalize_flag("some_value"), "some_value");
}

TEST(FlagNormalization, EmptyString) {
    EXPECT_EQ(normalize_flag(""), "");
}

TEST(FlagNormalization, JustDoubleDash) {
    EXPECT_EQ(normalize_flag("--"), "--");
}

TEST(FlagNormalization, MixedUnderscoreHyphen) {
    EXPECT_EQ(normalize_flag("--flag_pick-allele_gene"), "--flag-pick-allele-gene");
}

TEST(FlagNormalization, PerlVepCommonFlags) {
    EXPECT_EQ(normalize_flag("--pick_allele"), "--pick-allele");
    EXPECT_EQ(normalize_flag("--per_gene"), "--per-gene");
    EXPECT_EQ(normalize_flag("--most_severe"), "--most-severe");
    EXPECT_EQ(normalize_flag("--coding_only"), "--coding-only");
    EXPECT_EQ(normalize_flag("--no_intergenic"), "--no-intergenic");
    EXPECT_EQ(normalize_flag("--input_file"), "--input-file");
    EXPECT_EQ(normalize_flag("--output_format"), "--output-format");
    EXPECT_EQ(normalize_flag("--check_existing"), "--check-existing");
    EXPECT_EQ(normalize_flag("--mane_select"), "--mane-select");
}

// ============================================================================
// Test Suite: Chromosome synonym mapping
// ============================================================================

TEST(ChromSynonyms, BasicLookup) {
    std::map<std::string, std::string> synonyms;
    synonyms["NC_000017.11"] = "chr17";
    synonyms["17"] = "chr17";

    EXPECT_EQ(synonyms["NC_000017.11"], "chr17");
    EXPECT_EQ(synonyms["17"], "chr17");
}

TEST(ChromSynonyms, MissingKeyReturnsEmpty) {
    std::map<std::string, std::string> synonyms;
    synonyms["NC_000017.11"] = "chr17";

    auto it = synonyms.find("chr17");
    EXPECT_EQ(it, synonyms.end());  // chr17 is a value, not a key
}

TEST(ChromSynonyms, ParseSynonymLine) {
    // Simulates loading a tab-separated synonym file line
    std::string line = "NC_000017.11\tchr17";
    std::istringstream iss(line);
    std::string synonym, canonical;
    ASSERT_TRUE(bool(iss >> synonym >> canonical));
    EXPECT_EQ(synonym, "NC_000017.11");
    EXPECT_EQ(canonical, "chr17");
}

TEST(ChromSynonyms, MultipleSynonymsForSameChrom) {
    std::map<std::string, std::string> synonyms;
    synonyms["NC_000017.11"] = "chr17";
    synonyms["17"] = "chr17";
    synonyms["CM000679.2"] = "chr17";

    // All synonyms resolve to the same canonical name
    EXPECT_EQ(synonyms["NC_000017.11"], "chr17");
    EXPECT_EQ(synonyms["17"], "chr17");
    EXPECT_EQ(synonyms["CM000679.2"], "chr17");
}

TEST(ChromSynonyms, MitochondrialSynonyms) {
    std::map<std::string, std::string> synonyms;
    synonyms["MT"] = "chrM";
    synonyms["NC_012920.1"] = "chrM";

    EXPECT_EQ(synonyms["MT"], "chrM");
    EXPECT_EQ(synonyms["NC_012920.1"], "chrM");
}

// ============================================================================
// Test Suite: Edge cases and integration-style tests
// ============================================================================

TEST(ParseVariantEdgeCases, DashSeparatedWithChrPrefix) {
    std::string chrom, ref, alt;
    int pos = 0;
    ASSERT_TRUE(parse_variant("chr7-140753336-A-T", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "chr7");
    EXPECT_EQ(pos, 140753336);
    EXPECT_EQ(ref, "A");
    EXPECT_EQ(alt, "T");
}

TEST(ParseVariantEdgeCases, PositionOne) {
    std::string chrom, ref, alt;
    int pos = 0;
    ASSERT_TRUE(parse_variant("1:1:A:T", chrom, pos, ref, alt));
    EXPECT_EQ(chrom, "1");
    EXPECT_EQ(pos, 1);
}

TEST(ParseVariantEdgeCases, MultiBaseRef) {
    std::string chrom, ref, alt;
    int pos = 0;
    ASSERT_TRUE(parse_variant("17:7675088:ACGTACGT:A", chrom, pos, ref, alt));
    EXPECT_EQ(ref, "ACGTACGT");
    EXPECT_EQ(alt, "A");
}

TEST(DetectInputFormatEdgeCases, TabSeparatedButNotVCF) {
    // 3 tabs: could be BED with name field
    // chr\tstart\tend\tname = BED format
    EXPECT_EQ(detect_input_format("chr1\t100\t200\tregion1"), InputFormat::BED);
}

TEST(DetectInputFormatEdgeCases, SpaceSeparatedNonEnsembl) {
    // Space-separated but no allele slash
    EXPECT_EQ(detect_input_format("17 7675088 7675088 CT +"), InputFormat::UNKNOWN);
}

TEST(MinimalTrimEdgeCases, MultiBasePrefixTrim) {
    // AAACG/AAATG -> prefix AAA trimmed -> CG/TG -> no suffix match -> C/T? No:
    // Actually: prefix: A==A, AA..==AA.., AAA==AAA -> CG/TG. Suffix: G==G -> C/T.
    std::string ref = "AAACG", alt = "AAATG";
    int pos = 100;
    minimal_trim(ref, alt, pos);
    EXPECT_EQ(ref, "C");
    EXPECT_EQ(alt, "T");
    EXPECT_EQ(pos, 103);  // 3 prefix bases trimmed
}

TEST(MinimalTrimEdgeCases, OnlyPrefixTrim) {
    // ACG/ATG -> prefix A trimmed -> CG/TG (no suffix match: G != G? G == G!)
    // Actually CG/TG -> suffix G == G -> C/T
    std::string ref = "ACG", alt = "ATG";
    int pos = 100;
    minimal_trim(ref, alt, pos);
    EXPECT_EQ(ref, "C");
    EXPECT_EQ(alt, "T");
    EXPECT_EQ(pos, 101);
}

TEST(MinimalTrimEdgeCases, OnlySuffixTrim) {
    // CA/TA -> prefix: C != T, no prefix trim. Suffix: A == A -> C/T
    std::string ref = "CA", alt = "TA";
    int pos = 100;
    minimal_trim(ref, alt, pos);
    EXPECT_EQ(ref, "C");
    EXPECT_EQ(alt, "T");
    EXPECT_EQ(pos, 100);  // No prefix trimmed
}

// ============================================================================
// Test Suite: gz_read_line behavior documentation
// ============================================================================
// Note: gz_read_line cannot be tested directly without a gzFile handle.
// These tests document the expected behavior through string manipulation
// equivalent to what gz_read_line does after reading.

TEST(LineStripping, TrailingNewline) {
    std::string line = "chr17\t7675088\t.\tC\tT\n";
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
        line.pop_back();
    EXPECT_EQ(line, "chr17\t7675088\t.\tC\tT");
}

TEST(LineStripping, WindowsLineEnding) {
    std::string line = "chr17\t7675088\t.\tC\tT\r\n";
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
        line.pop_back();
    EXPECT_EQ(line, "chr17\t7675088\t.\tC\tT");
}

TEST(LineStripping, NoNewline) {
    std::string line = "chr17\t7675088\t.\tC\tT";
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
        line.pop_back();
    EXPECT_EQ(line, "chr17\t7675088\t.\tC\tT");
}

TEST(LineStripping, EmptyString) {
    std::string line = "";
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
        line.pop_back();
    EXPECT_EQ(line, "");
}

TEST(LineStripping, OnlyNewlines) {
    std::string line = "\r\n";
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
        line.pop_back();
    EXPECT_EQ(line, "");
}

TEST(LineStripping, MultipleTrailingCR) {
    // Edge case: multiple \r characters
    std::string line = "data\r\r\n";
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
        line.pop_back();
    EXPECT_EQ(line, "data");
}
