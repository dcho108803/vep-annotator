/**
 * Tests for codon translation table: standard codons, MT overrides, start/stop
 */

#include <gtest/gtest.h>
#include "vep_annotator.hpp"

using namespace vep;

// ============================================================================
// Standard codon table
// ============================================================================

TEST(CodonTable, StandardStartCodon) {
    EXPECT_EQ(CodonTable::translate("ATG"), 'M');
    EXPECT_TRUE(CodonTable::is_start_codon("ATG"));
}

TEST(CodonTable, StandardStopCodons) {
    EXPECT_EQ(CodonTable::translate("TAA"), '*');
    EXPECT_EQ(CodonTable::translate("TAG"), '*');
    EXPECT_EQ(CodonTable::translate("TGA"), '*');
    EXPECT_TRUE(CodonTable::is_stop_codon("TAA"));
    EXPECT_TRUE(CodonTable::is_stop_codon("TAG"));
    EXPECT_TRUE(CodonTable::is_stop_codon("TGA"));
}

TEST(CodonTable, NonStopCodons) {
    EXPECT_FALSE(CodonTable::is_stop_codon("ATG"));
    EXPECT_FALSE(CodonTable::is_stop_codon("GCA"));
    EXPECT_FALSE(CodonTable::is_stop_codon("TTT"));
}

TEST(CodonTable, AllAminoAcids) {
    // Test one codon for each amino acid
    EXPECT_EQ(CodonTable::translate("GCT"), 'A');  // Ala
    EXPECT_EQ(CodonTable::translate("TGT"), 'C');  // Cys
    EXPECT_EQ(CodonTable::translate("GAT"), 'D');  // Asp
    EXPECT_EQ(CodonTable::translate("GAA"), 'E');  // Glu
    EXPECT_EQ(CodonTable::translate("TTT"), 'F');  // Phe
    EXPECT_EQ(CodonTable::translate("GGT"), 'G');  // Gly
    EXPECT_EQ(CodonTable::translate("CAT"), 'H');  // His
    EXPECT_EQ(CodonTable::translate("ATT"), 'I');  // Ile
    EXPECT_EQ(CodonTable::translate("AAA"), 'K');  // Lys
    EXPECT_EQ(CodonTable::translate("TTA"), 'L');  // Leu
    EXPECT_EQ(CodonTable::translate("ATG"), 'M');  // Met
    EXPECT_EQ(CodonTable::translate("AAT"), 'N');  // Asn
    EXPECT_EQ(CodonTable::translate("CCT"), 'P');  // Pro
    EXPECT_EQ(CodonTable::translate("CAA"), 'Q');  // Gln
    EXPECT_EQ(CodonTable::translate("CGT"), 'R');  // Arg
    EXPECT_EQ(CodonTable::translate("TCT"), 'S');  // Ser
    EXPECT_EQ(CodonTable::translate("ACT"), 'T');  // Thr
    EXPECT_EQ(CodonTable::translate("GTT"), 'V');  // Val
    EXPECT_EQ(CodonTable::translate("TGG"), 'W');  // Trp
    EXPECT_EQ(CodonTable::translate("TAT"), 'Y');  // Tyr
}

TEST(CodonTable, DegenerateCodons) {
    // Alanine: GCN (4-fold degenerate)
    EXPECT_EQ(CodonTable::translate("GCA"), 'A');
    EXPECT_EQ(CodonTable::translate("GCC"), 'A');
    EXPECT_EQ(CodonTable::translate("GCG"), 'A');
    EXPECT_EQ(CodonTable::translate("GCT"), 'A');

    // Leucine: 6 codons
    EXPECT_EQ(CodonTable::translate("TTA"), 'L');
    EXPECT_EQ(CodonTable::translate("TTG"), 'L');
    EXPECT_EQ(CodonTable::translate("CTT"), 'L');
    EXPECT_EQ(CodonTable::translate("CTC"), 'L');
    EXPECT_EQ(CodonTable::translate("CTA"), 'L');
    EXPECT_EQ(CodonTable::translate("CTG"), 'L');

    // Serine: 6 codons
    EXPECT_EQ(CodonTable::translate("TCT"), 'S');
    EXPECT_EQ(CodonTable::translate("TCC"), 'S');
    EXPECT_EQ(CodonTable::translate("TCA"), 'S');
    EXPECT_EQ(CodonTable::translate("TCG"), 'S');
    EXPECT_EQ(CodonTable::translate("AGT"), 'S');
    EXPECT_EQ(CodonTable::translate("AGC"), 'S');

    // Arginine: 6 codons
    EXPECT_EQ(CodonTable::translate("CGT"), 'R');
    EXPECT_EQ(CodonTable::translate("CGC"), 'R');
    EXPECT_EQ(CodonTable::translate("CGA"), 'R');
    EXPECT_EQ(CodonTable::translate("CGG"), 'R');
    EXPECT_EQ(CodonTable::translate("AGA"), 'R');
    EXPECT_EQ(CodonTable::translate("AGG"), 'R');
}

TEST(CodonTable, UnknownCodon) {
    EXPECT_EQ(CodonTable::translate("NNN"), 'X');
    EXPECT_EQ(CodonTable::translate("XXX"), 'X');
    EXPECT_EQ(CodonTable::translate(""), 'X');
    EXPECT_EQ(CodonTable::translate("AT"), 'X');   // Too short
}

// ============================================================================
// Three-letter codes
// ============================================================================

TEST(CodonTable, ThreeLetterCodes) {
    EXPECT_EQ(CodonTable::get_three_letter('A'), "Ala");
    EXPECT_EQ(CodonTable::get_three_letter('C'), "Cys");
    EXPECT_EQ(CodonTable::get_three_letter('D'), "Asp");
    EXPECT_EQ(CodonTable::get_three_letter('E'), "Glu");
    EXPECT_EQ(CodonTable::get_three_letter('F'), "Phe");
    EXPECT_EQ(CodonTable::get_three_letter('G'), "Gly");
    EXPECT_EQ(CodonTable::get_three_letter('H'), "His");
    EXPECT_EQ(CodonTable::get_three_letter('I'), "Ile");
    EXPECT_EQ(CodonTable::get_three_letter('K'), "Lys");
    EXPECT_EQ(CodonTable::get_three_letter('L'), "Leu");
    EXPECT_EQ(CodonTable::get_three_letter('M'), "Met");
    EXPECT_EQ(CodonTable::get_three_letter('N'), "Asn");
    EXPECT_EQ(CodonTable::get_three_letter('P'), "Pro");
    EXPECT_EQ(CodonTable::get_three_letter('Q'), "Gln");
    EXPECT_EQ(CodonTable::get_three_letter('R'), "Arg");
    EXPECT_EQ(CodonTable::get_three_letter('S'), "Ser");
    EXPECT_EQ(CodonTable::get_three_letter('T'), "Thr");
    EXPECT_EQ(CodonTable::get_three_letter('V'), "Val");
    EXPECT_EQ(CodonTable::get_three_letter('W'), "Trp");
    EXPECT_EQ(CodonTable::get_three_letter('Y'), "Tyr");
    EXPECT_EQ(CodonTable::get_three_letter('*'), "Ter");
}

// ============================================================================
// Mitochondrial codon table
// ============================================================================

TEST(CodonTable, MitochondrialOverrides) {
    // AGA -> Stop (instead of Arg in standard table)
    EXPECT_EQ(CodonTable::translate_mt("AGA"), '*');
    // AGG -> Stop (instead of Arg)
    EXPECT_EQ(CodonTable::translate_mt("AGG"), '*');
    // ATA -> Met (instead of Ile)
    EXPECT_EQ(CodonTable::translate_mt("ATA"), 'M');
    // TGA -> Trp (instead of Stop)
    EXPECT_EQ(CodonTable::translate_mt("TGA"), 'W');
}

TEST(CodonTable, MitochondrialNonOverridden) {
    // Standard codons should still work
    EXPECT_EQ(CodonTable::translate_mt("ATG"), 'M');
    EXPECT_EQ(CodonTable::translate_mt("TAA"), '*');
    EXPECT_EQ(CodonTable::translate_mt("TAG"), '*');
    EXPECT_EQ(CodonTable::translate_mt("GCT"), 'A');
    EXPECT_EQ(CodonTable::translate_mt("TTT"), 'F');
}

// ============================================================================
// Chromosome-aware translation
// ============================================================================

TEST(CodonTable, ChromosomeAwareTranslation) {
    // Standard chromosome uses standard table
    EXPECT_EQ(CodonTable::translate("AGA", "7"), 'R');
    EXPECT_EQ(CodonTable::translate("TGA", "1"), '*');

    // MT chromosome uses mitochondrial table
    EXPECT_EQ(CodonTable::translate("AGA", "MT"), '*');
    EXPECT_EQ(CodonTable::translate("TGA", "MT"), 'W');
    EXPECT_EQ(CodonTable::translate("ATA", "MT"), 'M');

    // Alternative MT chromosome names
    EXPECT_EQ(CodonTable::translate("AGA", "M"), '*');
    EXPECT_EQ(CodonTable::translate("AGA", "chrM"), '*');
    EXPECT_EQ(CodonTable::translate("TGA", "chrMT"), 'W');
}
