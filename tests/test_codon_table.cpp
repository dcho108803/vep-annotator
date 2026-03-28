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

// ============================================================================
// Case insensitivity (3 tests)
// ============================================================================

TEST(CodonTable, CaseInsensitivityAllLower) {
    EXPECT_EQ(CodonTable::translate("atg"), 'M');
}

TEST(CodonTable, CaseInsensitivityMixedCase1) {
    EXPECT_EQ(CodonTable::translate("Atg"), 'M');
}

TEST(CodonTable, CaseInsensitivityMixedCase2) {
    EXPECT_EQ(CodonTable::translate("aTg"), 'M');
}

// ============================================================================
// All stop codons (2 tests)
// ============================================================================

TEST(CodonTable, StopCodonLowercase) {
    EXPECT_TRUE(CodonTable::is_stop_codon("taa"));
    EXPECT_TRUE(CodonTable::is_stop_codon("tag"));
    EXPECT_TRUE(CodonTable::is_stop_codon("tga"));
}

TEST(CodonTable, StopCodonChromosomeAware) {
    // In standard table, TGA is a stop codon
    EXPECT_TRUE(CodonTable::is_stop_codon("TGA", "1"));
    // In mitochondrial table, TGA is NOT a stop codon (codes for Trp)
    EXPECT_FALSE(CodonTable::is_stop_codon("TGA", "MT"));
}

// ============================================================================
// Mitochondrial stop codons (3 tests)
// ============================================================================

TEST(CodonTable, MitochondrialStopAGA) {
    EXPECT_TRUE(CodonTable::is_stop_codon("AGA", "MT"));
}

TEST(CodonTable, MitochondrialStopAGG) {
    EXPECT_TRUE(CodonTable::is_stop_codon("AGG", "MT"));
}

TEST(CodonTable, MitochondrialStopTAAStillStop) {
    // TAA is still a stop codon in mitochondrial code
    EXPECT_TRUE(CodonTable::is_stop_codon("TAA", "MT"));
    // TGA is NOT a stop in MT (it encodes Trp)
    EXPECT_FALSE(CodonTable::is_stop_codon("TGA", "MT"));
}

// ============================================================================
// Three-letter edge cases (3 tests)
// ============================================================================

TEST(CodonTable, ThreeLetterUnknownX) {
    EXPECT_EQ(CodonTable::get_three_letter('X'), "Xaa");
}

TEST(CodonTable, ThreeLetterSelenocysteine) {
    EXPECT_EQ(CodonTable::get_three_letter('U'), "Sec");
}

TEST(CodonTable, ThreeLetterInvalidChar) {
    // '?' is not in the table, should return default "Xaa"
    EXPECT_EQ(CodonTable::get_three_letter('?'), "Xaa");
}

// ============================================================================
// All remaining codons (5 tests)
// ============================================================================

TEST(CodonTable, RemainingIsoleucineAndThreonine) {
    // Isoleucine: ATC, ATA
    EXPECT_EQ(CodonTable::translate("ATC"), 'I');
    EXPECT_EQ(CodonTable::translate("ATA"), 'I');

    // Threonine: ACC, ACA, ACG
    EXPECT_EQ(CodonTable::translate("ACC"), 'T');
    EXPECT_EQ(CodonTable::translate("ACA"), 'T');
    EXPECT_EQ(CodonTable::translate("ACG"), 'T');
}

TEST(CodonTable, RemainingValineHistidineGlutamine) {
    // Valine: GTC, GTA, GTG
    EXPECT_EQ(CodonTable::translate("GTC"), 'V');
    EXPECT_EQ(CodonTable::translate("GTA"), 'V');
    EXPECT_EQ(CodonTable::translate("GTG"), 'V');

    // Histidine: CAC
    EXPECT_EQ(CodonTable::translate("CAC"), 'H');

    // Glutamine: CAG
    EXPECT_EQ(CodonTable::translate("CAG"), 'Q');
}

TEST(CodonTable, RemainingAsparagineLysineAspartateGlutamate) {
    // Asparagine: AAC
    EXPECT_EQ(CodonTable::translate("AAC"), 'N');

    // Lysine: AAG
    EXPECT_EQ(CodonTable::translate("AAG"), 'K');

    // Aspartate: GAC
    EXPECT_EQ(CodonTable::translate("GAC"), 'D');

    // Glutamate: GAG
    EXPECT_EQ(CodonTable::translate("GAG"), 'E');
}

TEST(CodonTable, RemainingPheTyrTrpCys) {
    // Phenylalanine: TTC
    EXPECT_EQ(CodonTable::translate("TTC"), 'F');

    // Tyrosine: TAC
    EXPECT_EQ(CodonTable::translate("TAC"), 'Y');

    // Tryptophan: TGG (only codon)
    EXPECT_EQ(CodonTable::translate("TGG"), 'W');

    // Cysteine: TGC
    EXPECT_EQ(CodonTable::translate("TGC"), 'C');
}

TEST(CodonTable, RemainingProlineGlycine) {
    // Proline: CCC, CCA, CCG
    EXPECT_EQ(CodonTable::translate("CCC"), 'P');
    EXPECT_EQ(CodonTable::translate("CCA"), 'P');
    EXPECT_EQ(CodonTable::translate("CCG"), 'P');

    // Glycine: GGC, GGA, GGG
    EXPECT_EQ(CodonTable::translate("GGC"), 'G');
    EXPECT_EQ(CodonTable::translate("GGA"), 'G');
    EXPECT_EQ(CodonTable::translate("GGG"), 'G');
}

// ============================================================================
// Codon table completeness (3 tests)
// ============================================================================

TEST(CodonTableCompleteness, All64CodonsReturnValidAminoAcid) {
    const std::string bases = "ACGT";
    int valid_count = 0;
    for (char b1 : bases) {
        for (char b2 : bases) {
            for (char b3 : bases) {
                std::string codon = {b1, b2, b3};
                char aa = CodonTable::translate(codon);
                EXPECT_NE(aa, 'X')
                    << "Codon " << codon << " returned 'X' (unknown)";
                valid_count++;
            }
        }
    }
    EXPECT_EQ(valid_count, 64);
}

TEST(CodonTableCompleteness, ExactlyThreeStopCodonsInStandardTable) {
    const std::string bases = "ACGT";
    int stop_count = 0;
    for (char b1 : bases) {
        for (char b2 : bases) {
            for (char b3 : bases) {
                std::string codon = {b1, b2, b3};
                if (CodonTable::translate(codon) == '*') {
                    stop_count++;
                }
            }
        }
    }
    EXPECT_EQ(stop_count, 3);
}

TEST(CodonTableCompleteness, ExactlyOneStartCodon) {
    const std::string bases = "ACGT";
    int start_count = 0;
    for (char b1 : bases) {
        for (char b2 : bases) {
            for (char b3 : bases) {
                std::string codon = {b1, b2, b3};
                if (CodonTable::is_start_codon(codon)) {
                    start_count++;
                    EXPECT_EQ(codon, "ATG");
                }
            }
        }
    }
    EXPECT_EQ(start_count, 1);
}

// ============================================================================
// Mitochondrial completeness (3 tests)
// ============================================================================

TEST(MitochondrialCompleteness, AllFourOverridesProduceExpectedResults) {
    EXPECT_EQ(CodonTable::translate_mt("AGA"), '*');  // Arg -> Stop
    EXPECT_EQ(CodonTable::translate_mt("AGG"), '*');  // Arg -> Stop
    EXPECT_EQ(CodonTable::translate_mt("ATA"), 'M');  // Ile -> Met
    EXPECT_EQ(CodonTable::translate_mt("TGA"), 'W');  // Stop -> Trp
}

TEST(MitochondrialCompleteness, NonOverriddenCodonsMatchStandardTable) {
    const std::string bases = "ACGT";
    // Check all 64 codons: non-overridden ones should match standard table
    std::set<std::string> overridden = {"AGA", "AGG", "ATA", "TGA"};
    for (char b1 : bases) {
        for (char b2 : bases) {
            for (char b3 : bases) {
                std::string codon = {b1, b2, b3};
                if (overridden.find(codon) == overridden.end()) {
                    EXPECT_EQ(CodonTable::translate_mt(codon), CodonTable::translate(codon))
                        << "MT translation for " << codon << " differs from standard table unexpectedly";
                }
            }
        }
    }
}

TEST(MitochondrialCompleteness, MTHasFourStopCodons) {
    const std::string bases = "ACGT";
    int mt_stop_count = 0;
    std::vector<std::string> mt_stops;
    for (char b1 : bases) {
        for (char b2 : bases) {
            for (char b3 : bases) {
                std::string codon = {b1, b2, b3};
                if (CodonTable::translate_mt(codon) == '*') {
                    mt_stop_count++;
                    mt_stops.push_back(codon);
                }
            }
        }
    }
    // MT stop codons: TAA, TAG, AGA, AGG (TGA is Trp in MT)
    EXPECT_EQ(mt_stop_count, 4);
    EXPECT_NE(std::find(mt_stops.begin(), mt_stops.end(), "TAA"), mt_stops.end());
    EXPECT_NE(std::find(mt_stops.begin(), mt_stops.end(), "TAG"), mt_stops.end());
    EXPECT_NE(std::find(mt_stops.begin(), mt_stops.end(), "AGA"), mt_stops.end());
    EXPECT_NE(std::find(mt_stops.begin(), mt_stops.end(), "AGG"), mt_stops.end());
}

// ============================================================================
// Invalid inputs (3 tests)
// ============================================================================

TEST(CodonTable, CodonWithNBases) {
    // N (ambiguous) bases should return 'X'
    EXPECT_EQ(CodonTable::translate("ANG"), 'X');
    EXPECT_EQ(CodonTable::translate("NAT"), 'X');
    EXPECT_EQ(CodonTable::translate("NNN"), 'X');
}

TEST(CodonTable, CodonTooLong) {
    // 4-character string is not a valid codon
    EXPECT_EQ(CodonTable::translate("ATGC"), 'X');
}

TEST(CodonTable, CodonTooShort) {
    // Single character is not a valid codon
    EXPECT_EQ(CodonTable::translate("A"), 'X');
}
