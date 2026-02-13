/**
 * Tests for HGVS parsing, generation, SPDI parsing, RefSeq mapping,
 * and Ensembl format parsing
 */

#include <gtest/gtest.h>
#include "hgvs_parser.hpp"

using namespace vep;

// ============================================================================
// HGVS notation detection
// ============================================================================

TEST(HGVSDetection, IsHGVS) {
    EXPECT_TRUE(is_hgvs_notation("NC_000007.14:g.140753336A>T"));
    EXPECT_TRUE(is_hgvs_notation("ENST00000366667:c.803C>T"));
    EXPECT_TRUE(is_hgvs_notation("BRAF:p.Val600Glu"));
    EXPECT_TRUE(is_hgvs_notation("NR_024540.1:n.1234A>G"));
    EXPECT_TRUE(is_hgvs_notation("NC_012920.1:m.1234A>G"));
}

TEST(HGVSDetection, NotHGVS) {
    EXPECT_FALSE(is_hgvs_notation("7:140753336:A:T"));
    EXPECT_FALSE(is_hgvs_notation("chr7:140753336"));
    EXPECT_FALSE(is_hgvs_notation("BRAF"));
    EXPECT_FALSE(is_hgvs_notation(""));
    EXPECT_FALSE(is_hgvs_notation(":c.803C>T"));
}

// ============================================================================
// Genomic HGVS parsing
// ============================================================================

TEST(HGVSParsing, GenomicSubstitution) {
    auto result = parse_hgvs("NC_000007.14:g.140753336A>T");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.hgvs_type, HGVSType::GENOMIC);
    EXPECT_EQ(result.variant_type, HGVSVariantType::SUBSTITUTION);
    EXPECT_EQ(result.reference_id, "NC_000007.14");
    EXPECT_EQ(result.start_pos, 140753336);
    EXPECT_EQ(result.ref_allele, "A");
    EXPECT_EQ(result.alt_allele, "T");
}

TEST(HGVSParsing, GenomicDeletion) {
    auto result = parse_hgvs("NC_000007.14:g.140753336_140753340del");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.variant_type, HGVSVariantType::DELETION);
    EXPECT_EQ(result.start_pos, 140753336);
    EXPECT_EQ(result.end_pos, 140753340);
}

TEST(HGVSParsing, GenomicInsertion) {
    auto result = parse_hgvs("NC_000007.14:g.140753336_140753337insACGT");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.variant_type, HGVSVariantType::INSERTION);
    EXPECT_EQ(result.start_pos, 140753336);
    EXPECT_EQ(result.end_pos, 140753337);
    EXPECT_EQ(result.alt_allele, "ACGT");
}

TEST(HGVSParsing, GenomicDuplication) {
    auto result = parse_hgvs("NC_000007.14:g.140753336_140753340dup");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.variant_type, HGVSVariantType::DUPLICATION);
    EXPECT_EQ(result.start_pos, 140753336);
    EXPECT_EQ(result.end_pos, 140753340);
}

TEST(HGVSParsing, GenomicDelins) {
    auto result = parse_hgvs("NC_000007.14:g.140753336_140753340delinsACGT");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.variant_type, HGVSVariantType::DELINS);
    EXPECT_EQ(result.start_pos, 140753336);
    EXPECT_EQ(result.end_pos, 140753340);
    EXPECT_EQ(result.alt_allele, "ACGT");
}

// ============================================================================
// Coding HGVS parsing
// ============================================================================

TEST(HGVSParsing, CodingSubstitution) {
    auto result = parse_hgvs("ENST00000366667:c.803C>T");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.hgvs_type, HGVSType::CODING);
    EXPECT_EQ(result.variant_type, HGVSVariantType::SUBSTITUTION);
    EXPECT_EQ(result.reference_id, "ENST00000366667");
    EXPECT_EQ(result.start_pos, 803);
    EXPECT_EQ(result.ref_allele, "C");
    EXPECT_EQ(result.alt_allele, "T");
}

TEST(HGVSParsing, CodingIntronicOffset) {
    auto result = parse_hgvs("NM_000546.5:c.123+5G>A");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.start_pos, 123);
    EXPECT_EQ(result.intron_offset, 5);
    EXPECT_EQ(result.ref_allele, "G");
    EXPECT_EQ(result.alt_allele, "A");
}

TEST(HGVSParsing, CodingNegativeIntronicOffset) {
    auto result = parse_hgvs("NM_000546.5:c.124-10G>A");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.start_pos, 124);
    EXPECT_EQ(result.intron_offset, -10);
}

TEST(HGVSParsing, CodingDeletion) {
    auto result = parse_hgvs("ENST00000366667:c.803del");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.variant_type, HGVSVariantType::DELETION);
    EXPECT_EQ(result.start_pos, 803);
}

TEST(HGVSParsing, CodingInsertion) {
    auto result = parse_hgvs("ENST00000366667:c.803_804insACG");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.variant_type, HGVSVariantType::INSERTION);
    EXPECT_EQ(result.start_pos, 803);
    EXPECT_EQ(result.end_pos, 804);
    EXPECT_EQ(result.alt_allele, "ACG");
}

// ============================================================================
// Protein HGVS parsing
// ============================================================================

TEST(HGVSParsing, ProteinMissense3Letter) {
    auto result = parse_hgvs("BRAF:p.Val600Glu");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.hgvs_type, HGVSType::PROTEIN);
    EXPECT_EQ(result.variant_type, HGVSVariantType::SUBSTITUTION);
    EXPECT_EQ(result.ref_aa, "Val");
    EXPECT_EQ(result.protein_pos, 600);
    EXPECT_EQ(result.alt_aa, "Glu");
}

TEST(HGVSParsing, ProteinMissense1Letter) {
    auto result = parse_hgvs("BRAF:p.V600E");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.protein_pos, 600);
    EXPECT_EQ(result.ref_aa, "Val");
    EXPECT_EQ(result.alt_aa, "Glu");
}

TEST(HGVSParsing, ProteinNonsense) {
    auto result = parse_hgvs("TP53:p.Arg234Ter");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.ref_aa, "Arg");
    EXPECT_EQ(result.protein_pos, 234);
    EXPECT_EQ(result.alt_aa, "Ter");
}

TEST(HGVSParsing, ProteinSynonymous) {
    auto result = parse_hgvs("BRAF:p.Val600=");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.ref_aa, "Val");
    EXPECT_EQ(result.alt_aa, "Val");
    EXPECT_EQ(result.protein_pos, 600);
}

TEST(HGVSParsing, ProteinFrameshift) {
    auto result = parse_hgvs("BRAF:p.Val600fs");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.ref_aa, "Val");
    EXPECT_EQ(result.protein_pos, 600);
    EXPECT_EQ(result.alt_aa, "fs");
}

// ============================================================================
// HGVS type parsing
// ============================================================================

TEST(HGVSType, ParseTypes) {
    EXPECT_EQ(parse_hgvs_type("g"), HGVSType::GENOMIC);
    EXPECT_EQ(parse_hgvs_type("c"), HGVSType::CODING);
    EXPECT_EQ(parse_hgvs_type("n"), HGVSType::NONCODING);
    EXPECT_EQ(parse_hgvs_type("p"), HGVSType::PROTEIN);
    EXPECT_EQ(parse_hgvs_type("m"), HGVSType::MITOCHONDRIAL);
    EXPECT_EQ(parse_hgvs_type("r"), HGVSType::RNA);
    EXPECT_EQ(parse_hgvs_type("x"), HGVSType::UNKNOWN);
}

TEST(HGVSType, TypeToString) {
    EXPECT_EQ(hgvs_type_to_string(HGVSType::GENOMIC), "g");
    EXPECT_EQ(hgvs_type_to_string(HGVSType::CODING), "c");
    EXPECT_EQ(hgvs_type_to_string(HGVSType::PROTEIN), "p");
}

// ============================================================================
// Amino acid conversion
// ============================================================================

TEST(AminoAcid, ThreeToOne) {
    EXPECT_EQ(aa_three_to_one("Ala"), "A");
    EXPECT_EQ(aa_three_to_one("Val"), "V");
    EXPECT_EQ(aa_three_to_one("Glu"), "E");
    EXPECT_EQ(aa_three_to_one("Ter"), "*");
    EXPECT_EQ(aa_three_to_one("ALA"), "A");  // Uppercase also works
}

TEST(AminoAcid, OneToThree) {
    EXPECT_EQ(aa_one_to_three('A'), "Ala");
    EXPECT_EQ(aa_one_to_three('V'), "Val");
    EXPECT_EQ(aa_one_to_three('E'), "Glu");
    EXPECT_EQ(aa_one_to_three('*'), "Ter");
}

// ============================================================================
// RefSeq to chromosome mapping
// ============================================================================

TEST(RefSeqMapping, GRCh38) {
    EXPECT_EQ(refseq_to_chromosome("NC_000007.14"), "7");
    EXPECT_EQ(refseq_to_chromosome("NC_000001.11"), "1");
    EXPECT_EQ(refseq_to_chromosome("NC_000023.11"), "X");
    EXPECT_EQ(refseq_to_chromosome("NC_000024.10"), "Y");
    EXPECT_EQ(refseq_to_chromosome("NC_012920.1"), "MT");
}

TEST(RefSeqMapping, GRCh37) {
    EXPECT_EQ(refseq_to_chromosome("NC_000007.13"), "7");
    EXPECT_EQ(refseq_to_chromosome("NC_000001.10"), "1");
}

TEST(RefSeqMapping, ChromosomeNames) {
    // Already chromosome names should pass through
    EXPECT_EQ(refseq_to_chromosome("1"), "1");
    EXPECT_EQ(refseq_to_chromosome("X"), "X");
    EXPECT_EQ(refseq_to_chromosome("MT"), "MT");
    EXPECT_EQ(refseq_to_chromosome("chr7"), "chr7");
}

TEST(RefSeqMapping, UnknownAccession) {
    EXPECT_EQ(refseq_to_chromosome("NC_999999.1"), "");
}

// ============================================================================
// HGVSg generation
// ============================================================================

TEST(HGVSGeneration, SubstitutionHGVSg) {
    std::string result = generate_hgvsg("7", 140753336, "A", "T");
    EXPECT_EQ(result, "NC_000007.14:g.140753336A>T");
}

TEST(HGVSGeneration, DeletionHGVSg) {
    std::string result = generate_hgvsg("7", 140753336, "ACG", "A");
    EXPECT_TRUE(result.find("del") != std::string::npos);
}

TEST(HGVSGeneration, InsertionHGVSg) {
    std::string result = generate_hgvsg("7", 140753336, "A", "ACGT");
    EXPECT_TRUE(result.find("ins") != std::string::npos);
}

// ============================================================================
// SPDI generation and parsing
// ============================================================================

TEST(SPDI, GenerateSPDI) {
    std::string result = generate_spdi("7", 140753336, "A", "T");
    // SPDI uses 0-based coordinates
    EXPECT_EQ(result, "NC_000007.14:140753335:A:T");
}

TEST(SPDI, ParseSPDI) {
    auto result = parse_spdi("NC_000007.14:140753335:A:T");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.chromosome, "7");
    EXPECT_EQ(result.position, 140753336);  // Converted to 1-based
    EXPECT_EQ(result.ref_allele, "A");
    EXPECT_EQ(result.alt_allele, "T");
}

TEST(SPDI, ParseSPDIWithChrom) {
    auto result = parse_spdi("NC_000001.11:100:A:T");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.chromosome, "1");
    EXPECT_EQ(result.position, 101);
}

TEST(SPDI, InvalidSPDI) {
    auto result = parse_spdi("invalid");
    EXPECT_FALSE(result.valid);

    auto result2 = parse_spdi("NC_000007.14:abc:A:T");
    EXPECT_FALSE(result2.valid);
}

TEST(SPDI, Detection) {
    EXPECT_TRUE(is_spdi_notation("NC_000007.14:140753335:A:T"));
    EXPECT_TRUE(is_spdi_notation("NC_000001.11:100:A:T"));
    EXPECT_FALSE(is_spdi_notation("7:140753336:A:T"));  // Regular variant format
    EXPECT_FALSE(is_spdi_notation("ENST00000366667:c.803C>T"));  // HGVS
}

// ============================================================================
// Ensembl format parsing
// ============================================================================

TEST(EnsemblFormat, BasicSNV) {
    auto result = parse_ensembl_format("1 100 100 A/T +");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.chromosome, "1");
    EXPECT_EQ(result.position, 100);
    EXPECT_EQ(result.ref_allele, "A");
    EXPECT_EQ(result.alt_allele, "T");
}

TEST(EnsemblFormat, Deletion) {
    auto result = parse_ensembl_format("X 100 102 ATG/- +");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.chromosome, "X");
    EXPECT_EQ(result.position, 100);
    EXPECT_EQ(result.ref_allele, "ATG");
    EXPECT_EQ(result.alt_allele, "");
}

TEST(EnsemblFormat, Insertion) {
    auto result = parse_ensembl_format("7 100 100 -/ACGT +");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.ref_allele, "");
    EXPECT_EQ(result.alt_allele, "ACGT");
}

TEST(EnsemblFormat, MinusStrand) {
    auto result = parse_ensembl_format("1 100 100 A/T -");
    EXPECT_TRUE(result.valid);
    // A complement = T, T complement = A
    EXPECT_EQ(result.ref_allele, "T");
    EXPECT_EQ(result.alt_allele, "A");
}

TEST(EnsemblFormat, Detection) {
    EXPECT_TRUE(is_ensembl_format("1 100 100 A/T +"));
    EXPECT_TRUE(is_ensembl_format("7 140753336 140753336 A/T 1"));
    EXPECT_FALSE(is_ensembl_format("7:140753336:A:T"));
    EXPECT_FALSE(is_ensembl_format("ENST00000366667:c.803C>T"));
    EXPECT_FALSE(is_ensembl_format(""));
}

TEST(EnsemblFormat, NoStrand) {
    auto result = parse_ensembl_format("1 100 100 A/T");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.ref_allele, "A");
    EXPECT_EQ(result.alt_allele, "T");
}

// ============================================================================
// Error handling
// ============================================================================

TEST(HGVSParsing, InvalidFormat) {
    auto result = parse_hgvs("invalid");
    EXPECT_FALSE(result.valid);
}

TEST(HGVSParsing, MissingColon) {
    auto result = parse_hgvs("BRAF.p.Val600Glu");
    EXPECT_FALSE(result.valid);
}

TEST(HGVSParsing, InvalidType) {
    auto result = parse_hgvs("BRAF:x.Val600Glu");
    EXPECT_FALSE(result.valid);
}

// ============================================================================
// REST-style region format
// ============================================================================

TEST(RESTRegion, Detection) {
    EXPECT_TRUE(is_rest_region_format("7:140753336-140753336:1/T"));
    EXPECT_TRUE(is_rest_region_format("1:12345-12346:-1/A"));
    EXPECT_FALSE(is_rest_region_format("7:140753336:A:T"));
    EXPECT_FALSE(is_rest_region_format("simple_string"));
}

TEST(RESTRegion, Parsing) {
    auto result = parse_rest_region("7:140753336-140753336:1/T");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.chromosome, "7");
    EXPECT_EQ(result.position, 140753336);
    EXPECT_EQ(result.end_position, 140753336);
    EXPECT_EQ(result.strand, 1);
    EXPECT_EQ(result.alt_allele, "T");
}

TEST(RESTRegion, NegativeStrand) {
    auto result = parse_rest_region("1:12345-12346:-1/A");
    EXPECT_TRUE(result.valid);
    EXPECT_EQ(result.chromosome, "1");
    EXPECT_EQ(result.position, 12345);
    EXPECT_EQ(result.end_position, 12346);
    EXPECT_EQ(result.strand, -1);
    EXPECT_EQ(result.alt_allele, "A");
}
