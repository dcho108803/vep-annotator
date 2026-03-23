/**
 * Tests for structural variant support: SV type parsing, consequence determination
 */

#include <gtest/gtest.h>
#include "structural_variant.hpp"

using namespace vep;

// ============================================================================
// SV Type string conversion
// ============================================================================

TEST(SVTypeString, AllTypes) {
    EXPECT_EQ(sv_type_to_string(SVType::SNV), "SNV");
    EXPECT_EQ(sv_type_to_string(SVType::INS), "INS");
    EXPECT_EQ(sv_type_to_string(SVType::DEL), "DEL");
    EXPECT_EQ(sv_type_to_string(SVType::DUP), "DUP");
    EXPECT_EQ(sv_type_to_string(SVType::TDUP), "TDUP");
    EXPECT_EQ(sv_type_to_string(SVType::INV), "INV");
    EXPECT_EQ(sv_type_to_string(SVType::CNV), "CNV");
    EXPECT_EQ(sv_type_to_string(SVType::BND), "BND");
    EXPECT_EQ(sv_type_to_string(SVType::UNKNOWN), "UNKNOWN");
}

// ============================================================================
// SV Type parsing
// ============================================================================

TEST(SVTypeParsing, PlainTypes) {
    EXPECT_EQ(parse_sv_type("DEL"), SVType::DEL);
    EXPECT_EQ(parse_sv_type("DUP"), SVType::DUP);
    EXPECT_EQ(parse_sv_type("INS"), SVType::INS);
    EXPECT_EQ(parse_sv_type("INV"), SVType::INV);
    EXPECT_EQ(parse_sv_type("CNV"), SVType::CNV);
    EXPECT_EQ(parse_sv_type("BND"), SVType::BND);
}

TEST(SVTypeParsing, CaseInsensitive) {
    EXPECT_EQ(parse_sv_type("del"), SVType::DEL);
    EXPECT_EQ(parse_sv_type("Del"), SVType::DEL);
    EXPECT_EQ(parse_sv_type("dup"), SVType::DUP);
    EXPECT_EQ(parse_sv_type("ins"), SVType::INS);
}

TEST(SVTypeParsing, SymbolicAlleles) {
    EXPECT_EQ(parse_sv_type("<DEL>"), SVType::DEL);
    EXPECT_EQ(parse_sv_type("<DUP>"), SVType::DUP);
    EXPECT_EQ(parse_sv_type("<INS>"), SVType::INS);
    EXPECT_EQ(parse_sv_type("<INV>"), SVType::INV);
    EXPECT_EQ(parse_sv_type("<CNV>"), SVType::CNV);
}

TEST(SVTypeParsing, TandemDuplication) {
    EXPECT_EQ(parse_sv_type("TDUP"), SVType::TDUP);
    EXPECT_EQ(parse_sv_type("DUP:TANDEM"), SVType::TDUP);
}

TEST(SVTypeParsing, CopyNumber) {
    EXPECT_EQ(parse_sv_type("CN0"), SVType::CNV);
    EXPECT_EQ(parse_sv_type("CN1"), SVType::CNV);
    EXPECT_EQ(parse_sv_type("CN4"), SVType::CNV);
    EXPECT_EQ(parse_sv_type("<CN0>"), SVType::CNV);
}

TEST(SVTypeParsing, Unknown) {
    EXPECT_EQ(parse_sv_type(""), SVType::UNKNOWN);
    EXPECT_EQ(parse_sv_type("FOO"), SVType::UNKNOWN);
    EXPECT_EQ(parse_sv_type("COMPLEX"), SVType::UNKNOWN);
}

TEST(SVTypeParsing, Roundtrip) {
    // parse -> string -> parse roundtrip for all named types
    EXPECT_EQ(parse_sv_type(sv_type_to_string(SVType::DEL)), SVType::DEL);
    EXPECT_EQ(parse_sv_type(sv_type_to_string(SVType::DUP)), SVType::DUP);
    EXPECT_EQ(parse_sv_type(sv_type_to_string(SVType::INS)), SVType::INS);
    EXPECT_EQ(parse_sv_type(sv_type_to_string(SVType::INV)), SVType::INV);
    EXPECT_EQ(parse_sv_type(sv_type_to_string(SVType::CNV)), SVType::CNV);
    EXPECT_EQ(parse_sv_type(sv_type_to_string(SVType::BND)), SVType::BND);
}

// ============================================================================
// StructuralVariant properties
// ============================================================================

TEST(StructuralVariantProps, Length) {
    StructuralVariant sv;
    sv.start = 100;
    sv.end = 200;
    sv.sv_len = 0;
    EXPECT_EQ(sv.length(), 101);  // end - start + 1

    sv.sv_len = -50;
    EXPECT_EQ(sv.length(), 50);   // abs(sv_len) takes precedence

    sv.sv_len = 300;
    EXPECT_EQ(sv.length(), 300);
}

TEST(StructuralVariantProps, IsSV) {
    StructuralVariant sv;
    sv.sv_type = SVType::SNV;
    EXPECT_FALSE(sv.is_sv());

    sv.sv_type = SVType::UNKNOWN;
    EXPECT_FALSE(sv.is_sv());

    sv.sv_type = SVType::DEL;
    EXPECT_TRUE(sv.is_sv());

    sv.sv_type = SVType::BND;
    EXPECT_TRUE(sv.is_sv());
}

TEST(StructuralVariantProps, Overlaps) {
    StructuralVariant sv;
    sv.start = 100;
    sv.end = 200;

    EXPECT_TRUE(sv.overlaps(150, 250));   // partial overlap
    EXPECT_TRUE(sv.overlaps(50, 150));    // partial overlap
    EXPECT_TRUE(sv.overlaps(120, 180));   // contained
    EXPECT_TRUE(sv.overlaps(50, 300));    // containing
    EXPECT_TRUE(sv.overlaps(100, 100));   // touches start
    EXPECT_TRUE(sv.overlaps(200, 200));   // touches end
    EXPECT_FALSE(sv.overlaps(201, 300));  // no overlap
    EXPECT_FALSE(sv.overlaps(50, 99));    // no overlap
}

TEST(StructuralVariantProps, Contains) {
    StructuralVariant sv;
    sv.start = 100;
    sv.end = 200;

    EXPECT_TRUE(sv.contains(100, 200));   // exact match
    EXPECT_TRUE(sv.contains(120, 180));   // strictly contained
    EXPECT_TRUE(sv.contains(150, 150));   // single position
    EXPECT_FALSE(sv.contains(50, 200));   // extends before
    EXPECT_FALSE(sv.contains(100, 300));  // extends after
    EXPECT_FALSE(sv.contains(50, 300));   // extends both ways
}

// ============================================================================
// SV Consequences
// ============================================================================

static Transcript make_coding_transcript() {
    Transcript t;
    t.id = "ENST00000001";
    t.chromosome = "chr1";
    t.start = 1000;
    t.end = 5000;
    t.strand = '+';
    t.cds_start = 1200;
    t.cds_end = 4800;
    t.biotype = "protein_coding";
    return t;
}

TEST(SVConsequences, NoOverlap) {
    Transcript t = make_coding_transcript();
    StructuralVariant sv;
    sv.sv_type = SVType::DEL;
    sv.start = 100;
    sv.end = 500;

    auto cons = get_sv_consequences(sv, t);
    ASSERT_FALSE(cons.empty());
    // Should be upstream for + strand (SV before transcript)
    EXPECT_EQ(cons[0], ConsequenceType::UPSTREAM_GENE_VARIANT);
}

TEST(SVConsequences, NoOverlapDownstream) {
    Transcript t = make_coding_transcript();
    StructuralVariant sv;
    sv.sv_type = SVType::DEL;
    sv.start = 5500;
    sv.end = 6000;

    auto cons = get_sv_consequences(sv, t);
    ASSERT_FALSE(cons.empty());
    EXPECT_EQ(cons[0], ConsequenceType::DOWNSTREAM_GENE_VARIANT);
}

TEST(SVConsequences, FarAway) {
    Transcript t = make_coding_transcript();
    StructuralVariant sv;
    sv.sv_type = SVType::DEL;
    sv.start = 50000;
    sv.end = 60000;

    auto cons = get_sv_consequences(sv, t);
    ASSERT_FALSE(cons.empty());
    EXPECT_EQ(cons[0], ConsequenceType::INTERGENIC_VARIANT);
}

TEST(SVConsequences, TranscriptAblation) {
    Transcript t = make_coding_transcript();
    StructuralVariant sv;
    sv.sv_type = SVType::DEL;
    sv.start = 900;
    sv.end = 5100;  // contains entire transcript

    auto cons = get_sv_consequences(sv, t);
    ASSERT_EQ(cons.size(), 1u);
    EXPECT_EQ(cons[0], ConsequenceType::TRANSCRIPT_ABLATION);
}

TEST(SVConsequences, TranscriptAmplification) {
    Transcript t = make_coding_transcript();
    StructuralVariant sv;
    sv.sv_type = SVType::DUP;
    sv.start = 900;
    sv.end = 5100;

    auto cons = get_sv_consequences(sv, t);
    ASSERT_EQ(cons.size(), 1u);
    EXPECT_EQ(cons[0], ConsequenceType::TRANSCRIPT_AMPLIFICATION);
}

TEST(SVConsequences, MinusStrandUpstream) {
    Transcript t = make_coding_transcript();
    t.strand = '-';

    StructuralVariant sv;
    sv.sv_type = SVType::DEL;
    sv.start = 100;
    sv.end = 500;

    auto cons = get_sv_consequences(sv, t);
    ASSERT_FALSE(cons.empty());
    // For minus strand, SV before transcript genomically is downstream
    EXPECT_EQ(cons[0], ConsequenceType::DOWNSTREAM_GENE_VARIANT);
}
