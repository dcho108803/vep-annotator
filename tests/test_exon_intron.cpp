/**
 * Tests for exon/intron numbering: position lookup, strand handling, formatting
 */

#include <gtest/gtest.h>
#include "exon_intron_numbers.hpp"

using namespace vep;

// ============================================================================
// Helper: make a 3-exon transcript layout
// Exons: [100-200], [400-500], [700-800]
// Introns: [201-399], [501-699]
// ============================================================================

static const std::vector<int> EXON_STARTS = {100, 400, 700};
static const std::vector<int> EXON_ENDS   = {200, 500, 800};

// ============================================================================
// Plus strand exon numbering
// ============================================================================

TEST(ExonIntronNumber, PlusStrandExon1) {
    auto info = get_exon_intron_number(150, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_TRUE(info.found);
    EXPECT_TRUE(info.is_exon);
    EXPECT_EQ(info.number, 1);
    EXPECT_EQ(info.total_exons, 3);
    EXPECT_EQ(info.total_introns, 2);
}

TEST(ExonIntronNumber, PlusStrandExon2) {
    auto info = get_exon_intron_number(450, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_TRUE(info.found);
    EXPECT_TRUE(info.is_exon);
    EXPECT_EQ(info.number, 2);
}

TEST(ExonIntronNumber, PlusStrandExon3) {
    auto info = get_exon_intron_number(750, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_TRUE(info.found);
    EXPECT_TRUE(info.is_exon);
    EXPECT_EQ(info.number, 3);
}

TEST(ExonIntronNumber, PlusStrandExonBoundaryStart) {
    auto info = get_exon_intron_number(100, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_TRUE(info.found);
    EXPECT_TRUE(info.is_exon);
    EXPECT_EQ(info.number, 1);
    EXPECT_EQ(info.position_in_feature, 1);
}

TEST(ExonIntronNumber, PlusStrandExonBoundaryEnd) {
    auto info = get_exon_intron_number(200, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_TRUE(info.found);
    EXPECT_TRUE(info.is_exon);
    EXPECT_EQ(info.number, 1);
    EXPECT_EQ(info.feature_length, 101);  // 200 - 100 + 1
}

// ============================================================================
// Plus strand intron numbering
// ============================================================================

TEST(ExonIntronNumber, PlusStrandIntron1) {
    auto info = get_exon_intron_number(300, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_TRUE(info.found);
    EXPECT_FALSE(info.is_exon);
    EXPECT_EQ(info.number, 1);
    EXPECT_EQ(info.total_introns, 2);
}

TEST(ExonIntronNumber, PlusStrandIntron2) {
    auto info = get_exon_intron_number(600, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_TRUE(info.found);
    EXPECT_FALSE(info.is_exon);
    EXPECT_EQ(info.number, 2);
}

TEST(ExonIntronNumber, IntronBoundaryStart) {
    auto info = get_exon_intron_number(201, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_TRUE(info.found);
    EXPECT_FALSE(info.is_exon);
    EXPECT_EQ(info.number, 1);
    EXPECT_EQ(info.position_in_feature, 1);
}

TEST(ExonIntronNumber, IntronBoundaryEnd) {
    auto info = get_exon_intron_number(399, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_TRUE(info.found);
    EXPECT_FALSE(info.is_exon);
    EXPECT_EQ(info.number, 1);
    EXPECT_EQ(info.feature_length, 199);  // 399 - 201 + 1
}

// ============================================================================
// Minus strand numbering (reversed)
// ============================================================================

TEST(ExonIntronNumber, MinusStrandExon1IsLast) {
    // First exon genomically (100-200) should be exon 3 on minus strand
    auto info = get_exon_intron_number(150, EXON_STARTS, EXON_ENDS, '-');
    EXPECT_TRUE(info.found);
    EXPECT_TRUE(info.is_exon);
    EXPECT_EQ(info.number, 3);
}

TEST(ExonIntronNumber, MinusStrandExon3IsFirst) {
    // Last exon genomically (700-800) should be exon 1 on minus strand
    auto info = get_exon_intron_number(750, EXON_STARTS, EXON_ENDS, '-');
    EXPECT_TRUE(info.found);
    EXPECT_TRUE(info.is_exon);
    EXPECT_EQ(info.number, 1);
}

TEST(ExonIntronNumber, MinusStrandIntron1) {
    // Second intron genomically (501-699) should be intron 1 on minus strand
    auto info = get_exon_intron_number(600, EXON_STARTS, EXON_ENDS, '-');
    EXPECT_TRUE(info.found);
    EXPECT_FALSE(info.is_exon);
    EXPECT_EQ(info.number, 1);
}

TEST(ExonIntronNumber, MinusStrandIntron2) {
    // First intron genomically (201-399) should be intron 2 on minus strand
    auto info = get_exon_intron_number(300, EXON_STARTS, EXON_ENDS, '-');
    EXPECT_TRUE(info.found);
    EXPECT_FALSE(info.is_exon);
    EXPECT_EQ(info.number, 2);
}

// ============================================================================
// Position outside transcript
// ============================================================================

TEST(ExonIntronNumber, BeforeTranscript) {
    auto info = get_exon_intron_number(50, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_FALSE(info.found);
}

TEST(ExonIntronNumber, AfterTranscript) {
    auto info = get_exon_intron_number(900, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_FALSE(info.found);
}

// ============================================================================
// Single exon transcript
// ============================================================================

TEST(ExonIntronNumber, SingleExon) {
    std::vector<int> starts = {100};
    std::vector<int> ends = {500};

    auto info = get_exon_intron_number(300, starts, ends, '+');
    EXPECT_TRUE(info.found);
    EXPECT_TRUE(info.is_exon);
    EXPECT_EQ(info.number, 1);
    EXPECT_EQ(info.total_exons, 1);
    EXPECT_EQ(info.total_introns, 0);
}

TEST(ExonIntronNumber, SingleExonOutside) {
    std::vector<int> starts = {100};
    std::vector<int> ends = {500};

    auto info = get_exon_intron_number(600, starts, ends, '+');
    EXPECT_FALSE(info.found);
}

// ============================================================================
// Empty input
// ============================================================================

TEST(ExonIntronNumber, EmptyExons) {
    std::vector<int> starts;
    std::vector<int> ends;

    auto info = get_exon_intron_number(100, starts, ends, '+');
    EXPECT_FALSE(info.found);
}

TEST(ExonIntronNumber, MismatchedSizes) {
    std::vector<int> starts = {100, 400};
    std::vector<int> ends = {200};

    auto info = get_exon_intron_number(150, starts, ends, '+');
    EXPECT_FALSE(info.found);
}

// ============================================================================
// Formatting
// ============================================================================

TEST(ExonIntronFormat, ExonToString) {
    auto info = get_exon_intron_number(450, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_EQ(info.to_string(), "2/3");
    EXPECT_EQ(format_exon_number(info), "2/3");
    EXPECT_EQ(format_intron_number(info), "");  // is_exon, so no intron string
}

TEST(ExonIntronFormat, IntronToString) {
    auto info = get_exon_intron_number(300, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_EQ(info.to_string(), "1/2");
    EXPECT_EQ(format_intron_number(info), "1/2");
    EXPECT_EQ(format_exon_number(info), "");  // is intron, so no exon string
}

TEST(ExonIntronFormat, NotFoundToString) {
    auto info = get_exon_intron_number(50, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_EQ(info.to_string(), "");
    EXPECT_EQ(format_exon_number(info), "");
    EXPECT_EQ(format_intron_number(info), "");
}

TEST(ExonIntronFormat, FeatureType) {
    auto exon_info = get_exon_intron_number(150, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_EQ(exon_info.feature_type(), "exon");

    auto intron_info = get_exon_intron_number(300, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_EQ(intron_info.feature_type(), "intron");

    ExonIntronInfo empty;
    EXPECT_EQ(empty.feature_type(), "");
}

// ============================================================================
// Position within feature
// ============================================================================

TEST(ExonIntronNumber, PositionInExon) {
    auto info = get_exon_intron_number(105, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_TRUE(info.found);
    EXPECT_EQ(info.position_in_feature, 6);  // 105 - 100 + 1
}

TEST(ExonIntronNumber, PositionInIntron) {
    auto info = get_exon_intron_number(210, EXON_STARTS, EXON_ENDS, '+');
    EXPECT_TRUE(info.found);
    EXPECT_EQ(info.position_in_feature, 10);  // 210 - 201 + 1
}
