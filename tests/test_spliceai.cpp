/**
 * Tests for SpliceAI annotation source: field names, factory functions, cutoff
 */

#include <gtest/gtest.h>
#include "annotation_sources.hpp"
#include <algorithm>

using namespace vep;

// ============================================================================
// Factory function tests
// ============================================================================

TEST(SpliceAI, LegacyFactoryCreatesSource) {
    auto source = create_spliceai_source("/nonexistent/path.vcf.gz");
    ASSERT_NE(source, nullptr);
    EXPECT_EQ(source->name(), "spliceai");
    EXPECT_EQ(source->type(), "splice");
}

TEST(SpliceAI, NewFactoryCreatesSource) {
    auto source = create_spliceai_source("snv.vcf.gz", "indel.vcf.gz", "", 0.5);
    ASSERT_NE(source, nullptr);
    EXPECT_EQ(source->name(), "spliceai");
}

TEST(SpliceAI, EmptyPathsFactory) {
    auto source = create_spliceai_source("", "", "unified.vcf.gz", -1.0);
    ASSERT_NE(source, nullptr);
}

// ============================================================================
// Field name tests (must match Perl VEP SpliceAI plugin)
// ============================================================================

TEST(SpliceAI, FieldNamesMatchPerlVEP) {
    auto source = create_spliceai_source("/nonexistent.vcf.gz");
    auto fields = source->get_fields();

    // Check all expected fields exist
    auto has_field = [&](const std::string& name) {
        return std::find(fields.begin(), fields.end(), name) != fields.end();
    };

    EXPECT_TRUE(has_field("SpliceAI_pred"));
    EXPECT_TRUE(has_field("SpliceAI_pred_SYMBOL"));
    EXPECT_TRUE(has_field("SpliceAI_pred_DS_AG"));
    EXPECT_TRUE(has_field("SpliceAI_pred_DS_AL"));
    EXPECT_TRUE(has_field("SpliceAI_pred_DS_DG"));
    EXPECT_TRUE(has_field("SpliceAI_pred_DS_DL"));
    EXPECT_TRUE(has_field("SpliceAI_pred_DP_AG"));
    EXPECT_TRUE(has_field("SpliceAI_pred_DP_AL"));
    EXPECT_TRUE(has_field("SpliceAI_pred_DP_DG"));
    EXPECT_TRUE(has_field("SpliceAI_pred_DP_DL"));
    EXPECT_TRUE(has_field("SpliceAI_pred_DS_max"));
}

TEST(SpliceAI, NoCutoffFieldByDefault) {
    auto source = create_spliceai_source("/nonexistent.vcf.gz");
    auto fields = source->get_fields();

    auto has_field = [&](const std::string& name) {
        return std::find(fields.begin(), fields.end(), name) != fields.end();
    };

    // No cutoff field when cutoff is disabled (default -1.0)
    EXPECT_FALSE(has_field("SpliceAI_cutoff"));
}

TEST(SpliceAI, CutoffFieldWhenEnabled) {
    auto source = create_spliceai_source("", "", "/nonexistent.vcf.gz", 0.5);
    auto fields = source->get_fields();

    auto has_field = [&](const std::string& name) {
        return std::find(fields.begin(), fields.end(), name) != fields.end();
    };

    EXPECT_TRUE(has_field("SpliceAI_cutoff"));
}

TEST(SpliceAI, FieldCount) {
    // Without cutoff: 11 fields (combined + SYMBOL + 4 DS + 4 DP + DS_max)
    auto source_no_cutoff = create_spliceai_source("/nonexistent.vcf.gz");
    EXPECT_EQ(source_no_cutoff->get_fields().size(), 11u);

    // With cutoff: 12 fields (+ SpliceAI_cutoff)
    auto source_with_cutoff = create_spliceai_source("", "", "/nonexistent.vcf.gz", 0.2);
    EXPECT_EQ(source_with_cutoff->get_fields().size(), 12u);
}

// ============================================================================
// Source properties
// ============================================================================

TEST(SpliceAI, ThreadSafe) {
    auto source = create_spliceai_source("/nonexistent.vcf.gz");
    EXPECT_TRUE(source->is_thread_safe());
}

TEST(SpliceAI, Description) {
    auto source = create_spliceai_source("/nonexistent.vcf.gz");
    EXPECT_FALSE(source->description().empty());
}

// ============================================================================
// Data path reporting
// ============================================================================

TEST(SpliceAI, DataPathUnified) {
    auto source = create_spliceai_source("/path/to/unified.vcf.gz");
    EXPECT_NE(source->get_data_path().find("unified.vcf.gz"), std::string::npos);
}

TEST(SpliceAI, DataPathSeparateFiles) {
    auto source = create_spliceai_source("snv.vcf.gz", "indel.vcf.gz", "", 0.5);
    auto path = source->get_data_path();
    EXPECT_NE(path.find("snv=snv.vcf.gz"), std::string::npos);
    EXPECT_NE(path.find("indel=indel.vcf.gz"), std::string::npos);
}

// ============================================================================
// Not ready without valid files
// ============================================================================

TEST(SpliceAI, NotReadyWithoutFiles) {
    auto source = create_spliceai_source("/nonexistent.vcf.gz");
    // Source is not ready until initialize() is called and files are valid
    EXPECT_FALSE(source->is_ready());
}
