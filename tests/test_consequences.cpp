/**
 * Tests for consequence type handling: enum, string conversion, impact, ranking
 */

#include <gtest/gtest.h>
#include "vep_annotator.hpp"
#include "transcript_filter.hpp"

using namespace vep;

// ============================================================================
// ConsequenceType to string conversion
// ============================================================================

TEST(ConsequenceString, HighImpactTypes) {
    EXPECT_EQ(consequence_to_string(ConsequenceType::TRANSCRIPT_ABLATION), "transcript_ablation");
    EXPECT_EQ(consequence_to_string(ConsequenceType::SPLICE_ACCEPTOR_VARIANT), "splice_acceptor_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::SPLICE_DONOR_VARIANT), "splice_donor_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::STOP_GAINED), "stop_gained");
    EXPECT_EQ(consequence_to_string(ConsequenceType::FRAMESHIFT_VARIANT), "frameshift_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::STOP_LOST), "stop_lost");
    EXPECT_EQ(consequence_to_string(ConsequenceType::START_LOST), "start_lost");
    EXPECT_EQ(consequence_to_string(ConsequenceType::TRANSCRIPT_AMPLIFICATION), "transcript_amplification");
    EXPECT_EQ(consequence_to_string(ConsequenceType::FEATURE_ELONGATION), "feature_elongation");
    EXPECT_EQ(consequence_to_string(ConsequenceType::FEATURE_TRUNCATION), "feature_truncation");
}

TEST(ConsequenceString, ModerateImpactTypes) {
    EXPECT_EQ(consequence_to_string(ConsequenceType::INFRAME_INSERTION), "inframe_insertion");
    EXPECT_EQ(consequence_to_string(ConsequenceType::INFRAME_DELETION), "inframe_deletion");
    EXPECT_EQ(consequence_to_string(ConsequenceType::MISSENSE_VARIANT), "missense_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::PROTEIN_ALTERING_VARIANT), "protein_altering_variant");
}

TEST(ConsequenceString, LowImpactTypes) {
    EXPECT_EQ(consequence_to_string(ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT), "splice_donor_5th_base_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::SPLICE_DONOR_REGION_VARIANT), "splice_donor_region_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT), "splice_polypyrimidine_tract_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::SPLICE_REGION_VARIANT), "splice_region_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::INCOMPLETE_TERMINAL_CODON_VARIANT), "incomplete_terminal_codon_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::START_RETAINED_VARIANT), "start_retained_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::STOP_RETAINED_VARIANT), "stop_retained_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::SYNONYMOUS_VARIANT), "synonymous_variant");
}

TEST(ConsequenceString, ModifierImpactTypes) {
    EXPECT_EQ(consequence_to_string(ConsequenceType::CODING_SEQUENCE_VARIANT), "coding_sequence_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::CODING_TRANSCRIPT_VARIANT), "coding_transcript_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::MATURE_MIRNA_VARIANT), "mature_miRNA_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::FIVE_PRIME_UTR_VARIANT), "5_prime_UTR_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::THREE_PRIME_UTR_VARIANT), "3_prime_UTR_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT), "non_coding_transcript_exon_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::INTRON_VARIANT), "intron_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::NMD_TRANSCRIPT_VARIANT), "NMD_transcript_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::NON_CODING_TRANSCRIPT_VARIANT), "non_coding_transcript_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::UPSTREAM_GENE_VARIANT), "upstream_gene_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::DOWNSTREAM_GENE_VARIANT), "downstream_gene_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::TFBS_ABLATION), "TFBS_ablation");
    EXPECT_EQ(consequence_to_string(ConsequenceType::TFBS_AMPLIFICATION), "TFBS_amplification");
    EXPECT_EQ(consequence_to_string(ConsequenceType::TF_BINDING_SITE_VARIANT), "TF_binding_site_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::REGULATORY_REGION_ABLATION), "regulatory_region_ablation");
    EXPECT_EQ(consequence_to_string(ConsequenceType::REGULATORY_REGION_AMPLIFICATION), "regulatory_region_amplification");
    EXPECT_EQ(consequence_to_string(ConsequenceType::REGULATORY_REGION_VARIANT), "regulatory_region_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::INTERGENIC_VARIANT), "intergenic_variant");
    EXPECT_EQ(consequence_to_string(ConsequenceType::SEQUENCE_VARIANT), "sequence_variant");
}

// ============================================================================
// Impact levels
// ============================================================================

TEST(ImpactLevel, HighImpact) {
    EXPECT_EQ(get_impact(ConsequenceType::TRANSCRIPT_ABLATION), Impact::HIGH);
    EXPECT_EQ(get_impact(ConsequenceType::SPLICE_ACCEPTOR_VARIANT), Impact::HIGH);
    EXPECT_EQ(get_impact(ConsequenceType::SPLICE_DONOR_VARIANT), Impact::HIGH);
    EXPECT_EQ(get_impact(ConsequenceType::STOP_GAINED), Impact::HIGH);
    EXPECT_EQ(get_impact(ConsequenceType::FRAMESHIFT_VARIANT), Impact::HIGH);
    EXPECT_EQ(get_impact(ConsequenceType::STOP_LOST), Impact::HIGH);
    EXPECT_EQ(get_impact(ConsequenceType::START_LOST), Impact::HIGH);
    EXPECT_EQ(get_impact(ConsequenceType::TRANSCRIPT_AMPLIFICATION), Impact::HIGH);
    EXPECT_EQ(get_impact(ConsequenceType::FEATURE_ELONGATION), Impact::HIGH);
    EXPECT_EQ(get_impact(ConsequenceType::FEATURE_TRUNCATION), Impact::HIGH);
}

TEST(ImpactLevel, ModerateImpact) {
    EXPECT_EQ(get_impact(ConsequenceType::INFRAME_INSERTION), Impact::MODERATE);
    EXPECT_EQ(get_impact(ConsequenceType::INFRAME_DELETION), Impact::MODERATE);
    EXPECT_EQ(get_impact(ConsequenceType::MISSENSE_VARIANT), Impact::MODERATE);
    EXPECT_EQ(get_impact(ConsequenceType::PROTEIN_ALTERING_VARIANT), Impact::MODERATE);
}

TEST(ImpactLevel, LowImpact) {
    EXPECT_EQ(get_impact(ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT), Impact::LOW);
    EXPECT_EQ(get_impact(ConsequenceType::SPLICE_DONOR_REGION_VARIANT), Impact::LOW);
    EXPECT_EQ(get_impact(ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT), Impact::LOW);
    EXPECT_EQ(get_impact(ConsequenceType::SPLICE_REGION_VARIANT), Impact::LOW);
    EXPECT_EQ(get_impact(ConsequenceType::SYNONYMOUS_VARIANT), Impact::LOW);
}

TEST(ImpactLevel, ModifierImpact) {
    EXPECT_EQ(get_impact(ConsequenceType::INTRON_VARIANT), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::UPSTREAM_GENE_VARIANT), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::DOWNSTREAM_GENE_VARIANT), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::INTERGENIC_VARIANT), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::REGULATORY_REGION_VARIANT), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::TF_BINDING_SITE_VARIANT), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::CODING_TRANSCRIPT_VARIANT), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::SEQUENCE_VARIANT), Impact::MODIFIER);
}

// ============================================================================
// Impact string conversion
// ============================================================================

TEST(ImpactString, AllLevels) {
    EXPECT_EQ(impact_to_string(Impact::HIGH), "HIGH");
    EXPECT_EQ(impact_to_string(Impact::MODERATE), "MODERATE");
    EXPECT_EQ(impact_to_string(Impact::LOW), "LOW");
    EXPECT_EQ(impact_to_string(Impact::MODIFIER), "MODIFIER");
}

// ============================================================================
// Consequence ranking
// ============================================================================

TEST(ConsequenceRanking, SeverityOrder) {
    // HIGH impact types should rank lower (more severe) than MODIFIER
    EXPECT_LT(get_consequence_rank(ConsequenceType::TRANSCRIPT_ABLATION),
              get_consequence_rank(ConsequenceType::INTERGENIC_VARIANT));

    EXPECT_LT(get_consequence_rank(ConsequenceType::SPLICE_DONOR_VARIANT),
              get_consequence_rank(ConsequenceType::MISSENSE_VARIANT));

    EXPECT_LT(get_consequence_rank(ConsequenceType::MISSENSE_VARIANT),
              get_consequence_rank(ConsequenceType::SYNONYMOUS_VARIANT));

    EXPECT_LT(get_consequence_rank(ConsequenceType::SYNONYMOUS_VARIANT),
              get_consequence_rank(ConsequenceType::INTRON_VARIANT));

    EXPECT_LT(get_consequence_rank(ConsequenceType::INTRON_VARIANT),
              get_consequence_rank(ConsequenceType::INTERGENIC_VARIANT));
}

TEST(ConsequenceRanking, NewTypesInCorrectPosition) {
    // TRANSCRIPT_AMPLIFICATION is HIGH but less severe than START_LOST
    EXPECT_LT(get_consequence_rank(ConsequenceType::START_LOST),
              get_consequence_rank(ConsequenceType::TRANSCRIPT_AMPLIFICATION));

    // Splice sub-types must match Perl VEP ranking exactly:
    // 15: splice_donor_5th_base_variant
    // 16: splice_region_variant
    // 17: splice_donor_region_variant
    // 18: splice_polypyrimidine_tract_variant
    EXPECT_LT(get_consequence_rank(ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT),
              get_consequence_rank(ConsequenceType::SPLICE_REGION_VARIANT));
    EXPECT_LT(get_consequence_rank(ConsequenceType::SPLICE_REGION_VARIANT),
              get_consequence_rank(ConsequenceType::SPLICE_DONOR_REGION_VARIANT));
    EXPECT_LT(get_consequence_rank(ConsequenceType::SPLICE_DONOR_REGION_VARIANT),
              get_consequence_rank(ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT));

    // CODING_TRANSCRIPT_VARIANT is after NON_CODING_TRANSCRIPT_VARIANT (rank 31 in Perl VEP)
    EXPECT_GT(get_consequence_rank(ConsequenceType::CODING_TRANSCRIPT_VARIANT),
              get_consequence_rank(ConsequenceType::NON_CODING_TRANSCRIPT_VARIANT));

    // MATURE_MIRNA_VARIANT comes right after CODING_SEQUENCE_VARIANT (rank 24)
    EXPECT_EQ(get_consequence_rank(ConsequenceType::MATURE_MIRNA_VARIANT),
              get_consequence_rank(ConsequenceType::CODING_SEQUENCE_VARIANT) + 1);

    // REGULATORY types
    EXPECT_GT(get_consequence_rank(ConsequenceType::REGULATORY_REGION_VARIANT),
              get_consequence_rank(ConsequenceType::INTRON_VARIANT));

    // SEQUENCE_VARIANT is the most generic
    EXPECT_GT(get_consequence_rank(ConsequenceType::SEQUENCE_VARIANT),
              get_consequence_rank(ConsequenceType::INTERGENIC_VARIANT));
}

TEST(ConsequenceRanking, AllTypesHaveRank) {
    // Ensure no type returns the default rank (999)
    std::vector<ConsequenceType> all_types = {
        ConsequenceType::TRANSCRIPT_ABLATION,
        ConsequenceType::SPLICE_ACCEPTOR_VARIANT,
        ConsequenceType::SPLICE_DONOR_VARIANT,
        ConsequenceType::STOP_GAINED,
        ConsequenceType::FRAMESHIFT_VARIANT,
        ConsequenceType::STOP_LOST,
        ConsequenceType::START_LOST,
        ConsequenceType::TRANSCRIPT_AMPLIFICATION,
        ConsequenceType::FEATURE_ELONGATION,
        ConsequenceType::FEATURE_TRUNCATION,
        ConsequenceType::INFRAME_INSERTION,
        ConsequenceType::INFRAME_DELETION,
        ConsequenceType::MISSENSE_VARIANT,
        ConsequenceType::PROTEIN_ALTERING_VARIANT,
        ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT,
        ConsequenceType::SPLICE_DONOR_REGION_VARIANT,
        ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT,
        ConsequenceType::SPLICE_REGION_VARIANT,
        ConsequenceType::INCOMPLETE_TERMINAL_CODON_VARIANT,
        ConsequenceType::START_RETAINED_VARIANT,
        ConsequenceType::STOP_RETAINED_VARIANT,
        ConsequenceType::SYNONYMOUS_VARIANT,
        ConsequenceType::CODING_SEQUENCE_VARIANT,
        ConsequenceType::CODING_TRANSCRIPT_VARIANT,
        ConsequenceType::MATURE_MIRNA_VARIANT,
        ConsequenceType::FIVE_PRIME_UTR_VARIANT,
        ConsequenceType::THREE_PRIME_UTR_VARIANT,
        ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT,
        ConsequenceType::INTRON_VARIANT,
        ConsequenceType::NMD_TRANSCRIPT_VARIANT,
        ConsequenceType::NON_CODING_TRANSCRIPT_VARIANT,
        ConsequenceType::UPSTREAM_GENE_VARIANT,
        ConsequenceType::DOWNSTREAM_GENE_VARIANT,
        ConsequenceType::TFBS_ABLATION,
        ConsequenceType::TFBS_AMPLIFICATION,
        ConsequenceType::TF_BINDING_SITE_VARIANT,
        ConsequenceType::REGULATORY_REGION_ABLATION,
        ConsequenceType::REGULATORY_REGION_AMPLIFICATION,
        ConsequenceType::REGULATORY_REGION_VARIANT,
        ConsequenceType::INTERGENIC_VARIANT,
        ConsequenceType::SEQUENCE_VARIANT,
    };

    for (const auto& type : all_types) {
        int rank = get_consequence_rank(type);
        EXPECT_LT(rank, 999) << "Type " << consequence_to_string(type) << " has no rank";
    }
}

// ============================================================================
// VariantAnnotation consequence string
// ============================================================================

TEST(VariantAnnotation, ConsequenceString) {
    VariantAnnotation ann;
    ann.consequences.push_back(ConsequenceType::MISSENSE_VARIANT);
    EXPECT_EQ(ann.get_consequence_string(), "missense_variant");
}

TEST(VariantAnnotation, MultipleConsequences) {
    VariantAnnotation ann;
    ann.consequences.push_back(ConsequenceType::SPLICE_REGION_VARIANT);
    ann.consequences.push_back(ConsequenceType::INTRON_VARIANT);
    std::string result = ann.get_consequence_string();
    // Should contain both consequences separated by &
    EXPECT_TRUE(result.find("splice_region_variant") != std::string::npos);
    EXPECT_TRUE(result.find("intron_variant") != std::string::npos);
}

// ============================================================================
// Variant class determination
// ============================================================================

TEST(VariantClass, SNV) {
    EXPECT_EQ(get_variant_class("A", "T"), "SNV");
    EXPECT_EQ(get_variant_class("C", "G"), "SNV");
}

TEST(VariantClass, Substitution) {
    EXPECT_EQ(get_variant_class("AT", "CG"), "substitution");
    EXPECT_EQ(get_variant_class("ACG", "TCA"), "substitution");
}

TEST(VariantClass, Insertion) {
    EXPECT_EQ(get_variant_class("A", "AT"), "insertion");
    EXPECT_EQ(get_variant_class("A", "ATCG"), "insertion");
}

TEST(VariantClass, Deletion) {
    EXPECT_EQ(get_variant_class("AT", "A"), "deletion");
    EXPECT_EQ(get_variant_class("ATCG", "A"), "deletion");
}

TEST(VariantClass, Indel) {
    EXPECT_EQ(get_variant_class("ATC", "AG"), "indel");
    EXPECT_EQ(get_variant_class("AG", "ATC"), "indel");
}

// ============================================================================
// Display term conversion
// ============================================================================

TEST(DisplayTerm, HighImpactTerms) {
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::STOP_GAINED), "STOP_GAINED");
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::FRAMESHIFT_VARIANT), "FRAMESHIFT_CODING");
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::SPLICE_ACCEPTOR_VARIANT), "ESSENTIAL_SPLICE_SITE");
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::SPLICE_DONOR_VARIANT), "ESSENTIAL_SPLICE_SITE");
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::START_LOST), "NON_SYNONYMOUS_CODING");
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::STOP_LOST), "STOP_LOST");
}

TEST(DisplayTerm, ModerateImpactTerms) {
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::MISSENSE_VARIANT), "NON_SYNONYMOUS_CODING");
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::INFRAME_INSERTION), "NON_SYNONYMOUS_CODING");
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::INFRAME_DELETION), "NON_SYNONYMOUS_CODING");
}

TEST(DisplayTerm, LowImpactTerms) {
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::SYNONYMOUS_VARIANT), "SYNONYMOUS_CODING");
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::SPLICE_REGION_VARIANT), "SPLICE_SITE");
}

TEST(DisplayTerm, ModifierImpactTerms) {
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::INTRON_VARIANT), "INTRONIC");
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::UPSTREAM_GENE_VARIANT), "UPSTREAM");
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::DOWNSTREAM_GENE_VARIANT), "DOWNSTREAM");
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::INTERGENIC_VARIANT), "INTERGENIC");
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::FIVE_PRIME_UTR_VARIANT), "5PRIME_UTR");
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::THREE_PRIME_UTR_VARIANT), "3PRIME_UTR");
}

// ============================================================================
// Exon/intron X/TOTAL format in VariantAnnotation
// ============================================================================

TEST(VariantAnnotation, ExonIntronTotalFields) {
    VariantAnnotation ann;
    ann.exon_number = 3;
    ann.total_exons = 10;
    ann.intron_number = 0;
    ann.total_introns = 9;

    EXPECT_EQ(ann.exon_number, 3);
    EXPECT_EQ(ann.total_exons, 10);
    EXPECT_EQ(ann.intron_number, 0);
    EXPECT_EQ(ann.total_introns, 9);
}

// ============================================================================
// Pick order default (MANE before CANONICAL)
// ============================================================================

TEST(PickOrderDefault, MANEBeforeCanonical) {
    TranscriptFilterConfig config;
    EXPECT_EQ(config.pick_order.size(), 9u);
    EXPECT_EQ(config.pick_order[0], PickCriteria::MANE_SELECT);
    EXPECT_EQ(config.pick_order[1], PickCriteria::MANE_PLUS);
    EXPECT_EQ(config.pick_order[2], PickCriteria::CANONICAL);
}
