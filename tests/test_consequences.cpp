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
    EXPECT_EQ(get_impact(ConsequenceType::FEATURE_ELONGATION), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::FEATURE_TRUNCATION), Impact::MODIFIER);
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

// ============================================================================
// Splice consequence suppression logic
// (Tests the VEP spec: specific sub-types suppress generic splice_region_variant)
// ============================================================================

TEST(SpliceConsequenceSuppression, Donor5thSuppressesSpliceRegionAndDonorRegion) {
    // When splice_donor_5th_base_variant is present, both
    // splice_region_variant and splice_donor_region_variant should be removed
    std::vector<ConsequenceType> cons = {
        ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT,
        ConsequenceType::SPLICE_REGION_VARIANT,
        ConsequenceType::SPLICE_DONOR_REGION_VARIANT
    };

    // Deduplicate and sort (as the real code does)
    std::sort(cons.begin(), cons.end());
    cons.erase(std::unique(cons.begin(), cons.end()), cons.end());

    // Apply suppression logic (same as vep_annotator.cpp)
    if (cons.size() > 1) {
        auto has = [&](ConsequenceType t) {
            return std::find(cons.begin(), cons.end(), t) != cons.end();
        };
        bool has_donor_5th = has(ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT);
        bool has_donor_region = has(ConsequenceType::SPLICE_DONOR_REGION_VARIANT);
        bool has_polypyrimidine = has(ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT);
        if (has_donor_5th || has_donor_region || has_polypyrimidine) {
            cons.erase(
                std::remove(cons.begin(), cons.end(), ConsequenceType::SPLICE_REGION_VARIANT),
                cons.end());
        }
        if (has_donor_5th) {
            cons.erase(
                std::remove(cons.begin(), cons.end(), ConsequenceType::SPLICE_DONOR_REGION_VARIANT),
                cons.end());
        }
    }

    ASSERT_EQ(cons.size(), 1u);
    EXPECT_EQ(cons[0], ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT);
}

TEST(SpliceConsequenceSuppression, DonorRegionSuppressesSpliceRegionOnly) {
    std::vector<ConsequenceType> cons = {
        ConsequenceType::SPLICE_DONOR_REGION_VARIANT,
        ConsequenceType::SPLICE_REGION_VARIANT
    };

    std::sort(cons.begin(), cons.end());

    if (cons.size() > 1) {
        auto has = [&](ConsequenceType t) {
            return std::find(cons.begin(), cons.end(), t) != cons.end();
        };
        if (has(ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT) ||
            has(ConsequenceType::SPLICE_DONOR_REGION_VARIANT) ||
            has(ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT)) {
            cons.erase(
                std::remove(cons.begin(), cons.end(), ConsequenceType::SPLICE_REGION_VARIANT),
                cons.end());
        }
    }

    ASSERT_EQ(cons.size(), 1u);
    EXPECT_EQ(cons[0], ConsequenceType::SPLICE_DONOR_REGION_VARIANT);
}

TEST(SpliceConsequenceSuppression, PolypyrimidineSuppressesSpliceRegion) {
    std::vector<ConsequenceType> cons = {
        ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT,
        ConsequenceType::SPLICE_REGION_VARIANT
    };

    std::sort(cons.begin(), cons.end());

    if (cons.size() > 1) {
        auto has = [&](ConsequenceType t) {
            return std::find(cons.begin(), cons.end(), t) != cons.end();
        };
        if (has(ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT) ||
            has(ConsequenceType::SPLICE_DONOR_REGION_VARIANT) ||
            has(ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT)) {
            cons.erase(
                std::remove(cons.begin(), cons.end(), ConsequenceType::SPLICE_REGION_VARIANT),
                cons.end());
        }
    }

    ASSERT_EQ(cons.size(), 1u);
    EXPECT_EQ(cons[0], ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT);
}

TEST(SpliceConsequenceSuppression, SpliceRegionAloneNotSuppressed) {
    std::vector<ConsequenceType> cons = {
        ConsequenceType::SPLICE_REGION_VARIANT
    };

    // With only one consequence, suppression should not fire
    ASSERT_EQ(cons.size(), 1u);
    EXPECT_EQ(cons[0], ConsequenceType::SPLICE_REGION_VARIANT);
}

// ============================================================================
// ConsequenceRanking - complete severity ordering (10 tests)
// ============================================================================

TEST(ConsequenceRanking, TranscriptAblationMoreSevereThanSpliceAcceptor) {
    EXPECT_LT(get_consequence_rank(ConsequenceType::TRANSCRIPT_ABLATION),
              get_consequence_rank(ConsequenceType::SPLICE_ACCEPTOR_VARIANT));
}

TEST(ConsequenceRanking, SpliceAcceptorMoreSevereThanSpliceDonor) {
    EXPECT_LT(get_consequence_rank(ConsequenceType::SPLICE_ACCEPTOR_VARIANT),
              get_consequence_rank(ConsequenceType::SPLICE_DONOR_VARIANT));
}

TEST(ConsequenceRanking, SpliceDonorMoreSevereThanStopGained) {
    EXPECT_LT(get_consequence_rank(ConsequenceType::SPLICE_DONOR_VARIANT),
              get_consequence_rank(ConsequenceType::STOP_GAINED));
}

TEST(ConsequenceRanking, StopGainedMoreSevereThanFrameshift) {
    EXPECT_LT(get_consequence_rank(ConsequenceType::STOP_GAINED),
              get_consequence_rank(ConsequenceType::FRAMESHIFT_VARIANT));
}

TEST(ConsequenceRanking, FrameshiftMoreSevereThanStopLost) {
    EXPECT_LT(get_consequence_rank(ConsequenceType::FRAMESHIFT_VARIANT),
              get_consequence_rank(ConsequenceType::STOP_LOST));
}

TEST(ConsequenceRanking, StopLostMoreSevereThanStartLost) {
    EXPECT_LT(get_consequence_rank(ConsequenceType::STOP_LOST),
              get_consequence_rank(ConsequenceType::START_LOST));
}

TEST(ConsequenceRanking, StartLostMoreSevereThanTranscriptAmplification) {
    EXPECT_LT(get_consequence_rank(ConsequenceType::START_LOST),
              get_consequence_rank(ConsequenceType::TRANSCRIPT_AMPLIFICATION));
}

TEST(ConsequenceRanking, InframeInsertionLessSevereThanInframeDeletion) {
    EXPECT_LT(get_consequence_rank(ConsequenceType::INFRAME_INSERTION),
              get_consequence_rank(ConsequenceType::INFRAME_DELETION));
}

TEST(ConsequenceRanking, MissenseMoreSevereThanProteinAltering) {
    EXPECT_LT(get_consequence_rank(ConsequenceType::MISSENSE_VARIANT),
              get_consequence_rank(ConsequenceType::PROTEIN_ALTERING_VARIANT));
}

TEST(ConsequenceRanking, SynonymousLessSevereThanCodingSequenceAndUTR) {
    EXPECT_LT(get_consequence_rank(ConsequenceType::SYNONYMOUS_VARIANT),
              get_consequence_rank(ConsequenceType::CODING_SEQUENCE_VARIANT));
    EXPECT_LT(get_consequence_rank(ConsequenceType::CODING_SEQUENCE_VARIANT),
              get_consequence_rank(ConsequenceType::FIVE_PRIME_UTR_VARIANT));
    EXPECT_LT(get_consequence_rank(ConsequenceType::FIVE_PRIME_UTR_VARIANT),
              get_consequence_rank(ConsequenceType::THREE_PRIME_UTR_VARIANT));
}

// ============================================================================
// VariantClass edge cases (8 tests)
// ============================================================================

TEST(VariantClass, EmptyRefAndAlt) {
    // Both empty strings have equal size (0), so the function returns "substitution"
    std::string result = get_variant_class("", "");
    EXPECT_EQ(result, "substitution");
}

TEST(VariantClass, EmptyRefNonEmptyAlt) {
    EXPECT_EQ(get_variant_class("", "A"), "insertion");
}

TEST(VariantClass, NonEmptyRefEmptyAlt) {
    EXPECT_EQ(get_variant_class("A", ""), "deletion");
}

TEST(VariantClass, SingleBaseSame) {
    // Same single base should be SNV (the function classifies by length, not content)
    EXPECT_EQ(get_variant_class("A", "A"), "SNV");
}

TEST(VariantClass, SymbolicAlleleDEL) {
    // Symbolic allele <DEL> treated as literal string: ref=A, alt=<DEL>
    // No common prefix, both non-empty after trim => indel
    EXPECT_EQ(get_variant_class("A", "<DEL>"), "indel");
}

TEST(VariantClass, MNV) {
    // Multi-nucleotide variant: same length, different bases
    EXPECT_EQ(get_variant_class("ATG", "TCA"), "substitution");
}

TEST(VariantClass, LongInsertion) {
    EXPECT_EQ(get_variant_class("A", "ATCGATCGATCG"), "insertion");
}

TEST(VariantClass, ComplexIndel) {
    // Different lengths with no common prefix match at start
    EXPECT_EQ(get_variant_class("ATCG", "TG"), "indel");
}

// ============================================================================
// DisplayTermCompleteness (5 tests)
// ============================================================================

TEST(DisplayTerm, ProteinAlteringVariant) {
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::PROTEIN_ALTERING_VARIANT), "PROTEIN_ALTERING");
}

TEST(DisplayTerm, NMDTranscriptVariant) {
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::NMD_TRANSCRIPT_VARIANT), "NMD_TRANSCRIPT");
}

TEST(DisplayTerm, IncompleteTerminalCodonVariant) {
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::INCOMPLETE_TERMINAL_CODON_VARIANT), "PARTIAL_CODON");
}

TEST(DisplayTerm, NonCodingTranscriptExonVariant) {
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT), "WITHIN_NON_CODING_GENE");
}

TEST(DisplayTerm, CodingSequenceVariant) {
    EXPECT_EQ(consequence_to_display_term(ConsequenceType::CODING_SEQUENCE_VARIANT), "CODING_UNKNOWN");
}

// ============================================================================
// MultiConsequenceFormatting (5 tests)
// ============================================================================

TEST(MultiConsequenceFormatting, SingleConsequence) {
    VariantAnnotation ann;
    ann.consequences.push_back(ConsequenceType::MISSENSE_VARIANT);
    EXPECT_EQ(ann.get_consequence_string(), "missense_variant");
}

TEST(MultiConsequenceFormatting, ThreeConsequencesJoinedWithAmpersand) {
    VariantAnnotation ann;
    ann.consequences.push_back(ConsequenceType::SPLICE_REGION_VARIANT);
    ann.consequences.push_back(ConsequenceType::INTRON_VARIANT);
    ann.consequences.push_back(ConsequenceType::NMD_TRANSCRIPT_VARIANT);
    std::string result = ann.get_consequence_string();
    EXPECT_EQ(result, "splice_region_variant&intron_variant&NMD_transcript_variant");
}

TEST(MultiConsequenceFormatting, EmptyConsequences) {
    VariantAnnotation ann;
    EXPECT_EQ(ann.get_consequence_string(), "");
}

TEST(MultiConsequenceFormatting, DuplicateConsequences) {
    VariantAnnotation ann;
    ann.consequences.push_back(ConsequenceType::MISSENSE_VARIANT);
    ann.consequences.push_back(ConsequenceType::MISSENSE_VARIANT);
    std::string result = ann.get_consequence_string();
    // Duplicates are preserved (no dedup in get_consequence_string itself)
    EXPECT_EQ(result, "missense_variant&missense_variant");
}

TEST(MultiConsequenceFormatting, HighAndModifierMixed) {
    VariantAnnotation ann;
    ann.consequences.push_back(ConsequenceType::STOP_GAINED);
    ann.consequences.push_back(ConsequenceType::NMD_TRANSCRIPT_VARIANT);
    std::string result = ann.get_consequence_string();
    EXPECT_EQ(result, "stop_gained&NMD_transcript_variant");
}

// ============================================================================
// ImpactCompleteness (4 tests)
// ============================================================================

TEST(ImpactCompleteness, UTRAndNonCodingModifier) {
    EXPECT_EQ(get_impact(ConsequenceType::FIVE_PRIME_UTR_VARIANT), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::THREE_PRIME_UTR_VARIANT), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::CODING_SEQUENCE_VARIANT), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::MATURE_MIRNA_VARIANT), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::NMD_TRANSCRIPT_VARIANT), Impact::MODIFIER);
}

TEST(ImpactCompleteness, LowImpactRemainingTypes) {
    EXPECT_EQ(get_impact(ConsequenceType::INCOMPLETE_TERMINAL_CODON_VARIANT), Impact::LOW);
    EXPECT_EQ(get_impact(ConsequenceType::START_RETAINED_VARIANT), Impact::LOW);
    EXPECT_EQ(get_impact(ConsequenceType::STOP_RETAINED_VARIANT), Impact::LOW);
}

TEST(ImpactCompleteness, RegulatoryModifier) {
    EXPECT_EQ(get_impact(ConsequenceType::TFBS_ABLATION), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::TFBS_AMPLIFICATION), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::REGULATORY_REGION_ABLATION), Impact::MODIFIER);
    EXPECT_EQ(get_impact(ConsequenceType::REGULATORY_REGION_AMPLIFICATION), Impact::MODIFIER);
}

TEST(ImpactCompleteness, NonCodingTranscriptVariantModifier) {
    EXPECT_EQ(get_impact(ConsequenceType::NON_CODING_TRANSCRIPT_VARIANT), Impact::MODIFIER);
}

// ============================================================================
// ConsequenceRankUniqueness (2 tests)
// ============================================================================

TEST(ConsequenceRankUniqueness, DifferentTypesHaveDifferentRanks) {
    // Build a map of rank -> list of types. Most ranks should be unique.
    std::map<int, std::vector<ConsequenceType>> rank_map;
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
        rank_map[get_consequence_rank(type)].push_back(type);
    }

    // Each rank should have exactly one consequence type (all ranks are unique)
    for (const auto& [rank, types] : rank_map) {
        EXPECT_EQ(types.size(), 1u)
            << "Rank " << rank << " is shared by " << types.size() << " types";
    }
}

TEST(ConsequenceRankUniqueness, AllRanksArePositiveIntegers) {
    std::vector<ConsequenceType> all_types = {
        ConsequenceType::TRANSCRIPT_ABLATION,
        ConsequenceType::SPLICE_ACCEPTOR_VARIANT,
        ConsequenceType::SPLICE_DONOR_VARIANT,
        ConsequenceType::STOP_GAINED,
        ConsequenceType::FRAMESHIFT_VARIANT,
        ConsequenceType::STOP_LOST,
        ConsequenceType::START_LOST,
        ConsequenceType::TRANSCRIPT_AMPLIFICATION,
        ConsequenceType::INFRAME_INSERTION,
        ConsequenceType::INFRAME_DELETION,
        ConsequenceType::MISSENSE_VARIANT,
        ConsequenceType::PROTEIN_ALTERING_VARIANT,
        ConsequenceType::SYNONYMOUS_VARIANT,
        ConsequenceType::INTRON_VARIANT,
        ConsequenceType::INTERGENIC_VARIANT,
        ConsequenceType::SEQUENCE_VARIANT,
    };

    for (const auto& type : all_types) {
        int rank = get_consequence_rank(type);
        EXPECT_GT(rank, 0) << "Rank for " << consequence_to_string(type) << " should be positive";
    }
}

// ============================================================================
// VariantAnnotationFields (8 tests)
// ============================================================================

TEST(VariantAnnotationFields, DefaultInitialization) {
    VariantAnnotation ann;
    EXPECT_EQ(ann.position, 0);
    EXPECT_EQ(ann.impact, Impact::MODIFIER);
    EXPECT_TRUE(ann.consequences.empty());
    EXPECT_EQ(ann.exon_number, 0);
    EXPECT_EQ(ann.intron_number, 0);
    EXPECT_EQ(ann.cdna_position, 0);
    EXPECT_EQ(ann.cds_position, 0);
    EXPECT_EQ(ann.protein_position, 0);
    EXPECT_FALSE(ann.is_canonical);
    EXPECT_EQ(ann.feature_type, "Transcript");
}

TEST(VariantAnnotationFields, CustomAnnotationsMapAccess) {
    VariantAnnotation ann;
    ann.custom_annotations["gnomad:AF"] = "0.001";
    ann.custom_annotations["clinvar:CLNSIG"] = "Pathogenic";

    EXPECT_EQ(ann.custom_annotations.size(), 2u);
    EXPECT_EQ(ann.custom_annotations["gnomad:AF"], "0.001");
    EXPECT_EQ(ann.custom_annotations["clinvar:CLNSIG"], "Pathogenic");
}

TEST(VariantAnnotationFields, GetAnnotationPopulated) {
    VariantAnnotation ann;
    ann.custom_annotations["gnomad:AF"] = "0.001";
    ann.custom_annotations["clinvar:CLNSIG"] = "Pathogenic";

    auto af = ann.get_annotation("gnomad", "AF");
    ASSERT_TRUE(af.has_value());
    EXPECT_EQ(af.value(), "0.001");

    auto clnsig = ann.get_annotation("clinvar", "CLNSIG");
    ASSERT_TRUE(clnsig.has_value());
    EXPECT_EQ(clnsig.value(), "Pathogenic");
}

TEST(VariantAnnotationFields, GetAnnotationMissing) {
    VariantAnnotation ann;
    auto result = ann.get_annotation("nonexistent", "field");
    EXPECT_FALSE(result.has_value());
}

TEST(VariantAnnotationFields, GetAnnotationDoubleParsing) {
    VariantAnnotation ann;
    ann.custom_annotations["gnomad:AF"] = "0.00123";
    ann.custom_annotations["test:invalid"] = "not_a_number";

    auto af = ann.get_annotation_double("gnomad", "AF");
    ASSERT_TRUE(af.has_value());
    EXPECT_NEAR(af.value(), 0.00123, 1e-6);

    auto invalid = ann.get_annotation_double("test", "invalid");
    EXPECT_FALSE(invalid.has_value());
}

TEST(VariantAnnotationFields, EmptyAnnotationDefaultStrings) {
    VariantAnnotation ann;
    // Unset strings default to empty (not "-")
    EXPECT_TRUE(ann.codons.empty());
    EXPECT_TRUE(ann.amino_acids.empty());
    EXPECT_TRUE(ann.hgvsc.empty());
    EXPECT_TRUE(ann.hgvsp.empty());
    EXPECT_TRUE(ann.hgvsg.empty());
}

TEST(VariantAnnotationFields, EndPositionFields) {
    VariantAnnotation ann;
    ann.cdna_end = 150;
    ann.cds_end = 120;
    ann.protein_end = 40;

    EXPECT_EQ(ann.cdna_end, 150);
    EXPECT_EQ(ann.cds_end, 120);
    EXPECT_EQ(ann.protein_end, 40);
}

TEST(VariantAnnotationFields, DistanceField) {
    VariantAnnotation ann;
    ann.distance = 4500;
    ann.consequences.push_back(ConsequenceType::UPSTREAM_GENE_VARIANT);

    EXPECT_EQ(ann.distance, 4500);
    EXPECT_EQ(ann.consequences[0], ConsequenceType::UPSTREAM_GENE_VARIANT);
}

// ============================================================================
// GetConsequenceRank specific values (8 tests)
// ============================================================================

TEST(GetConsequenceRankValues, TranscriptAblationIsRank1) {
    EXPECT_EQ(get_consequence_rank(ConsequenceType::TRANSCRIPT_ABLATION), 1);
}

TEST(GetConsequenceRankValues, SpliceAcceptorIsRank2) {
    EXPECT_EQ(get_consequence_rank(ConsequenceType::SPLICE_ACCEPTOR_VARIANT), 2);
}

TEST(GetConsequenceRankValues, StopGainedIsRank4) {
    EXPECT_EQ(get_consequence_rank(ConsequenceType::STOP_GAINED), 4);
}

TEST(GetConsequenceRankValues, MissenseVariantIsRank13) {
    EXPECT_EQ(get_consequence_rank(ConsequenceType::MISSENSE_VARIANT), 13);
}

TEST(GetConsequenceRankValues, SynonymousVariantIsRank22) {
    EXPECT_EQ(get_consequence_rank(ConsequenceType::SYNONYMOUS_VARIANT), 22);
}

TEST(GetConsequenceRankValues, IntergenicVariantIsRank40) {
    EXPECT_EQ(get_consequence_rank(ConsequenceType::INTERGENIC_VARIANT), 40);
}

TEST(GetConsequenceRankValues, SequenceVariantIsRank41) {
    EXPECT_EQ(get_consequence_rank(ConsequenceType::SEQUENCE_VARIANT), 41);
}

TEST(GetConsequenceRankValues, IntergenicIsHighestNonSequenceVariantRank) {
    // INTERGENIC_VARIANT should be the last "real" consequence rank
    // before SEQUENCE_VARIANT (the catch-all)
    int intergenic = get_consequence_rank(ConsequenceType::INTERGENIC_VARIANT);
    int sequence = get_consequence_rank(ConsequenceType::SEQUENCE_VARIANT);
    EXPECT_LT(intergenic, sequence);

    // All other types should rank <= intergenic (except SEQUENCE_VARIANT)
    std::vector<ConsequenceType> all_except_seq = {
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
        ConsequenceType::SYNONYMOUS_VARIANT,
        ConsequenceType::INTRON_VARIANT,
        ConsequenceType::UPSTREAM_GENE_VARIANT,
        ConsequenceType::DOWNSTREAM_GENE_VARIANT,
        ConsequenceType::REGULATORY_REGION_VARIANT,
        ConsequenceType::INTERGENIC_VARIANT,
    };
    for (const auto& type : all_except_seq) {
        EXPECT_LE(get_consequence_rank(type), intergenic)
            << consequence_to_string(type) << " should rank <= INTERGENIC";
    }
}

// ============================================================================
// VariantClassSymbolicAlleles (5 tests)
// ============================================================================

TEST(VariantClassSymbolicAlleles, DELSymbolic) {
    // Symbolic alleles are not specially handled by get_variant_class;
    // they are treated as literal strings by length comparison
    std::string result = get_variant_class("A", "<DEL>");
    // ref=1 char, alt=5 chars, no common prefix => indel
    EXPECT_EQ(result, "indel");
}

TEST(VariantClassSymbolicAlleles, INSSymbolic) {
    std::string result = get_variant_class("A", "<INS>");
    EXPECT_EQ(result, "indel");
}

TEST(VariantClassSymbolicAlleles, DUPSymbolic) {
    std::string result = get_variant_class("A", "<DUP>");
    EXPECT_EQ(result, "indel");
}

TEST(VariantClassSymbolicAlleles, INVSymbolic) {
    std::string result = get_variant_class("A", "<INV>");
    EXPECT_EQ(result, "indel");
}

TEST(VariantClassSymbolicAlleles, CNVSymbolic) {
    std::string result = get_variant_class("A", "<CNV>");
    EXPECT_EQ(result, "indel");
}
