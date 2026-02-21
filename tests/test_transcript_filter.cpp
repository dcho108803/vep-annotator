/**
 * Tests for TranscriptFilter: basic filtering, pick modes, pick ordering,
 * frequency filtering, and edge cases.
 */

#include <gtest/gtest.h>
#include "vep_annotator.hpp"
#include "transcript_filter.hpp"

using namespace vep;

// ============================================================================
// Helper functions to build VariantAnnotation objects for testing
// ============================================================================

/**
 * Build a minimal VariantAnnotation with the given parameters.
 */
static VariantAnnotation make_annotation(
    const std::string& transcript_id,
    const std::string& gene_symbol,
    const std::string& gene_id,
    const std::string& biotype,
    bool is_canonical,
    const std::vector<ConsequenceType>& consequences,
    Impact impact,
    const std::string& ref_allele = "A",
    const std::string& alt_allele = "T")
{
    VariantAnnotation ann;
    ann.transcript_id = transcript_id;
    ann.gene_symbol = gene_symbol;
    ann.gene_id = gene_id;
    ann.biotype = biotype;
    ann.is_canonical = is_canonical;
    ann.consequences = consequences;
    ann.impact = impact;
    ann.ref_allele = ref_allele;
    ann.alt_allele = alt_allele;
    ann.chromosome = "chr17";
    ann.position = 7675088;
    return ann;
}

/**
 * Build annotation with MANE_SELECT metadata.
 */
static VariantAnnotation make_mane_annotation(
    const std::string& transcript_id,
    const std::string& gene_symbol,
    const std::string& gene_id,
    const std::string& biotype,
    bool is_canonical,
    const std::vector<ConsequenceType>& consequences,
    Impact impact,
    const std::string& mane_select_value = "NM_000546.6")
{
    VariantAnnotation ann = make_annotation(
        transcript_id, gene_symbol, gene_id, biotype, is_canonical, consequences, impact);
    ann.custom_annotations["MANE_SELECT"] = mane_select_value;
    return ann;
}

/**
 * Build annotation with GENCODE_BASIC metadata.
 */
static VariantAnnotation make_gencode_basic_annotation(
    const std::string& transcript_id,
    const std::string& gene_symbol,
    const std::string& gene_id,
    const std::string& biotype,
    bool is_canonical,
    const std::vector<ConsequenceType>& consequences,
    Impact impact)
{
    VariantAnnotation ann = make_annotation(
        transcript_id, gene_symbol, gene_id, biotype, is_canonical, consequences, impact);
    ann.custom_annotations["GENCODE_BASIC"] = "YES";
    return ann;
}

/**
 * Build annotation with TSL, APPRIS, CCDS, TRANSCRIPT_LENGTH metadata.
 */
static VariantAnnotation make_rich_annotation(
    const std::string& transcript_id,
    const std::string& gene_symbol,
    const std::string& gene_id,
    const std::string& biotype,
    bool is_canonical,
    const std::vector<ConsequenceType>& consequences,
    Impact impact,
    int tsl,
    const std::string& appris,
    bool has_ccds,
    int transcript_length,
    const std::string& mane_select = "")
{
    VariantAnnotation ann = make_annotation(
        transcript_id, gene_symbol, gene_id, biotype, is_canonical, consequences, impact);
    if (tsl > 0) ann.custom_annotations["TSL"] = std::to_string(tsl);
    if (!appris.empty()) ann.custom_annotations["APPRIS"] = appris;
    if (has_ccds) ann.custom_annotations["CCDS"] = "CCDS12345";
    if (transcript_length > 0) ann.custom_annotations["TRANSCRIPT_LENGTH"] = std::to_string(transcript_length);
    if (!mane_select.empty()) ann.custom_annotations["MANE_SELECT"] = mane_select;
    return ann;
}

// ============================================================================
// 1. TranscriptFilter basic filtering
// ============================================================================

// --- canonical_only ---

TEST(TranscriptFilterBasic, CanonicalOnly_PassesCanonical) {
    TranscriptFilterConfig config;
    config.canonical_only = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);

    EXPECT_TRUE(filter.passes_filter(ann));
}

TEST(TranscriptFilterBasic, CanonicalOnly_FiltersNonCanonical) {
    TranscriptFilterConfig config;
    config.canonical_only = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);

    EXPECT_FALSE(filter.passes_filter(ann));
}

TEST(TranscriptFilterBasic, CanonicalOnly_FilterList) {
    TranscriptFilterConfig config;
    config.canonical_only = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW));
    anns.push_back(make_annotation(
        "ENST00000413465", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::INTRON_VARIANT}, Impact::MODIFIER));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].transcript_id, "ENST00000269305");
}

// --- coding_only ---

TEST(TranscriptFilterBasic, CodingOnly_PassesProteinCoding) {
    TranscriptFilterConfig config;
    config.coding_only = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);

    EXPECT_TRUE(filter.passes_filter(ann));
}

TEST(TranscriptFilterBasic, CodingOnly_FiltersNonCoding) {
    TranscriptFilterConfig config;
    config.coding_only = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "ENST00000500000", "LINC01234", "ENSG00000200000", "lncRNA",
        false, {ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT}, Impact::MODIFIER);

    EXPECT_FALSE(filter.passes_filter(ann));
}

TEST(TranscriptFilterBasic, CodingOnly_FiltersMiRNA) {
    TranscriptFilterConfig config;
    config.coding_only = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "ENST00000600000", "MIR1234", "ENSG00000300000", "miRNA",
        false, {ConsequenceType::MATURE_MIRNA_VARIANT}, Impact::MODIFIER);

    EXPECT_FALSE(filter.passes_filter(ann));
}

// --- mane_only ---

TEST(TranscriptFilterBasic, ManeOnly_PassesManeSelect) {
    TranscriptFilterConfig config;
    config.mane_only = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_mane_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);

    EXPECT_TRUE(filter.passes_filter(ann));
}

TEST(TranscriptFilterBasic, ManeOnly_FiltersNonMane) {
    TranscriptFilterConfig config;
    config.mane_only = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);

    EXPECT_FALSE(filter.passes_filter(ann));
}

TEST(TranscriptFilterBasic, ManeOnly_ManePlusClinicalDoesNotCount) {
    // MANE_ONLY only passes MANE_SELECT, not MANE_PLUS_CLINICAL (per Perl VEP)
    TranscriptFilterConfig config;
    config.mane_only = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    ann.custom_annotations["MANE_PLUS_CLINICAL"] = "NM_001276760.4";

    EXPECT_FALSE(filter.passes_filter(ann));
}

// --- gencode_basic ---

TEST(TranscriptFilterBasic, GencodeBasic_PassesGencodeBasic) {
    TranscriptFilterConfig config;
    config.gencode_basic = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_gencode_basic_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);

    EXPECT_TRUE(filter.passes_filter(ann));
}

TEST(TranscriptFilterBasic, GencodeBasic_FiltersNonGencodeBasic) {
    TranscriptFilterConfig config;
    config.gencode_basic = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);

    EXPECT_FALSE(filter.passes_filter(ann));
}

TEST(TranscriptFilterBasic, GencodeBasic_RejectsWrongValue) {
    TranscriptFilterConfig config;
    config.gencode_basic = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    ann.custom_annotations["GENCODE_BASIC"] = "NO";

    EXPECT_FALSE(filter.passes_filter(ann));
}

// --- exclude_predicted ---

TEST(TranscriptFilterBasic, ExcludePredicted_FiltersXM) {
    TranscriptFilterConfig config;
    config.exclude_predicted = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "XM_017024208.2", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);

    EXPECT_FALSE(filter.passes_filter(ann));
}

TEST(TranscriptFilterBasic, ExcludePredicted_FiltersXR) {
    TranscriptFilterConfig config;
    config.exclude_predicted = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "XR_001234567.1", "LOC12345", "ENSG00000200000", "lncRNA",
        false, {ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT}, Impact::MODIFIER);

    EXPECT_FALSE(filter.passes_filter(ann));
}

TEST(TranscriptFilterBasic, ExcludePredicted_PassesNM) {
    TranscriptFilterConfig config;
    config.exclude_predicted = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "NM_000546.6", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);

    EXPECT_TRUE(filter.passes_filter(ann));
}

TEST(TranscriptFilterBasic, ExcludePredicted_PassesENST) {
    TranscriptFilterConfig config;
    config.exclude_predicted = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);

    EXPECT_TRUE(filter.passes_filter(ann));
}

// --- no_intergenic ---

TEST(TranscriptFilterBasic, NoIntergenic_FiltersIntergenicVariants) {
    TranscriptFilterConfig config;
    config.no_intergenic = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "", "", "", "",
        false, {ConsequenceType::INTERGENIC_VARIANT}, Impact::MODIFIER);

    EXPECT_FALSE(filter.passes_filter(ann));
}

TEST(TranscriptFilterBasic, NoIntergenic_PassesNonIntergenicVariants) {
    TranscriptFilterConfig config;
    config.no_intergenic = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);

    EXPECT_TRUE(filter.passes_filter(ann));
}

TEST(TranscriptFilterBasic, NoIntergenic_FilterListRemovesIntergenic) {
    TranscriptFilterConfig config;
    config.no_intergenic = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));
    anns.push_back(make_annotation(
        "", "", "", "",
        false, {ConsequenceType::INTERGENIC_VARIANT}, Impact::MODIFIER));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].transcript_id, "ENST00000269305");
}

// --- biotype filter ---

TEST(TranscriptFilterBasic, BiotypeFilt_OnlyAllowedBiotypes) {
    TranscriptFilterConfig config;
    config.biotypes.insert("protein_coding");
    config.biotypes.insert("lncRNA");
    TranscriptFilter filter(config);

    VariantAnnotation pc = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    VariantAnnotation lnc = make_annotation(
        "ENST00000500000", "LINC01234", "ENSG00000200000", "lncRNA",
        false, {ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT}, Impact::MODIFIER);
    VariantAnnotation mirna = make_annotation(
        "ENST00000600000", "MIR1234", "ENSG00000300000", "miRNA",
        false, {ConsequenceType::MATURE_MIRNA_VARIANT}, Impact::MODIFIER);

    EXPECT_TRUE(filter.passes_filter(pc));
    EXPECT_TRUE(filter.passes_filter(lnc));
    EXPECT_FALSE(filter.passes_filter(mirna));
}

TEST(TranscriptFilterBasic, BiotypeFilt_EmptyMeansAllPass) {
    TranscriptFilterConfig config;
    // biotypes is empty by default
    TranscriptFilter filter(config);

    VariantAnnotation mirna = make_annotation(
        "ENST00000600000", "MIR1234", "ENSG00000300000", "miRNA",
        false, {ConsequenceType::MATURE_MIRNA_VARIANT}, Impact::MODIFIER);

    EXPECT_TRUE(filter.passes_filter(mirna));
}

// --- include_consequences ---

TEST(TranscriptFilterBasic, IncludeConsequences_OnlySpecifiedPass) {
    TranscriptFilterConfig config;
    config.include_consequences.insert("missense_variant");
    config.include_consequences.insert("stop_gained");
    TranscriptFilter filter(config);

    VariantAnnotation missense = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    VariantAnnotation syn = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW);
    VariantAnnotation stop = make_annotation(
        "ENST00000413465", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::STOP_GAINED}, Impact::HIGH);

    EXPECT_TRUE(filter.passes_filter(missense));
    EXPECT_FALSE(filter.passes_filter(syn));
    EXPECT_TRUE(filter.passes_filter(stop));
}

TEST(TranscriptFilterBasic, IncludeConsequences_CaseInsensitive) {
    // consequence_set_contains lowercases the search term and checks the set.
    // So the set must contain the lowercase form for the case-insensitive path to work.
    // When the set has "missense_variant" and the term is "Missense_Variant",
    // the function will lowercase the term to "missense_variant" and find it.
    TranscriptFilterConfig config;
    config.include_consequences.insert("missense_variant");
    TranscriptFilter filter(config);

    VariantAnnotation ann = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);

    // consequence_to_string returns "missense_variant" which matches directly
    EXPECT_TRUE(filter.passes_filter(ann));
}

TEST(TranscriptFilterBasic, IncludeConsequences_MultipleConsequencesOnAnnotation) {
    TranscriptFilterConfig config;
    config.include_consequences.insert("splice_region_variant");
    TranscriptFilter filter(config);

    // Annotation has both intron_variant and splice_region_variant
    VariantAnnotation ann = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::INTRON_VARIANT, ConsequenceType::SPLICE_REGION_VARIANT},
        Impact::LOW);

    // Should pass because splice_region_variant is in the include set
    EXPECT_TRUE(filter.passes_filter(ann));
}

// --- exclude_consequences ---

TEST(TranscriptFilterBasic, ExcludeConsequences_FiltersMostSevere) {
    TranscriptFilterConfig config;
    config.exclude_consequences.insert("intron_variant");
    TranscriptFilter filter(config);

    VariantAnnotation intron = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::INTRON_VARIANT}, Impact::MODIFIER);
    VariantAnnotation missense = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);

    EXPECT_FALSE(filter.passes_filter(intron));
    EXPECT_TRUE(filter.passes_filter(missense));
}

TEST(TranscriptFilterBasic, ExcludeConsequences_OnlyExcludesMostSevere) {
    // If an annotation has multiple consequences, only the most severe
    // is checked against the exclude set
    TranscriptFilterConfig config;
    config.exclude_consequences.insert("intron_variant");
    TranscriptFilter filter(config);

    // Has splice_region_variant (more severe) and intron_variant
    // Most severe is splice_region_variant, which is NOT in exclude set
    VariantAnnotation ann = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::INTRON_VARIANT, ConsequenceType::SPLICE_REGION_VARIANT},
        Impact::LOW);

    EXPECT_TRUE(filter.passes_filter(ann));
}

// --- impact filter ---

TEST(TranscriptFilterBasic, ImpactFilter_OnlySpecifiedImpacts) {
    TranscriptFilterConfig config;
    config.include_impacts.insert(Impact::HIGH);
    config.include_impacts.insert(Impact::MODERATE);
    TranscriptFilter filter(config);

    VariantAnnotation high = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::STOP_GAINED}, Impact::HIGH);
    VariantAnnotation moderate = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    VariantAnnotation low = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW);
    VariantAnnotation modifier = make_annotation(
        "ENST00000413465", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::INTRON_VARIANT}, Impact::MODIFIER);

    EXPECT_TRUE(filter.passes_filter(high));
    EXPECT_TRUE(filter.passes_filter(moderate));
    EXPECT_FALSE(filter.passes_filter(low));
    EXPECT_FALSE(filter.passes_filter(modifier));
}

// --- combined filters ---

TEST(TranscriptFilterBasic, CombinedFilters_CanonicalAndCodingOnly) {
    TranscriptFilterConfig config;
    config.canonical_only = true;
    config.coding_only = true;
    TranscriptFilter filter(config);

    // Canonical + protein_coding -> pass
    VariantAnnotation good = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    // Canonical + lncRNA -> fail
    VariantAnnotation bad1 = make_annotation(
        "ENST00000500000", "LINC01234", "ENSG00000200000", "lncRNA",
        true, {ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT}, Impact::MODIFIER);
    // Non-canonical + protein_coding -> fail
    VariantAnnotation bad2 = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);

    EXPECT_TRUE(filter.passes_filter(good));
    EXPECT_FALSE(filter.passes_filter(bad1));
    EXPECT_FALSE(filter.passes_filter(bad2));
}

// ============================================================================
// 2. Pick modes
// ============================================================================

// --- pick: single best per variant ---

TEST(TranscriptFilterPick, PickReturnsSingleAnnotation) {
    TranscriptFilterConfig config;
    config.pick = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW));
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));
    anns.push_back(make_annotation(
        "ENST00000413465", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::INTRON_VARIANT}, Impact::MODIFIER));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
}

// --- pick_allele: one per allele ---

TEST(TranscriptFilterPick, PickAllele_OnePerAllele) {
    TranscriptFilterConfig config;
    config.pick_allele = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // Allele A>T, two transcripts
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE, "A", "T"));
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW, "A", "T"));
    // Allele A>G, two transcripts
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::STOP_GAINED}, Impact::HIGH, "A", "G"));
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE, "A", "G"));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 2u);

    // Check we have one per allele
    std::set<std::string> alleles;
    for (const auto& r : result) {
        alleles.insert(r.ref_allele + ">" + r.alt_allele);
    }
    EXPECT_EQ(alleles.size(), 2u);
    EXPECT_TRUE(alleles.count("A>T"));
    EXPECT_TRUE(alleles.count("A>G"));
}

// --- per_gene: one per gene ---

TEST(TranscriptFilterPick, PerGene_OnePerGene) {
    TranscriptFilterConfig config;
    config.per_gene = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // Gene TP53, two transcripts
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW));
    // Gene BRCA1, one transcript
    anns.push_back(make_annotation(
        "ENST00000357654", "BRCA1", "ENSG00000012048", "protein_coding",
        true, {ConsequenceType::INTRON_VARIANT}, Impact::MODIFIER));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 2u);

    // Check we have one per gene_id
    std::set<std::string> genes;
    for (const auto& r : result) {
        genes.insert(r.gene_id);
    }
    EXPECT_EQ(genes.size(), 2u);
    EXPECT_TRUE(genes.count("ENSG00000141510"));
    EXPECT_TRUE(genes.count("ENSG00000012048"));
}

TEST(TranscriptFilterPick, PerGene_UsesGeneIdNotSymbol) {
    TranscriptFilterConfig config;
    config.per_gene = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // Two transcripts, same gene_symbol but different gene_id
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000999999", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW));

    auto result = filter.filter(anns);
    // Should be 2 because gene_ids differ
    ASSERT_EQ(result.size(), 2u);
}

// --- pick_allele_gene: one per allele+gene ---

TEST(TranscriptFilterPick, PickAlleleGene_OnePerAlleleGene) {
    TranscriptFilterConfig config;
    config.pick_allele_gene = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // Gene TP53, allele A>T, two transcripts
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE, "A", "T"));
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW, "A", "T"));
    // Gene TP53, allele A>G, one transcript
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::STOP_GAINED}, Impact::HIGH, "A", "G"));
    // Gene BRCA1, allele A>T, one transcript
    anns.push_back(make_annotation(
        "ENST00000357654", "BRCA1", "ENSG00000012048", "protein_coding",
        true, {ConsequenceType::INTRON_VARIANT}, Impact::MODIFIER, "A", "T"));

    auto result = filter.filter(anns);
    // 3 groups: TP53+A>T, TP53+A>G, BRCA1+A>T
    ASSERT_EQ(result.size(), 3u);
}

// --- flag_pick: flags best but keeps all ---

TEST(TranscriptFilterPick, FlagPick_KeepsAllAnnotations) {
    TranscriptFilterConfig config;
    config.flag_pick = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW));
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));
    anns.push_back(make_annotation(
        "ENST00000413465", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::INTRON_VARIANT}, Impact::MODIFIER));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 3u);

    // Exactly one should be flagged with PICK=1
    int pick_count = 0;
    for (const auto& r : result) {
        auto it = r.custom_annotations.find("PICK");
        if (it != r.custom_annotations.end() && it->second == "1") {
            ++pick_count;
        }
    }
    EXPECT_EQ(pick_count, 1);
}

TEST(TranscriptFilterPick, FlagPick_PickedHasPickAnnotation) {
    TranscriptFilterConfig config;
    config.flag_pick = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::INTRON_VARIANT}, Impact::MODIFIER));
    anns.push_back(make_mane_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 2u);

    // The one with PICK=1 should exist
    bool found_pick = false;
    for (const auto& r : result) {
        auto it = r.custom_annotations.find("PICK");
        if (it != r.custom_annotations.end() && it->second == "1") {
            found_pick = true;
        }
    }
    EXPECT_TRUE(found_pick);
}

// --- most_severe: only most severe consequence ---

TEST(TranscriptFilterPick, MostSevere_ReturnsMostSevere) {
    TranscriptFilterConfig config;
    config.most_severe = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW));
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::STOP_GAINED}, Impact::HIGH));
    anns.push_back(make_annotation(
        "ENST00000413465", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].consequences.size(), 1u);
    EXPECT_EQ(result[0].consequences[0], ConsequenceType::STOP_GAINED);
}

TEST(TranscriptFilterPick, MostSevere_TieBreaksToSingle) {
    TranscriptFilterConfig config;
    config.most_severe = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // Two annotations with the same most severe consequence
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));

    auto result = filter.filter(anns);
    // Should pick one out of the tied pair
    ASSERT_EQ(result.size(), 1u);
}

// ============================================================================
// 3. Pick ordering (pick_order criteria)
// ============================================================================

TEST(TranscriptFilterPickOrder, CanonicalPreferred) {
    TranscriptFilterConfig config;
    config.pick = true;
    // Use CANONICAL as the first (most important) criteria
    config.pick_order = {PickCriteria::CANONICAL, PickCriteria::RANK};
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::STOP_GAINED}, Impact::HIGH));
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].transcript_id, "ENST00000269305");
    EXPECT_TRUE(result[0].is_canonical);
}

TEST(TranscriptFilterPickOrder, ManePreferred) {
    TranscriptFilterConfig config;
    config.pick = true;
    config.pick_order = {PickCriteria::MANE_SELECT, PickCriteria::RANK};
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::STOP_GAINED}, Impact::HIGH));
    anns.push_back(make_mane_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].transcript_id, "ENST00000269305");
}

TEST(TranscriptFilterPickOrder, BiotypeRanking_ProteinCodingPreferred) {
    TranscriptFilterConfig config;
    config.pick = true;
    config.pick_order = {PickCriteria::BIOTYPE, PickCriteria::RANK};
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    anns.push_back(make_annotation(
        "ENST00000500000", "LINC01234", "ENSG00000200000", "lncRNA",
        false, {ConsequenceType::STOP_GAINED}, Impact::HIGH));
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].biotype, "protein_coding");
}

TEST(TranscriptFilterPickOrder, BiotypeRanking_LncRNAOverMiRNA) {
    TranscriptFilterConfig config;
    config.pick = true;
    config.pick_order = {PickCriteria::BIOTYPE};
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    anns.push_back(make_annotation(
        "ENST00000600000", "MIR1234", "ENSG00000300000", "miRNA",
        false, {ConsequenceType::MATURE_MIRNA_VARIANT}, Impact::MODIFIER));
    anns.push_back(make_annotation(
        "ENST00000500000", "LINC01234", "ENSG00000200000", "lncRNA",
        false, {ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT}, Impact::MODIFIER));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].biotype, "lncRNA");
}

TEST(TranscriptFilterPickOrder, ApprisRanking) {
    TranscriptFilterConfig config;
    config.pick = true;
    config.pick_order = {PickCriteria::APPRIS};
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    VariantAnnotation alt2 = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    alt2.custom_annotations["APPRIS"] = "alternative_2";
    anns.push_back(alt2);

    VariantAnnotation p1 = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    p1.custom_annotations["APPRIS"] = "principal_1";
    anns.push_back(p1);

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].transcript_id, "ENST00000269305");
}

TEST(TranscriptFilterPickOrder, ApprisRanking_PrincipalBeforeAlternative) {
    TranscriptFilterConfig config;
    config.pick = true;
    config.pick_order = {PickCriteria::APPRIS};
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    VariantAnnotation alt1 = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    alt1.custom_annotations["APPRIS"] = "alternative_1";
    anns.push_back(alt1);

    VariantAnnotation p5 = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    p5.custom_annotations["APPRIS"] = "principal_5";
    anns.push_back(p5);

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    // principal_5 (score 5) < alternative_1 (score 6)
    EXPECT_EQ(result[0].transcript_id, "ENST00000269305");
}

TEST(TranscriptFilterPickOrder, TslRanking) {
    TranscriptFilterConfig config;
    config.pick = true;
    config.pick_order = {PickCriteria::TSL};
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    VariantAnnotation tsl5 = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    tsl5.custom_annotations["TSL"] = "5";
    anns.push_back(tsl5);

    VariantAnnotation tsl1 = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    tsl1.custom_annotations["TSL"] = "1";
    anns.push_back(tsl1);

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].transcript_id, "ENST00000269305");
}

TEST(TranscriptFilterPickOrder, TslRanking_UnknownTslWorst) {
    TranscriptFilterConfig config;
    config.pick = true;
    config.pick_order = {PickCriteria::TSL};
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // No TSL set (treated as 0 -> score 10, worst)
    VariantAnnotation no_tsl = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    anns.push_back(no_tsl);

    VariantAnnotation tsl3 = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    tsl3.custom_annotations["TSL"] = "3";
    anns.push_back(tsl3);

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].transcript_id, "ENST00000269305");
}

TEST(TranscriptFilterPickOrder, RankBySeverity) {
    TranscriptFilterConfig config;
    config.pick = true;
    config.pick_order = {PickCriteria::RANK};
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW));
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::STOP_GAINED}, Impact::HIGH));
    anns.push_back(make_annotation(
        "ENST00000413465", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].consequences[0], ConsequenceType::STOP_GAINED);
}

TEST(TranscriptFilterPickOrder, LengthPreferred) {
    TranscriptFilterConfig config;
    config.pick = true;
    config.pick_order = {PickCriteria::LENGTH};
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    VariantAnnotation short_tx = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    short_tx.custom_annotations["TRANSCRIPT_LENGTH"] = "1000";
    anns.push_back(short_tx);

    VariantAnnotation long_tx = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    long_tx.custom_annotations["TRANSCRIPT_LENGTH"] = "5000";
    anns.push_back(long_tx);

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    // Longer transcript preferred
    EXPECT_EQ(result[0].transcript_id, "ENST00000269305");
}

TEST(TranscriptFilterPickOrder, CcdsPreferred) {
    TranscriptFilterConfig config;
    config.pick = true;
    config.pick_order = {PickCriteria::CCDS};
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    VariantAnnotation no_ccds = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    anns.push_back(no_ccds);

    VariantAnnotation with_ccds = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    with_ccds.custom_annotations["CCDS"] = "CCDS11118.2";
    anns.push_back(with_ccds);

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].transcript_id, "ENST00000269305");
}

TEST(TranscriptFilterPickOrder, DefaultPickOrder_ManeBeforeCanonical) {
    // Default pick order has MANE_SELECT before CANONICAL
    TranscriptFilterConfig config;
    config.pick = true;
    // Use default pick_order
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // Canonical but NOT MANE
    VariantAnnotation canonical_only = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    anns.push_back(canonical_only);

    // MANE but NOT canonical
    VariantAnnotation mane_only = make_mane_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    anns.push_back(mane_only);

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    // MANE_SELECT should be preferred over CANONICAL in default order
    EXPECT_EQ(result[0].transcript_id, "ENST00000269305");
}

TEST(TranscriptFilterPickOrder, EnsemblPreferred) {
    TranscriptFilterConfig config;
    config.pick = true;
    config.pick_order = {PickCriteria::ENSEMBL};
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    VariantAnnotation refseq_ann = make_annotation(
        "NM_000546.6", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    refseq_ann.source = "RefSeq";
    anns.push_back(refseq_ann);

    VariantAnnotation ensembl_ann = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    ensembl_ann.source = "Ensembl";
    anns.push_back(ensembl_ann);

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].source, "Ensembl");
}

TEST(TranscriptFilterPickOrder, RefSeqPreferred) {
    TranscriptFilterConfig config;
    config.pick = true;
    config.pick_order = {PickCriteria::REFSEQ};
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    VariantAnnotation ensembl_ann = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    ensembl_ann.source = "Ensembl";
    anns.push_back(ensembl_ann);

    VariantAnnotation refseq_ann = make_annotation(
        "NM_000546.6", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    refseq_ann.source = "RefSeq";
    anns.push_back(refseq_ann);

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].source, "RefSeq");
}

TEST(TranscriptFilterPickOrder, MultipleCriteriaTieBreak) {
    // When first criteria is tied, fall through to next criteria
    TranscriptFilterConfig config;
    config.pick = true;
    config.pick_order = {PickCriteria::CANONICAL, PickCriteria::RANK};
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // Both canonical, but different consequences
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW));
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::STOP_GAINED}, Impact::HIGH));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    // Both canonical (tied), so RANK (stop_gained more severe) wins
    EXPECT_EQ(result[0].consequences[0], ConsequenceType::STOP_GAINED);
}

// ============================================================================
// Rich annotation pick ordering test
// ============================================================================

TEST(TranscriptFilterPickOrder, ComplexPickWithAllCriteria) {
    TranscriptFilterConfig config;
    config.pick = true;
    // Use a complex pick order
    config.pick_order = {
        PickCriteria::MANE_SELECT,
        PickCriteria::CANONICAL,
        PickCriteria::APPRIS,
        PickCriteria::TSL,
        PickCriteria::BIOTYPE,
        PickCriteria::CCDS,
        PickCriteria::RANK,
        PickCriteria::LENGTH
    };
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;

    // Not MANE, not canonical, alternative_1, TSL 3
    anns.push_back(make_rich_annotation(
        "ENST00000413465", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::STOP_GAINED}, Impact::HIGH,
        3, "alternative_1", false, 2000));

    // Not MANE, canonical, principal_1, TSL 1
    anns.push_back(make_rich_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE,
        1, "principal_1", true, 5000));

    // MANE, not canonical, principal_2, TSL 1
    anns.push_back(make_rich_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE,
        1, "principal_2", true, 4000, "NM_000546.6"));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    // MANE_SELECT is first criterion -> ENST00000359597 wins
    EXPECT_EQ(result[0].transcript_id, "ENST00000359597");
}

// ============================================================================
// 4. Frequency filtering
// ============================================================================

TEST(TranscriptFilterFrequency, FreqLessThan_PassesRare) {
    TranscriptFilterConfig config;
    config.check_frequency = true;
    config.freq_pop = "MAX_AF";
    config.freq_threshold = 0.01;
    config.freq_gt = false; // keep variants with freq < threshold (rare)
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    VariantAnnotation rare = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    rare.custom_annotations["MAX_AF"] = "0.0001";
    anns.push_back(rare);

    VariantAnnotation common = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW);
    common.custom_annotations["MAX_AF"] = "0.15";
    anns.push_back(common);

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].transcript_id, "ENST00000269305");
}

TEST(TranscriptFilterFrequency, FreqGreaterThan_PassesCommon) {
    TranscriptFilterConfig config;
    config.check_frequency = true;
    config.freq_pop = "gnomAD_AF";
    config.freq_threshold = 0.05;
    config.freq_gt = true; // keep variants with freq > threshold (common)
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    VariantAnnotation rare = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    rare.custom_annotations["gnomAD_AF"] = "0.001";
    anns.push_back(rare);

    VariantAnnotation common = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW);
    common.custom_annotations["gnomAD_AF"] = "0.25";
    anns.push_back(common);

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].transcript_id, "ENST00000359597");
}

TEST(TranscriptFilterFrequency, FreqMissing_KeptByDefault) {
    // Annotations without the frequency field should NOT be removed
    TranscriptFilterConfig config;
    config.check_frequency = true;
    config.freq_pop = "MAX_AF";
    config.freq_threshold = 0.01;
    config.freq_gt = false;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    VariantAnnotation no_freq = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    // No MAX_AF set
    anns.push_back(no_freq);

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
}

TEST(TranscriptFilterFrequency, FreqAtThreshold_LessThanMode) {
    TranscriptFilterConfig config;
    config.check_frequency = true;
    config.freq_pop = "MAX_AF";
    config.freq_threshold = 0.01;
    config.freq_gt = false;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    VariantAnnotation at_threshold = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    at_threshold.custom_annotations["MAX_AF"] = "0.01";
    anns.push_back(at_threshold);

    auto result = filter.filter(anns);
    // freq (0.01) is NOT less than threshold (0.01), so it should be filtered
    ASSERT_EQ(result.size(), 0u);
}

TEST(TranscriptFilterFrequency, FreqAtThreshold_GreaterThanMode) {
    TranscriptFilterConfig config;
    config.check_frequency = true;
    config.freq_pop = "MAX_AF";
    config.freq_threshold = 0.01;
    config.freq_gt = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    VariantAnnotation at_threshold = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    at_threshold.custom_annotations["MAX_AF"] = "0.01";
    anns.push_back(at_threshold);

    auto result = filter.filter(anns);
    // freq (0.01) is NOT greater than threshold (0.01), so it should be filtered
    ASSERT_EQ(result.size(), 0u);
}

TEST(TranscriptFilterFrequency, FreqInvalidValue_NotFiltered) {
    TranscriptFilterConfig config;
    config.check_frequency = true;
    config.freq_pop = "MAX_AF";
    config.freq_threshold = 0.01;
    config.freq_gt = false;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    VariantAnnotation bad_freq = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    bad_freq.custom_annotations["MAX_AF"] = "not_a_number";
    anns.push_back(bad_freq);

    auto result = filter.filter(anns);
    // Invalid frequency value -> stod throws, caught by catch -> not filtered
    ASSERT_EQ(result.size(), 1u);
}

// ============================================================================
// 5. Edge cases
// ============================================================================

TEST(TranscriptFilterEdge, EmptyAnnotationsList) {
    TranscriptFilterConfig config;
    config.pick = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    auto result = filter.filter(anns);
    EXPECT_TRUE(result.empty());
}

TEST(TranscriptFilterEdge, SingleAnnotation_Pick) {
    TranscriptFilterConfig config;
    config.pick = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].transcript_id, "ENST00000269305");
}

TEST(TranscriptFilterEdge, SingleAnnotation_FlagPick) {
    TranscriptFilterConfig config;
    config.flag_pick = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
    auto it = result[0].custom_annotations.find("PICK");
    ASSERT_NE(it, result[0].custom_annotations.end());
    EXPECT_EQ(it->second, "1");
}

TEST(TranscriptFilterEdge, MultipleAnnotationsSameGene_PerGene) {
    TranscriptFilterConfig config;
    config.per_gene = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    for (int i = 0; i < 5; ++i) {
        anns.push_back(make_annotation(
            "ENST0000000000" + std::to_string(i), "TP53", "ENSG00000141510",
            "protein_coding", i == 0,
            {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));
    }

    auto result = filter.filter(anns);
    // All same gene -> only one should be returned
    ASSERT_EQ(result.size(), 1u);
}

TEST(TranscriptFilterEdge, AnnotationWithNoConsequences) {
    TranscriptFilterConfig config;
    config.pick = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // One with no consequences
    VariantAnnotation empty_csq;
    empty_csq.transcript_id = "ENST00000269305";
    empty_csq.gene_symbol = "TP53";
    empty_csq.gene_id = "ENSG00000141510";
    empty_csq.biotype = "protein_coding";
    empty_csq.is_canonical = true;
    empty_csq.ref_allele = "A";
    empty_csq.alt_allele = "T";
    // consequences is empty
    anns.push_back(empty_csq);

    // One with consequences
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));

    auto result = filter.filter(anns);
    ASSERT_EQ(result.size(), 1u);
}

TEST(TranscriptFilterEdge, AnnotationWithEmptyTranscriptId) {
    TranscriptFilterConfig config;
    TranscriptFilter filter(config);

    VariantAnnotation ann;
    ann.transcript_id = "";
    ann.biotype = "protein_coding";
    ann.is_canonical = false;
    ann.consequences = {ConsequenceType::INTERGENIC_VARIANT};
    ann.impact = Impact::MODIFIER;
    ann.ref_allele = "A";
    ann.alt_allele = "T";

    // Should still pass default filter (no constraints set)
    EXPECT_TRUE(filter.passes_filter(ann));
}

TEST(TranscriptFilterEdge, ExcludePredicted_EmptyTranscriptId) {
    TranscriptFilterConfig config;
    config.exclude_predicted = true;
    TranscriptFilter filter(config);

    VariantAnnotation ann;
    ann.transcript_id = "";
    ann.biotype = "";
    ann.consequences = {ConsequenceType::INTERGENIC_VARIANT};
    ann.impact = Impact::MODIFIER;
    ann.ref_allele = "A";
    ann.alt_allele = "T";

    // Empty transcript_id should not cause issues; exclude_predicted checks prefix
    // but the code guards with !ann.transcript_id.empty()
    EXPECT_TRUE(filter.passes_filter(ann));
}

TEST(TranscriptFilterEdge, PerGene_UsesGeneSymbolFallbackWhenIdEmpty) {
    TranscriptFilterConfig config;
    config.per_gene = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // gene_id empty, gene_symbol used as fallback
    VariantAnnotation ann1 = make_annotation(
        "ENST00000269305", "TP53", "", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    VariantAnnotation ann2 = make_annotation(
        "ENST00000359597", "TP53", "", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW);

    anns.push_back(ann1);
    anns.push_back(ann2);

    auto result = filter.filter(anns);
    // Both have same gene_symbol "TP53" and empty gene_id -> grouped together
    ASSERT_EQ(result.size(), 1u);
}

TEST(TranscriptFilterEdge, NoFilterApplied_AllPass) {
    TranscriptFilterConfig config;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));
    anns.push_back(make_annotation(
        "ENST00000500000", "LINC01234", "ENSG00000200000", "lncRNA",
        false, {ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT}, Impact::MODIFIER));
    anns.push_back(make_annotation(
        "", "", "", "",
        false, {ConsequenceType::INTERGENIC_VARIANT}, Impact::MODIFIER));

    auto result = filter.filter(anns);
    EXPECT_EQ(result.size(), 3u);
}

TEST(TranscriptFilterEdge, AllFilteredOut) {
    TranscriptFilterConfig config;
    config.canonical_only = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // None are canonical
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));
    anns.push_back(make_annotation(
        "ENST00000413465", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW));

    auto result = filter.filter(anns);
    EXPECT_TRUE(result.empty());
}

TEST(TranscriptFilterEdge, FlagPickAllele_FlagsOnePerAllele) {
    TranscriptFilterConfig config;
    config.flag_pick_allele = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // Allele A>T, two transcripts
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE, "A", "T"));
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW, "A", "T"));
    // Allele A>G, one transcript
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::STOP_GAINED}, Impact::HIGH, "A", "G"));

    auto result = filter.filter(anns);
    // All 3 kept
    ASSERT_EQ(result.size(), 3u);

    // Count PICK flags: should be 2 (one per allele)
    int pick_count = 0;
    for (const auto& r : result) {
        auto it = r.custom_annotations.find("PICK");
        if (it != r.custom_annotations.end() && it->second == "1") {
            ++pick_count;
        }
    }
    EXPECT_EQ(pick_count, 2);
}

TEST(TranscriptFilterEdge, FlagPickAlleleGene_FlagsOnePerAlleleGene) {
    TranscriptFilterConfig config;
    config.flag_pick_allele_gene = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // Gene TP53, allele A>T, two transcripts
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE, "A", "T"));
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW, "A", "T"));
    // Gene BRCA1, allele A>T, one transcript
    anns.push_back(make_annotation(
        "ENST00000357654", "BRCA1", "ENSG00000012048", "protein_coding",
        true, {ConsequenceType::INTRON_VARIANT}, Impact::MODIFIER, "A", "T"));

    auto result = filter.filter(anns);
    // All 3 kept
    ASSERT_EQ(result.size(), 3u);

    // Count PICK flags: should be 2 (TP53:A>T and BRCA1:A>T)
    int pick_count = 0;
    for (const auto& r : result) {
        auto it = r.custom_annotations.find("PICK");
        if (it != r.custom_annotations.end() && it->second == "1") {
            ++pick_count;
        }
    }
    EXPECT_EQ(pick_count, 2);
}

// ============================================================================
// 6. parse_pick_order and parse_biotypes utilities
// ============================================================================

TEST(TranscriptFilterParse, ParsePickOrder_BasicTokens) {
    auto order = parse_pick_order("canonical,mane,appris,tsl,biotype,ccds,rank,length");
    ASSERT_EQ(order.size(), 8u);
    EXPECT_EQ(order[0], PickCriteria::CANONICAL);
    EXPECT_EQ(order[1], PickCriteria::MANE_SELECT);
    EXPECT_EQ(order[2], PickCriteria::APPRIS);
    EXPECT_EQ(order[3], PickCriteria::TSL);
    EXPECT_EQ(order[4], PickCriteria::BIOTYPE);
    EXPECT_EQ(order[5], PickCriteria::CCDS);
    EXPECT_EQ(order[6], PickCriteria::RANK);
    EXPECT_EQ(order[7], PickCriteria::LENGTH);
}

TEST(TranscriptFilterParse, ParsePickOrder_CaseInsensitive) {
    auto order = parse_pick_order("CANONICAL,MANE_SELECT,APPRIS");
    ASSERT_EQ(order.size(), 3u);
    EXPECT_EQ(order[0], PickCriteria::CANONICAL);
    EXPECT_EQ(order[1], PickCriteria::MANE_SELECT);
    EXPECT_EQ(order[2], PickCriteria::APPRIS);
}

TEST(TranscriptFilterParse, ParsePickOrder_WithWhitespace) {
    auto order = parse_pick_order(" canonical , mane , rank ");
    ASSERT_EQ(order.size(), 3u);
    EXPECT_EQ(order[0], PickCriteria::CANONICAL);
    EXPECT_EQ(order[1], PickCriteria::MANE_SELECT);
    EXPECT_EQ(order[2], PickCriteria::RANK);
}

TEST(TranscriptFilterParse, ParsePickOrder_ManeAliases) {
    auto order = parse_pick_order("mane,mane_select,mane_plus,mane_plus_clinical");
    ASSERT_EQ(order.size(), 4u);
    EXPECT_EQ(order[0], PickCriteria::MANE_SELECT);
    EXPECT_EQ(order[1], PickCriteria::MANE_SELECT);
    EXPECT_EQ(order[2], PickCriteria::MANE_PLUS);
    EXPECT_EQ(order[3], PickCriteria::MANE_PLUS);
}

TEST(TranscriptFilterParse, ParsePickOrder_EnsemblRefseq) {
    auto order = parse_pick_order("ensembl,refseq");
    ASSERT_EQ(order.size(), 2u);
    EXPECT_EQ(order[0], PickCriteria::ENSEMBL);
    EXPECT_EQ(order[1], PickCriteria::REFSEQ);
}

TEST(TranscriptFilterParse, ParsePickOrder_UnknownTokensIgnored) {
    auto order = parse_pick_order("canonical,foobar,rank");
    ASSERT_EQ(order.size(), 2u);
    EXPECT_EQ(order[0], PickCriteria::CANONICAL);
    EXPECT_EQ(order[1], PickCriteria::RANK);
}

TEST(TranscriptFilterParse, ParsePickOrder_Empty) {
    auto order = parse_pick_order("");
    EXPECT_TRUE(order.empty());
}

TEST(TranscriptFilterParse, ParseBiotypes_Basic) {
    auto biotypes = parse_biotypes("protein_coding,lncRNA,miRNA");
    ASSERT_EQ(biotypes.size(), 3u);
    EXPECT_TRUE(biotypes.count("protein_coding"));
    EXPECT_TRUE(biotypes.count("lncRNA"));
    EXPECT_TRUE(biotypes.count("miRNA"));
}

TEST(TranscriptFilterParse, ParseBiotypes_WithWhitespace) {
    auto biotypes = parse_biotypes(" protein_coding , lncRNA ");
    ASSERT_EQ(biotypes.size(), 2u);
    EXPECT_TRUE(biotypes.count("protein_coding"));
    EXPECT_TRUE(biotypes.count("lncRNA"));
}

TEST(TranscriptFilterParse, ParseBiotypes_Empty) {
    auto biotypes = parse_biotypes("");
    EXPECT_TRUE(biotypes.empty());
}

// ============================================================================
// 7. Consequence rank utility
// ============================================================================

TEST(ConsequenceRank, SeverityOrdering) {
    // transcript_ablation should be most severe (lowest rank)
    EXPECT_LT(get_consequence_rank(ConsequenceType::TRANSCRIPT_ABLATION),
              get_consequence_rank(ConsequenceType::SPLICE_ACCEPTOR_VARIANT));
    EXPECT_LT(get_consequence_rank(ConsequenceType::SPLICE_ACCEPTOR_VARIANT),
              get_consequence_rank(ConsequenceType::STOP_GAINED));
    EXPECT_LT(get_consequence_rank(ConsequenceType::STOP_GAINED),
              get_consequence_rank(ConsequenceType::FRAMESHIFT_VARIANT));
    EXPECT_LT(get_consequence_rank(ConsequenceType::FRAMESHIFT_VARIANT),
              get_consequence_rank(ConsequenceType::MISSENSE_VARIANT));
    EXPECT_LT(get_consequence_rank(ConsequenceType::MISSENSE_VARIANT),
              get_consequence_rank(ConsequenceType::SYNONYMOUS_VARIANT));
    EXPECT_LT(get_consequence_rank(ConsequenceType::SYNONYMOUS_VARIANT),
              get_consequence_rank(ConsequenceType::INTRON_VARIANT));
    EXPECT_LT(get_consequence_rank(ConsequenceType::INTRON_VARIANT),
              get_consequence_rank(ConsequenceType::INTERGENIC_VARIANT));
}

TEST(ConsequenceRank, UnknownIsWorst) {
    EXPECT_GT(get_consequence_rank(ConsequenceType::UNKNOWN),
              get_consequence_rank(ConsequenceType::INTERGENIC_VARIANT));
}

// ============================================================================
// 8. APPRIS comparison
// ============================================================================

TEST(ApprisCompare, Principal1IsBest) {
    EXPECT_LT(compare_appris("principal_1", "principal_2"), 0);
    EXPECT_LT(compare_appris("principal_1", "alternative_1"), 0);
    EXPECT_LT(compare_appris("principal_1", ""), 0);
}

TEST(ApprisCompare, EmptyIsWorst) {
    EXPECT_GT(compare_appris("", "principal_1"), 0);
    EXPECT_GT(compare_appris("", "alternative_5"), 0);
}

TEST(ApprisCompare, EqualValues) {
    EXPECT_EQ(compare_appris("principal_1", "principal_1"), 0);
    EXPECT_EQ(compare_appris("alternative_2", "alternative_2"), 0);
    EXPECT_EQ(compare_appris("", ""), 0);
}

TEST(ApprisCompare, PrincipalBeforeAlternative) {
    EXPECT_LT(compare_appris("principal_5", "alternative_1"), 0);
}

// ============================================================================
// 9. AnnotationWithMeta::get_rank()
// ============================================================================

TEST(AnnotationWithMeta, GetRank_EmptyConsequences) {
    AnnotationWithMeta ext;
    ext.annotation.consequences = {};
    EXPECT_EQ(ext.get_rank(), 100);
}

TEST(AnnotationWithMeta, GetRank_SingleConsequence) {
    AnnotationWithMeta ext;
    ext.annotation.consequences = {ConsequenceType::MISSENSE_VARIANT};
    EXPECT_EQ(ext.get_rank(), get_consequence_rank(ConsequenceType::MISSENSE_VARIANT));
}

TEST(AnnotationWithMeta, GetRank_MultipleConsequences_ReturnsMostSevere) {
    AnnotationWithMeta ext;
    ext.annotation.consequences = {
        ConsequenceType::INTRON_VARIANT,
        ConsequenceType::SPLICE_REGION_VARIANT
    };
    // splice_region_variant is more severe than intron_variant
    EXPECT_EQ(ext.get_rank(), get_consequence_rank(ConsequenceType::SPLICE_REGION_VARIANT));
}

// ============================================================================
// 10. consequence_set_contains (case-insensitive)
// ============================================================================

TEST(ConsequenceSetContains, ExactMatch) {
    std::set<std::string> s = {"missense_variant", "stop_gained"};
    EXPECT_TRUE(consequence_set_contains(s, "missense_variant"));
    EXPECT_TRUE(consequence_set_contains(s, "stop_gained"));
    EXPECT_FALSE(consequence_set_contains(s, "synonymous_variant"));
}

TEST(ConsequenceSetContains, CaseInsensitiveMatch) {
    std::set<std::string> s = {"missense_variant"};
    EXPECT_TRUE(consequence_set_contains(s, "Missense_Variant"));
    EXPECT_TRUE(consequence_set_contains(s, "MISSENSE_VARIANT"));
    EXPECT_TRUE(consequence_set_contains(s, "missense_variant"));
}

TEST(ConsequenceSetContains, EmptySet) {
    std::set<std::string> s;
    EXPECT_FALSE(consequence_set_contains(s, "missense_variant"));
}

// ============================================================================
// 11. TranscriptFilterConfig defaults
// ============================================================================

TEST(TranscriptFilterConfig, DefaultPickOrderMatchesPerl) {
    TranscriptFilterConfig config;
    // Default order: MANE_SELECT, MANE_PLUS, CANONICAL, APPRIS, TSL, BIOTYPE, CCDS, RANK, LENGTH
    ASSERT_EQ(config.pick_order.size(), 9u);
    EXPECT_EQ(config.pick_order[0], PickCriteria::MANE_SELECT);
    EXPECT_EQ(config.pick_order[1], PickCriteria::MANE_PLUS);
    EXPECT_EQ(config.pick_order[2], PickCriteria::CANONICAL);
    EXPECT_EQ(config.pick_order[3], PickCriteria::APPRIS);
    EXPECT_EQ(config.pick_order[4], PickCriteria::TSL);
    EXPECT_EQ(config.pick_order[5], PickCriteria::BIOTYPE);
    EXPECT_EQ(config.pick_order[6], PickCriteria::CCDS);
    EXPECT_EQ(config.pick_order[7], PickCriteria::RANK);
    EXPECT_EQ(config.pick_order[8], PickCriteria::LENGTH);
}

TEST(TranscriptFilterConfig, DefaultFlagsFalse) {
    TranscriptFilterConfig config;
    EXPECT_FALSE(config.pick);
    EXPECT_FALSE(config.pick_allele);
    EXPECT_FALSE(config.pick_allele_gene);
    EXPECT_FALSE(config.per_gene);
    EXPECT_FALSE(config.most_severe);
    EXPECT_FALSE(config.flag_pick);
    EXPECT_FALSE(config.flag_pick_allele);
    EXPECT_FALSE(config.flag_pick_allele_gene);
    EXPECT_FALSE(config.canonical_only);
    EXPECT_FALSE(config.mane_only);
    EXPECT_FALSE(config.coding_only);
    EXPECT_FALSE(config.gencode_basic);
    EXPECT_FALSE(config.exclude_predicted);
    EXPECT_FALSE(config.no_intergenic);
    EXPECT_FALSE(config.check_frequency);
    EXPECT_FALSE(config.freq_gt);
}

TEST(TranscriptFilterConfig, DefaultFreqThreshold) {
    TranscriptFilterConfig config;
    EXPECT_DOUBLE_EQ(config.freq_threshold, 0.01);
}

// ============================================================================
// 12. set_config and config() accessors
// ============================================================================

TEST(TranscriptFilterAccessors, SetConfigUpdatesConfig) {
    TranscriptFilter filter;

    TranscriptFilterConfig new_config;
    new_config.canonical_only = true;
    new_config.coding_only = true;
    new_config.pick = true;
    filter.set_config(new_config);

    EXPECT_TRUE(filter.config().canonical_only);
    EXPECT_TRUE(filter.config().coding_only);
    EXPECT_TRUE(filter.config().pick);
}

// ============================================================================
// 13. Filter + Pick combined scenarios
// ============================================================================

TEST(TranscriptFilterCombined, CanonicalFilterThenPick) {
    TranscriptFilterConfig config;
    config.canonical_only = true;
    config.pick = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // Non-canonical with severe consequence
    anns.push_back(make_annotation(
        "ENST00000413465", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::STOP_GAINED}, Impact::HIGH));
    // Canonical with less severe consequence
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));
    // Another canonical
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW));

    auto result = filter.filter(anns);
    // First filter to canonical only (2 remain), then pick one
    ASSERT_EQ(result.size(), 1u);
    EXPECT_TRUE(result[0].is_canonical);
}

TEST(TranscriptFilterCombined, CodingOnlyThenPerGene) {
    TranscriptFilterConfig config;
    config.coding_only = true;
    config.per_gene = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // protein_coding TP53
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE));
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW));
    // lncRNA (should be filtered by coding_only)
    anns.push_back(make_annotation(
        "ENST00000500000", "LINC01234", "ENSG00000200000", "lncRNA",
        false, {ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT}, Impact::MODIFIER));

    auto result = filter.filter(anns);
    // lncRNA filtered out, then per_gene picks one for TP53
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].gene_id, "ENSG00000141510");
}

TEST(TranscriptFilterCombined, FrequencyFilterThenMostSevere) {
    TranscriptFilterConfig config;
    config.check_frequency = true;
    config.freq_pop = "MAX_AF";
    config.freq_threshold = 0.01;
    config.freq_gt = false;
    config.most_severe = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // Rare, stop_gained
    VariantAnnotation rare_severe = make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::STOP_GAINED}, Impact::HIGH);
    rare_severe.custom_annotations["MAX_AF"] = "0.0001";
    anns.push_back(rare_severe);

    // Common, missense
    VariantAnnotation common_moderate = make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE);
    common_moderate.custom_annotations["MAX_AF"] = "0.15";
    anns.push_back(common_moderate);

    // Rare, synonymous
    VariantAnnotation rare_low = make_annotation(
        "ENST00000413465", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW);
    rare_low.custom_annotations["MAX_AF"] = "0.005";
    anns.push_back(rare_low);

    auto result = filter.filter(anns);
    // Common one filtered out, then most_severe picks stop_gained
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].consequences[0], ConsequenceType::STOP_GAINED);
}

TEST(TranscriptFilterCombined, BiotypeFiltThenPickAlleleGene) {
    TranscriptFilterConfig config;
    config.biotypes.insert("protein_coding");
    config.pick_allele_gene = true;
    TranscriptFilter filter(config);

    std::vector<VariantAnnotation> anns;
    // protein_coding TP53 A>T - two transcripts
    anns.push_back(make_annotation(
        "ENST00000269305", "TP53", "ENSG00000141510", "protein_coding",
        true, {ConsequenceType::MISSENSE_VARIANT}, Impact::MODERATE, "A", "T"));
    anns.push_back(make_annotation(
        "ENST00000359597", "TP53", "ENSG00000141510", "protein_coding",
        false, {ConsequenceType::SYNONYMOUS_VARIANT}, Impact::LOW, "A", "T"));
    // lncRNA (filtered by biotype)
    anns.push_back(make_annotation(
        "ENST00000500000", "LINC01234", "ENSG00000200000", "lncRNA",
        false, {ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT}, Impact::MODIFIER, "A", "T"));

    auto result = filter.filter(anns);
    // lncRNA filtered, then pick_allele_gene: one for TP53:A>T
    ASSERT_EQ(result.size(), 1u);
    EXPECT_EQ(result[0].gene_id, "ENSG00000141510");
}
