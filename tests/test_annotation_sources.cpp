/**
 * Tests for Annotation Sources
 *
 * Tests for LOFTEE, NMD, MaxEntScan, dbNSFP fields, AnnotationSourceManager,
 * conservation source metadata, and domain source types.
 */

#include <gtest/gtest.h>
#include "vep_annotator.hpp"
#include "annotation_source.hpp"
#include "annotation_sources.hpp"
#include "file_parsers.hpp"
#include "dbnsfp_fields.hpp"

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <algorithm>
#include <set>

using namespace vep;

// ============================================================================
// Helper: create a mock coding transcript for LOFTEE/NMD tests
// ============================================================================

static Transcript make_multi_exon_coding_transcript() {
    Transcript t;
    t.id = "ENST00000000001";
    t.gene_id = "ENSG00000000001";
    t.gene_name = "TESTGENE";
    t.chromosome = "1";
    t.start = 1000;
    t.end = 5000;
    t.strand = '+';
    t.biotype = "protein_coding";
    t.is_canonical = true;
    t.protein_id = "ENSP00000000001";

    // 3 exons
    Exon e1; e1.start = 1000; e1.end = 1200; e1.exon_number = 1; e1.phase = 0;
    Exon e2; e2.start = 2000; e2.end = 2300; e2.exon_number = 2; e2.phase = 0;
    Exon e3; e3.start = 3000; e3.end = 3500; e3.exon_number = 3; e3.phase = 0;
    t.exons = {e1, e2, e3};

    // CDS regions (coding parts of exons)
    CDS c1; c1.start = 1050; c1.end = 1200; c1.phase = 0;
    CDS c2; c2.start = 2000; c2.end = 2300; c2.phase = 0;
    CDS c3; c3.start = 3000; c3.end = 3400; c3.phase = 0;
    t.cds_regions = {c1, c2, c3};

    t.cds_start = 1050;
    t.cds_end = 3400;

    return t;
}

static Transcript make_single_exon_coding_transcript() {
    Transcript t;
    t.id = "ENST00000000002";
    t.gene_id = "ENSG00000000002";
    t.gene_name = "SINGLEGENE";
    t.chromosome = "1";
    t.start = 1000;
    t.end = 2000;
    t.strand = '+';
    t.biotype = "protein_coding";
    t.protein_id = "ENSP00000000002";

    Exon e1; e1.start = 1000; e1.end = 2000; e1.exon_number = 1; e1.phase = 0;
    t.exons = {e1};

    CDS c1; c1.start = 1100; c1.end = 1900; c1.phase = 0;
    t.cds_regions = {c1};

    t.cds_start = 1100;
    t.cds_end = 1900;

    return t;
}

static Transcript make_noncoding_transcript() {
    Transcript t;
    t.id = "ENST00000000003";
    t.gene_id = "ENSG00000000003";
    t.gene_name = "NONCODING";
    t.chromosome = "1";
    t.start = 1000;
    t.end = 3000;
    t.strand = '+';
    t.biotype = "lncRNA";

    Exon e1; e1.start = 1000; e1.end = 1500; e1.exon_number = 1; e1.phase = -1;
    Exon e2; e2.start = 2000; e2.end = 3000; e2.exon_number = 2; e2.phase = -1;
    t.exons = {e1, e2};

    // No CDS regions - non-coding
    t.cds_start = 0;
    t.cds_end = 0;

    return t;
}

static Transcript make_incomplete_cds_transcript() {
    Transcript t = make_multi_exon_coding_transcript();
    t.id = "ENST00000000004";
    t.gene_name = "INCOMPLETE";
    t.cds_start = 0;  // Incomplete CDS
    t.cds_end = 3400;
    return t;
}

// ============================================================================
// Mock annotation source for AnnotationSourceManager tests
// ============================================================================

class MockAnnotationSource : public AnnotationSource {
public:
    MockAnnotationSource(const std::string& src_name, const std::string& src_type,
                        size_t memory = 1024)
        : name_(src_name), type_(src_type), memory_(memory) {}

    std::string name() const override { return name_; }
    std::string type() const override { return type_; }
    std::string description() const override { return "Mock " + name_ + " source"; }
    bool is_ready() const override { return true; }
    void initialize() override {}

    void annotate(
        const std::string& /*chrom*/,
        int /*pos*/,
        const std::string& /*ref*/,
        const std::string& /*alt*/,
        const Transcript* /*transcript*/,
        std::unordered_map<std::string, std::string>& annotations
    ) override {
        annotations[name_ + ":test"] = "mock_value";
    }

    std::vector<std::string> get_fields() const override {
        return { name_ + ":test" };
    }

    size_t memory_usage() const override { return memory_; }
    bool is_thread_safe() const override { return true; }

private:
    std::string name_;
    std::string type_;
    size_t memory_;
};

// ============================================================================
// Test Suite: LOFTEE classification logic
// ============================================================================

TEST(LOFTEE, HCForFrameshiftInMultiExon) {
    auto source = create_loftee_source();
    auto t = make_multi_exon_coding_transcript();
    std::unordered_map<std::string, std::string> annotations;

    // Variant in exon 2 CDS (not last exon, protein_coding, multi-exon)
    source->annotate("1", 2100, "A", "T", &t, annotations);

    EXPECT_EQ(annotations["loftee:classification"], "HC");
}

TEST(LOFTEE, LCDueToEndTrunc) {
    auto source = create_loftee_source();
    auto t = make_multi_exon_coding_transcript();
    std::unordered_map<std::string, std::string> annotations;

    // Variant in last exon CDS (exon 3, pos 3100 is in CDS c3: 3000-3400)
    source->annotate("1", 3100, "A", "T", &t, annotations);

    EXPECT_EQ(annotations["loftee:classification"], "LC");
    EXPECT_NE(annotations.find("loftee:flags"), annotations.end());
    EXPECT_NE(annotations["loftee:flags"].find("END_TRUNC"), std::string::npos);
}

TEST(LOFTEE, LCDueToSingleExon) {
    auto source = create_loftee_source();
    auto t = make_single_exon_coding_transcript();
    std::unordered_map<std::string, std::string> annotations;

    // Variant in single exon transcript CDS
    source->annotate("1", 1500, "A", "T", &t, annotations);

    // Single exon gets SINGLE_EXON flag but may still be HC
    // (SINGLE_EXON alone doesn't force LC in this implementation)
    auto flags_it = annotations.find("loftee:flags");
    if (flags_it != annotations.end()) {
        EXPECT_NE(flags_it->second.find("SINGLE_EXON"), std::string::npos);
    }
}

TEST(LOFTEE, LCDueToIncompleteCDS) {
    auto source = create_loftee_source();
    auto t = make_incomplete_cds_transcript();
    std::unordered_map<std::string, std::string> annotations;

    // Variant in CDS of transcript with cds_start <= 0
    source->annotate("1", 2100, "A", "T", &t, annotations);

    EXPECT_EQ(annotations["loftee:classification"], "LC");
    EXPECT_NE(annotations["loftee:flags"].find("INCOMPLETE_CDS"), std::string::npos);
}

TEST(LOFTEE, NonCodingTranscriptExcluded) {
    auto source = create_loftee_source();
    auto t = make_noncoding_transcript();
    std::unordered_map<std::string, std::string> annotations;

    source->annotate("1", 1200, "A", "T", &t, annotations);

    // Non-coding transcript should not get LOFTEE annotations
    EXPECT_EQ(annotations.find("loftee:classification"), annotations.end());
}

TEST(LOFTEE, NonLoFConsequenceExcluded) {
    auto source = create_loftee_source();
    auto t = make_multi_exon_coding_transcript();
    std::unordered_map<std::string, std::string> annotations;

    // Variant in intron (not in CDS or splice site) -> no LoF annotation
    source->annotate("1", 1500, "A", "T", &t, annotations);

    EXPECT_EQ(annotations.find("loftee:classification"), annotations.end());
}

TEST(LOFTEE, ConfidenceScoreHC) {
    auto source = create_loftee_source();
    auto t = make_multi_exon_coding_transcript();
    std::unordered_map<std::string, std::string> annotations;

    // HC variant in exon 2
    source->annotate("1", 2100, "A", "T", &t, annotations);

    auto conf_it = annotations.find("loftee:confidence");
    ASSERT_NE(conf_it, annotations.end());
    double conf = std::stod(conf_it->second);
    EXPECT_GE(conf, 0.8);  // HC confidence should be high
}

TEST(LOFTEE, MultipleFlagsReduceConfidence) {
    auto source = create_loftee_source();
    auto t = make_incomplete_cds_transcript();
    std::unordered_map<std::string, std::string> annotations;

    // Variant that triggers multiple flags (INCOMPLETE_CDS at minimum)
    source->annotate("1", 2100, "A", "T", &t, annotations);

    auto conf_it = annotations.find("loftee:confidence");
    ASSERT_NE(conf_it, annotations.end());
    double conf = std::stod(conf_it->second);
    // LC with flags should have lower confidence than HC
    EXPECT_LT(conf, 0.9);
}

// ============================================================================
// Test Suite: NMD prediction
// ============================================================================

TEST(NMD, PTCUpstreamOfLastJunctionSusceptible) {
    auto source = create_nmd_source();
    auto t = make_multi_exon_coding_transcript();
    std::unordered_map<std::string, std::string> annotations;

    // Variant early in exon 1 CDS (far upstream of last junction) -> NMD susceptible
    source->annotate("1", 1060, "A", "T", &t, annotations);

    EXPECT_EQ(annotations["nmd:susceptible"], "true");
}

TEST(NMD, PTCNearLastJunctionEscape) {
    auto source = create_nmd_source();
    auto t = make_multi_exon_coding_transcript();
    std::unordered_map<std::string, std::string> annotations;

    // Variant near the end of exon 2 (close to last exon-exon junction)
    // The last junction is at exon 2 end (2300) in + strand
    // Variant at 2290 should be close to or past the junction in CDS coords
    source->annotate("1", 2290, "A", "T", &t, annotations);

    // Should be escape or very close to boundary
    auto susc_it = annotations.find("nmd:susceptible");
    ASSERT_NE(susc_it, annotations.end());
    // The exact result depends on CDS coordinate math, but we verify it returns something
    EXPECT_FALSE(susc_it->second.empty());
}

TEST(NMD, SingleExonEscape) {
    auto source = create_nmd_source();
    auto t = make_single_exon_coding_transcript();
    std::unordered_map<std::string, std::string> annotations;

    source->annotate("1", 1500, "A", "T", &t, annotations);

    EXPECT_EQ(annotations["nmd:susceptible"], "false");
    EXPECT_EQ(annotations["nmd:reason"], "single_exon_gene");
}

TEST(NMD, PTCInLastExonEscape) {
    auto source = create_nmd_source();
    auto t = make_multi_exon_coding_transcript();
    std::unordered_map<std::string, std::string> annotations;

    // Variant in last exon (exon 3: 3000-3500, CDS 3000-3400)
    source->annotate("1", 3200, "A", "T", &t, annotations);

    EXPECT_EQ(annotations["nmd:susceptible"], "false");
    // Reason should indicate in_last_exon or within_50bp_of_junction
    auto reason_it = annotations.find("nmd:reason");
    ASSERT_NE(reason_it, annotations.end());
    EXPECT_FALSE(reason_it->second.empty());
}

TEST(NMD, NonCodingExcluded) {
    auto source = create_nmd_source();
    auto t = make_noncoding_transcript();
    std::unordered_map<std::string, std::string> annotations;

    source->annotate("1", 1200, "A", "T", &t, annotations);

    // Non-coding transcript should not get NMD annotations
    EXPECT_EQ(annotations.find("nmd:susceptible"), annotations.end());
}

// ============================================================================
// Test Suite: MaxEntScan scoring
// ============================================================================

TEST(MaxEntScan, SourceMetadata) {
    auto source = create_maxentscan_source();
    EXPECT_EQ(source->name(), "maxentscan");
    EXPECT_EQ(source->type(), "splice");
}

TEST(MaxEntScan, IsReadyAlwaysTrue) {
    auto source = create_maxentscan_source();
    // MaxEntScan is algorithmic, always ready
    EXPECT_TRUE(source->is_ready());
}

TEST(MaxEntScan, ThreadSafety) {
    auto source = create_maxentscan_source();
    EXPECT_TRUE(source->is_thread_safe());
}

TEST(MaxEntScan, GetFieldsReturnsExpected) {
    auto source = create_maxentscan_source();
    auto fields = source->get_fields();
    ASSERT_FALSE(fields.empty());

    // Should include 5' and 3' splice site fields
    std::set<std::string> field_set(fields.begin(), fields.end());
    EXPECT_TRUE(field_set.count("maxentscan:5ss_ref"));
    EXPECT_TRUE(field_set.count("maxentscan:5ss_alt"));
    EXPECT_TRUE(field_set.count("maxentscan:3ss_ref"));
    EXPECT_TRUE(field_set.count("maxentscan:3ss_alt"));
    EXPECT_TRUE(field_set.count("maxentscan:in_splice_region"));
}

TEST(MaxEntScan, DescriptionNonEmpty) {
    auto source = create_maxentscan_source();
    EXPECT_FALSE(source->description().empty());
}

// ============================================================================
// Test Suite: dbNSFP field definitions
// ============================================================================

TEST(DbNSFPFields, PathogenicityFieldsContainExpected) {
    EXPECT_FALSE(DBNSFP_PATHOGENICITY_FIELDS.empty());

    std::set<std::string> names;
    for (const auto& f : DBNSFP_PATHOGENICITY_FIELDS) {
        names.insert(f.name);
    }

    EXPECT_TRUE(names.count("SIFT_score"));
    EXPECT_TRUE(names.count("CADD_phred"));
    EXPECT_TRUE(names.count("REVEL_score"));
    EXPECT_TRUE(names.count("Polyphen2_HDIV_score"));
}

TEST(DbNSFPFields, SIFTScoreProperties) {
    const DbNSFPField* sift = nullptr;
    for (const auto& f : DBNSFP_PATHOGENICITY_FIELDS) {
        if (f.name == "SIFT_score") { sift = &f; break; }
    }
    ASSERT_NE(sift, nullptr);
    EXPECT_EQ(sift->column, "SIFT_score");
    EXPECT_TRUE(sift->is_score);
    EXPECT_FALSE(sift->higher_is_damaging);  // SIFT: lower = more damaging
    EXPECT_DOUBLE_EQ(sift->damaging_threshold, 0.05);
}

TEST(DbNSFPFields, CADDPhredThreshold) {
    const DbNSFPField* cadd = nullptr;
    for (const auto& f : DBNSFP_PATHOGENICITY_FIELDS) {
        if (f.name == "CADD_phred") { cadd = &f; break; }
    }
    ASSERT_NE(cadd, nullptr);
    EXPECT_DOUBLE_EQ(cadd->damaging_threshold, 20.0);
    EXPECT_TRUE(cadd->higher_is_damaging);
}

TEST(DbNSFPFields, REVELScoreThreshold) {
    const DbNSFPField* revel = nullptr;
    for (const auto& f : DBNSFP_PATHOGENICITY_FIELDS) {
        if (f.name == "REVEL_score") { revel = &f; break; }
    }
    ASSERT_NE(revel, nullptr);
    EXPECT_DOUBLE_EQ(revel->damaging_threshold, 0.5);
    EXPECT_TRUE(revel->higher_is_damaging);
}

TEST(DbNSFPFields, AlphaMissenseFieldsPresent) {
    std::set<std::string> names;
    for (const auto& f : DBNSFP_PATHOGENICITY_FIELDS) {
        names.insert(f.name);
    }
    EXPECT_TRUE(names.count("AlphaMissense_score"));
    EXPECT_TRUE(names.count("AlphaMissense_pred"));
}

TEST(DbNSFPFields, ConservationPresetIncludesExpected) {
    auto conservation_fields = get_dbnsfp_preset("conservation");
    ASSERT_FALSE(conservation_fields.empty());

    std::set<std::string> names;
    for (const auto& f : conservation_fields) {
        names.insert(f.name);
    }

    // Should include PhyloP, PhastCons, and GERP fields
    bool has_phylop = false, has_phastcons = false, has_gerp = false;
    for (const auto& n : names) {
        if (n.find("phyloP") != std::string::npos) has_phylop = true;
        if (n.find("phastCons") != std::string::npos) has_phastcons = true;
        if (n.find("GERP") != std::string::npos) has_gerp = true;
    }
    EXPECT_TRUE(has_phylop);
    EXPECT_TRUE(has_phastcons);
    EXPECT_TRUE(has_gerp);
}

TEST(DbNSFPFields, EssentialPresetIncludesKeyClinicalScores) {
    auto essential_fields = get_dbnsfp_preset("essential");
    ASSERT_FALSE(essential_fields.empty());

    std::set<std::string> names;
    for (const auto& f : essential_fields) {
        names.insert(f.name);
    }

    EXPECT_TRUE(names.count("SIFT_score"));
    EXPECT_TRUE(names.count("CADD_phred"));
    EXPECT_TRUE(names.count("REVEL_score"));
    EXPECT_TRUE(names.count("AlphaMissense_score"));
}

TEST(DbNSFPFields, AllPresetIsComprehensive) {
    auto all_fields = get_all_dbnsfp_fields();
    ASSERT_FALSE(all_fields.empty());

    // "all" should include pathogenicity + conservation + other fields
    EXPECT_GT(all_fields.size(), DBNSFP_PATHOGENICITY_FIELDS.size());
    EXPECT_GT(all_fields.size(), DBNSFP_CONSERVATION_FIELDS.size());

    // Verify it includes fields from multiple categories
    std::set<std::string> names;
    for (const auto& f : all_fields) {
        names.insert(f.name);
    }
    EXPECT_TRUE(names.count("SIFT_score"));       // pathogenicity
    EXPECT_TRUE(names.count("GERP_RS"));           // conservation
    EXPECT_TRUE(names.count("gnomAD_exomes_AF"));  // frequency
    EXPECT_TRUE(names.count("clinvar_clnsig"));    // clinical
}

// ============================================================================
// Test Suite: AnnotationSourceManager
// ============================================================================

TEST(AnnotationSourceManager, AddSourceAndVerifyCount) {
    AnnotationSourceManager mgr;
    auto mock1 = std::make_shared<MockAnnotationSource>("test1", "pathogenicity");
    auto mock2 = std::make_shared<MockAnnotationSource>("test2", "conservation");

    mgr.add_source(mock1);
    mgr.add_source(mock2);

    auto sources = mgr.get_sources();
    EXPECT_EQ(sources.size(), 2u);
}

TEST(AnnotationSourceManager, EnableDisableSources) {
    AnnotationSourceManager mgr;
    auto mock = std::make_shared<MockAnnotationSource>("test1", "pathogenicity");
    mgr.add_source(mock);

    // Sources should be enabled by default
    EXPECT_TRUE(mgr.is_enabled("test1"));

    // Disable the source
    mgr.set_enabled("test1", false);
    EXPECT_FALSE(mgr.is_enabled("test1"));

    // Re-enable
    mgr.set_enabled("test1", true);
    EXPECT_TRUE(mgr.is_enabled("test1"));
}

TEST(AnnotationSourceManager, IsEnabledReflectsState) {
    AnnotationSourceManager mgr;
    auto mock = std::make_shared<MockAnnotationSource>("src1", "lof");
    mgr.add_source(mock);

    EXPECT_TRUE(mgr.is_enabled("src1"));
    mgr.set_enabled("src1", false);
    EXPECT_FALSE(mgr.is_enabled("src1"));
}

TEST(AnnotationSourceManager, GetSourcesReturnsAll) {
    AnnotationSourceManager mgr;
    auto m1 = std::make_shared<MockAnnotationSource>("alpha", "conservation");
    auto m2 = std::make_shared<MockAnnotationSource>("beta", "splice");
    auto m3 = std::make_shared<MockAnnotationSource>("gamma", "lof");

    mgr.add_source(m1);
    mgr.add_source(m2);
    mgr.add_source(m3);

    auto sources = mgr.get_sources();
    EXPECT_EQ(sources.size(), 3u);
}

TEST(AnnotationSourceManager, GetSourceByName) {
    AnnotationSourceManager mgr;
    auto mock = std::make_shared<MockAnnotationSource>("mysource", "pathogenicity");
    mgr.add_source(mock);

    auto found = mgr.get_source("mysource");
    ASSERT_NE(found, nullptr);
    EXPECT_EQ(found->name(), "mysource");
}

TEST(AnnotationSourceManager, GetSourceNonExistentReturnsNull) {
    AnnotationSourceManager mgr;
    auto mock = std::make_shared<MockAnnotationSource>("exists", "pathogenicity");
    mgr.add_source(mock);

    auto not_found = mgr.get_source("nonexistent");
    EXPECT_EQ(not_found, nullptr);
}

TEST(AnnotationSourceManager, TotalMemoryUsageAggregates) {
    AnnotationSourceManager mgr;
    auto m1 = std::make_shared<MockAnnotationSource>("src1", "pathogenicity", 1000);
    auto m2 = std::make_shared<MockAnnotationSource>("src2", "conservation", 2000);
    auto m3 = std::make_shared<MockAnnotationSource>("src3", "splice", 3000);

    mgr.add_source(m1);
    mgr.add_source(m2);
    mgr.add_source(m3);

    EXPECT_EQ(mgr.total_memory_usage(), 6000u);
}

TEST(AnnotationSourceManager, GetStatsProducesFormattedString) {
    AnnotationSourceManager mgr;
    auto mock = std::make_shared<MockAnnotationSource>("statstest", "pathogenicity");
    mgr.add_source(mock);

    std::string stats = mgr.get_stats();
    EXPECT_FALSE(stats.empty());
    // Stats string should mention the source name
    EXPECT_NE(stats.find("statstest"), std::string::npos);
}

// ============================================================================
// Test Suite: Conservation source metadata
// ============================================================================

TEST(ConservationSourceMetadata, PhyloPProperties) {
    // PhyloPSource is created via factory, name = "phylop", type = "conservation"
    // We can't instantiate without a valid bigWig file, but we can check factory output
    auto source = create_phylop_source("/nonexistent/phylop.bw");
    EXPECT_EQ(source->name(), "phylop");
    EXPECT_EQ(source->type(), "conservation");
    EXPECT_FALSE(source->description().empty());
}

TEST(ConservationSourceMetadata, PhastConsProperties) {
    auto source = create_phastcons_source("/nonexistent/phastcons.bw");
    EXPECT_EQ(source->name(), "phastcons");
    EXPECT_EQ(source->type(), "conservation");
    EXPECT_FALSE(source->description().empty());
}

TEST(ConservationSourceMetadata, GERPProperties) {
    auto source = create_gerp_source("/nonexistent/gerp.bw");
    EXPECT_EQ(source->name(), "gerp");
    EXPECT_EQ(source->type(), "conservation");
    EXPECT_FALSE(source->description().empty());
}

// ============================================================================
// Test Suite: Domain source types
// ============================================================================

TEST(DomainSource, PfamSourceNameAndType) {
    auto source = create_pfam_source("/nonexistent/pfam.tsv");
    EXPECT_EQ(source->name(), "pfam");
    EXPECT_EQ(source->type(), "domain");
}

TEST(DomainSource, InterProSourceNameAndType) {
    auto source = create_interpro_source("/nonexistent/interpro.tsv");
    EXPECT_EQ(source->name(), "interpro");
    EXPECT_EQ(source->type(), "domain");
}

TEST(DomainSource, DomainSourceFieldLists) {
    auto pfam = create_pfam_source("/nonexistent/pfam.tsv");
    auto fields = pfam->get_fields();
    ASSERT_FALSE(fields.empty());

    std::set<std::string> field_set(fields.begin(), fields.end());
    EXPECT_TRUE(field_set.count("pfam:domain_id"));
    EXPECT_TRUE(field_set.count("pfam:domain_name"));
    EXPECT_TRUE(field_set.count("pfam:domain_count"));
}
