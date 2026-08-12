/**
 * CNV annotation feature tests: ClinGen dosage sensitivity (HI/TS scores),
 * population SV frequency, and ACMG CNV classification.
 */

#include <gtest/gtest.h>
#include <fstream>
#include <cstdio>
#include "gene_constraint.hpp"

using namespace vep;

namespace {

std::string write_temp(const std::string& name, const std::string& content) {
    std::string path = std::string(::testing::TempDir()) + name;
    std::ofstream out(path);
    out << content;
    return path;
}

} // namespace

// ============================================================================
// ClinGen dosage loading (--clingen-dosage)
// ============================================================================

TEST(ClinGenDosage, LoadsClinGenCurationFormat) {
    // Real ClinGen files carry #-prefixed prose plus a #-prefixed header row
    std::string path = write_temp("clingen_dosage.tsv",
        "#ClinGen Gene Curation Results\n"
        "#Downloaded from ftp.clinicalgenome.org\n"
        "#Gene Symbol\tGene ID\tcytoBand\tGenomic Location\t"
        "Haploinsufficiency Score\tHaploinsufficiency Description\t"
        "Triplosensitivity Score\tTriplosensitivity Description\n"
        "BRCA1\t672\t17q21.31\tchr17:43044295-43125483\t3\tSufficient evidence\t0\tNo evidence\n"
        "PMP22\t5376\t17p12\tchr17:15229777-15265326\t3\tSufficient evidence\t3\tSufficient evidence\n"
        "SHANK3\t85358\t22q13.33\tchr22:50674415-50733212\t3\tSufficient evidence\tNot yet evaluated\t\n"
        "AR\t367\tXq12\tchrX:67544021-67730619\t30\tGene associated with AR phenotype\t0\tNo evidence\n");

    GeneConstraintDB db;
    ASSERT_TRUE(db.load_clingen_dosage(path));
    EXPECT_TRUE(db.has_dosage_data());

    GeneConstraint brca1 = db.get_by_symbol("BRCA1");
    EXPECT_DOUBLE_EQ(brca1.hi_score, 3.0);
    EXPECT_DOUBLE_EQ(brca1.ts_score, 0.0);
    EXPECT_TRUE(brca1.has_data());

    GeneConstraint pmp22 = db.get_by_symbol("PMP22");
    EXPECT_DOUBLE_EQ(pmp22.hi_score, 3.0);
    EXPECT_DOUBLE_EQ(pmp22.ts_score, 3.0);

    // "Not yet evaluated" TS leaves the score unknown; numeric HI still loads
    GeneConstraint shank3 = db.get_by_symbol("SHANK3");
    EXPECT_DOUBLE_EQ(shank3.hi_score, 3.0);
    EXPECT_LT(shank3.ts_score, 0);

    // ClinGen code 30 (autosomal recessive) passes through
    EXPECT_DOUBLE_EQ(db.get_by_symbol("AR").hi_score, 30.0);

    std::remove(path.c_str());
}

TEST(ClinGenDosage, LoadsSimpleThreeColumnFormat) {
    std::string path = write_temp("dosage_simple.tsv",
        "GENE\tHI\tTS\n"
        "TESTGENE\t3\t0\n"
        "OTHERGENE\t2\t1\n");

    GeneConstraintDB db;
    ASSERT_TRUE(db.load_clingen_dosage(path));
    EXPECT_DOUBLE_EQ(db.get_by_symbol("TESTGENE").hi_score, 3.0);
    EXPECT_DOUBLE_EQ(db.get_by_symbol("OTHERGENE").ts_score, 1.0);

    std::remove(path.c_str());
}

TEST(ClinGenDosage, MergesIntoExistingConstraintEntry) {
    // pLI loaded first; dosage curations must merge, not replace
    std::string pli_path = write_temp("pli_for_merge.txt",
        "gene\tpLI\n"
        "TESTGENE\t0.98\n");
    std::string dosage_path = write_temp("dosage_for_merge.tsv",
        "TESTGENE\t3\t0\n");

    GeneConstraintDB db;
    ASSERT_TRUE(db.load_pli_scores(pli_path));
    ASSERT_TRUE(db.load_clingen_dosage(dosage_path));

    GeneConstraint gc = db.get_by_symbol("TESTGENE");
    EXPECT_DOUBLE_EQ(gc.pLI, 0.98);          // preserved
    EXPECT_DOUBLE_EQ(gc.hi_score, 3.0);      // merged
    EXPECT_DOUBLE_EQ(gc.ts_score, 0.0);

    std::remove(pli_path.c_str());
    std::remove(dosage_path.c_str());
}

TEST(ClinGenDosage, DosageOnlyGeneReportsHasData) {
    // has_data() must be true for a gene known only from dosage curations, so
    // the gene-ID -> symbol lookup fallback works for it
    std::string path = write_temp("dosage_only.tsv", "DOSAGEGENE\t2\t.\n");

    GeneConstraintDB db;
    ASSERT_TRUE(db.load_clingen_dosage(path));
    GeneConstraint gc = db.get_by_symbol("DOSAGEGENE");
    EXPECT_TRUE(gc.has_data());
    EXPECT_DOUBLE_EQ(gc.hi_score, 2.0);
    EXPECT_LT(gc.ts_score, 0);

    std::remove(path.c_str());
}

TEST(ClinGenDosage, NoDosageDataWithoutFile) {
    GeneConstraintDB db;
    EXPECT_FALSE(db.has_dosage_data());
    std::string pli_path = write_temp("pli_no_dosage.txt", "TESTGENE\t0.5\n");
    ASSERT_TRUE(db.load_pli_scores(pli_path));
    EXPECT_FALSE(db.has_dosage_data());
    std::remove(pli_path.c_str());
}

// ============================================================================
// Population SV frequency (--sv-frequency)
// ============================================================================

#include "sv_frequency.hpp"

namespace {

StructuralVariant make_sv(SVType type, int start, int end, double cn = -1) {
    StructuralVariant sv;
    sv.chromosome = "1";
    sv.sv_type = type;
    sv.start = start;
    sv.end = end;
    sv.copy_number = cn;
    return sv;
}

} // namespace

TEST(SVFrequency, LoadsTsvAndMatchesReciprocal) {
    std::string path = write_temp("sv_freq.tsv",
        "#CHROM\tSTART\tEND\tTYPE\tAF\tID\n"
        "1\t1000\t5000\tDEL\t0.012\tgnomAD-SV_v2_DEL_1\n"
        "1\t20000\t30000\tDUP\t0.34\tgnomAD-SV_v2_DUP_1\n"
        "2\t1000\t5000\tDEL\t0.5\tchr2_del\n");

    SVFrequencyDB db;
    ASSERT_TRUE(db.load(path));
    EXPECT_EQ(db.size(), 3u);

    // Near-identical DEL matches with high reciprocal overlap
    SVFrequencyDB::Match m;
    ASSERT_TRUE(db.best_match(make_sv(SVType::DEL, 1100, 5000), m, 0.7));
    EXPECT_EQ(m.record.id, "gnomAD-SV_v2_DEL_1");
    EXPECT_DOUBLE_EQ(m.record.af, 0.012);
    EXPECT_GT(m.overlap, 0.9);

    // Same interval as a DUP does not match the DEL record
    EXPECT_FALSE(db.best_match(make_sv(SVType::DUP, 1100, 5000), m, 0.7));

    // Insufficient reciprocal overlap fails the cutoff
    EXPECT_FALSE(db.best_match(make_sv(SVType::DEL, 4000, 20000), m, 0.7));

    std::remove(path.c_str());
}

TEST(SVFrequency, LoadsVcfFormat) {
    std::string path = write_temp("sv_freq.vcf",
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t999\tgnomAD_DEL_42\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=5000;AF=0.0123\n"
        "1\t9999\tgnomAD_INS_7\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=300;AF=0.002\n");

    SVFrequencyDB db;
    ASSERT_TRUE(db.load(path));
    EXPECT_EQ(db.size(), 2u);

    SVFrequencyDB::Match m;
    ASSERT_TRUE(db.best_match(make_sv(SVType::DEL, 1000, 5000), m, 0.7));
    EXPECT_EQ(m.record.id, "gnomAD_DEL_42");
    EXPECT_DOUBLE_EQ(m.record.af, 0.0123);

    std::remove(path.c_str());
}

TEST(SVFrequency, CnvQueryMatchesByCopyNumberClass) {
    std::string path = write_temp("sv_freq_cnv.tsv",
        "1\t1000\t5000\tDEL\t0.01\tdel_1\n"
        "1\t1000\t5000\tDUP\t0.02\tdup_1\n");

    SVFrequencyDB db;
    ASSERT_TRUE(db.load(path));

    // CNV with CN=1 is a loss: matches the DEL record
    SVFrequencyDB::Match m;
    ASSERT_TRUE(db.best_match(make_sv(SVType::CNV, 1000, 5000, 1), m, 0.7));
    EXPECT_EQ(m.record.id, "del_1");

    // CNV with CN=3 is a gain: matches the DUP record
    ASSERT_TRUE(db.best_match(make_sv(SVType::CNV, 1000, 5000, 3), m, 0.7));
    EXPECT_EQ(m.record.id, "dup_1");

    std::remove(path.c_str());
}

TEST(SVFrequency, BestOfMultipleCandidates) {
    std::string path = write_temp("sv_freq_multi.tsv",
        "1\t900\t5500\tDEL\t0.10\tloose\n"
        "1\t1000\t5000\tDEL\t0.01\ttight\n");

    SVFrequencyDB db;
    ASSERT_TRUE(db.load(path));

    // Exact interval prefers the highest-reciprocal-overlap record
    SVFrequencyDB::Match m;
    ASSERT_TRUE(db.best_match(make_sv(SVType::DEL, 1000, 5000), m, 0.5));
    EXPECT_EQ(m.record.id, "tight");
    EXPECT_DOUBLE_EQ(m.overlap, 1.0);

    std::remove(path.c_str());
}

TEST(SVFrequency, OverlapModes) {
    std::string path = write_temp("sv_freq_modes.tsv",
        "1\t2000\t3000\tDEL\t0.05\tinner\n");

    SVFrequencyDB db;
    ASSERT_TRUE(db.load(path));
    SVFrequencyDB::Match m;

    // SURROUNDING: query contains the record
    EXPECT_TRUE(db.best_match(make_sv(SVType::DEL, 1000, 5000), m, 0.7,
                              OverlapType::SURROUNDING));
    // WITHIN: query inside the record — not satisfied here
    EXPECT_FALSE(db.best_match(make_sv(SVType::DEL, 1000, 5000), m, 0.7,
                               OverlapType::WITHIN));
    // EXACT requires identical coordinates
    EXPECT_TRUE(db.best_match(make_sv(SVType::DEL, 2000, 3000), m, 0.7,
                              OverlapType::EXACT));
    EXPECT_FALSE(db.best_match(make_sv(SVType::DEL, 2000, 3001), m, 0.7,
                               OverlapType::EXACT));
    // ANY accepts a 1bp overlap
    EXPECT_TRUE(db.best_match(make_sv(SVType::DEL, 3000, 9000), m, 0.7,
                              OverlapType::ANY));

    std::remove(path.c_str());
}

TEST(SVFrequency, ChromosomeNamingNormalized) {
    std::string path = write_temp("sv_freq_chr.tsv",
        "chr1\t1000\t5000\tDEL\t0.01\twith_prefix\n");

    SVFrequencyDB db;
    ASSERT_TRUE(db.load(path));
    SVFrequencyDB::Match m;
    ASSERT_TRUE(db.best_match(make_sv(SVType::DEL, 1000, 5000), m, 0.7));
    EXPECT_EQ(m.record.id, "with_prefix");

    std::remove(path.c_str());
}

// ============================================================================
// ACMG CNV classification (--acmg-cnv)
// ============================================================================

#include "acmg_cnv.hpp"

namespace {

Transcript make_tx(const std::string& gene, int start, int end,
                   bool coding = true, const std::string& id_suffix = "T1") {
    Transcript t;
    t.id = gene + "_" + id_suffix;
    t.gene_id = "ENSG_" + gene;
    t.gene_name = gene;
    t.chromosome = "1";
    t.start = start;
    t.end = end;
    t.strand = '+';
    if (coding) {
        t.cds_start = start + 100;
        t.cds_end = end - 100;
        t.cds_regions.push_back({t.cds_start, t.cds_end});
    }
    t.biotype = coding ? "protein_coding" : "lincRNA";
    return t;
}

GeneConstraintDB make_dosage_db(const std::string& rows) {
    std::string path = write_temp("acmg_dosage.tsv", rows);
    GeneConstraintDB db;
    db.load_clingen_dosage(path);
    std::remove(path.c_str());
    return db;
}

} // namespace

TEST(AcmgCnv, LossContainingHiGeneIsPathogenic) {
    GeneConstraintDB db = make_dosage_db("HIGENE\t3\t0\n");
    Transcript tx = make_tx("HIGENE", 2000, 8000);
    std::vector<const Transcript*> txs = {&tx};

    auto sv = make_sv(SVType::DEL, 1000, 10000);
    auto res = classify_acmg_cnv(sv, txs, db);
    ASSERT_TRUE(res.applicable);
    EXPECT_DOUBLE_EQ(res.score, 1.0);   // 1A(0) + 2A(+1.0) + 3A(0)
    EXPECT_EQ(res.classification, "pathogenic");
    // Evidence trail names the gene
    bool has_2a = false;
    for (const auto& e : res.evidence) {
        if (e.find("2A:HIGENE") != std::string::npos) has_2a = true;
    }
    EXPECT_TRUE(has_2a);
}

TEST(AcmgCnv, PartialLossOfHiGeneCodingIsLikelyPathogenic) {
    GeneConstraintDB db = make_dosage_db("HIGENE\t3\t0\n");
    Transcript tx = make_tx("HIGENE", 2000, 8000);
    std::vector<const Transcript*> txs = {&tx};

    // Loss removes the 5' half of the gene including CDS
    auto sv = make_sv(SVType::DEL, 1000, 5000);
    auto res = classify_acmg_cnv(sv, txs, db);
    EXPECT_DOUBLE_EQ(res.score, 0.90);  // 2C-2E(+0.90)
    EXPECT_EQ(res.classification, "likely_pathogenic");
}

TEST(AcmgCnv, GeneFreeLossIsUncertain) {
    GeneConstraintDB db;
    std::vector<const Transcript*> txs;
    auto res = classify_acmg_cnv(make_sv(SVType::DEL, 1000, 5000), txs, db);
    ASSERT_TRUE(res.applicable);
    EXPECT_DOUBLE_EQ(res.score, -0.60);  // 1B
    EXPECT_EQ(res.classification, "uncertain_significance");
}

TEST(AcmgCnv, CommonLossIsBenign) {
    GeneConstraintDB db;
    Transcript tx = make_tx("SOMEGENE", 2000, 4000);
    std::vector<const Transcript*> txs = {&tx};

    // No dosage evidence, but the SV matches a 5% population variant
    auto res = classify_acmg_cnv(make_sv(SVType::DEL, 1000, 5000), txs, db, 0.05);
    EXPECT_DOUBLE_EQ(res.score, -1.0);   // 1A + 3A + 4O(-1.0)
    EXPECT_EQ(res.classification, "benign");
}

TEST(AcmgCnv, GainContainingTsGeneIsPathogenic) {
    GeneConstraintDB db = make_dosage_db("TSGENE\t0\t3\n");
    Transcript tx = make_tx("TSGENE", 2000, 8000);
    std::vector<const Transcript*> txs = {&tx};

    auto res = classify_acmg_cnv(make_sv(SVType::DUP, 1000, 10000), txs, db);
    EXPECT_DOUBLE_EQ(res.score, 1.0);
    EXPECT_EQ(res.classification, "pathogenic");
}

TEST(AcmgCnv, HiGeneDoesNotScoreForGains) {
    // Haploinsufficiency evidence applies to losses, not duplications
    GeneConstraintDB db = make_dosage_db("HIGENE\t3\t0\n");
    Transcript tx = make_tx("HIGENE", 2000, 8000);
    std::vector<const Transcript*> txs = {&tx};

    auto res = classify_acmg_cnv(make_sv(SVType::DUP, 1000, 10000), txs, db);
    EXPECT_DOUBLE_EQ(res.score, 0.0);
    EXPECT_EQ(res.classification, "uncertain_significance");
}

TEST(AcmgCnv, LossInsideDosageInsensitiveGeneIsBenign) {
    // ClinGen code 40 = dosage sensitivity unlikely; a CNV fully inside it
    // scores the benign 2F item
    GeneConstraintDB db = make_dosage_db("BENIGNGENE\t40\t40\n");
    Transcript tx = make_tx("BENIGNGENE", 1000, 20000);
    std::vector<const Transcript*> txs = {&tx};

    auto res = classify_acmg_cnv(make_sv(SVType::DEL, 5000, 6000), txs, db);
    EXPECT_DOUBLE_EQ(res.score, -1.0);   // 1A + 2F(-1.0) + 3A
    EXPECT_EQ(res.classification, "benign");
}

TEST(AcmgCnv, ManyGenesAddSection3Points) {
    GeneConstraintDB db;
    std::vector<Transcript> storage;
    storage.reserve(40);
    std::vector<const Transcript*> txs;
    for (int i = 0; i < 36; ++i) {
        storage.push_back(make_tx("G" + std::to_string(i),
                                  1000 + i * 100, 1000 + i * 100 + 50));
    }
    for (const auto& t : storage) txs.push_back(&t);

    // 36 genes deleted: 3C (+0.90) for losses
    auto res = classify_acmg_cnv(make_sv(SVType::DEL, 500, 10000), txs, db);
    EXPECT_DOUBLE_EQ(res.score, 0.90);
    EXPECT_EQ(res.classification, "likely_pathogenic");

    // The same 36 genes duplicated: gains need 35-49 for 3B (+0.45)
    auto res_dup = classify_acmg_cnv(make_sv(SVType::DUP, 500, 10000), txs, db);
    EXPECT_DOUBLE_EQ(res_dup.score, 0.45);
    EXPECT_EQ(res_dup.classification, "uncertain_significance");
}

TEST(AcmgCnv, CnvClassFromCopyNumber) {
    GeneConstraintDB db = make_dosage_db("HIGENE\t3\t0\n");
    Transcript tx = make_tx("HIGENE", 2000, 8000);
    std::vector<const Transcript*> txs = {&tx};

    // DRAGEN-style SVTYPE=CNV with CN=1 behaves as a loss
    auto res = classify_acmg_cnv(make_sv(SVType::CNV, 1000, 10000, 1), txs, db);
    EXPECT_EQ(res.classification, "pathogenic");

    // Unknown copy number: not classifiable
    auto res_unknown = classify_acmg_cnv(make_sv(SVType::CNV, 1000, 10000), txs, db);
    EXPECT_FALSE(res_unknown.applicable);
}

TEST(AcmgCnv, NotApplicableToInversionsAndBreakends) {
    GeneConstraintDB db;
    std::vector<const Transcript*> txs;
    EXPECT_FALSE(classify_acmg_cnv(make_sv(SVType::INV, 1000, 5000), txs, db).applicable);
    EXPECT_FALSE(classify_acmg_cnv(make_sv(SVType::BND, 1000, 1000), txs, db).applicable);
    EXPECT_FALSE(classify_acmg_cnv(make_sv(SVType::INS, 1000, 1000), txs, db).applicable);
}
