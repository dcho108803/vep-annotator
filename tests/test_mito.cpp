/**
 * Mitochondrial annotation tests
 *
 * End-to-end coverage for MT-specific behavior:
 * - chrM (UCSC) vs MT (Ensembl) contig naming interop through the full
 *   GTF/FASTA load + annotate path
 * - Mt_tRNA / Mt_rRNA biotype consequence calling
 * - Mitochondrial genetic code (table 2) dispatch in coding consequences
 * - m. HGVS output for MT variants
 */

#include <gtest/gtest.h>
#include <fstream>
#include <cstdio>
#include "vep_annotator.hpp"
#include "mito_annotation.hpp"

using namespace vep;

namespace {

// Synthetic mitochondrial genome (1-based, length 800), UCSC-named "chrM":
//   MT-TF   Mt_tRNA          gene/exon 100-171   (no CDS)
//   MT-RNR1 Mt_rRNA          gene/exon 200-350   (no CDS)
//   MT-ND1  protein_coding   gene/exon 400-450, CDS 400-450 (17 codons)
// CDS layout: ATG | TGG | GCT x14 | TAA  -> M W A...A *
// Everything outside the crafted CDS is 'A'.
struct MitoFixture {
    std::string gtf_path;
    std::string fasta_path;

    MitoFixture() {
        const std::string dir = ::testing::TempDir();
        gtf_path = dir + "mito_test.gtf";
        fasta_path = dir + "mito_test.fa";

        auto attrs = [](const char* gid, const char* tid, const char* name,
                        const char* biotype) {
            std::string a = "gene_id \"";
            a += gid;
            a += "\"; transcript_id \"";
            a += tid;
            a += "\"; gene_name \"";
            a += name;
            a += "\"; gene_biotype \"";
            a += biotype;
            a += "\"; transcript_biotype \"";
            a += biotype;
            a += "\";";
            return a;
        };

        {
            std::ofstream g(gtf_path);
            std::string trna = attrs("ENSG_MTTF", "ENST_MTTF", "MT-TF", "Mt_tRNA");
            g << "chrM\ttest\tgene\t100\t171\t.\t+\t.\t" << trna << "\n";
            g << "chrM\ttest\ttranscript\t100\t171\t.\t+\t.\t" << trna << "\n";
            g << "chrM\ttest\texon\t100\t171\t.\t+\t.\t" << trna << "\n";

            std::string rrna = attrs("ENSG_RNR1", "ENST_RNR1", "MT-RNR1", "Mt_rRNA");
            g << "chrM\ttest\tgene\t200\t350\t.\t+\t.\t" << rrna << "\n";
            g << "chrM\ttest\ttranscript\t200\t350\t.\t+\t.\t" << rrna << "\n";
            g << "chrM\ttest\texon\t200\t350\t.\t+\t.\t" << rrna << "\n";

            std::string nd1 = attrs("ENSG_ND1", "ENST_ND1", "MT-ND1", "protein_coding");
            g << "chrM\ttest\tgene\t400\t450\t.\t+\t.\t" << nd1 << "\n";
            g << "chrM\ttest\ttranscript\t400\t450\t.\t+\t.\t" << nd1 << "\n";
            g << "chrM\ttest\texon\t400\t450\t.\t+\t.\t" << nd1 << "\n";
            g << "chrM\ttest\tCDS\t400\t450\t.\t+\t0\t" << nd1 << "\n";
        }
        {
            std::string seq(800, 'A');
            auto put = [&](int pos1, const std::string& s) {
                for (size_t i = 0; i < s.size(); ++i) seq[pos1 - 1 + i] = s[i];
            };
            put(400, "ATG");
            put(403, "TGG");
            for (int c = 0; c < 14; ++c) put(406 + 3 * c, "GCT");
            put(448, "TAA");

            std::ofstream f(fasta_path);
            f << ">chrM\n";
            for (size_t i = 0; i < seq.size(); i += 70) {
                f << seq.substr(i, 70) << "\n";
            }
        }
    }
    ~MitoFixture() {
        std::remove(gtf_path.c_str());
        std::remove(fasta_path.c_str());
    }
};

const VariantAnnotation* find_ann(const std::vector<VariantAnnotation>& anns,
                                  const std::string& transcript_id) {
    for (const auto& a : anns) {
        if (a.transcript_id == transcript_id) return &a;
    }
    return nullptr;
}

bool ann_has(const VariantAnnotation& a, ConsequenceType c) {
    return std::find(a.consequences.begin(), a.consequences.end(), c) !=
           a.consequences.end();
}

} // namespace

TEST(MitoBiotypes, TrnaVariantIsNonCodingExon) {
    MitoFixture fx;
    VEPAnnotator annotator(fx.gtf_path, fx.fasta_path);

    auto anns = annotator.annotate("MT", 130, "A", "G");
    const VariantAnnotation* trna = find_ann(anns, "ENST_MTTF");
    ASSERT_NE(trna, nullptr) << "tRNA transcript not annotated";

    EXPECT_EQ(trna->biotype, "Mt_tRNA");
    EXPECT_TRUE(ann_has(*trna, ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT));
    // Mt_tRNA must not take the miRNA special case
    EXPECT_FALSE(ann_has(*trna, ConsequenceType::MATURE_MIRNA_VARIANT));
    EXPECT_EQ(trna->impact, Impact::MODIFIER);
}

TEST(MitoBiotypes, RrnaVariantIsNonCodingExon) {
    MitoFixture fx;
    VEPAnnotator annotator(fx.gtf_path, fx.fasta_path);

    auto anns = annotator.annotate("MT", 250, "A", "G");
    const VariantAnnotation* rrna = find_ann(anns, "ENST_RNR1");
    ASSERT_NE(rrna, nullptr) << "rRNA transcript not annotated";

    EXPECT_EQ(rrna->biotype, "Mt_rRNA");
    EXPECT_TRUE(ann_has(*rrna, ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT));
    EXPECT_FALSE(ann_has(*rrna, ConsequenceType::MATURE_MIRNA_VARIANT));
}

TEST(MitoNaming, UcscAndEnsemblQueriesMatch) {
    // GTF/FASTA are chrM-named; both MT and chrM queries must hit the
    // same transcripts and produce m. HGVS on the rCRS accession.
    MitoFixture fx;
    VEPAnnotator annotator(fx.gtf_path, fx.fasta_path);

    auto anns_mt = annotator.annotate("MT", 130, "A", "G");
    auto anns_chrm = annotator.annotate("chrM", 130, "A", "G");

    ASSERT_FALSE(anns_mt.empty());
    EXPECT_EQ(anns_mt.size(), anns_chrm.size());
    ASSERT_NE(find_ann(anns_mt, "ENST_MTTF"), nullptr);
    ASSERT_NE(find_ann(anns_chrm, "ENST_MTTF"), nullptr);

    EXPECT_EQ(find_ann(anns_mt, "ENST_MTTF")->hgvsg, "NC_012920.1:m.130A>G");
    EXPECT_EQ(find_ann(anns_chrm, "ENST_MTTF")->hgvsg, "NC_012920.1:m.130A>G");
}

TEST(MitoGeneticCode, TgaIsTrpNotStop) {
    // Codon 2 of MT-ND1 is TGG (Trp). m.405G>A turns it into TGA, which is
    // a stop in the nuclear code but tryptophan in the vertebrate mito code:
    // the end-to-end call must be synonymous, never stop_gained.
    MitoFixture fx;
    VEPAnnotator annotator(fx.gtf_path, fx.fasta_path);

    auto anns = annotator.annotate("MT", 405, "G", "A");
    const VariantAnnotation* nd1 = find_ann(anns, "ENST_ND1");
    ASSERT_NE(nd1, nullptr) << "MT protein-coding transcript not annotated";

    EXPECT_FALSE(ann_has(*nd1, ConsequenceType::STOP_GAINED))
        << "TGA on chrM was treated as a nuclear stop codon";
    EXPECT_TRUE(ann_has(*nd1, ConsequenceType::SYNONYMOUS_VARIANT));
}

TEST(MitoGeneticCode, MissenseSanity) {
    // m.404G>C: TGG (Trp) -> TCG (Ser), missense under both codes.
    MitoFixture fx;
    VEPAnnotator annotator(fx.gtf_path, fx.fasta_path);

    auto anns = annotator.annotate("MT", 404, "G", "C");
    const VariantAnnotation* nd1 = find_ann(anns, "ENST_ND1");
    ASSERT_NE(nd1, nullptr);

    EXPECT_TRUE(ann_has(*nd1, ConsequenceType::MISSENSE_VARIANT));
}

// ============================================================================
// Heteroplasmy fraction (compute_heteroplasmy)
// ============================================================================

TEST(Heteroplasmy, FormatAfPassthrough) {
    // Caller-reported AF is passed through verbatim
    EXPECT_EQ(compute_heteroplasmy("GT:AF:DP\t0/1:0.853:1200", 1), "0.853");
}

TEST(Heteroplasmy, MultiAllelicAfSelection) {
    // AF is Number=A: one entry per alt allele, chosen by allele index
    EXPECT_EQ(compute_heteroplasmy("GT:AF\t1/2:0.60,0.35", 1), "0.60");
    EXPECT_EQ(compute_heteroplasmy("GT:AF\t1/2:0.60,0.35", 2), "0.35");
    // Allele index past the AF entries -> unavailable
    EXPECT_EQ(compute_heteroplasmy("GT:AF\t1/2:0.60,0.35", 3), "");
}

TEST(Heteroplasmy, AdFallback) {
    // No AF: derived from AD (Number=R: ref,alt1,...) as alt/sum
    EXPECT_EQ(compute_heteroplasmy("GT:AD\t0/1:25,75", 1), "0.75");
    // Multi-allelic AD: entry allele_index (ref is entry 0)
    EXPECT_EQ(compute_heteroplasmy("GT:AD\t1/2:10,60,30", 2), "0.3");
}

TEST(Heteroplasmy, FirstSampleOnly) {
    // Mito calling is single-sample; values come from the first sample column
    EXPECT_EQ(compute_heteroplasmy("GT:AF\t0/1:0.90\t0/0:0.05", 1), "0.90");
}

TEST(Heteroplasmy, Unavailable) {
    EXPECT_EQ(compute_heteroplasmy("", 1), "");
    EXPECT_EQ(compute_heteroplasmy("GT:DP\t0/1:1200", 1), "");        // no AF/AD
    EXPECT_EQ(compute_heteroplasmy("GT:AF\t0/1:.", 1), "");           // missing value
    EXPECT_EQ(compute_heteroplasmy("GT:AF\t0/1:abc", 1), "");         // malformed
    EXPECT_EQ(compute_heteroplasmy("GT:AD\t0/1:0,0", 1), "");         // zero depth
    EXPECT_EQ(compute_heteroplasmy("GT", 1), "");                     // no sample
}

TEST(Heteroplasmy, MitoChromDetection) {
    EXPECT_TRUE(is_mito_chrom("MT"));
    EXPECT_TRUE(is_mito_chrom("chrM"));
    EXPECT_TRUE(is_mito_chrom("chrMT"));
    EXPECT_TRUE(is_mito_chrom("M"));
    EXPECT_FALSE(is_mito_chrom("1"));
    EXPECT_FALSE(is_mito_chrom("X"));
}

// ============================================================================
// Mito disease database (--mitomap)
// ============================================================================

namespace {

std::string write_mito_temp(const std::string& name, const std::string& content) {
    std::string path = std::string(::testing::TempDir()) + name;
    std::ofstream out(path);
    out << content;
    return path;
}

} // namespace

TEST(MitoDiseaseDB, LoadsHeaderedTsv) {
    std::string path = write_mito_temp("mitomap.tsv",
        "pos\tref\talt\tdisease\tstatus\n"
        "3243\tA\tG\tMELAS / Leigh Syndrome\tCfrm\n"
        "8993\tT\tG\tNARP / Leigh Syndrome\tCfrm\n");

    MitoDiseaseDB db;
    ASSERT_TRUE(db.load(path));
    EXPECT_EQ(db.size(), 2u);

    const auto* e = db.lookup(3243, "A", "G");
    ASSERT_NE(e, nullptr);
    // Spaces are sanitized to keep Extra-column values single-token
    EXPECT_EQ(e->disease, "MELAS_/_Leigh_Syndrome");
    EXPECT_EQ(e->status, "Cfrm");

    EXPECT_EQ(db.lookup(3243, "A", "T"), nullptr);  // different alt
    EXPECT_EQ(db.lookup(9999, "A", "G"), nullptr);  // absent position

    std::remove(path.c_str());
}

TEST(MitoDiseaseDB, LoadsMitomapAlleleColumn) {
    // MITOMAP exports carry a combined allele column
    std::string path = write_mito_temp("mitomap_allele.tsv",
        "Position\tLocus\tAllele\tDisease\tStatus\n"
        "3243\tMT-TL1\tA3243G\tMELAS\tCfrm\n"
        "11778\tMT-ND4\tm.11778G>A\tLHON\tCfrm\n");

    MitoDiseaseDB db;
    ASSERT_TRUE(db.load(path));
    ASSERT_NE(db.lookup(3243, "A", "G"), nullptr);
    ASSERT_NE(db.lookup(11778, "G", "A"), nullptr);
    EXPECT_EQ(db.lookup(11778, "G", "A")->disease, "LHON");

    std::remove(path.c_str());
}

TEST(MitoDiseaseDB, LoadsHeaderlessTsvAndVcf) {
    std::string tsv_path = write_mito_temp("mito_plain.tsv",
        "3243\tA\tG\tMELAS\tCfrm\n");
    MitoDiseaseDB tsv_db;
    ASSERT_TRUE(tsv_db.load(tsv_path));
    ASSERT_NE(tsv_db.lookup(3243, "A", "G"), nullptr);
    std::remove(tsv_path.c_str());

    std::string vcf_path = write_mito_temp("mito.vcf",
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "MT\t8993\t.\tT\tG,C\t.\t.\tDisease=NARP;CLNSIG=Pathogenic\n");
    MitoDiseaseDB vcf_db;
    ASSERT_TRUE(vcf_db.load(vcf_path));
    // Multi-allelic records key each alt separately
    ASSERT_NE(vcf_db.lookup(8993, "T", "G"), nullptr);
    ASSERT_NE(vcf_db.lookup(8993, "T", "C"), nullptr);
    EXPECT_EQ(vcf_db.lookup(8993, "T", "G")->disease, "NARP");
    EXPECT_EQ(vcf_db.lookup(8993, "T", "G")->status, "Pathogenic");
    std::remove(vcf_path.c_str());
}

TEST(MitoDiseaseDB, CaseInsensitiveAlleleLookup) {
    std::string path = write_mito_temp("mito_case.tsv", "3243\ta\tg\tMELAS\tCfrm\n");
    MitoDiseaseDB db;
    ASSERT_TRUE(db.load(path));
    EXPECT_NE(db.lookup(3243, "A", "G"), nullptr);
    std::remove(path.c_str());
}
