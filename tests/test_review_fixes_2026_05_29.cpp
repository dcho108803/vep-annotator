/**
 * Regression tests for the HIGH-severity fixes from the 2026-05-29 code review.
 *
 *  - #1/#2 Minus-strand multi-base variant spanning the 3'UTR/CDS boundary was
 *          mis-clamped to CDS position 1 and falsely reported as start_lost.
 *  - #3    generate_hgvsg dropped the inserted base for single-base-side delins
 *          (e.g. AGG>C emitted a plain deletion).
 *  - #6/#7 Transcript-selecting filters (--biotype, --coding-only, ...) ran
 *          against a single most-severe transcript; TranscriptFilterConfig::
 *          requires_all_transcripts() now drives "annotate all transcripts".
 *  - dead  Gene-constraint (pLI/LOEUF) was loaded but never queried; the query
 *          API the wiring depends on is exercised here.
 */

#include <gtest/gtest.h>
#include "vep_annotator.hpp"
#include "transcript_filter.hpp"
#include "hgvs_parser.hpp"
#include "gene_constraint.hpp"
#include "structural_variant.hpp"
#include "exon_intron_numbers.hpp"

#include <fstream>
#include <string>
#include <cstdio>
#include <algorithm>
#include <map>

using namespace vep;

// ============================================================================
// #3 - generate_hgvsg must not drop inserted bases for single-base-side delins
// ============================================================================

TEST(HGVSgDelinsFix, SingleBaseAltDelinsKeepsInsertedBase) {
    // REF=ATG (3 bases), ALT=C (single base, does NOT match anchor 'A').
    // This is a delins, not a deletion. The pre-fix code emitted "100_102del"
    // and silently dropped the inserted C.
    std::string result = generate_hgvsg("7", 100, "ATG", "C");
    EXPECT_NE(result.find("delins"), std::string::npos) << result;
    EXPECT_NE(result.find("100_102delinsC"), std::string::npos) << result;
}

TEST(HGVSgDelinsFix, SingleBaseRefDelinsKeepsAllInsertedBases) {
    // REF=A (single base), ALT=GT (does NOT match anchor) -> delins, not insertion.
    // Pre-fix emitted an insertion of "T" and dropped the A->G replacement.
    std::string result = generate_hgvsg("7", 100, "A", "GT");
    EXPECT_NE(result.find("delins"), std::string::npos) << result;
    EXPECT_NE(result.find("100delinsGT"), std::string::npos) << result;
}

TEST(HGVSgDelinsFix, AnchoredDeletionStillDeletion) {
    // REF=ATG, ALT=A (alt matches anchor) is still a pure deletion.
    std::string result = generate_hgvsg("7", 100, "ATG", "A");
    EXPECT_NE(result.find("del"), std::string::npos) << result;
    EXPECT_EQ(result.find("delins"), std::string::npos) << result;
    EXPECT_NE(result.find("101_102"), std::string::npos) << result;
}

TEST(HGVSgDelinsFix, AnchoredInsertionStillInsertion) {
    // REF=A, ALT=ATG (alt starts with anchor) is still a pure insertion.
    std::string result = generate_hgvsg("7", 100, "A", "ATG");
    EXPECT_NE(result.find("ins"), std::string::npos) << result;
    EXPECT_EQ(result.find("delins"), std::string::npos) << result;
}

// ============================================================================
// #6/#7 - TranscriptFilterConfig::requires_all_transcripts()
// ============================================================================

TEST(RequiresAllTranscripts, EmptyConfigDoesNotRequireAll) {
    TranscriptFilterConfig config;
    EXPECT_FALSE(config.requires_all_transcripts());
}

TEST(RequiresAllTranscripts, BiotypeFilterRequiresAll) {
    // The exact bug #7: --biotype without a pick flag.
    TranscriptFilterConfig config;
    config.biotypes.insert("protein_coding");
    EXPECT_TRUE(config.requires_all_transcripts());
}

TEST(RequiresAllTranscripts, TranscriptSelectingFlagsRequireAll) {
    {
        TranscriptFilterConfig c; c.coding_only = true;
        EXPECT_TRUE(c.requires_all_transcripts());
    }
    {
        TranscriptFilterConfig c; c.canonical_only = true;
        EXPECT_TRUE(c.requires_all_transcripts());
    }
    {
        TranscriptFilterConfig c; c.mane_only = true;
        EXPECT_TRUE(c.requires_all_transcripts());
    }
    {
        TranscriptFilterConfig c; c.no_intergenic = true;
        EXPECT_TRUE(c.requires_all_transcripts());
    }
    {
        TranscriptFilterConfig c; c.exclude_predicted = true;
        EXPECT_TRUE(c.requires_all_transcripts());
    }
    {
        TranscriptFilterConfig c; c.gencode_basic = true;
        EXPECT_TRUE(c.requires_all_transcripts());
    }
    {
        TranscriptFilterConfig c; c.pick = true;
        EXPECT_TRUE(c.requires_all_transcripts());
    }
    {
        TranscriptFilterConfig c; c.most_severe = true;
        EXPECT_TRUE(c.requires_all_transcripts());
    }
    {
        TranscriptFilterConfig c; c.include_consequences.insert("missense_variant");
        EXPECT_TRUE(c.requires_all_transcripts());
    }
    {
        TranscriptFilterConfig c; c.exclude_consequences.insert("intron_variant");
        EXPECT_TRUE(c.requires_all_transcripts());
    }
}

TEST(RequiresAllTranscripts, FrequencyOnlyDoesNotRequireAll) {
    // Frequency filtering is a variant-level decision with the same outcome
    // regardless of how many transcripts are present.
    TranscriptFilterConfig config;
    config.check_frequency = true;
    EXPECT_FALSE(config.requires_all_transcripts());
}

// ============================================================================
// Gene constraint (pLI / LOEUF) - query API the output wiring depends on
// ============================================================================

TEST(GeneConstraintQuery, LoadAndQueryBySymbol) {
    std::string path = std::string(::testing::TempDir()) + "test_pli_scores.txt";
    {
        std::ofstream out(path);
        out << "gene\tpLI\n";          // header (skipped)
        out << "TESTGENE\t0.95\n";
        out << "OTHERGENE\t0.01\n";
    }
    GeneConstraintDB db;
    ASSERT_TRUE(db.load_pli_scores(path));
    EXPECT_TRUE(db.is_loaded());

    GeneConstraint gc = db.get_by_symbol("TESTGENE");
    EXPECT_TRUE(gc.has_data());
    EXPECT_DOUBLE_EQ(gc.pLI, 0.95);
    EXPECT_EQ(format_constraint_score(gc.pLI), "0.95");

    // Absent gene -> no data, formats to "."
    GeneConstraint missing = db.get_by_symbol("NOSUCHGENE");
    EXPECT_FALSE(missing.has_data());
    EXPECT_EQ(format_constraint_score(missing.pLI), ".");

    std::remove(path.c_str());
}

TEST(GeneConstraintQuery, LoadLoeufAndFormat) {
    std::string path = std::string(::testing::TempDir()) + "test_loeuf_scores.txt";
    {
        std::ofstream out(path);
        out << "gene\toe_lof_upper\n";  // header (skipped)
        out << "TESTGENE\t0.25\n";
    }
    GeneConstraintDB db;
    ASSERT_TRUE(db.load_loeuf_scores(path));

    GeneConstraint gc = db.get_by_symbol("TESTGENE");
    EXPECT_TRUE(gc.has_data());
    EXPECT_DOUBLE_EQ(gc.oe_lof_upper, 0.25);
    EXPECT_EQ(format_constraint_score(gc.oe_lof_upper), "0.25");

    std::remove(path.c_str());
}

// ============================================================================
// #1/#2 - Minus-strand multi-base variant at the 3'UTR/CDS boundary
// ============================================================================
//
// Build a synthetic minus-strand coding transcript via a temp GTF + FASTA and
// run the full annotator. On the minus strand the 3'UTR sits at LOW genomic
// coordinates (below cds_start). A multi-base variant whose genomic-leftmost
// base is in the 3'UTR used to be mis-clamped to CDS position 1 and reported
// as start_lost; it must not be.

namespace {

// Deterministic reference base for a 1-based position (repeating ACGT pattern),
// so test variant REF alleles match the synthetic genome.
char base_at(int pos1) {
    static const char* P = "ACGT";
    return P[(pos1 - 1) % 4];
}

std::string ref_substr(int start1, int len) {
    std::string s;
    for (int i = 0; i < len; ++i) s += base_at(start1 + i);
    return s;
}

struct MinusStrandFixture {
    std::string gtf_path;
    std::string fasta_path;

    MinusStrandFixture() {
        const std::string dir = ::testing::TempDir();
        gtf_path = dir + "minus_strand_fix.gtf";
        fasta_path = dir + "minus_strand_fix.fa";

        // Transcript on chr20, minus strand, single exon 1000-2000,
        // CDS 1100-1900 (length 801, divisible by 3 -> no incomplete codon).
        //   3'UTR (minus strand) = genomic [1000,1099] (below cds_start)
        //   5'UTR (minus strand) = genomic [1901,2000] (above cds_end)
        {
            std::ofstream g(gtf_path);
            const char* attrs =
                "gene_id \"ENSGTEST\"; transcript_id \"ENSTTEST\"; "
                "gene_name \"MINUSGENE\"; gene_biotype \"protein_coding\"; "
                "transcript_biotype \"protein_coding\";";
            g << "20\ttest\tgene\t1000\t2000\t.\t-\t.\t" << attrs << "\n";
            g << "20\ttest\ttranscript\t1000\t2000\t.\t-\t.\t" << attrs << "\n";
            g << "20\ttest\texon\t1000\t2000\t.\t-\t.\t" << attrs << "\n";
            g << "20\ttest\tCDS\t1100\t1900\t.\t-\t0\t" << attrs << "\n";
        }
        {
            std::ofstream f(fasta_path);
            f << ">20\n";
            std::string seq;
            for (int p = 1; p <= 2100; ++p) seq += base_at(p);
            // wrap at 70 cols
            for (size_t i = 0; i < seq.size(); i += 70) {
                f << seq.substr(i, 70) << "\n";
            }
        }
    }
    ~MinusStrandFixture() {
        std::remove(gtf_path.c_str());
        std::remove(fasta_path.c_str());
    }
};

// Plus-strand single-exon coding transcript on chr20, CDS 1100-1900, ACGT-pattern
// genome (so codons are deterministic): CDS starts TAC GTA CGT ... = Tyr Val Arg ...
struct PlusStrandFixture {
    std::string gtf_path;
    std::string fasta_path;
    PlusStrandFixture() {
        const std::string dir = ::testing::TempDir();
        gtf_path = dir + "plus_strand_fix.gtf";
        fasta_path = dir + "plus_strand_fix.fa";
        {
            std::ofstream g(gtf_path);
            const char* attrs =
                "gene_id \"ENSGP\"; transcript_id \"ENSTP\"; gene_name \"PLUSGENE\"; "
                "gene_biotype \"protein_coding\"; transcript_biotype \"protein_coding\";";
            g << "20\ttest\tgene\t1000\t2000\t.\t+\t.\t" << attrs << "\n";
            g << "20\ttest\ttranscript\t1000\t2000\t.\t+\t.\t" << attrs << "\n";
            g << "20\ttest\texon\t1000\t2000\t.\t+\t.\t" << attrs << "\n";
            g << "20\ttest\tCDS\t1100\t1900\t.\t+\t0\t" << attrs << "\n";
        }
        {
            std::ofstream f(fasta_path);
            f << ">20\n";
            std::string seq;
            for (int p = 1; p <= 2100; ++p) seq += base_at(p);
            for (size_t i = 0; i < seq.size(); i += 70) f << seq.substr(i, 70) << "\n";
        }
    }
    ~PlusStrandFixture() {
        std::remove(gtf_path.c_str());
        std::remove(fasta_path.c_str());
    }
};

bool has_consequence(const std::vector<VariantAnnotation>& anns, ConsequenceType c) {
    for (const auto& a : anns) {
        if (std::find(a.consequences.begin(), a.consequences.end(), c) != a.consequences.end()) {
            return true;
        }
    }
    return false;
}

} // namespace

TEST(MinusStrandBoundary, VariantSpanning3UtrCdsBoundaryIsNotStartLost) {
    MinusStrandFixture fx;
    VEPAnnotator annotator(fx.gtf_path, fx.fasta_path);

    // Deletion spanning the 3'UTR/CDS boundary at the LOW genomic end.
    // pos=1098 (3'UTR), REF covers 1098..1100 (1100 is the first CDS base scanned,
    // which maps to the highest CDS position ~ the stop region), ALT=anchor.
    int pos = 1098;
    std::string ref = ref_substr(pos, 3);   // 1098,1099,1100
    std::string alt(1, ref[0]);             // anchored deletion of 1099,1100

    auto anns = annotator.annotate("20", pos, ref, alt);
    ASSERT_FALSE(anns.empty());

    // The regression: this must NOT be start_lost (it is near the stop, not start).
    EXPECT_FALSE(has_consequence(anns, ConsequenceType::START_LOST))
        << "minus-strand 3'UTR/CDS boundary variant was wrongly called start_lost";
}

TEST(MinusStrandBoundary, RealStartCodonVariantStillStartLost) {
    // Positive control: a deletion at the actual start codon (minus strand, near
    // cds_end) must still be detected as start_lost after the fix.
    MinusStrandFixture fx;
    VEPAnnotator annotator(fx.gtf_path, fx.fasta_path);

    int pos = 1898;                          // 1899,1900 are CDS positions 2,1 (start codon)
    std::string ref = ref_substr(pos, 3);    // 1898,1899,1900
    std::string alt(1, ref[0]);              // anchored deletion of 1899,1900

    auto anns = annotator.annotate("20", pos, ref, alt);
    ASSERT_FALSE(anns.empty());

    EXPECT_TRUE(has_consequence(anns, ConsequenceType::START_LOST))
        << "minus-strand start-codon deletion should still be start_lost";
}

// ============================================================================
// HGVS parser additions (inversion, protein delins/dup/ext, regex_match, flags)
// ============================================================================

TEST(HGVSParseReviewFixes, ProteinDelins) {
    auto r = parse_protein_hgvs("ENSP1", "Cys28delinsTrpVal");
    EXPECT_TRUE(r.valid);
    EXPECT_EQ(r.variant_type, HGVSVariantType::DELINS);
    EXPECT_EQ(r.protein_pos, 28);
    EXPECT_EQ(r.alt_aa, "TrpVal");
}

TEST(HGVSParseReviewFixes, ProteinDuplication) {
    auto r = parse_protein_hgvs("ENSP1", "Gly4_Gln6dup");
    EXPECT_TRUE(r.valid);
    EXPECT_EQ(r.variant_type, HGVSVariantType::DUPLICATION);
    EXPECT_EQ(r.start_pos, 4);
    EXPECT_EQ(r.end_pos, 6);
}

TEST(HGVSParseReviewFixes, ProteinExtension) {
    auto r = parse_protein_hgvs("ENSP1", "Ter110GlnextTer17");
    EXPECT_TRUE(r.valid);
    EXPECT_EQ(r.ref_aa, "Ter");
    EXPECT_EQ(r.protein_pos, 110);
}

TEST(HGVSParseReviewFixes, ProteinInsertionRejectsGarbage) {
    // regex_match (not search): well-formed accepted, leading garbage rejected.
    auto good = parse_protein_hgvs("ENSP1", "Gly12_Ala13insValTrp");
    EXPECT_TRUE(good.valid);
    EXPECT_EQ(good.alt_aa, "ValTrp");
    auto bad = parse_protein_hgvs("ENSP1", "garbage_Gly12_Ala13insVal");
    EXPECT_FALSE(bad.valid);
}

TEST(HGVSParseReviewFixes, CodingAndGenomicInversion) {
    auto c = parse_coding_hgvs("ENST1", "100_110inv");
    EXPECT_TRUE(c.valid);
    EXPECT_EQ(c.variant_type, HGVSVariantType::INVERSION);
    auto g = parse_genomic_hgvs("NC_000001.11", "100_110inv");
    EXPECT_TRUE(g.valid);
    EXPECT_EQ(g.variant_type, HGVSVariantType::INVERSION);
}

TEST(HGVSParseReviewFixes, ThreePrimeUtrStarFlag) {
    auto r = parse_coding_hgvs("ENST1", "*51A>G");
    EXPECT_TRUE(r.valid);
    EXPECT_TRUE(r.is_utr3);
    EXPECT_EQ(r.start_pos, 51);
    // A plain c.51 must NOT be flagged as 3'UTR.
    auto plain = parse_coding_hgvs("ENST1", "51A>G");
    EXPECT_FALSE(plain.is_utr3);
}

TEST(HGVSParseReviewFixes, EndIntronOffsetForInsertion) {
    auto r = parse_coding_hgvs("ENST1", "100+1_100+2insA");
    EXPECT_TRUE(r.valid);
    EXPECT_EQ(r.intron_offset, 1);
    EXPECT_EQ(r.end_intron_offset, 2);
}

// ============================================================================
// Transcript pick: deterministic tiebreaker on full ties
// ============================================================================

TEST(TranscriptPickReviewFixes, FullTieResolvedByTranscriptId) {
    TranscriptFilterConfig cfg;
    cfg.pick = true;
    TranscriptFilter filter(cfg);

    auto make = [](const std::string& tid) {
        VariantAnnotation a;
        a.transcript_id = tid;
        a.gene_id = "G1";
        a.biotype = "protein_coding";
        a.consequences = { ConsequenceType::MISSENSE_VARIANT };
        return a;
    };
    // Input order deliberately reversed; tie must resolve to the lexicographically
    // smaller transcript_id regardless of input order or std::sort instability.
    std::vector<VariantAnnotation> anns = { make("ENST00000002"), make("ENST00000001") };
    auto picked = filter.filter(anns);
    ASSERT_EQ(picked.size(), 1u);
    EXPECT_EQ(picked[0].transcript_id, "ENST00000001");
}

// ============================================================================
// Structural variants: strand-aware splice, CDS-overlap frame, BND bracket
// ============================================================================

namespace {
Transcript make_minus_sv_transcript() {
    Transcript t;
    t.id = "ENSTSV"; t.gene_id = "GSV"; t.gene_name = "SVGENE";
    t.chromosome = "1"; t.strand = '-'; t.biotype = "protein_coding";
    // Two exons: 1000-1100 and 2000-2100. Intron 1101-1999.
    Exon e1; e1.start = 1000; e1.end = 1100; e1.exon_number = 2; t.exons.push_back(e1);
    Exon e2; e2.start = 2000; e2.end = 2100; e2.exon_number = 1; t.exons.push_back(e2);
    CDS c1; c1.start = 1000; c1.end = 1100; t.cds_regions.push_back(c1);
    CDS c2; c2.start = 2000; c2.end = 2100; t.cds_regions.push_back(c2);
    t.start = 1000; t.end = 2100; t.cds_start = 1000; t.cds_end = 2100;
    return t;
}
}  // namespace

TEST(SVReviewFixes, MinusStrandSpliceDonorAcceptorSwapped) {
    Transcript t = make_minus_sv_transcript();
    // A small deletion hitting the high-coordinate boundary of exon e1 (exon.end+1..+2
    // = 1101..1102). On the MINUS strand that genomic side is the ACCEPTOR.
    std::map<std::string, std::string> info = {{"SVTYPE", "DEL"}, {"END", "1102"}};
    auto sv = parse_sv_from_vcf("1", 1100, "N", "<DEL>", info);
    auto cons = get_sv_consequences(sv, t);
    bool has_acceptor = std::find(cons.begin(), cons.end(),
        ConsequenceType::SPLICE_ACCEPTOR_VARIANT) != cons.end();
    bool has_donor = std::find(cons.begin(), cons.end(),
        ConsequenceType::SPLICE_DONOR_VARIANT) != cons.end();
    EXPECT_TRUE(has_acceptor);
    EXPECT_FALSE(has_donor);
}

TEST(SVReviewFixes, BndOrientationFromBracketChar) {
    std::map<std::string, std::string> info = {{"SVTYPE", "BND"}};
    // ']' bracket -> reverse join orientation (bnd_mate_forward == false), even
    // though the sequence precedes the bracket (old position-only logic said true).
    auto sv = parse_sv_from_vcf("1", 1000, "N", "N]2:5000]", info);
    EXPECT_EQ(sv.sv_type, SVType::BND);
    EXPECT_FALSE(sv.bnd_mate_forward);
    auto sv2 = parse_sv_from_vcf("1", 1000, "N", "N[2:5000[", info);
    EXPECT_TRUE(sv2.bnd_mate_forward);
}

// ============================================================================
// Exon/intron numbering: strand-aware position_in_feature
// ============================================================================

// ============================================================================
// Multi-codon MNV -> protein delins HGVSp (not a single truncated residue)
// ============================================================================

TEST(MultiCodonMNV, TwoCodonSubstitutionEmitsDelins) {
    PlusStrandFixture fx;
    VEPAnnotator annotator(fx.gtf_path, fx.fasta_path);

    // CDS codons 2 and 3 are GTA (Val) and CGT (Arg), at genomic 1103-1108.
    // Replace with AAA (Lys) and CCC (Pro): both residues change -> delins.
    std::string ref = ref_substr(1103, 6);   // "GTACGT"
    ASSERT_EQ(ref, "GTACGT");
    std::string alt = "AAACCC";

    auto anns = annotator.annotate("20", 1103, ref, alt);
    ASSERT_FALSE(anns.empty());
    const auto& a = anns[0];

    // amino_acids carries both changed residues on each side.
    EXPECT_EQ(a.amino_acids, "VR/KP");
    // HGVSp must be a delins over the two-residue range, not a single missense.
    EXPECT_NE(a.hgvsp.find("delins"), std::string::npos) << a.hgvsp;
    EXPECT_NE(a.hgvsp.find("Val2_Arg3delinsLysPro"), std::string::npos) << a.hgvsp;
}

TEST(MultiCodonMNV, SingleChangedResidueInMnvIsMissense) {
    PlusStrandFixture fx;
    VEPAnnotator annotator(fx.gtf_path, fx.fasta_path);

    // 6bp MNV spanning codons 2,3 but only changing codon 3 (CGT Arg -> CCC Pro):
    // codon 2 GTA is left intact. After trimming, this is a single missense at pos 3.
    std::string ref = ref_substr(1103, 6);   // "GTACGT"
    std::string alt = "GTACCC";              // codon2 unchanged, codon3 CGT->CCC

    auto anns = annotator.annotate("20", 1103, ref, alt);
    ASSERT_FALSE(anns.empty());
    const auto& a = anns[0];
    EXPECT_EQ(a.amino_acids, "VR/VP");
    // Trimmed to the single changed residue (Arg3 -> Pro), no delins.
    EXPECT_EQ(a.hgvsp.find("delins"), std::string::npos) << a.hgvsp;
    EXPECT_NE(a.hgvsp.find("Arg3Pro"), std::string::npos) << a.hgvsp;
}

TEST(MultiCodonMNV, StopGainingMnvEmitsTer) {
    PlusStrandFixture fx;
    VEPAnnotator annotator(fx.gtf_path, fx.fasta_path);

    // 6bp MNV over codons 2,3: codon2 GTA (Val) kept, codon3 CGT (Arg) -> TAA (stop).
    // HGVSp must be stop_gained at the stop position, not a missense/synonymous.
    std::string ref = ref_substr(1103, 6);   // "GTACGT"
    std::string alt = "GTATAA";

    auto anns = annotator.annotate("20", 1103, ref, alt);
    ASSERT_FALSE(anns.empty());
    const auto& a = anns[0];
    EXPECT_EQ(a.amino_acids, "VR/V*");
    EXPECT_NE(a.hgvsp.find("Arg3Ter"), std::string::npos) << a.hgvsp;
    EXPECT_EQ(a.hgvsp.find("delins"), std::string::npos) << a.hgvsp;
}

TEST(ExonIntronReviewFixes, MinusStrandPositionInFeature) {
    // Single exon 1000..1100 on the minus strand. position_in_feature is measured
    // from the high (transcript-5') coordinate: pos 1098 -> 1100-1098+1 = 3.
    std::vector<int> starts = {1000};
    std::vector<int> ends = {1100};
    auto info = get_exon_intron_number(1098, starts, ends, '-');
    ASSERT_TRUE(info.found);
    EXPECT_TRUE(info.is_exon);
    EXPECT_EQ(info.position_in_feature, 3);
    // Plus strand from the low coordinate: 1098-1000+1 = 99.
    auto info_plus = get_exon_intron_number(1098, starts, ends, '+');
    EXPECT_EQ(info_plus.position_in_feature, 99);
}
