/**
 * Tests for output_writer.hpp: TSVWriter, JSONWriter, VCFWriter,
 * AnnotationStats, utility functions, and factory function.
 */

#include <gtest/gtest.h>
#include "output_writer.hpp"

#include <fstream>
#include <sstream>
#include <filesystem>
#include <cstdio>

using namespace vep;

// ============================================================================
// Helper: create a minimal VariantAnnotation for testing
// ============================================================================

static VariantAnnotation make_basic_annotation() {
    VariantAnnotation ann;
    ann.input_variant = "chr17\t7675088\t.\tC\tT\t.\t.\t.";
    ann.chromosome = "chr17";
    ann.position = 7675088;
    ann.ref_allele = "C";
    ann.alt_allele = "T";
    ann.gene_symbol = "TP53";
    ann.gene_id = "ENSG00000141510";
    ann.transcript_id = "ENST00000269305";
    ann.biotype = "protein_coding";
    ann.feature_type = "Transcript";
    ann.is_canonical = true;
    ann.consequences = {ConsequenceType::MISSENSE_VARIANT};
    ann.impact = Impact::MODERATE;
    ann.strand = '-';
    ann.cdna_position = 782;
    ann.cds_position = 535;
    ann.protein_position = 179;
    ann.amino_acids = "H/R";
    ann.codons = "cAt/cGt";
    ann.hgvsc = "ENST00000269305.9:c.535A>G";
    ann.hgvsp = "ENSP00000269305.4:p.His179Arg";
    ann.hgvsg = "chr17:g.7675088C>T";
    ann.existing_variation = "rs28934576";
    ann.display_alt = "T";
    ann.exon_number = 5;
    ann.total_exons = 11;
    ann.source = "Ensembl";
    return ann;
}

static VariantAnnotation make_empty_annotation() {
    VariantAnnotation ann;
    ann.chromosome = "chr1";
    ann.position = 100;
    ann.ref_allele = "A";
    ann.alt_allele = "G";
    ann.consequences = {ConsequenceType::INTERGENIC_VARIANT};
    ann.impact = Impact::MODIFIER;
    ann.feature_type = "Intergenic";
    return ann;
}

static VariantAnnotation make_multibase_annotation() {
    VariantAnnotation ann;
    ann.chromosome = "chr17";
    ann.position = 7675080;
    ann.ref_allele = "ATCG";
    ann.alt_allele = "A";
    ann.gene_symbol = "TP53";
    ann.gene_id = "ENSG00000141510";
    ann.transcript_id = "ENST00000269305";
    ann.biotype = "protein_coding";
    ann.feature_type = "Transcript";
    ann.consequences = {ConsequenceType::FRAMESHIFT_VARIANT};
    ann.impact = Impact::HIGH;
    ann.strand = '-';
    ann.cdna_position = 790;
    ann.cdna_end = 792;
    ann.cds_position = 543;
    ann.cds_end = 545;
    ann.protein_position = 181;
    ann.protein_end = 182;
    ann.display_ref = "TCG";
    ann.display_alt = "-";
    ann.display_start = 7675081;
    ann.display_end = 7675083;
    ann.custom_annotations["CDS_LENGTH"] = "1182";
    ann.custom_annotations["TRANSCRIPT_LENGTH"] = "2591";
    return ann;
}

// Helper: read entire file contents
static std::string read_file(const std::string& path) {
    std::ifstream f(path);
    std::string content((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());
    return content;
}

// RAII temp file cleanup
class TempFile {
public:
    TempFile(const std::string& suffix = ".txt") {
        path_ = std::filesystem::temp_directory_path() / ("test_output_" + std::to_string(counter_++) + suffix);
    }
    ~TempFile() {
        std::filesystem::remove(path_);
    }
    std::string path() const { return path_.string(); }
private:
    std::filesystem::path path_;
    static inline int counter_ = 0;
};


// ============================================================================
// OutputFormat parsing
// ============================================================================

TEST(OutputFormat, ParseTSV) {
    EXPECT_EQ(parse_output_format("tsv"), OutputFormat::TSV);
    EXPECT_EQ(parse_output_format("TSV"), OutputFormat::TSV);
    EXPECT_EQ(parse_output_format("Tsv"), OutputFormat::TSV);
}

TEST(OutputFormat, ParseJSON) {
    EXPECT_EQ(parse_output_format("json"), OutputFormat::JSON);
    EXPECT_EQ(parse_output_format("JSON"), OutputFormat::JSON);
    EXPECT_EQ(parse_output_format("Json"), OutputFormat::JSON);
}

TEST(OutputFormat, ParseVCF) {
    EXPECT_EQ(parse_output_format("vcf"), OutputFormat::VCF);
    EXPECT_EQ(parse_output_format("VCF"), OutputFormat::VCF);
    EXPECT_EQ(parse_output_format("Vcf"), OutputFormat::VCF);
}

TEST(OutputFormat, ParseDefaultsTSV) {
    EXPECT_EQ(parse_output_format("tab"), OutputFormat::TSV);
    EXPECT_EQ(parse_output_format(""), OutputFormat::TSV);
    EXPECT_EQ(parse_output_format("unknown"), OutputFormat::TSV);
}


// ============================================================================
// AnnotationStats
// ============================================================================

TEST(AnnotationStats, InitialState) {
    AnnotationStats stats;
    EXPECT_EQ(stats.total_variants, 0);
    EXPECT_EQ(stats.annotated_variants, 0);
    EXPECT_TRUE(stats.consequence_counts.empty());
    EXPECT_TRUE(stats.impact_counts.empty());
    EXPECT_TRUE(stats.biotype_counts.empty());
}

TEST(AnnotationStats, AddAnnotatedVariant) {
    AnnotationStats stats;
    VariantAnnotation ann = make_basic_annotation();
    stats.add(ann);

    EXPECT_EQ(stats.total_variants, 1);
    EXPECT_EQ(stats.annotated_variants, 1);  // gene_symbol is non-empty
    EXPECT_EQ(stats.consequence_counts["missense_variant"], 1);
    EXPECT_EQ(stats.impact_counts["MODERATE"], 1);
    EXPECT_EQ(stats.biotype_counts["protein_coding"], 1);
}

TEST(AnnotationStats, AddUnannotatedVariant) {
    AnnotationStats stats;
    VariantAnnotation ann = make_empty_annotation();
    stats.add(ann);

    EXPECT_EQ(stats.total_variants, 1);
    EXPECT_EQ(stats.annotated_variants, 0);  // gene_symbol is empty
    EXPECT_EQ(stats.consequence_counts["intergenic_variant"], 1);
    EXPECT_EQ(stats.impact_counts["MODIFIER"], 1);
    EXPECT_TRUE(stats.biotype_counts.empty());  // biotype is empty
}

TEST(AnnotationStats, AddMultipleVariants) {
    AnnotationStats stats;

    VariantAnnotation ann1 = make_basic_annotation();
    VariantAnnotation ann2 = make_basic_annotation();
    ann2.consequences = {ConsequenceType::SYNONYMOUS_VARIANT};
    ann2.impact = Impact::LOW;
    ann2.biotype = "protein_coding";

    VariantAnnotation ann3 = make_empty_annotation();

    stats.add(ann1);
    stats.add(ann2);
    stats.add(ann3);

    EXPECT_EQ(stats.total_variants, 3);
    EXPECT_EQ(stats.annotated_variants, 2);
    EXPECT_EQ(stats.consequence_counts["missense_variant"], 1);
    EXPECT_EQ(stats.consequence_counts["synonymous_variant"], 1);
    EXPECT_EQ(stats.consequence_counts["intergenic_variant"], 1);
    EXPECT_EQ(stats.impact_counts["MODERATE"], 1);
    EXPECT_EQ(stats.impact_counts["LOW"], 1);
    EXPECT_EQ(stats.impact_counts["MODIFIER"], 1);
    EXPECT_EQ(stats.biotype_counts["protein_coding"], 2);
}

TEST(AnnotationStats, MultipleConsequencesPerVariant) {
    AnnotationStats stats;
    VariantAnnotation ann = make_basic_annotation();
    ann.consequences = {ConsequenceType::MISSENSE_VARIANT, ConsequenceType::SPLICE_REGION_VARIANT};
    stats.add(ann);

    EXPECT_EQ(stats.total_variants, 1);
    EXPECT_EQ(stats.consequence_counts["missense_variant"], 1);
    EXPECT_EQ(stats.consequence_counts["splice_region_variant"], 1);
}

TEST(AnnotationStats, ToStringFormat) {
    AnnotationStats stats;
    VariantAnnotation ann = make_basic_annotation();
    stats.add(ann);

    std::string s = stats.to_string();
    EXPECT_NE(s.find("=== Annotation Statistics ==="), std::string::npos);
    EXPECT_NE(s.find("Total variants: 1"), std::string::npos);
    EXPECT_NE(s.find("Annotated variants: 1"), std::string::npos);
    EXPECT_NE(s.find("missense_variant: 1"), std::string::npos);
    EXPECT_NE(s.find("MODERATE: 1"), std::string::npos);
}

TEST(AnnotationStats, ToJsonFormat) {
    AnnotationStats stats;
    VariantAnnotation ann = make_basic_annotation();
    stats.add(ann);

    std::string j = stats.to_json();
    EXPECT_NE(j.find("\"total_variants\": 1"), std::string::npos);
    EXPECT_NE(j.find("\"annotated_variants\": 1"), std::string::npos);
    EXPECT_NE(j.find("\"missense_variant\": 1"), std::string::npos);
    EXPECT_NE(j.find("\"MODERATE\": 1"), std::string::npos);
}


// ============================================================================
// Utility: format_position_with_total
// ============================================================================

TEST(FormatPosition, BasicPosition) {
    VariantAnnotation ann;
    EXPECT_EQ(format_position_with_total(100, 0, "CDS_LENGTH", ann, "-"), "100");
}

TEST(FormatPosition, PositionWithRange) {
    VariantAnnotation ann;
    EXPECT_EQ(format_position_with_total(100, 105, "CDS_LENGTH", ann, "-"), "100-105");
}

TEST(FormatPosition, PositionWithTotal) {
    VariantAnnotation ann;
    ann.custom_annotations["CDS_LENGTH"] = "1182";
    EXPECT_EQ(format_position_with_total(100, 0, "CDS_LENGTH", ann, "-"), "100/1182");
}

TEST(FormatPosition, PositionWithRangeAndTotal) {
    VariantAnnotation ann;
    ann.custom_annotations["CDS_LENGTH"] = "1182";
    EXPECT_EQ(format_position_with_total(100, 105, "CDS_LENGTH", ann, "-"), "100-105/1182");
}

TEST(FormatPosition, ZeroStartReturnsDash) {
    VariantAnnotation ann;
    EXPECT_EQ(format_position_with_total(0, 0, "CDS_LENGTH", ann, "-"), "-");
}

TEST(FormatPosition, NegativeStartReturnsDash) {
    VariantAnnotation ann;
    EXPECT_EQ(format_position_with_total(-5, 0, "CDS_LENGTH", ann, "-"), "-");
}

TEST(FormatPosition, ZeroStartReturnsEmptyForVCF) {
    VariantAnnotation ann;
    EXPECT_EQ(format_position_with_total(0, 0, "CDS_LENGTH", ann, ""), "");
}

TEST(FormatPosition, SameStartEnd) {
    VariantAnnotation ann;
    // When start == end, no range displayed
    EXPECT_EQ(format_position_with_total(100, 100, "CDS_LENGTH", ann, "-"), "100");
}


// ============================================================================
// Utility: format_protein_position
// ============================================================================

TEST(FormatProteinPosition, BasicPosition) {
    VariantAnnotation ann;
    ann.protein_position = 179;
    EXPECT_EQ(format_protein_position(ann, "-"), "179");
}

TEST(FormatProteinPosition, PositionWithRange) {
    VariantAnnotation ann;
    ann.protein_position = 181;
    ann.protein_end = 182;
    EXPECT_EQ(format_protein_position(ann, "-"), "181-182");
}

TEST(FormatProteinPosition, PositionWithCDSLength) {
    VariantAnnotation ann;
    ann.protein_position = 179;
    ann.custom_annotations["CDS_LENGTH"] = "1182";
    // CDS_LENGTH / 3 = 394
    EXPECT_EQ(format_protein_position(ann, "-"), "179/394");
}

TEST(FormatProteinPosition, ZeroReturnsDash) {
    VariantAnnotation ann;
    ann.protein_position = 0;
    EXPECT_EQ(format_protein_position(ann, "-"), "-");
}

TEST(FormatProteinPosition, ZeroReturnsEmptyForVCF) {
    VariantAnnotation ann;
    ann.protein_position = 0;
    EXPECT_EQ(format_protein_position(ann, ""), "");
}

TEST(FormatProteinPosition, InvalidCDSLength) {
    VariantAnnotation ann;
    ann.protein_position = 100;
    ann.custom_annotations["CDS_LENGTH"] = "abc";
    // Invalid CDS_LENGTH should be silently ignored (try-catch)
    EXPECT_EQ(format_protein_position(ann, "-"), "100");
}


// ============================================================================
// Utility: ends_with_gz
// ============================================================================

TEST(EndsWithGz, GzExtension) {
    EXPECT_TRUE(ends_with_gz("output.vcf.gz"));
    EXPECT_TRUE(ends_with_gz("test.tsv.gz"));
    EXPECT_TRUE(ends_with_gz("file.gz"));
}

TEST(EndsWithGz, NoGzExtension) {
    EXPECT_FALSE(ends_with_gz("output.vcf"));
    EXPECT_FALSE(ends_with_gz("test.tsv"));
    EXPECT_FALSE(ends_with_gz("file.txt"));
    EXPECT_FALSE(ends_with_gz(""));
    EXPECT_FALSE(ends_with_gz(".gz"));  // too short (size must be > 3)
    EXPECT_FALSE(ends_with_gz("gz"));
}


// ============================================================================
// Factory function
// ============================================================================

TEST(CreateOutputWriter, CreatesTSVWriter) {
    TempFile tmp(".tsv");
    auto writer = create_output_writer(tmp.path(), OutputFormat::TSV);
    ASSERT_NE(writer, nullptr);
    // Verify it's actually a TSVWriter by using it
    writer->write_header();
    writer->close();
    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("## ENSEMBL VARIANT EFFECT PREDICTOR"), std::string::npos);
}

TEST(CreateOutputWriter, CreatesJSONWriter) {
    TempFile tmp(".json");
    auto writer = create_output_writer(tmp.path(), OutputFormat::JSON);
    ASSERT_NE(writer, nullptr);
    writer->write_header();
    writer->write_footer();
    writer->close();
    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("["), std::string::npos);
    EXPECT_NE(content.find("]"), std::string::npos);
}

TEST(CreateOutputWriter, CreatesVCFWriter) {
    TempFile tmp(".vcf");
    auto writer = create_output_writer(tmp.path(), OutputFormat::VCF);
    ASSERT_NE(writer, nullptr);
    writer->write_header();
    writer->close();
    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("##fileformat=VCFv4.2"), std::string::npos);
}


// ============================================================================
// TSVWriter
// ============================================================================

class TSVWriterTest : public ::testing::Test {
protected:
    TempFile tmp{".tsv"};
};

TEST_F(TSVWriterTest, HeaderGeneration) {
    TSVWriter writer(tmp.path());
    writer.write_header({});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("## ENSEMBL VARIANT EFFECT PREDICTOR"), std::string::npos);
    EXPECT_NE(content.find("#Uploaded_variation\tLocation\tAllele\tGene\tFeature\t"
                           "Feature_type\tConsequence\tcDNA_position\tCDS_position\t"
                           "Protein_position\tAmino_acids\tCodons\tExisting_variation\tExtra"),
              std::string::npos);
}

TEST_F(TSVWriterTest, SkipHeader) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.write_header({});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_TRUE(content.empty());
}

TEST_F(TSVWriterTest, BasicAnnotation) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());

    // Uploaded_variation: should use rs ID when available
    EXPECT_NE(content.find("rs28934576\t"), std::string::npos);
    // Location
    EXPECT_NE(content.find("chr17:7675088\t"), std::string::npos);
    // Allele
    EXPECT_NE(content.find("\tT\t"), std::string::npos);
    // Gene
    EXPECT_NE(content.find("\tENSG00000141510\t"), std::string::npos);
    // Feature
    EXPECT_NE(content.find("\tENST00000269305\t"), std::string::npos);
    // Feature_type
    EXPECT_NE(content.find("\tTranscript\t"), std::string::npos);
    // Consequence
    EXPECT_NE(content.find("\tmissense_variant\t"), std::string::npos);
    // Amino acids
    EXPECT_NE(content.find("\tH/R\t"), std::string::npos);
    // Codons
    EXPECT_NE(content.find("\tcAt/cGt\t"), std::string::npos);
}

TEST_F(TSVWriterTest, UploadedVariationFallback) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    ann.vcf_id.clear();  // No VCF ID
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    // Should use CHR_POS_REF/ALT format
    EXPECT_NE(content.find("chr17_7675088_C/T\t"), std::string::npos);
}

TEST_F(TSVWriterTest, UploadedVariationDotId) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    ann.vcf_id = ".";  // VCF missing ID
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    // Dot ID should also fallback to CHR_POS_REF/ALT format
    EXPECT_NE(content.find("chr17_7675088_C/T\t"), std::string::npos);
}

TEST_F(TSVWriterTest, ExtraFieldOrdering) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    // Extra should contain IMPACT first
    size_t impact_pos = content.find("IMPACT=MODERATE");
    EXPECT_NE(impact_pos, std::string::npos);

    // STRAND should appear
    size_t strand_pos = content.find("STRAND=-1");
    EXPECT_NE(strand_pos, std::string::npos);

    // SYMBOL after STRAND
    size_t symbol_pos = content.find("SYMBOL=TP53");
    EXPECT_NE(symbol_pos, std::string::npos);
    EXPECT_GT(symbol_pos, strand_pos);

    // BIOTYPE after SYMBOL
    size_t biotype_pos = content.find("BIOTYPE=protein_coding");
    EXPECT_NE(biotype_pos, std::string::npos);
    EXPECT_GT(biotype_pos, symbol_pos);

    // CANONICAL=YES
    EXPECT_NE(content.find("CANONICAL=YES"), std::string::npos);

    // EXON
    EXPECT_NE(content.find("EXON=5/11"), std::string::npos);

    // HGVSc
    EXPECT_NE(content.find("HGVSc=ENST00000269305.9:c.535A>G"), std::string::npos);

    // HGVSp
    EXPECT_NE(content.find("HGVSp=ENSP00000269305.4:p.His179Arg"), std::string::npos);
}

TEST_F(TSVWriterTest, ExtraDistanceField) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    ann.distance = 4200;
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("DISTANCE=4200"), std::string::npos);
}

TEST_F(TSVWriterTest, ExtraNoDistanceWhenZero) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    ann.distance = 0;
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_EQ(content.find("DISTANCE="), std::string::npos);
}

TEST_F(TSVWriterTest, CustomColumnsInExtra) {
    TSVWriter writer(tmp.path());
    std::vector<std::string> custom_cols = {"SIFT", "PolyPhen", "gnomAD_AF"};
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["SIFT"] = "deleterious(0.01)";
    ann.custom_annotations["PolyPhen"] = "probably_damaging(0.99)";
    ann.custom_annotations["gnomAD_AF"] = "0.0001";
    writer.write_annotation(ann, custom_cols);
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("SIFT=deleterious(0.01)"), std::string::npos);
    EXPECT_NE(content.find("PolyPhen=probably_damaging(0.99)"), std::string::npos);
    EXPECT_NE(content.find("gnomAD_AF=0.0001"), std::string::npos);
}

TEST_F(TSVWriterTest, CustomColumnsDeduplicated) {
    // Custom columns that are already in phase1 (like BIOTYPE) should not be duplicated
    TSVWriter writer(tmp.path());
    std::vector<std::string> custom_cols = {"BIOTYPE", "CANONICAL", "IMPACT"};
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, custom_cols);
    writer.close();

    std::string content = read_file(tmp.path());
    // BIOTYPE should appear exactly once (via phase1 output)
    size_t first = content.find("BIOTYPE=protein_coding");
    size_t second = content.find("BIOTYPE=protein_coding", first + 1);
    EXPECT_NE(first, std::string::npos);
    EXPECT_EQ(second, std::string::npos);  // No second occurrence
}

TEST_F(TSVWriterTest, EmptyFieldsDash) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_empty_annotation();
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    // Multiple fields should show "-" for empty values
    // Amino acids, codons, existing variation, etc. should all be "-"
    // Count dashes as tab-separated fields
    std::istringstream iss(content);
    std::string line;
    std::getline(iss, line);
    // Fields: uploaded_var, location, allele, gene, feature, feature_type, consequence,
    // cdna_pos, cds_pos, protein_pos, amino_acids, codons, existing_variation, extra
    std::istringstream linestream(line);
    std::string field;
    std::vector<std::string> fields;
    while (std::getline(linestream, field, '\t')) {
        fields.push_back(field);
    }
    ASSERT_EQ(fields.size(), 14u);
    // cdna_position (7), cds_position (8), protein_position (9) should be "-"
    EXPECT_EQ(fields[7], "-");
    EXPECT_EQ(fields[8], "-");
    EXPECT_EQ(fields[9], "-");
    // amino_acids (10), codons (11), existing_variation (12) should be "-"
    EXPECT_EQ(fields[10], "-");
    EXPECT_EQ(fields[11], "-");
    EXPECT_EQ(fields[12], "-");
}

TEST_F(TSVWriterTest, MultiBaseLocationRange) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_multibase_annotation();
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    // Location should show range from display coords
    EXPECT_NE(content.find("chr17:7675081-7675083"), std::string::npos);
}

TEST_F(TSVWriterTest, MultiBasePositionRanges) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_multibase_annotation();
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    // cDNA position range with total
    EXPECT_NE(content.find("790-792/2591"), std::string::npos);
    // CDS position range with total
    EXPECT_NE(content.find("543-545/1182"), std::string::npos);
    // Protein position range with total (1182/3 = 394)
    EXPECT_NE(content.find("181-182/394"), std::string::npos);
}

TEST_F(TSVWriterTest, FlagsInExtra) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["FLAGS"] = "cds_start_NF";
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("FLAGS=cds_start_NF"), std::string::npos);
}

TEST_F(TSVWriterTest, IntronNumberInExtra) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    ann.exon_number = 0;    // clear exon
    ann.intron_number = 3;
    ann.total_introns = 10;
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("INTRON=3/10"), std::string::npos);
    EXPECT_EQ(content.find("EXON="), std::string::npos);
}

TEST_F(TSVWriterTest, HGVSOffsetInExtra) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["HGVS_OFFSET"] = "2";
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("HGVS_OFFSET=2"), std::string::npos);
}

TEST_F(TSVWriterTest, SymbolSourceAndHGNCInExtra) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["SYMBOL_SOURCE"] = "HGNC";
    ann.custom_annotations["HGNC_ID"] = "HGNC:11998";
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("SYMBOL_SOURCE=HGNC"), std::string::npos);
    EXPECT_NE(content.find("HGNC_ID=HGNC:11998"), std::string::npos);
}

TEST_F(TSVWriterTest, DisplayTermStyle) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_term_style("display");
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    // missense_variant display term is NON_SYNONYMOUS_CODING
    EXPECT_NE(content.find("NON_SYNONYMOUS_CODING"), std::string::npos);
    EXPECT_EQ(content.find("missense_variant"), std::string::npos);
}

TEST_F(TSVWriterTest, StatsTracking) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann1 = make_basic_annotation();
    VariantAnnotation ann2 = make_empty_annotation();
    writer.write_annotation(ann1, {});
    writer.write_annotation(ann2, {});
    writer.close();

    const auto& stats = writer.get_stats();
    EXPECT_EQ(stats.total_variants, 2);
    EXPECT_EQ(stats.annotated_variants, 1);
}

TEST_F(TSVWriterTest, WriteMultipleAnnotations) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    std::vector<VariantAnnotation> anns = {make_basic_annotation(), make_empty_annotation()};
    writer.write_annotations(anns, {});
    writer.close();

    std::string content = read_file(tmp.path());
    // Should contain two lines (each ending with \n)
    int line_count = 0;
    for (char c : content) {
        if (c == '\n') line_count++;
    }
    EXPECT_EQ(line_count, 2);
}

TEST_F(TSVWriterTest, SourceInExtra) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("SOURCE=Ensembl"), std::string::npos);
}

TEST_F(TSVWriterTest, HGVSgInExtra) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("HGVSg=chr17:g.7675088C>T"), std::string::npos);
}

TEST_F(TSVWriterTest, DisplayAltForDeletion) {
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_multibase_annotation();
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    // The allele column should show the display_alt ("-" for deletion)
    // Tab-separated: ...\t-\t...
    // After location and before gene
    EXPECT_NE(content.find("\t-\t"), std::string::npos);
}


// ============================================================================
// JSONWriter
// ============================================================================

class JSONWriterTest : public ::testing::Test {
protected:
    TempFile tmp{".json"};
};

TEST_F(JSONWriterTest, EmptyOutput) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_EQ(content, "[\n\n]\n");
}

TEST_F(JSONWriterTest, SkipHeader) {
    JSONWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // Should NOT start with "[" when skip_header is true
    EXPECT_TRUE(content.empty() || content[0] != '[');
}

TEST_F(JSONWriterTest, SingleAnnotation) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());

    // Check top-level structure
    EXPECT_NE(content.find("\"input\":"), std::string::npos);
    EXPECT_NE(content.find("\"assembly_name\": \"GRCh38\""), std::string::npos);
    EXPECT_NE(content.find("\"seq_region_name\": \"chr17\""), std::string::npos);
    EXPECT_NE(content.find("\"start\": 7675088"), std::string::npos);
    EXPECT_NE(content.find("\"end\": 7675088"), std::string::npos);
    EXPECT_NE(content.find("\"allele_string\": \"C/T\""), std::string::npos);
    EXPECT_NE(content.find("\"strand\": 1"), std::string::npos);
    EXPECT_NE(content.find("\"most_severe_consequence\": \"missense_variant\""), std::string::npos);
}

TEST_F(JSONWriterTest, TranscriptConsequenceFields) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());

    EXPECT_NE(content.find("\"transcript_consequences\""), std::string::npos);
    EXPECT_NE(content.find("\"gene_id\": \"ENSG00000141510\""), std::string::npos);
    EXPECT_NE(content.find("\"gene_symbol\": \"TP53\""), std::string::npos);
    EXPECT_NE(content.find("\"transcript_id\": \"ENST00000269305\""), std::string::npos);
    EXPECT_NE(content.find("\"biotype\": \"protein_coding\""), std::string::npos);
    EXPECT_NE(content.find("\"canonical\": 1"), std::string::npos);
    EXPECT_NE(content.find("\"variant_allele\": \"T\""), std::string::npos);
    EXPECT_NE(content.find("\"consequence_terms\": [\"missense_variant\"]"), std::string::npos);
    EXPECT_NE(content.find("\"impact\": \"MODERATE\""), std::string::npos);
    EXPECT_NE(content.find("\"strand\": -1"), std::string::npos);
}

TEST_F(JSONWriterTest, PositionFields) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());

    EXPECT_NE(content.find("\"cdna_start\": 782"), std::string::npos);
    EXPECT_NE(content.find("\"cdna_end\": 782"), std::string::npos);
    EXPECT_NE(content.find("\"cds_start\": 535"), std::string::npos);
    EXPECT_NE(content.find("\"cds_end\": 535"), std::string::npos);
    EXPECT_NE(content.find("\"protein_start\": 179"), std::string::npos);
    EXPECT_NE(content.find("\"protein_end\": 179"), std::string::npos);
}

TEST_F(JSONWriterTest, ExonIntronFields) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"exon\": \"5/11\""), std::string::npos);
}

TEST_F(JSONWriterTest, HGVSFields) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"hgvsc\": \"ENST00000269305.9:c.535A>G\""), std::string::npos);
    EXPECT_NE(content.find("\"hgvsp\": \"ENSP00000269305.4:p.His179Arg\""), std::string::npos);
    EXPECT_NE(content.find("\"hgvsg\": \"chr17:g.7675088C>T\""), std::string::npos);
}

TEST_F(JSONWriterTest, AminoAcidsAndCodons) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"amino_acids\": \"H/R\""), std::string::npos);
    EXPECT_NE(content.find("\"codons\": \"cAt/cGt\""), std::string::npos);
}

TEST_F(JSONWriterTest, ColocatedVariants) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"colocated_variants\""), std::string::npos);
    EXPECT_NE(content.find("\"id\": \"rs28934576\""), std::string::npos);
}

TEST_F(JSONWriterTest, MultipleColocatedVariants) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.existing_variation = "rs28934576,rs12345";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"id\": \"rs28934576\""), std::string::npos);
    EXPECT_NE(content.find("\"id\": \"rs12345\""), std::string::npos);
}

TEST_F(JSONWriterTest, NoColocatedVariantsWhenEmpty) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.existing_variation.clear();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_EQ(content.find("\"colocated_variants\""), std::string::npos);
}

TEST_F(JSONWriterTest, IntergenicConsequences) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_empty_annotation();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"intergenic_consequences\""), std::string::npos);
    EXPECT_NE(content.find("\"intergenic_variant\""), std::string::npos);
    EXPECT_EQ(content.find("\"transcript_consequences\""), std::string::npos);
}

TEST_F(JSONWriterTest, AssemblyNameCustom) {
    JSONWriter writer(tmp.path());
    writer.set_assembly_name("GRCh37");
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"assembly_name\": \"GRCh37\""), std::string::npos);
}

TEST_F(JSONWriterTest, MultipleVariantsSeparation) {
    JSONWriter writer(tmp.path());
    writer.write_header({});

    VariantAnnotation ann1 = make_basic_annotation();
    VariantAnnotation ann2 = make_basic_annotation();
    ann2.chromosome = "chr1";
    ann2.position = 100;
    ann2.ref_allele = "A";
    ann2.alt_allele = "G";

    writer.write_annotation(ann1, {});
    writer.write_annotation(ann2, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // Should have comma separating two variant objects
    // The structure is: [\n  {...},\n  {...}\n]\n
    EXPECT_NE(content.find("},\n"), std::string::npos);
    // Two seq_region_name entries
    size_t first_srn = content.find("\"seq_region_name\": \"chr17\"");
    size_t second_srn = content.find("\"seq_region_name\": \"chr1\"");
    EXPECT_NE(first_srn, std::string::npos);
    EXPECT_NE(second_srn, std::string::npos);
}

TEST_F(JSONWriterTest, GroupedTranscripts) {
    // Two annotations for the same variant (same chrom:pos:ref:alt) should be grouped
    JSONWriter writer(tmp.path());
    writer.write_header({});

    VariantAnnotation ann1 = make_basic_annotation();
    VariantAnnotation ann2 = make_basic_annotation();
    ann2.transcript_id = "ENST00000445888";
    ann2.is_canonical = false;

    writer.write_annotation(ann1, {});
    writer.write_annotation(ann2, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // Both transcripts should be in the same transcript_consequences array
    EXPECT_NE(content.find("ENST00000269305"), std::string::npos);
    EXPECT_NE(content.find("ENST00000445888"), std::string::npos);
    // Verify grouping: only one variant-level object
    // (both transcripts share the same "seq_region_name" block)
    size_t first_srn = content.find("\"seq_region_name\": \"chr17\"");
    EXPECT_NE(first_srn, std::string::npos);
    // A second occurrence is expected inside colocated_variants, but verify
    // there is no third occurrence (which would indicate a second variant object)
    size_t second_srn = content.find("\"seq_region_name\": \"chr17\"", first_srn + 1);
    if (second_srn != std::string::npos) {
        size_t third_srn = content.find("\"seq_region_name\": \"chr17\"", second_srn + 1);
        // At most two occurrences (variant-level + colocated_variants)
        EXPECT_EQ(third_srn, std::string::npos);
    }
}

TEST_F(JSONWriterTest, SIFTField) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["SIFT"] = "deleterious(0.01)";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"sift_prediction\": \"deleterious\""), std::string::npos);
    EXPECT_NE(content.find("\"sift_score\": 0.01"), std::string::npos);
}

TEST_F(JSONWriterTest, PolyPhenField) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["PolyPhen"] = "probably_damaging(0.997)";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"polyphen_prediction\": \"probably_damaging\""), std::string::npos);
    EXPECT_NE(content.find("\"polyphen_score\": 0.997"), std::string::npos);
}

TEST_F(JSONWriterTest, FlagsAsArray) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["FLAGS"] = "cds_start_NF,cds_end_NF";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"flags\": [\"cds_start_NF\", \"cds_end_NF\"]"), std::string::npos);
}

TEST_F(JSONWriterTest, DomainsField) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["pfam:domain_id"] = "PF00870";
    ann.custom_annotations["interpro:domain_id"] = "IPR011615";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"domains\""), std::string::npos);
    EXPECT_NE(content.find("\"db\": \"Pfam\""), std::string::npos);
    EXPECT_NE(content.find("\"name\": \"PF00870\""), std::string::npos);
    EXPECT_NE(content.find("\"db\": \"Interpro\""), std::string::npos);
    EXPECT_NE(content.find("\"name\": \"IPR011615\""), std::string::npos);
}

TEST_F(JSONWriterTest, MANEFields) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["MANE_SELECT"] = "NM_000546.6";
    ann.custom_annotations["MANE_PLUS_CLINICAL"] = "NM_001276695.2";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"mane_select\": \"NM_000546.6\""), std::string::npos);
    EXPECT_NE(content.find("\"mane_plus_clinical\": \"NM_001276695.2\""), std::string::npos);
}

TEST_F(JSONWriterTest, TSLField) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["TSL"] = "1";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // TSL is output as numeric
    EXPECT_NE(content.find("\"tsl\": 1"), std::string::npos);
}

TEST_F(JSONWriterTest, ProteinIdField) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["ENSP"] = "ENSP00000269305";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"protein_id\": \"ENSP00000269305\""), std::string::npos);
}

TEST_F(JSONWriterTest, CustomAnnotationNumericValue) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["my_score"] = "0.95";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // Numeric values should be output as numbers (no quotes)
    EXPECT_NE(content.find("\"my_score\": 0.95"), std::string::npos);
}

TEST_F(JSONWriterTest, CustomAnnotationStringValue) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["my_label"] = "important";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // String values should be quoted
    EXPECT_NE(content.find("\"my_label\": \"important\""), std::string::npos);
}

TEST_F(JSONWriterTest, JSONEscaping) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.input_variant = "has \"quotes\" and \\backslash";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("has \\\"quotes\\\" and \\\\backslash"), std::string::npos);
}

TEST_F(JSONWriterTest, IdWithVcfId) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.vcf_id = "rs28934576";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"id\": \"rs28934576\""), std::string::npos);
}

TEST_F(JSONWriterTest, IdFallbackWithoutVcfId) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.vcf_id.clear();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // Should use CHR_POS_REF/ALT format
    EXPECT_NE(content.find("\"id\": \"chr17_7675088_C/T\""), std::string::npos);
}

TEST_F(JSONWriterTest, EmptyAllelesDashInId) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.vcf_id.clear();
    ann.ref_allele = "";
    ann.alt_allele = "A";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // Empty ref allele should display as "-"
    EXPECT_NE(content.find("-/A"), std::string::npos);
}

TEST_F(JSONWriterTest, VariantClassField) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["VARIANT_CLASS"] = "SNV";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"variant_class\": \"SNV\""), std::string::npos);
}

TEST_F(JSONWriterTest, ClinSigInColocated) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["CLIN_SIG"] = "pathogenic,likely_pathogenic";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"clin_sig\""), std::string::npos);
    EXPECT_NE(content.find("\"pathogenic\""), std::string::npos);
    EXPECT_NE(content.find("\"likely_pathogenic\""), std::string::npos);
}

TEST_F(JSONWriterTest, DisplayTermStyle) {
    JSONWriter writer(tmp.path());
    writer.set_term_style("display");
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"most_severe_consequence\": \"NON_SYNONYMOUS_CODING\""), std::string::npos);
    EXPECT_NE(content.find("\"NON_SYNONYMOUS_CODING\""), std::string::npos);
}

TEST_F(JSONWriterTest, RegulatoryConsequences) {
    JSONWriter writer(tmp.path());
    writer.write_header({});

    VariantAnnotation ann;
    ann.chromosome = "chr1";
    ann.position = 1000;
    ann.ref_allele = "A";
    ann.alt_allele = "G";
    ann.consequences = {ConsequenceType::REGULATORY_REGION_VARIANT};
    ann.impact = Impact::MODIFIER;
    ann.gene_id = "ENSR00000012345";
    ann.transcript_id = "";  // Empty transcript_id signals regulatory
    ann.custom_annotations["regulatory:feature_type"] = "promoter";
    ann.display_alt = "G";

    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"regulatory_feature_consequences\""), std::string::npos);
    EXPECT_NE(content.find("\"biotype\": \"promoter\""), std::string::npos);
    EXPECT_NE(content.find("\"regulatory_feature_id\": \"ENSR00000012345\""), std::string::npos);
}

TEST_F(JSONWriterTest, SourceField) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"source\": \"Ensembl\""), std::string::npos);
}

TEST_F(JSONWriterTest, NonCanonicalOmitsCanonicalField) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.is_canonical = false;
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_EQ(content.find("\"canonical\""), std::string::npos);
}

TEST_F(JSONWriterTest, DistanceField) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.distance = 3500;
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"distance\": 3500"), std::string::npos);
}

TEST_F(JSONWriterTest, NoDistanceWhenZero) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    ann.distance = 0;
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // "distance" should not appear in transcript_consequences when 0
    // (but it may appear in other contexts - check carefully)
    std::string tc_section;
    auto tc_start = content.find("\"transcript_consequences\"");
    if (tc_start != std::string::npos) {
        tc_section = content.substr(tc_start);
    }
    EXPECT_EQ(tc_section.find("\"distance\":"), std::string::npos);
}

TEST_F(JSONWriterTest, StatsTracking) {
    JSONWriter writer(tmp.path());
    writer.write_header({});
    writer.write_annotation(make_basic_annotation(), {});
    writer.write_annotation(make_empty_annotation(), {});
    writer.write_footer();
    writer.close();

    const auto& stats = writer.get_stats();
    EXPECT_EQ(stats.total_variants, 2);
    EXPECT_EQ(stats.annotated_variants, 1);
}


// ============================================================================
// VCFWriter
// ============================================================================

class VCFWriterTest : public ::testing::Test {
protected:
    TempFile tmp{".vcf"};
};

TEST_F(VCFWriterTest, DefaultHeader) {
    VCFWriter writer(tmp.path());
    writer.write_header({});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("##fileformat=VCFv4.2"), std::string::npos);
    EXPECT_NE(content.find("##INFO=<ID=CSQ,"), std::string::npos);
    EXPECT_NE(content.find("##FILTER=<ID=PASS,"), std::string::npos);
    EXPECT_NE(content.find("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), std::string::npos);
}

TEST_F(VCFWriterTest, CSQFormatInHeader) {
    VCFWriter writer(tmp.path());
    writer.write_header({});
    writer.close();

    std::string content = read_file(tmp.path());
    // Default CSQ format should contain standard fields
    EXPECT_NE(content.find("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|"
                           "EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|"
                           "Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS"),
              std::string::npos);
}

TEST_F(VCFWriterTest, CSQFormatWithCustomColumns) {
    VCFWriter writer(tmp.path());
    std::vector<std::string> custom_cols = {"SIFT", "PolyPhen"};
    writer.write_header(custom_cols);
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("|SIFT|PolyPhen"), std::string::npos);
}

TEST_F(VCFWriterTest, CustomInfoFieldName) {
    VCFWriter writer(tmp.path());
    writer.set_info_field_name("ANN");
    writer.write_header({});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("##INFO=<ID=ANN,"), std::string::npos);
    EXPECT_EQ(content.find("##INFO=<ID=CSQ,"), std::string::npos);
}

TEST_F(VCFWriterTest, CustomFieldOrder) {
    VCFWriter writer(tmp.path());
    writer.set_field_order({"Allele", "Consequence", "Gene"});
    writer.write_header({});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("Format: Allele|Consequence|Gene"), std::string::npos);
}

TEST_F(VCFWriterTest, SkipHeader) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.write_header({});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_TRUE(content.empty());
}

TEST_F(VCFWriterTest, PassthroughHeaders) {
    VCFWriter writer(tmp.path());
    writer.add_passthrough_header("##fileformat=VCFv4.1");
    writer.add_passthrough_header("##contig=<ID=chr17,length=83257441>");
    writer.add_passthrough_header("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">");
    writer.set_column_header("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1");
    writer.write_header({});
    writer.close();

    std::string content = read_file(tmp.path());
    // Passthrough headers preserved
    EXPECT_NE(content.find("##fileformat=VCFv4.1"), std::string::npos);
    EXPECT_NE(content.find("##contig=<ID=chr17"), std::string::npos);
    EXPECT_NE(content.find("##INFO=<ID=DP,"), std::string::npos);
    // CSQ header added
    EXPECT_NE(content.find("##INFO=<ID=CSQ,"), std::string::npos);
    // No duplicate fileformat or FILTER lines from default
    EXPECT_EQ(content.find("##fileformat=VCFv4.2"), std::string::npos);
    EXPECT_EQ(content.find("##FILTER=<ID=PASS"), std::string::npos);
    // Original column header preserved with sample columns
    EXPECT_NE(content.find("SAMPLE1"), std::string::npos);
}

TEST_F(VCFWriterTest, PassthroughExcludesExistingCSQ) {
    VCFWriter writer(tmp.path());
    writer.add_passthrough_header("##fileformat=VCFv4.2");
    writer.add_passthrough_header("##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Old CSQ\">");
    writer.add_passthrough_header("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">");
    writer.write_header({});
    writer.close();

    std::string content = read_file(tmp.path());
    // Old CSQ header should be filtered out
    EXPECT_EQ(content.find("Old CSQ"), std::string::npos);
    // New CSQ header should be present
    EXPECT_NE(content.find("##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations"), std::string::npos);
    // Other headers preserved
    EXPECT_NE(content.find("##INFO=<ID=DP,"), std::string::npos);
}

TEST_F(VCFWriterTest, BasicCSQOutput) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Should contain CSQ= in INFO column
    EXPECT_NE(content.find("CSQ="), std::string::npos);
    // Check essential CSQ fields
    EXPECT_NE(content.find("missense_variant"), std::string::npos);
    EXPECT_NE(content.find("MODERATE"), std::string::npos);
    EXPECT_NE(content.find("TP53"), std::string::npos);
    EXPECT_NE(content.find("ENSG00000141510"), std::string::npos);
    EXPECT_NE(content.find("ENST00000269305"), std::string::npos);
    EXPECT_NE(content.find("protein_coding"), std::string::npos);
}

TEST_F(VCFWriterTest, CSQFieldSeparation) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "Consequence", "IMPACT"});
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // With custom field order, CSQ should be: T|missense_variant|MODERATE
    EXPECT_NE(content.find("T|missense_variant|MODERATE"), std::string::npos);
}

TEST_F(VCFWriterTest, VCFEscapePipe) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "Existing_variation"});
    VariantAnnotation ann = make_basic_annotation();
    ann.existing_variation = "rs123|rs456";  // Pipe in existing variation
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Pipe should be escaped to &
    EXPECT_NE(content.find("rs123&rs456"), std::string::npos);
    // No raw pipe in the CSQ value itself (besides field separators)
}

TEST_F(VCFWriterTest, VCFEscapeComma) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "Existing_variation"});
    VariantAnnotation ann = make_basic_annotation();
    ann.existing_variation = "rs123,rs456";  // Comma in existing variation
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Comma should be escaped to &
    EXPECT_NE(content.find("rs123&rs456"), std::string::npos);
}

TEST_F(VCFWriterTest, VCFEscapeSemicolon) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "SYMBOL"});
    VariantAnnotation ann = make_basic_annotation();
    ann.gene_symbol = "gene;name";
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("gene%3Bname"), std::string::npos);
}

TEST_F(VCFWriterTest, VCFEscapeEquals) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "SYMBOL"});
    VariantAnnotation ann = make_basic_annotation();
    ann.gene_symbol = "gene=name";
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("gene%3Dname"), std::string::npos);
}

TEST_F(VCFWriterTest, VCFEscapePercent) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "SYMBOL"});
    VariantAnnotation ann = make_basic_annotation();
    ann.gene_symbol = "100%gene";
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("100%25gene"), std::string::npos);
}

TEST_F(VCFWriterTest, VCFEscapeSpace) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "SYMBOL"});
    VariantAnnotation ann = make_basic_annotation();
    ann.gene_symbol = "gene name";
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("gene_name"), std::string::npos);
}

TEST_F(VCFWriterTest, NoEscapeMode) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_no_escape(true);
    writer.set_field_order({"Allele", "SYMBOL"});
    VariantAnnotation ann = make_basic_annotation();
    ann.gene_symbol = "gene;name";
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // No escaping should preserve raw characters
    EXPECT_NE(content.find("gene;name"), std::string::npos);
}

TEST_F(VCFWriterTest, MultiAlleleCSQMerging) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann1 = make_basic_annotation();
    VariantAnnotation ann2 = make_basic_annotation();
    ann2.transcript_id = "ENST00000445888";
    ann2.consequences = {ConsequenceType::SYNONYMOUS_VARIANT};
    ann2.impact = Impact::LOW;

    writer.write_annotation(ann1, {});
    writer.write_annotation(ann2, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Both CSQ entries should be comma-separated on same line
    EXPECT_NE(content.find("missense_variant"), std::string::npos);
    EXPECT_NE(content.find("synonymous_variant"), std::string::npos);
    // Only one line (one VCF data row)
    int line_count = 0;
    for (char c : content) {
        if (c == '\n') line_count++;
    }
    EXPECT_EQ(line_count, 1);
}

TEST_F(VCFWriterTest, VCFFieldPassthrough) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.vcf_id = "rs28934576";
    ann.vcf_qual = "99";
    ann.vcf_filter = "PASS";
    ann.vcf_info = "DP=100;AF=0.5";
    ann.vcf_sample_columns = "GT:DP\t0/1:50";

    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("rs28934576"), std::string::npos);
    EXPECT_NE(content.find("\t99\t"), std::string::npos);
    EXPECT_NE(content.find("\tPASS\t"), std::string::npos);
    EXPECT_NE(content.find("DP=100;AF=0.5"), std::string::npos);
    EXPECT_NE(content.find("GT:DP\t0/1:50"), std::string::npos);
}

TEST_F(VCFWriterTest, VCFInfoPreservation) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.vcf_info = "DP=100;AF=0.5";
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Original INFO fields should be preserved with CSQ appended
    EXPECT_NE(content.find("DP=100;AF=0.5;CSQ="), std::string::npos);
}

TEST_F(VCFWriterTest, VCFInfoReplacesExistingCSQ) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.vcf_info = "DP=100;CSQ=old_data;AF=0.5";
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Old CSQ should be removed, new one added
    EXPECT_EQ(content.find("old_data"), std::string::npos);
    EXPECT_NE(content.find("DP=100;AF=0.5;CSQ="), std::string::npos);
}

TEST_F(VCFWriterTest, KeepCSQMode) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_keep_csq(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["EXISTING_CSQ"] = "old_allele|old_consequence|HIGH";
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Existing CSQ should be prepended
    size_t csq_pos = content.find("CSQ=");
    ASSERT_NE(csq_pos, std::string::npos);
    std::string after_csq = content.substr(csq_pos + 4);
    // Old CSQ should come before new CSQ
    size_t old_pos = after_csq.find("old_allele");
    size_t new_pos = after_csq.find("missense_variant");
    EXPECT_NE(old_pos, std::string::npos);
    EXPECT_NE(new_pos, std::string::npos);
    EXPECT_LT(old_pos, new_pos);
}

TEST_F(VCFWriterTest, DefaultCSQFormat) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Default format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|
    // EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|
    // Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS

    // Extract CSQ value
    size_t csq_start = content.find("CSQ=");
    ASSERT_NE(csq_start, std::string::npos);
    std::string csq = content.substr(csq_start + 4);
    size_t csq_end = csq.find_first_of("\t\n");
    if (csq_end != std::string::npos) csq = csq.substr(0, csq_end);

    // Split by | to verify field count
    int pipe_count = 0;
    for (char c : csq) {
        if (c == '|') pipe_count++;
    }
    // 21 fields = 20 pipes (Allele through FLAGS)
    EXPECT_GE(pipe_count, 20);
}

TEST_F(VCFWriterTest, ExonIntronInCSQ) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Exon should appear as "5/11"
    EXPECT_NE(content.find("5/11"), std::string::npos);
}

TEST_F(VCFWriterTest, CanonicalInCSQ) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "CANONICAL"});

    VariantAnnotation ann = make_basic_annotation();
    ann.is_canonical = true;
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("T|YES"), std::string::npos);
}

TEST_F(VCFWriterTest, CanonicalEmptyWhenFalse) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "CANONICAL"});

    VariantAnnotation ann = make_basic_annotation();
    ann.is_canonical = false;
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // When not canonical, the field should be empty (T|)
    EXPECT_NE(content.find("T|"), std::string::npos);
    EXPECT_EQ(content.find("YES"), std::string::npos);
}

TEST_F(VCFWriterTest, StrandInCSQ) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "STRAND"});

    VariantAnnotation ann = make_basic_annotation();
    ann.strand = '-';
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("T|-1"), std::string::npos);
}

TEST_F(VCFWriterTest, DistanceInCSQ) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "DISTANCE"});

    VariantAnnotation ann = make_basic_annotation();
    ann.distance = 3000;
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("T|3000"), std::string::npos);
}

TEST_F(VCFWriterTest, DistanceEmptyWhenZero) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "DISTANCE"});

    VariantAnnotation ann = make_basic_annotation();
    ann.distance = 0;
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Field should be empty when distance is 0
    size_t csq_pos = content.find("CSQ=");
    ASSERT_NE(csq_pos, std::string::npos);
    std::string csq_val = content.substr(csq_pos + 4);
    size_t end = csq_val.find_first_of("\t\n");
    if (end != std::string::npos) csq_val = csq_val.substr(0, end);
    // Should be "T|" (empty after pipe)
    EXPECT_EQ(csq_val, "T|");
}

TEST_F(VCFWriterTest, CustomAnnotationsInCSQ) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    std::vector<std::string> custom_cols = {"SIFT", "PolyPhen"};
    writer.write_header(custom_cols);

    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["SIFT"] = "deleterious(0.01)";
    ann.custom_annotations["PolyPhen"] = "probably_damaging(0.99)";
    writer.write_annotation(ann, custom_cols);
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("deleterious(0.01)"), std::string::npos);
    EXPECT_NE(content.find("probably_damaging(0.99)"), std::string::npos);
}

TEST_F(VCFWriterTest, FlagsInCSQ) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["FLAGS"] = "cds_start_NF";
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("cds_start_NF"), std::string::npos);
}

TEST_F(VCFWriterTest, VcfRefAltPreserved) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.vcf_ref = "CGT";
    ann.vcf_alt = "CAGT,CA";
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\tCGT\t"), std::string::npos);
    EXPECT_NE(content.find("\tCAGT,CA\t"), std::string::npos);
}

TEST_F(VCFWriterTest, EmptyRefAltDash) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.ref_allele = "";
    ann.alt_allele = "A";
    ann.vcf_ref.clear();
    ann.vcf_alt.clear();
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Empty ref should be displayed as "-"
    EXPECT_NE(content.find("\t-\t"), std::string::npos);
}

TEST_F(VCFWriterTest, FlushOnWriteFooter) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    // Don't call flush_variant() - write_footer() should flush
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("CSQ="), std::string::npos);
}

TEST_F(VCFWriterTest, FlushOnClose) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    // Don't call flush_variant() or write_footer() - close() should flush
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("CSQ="), std::string::npos);
}

TEST_F(VCFWriterTest, MissingVCFFieldsDefaultToDot) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.vcf_id.clear();
    ann.vcf_qual.clear();
    ann.vcf_filter.clear();
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Tab-separated: CHROM POS ID REF ALT QUAL FILTER INFO
    // ID, QUAL, FILTER should all be "."
    std::istringstream iss(content);
    std::string line;
    std::getline(iss, line);
    std::istringstream linestream(line);
    std::string field;
    std::vector<std::string> fields;
    while (std::getline(linestream, field, '\t')) {
        fields.push_back(field);
    }
    ASSERT_GE(fields.size(), 8u);
    EXPECT_EQ(fields[2], ".");      // ID
    EXPECT_EQ(fields[5], ".");      // QUAL
    EXPECT_EQ(fields[6], ".");      // FILTER
}

TEST_F(VCFWriterTest, StatsTracking) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    VariantAnnotation ann1 = make_basic_annotation();
    VariantAnnotation ann2 = make_empty_annotation();
    writer.write_annotation(ann1, {});
    writer.write_annotation(ann2, {});
    writer.flush_variant();
    writer.close();

    const auto& stats = writer.get_stats();
    EXPECT_EQ(stats.total_variants, 2);
    EXPECT_EQ(stats.annotated_variants, 1);
}

TEST_F(VCFWriterTest, PositionRangesInCSQ) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "cDNA_position", "CDS_position", "Protein_position"});

    VariantAnnotation ann = make_multibase_annotation();
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // cDNA position with range and total
    EXPECT_NE(content.find("790-792/2591"), std::string::npos);
    // CDS position with range and total
    EXPECT_NE(content.find("543-545/1182"), std::string::npos);
    // Protein position with range and total (1182/3 = 394)
    EXPECT_NE(content.find("181-182/394"), std::string::npos);
}

TEST_F(VCFWriterTest, EmptyPositionsInCSQ) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "cDNA_position", "CDS_position", "Protein_position"});

    VariantAnnotation ann = make_empty_annotation();
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // All position fields should be empty (not "-" as in TSV)
    size_t csq_pos = content.find("CSQ=");
    ASSERT_NE(csq_pos, std::string::npos);
    std::string csq_val = content.substr(csq_pos + 4);
    size_t end = csq_val.find_first_of("\t\n");
    if (end != std::string::npos) csq_val = csq_val.substr(0, end);
    // Should be "G|||" (allele, then three empty position fields)
    EXPECT_EQ(csq_val, "G|||");
}

TEST_F(VCFWriterTest, SourceFieldInCSQ) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "SOURCE"});

    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("T|Ensembl"), std::string::npos);
}

TEST_F(VCFWriterTest, HGVSFieldsInCSQ) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "HGVSc", "HGVSp", "HGVSg"});

    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("ENST00000269305.9:c.535A>G"), std::string::npos);
    EXPECT_NE(content.find("ENSP00000269305.4:p.His179Arg"), std::string::npos);
    EXPECT_NE(content.find("chr17:g.7675088C>T"), std::string::npos);
}

TEST_F(VCFWriterTest, CustomAnnotationFallthrough) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "MY_CUSTOM_FIELD"});

    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["MY_CUSTOM_FIELD"] = "custom_value";
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("T|custom_value"), std::string::npos);
}

TEST_F(VCFWriterTest, DisplayAltUsedForAllele) {
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele"});

    VariantAnnotation ann = make_multibase_annotation();
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Allele field should use display_alt ("-" for deletion)
    EXPECT_NE(content.find("CSQ=-"), std::string::npos);
}


// ============================================================================
// Edge cases
// ============================================================================

TEST(OutputEdgeCases, MultipleConsequencesJoinedWithAmpersand) {
    TempFile tmp(".tsv");
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.consequences = {ConsequenceType::MISSENSE_VARIANT, ConsequenceType::SPLICE_REGION_VARIANT};
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("missense_variant&splice_region_variant"), std::string::npos);
}

TEST(OutputEdgeCases, SpecialCharsInInputVariant) {
    TempFile tmp(".json");
    JSONWriter writer(tmp.path());
    writer.write_header({});

    VariantAnnotation ann = make_basic_annotation();
    ann.input_variant = "chr17\t7675088\t.\tC\tT\t.\t.\t.";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // Tab should be escaped in JSON
    EXPECT_NE(content.find("\\t"), std::string::npos);
}

TEST(OutputEdgeCases, EmptyExistingVariation) {
    TempFile tmp(".tsv");
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.existing_variation.clear();
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    // Should show "-" for empty existing_variation
    // Find the 13th field (0-indexed 12th)
    std::istringstream iss(content);
    std::string line;
    std::getline(iss, line);
    std::istringstream linestream(line);
    std::string field;
    std::vector<std::string> fields;
    while (std::getline(linestream, field, '\t')) {
        fields.push_back(field);
    }
    ASSERT_GE(fields.size(), 13u);
    EXPECT_EQ(fields[12], "-");
}

TEST(OutputEdgeCases, NoConsequencesAnnotation) {
    TempFile tmp(".tsv");
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.consequences.clear();
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    // Should have an empty consequence field (or at least not crash)
    EXPECT_FALSE(content.empty());
}

TEST(OutputEdgeCases, VCFWriterFlushEmptyDoesNothing) {
    TempFile tmp(".vcf");
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.flush_variant();  // Nothing buffered
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_TRUE(content.empty());
}

TEST(OutputEdgeCases, JSONEmptyAlleleStringDash) {
    TempFile tmp(".json");
    JSONWriter writer(tmp.path());
    writer.write_header({});

    VariantAnnotation ann = make_basic_annotation();
    ann.ref_allele = "";
    ann.alt_allele = "A";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // allele_string should show "-/A" for empty ref
    EXPECT_NE(content.find("\"allele_string\": \"-/A\""), std::string::npos);
}

TEST(OutputEdgeCases, TSVWriterToFile) {
    TempFile tmp(".tsv");
    TSVWriter writer(tmp.path());
    writer.write_header({});
    VariantAnnotation ann = make_basic_annotation();
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("## ENSEMBL VARIANT EFFECT PREDICTOR"), std::string::npos);
    EXPECT_NE(content.find("missense_variant"), std::string::npos);
}

TEST(OutputEdgeCases, JSONWriterEndPosCalculation) {
    TempFile tmp(".json");
    JSONWriter writer(tmp.path());
    writer.write_header({});

    VariantAnnotation ann = make_basic_annotation();
    // Multi-base ref: end should be position + ref_allele.size() - 1
    ann.ref_allele = "ATG";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // end = 7675088 + 3 - 1 = 7675090
    EXPECT_NE(content.find("\"end\": 7675090"), std::string::npos);
}

TEST(OutputEdgeCases, JSONMostSevereFromMultipleTranscripts) {
    TempFile tmp(".json");
    JSONWriter writer(tmp.path());
    writer.write_header({});

    // Two transcripts for same variant with different consequences
    VariantAnnotation ann1 = make_basic_annotation();
    ann1.consequences = {ConsequenceType::SYNONYMOUS_VARIANT};

    VariantAnnotation ann2 = make_basic_annotation();
    ann2.transcript_id = "ENST00000445888";
    ann2.consequences = {ConsequenceType::STOP_GAINED};

    writer.write_annotation(ann1, {});
    writer.write_annotation(ann2, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // stop_gained is more severe than synonymous_variant
    EXPECT_NE(content.find("\"most_severe_consequence\": \"stop_gained\""), std::string::npos);
}

TEST(OutputEdgeCases, VCFWriterStrandPlusStrand) {
    TempFile tmp(".vcf");
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "STRAND"});

    VariantAnnotation ann = make_basic_annotation();
    ann.strand = '+';
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("T|1"), std::string::npos);
}

TEST(OutputEdgeCases, VCFWriterStrandNull) {
    TempFile tmp(".vcf");
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);
    writer.set_field_order({"Allele", "STRAND"});

    VariantAnnotation ann = make_basic_annotation();
    ann.strand = '\0';
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Strand should be empty when null
    size_t csq_pos = content.find("CSQ=");
    ASSERT_NE(csq_pos, std::string::npos);
    std::string csq_val = content.substr(csq_pos + 4);
    size_t end = csq_val.find_first_of("\t\n");
    if (end != std::string::npos) csq_val = csq_val.substr(0, end);
    EXPECT_EQ(csq_val, "T|");
}

TEST(OutputEdgeCases, TSVNoStrandWhenNull) {
    TempFile tmp(".tsv");
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.strand = '\0';
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_EQ(content.find("STRAND="), std::string::npos);
}

TEST(OutputEdgeCases, TSVNonCanonicalOmitsCanonical) {
    TempFile tmp(".tsv");
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.is_canonical = false;
    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_EQ(content.find("CANONICAL="), std::string::npos);
}

TEST(OutputEdgeCases, JSONIntronField) {
    TempFile tmp(".json");
    JSONWriter writer(tmp.path());
    writer.write_header({});

    VariantAnnotation ann = make_basic_annotation();
    ann.exon_number = 0;
    ann.intron_number = 4;
    ann.total_introns = 10;
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"intron\": \"4/10\""), std::string::npos);
    EXPECT_EQ(content.find("\"exon\""), std::string::npos);
}

TEST(OutputEdgeCases, JSONNoPositionsForIntergenic) {
    TempFile tmp(".json");
    JSONWriter writer(tmp.path());
    writer.write_header({});

    VariantAnnotation ann = make_empty_annotation();
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // Should NOT have cdna_start, cds_start, protein_start for intergenic
    EXPECT_EQ(content.find("\"cdna_start\""), std::string::npos);
    EXPECT_EQ(content.find("\"cds_start\""), std::string::npos);
    EXPECT_EQ(content.find("\"protein_start\""), std::string::npos);
}

TEST(OutputEdgeCases, ExtraFieldDashWhenAllEmpty) {
    TempFile tmp(".tsv");
    TSVWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann;
    ann.chromosome = "chr1";
    ann.position = 100;
    ann.ref_allele = "A";
    ann.alt_allele = "G";
    ann.consequences.clear();  // No consequences
    ann.impact = Impact::MODIFIER;
    ann.feature_type = "Transcript";
    // All optional fields empty, strand null, no canonical, no distance

    writer.write_annotation(ann, {});
    writer.close();

    std::string content = read_file(tmp.path());
    // Extra field should at minimum have IMPACT
    EXPECT_NE(content.find("IMPACT=MODIFIER"), std::string::npos);
}

TEST(OutputEdgeCases, JSONSkippedCustomFields) {
    TempFile tmp(".json");
    JSONWriter writer(tmp.path());
    writer.write_header({});

    VariantAnnotation ann = make_basic_annotation();
    // These should be skipped in custom annotations output (already output elsewhere)
    ann.custom_annotations["CANONICAL"] = "YES";
    ann.custom_annotations["BIOTYPE"] = "protein_coding";
    ann.custom_annotations["STRAND"] = "-1";
    ann.custom_annotations["_consequences"] = "missense_variant";
    ann.custom_annotations["VARIANT_CLASS"] = "SNV";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    // These should NOT appear as lowercased custom keys in transcript_consequences
    // (but may appear as their proper explicit field names)
    // Check that "canonical" appears exactly as the integer field, not as string
    // Find transcript_consequences section
    auto tc_pos = content.find("\"transcript_consequences\"");
    ASSERT_NE(tc_pos, std::string::npos);
    std::string tc_section = content.substr(tc_pos);

    // _consequences should not leak through
    EXPECT_EQ(tc_section.find("\"_consequences\""), std::string::npos);
}

TEST(OutputEdgeCases, JSONSymbolSourceAndHGNC) {
    TempFile tmp(".json");
    JSONWriter writer(tmp.path());
    writer.write_header({});

    VariantAnnotation ann = make_basic_annotation();
    ann.custom_annotations["SYMBOL_SOURCE"] = "HGNC";
    ann.custom_annotations["HGNC_ID"] = "HGNC:11998";
    writer.write_annotation(ann, {});
    writer.write_footer();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("\"gene_symbol_source\": \"HGNC\""), std::string::npos);
    EXPECT_NE(content.find("\"hgnc_id\": \"HGNC:11998\""), std::string::npos);
}

TEST(OutputEdgeCases, VCFDotInfoField) {
    TempFile tmp(".vcf");
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.vcf_info = ".";
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    // Dot INFO should result in just CSQ= (not ".;CSQ=")
    EXPECT_NE(content.find("CSQ="), std::string::npos);
    EXPECT_EQ(content.find(".;CSQ="), std::string::npos);
}

TEST(OutputEdgeCases, VCFEmptyInfoField) {
    TempFile tmp(".vcf");
    VCFWriter writer(tmp.path());
    writer.set_skip_header(true);

    VariantAnnotation ann = make_basic_annotation();
    ann.vcf_info.clear();
    writer.write_annotation(ann, {});
    writer.flush_variant();
    writer.close();

    std::string content = read_file(tmp.path());
    EXPECT_NE(content.find("CSQ="), std::string::npos);
}
