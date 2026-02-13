/**
 * VEP Variant Annotator - Pure C++ Local Implementation
 *
 * A C++ library to annotate genetic variants using local annotation files.
 * No external API calls - everything runs locally.
 *
 * Required data files:
 * - GTF/GFF3 annotation file (gene/transcript definitions)
 * - FASTA reference genome file
 * - Optional: ClinVar VCF, gnomAD VCF for additional annotations
 */

#ifndef VEP_ANNOTATOR_HPP
#define VEP_ANNOTATOR_HPP

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <optional>
#include <memory>
#include <fstream>
#include <functional>
#include <shared_mutex>

namespace vep {

// Forward declarations
class TranscriptDatabase;
class ReferenceGenome;
class VCFAnnotationDatabase;
class AnnotationSource;
class AnnotationSourceManager;

/**
 * Sequence Ontology consequence types (SO terms)
 */
enum class ConsequenceType {
    // High impact
    TRANSCRIPT_ABLATION,
    SPLICE_ACCEPTOR_VARIANT,
    SPLICE_DONOR_VARIANT,
    STOP_GAINED,
    FRAMESHIFT_VARIANT,
    STOP_LOST,
    START_LOST,
    TRANSCRIPT_AMPLIFICATION,
    FEATURE_ELONGATION,
    FEATURE_TRUNCATION,

    // Moderate impact
    INFRAME_INSERTION,
    INFRAME_DELETION,
    MISSENSE_VARIANT,
    PROTEIN_ALTERING_VARIANT,

    // Low impact
    SPLICE_DONOR_5TH_BASE_VARIANT,
    SPLICE_DONOR_REGION_VARIANT,
    SPLICE_POLYPYRIMIDINE_TRACT_VARIANT,
    SPLICE_REGION_VARIANT,
    INCOMPLETE_TERMINAL_CODON_VARIANT,
    START_RETAINED_VARIANT,
    STOP_RETAINED_VARIANT,
    SYNONYMOUS_VARIANT,

    // Modifier
    CODING_SEQUENCE_VARIANT,
    CODING_TRANSCRIPT_VARIANT,
    MATURE_MIRNA_VARIANT,
    FIVE_PRIME_UTR_VARIANT,
    THREE_PRIME_UTR_VARIANT,
    NON_CODING_TRANSCRIPT_EXON_VARIANT,
    INTRON_VARIANT,
    NMD_TRANSCRIPT_VARIANT,
    NON_CODING_TRANSCRIPT_VARIANT,
    UPSTREAM_GENE_VARIANT,
    DOWNSTREAM_GENE_VARIANT,
    TFBS_ABLATION,
    TFBS_AMPLIFICATION,
    TF_BINDING_SITE_VARIANT,
    REGULATORY_REGION_ABLATION,
    REGULATORY_REGION_AMPLIFICATION,
    REGULATORY_REGION_VARIANT,
    INTERGENIC_VARIANT,
    SEQUENCE_VARIANT,

    UNKNOWN
};

/**
 * Impact severity levels
 */
enum class Impact {
    HIGH,
    MODERATE,
    LOW,
    MODIFIER
};

/**
 * Get string representation of consequence type
 */
std::string consequence_to_string(ConsequenceType type);

/**
 * Get impact level for a consequence type
 */
Impact get_impact(ConsequenceType type);

/**
 * Get string representation of impact
 */
std::string impact_to_string(Impact impact);

/**
 * Get SO variant class string from ref/alt alleles
 */
inline std::string get_variant_class(const std::string& ref, const std::string& alt) {
    if (ref.size() == 1 && alt.size() == 1) return "SNV";
    if (ref.size() == alt.size()) return "substitution";
    // Trim common prefix to determine insertion/deletion/indel
    size_t prefix = 0;
    while (prefix < ref.size() && prefix < alt.size() && ref[prefix] == alt[prefix]) {
        ++prefix;
    }
    std::string trimmed_ref = ref.substr(prefix);
    std::string trimmed_alt = alt.substr(prefix);
    if (trimmed_ref.empty() && !trimmed_alt.empty()) return "insertion";
    if (!trimmed_ref.empty() && trimmed_alt.empty()) return "deletion";
    if (!trimmed_ref.empty() && !trimmed_alt.empty()) return "indel";
    return "sequence_alteration";
}

/**
 * Get display term for a consequence type (Perl VEP Constants.pm format)
 */
inline std::string consequence_to_display_term(ConsequenceType type) {
    switch (type) {
        case ConsequenceType::TRANSCRIPT_ABLATION: return "TRANSCRIPT_ABLATION";
        case ConsequenceType::SPLICE_ACCEPTOR_VARIANT: return "ESSENTIAL_SPLICE_SITE";
        case ConsequenceType::SPLICE_DONOR_VARIANT: return "ESSENTIAL_SPLICE_SITE";
        case ConsequenceType::STOP_GAINED: return "STOP_GAINED";
        case ConsequenceType::FRAMESHIFT_VARIANT: return "FRAMESHIFT_CODING";
        case ConsequenceType::STOP_LOST: return "STOP_LOST";
        case ConsequenceType::START_LOST: return "NON_SYNONYMOUS_CODING";
        case ConsequenceType::TRANSCRIPT_AMPLIFICATION: return "TRANSCRIPT_AMPLIFICATION";
        case ConsequenceType::FEATURE_ELONGATION: return "FEATURE_ELONGATION";
        case ConsequenceType::FEATURE_TRUNCATION: return "FEATURE_TRUNCATION";
        case ConsequenceType::INFRAME_INSERTION: return "NON_SYNONYMOUS_CODING";
        case ConsequenceType::INFRAME_DELETION: return "NON_SYNONYMOUS_CODING";
        case ConsequenceType::MISSENSE_VARIANT: return "NON_SYNONYMOUS_CODING";
        case ConsequenceType::PROTEIN_ALTERING_VARIANT: return "PROTEIN_ALTERING";
        case ConsequenceType::SPLICE_DONOR_5TH_BASE_VARIANT: return "SPLICE_SITE";
        case ConsequenceType::SPLICE_REGION_VARIANT: return "SPLICE_SITE";
        case ConsequenceType::SPLICE_DONOR_REGION_VARIANT: return "SPLICE_SITE";
        case ConsequenceType::SPLICE_POLYPYRIMIDINE_TRACT_VARIANT: return "SPLICE_SITE";
        case ConsequenceType::INCOMPLETE_TERMINAL_CODON_VARIANT: return "PARTIAL_CODON";
        case ConsequenceType::START_RETAINED_VARIANT: return "SYNONYMOUS_CODING";
        case ConsequenceType::STOP_RETAINED_VARIANT: return "SYNONYMOUS_CODING";
        case ConsequenceType::SYNONYMOUS_VARIANT: return "SYNONYMOUS_CODING";
        case ConsequenceType::CODING_SEQUENCE_VARIANT: return "CODING_UNKNOWN";
        case ConsequenceType::CODING_TRANSCRIPT_VARIANT: return "CODING_TRANSCRIPT_VARIANT";
        case ConsequenceType::MATURE_MIRNA_VARIANT: return "WITHIN_MATURE_miRNA";
        case ConsequenceType::FIVE_PRIME_UTR_VARIANT: return "5PRIME_UTR";
        case ConsequenceType::THREE_PRIME_UTR_VARIANT: return "3PRIME_UTR";
        case ConsequenceType::NON_CODING_TRANSCRIPT_EXON_VARIANT: return "WITHIN_NON_CODING_GENE";
        case ConsequenceType::INTRON_VARIANT: return "INTRONIC";
        case ConsequenceType::NMD_TRANSCRIPT_VARIANT: return "NMD_TRANSCRIPT";
        case ConsequenceType::NON_CODING_TRANSCRIPT_VARIANT: return "WITHIN_NON_CODING_GENE";
        case ConsequenceType::UPSTREAM_GENE_VARIANT: return "UPSTREAM";
        case ConsequenceType::DOWNSTREAM_GENE_VARIANT: return "DOWNSTREAM";
        case ConsequenceType::TFBS_ABLATION: return "TFBS_ABLATION";
        case ConsequenceType::TFBS_AMPLIFICATION: return "TFBS_AMPLIFICATION";
        case ConsequenceType::TF_BINDING_SITE_VARIANT: return "REGULATORY_REGION";
        case ConsequenceType::REGULATORY_REGION_ABLATION: return "REGULATORY_ABLATION";
        case ConsequenceType::REGULATORY_REGION_AMPLIFICATION: return "REGULATORY_AMPLIFICATION";
        case ConsequenceType::REGULATORY_REGION_VARIANT: return "REGULATORY_REGION";
        case ConsequenceType::INTERGENIC_VARIANT: return "INTERGENIC";
        case ConsequenceType::SEQUENCE_VARIANT: return "SEQUENCE_VARIANT";
        default: return "UNKNOWN";
    }
}

/**
 * Represents an exon within a transcript
 */
struct Exon {
    int start;          // 1-based, inclusive
    int end;            // 1-based, inclusive
    int exon_number;
    int phase = -1;     // Reading frame phase (0, 1, 2, or -1 if unknown)
};

/**
 * Represents a CDS (coding sequence) region
 */
struct CDS {
    int start;
    int end;
    int phase = 0;
};

/**
 * Represents a transcript (mRNA)
 */
struct Transcript {
    std::string id;
    std::string gene_id;
    std::string gene_name;
    std::string chromosome;
    int start;
    int end;
    char strand;                    // '+' or '-'
    std::string biotype;            // protein_coding, lncRNA, etc.
    bool is_canonical = false;
    std::string ccds_id;
    std::string protein_id;         // Ensembl protein ID (ENSP)
    std::string version;            // Transcript version (e.g., "5")
    std::string mane_select;        // MANE Select tag value (e.g., "NM_000546.6")
    std::string mane_plus_clinical; // MANE Plus Clinical tag value
    int tsl = 0;                    // Transcript Support Level (1-5, 0=unknown)
    std::string appris;             // APPRIS annotation (e.g., "principal1")
    bool gencode_basic = false;     // Part of GENCODE basic set
    bool cds_start_NF = false;      // CDS start not found (incomplete 5' end)
    bool cds_end_NF = false;        // CDS end not found (incomplete 3' end)

    std::vector<Exon> exons;
    std::vector<CDS> cds_regions;

    int cds_start = 0;              // Translation start
    int cds_end = 0;                // Translation end

    // Computed properties
    int get_cds_length() const;
    bool is_coding() const { return !cds_regions.empty(); }
    bool overlaps(int pos) const { return pos >= start && pos <= end; }
    bool overlaps(int s, int e) const { return s <= end && e >= start; }
};

/**
 * Represents a gene
 */
struct Gene {
    std::string id;
    std::string name;
    std::string chromosome;
    int start;
    int end;
    char strand;
    std::string biotype;
    std::string source;             // Gene source (e.g., "ensembl_havana", "havana")
    std::string hgnc_id;            // HGNC gene ID (e.g., "HGNC:1097")

    std::vector<std::string> transcript_ids;
};

/**
 * Variant annotation result
 */
struct VariantAnnotation {
    std::string input_variant;
    std::string chromosome;
    int position = 0;
    std::string ref_allele;
    std::string alt_allele;

    // Gene/transcript info
    std::string gene_symbol;
    std::string gene_id;
    std::string transcript_id;
    std::string biotype;
    std::string feature_type = "Transcript";  // "Transcript", "RegulatoryFeature", "MotifFeature"
    bool is_canonical = false;

    // Consequence
    std::vector<ConsequenceType> consequences;
    Impact impact = Impact::MODIFIER;

    // Position details
    int exon_number = 0;
    int intron_number = 0;
    int total_exons = 0;
    int total_introns = 0;
    int cdna_position = 0;
    int cdna_end = 0;           // End position for multi-base variants
    int cds_position = 0;
    int cds_end = 0;            // End position for multi-base variants
    int protein_position = 0;
    int protein_end = 0;        // End position for multi-base variants
    int distance = 0;           // Distance to transcript for upstream/downstream variants
    char strand = '\0';         // Transcript strand ('+' or '-')
    std::string source;         // Transcript source ("Ensembl" or "RefSeq")

    // Display alleles (Ensembl format: stripped shared leading base, "-" for empty)
    std::string display_ref;            // Ensembl-format REF (e.g., "TG" for VCF REF="ATG")
    std::string display_alt;            // Ensembl-format ALT (e.g., "-" for VCF deletion)
    int display_start = 0;              // Ensembl-format start position
    int display_end = 0;                // Ensembl-format end position

    // Sequence changes
    std::string codons;             // REF/ALT codons (e.g., "Gcc/Acc")
    std::string amino_acids;        // REF/ALT amino acids (e.g., "A/T")
    std::string hgvsc;              // HGVS coding notation
    std::string hgvsp;              // HGVS protein notation
    std::string hgvsg;              // HGVS genomic notation

    // Co-located known variants (e.g., rs IDs from dbSNP)
    std::string existing_variation;

    // Original VCF fields for passthrough
    std::string vcf_id;              // VCF ID column (col 2)
    std::string vcf_qual;            // VCF QUAL column (col 5)
    std::string vcf_filter;          // VCF FILTER column (col 6)
    std::string vcf_info;            // VCF INFO column (col 7) - original INFO for passthrough
    std::string vcf_sample_columns;  // FORMAT + sample columns (col 8 onwards)
    std::string vcf_ref;             // Original VCF REF column (for multi-allele output)
    std::string vcf_alt;             // Original VCF ALT column (comma-separated for multi-allele)

    // Custom annotations from VCF files
    // Key: "source:field" (e.g., "gnomad:AF", "clinvar:CLNSIG")
    // Value: annotation value as string
    std::unordered_map<std::string, std::string> custom_annotations;

    /**
     * Get a custom annotation value
     * @param source Source name (e.g., "gnomad")
     * @param field Field name (e.g., "AF")
     * @return Optional value if found
     */
    std::optional<std::string> get_annotation(const std::string& source, const std::string& field) const;

    /**
     * Get a custom annotation as double
     */
    std::optional<double> get_annotation_double(const std::string& source, const std::string& field) const;

    /**
     * Get consequence string
     */
    std::string get_consequence_string() const;

    /**
     * Get summary map
     */
    std::map<std::string, std::string> get_summary() const;
};

/**
 * FASTA reference genome reader
 */
class ReferenceGenome {
public:
    /**
     * Load reference genome from FASTA file
     * @param fasta_path Path to .fa or .fa.gz file
     * @param load_all If true, load all sequences into memory; if false, use indexed access
     */
    explicit ReferenceGenome(const std::string& fasta_path, bool load_all = false);
    ~ReferenceGenome();

    /**
     * Get sequence at a position
     * @param chrom Chromosome name
     * @param start 1-based start position
     * @param end 1-based end position (inclusive)
     * @return Uppercase sequence string
     */
    std::string get_sequence(const std::string& chrom, int start, int end) const;

    /**
     * Get single base at position
     */
    char get_base(const std::string& chrom, int pos) const;

    /**
     * Check if chromosome exists
     */
    bool has_chromosome(const std::string& chrom) const;

    /**
     * Get chromosome length
     */
    int get_chromosome_length(const std::string& chrom) const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

/**
 * GTF/GFF annotation file parser and database
 */
class TranscriptDatabase {
public:
    /**
     * Load annotations from GTF/GFF file
     * @param gtf_path Path to .gtf or .gtf.gz file
     */
    explicit TranscriptDatabase(const std::string& gtf_path);
    ~TranscriptDatabase();

    /**
     * Get all transcripts overlapping a position
     */
    std::vector<const Transcript*> get_transcripts_at(const std::string& chrom, int pos) const;

    /**
     * Get all transcripts overlapping a region
     */
    std::vector<const Transcript*> get_transcripts_in_region(
        const std::string& chrom, int start, int end) const;

    /**
     * Get transcript by ID
     */
    const Transcript* get_transcript(const std::string& transcript_id) const;

    /**
     * Get gene by ID
     */
    const Gene* get_gene(const std::string& gene_id) const;

    /**
     * Get all transcripts for a gene
     */
    std::vector<const Transcript*> get_gene_transcripts(const std::string& gene_id) const;

    /**
     * Get nearby genes (for upstream/downstream annotation)
     */
    std::vector<const Gene*> get_nearby_genes(
        const std::string& chrom, int pos, int distance = 5000) const;

    /**
     * Get total number of transcripts loaded
     */
    size_t transcript_count() const;

    /**
     * Get total number of genes loaded
     */
    size_t gene_count() const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

/**
 * Codon translation table
 */
class CodonTable {
public:
    /**
     * Translate a codon to amino acid
     * @param codon 3-letter DNA codon (uppercase)
     * @return Single-letter amino acid code, or '*' for stop, or 'X' for unknown
     */
    static char translate(const std::string& codon);

    /**
     * Get three-letter amino acid code
     */
    static std::string get_three_letter(char aa);

    /**
     * Check if codon is a start codon
     */
    static bool is_start_codon(const std::string& codon);

    /**
     * Check if codon is a stop codon
     */
    static bool is_stop_codon(const std::string& codon);
    static bool is_stop_codon(const std::string& codon, const std::string& chromosome);

    /**
     * Translate using vertebrate mitochondrial codon table
     */
    static char translate_mt(const std::string& codon);

    /**
     * Translate with automatic table selection based on chromosome
     */
    static char translate(const std::string& codon, const std::string& chromosome);
};

/**
 * Configuration for a custom VCF annotation source
 */
struct VCFAnnotationConfig {
    std::string name;                       // Source name (e.g., "gnomad", "clinvar")
    std::string vcf_path;                   // Path to VCF file
    std::vector<std::string> info_fields;   // INFO fields to extract (empty = all)
    bool match_allele = true;               // Require REF/ALT match (not just position)
    bool use_tabix = false;                 // Use tabix for on-disk queries (for large files)
};

/**
 * Represents a single VCF record for annotation lookup
 */
struct VCFRecord {
    std::string chrom;
    int pos;
    std::string id;
    std::string ref;
    std::vector<std::string> alts;
    std::map<std::string, std::string> info;  // Parsed INFO fields
};

/**
 * Database for custom VCF annotations
 * Loads VCF files and provides fast position-based lookups
 */
class VCFAnnotationDatabase {
public:
    VCFAnnotationDatabase();
    ~VCFAnnotationDatabase();

    /**
     * Add a VCF annotation source
     * @param config Configuration for the VCF source
     */
    void add_source(const VCFAnnotationConfig& config);

    /**
     * Add a VCF annotation source with simple parameters
     * @param name Source name
     * @param vcf_path Path to VCF file
     * @param info_fields INFO fields to extract (comma-separated, or empty for all)
     */
    void add_source(const std::string& name, const std::string& vcf_path,
                    const std::string& info_fields = "");

    /**
     * Get annotations for a variant
     * @param chrom Chromosome
     * @param pos Position
     * @param ref Reference allele
     * @param alt Alternative allele
     * @return Map of "source:field" -> value
     */
    std::map<std::string, std::string> get_annotations(
        const std::string& chrom, int pos,
        const std::string& ref, const std::string& alt) const;

    /**
     * Get all records at a position (regardless of allele)
     */
    std::vector<const VCFRecord*> get_records_at(const std::string& chrom, int pos) const;

    /**
     * Get list of loaded sources
     */
    std::vector<std::string> get_sources() const;

    /**
     * Get fields available for a source
     */
    std::vector<std::string> get_fields(const std::string& source) const;

    /**
     * Get statistics about loaded data
     */
    std::string get_stats() const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

/**
 * Main VEP annotator class
 */
class VEPAnnotator {
public:
    /**
     * Initialize annotator with required data files
     * @param gtf_path Path to GTF/GFF annotation file
     * @param fasta_path Path to reference FASTA file
     */
    VEPAnnotator(const std::string& gtf_path, const std::string& fasta_path);
    ~VEPAnnotator();

    // Prevent copying
    VEPAnnotator(const VEPAnnotator&) = delete;
    VEPAnnotator& operator=(const VEPAnnotator&) = delete;

    /**
     * Add a custom VCF annotation source
     * @param name Source name (e.g., "gnomad", "clinvar", "mydata")
     * @param vcf_path Path to VCF file (.vcf or .vcf.gz)
     * @param info_fields Comma-separated INFO fields to extract (empty = all)
     *
     * Example:
     *   annotator.add_annotation_source("gnomad", "gnomad.vcf.gz", "AF,AF_popmax");
     *   annotator.add_annotation_source("clinvar", "clinvar.vcf.gz", "CLNSIG,CLNDN");
     */
    void add_annotation_source(const std::string& name, const std::string& vcf_path,
                               const std::string& info_fields = "");

    /**
     * Add a custom VCF annotation source with full configuration
     */
    void add_annotation_source(const VCFAnnotationConfig& config);

    /**
     * Get list of loaded annotation sources
     */
    std::vector<std::string> get_annotation_sources() const;

    /**
     * Annotate a variant
     * @param chrom Chromosome
     * @param pos 1-based position
     * @param ref Reference allele
     * @param alt Alternative allele
     * @return Vector of annotations (one per overlapping transcript)
     */
    std::vector<VariantAnnotation> annotate(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt
    );

    /**
     * Annotate a variant, returning only the most severe consequence
     */
    VariantAnnotation annotate_most_severe(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt
    );

    /**
     * Get annotation database stats
     */
    std::string get_stats() const;

    // =========================================================================
    // Extended Annotation Source System
    // =========================================================================

    /**
     * Add a custom annotation source (pathogenicity, conservation, etc.)
     * @param source Shared pointer to annotation source
     */
    void add_source(std::shared_ptr<AnnotationSource> source);

    /**
     * Get all registered annotation sources
     */
    std::vector<std::shared_ptr<AnnotationSource>> get_sources() const;

    /**
     * Enable/disable an annotation source by name
     */
    void set_source_enabled(const std::string& name, bool enabled);

    /**
     * Get list of all available annotation fields
     * @return Vector of (field_name, source_type) pairs
     */
    std::vector<std::pair<std::string, std::string>> get_available_fields() const;

    /**
     * Initialize all annotation sources (call before annotating)
     */
    void initialize_sources();

    /**
     * Get transcript by ID (public access for HGVS resolution)
     */
    const Transcript* get_transcript(const std::string& transcript_id) const;

    /**
     * Map CDS position back to genomic coordinate
     */
    int map_cds_to_genomic(int cds_pos, const Transcript& transcript) const;

    void set_upstream_distance(int d) { upstream_distance_ = d; }
    void set_downstream_distance(int d) { downstream_distance_ = d; }
    void set_shift_3prime(bool v) { shift_3prime_ = v; }
    void set_shift_genomic(bool v) { shift_genomic_ = v; }
    void set_transcript_version(bool v) { transcript_version_ = v; }
    void set_include_numbers(bool v) { include_numbers_ = v; }
    void set_include_total_length(bool v) { include_total_length_ = v; }

    const ReferenceGenome* get_reference() const { return reference_.get(); }

private:
    std::unique_ptr<TranscriptDatabase> transcript_db_;
    std::unique_ptr<ReferenceGenome> reference_;
    std::unique_ptr<VCFAnnotationDatabase> vcf_annotations_;
    std::unique_ptr<AnnotationSourceManager> source_manager_;

    int upstream_distance_ = 5000;
    int downstream_distance_ = 5000;
    bool shift_3prime_ = false;
    bool shift_genomic_ = false;
    bool transcript_version_ = false;
    bool include_numbers_ = false;
    bool include_total_length_ = false;

    /**
     * Annotate variant against a single transcript
     */
    VariantAnnotation annotate_transcript(
        const std::string& chrom,
        int pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript& transcript
    );

    /**
     * Determine consequence type for a variant in a transcript
     */
    std::vector<ConsequenceType> determine_consequences(
        int pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript& transcript
    );

    /**
     * Calculate CDS position for a genomic position
     */
    int calculate_cds_position(int genomic_pos, const Transcript& transcript) const;

    /**
     * Build the CDS sequence for a transcript (builds once, reuse everywhere)
     */
    std::string build_cds_sequence(const std::string& chrom, const Transcript& transcript) const;

    /**
     * Get the codon affected by a variant
     */
    std::pair<std::string, std::string> get_affected_codons(
        int cds_pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript& transcript
    ) const;

    /** Overload accepting pre-built CDS sequence */
    std::pair<std::string, std::string> get_affected_codons(
        int cds_pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript& transcript,
        const std::string& cds_seq
    ) const;

    /**
     * Annotate coding region: CDS position, codons, amino acids, HGVSc, HGVSp
     */
    void annotate_coding_region(
        const std::string& chrom, int pos,
        const std::string& ref, const std::string& alt,
        const Transcript& transcript,
        const std::string& cached_cds,
        VariantAnnotation& ann);

    /**
     * Annotate non-CDS HGVSc: intronic (c.N+/-offset) and UTR (c.-N, c.*N) notation
     */
    void annotate_noncds_hgvsc(
        int pos,
        const std::string& ref, const std::string& alt,
        const Transcript& transcript,
        VariantAnnotation& ann);

    /**
     * Populate transcript metadata: CCDS, ENSP, CANONICAL, MANE, TSL, APPRIS, etc.
     */
    void populate_transcript_metadata(
        const std::string& chrom, int pos,
        const std::string& ref, const std::string& alt,
        const Transcript& transcript,
        int hgvs_offset,
        VariantAnnotation& ann);

    /**
     * Append regulatory consequences based on annotation source flags
     */
    void append_regulatory_consequences(VariantAnnotation& ann);

    /**
     * Right-normalize (3' shift) an indel for HGVS generation
     */
    std::tuple<int, std::string, std::string> right_normalize(
        const std::string& chrom, int pos,
        const std::string& ref, const std::string& alt,
        const Transcript& transcript) const;

    /**
     * Calculate cDNA position (exonic position including UTRs)
     */
    int calculate_cdna_position(int genomic_pos, const Transcript& transcript) const;

    /**
     * Generate HGVS notation
     */
    std::string generate_hgvsc(
        int cds_pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript& transcript
    ) const;

    /** Overload accepting pre-built CDS sequence */
    std::string generate_hgvsc(
        int cds_pos,
        const std::string& ref,
        const std::string& alt,
        const Transcript& transcript,
        const std::string& cds_seq
    ) const;

    std::string generate_hgvsp(
        const std::string& ref_aa,
        const std::string& alt_aa,
        int protein_pos,
        const Transcript& transcript,
        const std::vector<ConsequenceType>& consequences = {},
        int protein_end = 0,
        int fs_ter_distance = 0,
        const std::string& end_ref_aa = ""
    ) const;
};

/**
 * Annotate variants from a VCF file
 * @param vcf_input Input VCF file path
 * @param output_path Output TSV file path
 * @param gtf_path GTF annotation file path
 * @param fasta_path Reference FASTA file path
 * @param annotation_vcfs Optional vector of (name, vcf_path, info_fields, use_tabix) tuples
 */
void annotate_vcf_file(
    const std::string& vcf_input,
    const std::string& output_path,
    const std::string& gtf_path,
    const std::string& fasta_path,
    const std::vector<std::tuple<std::string, std::string, std::string, bool>>& annotation_vcfs = {}
);

/**
 * Logging utilities
 */
enum class LogLevel { DEBUG, INFO, WARNING, ERROR };
void set_log_level(LogLevel level);
void log(LogLevel level, const std::string& message);

} // namespace vep

#endif // VEP_ANNOTATOR_HPP
