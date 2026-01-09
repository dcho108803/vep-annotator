/**
 * VEP Variant Annotator - Main Entry Point
 *
 * Pure C++ local implementation for annotating genetic variants.
 * Requires local GTF and FASTA files - no external API calls.
 */

#include "vep_annotator.hpp"
#include "annotation_sources.hpp"
#include "plugin.hpp"
#include "output_writer.hpp"
#include "transcript_filter.hpp"
#include "structural_variant.hpp"
#include "hgvs_parser.hpp"
#include "gene_constraint.hpp"
#include "exon_intron_numbers.hpp"
#include "filter_vep.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <set>
#include <cctype>
#include <memory>

void print_usage(const char* program_name) {
    std::cout << "VEP Variant Annotator - Pure C++ Local Implementation\n"
              << "======================================================\n\n"
              << "Usage: " << program_name << " [OPTIONS]\n\n"
              << "Required Data Files:\n"
              << "  --gtf FILE              GTF/GFF annotation file (genes/transcripts)\n"
              << "  --fasta FILE            Reference genome FASTA file\n\n"
              << "Variant Input (choose one):\n"
              << "  -v, --variant CHR:POS:REF:ALT   Single variant to annotate\n"
              << "  --hgvs NOTATION         HGVS notation (e.g., ENST00000366667:c.803C>T)\n"
              << "  --vcf FILE              VCF file for batch annotation\n\n"
              << "Output Format Options:\n"
              << "  -o, --output FILE       Output file path\n"
              << "  --output-format FORMAT  Output format: tsv (default), json, vcf\n"
              << "  --compress              Compress output with gzip (.gz)\n"
              << "  --stats                 Print statistics summary at end\n"
              << "  --hgvsg                 Include genomic HGVS (g.) notation in output\n"
              << "  --spdi                  Include SPDI format in output\n"
              << "  --numbers               Include exon/intron numbers in output\n\n"
              << "Transcript Selection:\n"
              << "  --all-transcripts       Output all transcript annotations\n"
              << "  --pick                  Pick one consequence per variant\n"
              << "  --pick-allele           Pick one consequence per allele\n"
              << "  --per-gene              Pick one consequence per gene\n"
              << "  --most-severe           Only output most severe consequence\n"
              << "  --flag-pick             Flag picked transcript without filtering\n"
              << "  --pick-order ORDER      Comma-separated pick order criteria:\n"
              << "                          canonical,mane,appris,tsl,biotype,ccds,rank,length\n\n"
              << "Transcript Filtering:\n"
              << "  --canonical             Only annotate canonical transcripts\n"
              << "  --mane                  Only annotate MANE Select transcripts\n"
              << "  --coding-only           Only annotate protein-coding transcripts\n"
              << "  --biotype LIST          Only annotate these biotypes (comma-separated)\n"
              << "  --no-intergenic         Skip intergenic variants\n\n"
              << "Frequency Filtering:\n"
              << "  --check-frequency       Enable frequency-based filtering\n"
              << "  --freq-pop FIELD        Frequency field to check (default: gnomAD_AF)\n"
              << "  --freq-threshold VAL    Maximum frequency threshold (default: 0.01)\n"
              << "  --freq-gt               Filter for frequency greater than threshold\n\n"
              << "Custom Annotations (can be used multiple times):\n"
              << "  --annotation NAME:VCF_PATH[:FIELDS]\n"
              << "                          Add custom VCF annotation source (loads into memory)\n"
              << "  --annotation-tabix NAME:VCF_PATH[:FIELDS]\n"
              << "                          Add tabix-indexed VCF source (on-disk queries)\n\n"
              << "Pathogenicity Predictions:\n"
              << "  --dbnsfp FILE           dbNSFP database (tabix-indexed .txt.gz)\n"
              << "  --dbnsfp-fields FIELDS  Fields to extract (default: essential)\n"
              << "                          Presets: essential, pathogenicity, conservation, all\n\n"
              << "Splice Predictions:\n"
              << "  --spliceai FILE         SpliceAI VCF file (tabix-indexed)\n"
              << "  --maxentscan            Enable MaxEntScan splice site scoring\n"
              << "  --dbscsnv FILE          dbscSNV file (tabix-indexed)\n\n"
              << "Conservation Scores:\n"
              << "  --phylop FILE           PhyloP bigWig file (requires libBigWig)\n"
              << "  --phastcons FILE        PhastCons bigWig file (requires libBigWig)\n"
              << "  --gerp FILE             GERP++ bigWig file (requires libBigWig)\n\n"
              << "Regulatory Annotations:\n"
              << "  --regulatory FILE       Ensembl Regulatory Build GFF3 file\n\n"
              << "Protein Domains:\n"
              << "  --pfam FILE             Pfam domain annotations TSV\n"
              << "  --interpro FILE         InterPro domain annotations TSV\n\n"
              << "Loss-of-Function:\n"
              << "  --loftee                Enable LOFTEE-style LoF classification\n"
              << "  --nmd                   Enable NMD prediction\n"
              << "  --loftool FILE          LoFtool gene constraint scores\n\n"
              << "Gene Constraint Scores:\n"
              << "  --constraint FILE       gnomAD gene constraint file (TSV)\n"
              << "  --pli FILE              pLI scores file (GENE\\tpLI)\n"
              << "  --loeuf FILE            LOEUF scores file (GENE\\tLOEUF)\n\n"
              << "Plugins:\n"
              << "  --plugin PATH[:CONFIG]  Load plugin from shared library\n"
              << "  --plugin-dir DIR        Directory to search for plugins\n\n"
              << "Shortcut Flags:\n"
              << "  --everything            Enable all output fields (equivalent to multiple flags)\n\n"
              << "Identifiers:\n"
              << "  --ccds                  Add CCDS transcript identifiers\n"
              << "  --uniprot               Add UniProt identifiers\n"
              << "  --xref-refseq           Add RefSeq cross-references\n\n"
              << "Reference Checking:\n"
              << "  --check-ref             Check reference allele matches FASTA\n"
              << "  --allow-non-variant     Allow non-variant positions (e.g., REF=A, ALT=A)\n\n"
              << "Species/Assembly:\n"
              << "  --species NAME          Species name (default: homo_sapiens)\n"
              << "  --assembly NAME         Assembly name (default: GRCh38)\n\n"
              << "Additional Pick Options:\n"
              << "  --pick-allele-gene      Pick one consequence per allele and gene\n"
              << "  --flag-pick-allele      Flag picked allele without filtering\n"
              << "  --flag-pick-allele-gene Flag picked allele+gene without filtering\n"
              << "  --summary               Output summary of all consequences\n"
              << "  --total-length          Include total transcript length in output\n\n"
              << "Additional Filtering:\n"
              << "  --gencode-basic         Only use GENCODE basic transcripts\n"
              << "  --all-refseq            Include all RefSeq transcripts\n"
              << "  --filter-common         Filter out common variants (AF > 0.01)\n"
              << "  --max-af                Report maximum AF across populations\n"
              << "  --overlaps TYPE         SV overlap type: any, within, surrounding, exact\n\n"
              << "Performance Options:\n"
              << "  --buffer-size N         Number of variants to buffer (default: 5000)\n"
              << "  --quiet                 Suppress progress messages\n"
              << "  --no-progress           Don't show progress bar\n"
              << "  --minimal               Only output most essential fields\n\n"
              << "Other Options:\n"
              << "  -h, --help              Show this help message\n"
              << "  --debug                 Enable debug logging\n\n"
              << "Examples:\n"
              << "  # Annotate a single SNV\n"
              << "  " << program_name << " --gtf genes.gtf --fasta genome.fa -v 7:140753336:A:T\n\n"
              << "  # Annotate with JSON output and transcript picking\n"
              << "  " << program_name << " --gtf genes.gtf --fasta genome.fa \\\n"
              << "      --output-format json --pick --canonical \\\n"
              << "      --vcf variants.vcf -o results.json\n\n"
              << "  # Annotate with VCF output (CSQ INFO field)\n"
              << "  " << program_name << " --gtf genes.gtf --fasta genome.fa \\\n"
              << "      --output-format vcf --compress \\\n"
              << "      --vcf variants.vcf -o annotated.vcf.gz\n\n"
              << "  # Full annotation pipeline with filtering\n"
              << "  " << program_name << " --gtf genes.gtf --fasta genome.fa \\\n"
              << "      --dbnsfp dbNSFP.txt.gz --spliceai spliceai.vcf.gz \\\n"
              << "      --pick --coding-only --check-frequency --freq-threshold 0.001 \\\n"
              << "      --annotation-tabix gnomad:gnomad.vcf.bgz:AF \\\n"
              << "      --vcf variants.vcf -o results.tsv --stats\n"
              << std::endl;
}

void print_annotation(const vep::VariantAnnotation& ann) {
    std::cout << "\n=== Variant Annotation ===" << std::endl;

    auto summary = ann.get_summary();
    for (const auto& kv : summary) {
        if (!kv.second.empty()) {
            std::cout << kv.first << ": " << kv.second << std::endl;
        }
    }
}

// Parse variant string in format CHR:POS:REF:ALT or CHR-POS-REF-ALT
bool parse_variant(const std::string& variant, std::string& chrom, int& pos,
                   std::string& ref, std::string& alt) {
    std::string normalized = variant;

    // Replace - with : for uniform parsing
    for (size_t i = 0; i < normalized.size(); ++i) {
        if (normalized[i] == '-') normalized[i] = ':';
    }

    size_t p1 = normalized.find(':');
    if (p1 == std::string::npos) return false;

    size_t p2 = normalized.find(':', p1 + 1);
    if (p2 == std::string::npos) return false;

    size_t p3 = normalized.find(':', p2 + 1);
    if (p3 == std::string::npos) return false;

    try {
        chrom = normalized.substr(0, p1);
        pos = std::stoi(normalized.substr(p1 + 1, p2 - p1 - 1));
        ref = normalized.substr(p2 + 1, p3 - p2 - 1);
        alt = normalized.substr(p3 + 1);
        return true;
    } catch (...) {
        return false;
    }
}

// Parse annotation source string: NAME:VCF_PATH[:FIELDS]
bool parse_annotation_source(const std::string& source,
                             std::string& name, std::string& vcf_path, std::string& fields) {
    size_t p1 = source.find(':');
    if (p1 == std::string::npos) return false;

    name = source.substr(0, p1);

    size_t p2 = source.find(':', p1 + 1);
    if (p2 == std::string::npos) {
        vcf_path = source.substr(p1 + 1);
        fields = "";
    } else {
        vcf_path = source.substr(p1 + 1, p2 - p1 - 1);
        fields = source.substr(p2 + 1);
    }

    return !name.empty() && !vcf_path.empty();
}

// Helper function to verify reference allele against FASTA
bool verify_reference(const vep::ReferenceGenome& reference, const std::string& chrom,
                     int pos, const std::string& ref, bool quiet) {
    std::string fasta_ref = reference.get_sequence(chrom, pos, pos + static_cast<int>(ref.length()) - 1);

    // Convert to uppercase for comparison
    std::string ref_upper = ref;
    std::string fasta_upper = fasta_ref;
    for (size_t i = 0; i < ref_upper.size(); ++i) {
        ref_upper[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(ref_upper[i])));
    }
    for (size_t i = 0; i < fasta_upper.size(); ++i) {
        fasta_upper[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(fasta_upper[i])));
    }

    if (ref_upper != fasta_upper) {
        if (!quiet) {
            std::cerr << "Warning: Reference allele mismatch at " << chrom << ":" << pos
                      << " - VCF: " << ref << ", FASTA: " << fasta_ref << std::endl;
        }
        return false;
    }
    return true;
}

int main(int argc, char* argv[]) {
    // Basic options
    std::string gtf_path;
    std::string fasta_path;
    std::string variant;
    std::string vcf_path;
    std::string output_path;
    bool debug = false;

    // Output format options
    vep::OutputFormat output_format = vep::OutputFormat::TSV;
    bool compress_output = false;
    bool show_stats = false;
    bool include_hgvsg = false;
    bool include_spdi = false;
    bool include_numbers = false;

    // HGVS input
    std::string hgvs_input;

    // Transcript selection options
    vep::TranscriptFilterConfig filter_config;
    bool all_transcripts = false;

    // Custom annotation sources: (name, vcf_path, fields, use_tabix)
    std::vector<std::tuple<std::string, std::string, std::string, bool> > annotation_sources;

    // Pathogenicity predictions
    std::string dbnsfp_path;
    std::string dbnsfp_fields = "essential";

    // Splice predictions
    std::string spliceai_path;
    bool use_maxentscan = false;
    std::string dbscsnv_path;

    // Conservation scores
    std::string phylop_path;
    std::string phastcons_path;
    std::string gerp_path;

    // Regulatory annotations
    std::string regulatory_path;

    // Protein domains
    std::string pfam_path;
    std::string interpro_path;

    // LoF annotations
    bool use_loftee = false;
    bool use_nmd = false;
    std::string loftool_path;

    // Gene constraint scores
    std::string constraint_path;
    std::string pli_path;
    std::string loeuf_path;

    // Plugins
    std::vector<std::pair<std::string, std::string> > plugins;
    std::vector<std::string> plugin_dirs;

    // New options - Shortcut flags
    bool use_everything = false;

    // Identifier options
    bool include_ccds = false;
    bool include_uniprot = false;
    bool include_xref_refseq = false;

    // Reference checking
    bool check_ref = false;
    bool allow_non_variant = false;

    // Species/Assembly
    std::string species = "homo_sapiens";
    std::string assembly = "GRCh38";

    // Additional pick options
    bool pick_allele_gene = false;
    bool flag_pick_allele = false;
    bool flag_pick_allele_gene = false;
    bool show_summary = false;
    bool include_total_length = false;

    // Additional filtering
    bool gencode_basic = false;
    bool all_refseq = false;
    bool filter_common = false;
    bool show_max_af = false;
    std::string overlaps_type = "any";

    // Performance options
    int buffer_size = 5000;
    bool quiet_mode = false;
    bool no_progress = false;
    bool minimal_output = false;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        }
        // Basic options
        else if (arg == "--gtf" && i + 1 < argc) {
            gtf_path = argv[++i];
        } else if (arg == "--fasta" && i + 1 < argc) {
            fasta_path = argv[++i];
        } else if ((arg == "-v" || arg == "--variant") && i + 1 < argc) {
            variant = argv[++i];
        } else if (arg == "--vcf" && i + 1 < argc) {
            vcf_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            output_path = argv[++i];
        } else if (arg == "--debug") {
            debug = true;
        }
        // Output format options
        else if (arg == "--output-format" && i + 1 < argc) {
            output_format = vep::parse_output_format(argv[++i]);
        } else if (arg == "--compress") {
            compress_output = true;
        } else if (arg == "--stats") {
            show_stats = true;
        } else if (arg == "--hgvsg") {
            include_hgvsg = true;
        } else if (arg == "--spdi") {
            include_spdi = true;
        } else if (arg == "--numbers") {
            include_numbers = true;
        } else if (arg == "--hgvs" && i + 1 < argc) {
            hgvs_input = argv[++i];
        }
        // Transcript selection options
        else if (arg == "--all-transcripts") {
            all_transcripts = true;
        } else if (arg == "--pick") {
            filter_config.pick = true;
        } else if (arg == "--pick-allele") {
            filter_config.pick_allele = true;
        } else if (arg == "--per-gene") {
            filter_config.per_gene = true;
        } else if (arg == "--most-severe") {
            filter_config.most_severe = true;
        } else if (arg == "--flag-pick") {
            filter_config.flag_pick = true;
        } else if (arg == "--pick-order" && i + 1 < argc) {
            filter_config.pick_order = vep::parse_pick_order(argv[++i]);
        }
        // Transcript filtering options
        else if (arg == "--canonical") {
            filter_config.canonical_only = true;
        } else if (arg == "--mane") {
            filter_config.mane_only = true;
        } else if (arg == "--coding-only") {
            filter_config.coding_only = true;
        } else if (arg == "--biotype" && i + 1 < argc) {
            filter_config.biotypes = vep::parse_biotypes(argv[++i]);
        } else if (arg == "--no-intergenic") {
            filter_config.no_intergenic = true;
        }
        // Frequency filtering options
        else if (arg == "--check-frequency") {
            filter_config.check_frequency = true;
        } else if (arg == "--freq-pop" && i + 1 < argc) {
            filter_config.freq_pop = argv[++i];
        } else if (arg == "--freq-threshold" && i + 1 < argc) {
            filter_config.freq_threshold = std::stod(argv[++i]);
        } else if (arg == "--freq-gt") {
            filter_config.freq_gt = true;
        }
        // Custom annotation sources
        else if (arg == "--annotation" && i + 1 < argc) {
            std::string source_str = argv[++i];
            std::string name, path, fields;
            if (parse_annotation_source(source_str, name, path, fields)) {
                annotation_sources.push_back(std::make_tuple(name, path, fields, false));
            } else {
                std::cerr << "Error: Invalid annotation format: " << source_str << "\n"
                          << "Expected: NAME:VCF_PATH[:FIELDS]" << std::endl;
                return 1;
            }
        } else if (arg == "--annotation-tabix" && i + 1 < argc) {
            std::string source_str = argv[++i];
            std::string name, path, fields;
            if (parse_annotation_source(source_str, name, path, fields)) {
                annotation_sources.push_back(std::make_tuple(name, path, fields, true));
            } else {
                std::cerr << "Error: Invalid annotation format: " << source_str << "\n"
                          << "Expected: NAME:VCF_PATH[:FIELDS]" << std::endl;
                return 1;
            }
        }
        // Pathogenicity predictions
        else if (arg == "--dbnsfp" && i + 1 < argc) {
            dbnsfp_path = argv[++i];
        } else if (arg == "--dbnsfp-fields" && i + 1 < argc) {
            dbnsfp_fields = argv[++i];
        }
        // Splice predictions
        else if (arg == "--spliceai" && i + 1 < argc) {
            spliceai_path = argv[++i];
        } else if (arg == "--maxentscan") {
            use_maxentscan = true;
        } else if (arg == "--dbscsnv" && i + 1 < argc) {
            dbscsnv_path = argv[++i];
        }
        // Conservation scores
        else if (arg == "--phylop" && i + 1 < argc) {
            phylop_path = argv[++i];
        } else if (arg == "--phastcons" && i + 1 < argc) {
            phastcons_path = argv[++i];
        } else if (arg == "--gerp" && i + 1 < argc) {
            gerp_path = argv[++i];
        }
        // Regulatory annotations
        else if (arg == "--regulatory" && i + 1 < argc) {
            regulatory_path = argv[++i];
        }
        // Protein domains
        else if (arg == "--pfam" && i + 1 < argc) {
            pfam_path = argv[++i];
        } else if (arg == "--interpro" && i + 1 < argc) {
            interpro_path = argv[++i];
        }
        // LoF annotations
        else if (arg == "--loftee") {
            use_loftee = true;
        } else if (arg == "--nmd") {
            use_nmd = true;
        } else if (arg == "--loftool" && i + 1 < argc) {
            loftool_path = argv[++i];
        }
        // Gene constraint scores
        else if (arg == "--constraint" && i + 1 < argc) {
            constraint_path = argv[++i];
        } else if (arg == "--pli" && i + 1 < argc) {
            pli_path = argv[++i];
        } else if (arg == "--loeuf" && i + 1 < argc) {
            loeuf_path = argv[++i];
        }
        // Plugins
        else if (arg == "--plugin" && i + 1 < argc) {
            std::string plugin_arg = argv[++i];
            size_t colon = plugin_arg.find(':');
            if (colon != std::string::npos) {
                plugins.push_back(std::make_pair(plugin_arg.substr(0, colon), plugin_arg.substr(colon + 1)));
            } else {
                plugins.push_back(std::make_pair(plugin_arg, std::string("")));
            }
        } else if (arg == "--plugin-dir" && i + 1 < argc) {
            plugin_dirs.push_back(argv[++i]);
        }
        // Shortcut flags
        else if (arg == "--everything") {
            use_everything = true;
            // Enable all output fields
            include_hgvsg = true;
            include_spdi = true;
            include_numbers = true;
            include_ccds = true;
            include_uniprot = true;
            include_xref_refseq = true;
            include_total_length = true;
            show_max_af = true;
            show_stats = true;
        }
        // Identifier options
        else if (arg == "--ccds") {
            include_ccds = true;
        } else if (arg == "--uniprot") {
            include_uniprot = true;
        } else if (arg == "--xref-refseq") {
            include_xref_refseq = true;
        }
        // Reference checking
        else if (arg == "--check-ref") {
            check_ref = true;
        } else if (arg == "--allow-non-variant") {
            allow_non_variant = true;
        }
        // Species/Assembly
        else if (arg == "--species" && i + 1 < argc) {
            species = argv[++i];
        } else if (arg == "--assembly" && i + 1 < argc) {
            assembly = argv[++i];
        }
        // Additional pick options
        else if (arg == "--pick-allele-gene") {
            pick_allele_gene = true;
            filter_config.pick_allele_gene = true;
        } else if (arg == "--flag-pick-allele") {
            flag_pick_allele = true;
            filter_config.flag_pick_allele = true;
        } else if (arg == "--flag-pick-allele-gene") {
            flag_pick_allele_gene = true;
            filter_config.flag_pick_allele_gene = true;
        } else if (arg == "--summary") {
            show_summary = true;
        } else if (arg == "--total-length") {
            include_total_length = true;
        }
        // Additional filtering
        else if (arg == "--gencode-basic") {
            gencode_basic = true;
            filter_config.gencode_basic = true;
        } else if (arg == "--all-refseq") {
            all_refseq = true;
            filter_config.all_refseq = true;
        } else if (arg == "--filter-common") {
            filter_common = true;
            filter_config.check_frequency = true;
            filter_config.freq_threshold = 0.01;
        } else if (arg == "--max-af") {
            show_max_af = true;
        } else if (arg == "--overlaps" && i + 1 < argc) {
            overlaps_type = argv[++i];
        }
        // Performance options
        else if (arg == "--buffer-size" && i + 1 < argc) {
            buffer_size = std::stoi(argv[++i]);
        } else if (arg == "--quiet") {
            quiet_mode = true;
        } else if (arg == "--no-progress") {
            no_progress = true;
        } else if (arg == "--minimal") {
            minimal_output = true;
        }
        // Unknown option
        else {
            std::cerr << "Unknown option: " << arg << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    }

    // Set log level
    if (debug) {
        vep::set_log_level(vep::LogLevel::DEBUG);
    }

    // Validate required arguments
    if (gtf_path.empty() || fasta_path.empty()) {
        std::cerr << "Error: Both --gtf and --fasta are required.\n" << std::endl;
        print_usage(argv[0]);
        return 1;
    }

    if (variant.empty() && vcf_path.empty() && hgvs_input.empty()) {
        std::cerr << "Error: Either --variant, --hgvs, or --vcf is required.\n" << std::endl;
        print_usage(argv[0]);
        return 1;
    }

    // Load gene constraint database if specified
    vep::GeneConstraintDB& constraint_db = vep::get_gene_constraint_db();
    if (!constraint_path.empty()) {
        if (!constraint_db.load_gnomad_constraint(constraint_path)) {
            std::cerr << "Warning: Failed to load constraint file: " << constraint_path << std::endl;
        } else {
            vep::log(vep::LogLevel::INFO, "Loaded " + std::to_string(constraint_db.size()) + " genes from constraint file");
        }
    }
    if (!pli_path.empty()) {
        if (!constraint_db.load_pli_scores(pli_path)) {
            std::cerr << "Warning: Failed to load pLI file: " << pli_path << std::endl;
        }
    }
    if (!loeuf_path.empty()) {
        if (!constraint_db.load_loeuf_scores(loeuf_path)) {
            std::cerr << "Warning: Failed to load LOEUF file: " << loeuf_path << std::endl;
        }
    }

    // Create transcript filter
    vep::TranscriptFilter transcript_filter(filter_config);

    try {
        if (!vcf_path.empty()) {
            // Batch annotation from VCF
            if (output_path.empty()) {
                size_t dot_pos = vcf_path.rfind('.');
                std::string base = (dot_pos != std::string::npos) ?
                    vcf_path.substr(0, dot_pos) : vcf_path;

                if (output_format == vep::OutputFormat::JSON) {
                    output_path = base + "_annotated.json";
                } else if (output_format == vep::OutputFormat::VCF) {
                    output_path = base + "_annotated.vcf";
                } else {
                    output_path = base + "_annotated.tsv";
                }

                if (compress_output) {
                    output_path += ".gz";
                }
            }

            // Create annotator
            vep::VEPAnnotator annotator(gtf_path, fasta_path);

            // Create reference genome for check_ref if needed
            std::unique_ptr<vep::ReferenceGenome> check_ref_genome;
            if (check_ref) {
                check_ref_genome = std::make_unique<vep::ReferenceGenome>(fasta_path, false);
            }

            // Add custom annotation sources
            std::vector<std::string> custom_columns;
            for (size_t idx = 0; idx < annotation_sources.size(); ++idx) {
                const std::string& name = std::get<0>(annotation_sources[idx]);
                const std::string& path = std::get<1>(annotation_sources[idx]);
                const std::string& info_fields = std::get<2>(annotation_sources[idx]);
                bool use_tabix = std::get<3>(annotation_sources[idx]);

                vep::VCFAnnotationConfig config;
                config.name = name;
                config.vcf_path = path;
                config.use_tabix = use_tabix;

                if (!info_fields.empty()) {
                    std::istringstream iss(info_fields);
                    std::string field;
                    while (std::getline(iss, field, ',')) {
                        size_t start = field.find_first_not_of(" \t");
                        size_t end = field.find_last_not_of(" \t");
                        if (start != std::string::npos && end != std::string::npos) {
                            field = field.substr(start, end - start + 1);
                            if (!field.empty()) {
                                config.info_fields.push_back(field);
                                custom_columns.push_back(name + ":" + field);
                            }
                        }
                    }
                }

                annotator.add_annotation_source(config);
            }

            // Add pathogenicity prediction sources
            if (!dbnsfp_path.empty()) {
                auto dbnsfp_source = vep::create_dbnsfp_source(dbnsfp_path, dbnsfp_fields);
                annotator.add_source(dbnsfp_source);
            }

            // Add splice prediction sources
            if (!spliceai_path.empty()) {
                auto spliceai_source = vep::create_spliceai_source(spliceai_path);
                annotator.add_source(spliceai_source);
            }
            if (use_maxentscan) {
                auto mes_source = vep::create_maxentscan_source();
                annotator.add_source(mes_source);
            }
            if (!dbscsnv_path.empty()) {
                auto dbscsnv_source = vep::create_dbscsnv_source(dbscsnv_path);
                annotator.add_source(dbscsnv_source);
            }

            // Add conservation score sources
            if (!phylop_path.empty()) {
                auto phylop_source = vep::create_phylop_source(phylop_path);
                annotator.add_source(phylop_source);
            }
            if (!phastcons_path.empty()) {
                auto phastcons_source = vep::create_phastcons_source(phastcons_path);
                annotator.add_source(phastcons_source);
            }
            if (!gerp_path.empty()) {
                auto gerp_source = vep::create_gerp_source(gerp_path);
                annotator.add_source(gerp_source);
            }

            // Add regulatory annotation source
            if (!regulatory_path.empty()) {
                auto regulatory_source = vep::create_regulatory_source(regulatory_path, std::set<std::string>());
                annotator.add_source(regulatory_source);
            }

            // Add protein domain sources
            if (!pfam_path.empty()) {
                auto pfam_source = vep::create_pfam_source(pfam_path);
                annotator.add_source(pfam_source);
            }
            if (!interpro_path.empty()) {
                auto interpro_source = vep::create_interpro_source(interpro_path);
                annotator.add_source(interpro_source);
            }

            // Add LoF annotation sources
            if (use_loftee) {
                auto loftee_source = vep::create_loftee_source();
                annotator.add_source(loftee_source);
            }
            if (use_nmd) {
                auto nmd_source = vep::create_nmd_source();
                annotator.add_source(nmd_source);
            }
            if (!loftool_path.empty()) {
                auto loftool_source = vep::create_loftool_source(loftool_path);
                annotator.add_source(loftool_source);
            }

            // Initialize all annotation sources
            annotator.initialize_sources();

            // Create output writer
            std::unique_ptr<vep::OutputWriter> writer =
                vep::create_output_writer(output_path, output_format, compress_output);

            writer->write_header(custom_columns);

            // Read and annotate VCF
            bool is_gzipped = (vcf_path.size() > 3 &&
                              vcf_path.substr(vcf_path.size() - 3) == ".gz");

            gzFile gz_file = nullptr;
            std::ifstream plain_file;

            if (is_gzipped) {
                gz_file = gzopen(vcf_path.c_str(), "rb");
                if (!gz_file) {
                    throw std::runtime_error("Cannot open gzipped VCF file: " + vcf_path);
                }
            } else {
                plain_file.open(vcf_path.c_str());
                if (!plain_file.is_open()) {
                    throw std::runtime_error("Cannot open VCF file: " + vcf_path);
                }
            }

            std::string line;
            int variant_count = 0;
            char gz_buffer[65536];

            while (true) {
                // Read next line
                if (is_gzipped) {
                    if (gzgets(gz_file, gz_buffer, sizeof(gz_buffer)) == nullptr) {
                        break;
                    }
                    line = gz_buffer;
                    while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
                        line.pop_back();
                    }
                } else {
                    if (!std::getline(plain_file, line)) {
                        break;
                    }
                }

                if (line.empty() || line[0] == '#') continue;

                std::istringstream iss(line);
                std::string chrom, id, ref, alt, qual, filter_str, info;
                int pos;

                if (!(iss >> chrom >> pos >> id >> ref >> alt)) continue;

                // Check reference allele if requested
                if (check_ref && check_ref_genome) {
                    if (!verify_reference(*check_ref_genome, chrom, pos, ref, quiet_mode)) {
                        continue; // Skip variant with mismatched reference
                    }
                }

                // Skip non-variant sites unless allowed
                if (!allow_non_variant && ref == alt) {
                    continue;
                }

                // Handle multiple alt alleles
                std::istringstream alt_iss(alt);
                std::string single_alt;

                while (std::getline(alt_iss, single_alt, ',')) {
                    // Skip if alt equals ref (non-variant) unless allowed
                    if (!allow_non_variant && single_alt == ref) {
                        continue;
                    }

                    std::vector<vep::VariantAnnotation> annotations;

                    if (all_transcripts) {
                        annotations = annotator.annotate(chrom, pos, ref, single_alt);
                    } else {
                        auto ann = annotator.annotate_most_severe(chrom, pos, ref, single_alt);
                        annotations.push_back(ann);
                    }

                    // Apply transcript filter
                    if (filter_config.pick || filter_config.pick_allele ||
                        filter_config.pick_allele_gene || filter_config.per_gene ||
                        filter_config.most_severe || filter_config.flag_pick ||
                        filter_config.flag_pick_allele || filter_config.flag_pick_allele_gene ||
                        filter_config.canonical_only || filter_config.mane_only ||
                        filter_config.coding_only || filter_config.gencode_basic ||
                        !filter_config.biotypes.empty() || filter_config.no_intergenic ||
                        filter_config.check_frequency) {
                        annotations = transcript_filter.filter(annotations);
                    }

                    // Write annotations
                    for (const auto& ann : annotations) {
                        writer->write_annotation(ann, custom_columns);
                    }

                    // For VCF output, flush after each variant position
                    if (output_format == vep::OutputFormat::VCF) {
                        vep::VCFWriter* vcf_writer = dynamic_cast<vep::VCFWriter*>(writer.get());
                        if (vcf_writer) {
                            vcf_writer->flush_variant();
                        }
                    }

                    variant_count++;
                }

                if (!quiet_mode && !no_progress && variant_count % 1000 == 0) {
                    vep::log(vep::LogLevel::INFO, "Processed " + std::to_string(variant_count) + " variants...");
                }
            }

            // Cleanup
            if (gz_file) {
                gzclose(gz_file);
            }

            writer->write_footer();
            writer->close();

            if (!quiet_mode) {
                std::cout << "Annotation complete. " << variant_count << " variants written to: " << output_path << std::endl;
            }

            // Print statistics if requested
            if (show_stats && !quiet_mode) {
                std::cout << "\n" << writer->get_stats().to_string() << std::endl;
            }

        } else {
            // Single variant annotation
            std::string chrom, ref, alt;
            int pos;
            bool parsed = false;

            // Try HGVS notation first if provided
            if (!hgvs_input.empty()) {
                if (vep::is_hgvs_notation(hgvs_input)) {
                    vep::HGVSParseResult hgvs_result = vep::parse_hgvs(hgvs_input);
                    if (hgvs_result.valid) {
                        if (hgvs_result.hgvs_type == vep::HGVSType::GENOMIC) {
                            // Genomic HGVS can be directly converted
                            chrom = vep::refseq_to_chromosome(hgvs_result.reference_id);
                            if (chrom.empty()) chrom = hgvs_result.reference_id;
                            pos = hgvs_result.start_pos;
                            ref = hgvs_result.ref_allele;
                            alt = hgvs_result.alt_allele;
                            parsed = true;
                            std::cout << "Parsed HGVS: " << chrom << ":" << pos << ":" << ref << ">" << alt << std::endl;
                        } else {
                            // Coding/protein HGVS requires transcript mapping
                            std::cerr << "Note: Coding/protein HGVS requires transcript lookup (not implemented in standalone mode).\n"
                                      << "Please use genomic HGVS (g.) or CHR:POS:REF:ALT format." << std::endl;
                            return 1;
                        }
                    } else {
                        std::cerr << "Error: Invalid HGVS notation: " << hgvs_result.error_message << std::endl;
                        return 1;
                    }
                } else {
                    std::cerr << "Error: Invalid HGVS notation format: " << hgvs_input << std::endl;
                    return 1;
                }
            }

            // Try regular variant format
            if (!parsed && !variant.empty()) {
                if (!parse_variant(variant, chrom, pos, ref, alt)) {
                    std::cerr << "Error: Invalid variant format. Use CHR:POS:REF:ALT\n"
                              << "Example: 7:140753336:A:T" << std::endl;
                    return 1;
                }
                parsed = true;
            }

            if (!parsed) {
                std::cerr << "Error: No variant specified." << std::endl;
                return 1;
            }

            vep::VEPAnnotator annotator(gtf_path, fasta_path);

            // Add custom annotation sources
            for (size_t idx = 0; idx < annotation_sources.size(); ++idx) {
                const std::string& name = std::get<0>(annotation_sources[idx]);
                const std::string& path = std::get<1>(annotation_sources[idx]);
                const std::string& fields = std::get<2>(annotation_sources[idx]);
                bool use_tabix = std::get<3>(annotation_sources[idx]);

                vep::VCFAnnotationConfig config;
                config.name = name;
                config.vcf_path = path;
                config.use_tabix = use_tabix;

                if (!fields.empty()) {
                    std::istringstream iss(fields);
                    std::string field;
                    while (std::getline(iss, field, ',')) {
                        size_t start = field.find_first_not_of(" \t");
                        size_t end = field.find_last_not_of(" \t");
                        if (start != std::string::npos && end != std::string::npos) {
                            field = field.substr(start, end - start + 1);
                            if (!field.empty()) {
                                config.info_fields.push_back(field);
                            }
                        }
                    }
                }

                annotator.add_annotation_source(config);
            }

            // Add pathogenicity prediction sources
            if (!dbnsfp_path.empty()) {
                auto dbnsfp_source = vep::create_dbnsfp_source(dbnsfp_path, dbnsfp_fields);
                annotator.add_source(dbnsfp_source);
            }

            // Add splice prediction sources
            if (!spliceai_path.empty()) {
                auto spliceai_source = vep::create_spliceai_source(spliceai_path);
                annotator.add_source(spliceai_source);
            }
            if (use_maxentscan) {
                auto mes_source = vep::create_maxentscan_source();
                annotator.add_source(mes_source);
            }
            if (!dbscsnv_path.empty()) {
                auto dbscsnv_source = vep::create_dbscsnv_source(dbscsnv_path);
                annotator.add_source(dbscsnv_source);
            }

            // Add conservation score sources
            if (!phylop_path.empty()) {
                auto phylop_source = vep::create_phylop_source(phylop_path);
                annotator.add_source(phylop_source);
            }
            if (!phastcons_path.empty()) {
                auto phastcons_source = vep::create_phastcons_source(phastcons_path);
                annotator.add_source(phastcons_source);
            }
            if (!gerp_path.empty()) {
                auto gerp_source = vep::create_gerp_source(gerp_path);
                annotator.add_source(gerp_source);
            }

            // Add regulatory annotation source
            if (!regulatory_path.empty()) {
                auto regulatory_source = vep::create_regulatory_source(regulatory_path, std::set<std::string>());
                annotator.add_source(regulatory_source);
            }

            // Add protein domain sources
            if (!pfam_path.empty()) {
                auto pfam_source = vep::create_pfam_source(pfam_path);
                annotator.add_source(pfam_source);
            }
            if (!interpro_path.empty()) {
                auto interpro_source = vep::create_interpro_source(interpro_path);
                annotator.add_source(interpro_source);
            }

            // Add LoF annotation sources
            if (use_loftee) {
                auto loftee_source = vep::create_loftee_source();
                annotator.add_source(loftee_source);
            }
            if (use_nmd) {
                auto nmd_source = vep::create_nmd_source();
                annotator.add_source(nmd_source);
            }
            if (!loftool_path.empty()) {
                auto loftool_source = vep::create_loftool_source(loftool_path);
                annotator.add_source(loftool_source);
            }

            // Load plugins
            vep::PluginLoader plugin_loader;
            for (size_t i = 0; i < plugin_dirs.size(); ++i) {
                plugin_loader.add_plugin_dir(plugin_dirs[i]);
            }
            plugin_loader.load_plugins_from_dirs();

            for (size_t i = 0; i < plugins.size(); ++i) {
                const std::string& path = plugins[i].first;
                const std::string& config = plugins[i].second;
                if (!plugin_loader.load_plugin(path, config)) {
                    std::cerr << "Warning: Failed to load plugin: " << path << std::endl;
                }
            }

            // Add plugin sources
            auto plugin_sources = plugin_loader.get_all_sources();
            for (size_t i = 0; i < plugin_sources.size(); ++i) {
                annotator.add_source(plugin_sources[i]);
            }

            // Initialize all annotation sources
            annotator.initialize_sources();

            std::cout << "\nDatabase loaded:\n" << annotator.get_stats() << std::endl;

            std::vector<vep::VariantAnnotation> annotations;

            if (all_transcripts) {
                annotations = annotator.annotate(chrom, pos, ref, alt);
            } else {
                auto ann = annotator.annotate_most_severe(chrom, pos, ref, alt);
                annotations.push_back(ann);
            }

            // Apply transcript filter
            if (filter_config.pick || filter_config.pick_allele ||
                filter_config.pick_allele_gene || filter_config.per_gene ||
                filter_config.most_severe || filter_config.flag_pick ||
                filter_config.flag_pick_allele || filter_config.flag_pick_allele_gene ||
                filter_config.canonical_only || filter_config.mane_only ||
                filter_config.coding_only || filter_config.gencode_basic ||
                !filter_config.biotypes.empty() || filter_config.no_intergenic ||
                filter_config.check_frequency) {
                annotations = transcript_filter.filter(annotations);
            }

            if (!quiet_mode) {
                std::cout << "\nFound " << annotations.size() << " annotation(s):" << std::endl;
            }

            // Output HGVSg and SPDI if requested
            if (include_hgvsg || include_spdi) {
                std::cout << "\nAdditional identifiers:" << std::endl;
                if (include_hgvsg) {
                    std::string hgvsg = vep::generate_hgvsg(chrom, pos, ref, alt);
                    std::cout << "  HGVSg: " << hgvsg << std::endl;
                }
                if (include_spdi) {
                    std::string spdi = vep::generate_spdi(chrom, pos, ref, alt);
                    std::cout << "  SPDI:  " << spdi << std::endl;
                }
            }

            // Output gene constraint info if available
            if (constraint_db.is_loaded() && !annotations.empty()) {
                std::set<std::string> genes_shown;
                std::cout << "\nGene constraint scores:" << std::endl;
                for (size_t i = 0; i < annotations.size(); ++i) {
                    const std::string& gene = annotations[i].gene_symbol;
                    if (!gene.empty() && genes_shown.count(gene) == 0) {
                        vep::GeneConstraint gc = constraint_db.get_by_symbol(gene);
                        if (gc.has_data()) {
                            std::cout << "  " << gene << ": pLI=" << vep::format_constraint_score(gc.pLI)
                                      << " LOEUF=" << vep::format_constraint_score(gc.oe_lof_upper)
                                      << " (" << gc.get_constraint_level() << ")" << std::endl;
                        }
                        genes_shown.insert(gene);
                    }
                }
            }

            for (size_t i = 0; i < annotations.size(); ++i) {
                print_annotation(annotations[i]);
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
