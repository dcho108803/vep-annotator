/**
 * VEP Variant Annotator - Main Entry Point
 *
 * Pure C++ local implementation for annotating genetic variants.
 * Requires local GTF and FASTA files - no external API calls.
 */

#include "vep_annotator.hpp"
#include "annotation_sources.hpp"
#include "plugin.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>

void print_usage(const char* program_name) {
    std::cout << "VEP Variant Annotator - Pure C++ Local Implementation\n"
              << "======================================================\n\n"
              << "Usage: " << program_name << " [OPTIONS]\n\n"
              << "Required Data Files:\n"
              << "  --gtf FILE              GTF/GFF annotation file (genes/transcripts)\n"
              << "  --fasta FILE            Reference genome FASTA file\n\n"
              << "Variant Input (choose one):\n"
              << "  -v, --variant CHR:POS:REF:ALT   Single variant to annotate\n"
              << "  --vcf FILE              VCF file for batch annotation\n\n"
              << "Custom Annotations (can be used multiple times):\n"
              << "  --annotation NAME:VCF_PATH[:FIELDS]\n"
              << "                          Add custom VCF annotation source (loads into memory)\n"
              << "                          NAME: source name (e.g., gnomad, clinvar)\n"
              << "                          VCF_PATH: path to VCF file\n"
              << "                          FIELDS: comma-separated INFO fields (optional)\n\n"
              << "  --annotation-tabix NAME:VCF_PATH[:FIELDS]\n"
              << "                          Add tabix-indexed VCF source (on-disk queries)\n"
              << "                          Use for large files like gnomAD (requires htslib)\n"
              << "                          VCF must be bgzip compressed with .tbi index\n\n"
              << "Pathogenicity Predictions:\n"
              << "  --dbnsfp FILE           dbNSFP database (tabix-indexed .txt.gz)\n"
              << "  --dbnsfp-fields FIELDS  Fields to extract (default: essential)\n"
              << "                          Presets: essential, pathogenicity, conservation, all\n"
              << "                          Or comma-separated: SIFT_score,CADD_phred,REVEL_score\n\n"
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
              << "Plugins:\n"
              << "  --plugin PATH[:CONFIG]  Load plugin from shared library\n"
              << "                          CONFIG: key1=value1;key2=value2\n"
              << "  --plugin-dir DIR        Directory to search for plugins\n\n"
              << "Output Options:\n"
              << "  -o, --output FILE       Output file path (for VCF annotation)\n"
              << "  --all-transcripts       Output all transcript annotations (default: most severe)\n\n"
              << "Other Options:\n"
              << "  -h, --help              Show this help message\n"
              << "  --debug                 Enable debug logging\n\n"
              << "Examples:\n"
              << "  # Annotate a single SNV\n"
              << "  " << program_name << " --gtf genes.gtf --fasta genome.fa -v 7:140753336:A:T\n\n"
              << "  # Annotate with gnomAD frequencies (in-memory, for small files)\n"
              << "  " << program_name << " --gtf genes.gtf --fasta genome.fa \\\n"
              << "      --annotation gnomad:gnomad.vcf.gz:AF,AF_popmax \\\n"
              << "      -v 7:140753336:A:T\n\n"
              << "  # Annotate with gnomAD using tabix (on-disk, for large files)\n"
              << "  " << program_name << " --gtf genes.gtf --fasta genome.fa \\\n"
              << "      --annotation-tabix gnomad:gnomad.exomes.vcf.bgz:AF,AF_popmax \\\n"
              << "      -v 7:140753336:A:T\n\n"
              << "  # Annotate VCF with multiple annotation sources\n"
              << "  " << program_name << " --gtf genes.gtf --fasta genome.fa \\\n"
              << "      --annotation-tabix gnomad:gnomad.vcf.bgz:AF \\\n"
              << "      --annotation clinvar:clinvar.vcf.gz:CLNSIG,CLNDN \\\n"
              << "      --vcf variants.vcf -o results.tsv\n\n"
              << "Data Files:\n"
              << "  GTF: https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/\n"
              << "  FASTA: https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/\n"
              << "  gnomAD: https://gnomad.broadinstitute.org/downloads\n"
              << "  ClinVar: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/\n"
              << std::endl;
}

void print_annotation(const vep::VariantAnnotation& ann) {
    std::cout << "\n=== Variant Annotation ===" << std::endl;

    auto summary = ann.get_summary();
    for (const auto& [key, value] : summary) {
        if (!value.empty()) {
            std::cout << key << ": " << value << std::endl;
        }
    }
}

// Parse variant string in format CHR:POS:REF:ALT or CHR-POS-REF-ALT
bool parse_variant(const std::string& variant, std::string& chrom, int& pos,
                   std::string& ref, std::string& alt) {
    std::string normalized = variant;

    // Replace - with : for uniform parsing
    for (char& c : normalized) {
        if (c == '-') c = ':';
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
        // No fields specified
        vcf_path = source.substr(p1 + 1);
        fields = "";
    } else {
        vcf_path = source.substr(p1 + 1, p2 - p1 - 1);
        fields = source.substr(p2 + 1);
    }

    return !name.empty() && !vcf_path.empty();
}

int main(int argc, char* argv[]) {
    std::string gtf_path;
    std::string fasta_path;
    std::string variant;
    std::string vcf_path;
    std::string output_path;
    bool all_transcripts = false;
    bool debug = false;

    // Custom annotation sources: (name, vcf_path, fields, use_tabix)
    std::vector<std::tuple<std::string, std::string, std::string, bool>> annotation_sources;

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

    // Plugins
    std::vector<std::pair<std::string, std::string>> plugins;  // (path, config)
    std::vector<std::string> plugin_dirs;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--gtf" && i + 1 < argc) {
            gtf_path = argv[++i];
        } else if (arg == "--fasta" && i + 1 < argc) {
            fasta_path = argv[++i];
        } else if ((arg == "-v" || arg == "--variant") && i + 1 < argc) {
            variant = argv[++i];
        } else if (arg == "--vcf" && i + 1 < argc) {
            vcf_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            output_path = argv[++i];
        } else if (arg == "--annotation" && i + 1 < argc) {
            std::string source_str = argv[++i];
            std::string name, path, fields;
            if (parse_annotation_source(source_str, name, path, fields)) {
                annotation_sources.emplace_back(name, path, fields, false);  // use_tabix = false
            } else {
                std::cerr << "Error: Invalid annotation format: " << source_str << "\n"
                          << "Expected: NAME:VCF_PATH[:FIELDS]" << std::endl;
                return 1;
            }
        } else if (arg == "--annotation-tabix" && i + 1 < argc) {
            std::string source_str = argv[++i];
            std::string name, path, fields;
            if (parse_annotation_source(source_str, name, path, fields)) {
                annotation_sources.emplace_back(name, path, fields, true);  // use_tabix = true
            } else {
                std::cerr << "Error: Invalid annotation format: " << source_str << "\n"
                          << "Expected: NAME:VCF_PATH[:FIELDS]" << std::endl;
                return 1;
            }
        } else if (arg == "--all-transcripts") {
            all_transcripts = true;
        } else if (arg == "--debug") {
            debug = true;
        } else if (arg == "--dbnsfp" && i + 1 < argc) {
            dbnsfp_path = argv[++i];
        } else if (arg == "--dbnsfp-fields" && i + 1 < argc) {
            dbnsfp_fields = argv[++i];
        } else if (arg == "--spliceai" && i + 1 < argc) {
            spliceai_path = argv[++i];
        } else if (arg == "--maxentscan") {
            use_maxentscan = true;
        } else if (arg == "--dbscsnv" && i + 1 < argc) {
            dbscsnv_path = argv[++i];
        } else if (arg == "--phylop" && i + 1 < argc) {
            phylop_path = argv[++i];
        } else if (arg == "--phastcons" && i + 1 < argc) {
            phastcons_path = argv[++i];
        } else if (arg == "--gerp" && i + 1 < argc) {
            gerp_path = argv[++i];
        } else if (arg == "--regulatory" && i + 1 < argc) {
            regulatory_path = argv[++i];
        } else if (arg == "--pfam" && i + 1 < argc) {
            pfam_path = argv[++i];
        } else if (arg == "--interpro" && i + 1 < argc) {
            interpro_path = argv[++i];
        } else if (arg == "--loftee") {
            use_loftee = true;
        } else if (arg == "--nmd") {
            use_nmd = true;
        } else if (arg == "--loftool" && i + 1 < argc) {
            loftool_path = argv[++i];
        } else if (arg == "--plugin" && i + 1 < argc) {
            std::string plugin_arg = argv[++i];
            size_t colon = plugin_arg.find(':');
            if (colon != std::string::npos) {
                plugins.emplace_back(plugin_arg.substr(0, colon), plugin_arg.substr(colon + 1));
            } else {
                plugins.emplace_back(plugin_arg, "");
            }
        } else if (arg == "--plugin-dir" && i + 1 < argc) {
            plugin_dirs.push_back(argv[++i]);
        } else {
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

    if (variant.empty() && vcf_path.empty()) {
        std::cerr << "Error: Either --variant or --vcf is required.\n" << std::endl;
        print_usage(argv[0]);
        return 1;
    }

    try {
        if (!vcf_path.empty()) {
            // Batch annotation from VCF
            if (output_path.empty()) {
                size_t dot_pos = vcf_path.rfind('.');
                if (dot_pos != std::string::npos) {
                    output_path = vcf_path.substr(0, dot_pos) + "_annotated.tsv";
                } else {
                    output_path = vcf_path + "_annotated.tsv";
                }
            }

            vep::annotate_vcf_file(vcf_path, output_path, gtf_path, fasta_path, annotation_sources);
            std::cout << "Annotation complete. Results saved to: " << output_path << std::endl;

        } else {
            // Single variant annotation
            std::string chrom, ref, alt;
            int pos;

            if (!parse_variant(variant, chrom, pos, ref, alt)) {
                std::cerr << "Error: Invalid variant format. Use CHR:POS:REF:ALT\n"
                          << "Example: 7:140753336:A:T" << std::endl;
                return 1;
            }

            vep::VEPAnnotator annotator(gtf_path, fasta_path);

            // Add custom annotation sources
            for (const auto& [name, path, fields, use_tabix] : annotation_sources) {
                vep::VCFAnnotationConfig config;
                config.name = name;
                config.vcf_path = path;
                config.use_tabix = use_tabix;

                // Parse fields
                if (!fields.empty()) {
                    std::istringstream iss(fields);
                    std::string field;
                    while (std::getline(iss, field, ',')) {
                        field.erase(0, field.find_first_not_of(" \t"));
                        field.erase(field.find_last_not_of(" \t") + 1);
                        if (!field.empty()) {
                            config.info_fields.push_back(field);
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
                auto regulatory_source = vep::create_regulatory_source(regulatory_path, {});
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
            for (const auto& dir : plugin_dirs) {
                plugin_loader.add_plugin_dir(dir);
            }
            plugin_loader.load_plugins_from_dirs();

            for (const auto& [path, config] : plugins) {
                if (!plugin_loader.load_plugin(path, config)) {
                    std::cerr << "Warning: Failed to load plugin: " << path << std::endl;
                }
            }

            // Add plugin sources
            auto plugin_sources = plugin_loader.get_all_sources();
            for (auto& source : plugin_sources) {
                annotator.add_source(source);
            }

            // Initialize all annotation sources
            annotator.initialize_sources();

            std::cout << "\nDatabase loaded:\n" << annotator.get_stats() << std::endl;

            if (all_transcripts) {
                auto annotations = annotator.annotate(chrom, pos, ref, alt);

                std::cout << "\nFound " << annotations.size() << " transcript annotations:" << std::endl;

                for (const auto& ann : annotations) {
                    print_annotation(ann);
                }
            } else {
                auto annotation = annotator.annotate_most_severe(chrom, pos, ref, alt);
                print_annotation(annotation);
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
