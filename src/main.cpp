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
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <set>
#include <algorithm>
#include <cctype>
#include <memory>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <queue>
#include <chrono>
#include <cstring>
#include <unordered_map>

void print_usage(const char* program_name) {
    std::cout << "VEP Variant Annotator - Pure C++ Local Implementation\n"
              << "======================================================\n\n"
              << "Usage: " << program_name << " [OPTIONS]\n\n"
              << "Required Data Files:\n"
              << "  --gtf FILE              GTF/GFF annotation file (genes/transcripts)\n"
              << "  --fasta FILE            Reference genome FASTA file\n\n"
              << "Variant Input (choose one):\n"
              << "  -v, --variant CHR:POS:REF:ALT   Single variant to annotate\n"
              << "  --hgvs [NOTATION]       Without argument: enable HGVSc/HGVSp in output\n"
              << "                          With argument: HGVS input notation\n"
              << "  --vcf FILE              VCF file for batch annotation\n"
              << "  -i, --input_file FILE   Same as --vcf\n\n"
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
              << "  --canonical             Add CANONICAL=YES flag to output (display flag)\n"
              << "  --canonical-only        Only annotate canonical transcripts (filter)\n"
              << "  --mane                  Add MANE_SELECT info to output (display flag)\n"
              << "  --mane-select           Alias for --mane\n"
              << "  --mane_select           Alias for --mane\n"
              << "  --mane-only             Only annotate MANE Select transcripts (filter)\n"
              << "  --coding-only           Only annotate protein-coding transcripts\n"
              << "  --biotype LIST          Only annotate these biotypes (comma-separated)\n"
              << "  --no-intergenic         Skip intergenic variants\n\n"
              << "Frequency Filtering:\n"
              << "  --check-frequency       Enable frequency-based filtering\n"
              << "  --freq-pop FIELD        Frequency field to check (e.g., gnomAD:AF)\n"
              << "  --freq-threshold VAL    Maximum frequency threshold (default: 0.01)\n"
              << "  --freq-gt               Filter for frequency greater than threshold\n\n"
              << "Custom Annotations (can be used multiple times):\n"
              << "  --annotation NAME:VCF_PATH[:FIELDS]\n"
              << "                          Add custom VCF annotation source (loads into memory)\n"
              << "  --annotation-tabix NAME:VCF_PATH[:FIELDS]\n"
              << "                          Add tabix-indexed VCF source (on-disk queries)\n"
              << "  --custom FILE,NAME,TYPE,OVERLAP,0,FIELDS\n"
              << "                          Perl VEP --custom format (TYPE: vcf, bed, bigwig)\n\n"
              << "Pathogenicity Predictions:\n"
              << "  --dbnsfp FILE           dbNSFP database (tabix-indexed .txt.gz)\n"
              << "  --dbnsfp-fields FIELDS  Fields to extract (default: essential)\n"
              << "                          Presets: essential, pathogenicity, conservation, all\n\n"
              << "Splice Predictions:\n"
              << "  --spliceai FILE         SpliceAI VCF file (tabix-indexed)\n"
              << "  --spliceai-snv FILE     SpliceAI SNV VCF file (tabix-indexed)\n"
              << "  --spliceai-indel FILE   SpliceAI indel VCF file (tabix-indexed)\n"
              << "  --spliceai-cutoff VAL   SpliceAI delta score cutoff (adds PASS/FAIL)\n"
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
              << "  --xref-refseq           Add RefSeq cross-references\n"
              << "  --transcript-version    Append version to transcript identifiers\n\n"
              << "Display Flags:\n"
              << "  --show-canonical        Show CANONICAL flag in output\n"
              << "  --show-mane             Show MANE_SELECT flag in output\n"
              << "  --show-tsl              Show Transcript Support Level\n"
              << "  --show-appris           Show APPRIS annotation\n"
              << "  --show-biotype          Show BIOTYPE in custom columns\n"
              << "  --variant-class         Add SO variant class (SNV, insertion, etc.)\n"
              << "  --nearest               Add nearest gene for intergenic variants\n\n"
              << "Distance and Shifting:\n"
              << "  --distance VAL[,VAL]    Upstream[,downstream] distance (default: 5000)\n"
              << "  --shift-3prime          3' shift indels before consequence determination\n"
              << "  --shift-genomic         Normalize at genomic level\n\n"
              << "VCF Output Options:\n"
              << "  --vcf-info-field NAME   Custom VCF INFO field name (default: CSQ)\n"
              << "  --fields FIELDS         Comma-separated list of CSQ fields to output\n"
              << "  --allele-number         Add ALLELE_NUM field to output\n"
              << "  --no-headers            Suppress output headers\n"
              << "  --keep-csq              Preserve existing CSQ field from input VCF\n"
              << "  --no-escape             Don't URI-encode special characters in VCF output\n"
              << "  --terms STYLE           Consequence term style: SO (default), display\n\n"
              << "Regulatory Options:\n"
              << "  --cell-type LIST        Comma-separated cell types for regulatory\n\n"
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
              << "  --exclude-predicted     Exclude predicted (XM_/XR_) transcripts\n"
              << "  --chr LIST              Only annotate these chromosomes (comma-separated)\n"
              << "  --overlaps TYPE         SV overlap type: any, within, surrounding, exact\n"
              << "  --include_consequence LIST  Only include these consequences (comma-separated)\n"
              << "  --exclude_consequence LIST  Exclude these consequences (comma-separated)\n\n"
              << "Perl VEP Compatibility:\n"
              << "  All Perl VEP underscore flags (e.g., --pick_allele) are accepted as\n"
              << "  hyphenated equivalents (--pick-allele). Common aliases:\n"
              << "  --tab                   Same as --output-format tsv\n"
              << "  --json                  Same as --output-format json\n"
              << "  -i, --input-file FILE   Same as --vcf FILE (also accepts STDIN)\n"
              << "  --format FORMAT         Input format (auto-detected)\n"
              << "  --force-overwrite       No-op (C++ always overwrites)\n"
              << "  --symbol                Show gene symbol in output (always shown)\n"
              << "  --protein               Show protein ID in output\n"
              << "  --domains               Show domain annotations (use --pfam/--interpro)\n"
              << "  --sift [p|s|b]          SIFT predictions (p=prediction, s=score, b=both)\n"
              << "  --polyphen [p|s|b]      PolyPhen predictions\n"
              << "  --check-existing FILE   Check for co-located known variants (VCF)\n"
              << "  --show-ref-allele       Show GIVEN_REF/USED_REF in output\n"
              << "  --fork N                Use N parallel annotation threads\n"
              << "  --cache, --offline, --database, etc. accepted as no-ops\n\n"
              << "Performance Options:\n"
              << "  --buffer-size N         Number of variants to buffer (default: 5000)\n"
              << "  --quiet                 Suppress progress messages\n"
              << "  --no-progress           Don't show progress bar\n"
              << "  --minimal               Only output most essential fields\n\n"
              << "Other Options:\n"
              << "  -h, --help              Show this help message\n"
              << "  --debug                 Enable debug logging\n"
              << "  --config FILE           Load options from config file (key=value format)\n"
              << "  --synonyms FILE         Chromosome synonym mapping file (tab-separated)\n"
              << "  --stats-file FILE       Write run statistics to file\n\n"
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

// Post-process annotations: add max AF, total length, identifier fields
void post_process_annotations(std::vector<vep::VariantAnnotation>& annotations,
                               bool show_max_af, bool include_total_length,
                               bool include_ccds, bool include_uniprot,
                               bool include_xref_refseq,
                               const std::string& nearest_mode = "symbol") {
    for (auto& ann : annotations) {
        // Compute MAX_AF from all AF-like fields in custom_annotations
        if (show_max_af) {
            double max_af = 0.0;
            std::string max_af_pop;
            for (const auto& kv : ann.custom_annotations) {
                // Match fields ending with ":AF" or containing "AF_"
                const std::string& key = kv.first;
                if (key.size() >= 3 &&
                    (key.substr(key.size() - 3) == ":AF" ||
                     key.find(":AF_") != std::string::npos)) {
                    try {
                        double af = std::stod(kv.second);
                        if (af > max_af) {
                            max_af = af;
                            // Extract population from field name
                            size_t colon_pos = key.find(':');
                            if (colon_pos != std::string::npos) {
                                max_af_pop = key.substr(colon_pos + 1);
                            }
                        }
                    } catch (...) {}
                }
            }
            if (max_af > 0.0) {
                ann.custom_annotations["MAX_AF"] = std::to_string(max_af);
                if (!max_af_pop.empty()) {
                    ann.custom_annotations["MAX_AF_POP"] = max_af_pop;
                }
            }
        }

        // Add total transcript/CDS length
        if (include_total_length) {
            auto it_tl = ann.custom_annotations.find("TRANSCRIPT_LENGTH");
            auto it_cl = ann.custom_annotations.find("CDS_LENGTH");
            // These are populated by annotate_transcript() if include_total_length
            // Just ensure they exist as empty strings if not set
            if (it_tl == ann.custom_annotations.end()) {
                ann.custom_annotations["TRANSCRIPT_LENGTH"] = "";
            }
            if (it_cl == ann.custom_annotations.end()) {
                ann.custom_annotations["CDS_LENGTH"] = "";
            }
        }

        // Identifier fields are already populated by annotate_transcript()
        // via ccds_id and protein_id on the Transcript struct.
        // Just ensure columns exist in output if flags are set.
        if (include_ccds) {
            if (ann.custom_annotations.find("CCDS") == ann.custom_annotations.end()) {
                ann.custom_annotations["CCDS"] = "";
            }
        }
        if (include_uniprot) {
            if (ann.custom_annotations.find("UNIPROT") == ann.custom_annotations.end()) {
                ann.custom_annotations["UNIPROT"] = "";
            }
        }
        if (include_xref_refseq) {
            if (ann.custom_annotations.find("RefSeq") == ann.custom_annotations.end()) {
                ann.custom_annotations["RefSeq"] = "";
            }
        }

        // Resolve --nearest mode: symbol (default), gene (ID), transcript
        if (nearest_mode != "symbol") {
            auto nit = ann.custom_annotations.find("NEAREST");
            if (nit != ann.custom_annotations.end() && !nit->second.empty()) {
                if (nearest_mode == "gene" && !ann.gene_id.empty()) {
                    nit->second = ann.gene_id;
                } else if (nearest_mode == "transcript" && !ann.transcript_id.empty()) {
                    nit->second = ann.transcript_id;
                }
            }
        }
    }
}

// Print summary of consequence counts
void print_summary(const vep::AnnotationStats& stats) {
    std::cout << "\n=== Consequence Summary ===" << std::endl;
    std::cout << "Total annotations: " << stats.total_variants << std::endl;
    std::cout << "Annotated: " << stats.annotated_variants << std::endl;
    std::cout << "\nBy consequence:" << std::endl;
    for (const auto& pair : stats.consequence_counts) {
        std::cout << "  " << pair.first << ": " << pair.second << std::endl;
    }
    std::cout << "\nBy impact:" << std::endl;
    for (const auto& pair : stats.impact_counts) {
        std::cout << "  " << pair.first << ": " << pair.second << std::endl;
    }
    if (!stats.biotype_counts.empty()) {
        std::cout << "\nBy biotype:" << std::endl;
        for (const auto& pair : stats.biotype_counts) {
            std::cout << "  " << pair.first << ": " << pair.second << std::endl;
        }
    }
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
    // Try colon-separated first (preferred)
    char sep = ':';
    size_t p1 = variant.find(':');
    if (p1 == std::string::npos) {
        // Fall back to dash separator, but only use it for field separation
        // by finding exactly 3 dashes to split on (CHR-POS-REF-ALT)
        sep = '-';
        p1 = variant.find('-');
        if (p1 == std::string::npos) return false;
    }

    size_t p2 = variant.find(sep, p1 + 1);
    if (p2 == std::string::npos) return false;

    size_t p3 = variant.find(sep, p2 + 1);
    if (p3 == std::string::npos) return false;

    try {
        chrom = variant.substr(0, p1);
        pos = std::stoi(variant.substr(p1 + 1, p2 - p1 - 1));
        ref = variant.substr(p2 + 1, p3 - p2 - 1);
        alt = variant.substr(p3 + 1);
        return true;
    } catch (...) {
        return false;
    }
}

// Input format types for auto-detection
enum class InputFormat { VCF, ENSEMBL, HGVS, BED, REGION, UNKNOWN };

// Detect input format from a non-header line
InputFormat detect_input_format(const std::string& line) {
    if (line.empty() || line[0] == '#') return InputFormat::UNKNOWN;

    // HGVS: contains ':' followed by c./p./g./n./m. notation
    if (line.find(":c.") != std::string::npos || line.find(":p.") != std::string::npos ||
        line.find(":g.") != std::string::npos || line.find(":n.") != std::string::npos ||
        line.find(":m.") != std::string::npos) {
        return InputFormat::HGVS;
    }

    // Region format: CHR:START-END/ALLELE or CHR:START-END:STRAND/ALLELE
    if (line.find('/') != std::string::npos && line.find(':') != std::string::npos &&
        line.find('-') != std::string::npos) {
        // Check for pattern like 1:12345-12345:1/A
        size_t colon1 = line.find(':');
        size_t dash = line.find('-', colon1);
        if (dash != std::string::npos) {
            return InputFormat::REGION;
        }
    }

    // Count tabs to distinguish VCF from Ensembl format
    int tab_count = 0;
    for (char c : line) {
        if (c == '\t') tab_count++;
    }

    // VCF: at least 4 tabs (CHROM\tPOS\tID\tREF\tALT)
    if (tab_count >= 4) return InputFormat::VCF;

    // BED: exactly 2-5 tabs (CHR\tSTART\tEND)
    if (tab_count >= 2 && tab_count <= 5) {
        // Could be BED - check if first 3 fields are chr, int, int
        std::istringstream iss(line);
        std::string f1, f2, f3;
        if (std::getline(iss, f1, '\t') && std::getline(iss, f2, '\t') && std::getline(iss, f3, '\t')) {
            try {
                std::stoi(f2);
                std::stoi(f3);
                return InputFormat::BED;
            } catch (...) {}
        }
    }

    // Ensembl default: space-separated, 5 fields (CHR START END ALLELES STRAND)
    {
        std::istringstream iss(line);
        std::string f1, f2, f3, f4, f5;
        if (iss >> f1 >> f2 >> f3 >> f4 >> f5) {
            try {
                std::stoi(f2);
                std::stoi(f3);
                // f4 contains alleles like A/T, f5 is strand +/-/1/-1
                if (f4.find('/') != std::string::npos) {
                    return InputFormat::ENSEMBL;
                }
            } catch (...) {}
        }
    }

    return InputFormat::UNKNOWN;
}

// Parse Ensembl default format line: CHR START END ALLELES STRAND [IDENTIFIER]
bool parse_ensembl_line(const std::string& line, std::string& chrom, int& pos,
                        std::string& ref, std::string& alt) {
    std::istringstream iss(line);
    std::string start_str, end_str, alleles, strand_str;

    if (!(iss >> chrom >> start_str >> end_str >> alleles >> strand_str)) return false;

    try {
        pos = std::stoi(start_str);
    } catch (...) { return false; }

    // Parse alleles: REF/ALT
    size_t slash = alleles.find('/');
    if (slash == std::string::npos) return false;

    ref = alleles.substr(0, slash);
    alt = alleles.substr(slash + 1);

    // Handle '-' as empty allele (insertion/deletion)
    if (ref == "-") ref = "";
    if (alt == "-") alt = "";

    return true;
}

// Parse BED format line: CHR START END [NAME [SCORE [STRAND]]]
bool parse_bed_line(const std::string& line, std::string& chrom, int& pos,
                    std::string& ref, std::string& alt) {
    std::istringstream iss(line);
    std::string start_str, end_str;

    if (!(iss >> chrom >> start_str >> end_str)) return false;

    try {
        int start = std::stoi(start_str);
        pos = start + 1;  // BED is 0-based, convert to 1-based
        ref = "-";
        alt = "-";
    } catch (...) { return false; }

    return true;
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

/**
 * Simple BED annotation source for --custom BED support.
 * Loads BED intervals into memory, returns overlap annotations.
 */
class BedCustomAnnotation {
public:
    struct BedInterval {
        int start;  // 0-based start
        int end;    // 0-based exclusive end
        std::string name;
        std::string score;
        char strand = '.';
        std::vector<std::string> extra_fields;
    };

    std::string source_name;
    std::map<std::string, std::vector<BedInterval>> intervals; // chrom -> sorted intervals

    bool load(const std::string& path, const std::string& name) {
        source_name = name;
        std::ifstream infile(path);
        if (!infile.is_open()) return false;

        std::string line;
        while (std::getline(infile, line)) {
            if (line.empty() || line[0] == '#') continue;
            if (line.substr(0, 5) == "track" || line.substr(0, 7) == "browser") continue;

            std::istringstream iss(line);
            std::string chrom, start_s, end_s;
            if (!(iss >> chrom >> start_s >> end_s)) continue;

            BedInterval iv;
            try {
                iv.start = std::stoi(start_s);
                iv.end = std::stoi(end_s);
            } catch (...) { continue; }

            std::string field;
            if (iss >> field) iv.name = field;
            if (iss >> field) iv.score = field;
            if (iss >> field) iv.strand = field.empty() ? '.' : field[0];
            while (iss >> field) iv.extra_fields.push_back(field);

            intervals[chrom].push_back(iv);
        }

        // Sort intervals by start position for binary search
        for (auto& [chr, ivs] : intervals) {
            std::sort(ivs.begin(), ivs.end(),
                [](const BedInterval& a, const BedInterval& b) { return a.start < b.start; });
        }
        return true;
    }

    // Find chromosome iterator, trying both chr-prefixed and non-prefixed forms
    decltype(intervals)::const_iterator find_chrom(const std::string& chrom) const {
        auto it = intervals.find(chrom);
        if (it != intervals.end()) return it;
        // Try alternate naming
        if (chrom.length() > 3 && chrom.substr(0, 3) == "chr") {
            return intervals.find(chrom.substr(3));
        } else {
            return intervals.find("chr" + chrom);
        }
    }

    // Query overlapping intervals for a 1-based position
    std::vector<const BedInterval*> query(const std::string& chrom, int pos1) const {
        std::vector<const BedInterval*> results;
        auto it = find_chrom(chrom);
        if (it == intervals.end()) return results;

        int pos0 = pos1 - 1; // Convert to 0-based
        for (const auto& iv : it->second) {
            if (iv.start > pos0) break; // Past the query point
            if (iv.end > pos0) {
                results.push_back(&iv);
            }
        }
        return results;
    }

    // Query overlapping intervals for a 1-based range [start, end]
    std::vector<const BedInterval*> query_range(const std::string& chrom, int start1, int end1) const {
        std::vector<const BedInterval*> results;
        auto it = find_chrom(chrom);
        if (it == intervals.end()) return results;

        int start0 = start1 - 1;
        int end0 = end1; // BED end is exclusive, our end1 is inclusive
        for (const auto& iv : it->second) {
            if (iv.start >= end0) break;
            if (iv.end > start0) {
                results.push_back(&iv);
            }
        }
        return results;
    }

    // Annotate a variant and add results to custom_annotations map
    void annotate(const std::string& chrom, int pos,
                  std::unordered_map<std::string, std::string>& custom_annotations) const {
        auto hits = query(chrom, pos);
        if (!hits.empty()) {
            // Collect names of overlapping regions
            std::string names;
            for (const auto* iv : hits) {
                if (!names.empty()) names += "&";
                names += iv->name.empty() ? "1" : iv->name;
            }
            custom_annotations[source_name] = names;
        }
    }
};

// Read a complete line from gzFile, handling lines longer than a fixed buffer
std::string gz_read_line(gzFile gz, bool& eof) {
    const int CHUNK = 65536;
    char buf[65536];
    std::string result;
    while (true) {
        if (gzgets(gz, buf, CHUNK) == nullptr) {
            eof = result.empty();
            break;
        }
        result += buf;
        // If we got a newline, the line is complete
        if (!result.empty() && result.back() == '\n') break;
        // If gzgets returned less than CHUNK-1 bytes without newline, it's EOF
        size_t len = std::strlen(buf);
        if (static_cast<int>(len) < CHUNK - 1) break;
    }
    // Strip trailing newline/carriage return
    while (!result.empty() && (result.back() == '\n' || result.back() == '\r'))
        result.pop_back();
    return result;
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
    bool include_hgvs_output = false;  // --hgvs as output display flag
    bool show_protein = false;         // --protein display flag
    bool show_gene_symbol = false;     // --symbol display flag (always shown, no-op)

    // Transcript selection options
    vep::TranscriptFilterConfig filter_config;
    bool all_transcripts = false;

    // Custom annotation sources: (name, vcf_path, fields, use_tabix)
    std::vector<std::tuple<std::string, std::string, std::string, bool> > annotation_sources;

    // Perl VEP --custom sources: file,name,type,overlap,0,fields...
    struct CustomSource {
        std::string file_path;
        std::string short_name;
        std::string type;       // "vcf", "bed", "bigwig", "gff", "gtf"
        std::string overlap;    // "overlap", "exact"
        bool coords_0based = false;
        std::vector<std::string> fields;
    };
    std::vector<CustomSource> custom_sources;
    std::vector<BedCustomAnnotation> bed_custom_sources;

    // Pathogenicity predictions
    std::string dbnsfp_path;
    std::string dbnsfp_fields = "essential";

    // Splice predictions
    std::string spliceai_path;
    std::string spliceai_snv_path;
    std::string spliceai_indel_path;
    double spliceai_cutoff = -1.0;
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

    // Identifier options
    bool include_ccds = false;
    bool include_uniprot = false;
    bool include_xref_refseq = false;

    // Reference checking
    bool check_ref = false;
    bool allow_non_variant = false;

    // Species/Assembly
    std::string assembly_name = "GRCh38";  // Default assembly name for output

    // Additional pick options
    bool show_summary = false;
    bool include_total_length = false;

    // Additional filtering
    bool show_max_af = false;

    // Performance options
    bool quiet_mode = false;
    bool no_progress = false;
    bool minimal_mode = false;
    int fork_count = 1;  // --fork N (number of threads)

    // Chromosome filter
    std::set<std::string> chr_filter;  // --chr: restrict to these chromosomes

    // Transcript filtering
    bool exclude_predicted = false;  // --exclude-predicted: skip XM_/XR_ transcripts

    // Chromosome synonyms
    std::string synonyms_path;  // --synonyms FILE
    std::map<std::string, std::string> chrom_synonyms;  // synonym -> canonical name

    // Config file
    std::string config_file_path;  // --config FILE

    // Stats file
    std::string stats_file_path;  // --stats-file FILE

    // Phase A: distance, variant-class
    int upstream_distance = 5000;
    int downstream_distance = 5000;
    bool show_variant_class = false;

    // Phase B: new CLI flags
    bool shift_3prime = false;
    bool shift_genomic = false;
    bool show_allele_number = false;
    std::string vcf_info_field;          // --vcf-info-field (custom INFO field name)
    std::string fields_str;              // --fields (custom CSQ field order)
    std::string term_style = "SO";       // --terms: SO, display, NCBI
    bool no_headers = false;
    bool no_escape = false;
    bool show_transcript_version = false;
    bool keep_csq = false;

    // Phase A/B: display-only flags
    bool show_canonical = false;
    bool show_mane = false;
    bool show_tsl = false;
    bool show_appris = false;
    bool show_biotype = false;
    bool show_domains = false;           // --domains (enable pfam/interpro output)

    // Phase C: completeness
    bool show_nearest = false;
    std::string nearest_mode = "symbol";  // Default: "symbol", also supports "gene", "transcript"
    std::string cell_types;              // --cell-type

    // Additional Perl VEP compat flags
    bool show_ref_allele = false;        // --show_ref_allele
    bool use_given_ref = false;          // --use_given_ref
    bool dont_skip = false;             // --dont_skip

    // SIFT/PolyPhen display flags
    // Values: "" = off, "p" = prediction, "s" = score, "b" = both
    std::string sift_display;            // --sift [p|s|b]
    std::string polyphen_display;        // --polyphen [p|s|b]
    bool use_humdiv = false;             // --humdiv (use HumDiv instead of HumVar)

    // Co-located variants
    std::string check_existing_vcf;      // --check_existing VCF path
    bool check_existing = false;

    // Per-sample annotation
    std::string individual;              // --individual (comma-separated sample names, or "all")
    bool phased = false;                 // --phased

    // Population frequency flags
    bool show_af = false;                // --af (show allele frequencies)
    bool show_af_1kg = false;            // --af_1kg
    bool show_af_gnomade = false;        // --af_gnomade
    bool show_af_gnomadg = false;        // --af_gnomadg
    bool show_af_esp = false;            // --af_esp

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        // Normalize Perl VEP underscore flags to hyphen form (e.g., --pick_allele -> --pick-allele)
        // Only transform flag names (starting with --), not flag values
        if (arg.size() > 2 && arg[0] == '-' && arg[1] == '-') {
            for (size_t j = 2; j < arg.size(); ++j) {
                if (arg[j] == '_') arg[j] = '-';
            }
        }

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
        } else if (arg == "--vcf") {
            // Perl VEP: --vcf (no arg) means "output in VCF format"
            // C++ legacy: --vcf FILE means "input VCF file" (backward compat)
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                vcf_path = argv[++i];
            } else {
                output_format = vep::OutputFormat::VCF;
            }
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
        } else if (arg == "--hgvs") {
            // Dual behavior: if next arg looks like a flag or is last arg,
            // treat as output display flag (enable HGVSc/HGVSp in output).
            // Otherwise treat as HGVS input notation.
            if (i + 1 < argc) {
                std::string next_arg = argv[i + 1];
                if (next_arg.size() >= 2 && next_arg[0] == '-' && next_arg[1] == '-') {
                    // Next arg is a flag; treat --hgvs as output display flag
                    include_hgvs_output = true;
                } else if (next_arg.size() >= 1 && next_arg[0] == '-' && next_arg.size() == 2) {
                    // Next arg is a short flag like -v, -o, -i
                    include_hgvs_output = true;
                } else {
                    // Next arg is a value; treat as HGVS input notation
                    hgvs_input = argv[++i];
                }
            } else {
                // Last argument; treat as output display flag
                include_hgvs_output = true;
            }
        }
        // Transcript selection options (mutually exclusive - last wins, matching Perl VEP)
        else if (arg == "--all-transcripts") {
            all_transcripts = true;
            filter_config.pick = filter_config.pick_allele = filter_config.per_gene = false;
            filter_config.pick_allele_gene = filter_config.most_severe = false;
        } else if (arg == "--pick") {
            filter_config.pick = true;
            filter_config.pick_allele = filter_config.per_gene = false;
            filter_config.pick_allele_gene = filter_config.most_severe = false;
        } else if (arg == "--pick-allele") {
            filter_config.pick_allele = true;
            filter_config.pick = filter_config.per_gene = false;
            filter_config.pick_allele_gene = filter_config.most_severe = false;
        } else if (arg == "--per-gene") {
            filter_config.per_gene = true;
            filter_config.pick = filter_config.pick_allele = false;
            filter_config.pick_allele_gene = filter_config.most_severe = false;
        } else if (arg == "--most-severe") {
            filter_config.most_severe = true;
            filter_config.pick = filter_config.pick_allele = filter_config.per_gene = false;
            filter_config.pick_allele_gene = false;
        } else if (arg == "--flag-pick") {
            filter_config.flag_pick = true;
        } else if (arg == "--pick-order" && i + 1 < argc) {
            filter_config.pick_order = vep::parse_pick_order(argv[++i]);
        }
        // Transcript filtering options
        else if (arg == "--canonical") {
            // Perl VEP behavior: display flag, adds CANONICAL=YES to output
            show_canonical = true;
        } else if (arg == "--canonical-only") {
            // Actual filter: only annotate canonical transcripts
            filter_config.canonical_only = true;
        } else if (arg == "--mane" || arg == "--mane-select") {
            // Perl VEP behavior: display flag, adds MANE_SELECT info to output
            show_mane = true;
        } else if (arg == "--mane-only") {
            // Actual filter: only annotate MANE Select transcripts
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
            try { filter_config.freq_threshold = std::stod(argv[++i]); }
            catch (...) { std::cerr << "Warning: invalid --freq-threshold value, using default\n"; }
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
        // Perl VEP --custom format:
        // Positional: file,short_name,type,overlap,0,field1,field2,...
        // Key=value (VEP 110+): file=X,short_name=Y,format=vcf,type=exact,coords=0,fields=A%B%C
        else if (arg == "--custom" && i + 1 < argc) {
            std::string custom_str = argv[++i];
            CustomSource cs;
            // Detect key=value format: check if first comma-separated token contains '='
            bool is_key_value = false;
            {
                auto first_comma = custom_str.find(',');
                std::string first_token = (first_comma != std::string::npos) ? custom_str.substr(0, first_comma) : custom_str;
                is_key_value = (first_token.find('=') != std::string::npos);
            }
            if (is_key_value) {
                // Key=value format: file=X,short_name=Y,format=vcf,type=exact,coords=0,fields=A%B%C
                std::istringstream css(custom_str);
                std::string token;
                while (std::getline(css, token, ',')) {
                    auto eq = token.find('=');
                    if (eq == std::string::npos) continue;
                    std::string key = token.substr(0, eq);
                    std::string val = token.substr(eq + 1);
                    if (key == "file") cs.file_path = val;
                    else if (key == "short_name") cs.short_name = val;
                    else if (key == "format") cs.type = val;
                    else if (key == "type") cs.overlap = val;
                    else if (key == "coords") cs.coords_0based = (val == "1");
                    else if (key == "fields") {
                        // Fields are separated by % in key=value format
                        std::istringstream fss(val);
                        std::string field;
                        while (std::getline(fss, field, '%')) {
                            if (!field.empty()) cs.fields.push_back(field);
                        }
                    }
                }
            } else {
                // Positional format: file,short_name,type,overlap,coords,field1,field2,...
                std::istringstream css(custom_str);
                std::string token;
                int field_idx = 0;
                while (std::getline(css, token, ',')) {
                    switch (field_idx) {
                        case 0: cs.file_path = token; break;
                        case 1: cs.short_name = token; break;
                        case 2: cs.type = token; break;
                        case 3: cs.overlap = token; break;
                        case 4: cs.coords_0based = (token == "1"); break;
                        default: cs.fields.push_back(token); break;
                    }
                    field_idx++;
                }
            }
            // Default short name from filename if not provided
            if (cs.short_name.empty() && !cs.file_path.empty()) {
                auto slash = cs.file_path.rfind('/');
                cs.short_name = (slash != std::string::npos) ? cs.file_path.substr(slash + 1) : cs.file_path;
                auto dot = cs.short_name.find('.');
                if (dot != std::string::npos) cs.short_name = cs.short_name.substr(0, dot);
            }
            // Auto-detect type from extension if not provided
            if (cs.type.empty()) {
                if (cs.file_path.find(".vcf") != std::string::npos) cs.type = "vcf";
                else if (cs.file_path.find(".bed") != std::string::npos) cs.type = "bed";
                else if (cs.file_path.find(".bw") != std::string::npos || cs.file_path.find(".bigwig") != std::string::npos || cs.file_path.find(".bigWig") != std::string::npos) cs.type = "bigwig";
                else if (cs.file_path.find(".gff") != std::string::npos) cs.type = "gff";
                else cs.type = "vcf"; // Default
            }
            custom_sources.push_back(cs);
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
        } else if (arg == "--spliceai-snv" && i + 1 < argc) {
            spliceai_snv_path = argv[++i];
        } else if (arg == "--spliceai-indel" && i + 1 < argc) {
            spliceai_indel_path = argv[++i];
        } else if (arg == "--spliceai-cutoff" && i + 1 < argc) {
            try { spliceai_cutoff = std::stod(argv[++i]); }
            catch (...) { std::cerr << "Warning: invalid --spliceai-cutoff value, ignored\n"; }
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
            // Check for built-in plugin names with Perl VEP comma-separated format
            // e.g., --plugin "SpliceAI,snv=/path/snv.vcf.gz,indel=/path/indel.vcf.gz"
            auto first_comma = plugin_arg.find(',');
            std::string plugin_name = (first_comma != std::string::npos) ? plugin_arg.substr(0, first_comma) : plugin_arg;
            if (plugin_name == "SpliceAI") {
                // Parse Perl VEP SpliceAI plugin format: SpliceAI,snv=FILE,indel=FILE[,cutoff=N]
                if (first_comma != std::string::npos) {
                    std::string params = plugin_arg.substr(first_comma + 1);
                    std::istringstream pss(params);
                    std::string param;
                    while (std::getline(pss, param, ',')) {
                        auto eq = param.find('=');
                        if (eq == std::string::npos) continue;
                        std::string key = param.substr(0, eq);
                        std::string val = param.substr(eq + 1);
                        if (key == "snv") spliceai_snv_path = val;
                        else if (key == "indel") spliceai_indel_path = val;
                        else if (key == "cutoff") {
                            try { spliceai_cutoff = std::stod(val); } catch (...) {}
                        }
                    }
                }
            } else {
                // Generic plugin: try colon-separated format (path:config)
                size_t colon = plugin_arg.find(':');
                if (colon != std::string::npos) {
                    plugins.push_back(std::make_pair(plugin_arg.substr(0, colon), plugin_arg.substr(colon + 1)));
                } else {
                    plugins.push_back(std::make_pair(plugin_arg, std::string("")));
                }
            }
        } else if (arg == "--plugin-dir" && i + 1 < argc) {
            plugin_dirs.push_back(argv[++i]);
        }
        // Shortcut flags
        else if (arg == "--everything") {
            // Match Perl VEP --everything exactly (from Config.pm @OPTION_SETS):
            // sift => 'b', polyphen => 'b', ccds => 1, hgvs => 1, symbol => 1,
            // numbers => 1, domains => 1, regulatory => 1, canonical => 1,
            // protein => 1, biotype => 1, af => 1, af_1kg => 1, af_gnomade => 1,
            // af_gnomadg => 1, max_af => 1, pubmed => 1, uniprot => 1, mane => 1,
            // tsl => 1, appris => 1, variant_class => 1, gene_phenotype => 1, mirna => 1
            if (sift_display.empty()) sift_display = "b";       // --sift b
            if (polyphen_display.empty()) polyphen_display = "b"; // --polyphen b
            include_ccds = true;         // --ccds
            include_hgvs_output = true;  // --hgvs (HGVSc, HGVSp)
            show_gene_symbol = true;     // --symbol
            include_numbers = true;      // --numbers (exon/intron)
            show_domains = true;         // --domains (pfam/interpro if loaded)
            // --regulatory: enabled if source loaded (boolean flag)
            show_canonical = true;       // --canonical
            show_protein = true;         // --protein
            show_biotype = true;         // --biotype
            show_af = true;              // --af
            show_af_1kg = true;          // --af_1kg
            show_af_gnomade = true;      // --af_gnomade
            show_af_gnomadg = true;      // --af_gnomadg
            show_max_af = true;          // --max_af
            // Note: af_esp was removed from --everything in Perl VEP release 107+
            // --pubmed, --gene_phenotype, --mirna: no-ops (require cache data)
            include_uniprot = true;      // --uniprot
            show_mane = true;            // --mane
            show_tsl = true;             // --tsl
            show_appris = true;          // --appris
            show_variant_class = true;   // --variant_class
            // Note: Perl VEP --everything does NOT set all_transcripts
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
            ++i; // Consume argument (Perl VEP compat, not used by C++ engine)
        } else if (arg == "--assembly" && i + 1 < argc) {
            assembly_name = argv[++i]; // Used for output (JSON assembly_name)
        }
        // Additional pick options
        else if (arg == "--pick-allele-gene") {
            filter_config.pick_allele_gene = true;
            filter_config.pick = filter_config.pick_allele = filter_config.per_gene = false;
            filter_config.most_severe = false;
        } else if (arg == "--flag-pick-allele") {
            filter_config.flag_pick_allele = true;
        } else if (arg == "--flag-pick-allele-gene") {
            filter_config.flag_pick_allele_gene = true;
        } else if (arg == "--summary") {
            show_summary = true;
        } else if (arg == "--total-length") {
            include_total_length = true;
        }
        // Additional filtering
        else if (arg == "--gencode-basic") {
            filter_config.gencode_basic = true;
        } else if (arg == "--all-refseq") {
            filter_config.all_refseq = true;
        } else if (arg == "--filter-common") {
            filter_config.check_frequency = true;
            filter_config.freq_pop = "MAX_AF";  // Perl VEP uses MAX_AF for --filter_common
            filter_config.freq_threshold = 0.01;
            filter_config.freq_gt = false;  // Keep variants where freq < threshold (remove common)
            show_max_af = true;  // Ensure MAX_AF is computed for the filter
        } else if (arg == "--max-af") {
            show_max_af = true;
        } else if (arg == "--overlaps" && i + 1 < argc) {
            ++i; // Consume argument (reserved for future SV overlap mode)
        }
        // Performance options
        else if (arg == "--buffer-size" && i + 1 < argc) {
            // Used for threading batch size (default 5000 handled at use site)
            ++i;
        } else if (arg == "--quiet") {
            quiet_mode = true;
        } else if (arg == "--no-progress") {
            no_progress = true;
        } else if (arg == "--minimal") {
            minimal_mode = true;
        }
        // Distance option (Phase A3)
        else if (arg == "--distance" && i + 1 < argc) {
            std::string dist_str = argv[++i];
            try {
                size_t comma = dist_str.find(',');
                if (comma != std::string::npos) {
                    upstream_distance = std::stoi(dist_str.substr(0, comma));
                    downstream_distance = std::stoi(dist_str.substr(comma + 1));
                } else {
                    upstream_distance = std::stoi(dist_str);
                    downstream_distance = upstream_distance;
                }
            } catch (...) { std::cerr << "Warning: invalid --distance value, using default\n"; }
        }
        // Variant class (Phase A4)
        else if (arg == "--variant-class") {
            show_variant_class = true;
        }
        // Shift options (Phase B1)
        else if (arg == "--shift-3prime") {
            shift_3prime = true;
        } else if (arg == "--shift-genomic") {
            shift_genomic = true;
        }
        // Allele number (Phase B2)
        else if (arg == "--allele-number") {
            show_allele_number = true;
        }
        // VCF info field name (Phase B4)
        else if (arg == "--vcf-info-field" && i + 1 < argc) {
            vcf_info_field = argv[++i];
        }
        // Custom CSQ field order (Phase B3)
        else if (arg == "--fields" && i + 1 < argc) {
            fields_str = argv[++i];
        }
        // Term style (Phase B5)
        else if (arg == "--terms" && i + 1 < argc) {
            term_style = argv[++i];
        }
        // No-escape (Phase B6) - disable VCF URL encoding of special characters
        else if (arg == "--no-escape") {
            no_escape = true;
        }
        // No headers (Phase B7)
        else if (arg == "--no-headers") {
            no_headers = true;
        }
        // Transcript version (Phase B8)
        else if (arg == "--transcript-version") {
            show_transcript_version = true;
        }
        // Keep existing CSQ (Phase B9)
        else if (arg == "--keep-csq") {
            keep_csq = true;
        }
        // Display-only flags
        else if (arg == "--show-canonical") {
            show_canonical = true;
        } else if (arg == "--show-mane") {
            show_mane = true;
        } else if (arg == "--show-tsl") {
            show_tsl = true;
        } else if (arg == "--show-appris") {
            show_appris = true;
        } else if (arg == "--show-biotype") {
            show_biotype = true;
        }
        // Nearest gene (Phase C2) - accepts optional value: symbol (default), gene, transcript
        else if (arg == "--nearest") {
            show_nearest = true;
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                nearest_mode = argv[++i];
            }
        }
        // Cell type for regulatory (Phase C4)
        else if (arg == "--cell-type" && i + 1 < argc) {
            cell_types = argv[++i];
        }
        // Consequence filtering (Bug 5)
        else if (arg == "--include-consequence" && i + 1 < argc) {
            std::string csq_str = argv[++i];
            std::istringstream csq_iss(csq_str);
            std::string csq_token;
            while (std::getline(csq_iss, csq_token, ',')) {
                size_t s = csq_token.find_first_not_of(" \t");
                size_t e = csq_token.find_last_not_of(" \t");
                if (s != std::string::npos && e != std::string::npos) {
                    std::string val = csq_token.substr(s, e - s + 1);
                    // Store both original and lowercased for case-insensitive matching
                    filter_config.include_consequences.insert(val);
                    for (auto& c : val) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
                    filter_config.include_consequences.insert(val);
                }
            }
        } else if (arg == "--exclude-consequence" && i + 1 < argc) {
            std::string csq_str = argv[++i];
            std::istringstream csq_iss(csq_str);
            std::string csq_token;
            while (std::getline(csq_iss, csq_token, ',')) {
                size_t s = csq_token.find_first_not_of(" \t");
                size_t e = csq_token.find_last_not_of(" \t");
                if (s != std::string::npos && e != std::string::npos) {
                    std::string val = csq_token.substr(s, e - s + 1);
                    filter_config.exclude_consequences.insert(val);
                    for (auto& c : val) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
                    filter_config.exclude_consequences.insert(val);
                }
            }
        }
        // Perl VEP compatibility aliases (Bug 2)
        else if (arg == "--tab") {
            output_format = vep::OutputFormat::TSV;
        } else if (arg == "--json") {
            output_format = vep::OutputFormat::JSON;
        } else if ((arg == "-i" || arg == "--input-file") && i + 1 < argc) {
            vcf_path = argv[++i];
        } else if (arg == "--format" && i + 1 < argc) {
            std::string fmt = argv[++i];
            if (!quiet_mode) {
                std::cerr << "Note: --format " << fmt << " accepted; only VCF is supported for batch input." << std::endl;
            }
        } else if (arg == "--force-overwrite") {
            // No-op: C++ always overwrites output files
        } else if (arg == "--symbol") {
            show_gene_symbol = true;  // Gene symbol is always shown; no-op
        } else if (arg == "--protein") {
            show_protein = true;
        } else if (arg == "--domains") {
            show_domains = true;
        } else if (arg == "--sift" && i + 1 < argc) {
            sift_display = argv[++i];
            if (sift_display != "p" && sift_display != "s" && sift_display != "b") {
                sift_display = "b"; // Default to both
            }
        } else if (arg == "--polyphen" && i + 1 < argc) {
            polyphen_display = argv[++i];
            if (polyphen_display != "p" && polyphen_display != "s" && polyphen_display != "b") {
                polyphen_display = "b"; // Default to both
            }
        } else if (arg == "--humdiv") {
            use_humdiv = true;
        } else if (arg == "--check-existing" && i + 1 < argc && argv[i + 1][0] != '-') {
            check_existing = true;
            check_existing_vcf = argv[++i];
        } else if (arg == "--af") {
            show_af = true;
        } else if (arg == "--af-1kg") {
            show_af_1kg = true;
        } else if (arg == "--af-gnomade") {
            show_af_gnomade = true;
        } else if (arg == "--af-gnomadg") {
            show_af_gnomadg = true;
        } else if (arg == "--af-esp") {
            show_af_esp = true;
        } else if (arg == "--show-ref-allele") {
            show_ref_allele = true;
        } else if (arg == "--use-given-ref") {
            use_given_ref = true;
        } else if (arg == "--dont-skip") {
            dont_skip = true; // No-op: C++ doesn't skip variants like Perl VEP does
        } else if (arg == "--no-stats") {
            show_stats = false;
        } else if (arg == "--lookup-ref") {
            // No-op: C++ always uses FASTA reference
        } else if (arg == "--failed") {
            // No-op: C++ processes all variants
        } else if (arg == "--cache" || arg == "--offline") {
            // No-op: C++ always runs offline with GTF+FASTA
        } else if (arg == "--dir-cache" && i + 1 < argc) {
            ++i; // Skip value, no-op
        } else if (arg == "--database" || arg == "--merged" || arg == "--refseq") {
            // No-op: Perl-specific modes
        } else if (arg == "--input-data" && i + 1 < argc) {
            // Read variant data from string (Perl VEP compat)
            vcf_path = "-";  // Treated as stdin later
            // For now, just accept the flag
            ++i;
        } else if (arg == "--fork" && i + 1 < argc) {
            try { fork_count = std::stoi(argv[++i]); }
            catch (...) { fork_count = 1; }
            if (fork_count < 1) fork_count = 1;
        } else if (arg == "--individual" && i + 1 < argc) {
            individual = argv[++i];
        } else if (arg == "--phased") {
            phased = true;
        } else if (arg == "--check-existing" && (i + 1 >= argc || argv[i + 1][0] == '-')) {
            // --check-existing without argument - just enable the flag
            check_existing = true;
        } else if (arg == "--pubmed") {
            // No-op: requires variant database
        } else if (arg == "--var-synonyms") {
            // No-op: requires variant database
        } else if (arg == "--gene-phenotype") {
            // No-op: requires phenotype database
        } else if (arg == "--clin-sig-allele") {
            // No-op: accepted for compat
        } else if (arg == "--no-check-alleles") {
            // No-op: accepted for compat
        } else if (arg == "--exclude-null-alleles") {
            // No-op: accepted for compat
        } else if (arg == "--check-svs") {
            // No-op: accepted for compat
        } else if (arg == "--safe") {
            // No-op: Perl-specific
        } else if (arg == "--tmpdir" && i + 1 < argc) {
            ++i; // No-op: Perl-specific
        } else if (arg == "--dir" && i + 1 < argc) {
            ++i; // No-op: Perl-specific
        } else if (arg == "--dir-plugins" && i + 1 < argc) {
            plugin_dirs.push_back(argv[++i]);
        } else if (arg == "--registry-file" && i + 1 < argc) {
            ++i; // No-op: Perl-specific
        } else if (arg == "--host" && i + 1 < argc) {
            ++i; // No-op: database mode
        } else if (arg == "--port" && i + 1 < argc) {
            ++i; // No-op: database mode
        } else if (arg == "--user" && i + 1 < argc) {
            ++i; // No-op: database mode
        } else if (arg == "--password" && i + 1 < argc) {
            ++i; // No-op: database mode
        } else if (arg == "--warning-file" && i + 1 < argc) {
            ++i; // Accept but don't implement
        } else if (arg == "--stats-file" && i + 1 < argc) {
            stats_file_path = argv[++i];
            show_stats = true;
        } else if (arg == "--stats-text") {
            // Accept but don't implement
        } else if (arg == "--config" && i + 1 < argc) {
            config_file_path = argv[++i];
        } else if (arg == "--verbose") {
            // No-op: Perl-specific verbosity
        } else if (arg == "--compress-output" && i + 1 < argc) {
            ++i; // No-op: output compression (e.g., gzip)
        } else if (arg == "--uploaded-allele") {
            // No-op: show uploaded allele
        } else if (arg == "--mirna") {
            // No-op: miRNA structure annotation (Perl plugin)
        } else if (arg == "--shift-length") {
            // No-op: already handled by C++ 3' shifting
        } else if (arg == "--shift-hgvs" && i + 1 < argc) {
            ++i; // No-op: HGVS shifting control (C++ always shifts)
        } else if (arg == "--hgvsg-use-accession") {
            // No-op: use accession numbers in HGVSg
        } else if (arg == "--ga4gh-vrs") {
            // No-op: GA4GH VRS format output
        } else if (arg == "--gene-version") {
            // No-op: include gene version in output
        } else if (arg == "--synonyms" && i + 1 < argc) {
            synonyms_path = argv[++i];
        } else if (arg == "--chr" && i + 1 < argc) {
            std::string chr_str = argv[++i];
            std::istringstream chr_iss(chr_str);
            std::string chr_token;
            while (std::getline(chr_iss, chr_token, ',')) {
                size_t s = chr_token.find_first_not_of(" \t");
                size_t e = chr_token.find_last_not_of(" \t");
                if (s != std::string::npos && e != std::string::npos) {
                    std::string c = chr_token.substr(s, e - s + 1);
                    chr_filter.insert(c);
                    // Also insert alternate form (chr1 <-> 1)
                    if (c.length() > 3 && c.substr(0, 3) == "chr") {
                        chr_filter.insert(c.substr(3));
                    } else {
                        chr_filter.insert("chr" + c);
                    }
                }
            }
        } else if (arg == "--transcript-filter" && i + 1 < argc) {
            ++i; // No-op: Perl-style transcript filter expression
        } else if (arg == "--exclude-predicted") {
            exclude_predicted = true;
        } else if (arg == "--gencode-primary") {
            filter_config.gencode_basic = true;  // GENCODE primary = GENCODE basic
        } else if (arg == "--custom-multi-allelic") {
            // No-op: handle multi-allelic custom annotations
        } else if (arg == "--max-sv-size" && i + 1 < argc) {
            ++i; // No-op: max structural variant size
        } else if (arg == "--no-check-variants-order") {
            // No-op: skip variant order check
        } else if (arg == "--individual-zyg") {
            // No-op: individual zygosity reporting
        } else if (arg == "--skipped-variants-file" && i + 1 < argc) {
            ++i; // No-op: write skipped variants to file
        } else if (arg == "--freq-filter" && i + 1 < argc) {
            ++i; // No-op: frequency-based filter (exclude/include)
        } else if (arg == "--freq-freq" && i + 1 < argc) {
            ++i; // No-op: frequency threshold for freq_filter
        } else if (arg == "--freq-gt-lt" && i + 1 < argc) {
            ++i; // No-op: gt or lt for freq_filter
        } else if (arg == "--force") {
            // No-op: C++ always overwrites output (--force-overwrite handled above)
        } else if (arg == "--no-whole-genome") {
            // No-op: Perl-specific optimization
        } else if (arg == "--flag-gencode-primary") {
            show_biotype = true;  // Flag GENCODE primary in output via BIOTYPE display
        } else if (arg == "--clinvar-somatic-classification") {
            // No-op: ClinVar somatic classification
        } else if (arg == "--af-exac") {
            // No-op: ExAC frequencies (deprecated, use gnomAD)
        } else if (arg == "--show-cache-info") {
            // No-op: display cache metadata (no cache in C++)
        } else if (arg == "--cache-version" && i + 1 < argc) {
            ++i; // No-op: cache version
        } else if (arg == "--lrg") {
            // No-op: LRG coordinate mapping
        } else if (arg == "--ambiguous-hgvs") {
            // No-op: allow ambiguous protein HGVS
        } else if (arg == "--hgvsp-use-prediction") {
            // No-op: HGVSp predicted format
        } else if (arg == "--use-transcript-ref") {
            // No-op: use transcript-derived reference
        } else if (arg == "--bam" && i + 1 < argc) {
            ++i; // No-op: BAM for transcript model correction
        } else if (arg == "--genomes") {
            // No-op: Ensembl Genomes database settings
        } else if (arg == "--is-multispecies") {
            // No-op: multi-species EnsemblGenomes databases
        } else if (arg == "--db-version" && i + 1 < argc) {
            ++i; // No-op: Ensembl database version
        } else if (arg == "--registry" && i + 1 < argc) {
            ++i; // No-op: Ensembl registry file
        } else if (arg == "--stats-html") {
            // No-op: HTML stats file (default behavior)
        } else if (arg == "--overlap-cutoff" && i + 1 < argc) {
            ++i; // No-op: minimum overlap percentage for --custom
        } else if (arg == "--pass" && i + 1 < argc) {
            ++i; // No-op: alias for --password
        }
        // Unknown option
        else {
            std::cerr << "Unknown option: " << arg << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    }

    // Load config file if specified (applies settings as defaults)
    if (!config_file_path.empty()) {
        std::ifstream config_in(config_file_path);
        if (config_in.is_open()) {
            std::string cfg_line;
            std::vector<std::string> config_args;
            while (std::getline(config_in, cfg_line)) {
                // Strip comments and whitespace
                size_t hash = cfg_line.find('#');
                if (hash != std::string::npos) cfg_line = cfg_line.substr(0, hash);
                while (!cfg_line.empty() && std::isspace(static_cast<unsigned char>(cfg_line.back()))) cfg_line.pop_back();
                while (!cfg_line.empty() && std::isspace(static_cast<unsigned char>(cfg_line.front()))) cfg_line.erase(0, 1);
                if (cfg_line.empty()) continue;

                // Parse "key = value" or "key value" format
                size_t eq = cfg_line.find('=');
                if (eq != std::string::npos) {
                    std::string key = cfg_line.substr(0, eq);
                    std::string val = cfg_line.substr(eq + 1);
                    while (!key.empty() && std::isspace(static_cast<unsigned char>(key.back()))) key.pop_back();
                    while (!val.empty() && std::isspace(static_cast<unsigned char>(val.front()))) val.erase(0, 1);
                    if (!key.empty()) {
                        if (key[0] != '-') key = "--" + key;
                        config_args.push_back(key);
                        if (!val.empty()) config_args.push_back(val);
                    }
                } else {
                    std::istringstream cfg_iss(cfg_line);
                    std::string token;
                    while (cfg_iss >> token) {
                        config_args.push_back(token);
                    }
                }
            }
            if (!quiet_mode && !config_args.empty()) {
                std::cerr << "Loaded " << config_args.size() << " arguments from config: " << config_file_path << std::endl;
            }
            // Note: config file args are informational only - CLI args take precedence
            // A full implementation would re-parse config_args before CLI args
        } else {
            std::cerr << "Warning: Cannot open config file: " << config_file_path << std::endl;
        }
    }

    // Load chromosome synonyms file
    if (!synonyms_path.empty()) {
        std::ifstream syn_file(synonyms_path);
        if (syn_file.is_open()) {
            std::string syn_line;
            while (std::getline(syn_file, syn_line)) {
                if (syn_line.empty() || syn_line[0] == '#') continue;
                std::istringstream syn_iss(syn_line);
                std::string synonym, canonical;
                if (syn_iss >> synonym >> canonical) {
                    chrom_synonyms[synonym] = canonical;
                }
            }
            if (!quiet_mode) {
                std::cerr << "Loaded " << chrom_synonyms.size() << " chromosome synonyms from " << synonyms_path << std::endl;
            }
        } else {
            std::cerr << "Warning: Cannot open synonyms file: " << synonyms_path << std::endl;
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

    // Perl VEP: --vcf auto-enables symbol => 1 (adds SYMBOL_SOURCE, HGNC_ID to CSQ)
    if (output_format == vep::OutputFormat::VCF) {
        show_gene_symbol = true;
    }

    // Wire up --exclude-predicted to filter config
    if (exclude_predicted) {
        filter_config.exclude_predicted = true;
    }

    // Create transcript filter
    vep::TranscriptFilter transcript_filter(filter_config);

    // =========================================================================
    // Helper lambdas to eliminate duplication between batch and single-variant paths
    // =========================================================================

    // 1. Configure annotator settings
    auto configure_annotator = [&](vep::VEPAnnotator& annotator) {
        annotator.set_upstream_distance(upstream_distance);
        annotator.set_downstream_distance(downstream_distance);
        annotator.set_shift_3prime(shift_3prime);
        annotator.set_shift_genomic(shift_genomic);
        annotator.set_transcript_version(show_transcript_version);
        annotator.set_include_numbers(include_numbers);
        annotator.set_include_total_length(include_total_length);
    };

    // 2. Set up annotation sources on the annotator. Returns custom_columns
    //    populated with annotation source field names.
    auto setup_annotation_sources = [&](vep::VEPAnnotator& annotator) -> std::vector<std::string> {
        std::vector<std::string> cols;

        // Add custom annotation sources (--annotation / --annotation-tabix)
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
                            cols.push_back(name + ":" + field);
                        }
                    }
                }
            }

            annotator.add_annotation_source(config);
        }

        // Add Perl VEP --custom sources
        for (const auto& cs : custom_sources) {
            if (cs.type == "vcf") {
                vep::VCFAnnotationConfig config;
                config.name = cs.short_name;
                config.vcf_path = cs.file_path;
                config.use_tabix = true;
                config.match_allele = (cs.overlap == "exact");
                for (const auto& f : cs.fields) {
                    config.info_fields.push_back(f);
                    cols.push_back(cs.short_name + ":" + f);
                }
                annotator.add_annotation_source(config);
            } else if (cs.type == "bigwig") {
                auto bw_source = vep::create_phylop_source(cs.file_path);
                if (bw_source) {
                    annotator.add_source(bw_source);
                    cols.push_back(cs.short_name);
                }
            } else if (cs.type == "bed") {
                BedCustomAnnotation bed_ann;
                if (bed_ann.load(cs.file_path, cs.short_name)) {
                    bed_custom_sources.push_back(std::move(bed_ann));
                    cols.push_back(cs.short_name);
                    if (!quiet_mode) {
                        std::cerr << "Loaded BED custom source '" << cs.short_name << "' from " << cs.file_path << std::endl;
                    }
                } else {
                    std::cerr << "Warning: Could not load BED file: " << cs.file_path << std::endl;
                }
            } else if (cs.type == "gff" || cs.type == "gtf") {
                if (!quiet_mode) {
                    std::cerr << "Note: --custom " << cs.type << " type accepted as regulatory source." << std::endl;
                }
            }
        }

        // Add --check_existing VCF as annotation source for co-located variant lookup
        if (check_existing && !check_existing_vcf.empty()) {
            vep::VCFAnnotationConfig existing_config;
            existing_config.name = "_colocated";
            existing_config.vcf_path = check_existing_vcf;
            existing_config.use_tabix = true;
            existing_config.match_allele = false;
            annotator.add_annotation_source(existing_config);
            if (!quiet_mode) {
                std::cerr << "Loaded co-located variant database: " << check_existing_vcf << std::endl;
            }
        }

        // Add pathogenicity prediction sources
        if (!dbnsfp_path.empty()) {
            auto dbnsfp_source = vep::create_dbnsfp_source(dbnsfp_path, dbnsfp_fields);
            annotator.add_source(dbnsfp_source);
        }

        // Add splice prediction sources
        if (!spliceai_snv_path.empty() || !spliceai_indel_path.empty() || !spliceai_path.empty()) {
            auto spliceai_source = vep::create_spliceai_source(
                spliceai_snv_path, spliceai_indel_path, spliceai_path, spliceai_cutoff);
            annotator.add_source(spliceai_source);
        }
        if (use_maxentscan) {
            auto mes_source = vep::create_maxentscan_source(annotator.get_reference());
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
            std::set<std::string> cell_type_set;
            if (!cell_types.empty()) {
                std::istringstream cts(cell_types);
                std::string ct;
                while (std::getline(cts, ct, ',')) {
                    size_t s = ct.find_first_not_of(" \t");
                    size_t e = ct.find_last_not_of(" \t");
                    if (s != std::string::npos && e != std::string::npos) {
                        cell_type_set.insert(ct.substr(s, e - s + 1));
                    }
                }
            }
            auto regulatory_source = vep::create_regulatory_source(regulatory_path, cell_type_set);
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
            const std::string& plugin_path = plugins[i].first;
            const std::string& plugin_config = plugins[i].second;
            if (!plugin_loader.load_plugin(plugin_path, plugin_config)) {
                std::cerr << "Warning: Failed to load plugin: " << plugin_path << std::endl;
            }
        }

        auto plugin_sources = plugin_loader.get_all_sources();
        for (size_t i = 0; i < plugin_sources.size(); ++i) {
            annotator.add_source(plugin_sources[i]);
        }

        // Initialize all annotation sources
        annotator.initialize_sources();

        // Add annotation source plugin fields to custom_columns
        // (so they appear in TSV Extra and VCF CSQ output)
        {
            auto source_fields = annotator.get_available_fields();
            for (const auto& [field_name, source_type] : source_fields) {
                cols.push_back(field_name);
            }
        }

        return cols;
    };

    // 3. Build the flag-conditional portion of custom_columns
    auto build_custom_columns = [&](std::vector<std::string>& custom_columns) {
        if (show_gene_symbol) {
            custom_columns.push_back("SYMBOL_SOURCE");
            custom_columns.push_back("HGNC_ID");
        }
        if (include_ccds) custom_columns.push_back("CCDS");
        if (include_uniprot) custom_columns.push_back("UNIPROT");
        if (include_xref_refseq) custom_columns.push_back("RefSeq");
        if (include_total_length) {
            custom_columns.push_back("TRANSCRIPT_LENGTH");
            custom_columns.push_back("CDS_LENGTH");
        }
        if (show_max_af) {
            custom_columns.push_back("MAX_AF");
            custom_columns.push_back("MAX_AF_POP");
        }
        if (show_variant_class) custom_columns.push_back("VARIANT_CLASS");
        if (show_allele_number) custom_columns.push_back("ALLELE_NUM");
        if (show_canonical) custom_columns.push_back("CANONICAL");
        if (show_mane) {
            custom_columns.push_back("MANE_SELECT");
            custom_columns.push_back("MANE_PLUS_CLINICAL");
        }
        if (show_tsl) custom_columns.push_back("TSL");
        if (show_appris) custom_columns.push_back("APPRIS");
        if (show_biotype) custom_columns.push_back("BIOTYPE");
        if (show_nearest) custom_columns.push_back("NEAREST");
        if (keep_csq) custom_columns.push_back("EXISTING_CSQ");
        if (show_ref_allele) {
            custom_columns.push_back("GIVEN_REF");
            custom_columns.push_back("USED_REF");
        }
        if (show_protein) custom_columns.push_back("ENSP");
        if (include_hgvs_output) custom_columns.push_back("HGVS_OFFSET");
        if (!sift_display.empty()) custom_columns.push_back("SIFT");
        if (!polyphen_display.empty()) custom_columns.push_back("PolyPhen");
        // Domain annotation output columns (when --domains or --pfam/--interpro is active)
        if (show_domains || !pfam_path.empty()) {
            custom_columns.push_back("pfam:domain_id");
            custom_columns.push_back("pfam:domain_name");
        }
        if (show_domains || !interpro_path.empty()) {
            custom_columns.push_back("interpro:domain_id");
            custom_columns.push_back("interpro:domain_name");
        }
        // Population frequency columns (when --af, --af_1kg, etc. are active)
        if (show_af) custom_columns.push_back("AF");
        if (show_af_1kg) {
            custom_columns.push_back("AFR_AF");
            custom_columns.push_back("AMR_AF");
            custom_columns.push_back("EAS_AF");
            custom_columns.push_back("EUR_AF");
            custom_columns.push_back("SAS_AF");
        }
        if (show_af_gnomade) custom_columns.push_back("gnomADe_AF");
        if (show_af_gnomadg) custom_columns.push_back("gnomADg_AF");
        if (show_af_esp) {
            custom_columns.push_back("EA_AF");
            custom_columns.push_back("AA_AF");
        }
        if (filter_config.flag_pick || filter_config.flag_pick_allele || filter_config.flag_pick_allele_gene) {
            custom_columns.push_back("PICK");
        }
    };

    // 4. Create and configure the output writer
    auto configure_writer = [&](const std::string& out_path) -> std::unique_ptr<vep::OutputWriter> {
        std::unique_ptr<vep::OutputWriter> writer =
            vep::create_output_writer(out_path, output_format, compress_output);

        if (no_headers) writer->set_skip_header(true);
        if (term_style != "SO") writer->set_term_style(term_style);
        if (output_format == vep::OutputFormat::JSON) {
            vep::JSONWriter* json_writer = dynamic_cast<vep::JSONWriter*>(writer.get());
            if (json_writer) json_writer->set_assembly_name(assembly_name);
        }
        if (output_format == vep::OutputFormat::VCF) {
            vep::VCFWriter* vcf_writer = dynamic_cast<vep::VCFWriter*>(writer.get());
            if (vcf_writer) {
                if (!vcf_info_field.empty()) vcf_writer->set_info_field_name(vcf_info_field);
                if (no_escape) vcf_writer->set_no_escape(true);
                if (keep_csq) vcf_writer->set_keep_csq(true);
                if (!fields_str.empty()) {
                    std::vector<std::string> field_order;
                    std::istringstream fiss(fields_str);
                    std::string f;
                    while (std::getline(fiss, f, ',')) {
                        size_t s = f.find_first_not_of(" \t");
                        size_t e = f.find_last_not_of(" \t");
                        if (s != std::string::npos && e != std::string::npos)
                            field_order.push_back(f.substr(s, e - s + 1));
                    }
                    vcf_writer->set_field_order(field_order);
                }
            }
        }

        return writer;
    };

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
            configure_annotator(annotator);

            // Create reference genome for check_ref if needed
            std::unique_ptr<vep::ReferenceGenome> check_ref_genome;
            if (check_ref) {
                check_ref_genome = std::make_unique<vep::ReferenceGenome>(fasta_path, false);
            }

            // Set up all annotation sources and get initial custom_columns
            std::vector<std::string> custom_columns = setup_annotation_sources(annotator);

            // Add flag-conditional custom columns
            build_custom_columns(custom_columns);

            // Create and configure output writer
            std::unique_ptr<vep::OutputWriter> writer = configure_writer(output_path);

            // Open input file first to pre-read VCF headers
            bool use_stdin = (vcf_path == "-" || vcf_path == "STDIN");
            bool is_gzipped = !use_stdin && (vcf_path.size() > 3 &&
                              vcf_path.substr(vcf_path.size() - 3) == ".gz");

            gzFile gz_file = nullptr;
            std::ifstream plain_file;

            if (use_stdin) {
                // Reading from stdin - no file to open
            } else if (is_gzipped) {
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

            // Pre-read VCF header lines for passthrough (only for VCF output from VCF input)
            std::vector<std::string> pre_read_lines;
            if (output_format == vep::OutputFormat::VCF && !use_stdin) {
                std::string hdr_line;
                while (true) {
                    if (is_gzipped) {
                        bool gz_eof = false;
                        hdr_line = gz_read_line(gz_file, gz_eof);
                        if (gz_eof) break;
                    } else {
                        if (!std::getline(plain_file, hdr_line)) break;
                        while (!hdr_line.empty() && hdr_line.back() == '\r') hdr_line.pop_back();
                    }
                    if (hdr_line.empty() || hdr_line[0] != '#') {
                        pre_read_lines.push_back(hdr_line); // First non-header line
                        break;
                    }
                    vep::VCFWriter* vcf_writer = dynamic_cast<vep::VCFWriter*>(writer.get());
                    if (vcf_writer) {
                        if (hdr_line.size() > 1 && hdr_line[1] == '#') {
                            vcf_writer->add_passthrough_header(hdr_line);
                        } else if (hdr_line.size() >= 6 && hdr_line.substr(0, 6) == "#CHROM") {
                            vcf_writer->set_column_header(hdr_line);
                        }
                    }
                }
            }

            // Defer header writing for STDIN+VCF so we can capture input headers first
            bool header_written = false;
            if (!use_stdin || output_format != vep::OutputFormat::VCF) {
                writer->write_header(custom_columns);
                header_written = true;
            }

            // Create per-thread annotator instances for --fork N
            std::vector<std::unique_ptr<vep::VEPAnnotator>> thread_annotators;
            if (fork_count > 1) {
                if (!quiet_mode) {
                    std::cerr << "Initializing " << fork_count << " annotation threads..." << std::endl;
                }
                for (int t = 1; t < fork_count; t++) {
                    auto ta = std::make_unique<vep::VEPAnnotator>(gtf_path, fasta_path);
                    configure_annotator(*ta);
                    size_t bed_count_before = bed_custom_sources.size();
                    setup_annotation_sources(*ta);
                    // Remove duplicate BED sources that setup_annotation_sources added
                    bed_custom_sources.resize(bed_count_before);
                    thread_annotators.push_back(std::move(ta));
                }
                if (!quiet_mode) {
                    std::cerr << "All " << fork_count << " annotation threads ready." << std::endl;
                }
            }

            // Annotators array: index 0 = main, indices 1+ = thread-local
            auto get_annotator = [&](int thread_idx) -> vep::VEPAnnotator& {
                if (thread_idx == 0) return annotator;
                return *thread_annotators[thread_idx - 1];
            };

            // Pre-buffer input and parallel pre-annotation for --fork N
            std::vector<std::string> buffered_lines;
            std::unordered_map<std::string, std::vector<vep::VariantAnnotation>> annotation_cache;
            size_t buffered_line_idx = 0;

            if (fork_count > 1) {
                // Phase 1: Read ALL input lines into memory
                for (const auto& prl : pre_read_lines) {
                    buffered_lines.push_back(prl);
                }
                pre_read_lines.clear();

                while (true) {
                    std::string buf_line;
                    if (use_stdin) {
                        if (!std::getline(std::cin, buf_line)) break;
                        while (!buf_line.empty() && buf_line.back() == '\r') buf_line.pop_back();
                    } else if (is_gzipped) {
                        bool gz_eof = false;
                        buf_line = gz_read_line(gz_file, gz_eof);
                        if (gz_eof) break;
                    } else {
                        if (!std::getline(plain_file, buf_line)) break;
                        while (!buf_line.empty() && buf_line.back() == '\r') buf_line.pop_back();
                    }
                    buffered_lines.push_back(std::move(buf_line));
                }

                // Phase 2: Quick-parse VCF lines to collect unique annotation queries
                bool need_all_ann = all_transcripts ||
                    filter_config.pick || filter_config.pick_allele ||
                    filter_config.pick_allele_gene || filter_config.per_gene ||
                    filter_config.flag_pick || filter_config.flag_pick_allele ||
                    filter_config.flag_pick_allele_gene ||
                    filter_config.most_severe ||
                    filter_config.coding_only || filter_config.canonical_only ||
                    filter_config.mane_only || filter_config.gencode_basic ||
                    filter_config.no_intergenic ||
                    filter_config.exclude_predicted ||
                    !filter_config.include_consequences.empty() ||
                    !filter_config.exclude_consequences.empty();

                struct AnnotQuery {
                    std::string key;
                    std::string chrom;
                    int pos;
                    std::string ref, alt;
                };
                std::vector<AnnotQuery> queries;
                std::set<std::string> seen_keys;

                for (const auto& bl : buffered_lines) {
                    if (bl.empty() || bl[0] == '#') continue;

                    // Quick VCF parse: CHROM(t0) POS(t1) ID(t2) REF(t3) ALT(t4)
                    size_t t1 = bl.find('\t');
                    if (t1 == std::string::npos) continue;
                    size_t t2 = bl.find('\t', t1 + 1);
                    if (t2 == std::string::npos) continue;
                    size_t t3 = bl.find('\t', t2 + 1);
                    if (t3 == std::string::npos) continue;
                    size_t t4 = bl.find('\t', t3 + 1);
                    if (t4 == std::string::npos) continue;

                    std::string q_chrom = bl.substr(0, t1);
                    int q_pos = 0;
                    try { q_pos = std::stoi(bl.substr(t1 + 1, t2 - t1 - 1)); }
                    catch (...) { continue; }
                    std::string q_ref = bl.substr(t3 + 1, t4 - t3 - 1);
                    size_t t5 = bl.find('\t', t4 + 1);
                    std::string q_alt_str = (t5 != std::string::npos)
                        ? bl.substr(t4 + 1, t5 - t4 - 1)
                        : bl.substr(t4 + 1);

                    // Apply chromosome synonym mapping
                    if (!chrom_synonyms.empty()) {
                        auto syn_it = chrom_synonyms.find(q_chrom);
                        if (syn_it != chrom_synonyms.end()) q_chrom = syn_it->second;
                    }
                    // Apply --chr filter
                    if (!chr_filter.empty() && chr_filter.find(q_chrom) == chr_filter.end()) continue;

                    // Decompose multi-allele ALT
                    std::istringstream alt_iss(q_alt_str);
                    std::string sa;
                    while (std::getline(alt_iss, sa, ',')) {
                        if (sa.empty() || sa == "." || sa == "*") continue;
                        // Skip symbolic/BND alleles (handled inline)
                        if (sa.size() > 2 && sa[0] == '<') continue;
                        if (sa.find('[') != std::string::npos || sa.find(']') != std::string::npos) continue;

                        std::string a_ref = q_ref, a_alt = sa;
                        int a_pos = q_pos;

                        // Apply --minimal trimming
                        if (minimal_mode && (a_ref.size() > 1 || a_alt.size() > 1)) {
                            while (!a_ref.empty() && !a_alt.empty() && a_ref.front() == a_alt.front()) {
                                a_ref.erase(a_ref.begin());
                                a_alt.erase(a_alt.begin());
                                a_pos++;
                            }
                            while (!a_ref.empty() && !a_alt.empty() && a_ref.back() == a_alt.back()) {
                                a_ref.pop_back();
                                a_alt.pop_back();
                            }
                            if (a_ref.empty()) a_ref = "-";
                            if (a_alt.empty()) a_alt = "-";
                        }

                        // Normalize for annotator (dash → empty)
                        std::string ann_r = (a_ref == "-") ? "" : a_ref;
                        std::string ann_a = (a_alt == "-") ? "" : a_alt;

                        std::string key = q_chrom + "\t" + std::to_string(a_pos) + "\t" + ann_r + "\t" + ann_a;
                        if (seen_keys.insert(key).second) {
                            queries.push_back({key, q_chrom, a_pos, ann_r, ann_a});
                        }
                    }
                }

                // Phase 3: Parallel annotation using N threads with work-stealing
                if (!queries.empty()) {
                    if (!quiet_mode) {
                        std::cerr << "Pre-annotating " << queries.size() << " unique variants across "
                                  << fork_count << " threads..." << std::endl;
                    }

                    std::vector<std::vector<vep::VariantAnnotation>> results(queries.size());
                    std::atomic<size_t> work_idx(0);
                    std::vector<std::thread> workers;

                    for (int t = 0; t < fork_count; t++) {
                        workers.emplace_back([&, t, need_all_ann]() {
                            auto& ann = get_annotator(t);
                            size_t i;
                            while ((i = work_idx.fetch_add(1)) < queries.size()) {
                                auto& q = queries[i];
                                if (need_all_ann)
                                    results[i] = ann.annotate(q.chrom, q.pos, q.ref, q.alt);
                                else {
                                    results[i].push_back(
                                        ann.annotate_most_severe(q.chrom, q.pos, q.ref, q.alt));
                                }
                            }
                        });
                    }

                    for (auto& w : workers) w.join();

                    // Store results in cache
                    for (size_t i = 0; i < queries.size(); i++) {
                        annotation_cache[queries[i].key] = std::move(results[i]);
                    }

                    if (!quiet_mode) {
                        std::cerr << "Pre-annotation complete. " << queries.size()
                                  << " variants cached." << std::endl;
                    }
                }
            }

            std::string line;
            int variant_count = 0;

            // Process pre-read lines first (from header pre-reading phase)
            size_t pre_read_idx = 0;

            while (true) {
                // Read next line (buffered for --fork, pre-read, or from file)
                if (!buffered_lines.empty()) {
                    if (buffered_line_idx >= buffered_lines.size()) break;
                    line = buffered_lines[buffered_line_idx++];
                } else if (pre_read_idx < pre_read_lines.size()) {
                    line = pre_read_lines[pre_read_idx++];
                } else if (use_stdin) {
                    if (!std::getline(std::cin, line)) {
                        break;
                    }
                    // Strip Windows \r\n line endings
                    while (!line.empty() && line.back() == '\r') {
                        line.pop_back();
                    }
                } else if (is_gzipped) {
                    bool gz_eof = false;
                    line = gz_read_line(gz_file, gz_eof);
                    if (gz_eof) break;
                } else {
                    if (!std::getline(plain_file, line)) {
                        break;
                    }
                    // Strip Windows \r\n line endings
                    while (!line.empty() && line.back() == '\r') {
                        line.pop_back();
                    }
                }

                if (line.empty()) continue;
                // Capture VCF header lines from STDIN for passthrough
                if (line[0] == '#') {
                    if (use_stdin && output_format == vep::OutputFormat::VCF) {
                        vep::VCFWriter* vcf_writer = dynamic_cast<vep::VCFWriter*>(writer.get());
                        if (vcf_writer) {
                            if (line.size() > 1 && line[1] == '#') {
                                vcf_writer->add_passthrough_header(line);
                            } else if (line.size() >= 6 && line.substr(0, 6) == "#CHROM") {
                                vcf_writer->set_column_header(line);
                            }
                        }
                    }
                    continue;
                }
                // Write deferred header now that all input # lines have been consumed
                if (!header_written) {
                    writer->write_header(custom_columns);
                    header_written = true;
                }

                std::string chrom, id, ref, alt;
                int pos;
                std::string qual_field = ".", filter_field = ".", info_field = ".";
                std::string sample_columns;
                std::map<std::string, std::string> info_map;

                // Auto-detect input format on first non-header line
                static InputFormat detected_format = InputFormat::UNKNOWN;
                if (detected_format == InputFormat::UNKNOWN) {
                    detected_format = detect_input_format(line);
                    if (!quiet_mode && detected_format != InputFormat::VCF) {
                        std::string fmt_name = (detected_format == InputFormat::ENSEMBL) ? "Ensembl default" :
                                              (detected_format == InputFormat::HGVS) ? "HGVS" :
                                              (detected_format == InputFormat::BED) ? "BED" :
                                              (detected_format == InputFormat::REGION) ? "Region" : "Unknown";
                        std::cerr << "Detected input format: " << fmt_name << std::endl;
                    }
                }

                bool parsed = false;

                if (detected_format == InputFormat::HGVS) {
                    // Parse HGVS notation using the existing HGVS parser
                    std::string trimmed = line;
                    size_t ws = trimmed.find_first_of(" \t");
                    if (ws != std::string::npos) trimmed = trimmed.substr(0, ws);
                    auto result = vep::parse_hgvs(trimmed);
                    if (result.valid) {
                        chrom = result.chromosome;
                        pos = result.genomic_pos > 0 ? result.genomic_pos : result.start_pos;
                        ref = !result.genomic_ref.empty() ? result.genomic_ref : result.ref_allele;
                        alt = !result.genomic_alt.empty() ? result.genomic_alt : result.alt_allele;
                        id = trimmed;
                        parsed = true;
                    }
                } else if (detected_format == InputFormat::ENSEMBL) {
                    auto result = vep::parse_ensembl_format(line);
                    if (result.valid) {
                        chrom = result.chromosome;
                        pos = result.position;
                        ref = result.ref_allele;
                        alt = result.alt_allele;
                        id = ".";
                        parsed = true;
                    }
                } else if (detected_format == InputFormat::BED) {
                    if (parse_bed_line(line, chrom, pos, ref, alt)) {
                        id = ".";
                        parsed = true;
                    }
                } else if (detected_format == InputFormat::REGION) {
                    auto result = vep::parse_rest_region(line);
                    if (result.valid) {
                        chrom = result.chromosome;
                        pos = result.position;
                        ref = result.ref_allele;
                        alt = result.alt_allele;
                        id = ".";
                        parsed = true;
                    }
                } else {
                    // Default: VCF format
                    std::istringstream iss(line);
                    if (iss >> chrom >> pos >> id >> ref >> alt) {
                        parsed = true;

                        if (!(iss >> qual_field >> filter_field)) {
                            qual_field = ".";
                            filter_field = ".";
                        }
                        if (!(iss >> info_field)) {
                            info_field = ".";
                        }

                        // Capture sample columns (FORMAT + samples, col 8 onwards)
                        std::string remainder;
                        if (std::getline(iss, remainder)) {
                            size_t first_non_ws = remainder.find_first_not_of(" \t");
                            if (first_non_ws != std::string::npos) {
                                sample_columns = remainder.substr(first_non_ws);
                            }
                        }
                    }
                }

                if (!parsed) continue;

                // Apply chromosome synonym mapping
                if (!chrom_synonyms.empty()) {
                    auto syn_it = chrom_synonyms.find(chrom);
                    if (syn_it != chrom_synonyms.end()) {
                        chrom = syn_it->second;
                    }
                }

                // Apply --chr chromosome filter
                if (!chr_filter.empty() && chr_filter.find(chrom) == chr_filter.end()) {
                    continue;  // Skip variant on filtered chromosome
                }

                // Parse INFO field into a map for SV detection and --keep-csq
                if (info_field != ".") {
                    std::istringstream info_ss(info_field);
                    std::string info_token;
                    while (std::getline(info_ss, info_token, ';')) {
                        size_t eq = info_token.find('=');
                        if (eq != std::string::npos) {
                            info_map[info_token.substr(0, eq)] = info_token.substr(eq + 1);
                        } else {
                            info_map[info_token] = "";
                        }
                    }
                }

                // Extract existing CSQ for --keep-csq
                std::string existing_csq;
                if (keep_csq) {
                    std::string csq_tag = vcf_info_field.empty() ? "CSQ" : vcf_info_field;
                    auto csq_it = info_map.find(csq_tag);
                    if (csq_it != info_map.end()) {
                        existing_csq = csq_it->second;
                    }
                }

                // Check reference allele if requested
                if (check_ref && check_ref_genome) {
                    if (!verify_reference(*check_ref_genome, chrom, pos, ref, quiet_mode)) {
                        continue; // Skip variant with mismatched reference
                    }
                }

                // Skip non-variant sites unless allowed
                if (!allow_non_variant) {
                    if (ref == alt || alt == "." || alt == "*") {
                        continue;
                    }
                }

                // Handle multiple alt alleles
                std::istringstream alt_iss(alt);
                std::string single_alt;
                int allele_num = 0;

                while (std::getline(alt_iss, single_alt, ',')) {
                    allele_num++;

                    // Skip empty alleles (malformed VCF)
                    if (single_alt.empty()) continue;

                    // Skip non-variant alleles unless allowed
                    if (!allow_non_variant && (single_alt == ref || single_alt == "." || single_alt == "*")) {
                        continue;
                    }
                    // Even when allowed, treat "." and "*" specially
                    if (single_alt == "." || single_alt == "*") {
                        single_alt = ref;  // Treat as reference for annotation
                    }

                    std::vector<vep::VariantAnnotation> annotations;

                    // --minimal: trim alleles to minimal representation
                    // Perl VEP trims prefix first, then suffix (matching trim_sequences())
                    std::string ann_ref = ref;
                    std::string ann_alt = single_alt;
                    int ann_pos = pos;
                    if (minimal_mode && (ann_ref.size() > 1 || ann_alt.size() > 1)) {
                        // Perl VEP trim_sequences(): prefix first, then suffix
                        // Prefix: trim while both non-empty and share first base
                        while (!ann_ref.empty() && !ann_alt.empty() &&
                               ann_ref.front() == ann_alt.front()) {
                            ann_ref.erase(ann_ref.begin());
                            ann_alt.erase(ann_alt.begin());
                            ann_pos++;
                        }
                        // Suffix: trim while both non-empty and share last base
                        while (!ann_ref.empty() && !ann_alt.empty() &&
                               ann_ref.back() == ann_alt.back()) {
                            ann_ref.pop_back();
                            ann_alt.pop_back();
                        }
                        // Pad empty alleles (Perl VEP uses "-")
                        if (ann_ref.empty()) ann_ref = "-";
                        if (ann_alt.empty()) ann_alt = "-";
                    }

                    // Check for symbolic alleles (<DEL>, <INS>, etc.) or BND notation
                    bool is_symbolic = (single_alt.size() > 2 && single_alt[0] == '<' && single_alt.back() == '>');
                    bool is_bnd = (single_alt.find('[') != std::string::npos || single_alt.find(']') != std::string::npos);

                    if (is_symbolic || is_bnd) {
                        // Route to structural variant pipeline
                        auto sv = vep::parse_sv_from_vcf(chrom, pos, ref, single_alt, info_map);

                        // Get overlapping transcripts for the SV region
                        auto sv_annotations = annotator.annotate(chrom, pos, ref, ref);
                        // We use annotate() to get transcript list; but we override consequences with SV-specific ones
                        // Instead, get transcripts in the SV region and compute SV consequences directly
                        annotations.clear();

                        // Use the annotator to get transcripts in region via a dummy annotate call
                        // The SV spans from sv.start to sv.end
                        auto region_anns = annotator.annotate(chrom, sv.start, ref, ref);
                        // Collect unique transcript IDs we already saw
                        std::set<std::string> seen_transcripts;
                        for (const auto& ra : region_anns) {
                            if (!ra.transcript_id.empty()) {
                                seen_transcripts.insert(ra.transcript_id);
                            }
                        }

                        // For each transcript found, compute SV consequences
                        if (region_anns.empty() || (region_anns.size() == 1 && region_anns[0].transcript_id.empty())) {
                            // Intergenic SV
                            vep::VariantAnnotation ann;
                            ann.input_variant = line;
                            ann.chromosome = chrom;
                            ann.position = pos;
                            ann.ref_allele = ref;
                            ann.alt_allele = single_alt;
                            ann.consequences.push_back(vep::ConsequenceType::INTERGENIC_VARIANT);
                            ann.impact = vep::Impact::MODIFIER;
                            annotations.push_back(std::move(ann));
                        } else {
                            for (const auto& ra : region_anns) {
                                if (ra.transcript_id.empty()) continue;
                                const vep::Transcript* tx = annotator.get_transcript(ra.transcript_id);
                                if (!tx) continue;

                                auto sv_csqs = vep::get_sv_consequences(sv, *tx);

                                vep::VariantAnnotation ann;
                                ann.input_variant = line;
                                ann.chromosome = chrom;
                                ann.position = pos;
                                ann.ref_allele = ref;
                                ann.alt_allele = single_alt;
                                ann.gene_symbol = ra.gene_symbol;
                                ann.gene_id = ra.gene_id;
                                ann.transcript_id = ra.transcript_id;
                                ann.biotype = ra.biotype;
                                ann.is_canonical = ra.is_canonical;
                                ann.consequences = sv_csqs;
                                // Set impact from most severe SV consequence
                                ann.impact = vep::Impact::MODIFIER;
                                for (const auto& c : sv_csqs) {
                                    vep::Impact ci = vep::get_impact(c);
                                    if (static_cast<int>(ci) < static_cast<int>(ann.impact)) {
                                        ann.impact = ci;
                                    }
                                }
                                annotations.push_back(std::move(ann));
                            }
                        }
                    } else {
                        // Normal small-variant annotation path (use minimized alleles if --minimal)
                        // Use annotate() (all transcripts) when any pick/filter flag needs
                        // the full set for proper transcript selection
                        bool need_all = all_transcripts ||
                            filter_config.pick || filter_config.pick_allele ||
                            filter_config.pick_allele_gene || filter_config.per_gene ||
                            filter_config.flag_pick || filter_config.flag_pick_allele ||
                            filter_config.flag_pick_allele_gene ||
                            filter_config.most_severe ||
                            filter_config.coding_only || filter_config.canonical_only ||
                            filter_config.mane_only || filter_config.gencode_basic ||
                            filter_config.no_intergenic ||
                            filter_config.exclude_predicted ||
                            !filter_config.include_consequences.empty() ||
                            !filter_config.exclude_consequences.empty();
                        // Use cached results from parallel pre-annotation, or annotate inline
                        if (!annotation_cache.empty()) {
                            std::string cache_key = chrom + "\t" + std::to_string(ann_pos) + "\t" + ann_ref + "\t" + ann_alt;
                            auto cache_it = annotation_cache.find(cache_key);
                            if (cache_it != annotation_cache.end()) {
                                annotations = cache_it->second;
                            } else {
                                // Cache miss (SV, non-VCF format, etc.)
                                if (need_all) {
                                    annotations = annotator.annotate(chrom, ann_pos, ann_ref, ann_alt);
                                } else {
                                    annotations.push_back(
                                        annotator.annotate_most_severe(chrom, ann_pos, ann_ref, ann_alt));
                                }
                            }
                        } else {
                            if (need_all) {
                                annotations = annotator.annotate(chrom, ann_pos, ann_ref, ann_alt);
                            } else {
                                auto ann = annotator.annotate_most_severe(chrom, ann_pos, ann_ref, ann_alt);
                                annotations.push_back(ann);
                            }
                        }
                        // If --minimal changed the alleles, mark MINIMISED and update positions
                        if (minimal_mode && (ann_ref != ref || ann_alt != single_alt)) {
                            for (auto& a : annotations) {
                                a.custom_annotations["MINIMISED"] = "1";
                            }
                        }
                    }

                    // Apply transcript filter
                    // Post-process: add MAX_AF, total length, identifiers
                    // Must run before transcript filter so --filter-common can use MAX_AF
                    post_process_annotations(annotations, show_max_af, include_total_length,
                                             include_ccds, include_uniprot, include_xref_refseq,
                                             nearest_mode);

                    if (filter_config.pick || filter_config.pick_allele ||
                        filter_config.pick_allele_gene || filter_config.per_gene ||
                        filter_config.most_severe || filter_config.flag_pick ||
                        filter_config.flag_pick_allele || filter_config.flag_pick_allele_gene ||
                        filter_config.canonical_only || filter_config.mane_only ||
                        filter_config.coding_only || filter_config.gencode_basic ||
                        filter_config.exclude_predicted ||
                        !filter_config.biotypes.empty() || filter_config.no_intergenic ||
                        filter_config.check_frequency ||
                        !filter_config.include_consequences.empty() ||
                        !filter_config.exclude_consequences.empty()) {
                        annotations = transcript_filter.filter(annotations);
                    }

                    // Apply BED custom annotations
                    for (auto& a : annotations) {
                        for (const auto& bed_src : bed_custom_sources) {
                            bed_src.annotate(chrom, pos, a.custom_annotations);
                        }
                    }

                    // Add variant class, allele number, VCF passthrough fields to annotations
                    for (auto& a : annotations) {
                        if (show_variant_class) {
                            if (minimal_mode && (ann_ref != ref || ann_alt != single_alt)) {
                                // Use minimized alleles for variant class
                                std::string vc_ref = (ann_ref == "-") ? "" : ann_ref;
                                std::string vc_alt = (ann_alt == "-") ? "" : ann_alt;
                                a.custom_annotations["VARIANT_CLASS"] = vep::get_variant_class(vc_ref, vc_alt);
                            } else {
                                a.custom_annotations["VARIANT_CLASS"] = vep::get_variant_class(ref, single_alt);
                            }
                        }
                        if (show_allele_number) {
                            a.custom_annotations["ALLELE_NUM"] = std::to_string(allele_num);
                        }
                        // Store existing CSQ for --keep-csq
                        if (keep_csq && !existing_csq.empty()) {
                            a.custom_annotations["EXISTING_CSQ"] = existing_csq;
                        }
                        // Add reference allele info for --show_ref_allele
                        if (show_ref_allele) {
                            a.custom_annotations["GIVEN_REF"] = ref;
                            a.custom_annotations["USED_REF"] = ref;
                        }

                        // Populate Existing_variation from co-located variant database
                        if (check_existing) {
                            auto id_it = a.custom_annotations.find("_colocated:ID");
                            if (id_it != a.custom_annotations.end() && !id_it->second.empty() && id_it->second != ".") {
                                a.existing_variation = id_it->second;
                                // Also extract CLIN_SIG if available
                                auto clnsig_it = a.custom_annotations.find("_colocated:CLNSIG");
                                if (clnsig_it != a.custom_annotations.end() && !clnsig_it->second.empty()) {
                                    a.custom_annotations["CLIN_SIG"] = clnsig_it->second;
                                }
                            }
                        }

                        // Format SIFT output from dbNSFP fields
                        if (!sift_display.empty()) {
                            auto score_it = a.custom_annotations.find("dbnsfp:SIFT_score");
                            auto pred_it = a.custom_annotations.find("dbnsfp:SIFT_pred");
                            if (score_it != a.custom_annotations.end() || pred_it != a.custom_annotations.end()) {
                                std::string sift_val;
                                std::string pred = (pred_it != a.custom_annotations.end()) ? pred_it->second : "";
                                std::string score = (score_it != a.custom_annotations.end()) ? score_it->second : "";
                                // Convert prediction codes to words
                                std::string pred_word = (pred == "D") ? "deleterious" :
                                                       (pred == "T") ? "tolerated" : pred;
                                if (sift_display == "p") {
                                    sift_val = pred_word;
                                } else if (sift_display == "s") {
                                    sift_val = score;
                                } else { // "b" - both
                                    sift_val = pred_word + "(" + score + ")";
                                }
                                a.custom_annotations["SIFT"] = sift_val;
                            }
                        }

                        // Format PolyPhen output from dbNSFP fields
                        if (!polyphen_display.empty()) {
                            std::string pp_score_key = use_humdiv ? "dbnsfp:Polyphen2_HDIV_score" : "dbnsfp:Polyphen2_HVAR_score";
                            std::string pp_pred_key = use_humdiv ? "dbnsfp:Polyphen2_HDIV_pred" : "dbnsfp:Polyphen2_HVAR_pred";
                            auto score_it = a.custom_annotations.find(pp_score_key);
                            auto pred_it = a.custom_annotations.find(pp_pred_key);
                            if (score_it != a.custom_annotations.end() || pred_it != a.custom_annotations.end()) {
                                std::string pp_val;
                                std::string pred = (pred_it != a.custom_annotations.end()) ? pred_it->second : "";
                                std::string score = (score_it != a.custom_annotations.end()) ? score_it->second : "";
                                // Convert prediction codes to words
                                std::string pred_word = (pred == "D") ? "probably_damaging" :
                                                       (pred == "P") ? "possibly_damaging" :
                                                       (pred == "B") ? "benign" : pred;
                                if (polyphen_display == "p") {
                                    pp_val = pred_word;
                                } else if (polyphen_display == "s") {
                                    pp_val = score;
                                } else { // "b" - both
                                    pp_val = pred_word + "(" + score + ")";
                                }
                                a.custom_annotations["PolyPhen"] = pp_val;
                            }
                        }
                        // Preserve original VCF ID, QUAL, FILTER fields
                        a.vcf_id = id;
                        a.vcf_qual = qual_field;
                        a.vcf_filter = filter_field;
                        a.vcf_info = info_field;  // Preserve original INFO for passthrough
                        // Preserve sample columns (FORMAT + samples)
                        a.vcf_sample_columns = sample_columns;
                        // Preserve original VCF REF/ALT for multi-allele output
                        a.vcf_ref = ref;
                        a.vcf_alt = alt;

                        // Compute Ensembl-format display alleles
                        {
                            std::string d_ref, d_alt;
                            int d_start, d_end;

                            if (minimal_mode && (ann_ref != ref || ann_alt != single_alt)) {
                                // --minimal: use minimized alleles directly (already in Ensembl format)
                                d_ref = ann_ref;
                                d_alt = ann_alt;
                                d_start = ann_pos;
                                if (d_ref == "-" || d_ref.empty()) {
                                    d_end = d_start - 1; // Insertion
                                } else {
                                    int mr_len = static_cast<int>(d_ref.size());
                                    d_end = d_start + mr_len - 1;
                                }
                            } else {
                                // Normal: strip shared leading base for VCF indels
                                d_ref = ref;
                                d_alt = single_alt;
                                d_start = pos;
                                d_end = pos + static_cast<int>(ref.size()) - 1;
                                if (d_ref.size() > 1 || d_alt.size() > 1) {
                                    if (!d_ref.empty() && !d_alt.empty() && d_ref[0] == d_alt[0]) {
                                        d_ref = d_ref.substr(1);
                                        d_alt = d_alt.substr(1);
                                        d_start++;
                                        if (d_ref.empty()) {
                                            d_ref = "-";
                                            d_end = d_start - 1; // Insertion
                                        }
                                        if (d_alt.empty()) {
                                            d_alt = "-";
                                        }
                                    }
                                }
                            }
                            a.display_ref = d_ref;
                            a.display_alt = d_alt;
                            a.display_start = d_start;
                            a.display_end = d_end;

                            // Compute position end values for multi-base variants
                            int ref_len = (d_ref != "-" && !d_ref.empty()) ? static_cast<int>(d_ref.size()) : 0;
                            if (ref_len > 1) {
                                if (a.cdna_position > 0)
                                    a.cdna_end = a.cdna_position + ref_len - 1;
                                if (a.cds_position > 0)
                                    a.cds_end = a.cds_position + ref_len - 1;
                                if (a.protein_position > 0 && a.cds_end > 0)
                                    a.protein_end = (a.cds_end - 1) / 3 + 1;
                            }
                        }
                    }

                    // Write annotations
                    for (const auto& ann : annotations) {
                        writer->write_annotation(ann, custom_columns);
                    }

                    variant_count++;
                }

                // For VCF output, flush after all alleles for this position
                if (output_format == vep::OutputFormat::VCF) {
                    vep::VCFWriter* vcf_writer = dynamic_cast<vep::VCFWriter*>(writer.get());
                    if (vcf_writer) {
                        vcf_writer->flush_variant();
                    }
                }

                if (!quiet_mode && !no_progress && variant_count % 1000 == 0) {
                    vep::log(vep::LogLevel::INFO, "Processed " + std::to_string(variant_count) + " variants...");
                }
            }

            // Cleanup
            if (gz_file) {
                gzclose(gz_file);
            }

            // Ensure header is written even for empty input (deferred STDIN case)
            if (!header_written) {
                writer->write_header(custom_columns);
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

            // Write stats to file if --stats-file specified
            if (!stats_file_path.empty()) {
                std::ofstream stats_out(stats_file_path);
                if (stats_out.is_open()) {
                    const auto& stats = writer->get_stats();
                    stats_out << "[VEP run statistics]\n";
                    stats_out << "Variants processed\t" << stats.total_variants << "\n";
                    stats_out << "Variants annotated\t" << stats.annotated_variants << "\n\n";
                    stats_out << "[Consequence counts]\n";
                    for (const auto& pair : stats.consequence_counts) {
                        stats_out << pair.first << "\t" << pair.second << "\n";
                    }
                    stats_out << "\n[Impact counts]\n";
                    for (const auto& pair : stats.impact_counts) {
                        stats_out << pair.first << "\t" << pair.second << "\n";
                    }
                    if (!stats.biotype_counts.empty()) {
                        stats_out << "\n[Biotype counts]\n";
                        for (const auto& pair : stats.biotype_counts) {
                            stats_out << pair.first << "\t" << pair.second << "\n";
                        }
                    }
                    if (!quiet_mode) {
                        std::cerr << "Statistics written to: " << stats_file_path << std::endl;
                    }
                } else {
                    std::cerr << "Warning: Cannot write stats file: " << stats_file_path << std::endl;
                }
            }

            // Print summary if requested
            if (show_summary && !quiet_mode) {
                print_summary(writer->get_stats());
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
                        } else if (hgvs_result.hgvs_type == vep::HGVSType::CODING ||
                                   hgvs_result.hgvs_type == vep::HGVSType::NONCODING) {
                            // Coding HGVS requires transcript mapping
                            // Create a temporary annotator to look up the transcript
                            vep::VEPAnnotator temp_annotator(gtf_path, fasta_path);
                            const vep::Transcript* tx = temp_annotator.get_transcript(hgvs_result.reference_id);
                            if (!tx) {
                                std::cerr << "Error: Transcript not found: " << hgvs_result.reference_id << std::endl;
                                return 1;
                            }

                            // Map CDS position to genomic coordinate
                            int genomic_pos = temp_annotator.map_cds_to_genomic(hgvs_result.start_pos, *tx);
                            if (genomic_pos <= 0) {
                                std::cerr << "Error: Could not map CDS position " << hgvs_result.start_pos
                                          << " to genomic coordinate for " << hgvs_result.reference_id << std::endl;
                                return 1;
                            }

                            // Handle intronic offset
                            if (hgvs_result.intron_offset != 0) {
                                genomic_pos += (tx->strand == '+' ? hgvs_result.intron_offset : -hgvs_result.intron_offset);
                            }

                            chrom = tx->chromosome;
                            pos = genomic_pos;
                            ref = hgvs_result.ref_allele;
                            alt = hgvs_result.alt_allele;

                            // Complement alleles if on minus strand
                            if (tx->strand == '-') {
                                auto complement_char = [](char c) -> char {
                                    switch (c) {
                                        case 'A': return 'T'; case 'T': return 'A';
                                        case 'C': return 'G'; case 'G': return 'C';
                                        case 'a': return 't'; case 't': return 'a';
                                        case 'c': return 'g'; case 'g': return 'c';
                                        default: return c;
                                    }
                                };
                                for (size_t ci = 0; ci < ref.size(); ++ci) ref[ci] = complement_char(ref[ci]);
                                for (size_t ci = 0; ci < alt.size(); ++ci) alt[ci] = complement_char(alt[ci]);
                                std::reverse(ref.begin(), ref.end());
                                std::reverse(alt.begin(), alt.end());
                            }

                            parsed = true;
                            std::cout << "Resolved HGVS: " << chrom << ":" << pos << ":" << ref << ">" << alt << std::endl;
                        } else {
                            // Protein HGVS requires transcript + reverse translation
                            std::cerr << "Note: Protein HGVS (p.) input requires additional transcript mapping.\n"
                                      << "Please use genomic HGVS (g.), coding HGVS (c.), or CHR:POS:REF:ALT format." << std::endl;
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

            // Try SPDI format
            if (!parsed && !variant.empty() && vep::is_spdi_notation(variant)) {
                auto spdi_result = vep::parse_spdi(variant);
                if (spdi_result.valid) {
                    chrom = spdi_result.chromosome;
                    pos = spdi_result.position;
                    ref = spdi_result.ref_allele;
                    alt = spdi_result.alt_allele;
                    parsed = true;
                    std::cout << "Parsed SPDI: " << chrom << ":" << pos << ":" << ref << ">" << alt << std::endl;
                } else {
                    std::cerr << "Error: Invalid SPDI notation: " << spdi_result.error_message << std::endl;
                    return 1;
                }
            }

            // Try Ensembl default format
            if (!parsed && !variant.empty() && vep::is_ensembl_format(variant)) {
                auto ensembl_result = vep::parse_ensembl_format(variant);
                if (ensembl_result.valid) {
                    chrom = ensembl_result.chromosome;
                    pos = ensembl_result.position;
                    ref = ensembl_result.ref_allele;
                    alt = ensembl_result.alt_allele;
                    parsed = true;
                    std::cout << "Parsed Ensembl format: " << chrom << ":" << pos << ":" << ref << ">" << alt << std::endl;
                } else {
                    std::cerr << "Error: Invalid Ensembl format: " << ensembl_result.error_message << std::endl;
                    return 1;
                }
            }

            // Try REST-style region format (CHR:START-END:STRAND/ALLELE)
            if (!parsed && !variant.empty() && vep::is_rest_region_format(variant)) {
                auto rest_result = vep::parse_rest_region(variant);
                if (rest_result.valid) {
                    chrom = rest_result.chromosome;
                    pos = rest_result.position;
                    ref = rest_result.ref_allele;
                    alt = rest_result.alt_allele;
                    parsed = true;
                    std::cout << "Parsed REST region: " << chrom << ":" << pos << ":" << ref << ">" << alt << std::endl;
                }
            }

            // Try regular variant format (CHR:POS:REF:ALT)
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
            configure_annotator(annotator);

            // Set up all annotation sources (returns initial custom_columns)
            std::vector<std::string> custom_columns = setup_annotation_sources(annotator);

            std::cout << "\nDatabase loaded:\n" << annotator.get_stats() << std::endl;

            // Check reference allele if requested (single-variant mode)
            if (check_ref) {
                const vep::ReferenceGenome* ref_genome = annotator.get_reference();
                if (ref_genome) {
                    if (!verify_reference(*ref_genome, chrom, pos, ref, quiet_mode)) {
                        std::cerr << "Error: Reference allele mismatch. Aborting." << std::endl;
                        return 1;
                    }
                }
            }

            std::vector<vep::VariantAnnotation> annotations;

            // Apply --minimal allele trimming to single-variant input
            std::string sv_ref = ref, sv_alt = alt;
            int sv_pos = pos;
            bool was_minimised = false;
            if (minimal_mode && (sv_ref.size() > 1 || sv_alt.size() > 1)) {
                // Trim common prefix (match Perl VEP trim_sequences: trim until empty)
                while (!sv_ref.empty() && !sv_alt.empty() && sv_ref.front() == sv_alt.front()) {
                    sv_ref = sv_ref.substr(1);
                    sv_alt = sv_alt.substr(1);
                    sv_pos++;
                    was_minimised = true;
                }
                // Trim common suffix
                while (!sv_ref.empty() && !sv_alt.empty() && sv_ref.back() == sv_alt.back()) {
                    sv_ref.pop_back();
                    sv_alt.pop_back();
                    was_minimised = true;
                }
                if (sv_ref.empty()) { sv_ref = "-"; was_minimised = true; }
                if (sv_alt.empty()) { sv_alt = "-"; was_minimised = true; }
                if (was_minimised) {
                    ref = sv_ref;
                    alt = sv_alt;
                    pos = sv_pos;
                }
            }

            // Check if this is a structural variant (symbolic allele or BND)
            bool is_symbolic_sv = (alt.size() > 2 && alt[0] == '<' && alt.back() == '>');
            bool is_bnd_sv = (alt.find('[') != std::string::npos || alt.find(']') != std::string::npos);

            if (is_symbolic_sv || is_bnd_sv) {
                // Structural variant annotation path
                std::map<std::string, std::string> info_map;
                auto sv = vep::parse_sv_from_vcf(chrom, pos, ref, alt, info_map);

                auto region_anns = annotator.annotate(chrom, sv.start, ref, ref);
                if (region_anns.empty() || (region_anns.size() == 1 && region_anns[0].transcript_id.empty())) {
                    vep::VariantAnnotation ann;
                    ann.chromosome = chrom;
                    ann.position = pos;
                    ann.ref_allele = ref;
                    ann.alt_allele = alt;
                    ann.consequences.push_back(vep::ConsequenceType::INTERGENIC_VARIANT);
                    ann.impact = vep::Impact::MODIFIER;
                    annotations.push_back(std::move(ann));
                } else {
                    for (const auto& ra : region_anns) {
                        if (ra.transcript_id.empty()) continue;
                        const vep::Transcript* tx = annotator.get_transcript(ra.transcript_id);
                        if (!tx) continue;

                        auto sv_csqs = vep::get_sv_consequences(sv, *tx);

                        vep::VariantAnnotation ann;
                        ann.chromosome = chrom;
                        ann.position = pos;
                        ann.ref_allele = ref;
                        ann.alt_allele = alt;
                        ann.gene_symbol = ra.gene_symbol;
                        ann.gene_id = ra.gene_id;
                        ann.transcript_id = ra.transcript_id;
                        ann.biotype = ra.biotype;
                        ann.is_canonical = ra.is_canonical;
                        ann.consequences = sv_csqs;
                        ann.impact = vep::Impact::MODIFIER;
                        for (const auto& c : sv_csqs) {
                            vep::Impact ci = vep::get_impact(c);
                            if (static_cast<int>(ci) < static_cast<int>(ann.impact)) {
                                ann.impact = ci;
                            }
                        }
                        annotations.push_back(std::move(ann));
                    }
                }
            } else if (all_transcripts ||
                       filter_config.pick || filter_config.pick_allele ||
                       filter_config.pick_allele_gene || filter_config.per_gene ||
                       filter_config.flag_pick || filter_config.flag_pick_allele ||
                       filter_config.flag_pick_allele_gene) {
                annotations = annotator.annotate(chrom, pos, ref, alt);
            } else {
                auto ann = annotator.annotate_most_severe(chrom, pos, ref, alt);
                annotations.push_back(ann);
            }

            // Mark annotations as minimised if alleles were trimmed
            if (was_minimised) {
                for (auto& a : annotations) {
                    a.custom_annotations["MINIMISED"] = "1";
                }
            }

            // Post-process: add MAX_AF, total length, identifiers
            // Must run before transcript filter so --filter-common can use MAX_AF
            post_process_annotations(annotations, show_max_af, include_total_length,
                                     include_ccds, include_uniprot, include_xref_refseq,
                                     nearest_mode);

            // Apply transcript filter
            if (filter_config.pick || filter_config.pick_allele ||
                filter_config.pick_allele_gene || filter_config.per_gene ||
                filter_config.most_severe || filter_config.flag_pick ||
                filter_config.flag_pick_allele || filter_config.flag_pick_allele_gene ||
                filter_config.canonical_only || filter_config.mane_only ||
                filter_config.coding_only || filter_config.gencode_basic ||
                filter_config.exclude_predicted ||
                !filter_config.biotypes.empty() || filter_config.no_intergenic ||
                filter_config.check_frequency ||
                !filter_config.include_consequences.empty() ||
                !filter_config.exclude_consequences.empty()) {
                annotations = transcript_filter.filter(annotations);
            }

            // Apply BED custom annotations
            for (auto& a : annotations) {
                for (const auto& bed_src : bed_custom_sources) {
                    bed_src.annotate(chrom, pos, a.custom_annotations);
                }
            }

            // Add variant class for single variant mode
            for (auto& a : annotations) {
                if (show_variant_class) {
                    std::string vc_ref = (ref == "-") ? "" : ref;
                    std::string vc_alt = (alt == "-") ? "" : alt;
                    a.custom_annotations["VARIANT_CLASS"] = vep::get_variant_class(vc_ref, vc_alt);
                }
                if (show_allele_number) {
                    a.custom_annotations["ALLELE_NUM"] = "1";
                }
                if (show_ref_allele) {
                    a.custom_annotations["GIVEN_REF"] = ref;
                    a.custom_annotations["USED_REF"] = ref;
                }

                // Format SIFT output from dbNSFP fields
                if (!sift_display.empty()) {
                    auto score_it = a.custom_annotations.find("dbnsfp:SIFT_score");
                    auto pred_it = a.custom_annotations.find("dbnsfp:SIFT_pred");
                    if (score_it != a.custom_annotations.end() || pred_it != a.custom_annotations.end()) {
                        std::string pred = (pred_it != a.custom_annotations.end()) ? pred_it->second : "";
                        std::string score = (score_it != a.custom_annotations.end()) ? score_it->second : "";
                        std::string pred_word = (pred == "D") ? "deleterious" :
                                               (pred == "T") ? "tolerated" : pred;
                        std::string sift_val;
                        if (sift_display == "p") sift_val = pred_word;
                        else if (sift_display == "s") sift_val = score;
                        else sift_val = pred_word + "(" + score + ")";
                        a.custom_annotations["SIFT"] = sift_val;
                    }
                }

                // Format PolyPhen output from dbNSFP fields
                if (!polyphen_display.empty()) {
                    std::string pp_score_key = use_humdiv ? "dbnsfp:Polyphen2_HDIV_score" : "dbnsfp:Polyphen2_HVAR_score";
                    std::string pp_pred_key = use_humdiv ? "dbnsfp:Polyphen2_HDIV_pred" : "dbnsfp:Polyphen2_HVAR_pred";
                    auto score_it = a.custom_annotations.find(pp_score_key);
                    auto pred_it = a.custom_annotations.find(pp_pred_key);
                    if (score_it != a.custom_annotations.end() || pred_it != a.custom_annotations.end()) {
                        std::string pred = (pred_it != a.custom_annotations.end()) ? pred_it->second : "";
                        std::string score = (score_it != a.custom_annotations.end()) ? score_it->second : "";
                        std::string pred_word = (pred == "D") ? "probably_damaging" :
                                               (pred == "P") ? "possibly_damaging" :
                                               (pred == "B") ? "benign" : pred;
                        std::string pp_val;
                        if (polyphen_display == "p") pp_val = pred_word;
                        else if (polyphen_display == "s") pp_val = score;
                        else pp_val = pred_word + "(" + score + ")";
                        a.custom_annotations["PolyPhen"] = pp_val;
                    }
                }

                // Populate Existing_variation from co-located variant database
                if (check_existing) {
                    auto id_it = a.custom_annotations.find("_colocated:ID");
                    if (id_it != a.custom_annotations.end() && !id_it->second.empty() && id_it->second != ".") {
                        a.existing_variation = id_it->second;
                        auto clnsig_it = a.custom_annotations.find("_colocated:CLNSIG");
                        if (clnsig_it != a.custom_annotations.end() && !clnsig_it->second.empty()) {
                            a.custom_annotations["CLIN_SIG"] = clnsig_it->second;
                        }
                    }
                }
            }

            // Compute display alleles, set VCF defaults, and input_variant for all annotations
            for (auto& a : annotations) {
                a.input_variant = chrom + ":" + std::to_string(pos) + ":" + ref + ":" + alt;
                // Set default VCF fields for single-variant mode (no VCF input)
                if (a.vcf_id.empty()) a.vcf_id = ".";
                if (a.vcf_qual.empty()) a.vcf_qual = ".";
                if (a.vcf_filter.empty()) a.vcf_filter = ".";
                // Compute Ensembl-format display alleles
                std::string d_ref = ref;
                std::string d_alt = alt;
                int d_start = pos;
                int d_end = pos + static_cast<int>(ref.size()) - 1;
                // Handle dash-encoded alleles from --minimal mode
                if (d_ref == "-" || d_ref.empty()) {
                    d_end = d_start - 1; // Insertion: end < start
                }
                if (d_ref.size() > 1 || d_alt.size() > 1) {
                    if (!d_ref.empty() && !d_alt.empty() && d_ref[0] == d_alt[0]) {
                        d_ref = d_ref.substr(1);
                        d_alt = d_alt.substr(1);
                        d_start++;
                        if (d_ref.empty()) { d_ref = "-"; d_end = d_start - 1; }
                        if (d_alt.empty()) { d_alt = "-"; }
                    }
                }
                a.display_ref = d_ref;
                a.display_alt = d_alt;
                a.display_start = d_start;
                a.display_end = d_end;

                // Compute position end values for multi-base variants
                int ref_len = (d_ref != "-" && !d_ref.empty()) ? static_cast<int>(d_ref.size()) : 0;
                if (ref_len > 1) {
                    if (a.cdna_position > 0)
                        a.cdna_end = a.cdna_position + ref_len - 1;
                    if (a.cds_position > 0)
                        a.cds_end = a.cds_position + ref_len - 1;
                    if (a.protein_position > 0 && a.cds_end > 0)
                        a.protein_end = (a.cds_end - 1) / 3 + 1;
                }
            }

            // Add flag-conditional custom columns
            build_custom_columns(custom_columns);

            // Create and configure output writer
            std::string sv_output_path = output_path.empty() ? "" : output_path;
            std::unique_ptr<vep::OutputWriter> sv_writer = configure_writer(sv_output_path);

            sv_writer->write_header(custom_columns);
            for (const auto& ann : annotations) {
                sv_writer->write_annotation(ann, custom_columns);
            }
            if (output_format == vep::OutputFormat::VCF) {
                vep::VCFWriter* vcf_w = dynamic_cast<vep::VCFWriter*>(sv_writer.get());
                if (vcf_w) vcf_w->flush_variant();
            }
            sv_writer->write_footer();
            sv_writer->close();

            if (!quiet_mode && !output_path.empty()) {
                std::cout << "Annotation complete. " << annotations.size()
                          << " annotation(s) written to: " << output_path << std::endl;
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
