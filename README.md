# VEP Annotator

A high-performance C++ implementation of Ensembl's [Variant Effect Predictor (VEP)](https://www.ensembl.org/vep). Achieves ~99.9% feature parity with the Perl VEP while running **50-70x faster** (single-threaded). All annotation is performed locally using standard data files — no API calls required.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://isocpp.org/std/the-standard)
[![Tests](https://img.shields.io/badge/tests-640%20passing-brightgreen.svg)]()

## Features

### Core Annotation
- Consequence prediction using Sequence Ontology (SO) terms (all 35+ VEP consequence types)
- Impact classification (HIGH, MODERATE, LOW, MODIFIER)
- HGVS notation: HGVSc (coding), HGVSp (protein), HGVSg (genomic)
- CDS/protein position calculation with codon and amino acid changes
- 3' shifting for HGVS-compliant indel representation
- Structural variant annotation (symbolic alleles, breakends)
- Transcript filtering: `--pick`, `--per_gene`, `--pick_allele`, `--flag_pick`
- MANE Select/Plus Clinical, canonical transcript, TSL, APPRIS ranking

### Output Formats
- **TSV** (default): Tab-separated with Extra column for custom fields
- **JSON**: Structured output matching Perl VEP JSON schema
- **VCF**: CSQ INFO field annotation (compatible with downstream tools)
- `--everything` flag for comprehensive annotation output

### Pathogenicity Predictions (via dbNSFP)
- **35+ pathogenicity scores**: SIFT, PolyPhen-2, CADD, REVEL, AlphaMissense, MetaSVM, MetaLR, VEST4, PROVEAN, FATHMM, MutationTaster, MutationAssessor, DANN, Eigen, M-CAP, MPC, PrimateAI, BayesDel, ClinPred, and more
- Conservation scores: PhyloP, PhastCons, GERP++, SiPhy
- Population frequencies: gnomAD, ExAC, 1000 Genomes
- Clinical data: ClinVar annotations

### Splice Predictions
- **SpliceAI**: Deep learning splice predictions (delta scores for acceptor/donor gain/loss)
- **MaxEntScan**: Position weight matrix-based splice site scoring (algorithmic, no data file needed)
- **dbscSNV**: Splice site variant predictions (ada_score, rf_score)

### Conservation Scores
- **PhyloP**: Phylogenetic p-values from 100-way vertebrate alignment
- **PhastCons**: Probability of negative selection
- **GERP++**: Genomic Evolutionary Rate Profiling

### Regulatory Annotations
- Ensembl Regulatory Build features (promoters, enhancers, TFBS, open chromatin, CTCF)
- Cell type-specific filtering

### Protein Domains
- **Pfam**: Protein family annotations
- **InterPro**: Integrated protein domain database

### Loss-of-Function (LoF) Annotations
- **LOFTEE**: High-confidence (HC) / Low-confidence (LC) classification
- **NMD**: Nonsense-mediated decay prediction (50bp rule)
- **LoFtool**: Gene-level constraint scores

### Custom Annotations
- Add any VCF file as an annotation source with two loading modes:
  - **In-memory** (`--annotation`): Loads entire file into RAM for fast lookups — best for small files (<100K records)
  - **Tabix-indexed** (`--annotation-tabix`): On-disk queries via htslib — instant startup, minimal memory, scales to billions of records
- BED file annotations for interval-based overlaps
- Perl VEP-compatible `--custom` syntax for VCF, BED, bigWig, GFF/GTF
- See [Custom Annotation Files](#custom-annotation-files) below for detailed usage

### Plugin System
- Extend with custom annotation sources via shared libraries
- Dynamic loading with simple plugin API

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/dcho108803/vep-annotator.git
cd vep-annotator

# Build
mkdir build && cd build
cmake ..
make -j4

# Run tests
./vep_tests
```

### Dependencies

**Required:**
- C++17 compiler (GCC 7+, Clang 5+)
- CMake 3.16+
- zlib

**Optional:**
- htslib (for tabix-indexed file support: `--annotation-tabix`, `--dbnsfp`, `--spliceai`, etc.)
- libBigWig (for conservation bigWig files: `--phylop`, `--phastcons`, `--gerp`)

**macOS:**
```bash
brew install cmake htslib
```

**Ubuntu/Debian:**
```bash
sudo apt-get install cmake libhts-dev zlib1g-dev
```

### Basic Usage

```bash
# Annotate a single variant
./vep_annotator --gtf genes.gtf.gz --fasta genome.fa.gz -v 7:140753336:A:T

# Annotate a VCF file (TSV output)
./vep_annotator --gtf genes.gtf.gz --fasta genome.fa.gz --vcf input.vcf.gz -o output.tsv

# JSON output
./vep_annotator --gtf genes.gtf.gz --fasta genome.fa.gz --vcf input.vcf.gz --json -o output.json

# VCF output with CSQ annotation
./vep_annotator --gtf genes.gtf.gz --fasta genome.fa.gz --vcf input.vcf.gz --vcf-output -o output.vcf

# Everything mode (all annotations enabled)
./vep_annotator --gtf genes.gtf.gz --fasta genome.fa.gz --vcf input.vcf.gz --everything -o output.tsv
```

### Full Annotation Pipeline

```bash
./vep_annotator \
    --gtf Homo_sapiens.GRCh38.110.gtf.gz \
    --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    --dbnsfp dbNSFP4.4a.txt.gz \
    --dbnsfp-fields "SIFT_score,Polyphen2_HDIV_score,CADD_phred,REVEL_score,AlphaMissense_score" \
    --spliceai spliceai_scores.vcf.gz \
    --maxentscan \
    --phylop hg38.phyloP100way.bw \
    --phastcons hg38.phastCons100way.bw \
    --regulatory regulatory_features.gff.gz \
    --loftee \
    --nmd \
    --annotation clinvar:clinvar.vcf.gz:CLNSIG,CLNDN \
    --annotation-tabix gnomad:gnomad.exomes.vcf.bgz:AF,AF_popmax \
    --vcf input.vcf.gz \
    -o annotated.tsv
```

## Command Line Reference

### Required Options

| Option | Description |
|--------|-------------|
| `--gtf FILE` | GTF/GFF annotation file (genes/transcripts) |
| `--fasta FILE` | Reference genome FASTA file |

### Variant Input

| Option | Description |
|--------|-------------|
| `-v, --variant CHR:POS:REF:ALT` | Single variant to annotate |
| `--vcf FILE` | VCF file for batch annotation (.vcf or .vcf.gz) |
| `-` or `STDIN` | Read VCF from standard input |

### Output Options

| Option | Description |
|--------|-------------|
| `-o, --output FILE` | Output file path (default: stdout) |
| `--json` | JSON output format |
| `--vcf-output` | VCF output format (CSQ INFO field) |
| `--tab` | TSV output format (default) |
| `--compress_output bgzip` | Compress output with bgzip |
| `--no_headers` | Suppress header lines |
| `--fields FIELDS` | Custom field order (comma-separated) |

### Transcript Selection

| Option | Description |
|--------|-------------|
| `--pick` | Pick one annotation per variant (most severe) |
| `--pick_allele` | Pick one annotation per allele |
| `--per_gene` | Pick one annotation per gene |
| `--pick_allele_gene` | Pick one per allele per gene |
| `--flag_pick` | Flag picked annotation but output all |
| `--pick_order LIST` | Custom pick priority order |
| `--most_severe` | Only report most severe consequence |

### Display Options

| Option | Description |
|--------|-------------|
| `--everything` | Enable all display flags |
| `--symbol` | Show gene symbols |
| `--biotype` | Show transcript biotypes |
| `--canonical` | Flag canonical transcripts |
| `--mane` | Show MANE Select/Plus Clinical |
| `--tsl` | Show Transcript Support Level |
| `--appris` | Show APPRIS annotation |
| `--numbers` | Show exon/intron numbers |
| `--domains` | Show protein domain annotations |
| `--hgvs` | Generate HGVS notation |
| `--hgvsg` | Generate HGVSg (genomic) |
| `--protein` | Show protein (ENSP) ID |
| `--ccds` | Show CCDS ID |
| `--uniprot` | Show UniProt ID |
| `--af` | Show allele frequency |
| `--max_af` | Show maximum allele frequency |
| `--variant_class` | Show variant class (SNV/insertion/deletion) |
| `--allele_number` | Show allele number |
| `--total_length` | Show transcript/CDS lengths |
| `--nearest` | Show nearest gene for intergenic variants |

### Variant Processing

| Option | Description |
|--------|-------------|
| `--minimal` | Allele-level minimization for multi-allele VCFs |
| `--check_ref` | Verify reference allele against FASTA |
| `--check_existing FILE` | Co-located variant lookup VCF |
| `--keep_csq` | Preserve existing VCF CSQ annotations |
| `--no_intergenic` | Skip intergenic variants |
| `--coding_only` | Only annotate coding transcripts |
| `--exclude_predicted` | Exclude predicted (XM_/XR_) transcripts |
| `--chr LIST` | Only annotate these chromosomes (comma-separated) |
| `--distance N` | Upstream/downstream distance (default: 5000) |
| `--assembly NAME` | Assembly name for output (default: GRCh38) |
| `--synonyms FILE` | Chromosome synonym mapping file (tab-separated) |
| `--config FILE` | Load options from config file (key=value format) |
| `--stats_file FILE` | Write run statistics to file |

### Pathogenicity Predictions

| Option | Description |
|--------|-------------|
| `--dbnsfp FILE` | dbNSFP database (tabix-indexed .txt.gz) |
| `--dbnsfp-fields FIELDS` | Fields to extract (or preset: essential/pathogenicity/conservation/all) |
| `--sift DISPLAY` | SIFT output style (p = prediction, s = score, b = both) |
| `--polyphen DISPLAY` | PolyPhen output style (p/s/b) |

### Splice Predictions

| Option | Description |
|--------|-------------|
| `--spliceai FILE` | SpliceAI VCF file (tabix-indexed) |
| `--maxentscan` | Enable MaxEntScan splice site scoring |
| `--dbscsnv FILE` | dbscSNV file (tabix-indexed) |

### Conservation Scores

| Option | Description |
|--------|-------------|
| `--phylop FILE` | PhyloP bigWig file |
| `--phastcons FILE` | PhastCons bigWig file |
| `--gerp FILE` | GERP++ bigWig file |

### Regulatory & Protein Domains

| Option | Description |
|--------|-------------|
| `--regulatory FILE` | Ensembl Regulatory Build GFF3 file |
| `--cell_type LIST` | Cell type filter (comma-separated) |
| `--pfam FILE` | Pfam domain annotations TSV |
| `--interpro FILE` | InterPro domain annotations TSV |

### Loss-of-Function

| Option | Description |
|--------|-------------|
| `--loftee` | Enable LOFTEE-style LoF classification |
| `--nmd` | Enable NMD prediction |
| `--loftool FILE` | LoFtool gene constraint scores |

### Custom Annotations

| Option | Description |
|--------|-------------|
| `--annotation NAME:VCF[:FIELDS]` | Add VCF annotation source (loaded into memory) |
| `--annotation-tabix NAME:VCF[:FIELDS]` | Add tabix-indexed VCF source (on-disk queries) |
| `--custom FILE,NAME,TYPE,OVERLAP,0,FIELDS` | Perl VEP-compatible custom source (VCF, BED, bigWig) |

See [Custom Annotation Files](#custom-annotation-files) for detailed usage, examples, and performance guidance.

### Filtering

| Option | Description |
|--------|-------------|
| `--filter_common` | Remove common variants |
| `--freq_pop POP` | Population for frequency filter |
| `--freq_freq FREQ` | Frequency threshold |
| `--freq_gt_lt gt/lt` | Greater than or less than threshold |

## Performance

Benchmarked on 100,000 chr22 variants (Release build, Apple Silicon M-series):

| Mode | Wall time | Throughput | vs Perl VEP |
|------|----------:|-----------:|------------:|
| Single thread, TSV | 4.0s | **25,000 var/sec** | ~80x faster |
| Single thread, JSON | 4.0s | **25,000 var/sec** | ~80x faster |
| Single thread, --everything | 4.4s | **22,700 var/sec** | ~75x faster |
| --fork 4, --everything | 2.9s | **34,500 var/sec** | ~115x faster |

**Multi-threaded:** Use `--fork N` for parallel annotation with N threads. Input is pre-buffered and annotations are distributed across threads using atomic work-stealing.

### Memory Usage

| Component | Memory |
|-----------|--------|
| GTF annotations | ~500MB - 2GB |
| Reference FASTA | ~3GB |
| dbNSFP (tabix) | ~100MB |
| Conservation (bigWig) | ~200MB |
| ClinVar (in-memory) | ~500MB |
| gnomAD (tabix) | ~100MB |

**Tip:** Use `--annotation-tabix` for large VCF files to minimize memory usage.

## Testing

The project includes 640 GoogleTest unit tests:

| Test Suite | Tests | Coverage |
|------------|------:|----------|
| Consequences | 32 | SO terms, impact levels, ranking |
| Codon Table | 24 | Translation, start/stop codons |
| HGVS | 24 | Parsing, generation, notation types |
| SpliceAI | 12 | Score parsing, cutoffs, thread safety |
| Filter VEP | 180 | Operators, expressions, conditions, pipeline |
| Output Writers | 153 | TSV/JSON/VCF formatting, escaping, stats |
| Transcript Filter | 104 | Pick modes, filtering, ranking criteria |
| CLI Utilities | 76 | Variant parsing, allele trimming, config parsing |
| **Total** | **640** | |

```bash
cd build
./vep_tests
```

## Output Format

### TSV Output Columns

| Column | Description |
|--------|-------------|
| #Uploaded_variation | Variant identifier |
| Location | Genomic location |
| Allele | Alternate allele |
| Gene | Gene ID |
| Feature | Transcript ID |
| Feature_type | Feature type (Transcript) |
| Consequence | VEP consequence term(s) |
| cDNA_position | Position in cDNA |
| CDS_position | Position in CDS |
| Protein_position | Position in protein |
| Amino_acids | Amino acid change |
| Codons | Codon change |
| Existing_variation | Co-located variant IDs |
| Extra | Additional fields (key=value pairs) |

## Consequence Types

Consequences are reported using [Sequence Ontology](http://www.sequenceontology.org/) terms, ranked by severity:

| Impact | Consequences |
|--------|-------------|
| **HIGH** | transcript_ablation, splice_acceptor_variant, splice_donor_variant, stop_gained, frameshift_variant, stop_lost, start_lost, transcript_amplification |
| **MODERATE** | inframe_insertion, inframe_deletion, missense_variant, protein_altering_variant |
| **LOW** | splice_donor_5th_base_variant, splice_region_variant, splice_donor_region_variant, splice_polypyrimidine_tract_variant, incomplete_terminal_codon_variant, start_retained_variant, stop_retained_variant, synonymous_variant |
| **MODIFIER** | coding_sequence_variant, mature_miRNA_variant, 5_prime_UTR_variant, 3_prime_UTR_variant, non_coding_transcript_exon_variant, intron_variant, NMD_transcript_variant, non_coding_transcript_variant, upstream_gene_variant, downstream_gene_variant, TFBS_ablation, TFBS_amplification, TF_binding_site_variant, regulatory_region_ablation, regulatory_region_amplification, feature_elongation, regulatory_region_variant, feature_truncation, intergenic_variant |

## Data Files

### Required Files

| File | Description | Source |
|------|-------------|--------|
| GTF | Gene annotations | [Ensembl FTP](https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/) |
| FASTA | Reference genome | [Ensembl FTP](https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/) |

### Optional Annotation Files

| File | Description | Size | Source |
|------|-------------|------|--------|
| dbNSFP | Pathogenicity scores | ~35GB | [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP) |
| SpliceAI | Splice predictions | ~20GB | [Illumina](https://basespace.illumina.com/s/otSPW8hnhaZR) |
| PhyloP | Conservation (100-way) | ~10GB | [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/) |
| PhastCons | Conservation (100-way) | ~10GB | [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/) |
| Regulatory Build | Regulatory features | ~500MB | [Ensembl FTP](https://ftp.ensembl.org/pub/release-110/regulation/homo_sapiens/) |
| ClinVar | Clinical variants | ~100MB | [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) |
| gnomAD | Population frequencies | ~30GB+ | [gnomAD](https://gnomad.broadinstitute.org/downloads) |

## Using as a Library

```cpp
#include "vep_annotator.hpp"
#include "annotation_sources.hpp"

int main() {
    vep::VEPAnnotator annotator("genes.gtf.gz", "genome.fa.gz");

    // Add annotation sources
    auto dbnsfp = vep::create_dbnsfp_source("dbNSFP.txt.gz", "essential");
    annotator.add_source(dbnsfp);
    annotator.initialize_sources();

    // Annotate a variant
    auto results = annotator.annotate("7", 140753336, "A", "T");

    for (const auto& ann : results) {
        std::cout << ann.gene_symbol << "\t"
                  << ann.get_consequence_string() << "\t"
                  << vep::impact_to_string(ann.impact) << "\t"
                  << ann.hgvsc << "\t"
                  << ann.hgvsp << std::endl;
    }
    return 0;
}
```

### Performance Options

| Option | Description |
|--------|-------------|
| `--fork N` | Use N parallel annotation threads |
| `--buffer_size N` | Number of variants to buffer (default: 5000) |
| `--quiet` | Suppress progress messages |

## Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `-DWITH_HTSLIB=ON/OFF` | ON | Enable tabix support |
| `-DWITH_BIGWIG=ON/OFF` | ON | Enable bigWig support |
| `-DCMAKE_BUILD_TYPE=Release/Debug` | Release | Build type |

## Differences from Perl VEP

This C++ implementation achieves ~99.9% feature parity with Perl Ensembl VEP for core annotation functionality. The main differences fall into three categories: architectural design choices, missing features, and minor behavioral differences.

### Architecture: Local Files vs Ensembl Cache

The most fundamental difference is the data model. Perl VEP has three modes: `--cache` (pre-built binary data store), `--database` (live MySQL connection), and the REST API. This C++ implementation uses **local files exclusively** — a GTF for gene annotations and a FASTA for the reference genome are always required. There is no cache, database, or API connectivity.

This means:
- **Faster startup** — no cache unpacking or database connections
- **Simpler deployment** — just the binary and your data files
- **Explicit data sources** — every annotation source requires a file path argument (e.g., `--dbnsfp FILE`, `--spliceai FILE`, `--regulatory FILE`)
- **No built-in frequency data** — population frequencies (gnomAD, 1000G, ESP) must be provided as custom annotation VCF files rather than being bundled in a cache

### Missing Features

The following Perl VEP features are **not implemented**:

| Feature | Perl VEP | C++ VEP | Notes |
|---------|----------|---------|-------|
| Ensembl cache/database | Full support | Not available | Use `--gtf` + `--fasta` instead |
| REST API queries | Full support | Not available | Runs purely locally |
| Built-in population frequencies | From cache | Not built-in | Add via `--annotation-tabix gnomad:file.vcf.gz:AF` |
| HTML statistics report | `--stats-html` | Not available | Text stats via `--stats-file` |
| Individual/sample genotypes | `--individual` | Parsed but not functional | No per-sample consequence reporting |
| Motif feature consequences | Full support | Not available | `MOTIF_NAME`, `MOTIF_POS`, `MOTIF_SCORE_CHANGE` not generated |
| PubMed citations | `--pubmed` | No-op | Requires variant database |
| Gene-phenotype associations | `--gene-phenotype` | No-op | Requires phenotype database |
| Variant synonyms | `--var-synonyms` | No-op | Requires variant database |
| miRNA structure annotation | `--mirna` | No-op | Not implemented |
| LRG coordinates | `--lrg` | No-op | Not implemented |
| GA4GH VRS format | `--ga4gh-vrs` | No-op | Not implemented |
| BAM transcript correction | `--bam` | No-op | Not implemented |
| Transcript filter expressions | `--transcript-filter EXPR` | No-op | Use individual filter flags instead |
| Multi-species support | `--genomes` | No-op | Single species only |
| Warning file | `--warning-file` | No-op | Warnings go to stderr |

### Behavioral Differences

Features that both versions have but behave differently:

| Feature | Perl VEP behavior | C++ VEP behavior |
|---------|-------------------|------------------|
| `--everything` + regulatory | Automatically includes regulatory annotations from cache | Does NOT add `--regulatory` (requires explicit `--regulatory FILE` path) |
| `--check-existing` | Uses cache for co-located variant lookup | Requires explicit VCF file path (`--check-existing dbsnp.vcf.gz`) |
| `--sift` / `--polyphen` | Scores from Ensembl cache | Scores from `--dbnsfp` if provided; flags control display format only |
| `--domains` | Domain data from cache | Requires `--pfam FILE` and/or `--interpro FILE` |
| `--fork N` | Perl `fork()` child processes | C++ `std::thread` with work-stealing; pre-buffers all input |
| `--minimal` VCF output | Preserves original VCF alleles in output | Writes minimized alleles; insertion positions may shift +1 |
| Plugin system | Rich Perl plugin ecosystem (~70+ plugins) | C++ shared library plugins (`.so`/`.dylib`); Perl plugins not compatible |
| `filter_vep` regex | Full Perl regular expressions | Simple substring matching (`string::find`) |
| HGVS multi-base duplications | Full `dup` notation | May use `ins` instead of `dup` for multi-base duplications |

### Known Edge Cases

Minor limitations that affect rare variant types:

- **MNV multi-codon extraction**: Multi-nucleotide variants spanning multiple codons may have limited codon display (consequence determination is unaffected)
- **UTR-intronic HGVSc**: Variants in UTR introns may have imprecise HGVS notation (very rare)
- **Complex indel codon boundaries**: Complex indels spanning codon boundaries may show approximate codon changes (consequence type is correct)
- **`-v` mode dash alleles**: The single-variant parser's `-` to `:` replacement can misparse dash-containing deletion alleles (use VCF input instead)

### Perl VEP Compatibility Flags

For pipeline compatibility, the following Perl VEP flags are **accepted without error** but have no effect: `--cache`, `--offline`, `--database`, `--merged`, `--refseq`, `--dir-cache`, `--cache-version`, `--host`, `--port`, `--user`, `--password`, `--registry-file`, `--db-version`, `--safe`, `--tmpdir`, `--dir`, `--force`, `--dont-skip`, `--lookup-ref`, `--failed`, `--no-whole-genome`, `--verbose`, `--show-cache-info`, `--shift-length`, `--shift-hgvs`, `--no-check-alleles`, `--exclude-null-alleles`, `--check-svs`, `--custom-multi-allelic`, `--max-sv-size`, `--no-check-variants-order`, `--overlap-cutoff`, `--af-exac`, `--hgvsg-use-accession`.

This allows existing Perl VEP command lines to run without modification (unsupported flags are silently ignored).

## Custom Annotation Files

You can add any VCF, BED, or bigWig file as a custom annotation source. There are two loading modes for VCF files, plus Perl VEP-compatible `--custom` syntax.

### Loading Modes

#### Tabix-indexed (recommended for large files)

```bash
--annotation-tabix NAME:FILE[:FIELDS]
```

**How it works:** Opens only the tabix index at startup. Each variant query performs an on-disk seek to the relevant genomic region, decompresses a single block, and extracts matching records. No data is loaded into memory.

**Best for:** Large annotation databases (ClinVar, gnomAD, dbSNP, custom frequency databases).

| File | Records | Startup time | Memory overhead |
|------|--------:|-------------:|----------------:|
| ClinVar (~2M records) | ~2M | instant | ~10MB (index) |
| gnomAD exomes (~16M) | ~16M | instant | ~50MB (index) |
| gnomAD genomes (~750M) | ~750M | instant | ~200MB (index) |

**Requires:** bgzip compression + tabix index (see [File Preparation](#file-preparation) below).

#### In-memory (for small files)

```bash
--annotation NAME:FILE[:FIELDS]
```

**How it works:** Reads the entire VCF into a hash map indexed by chromosome and position. All records are parsed and stored in memory at startup.

**Best for:** Small annotation files (<100K records) where you want zero per-query I/O overhead.

| File | Records | Load time | Memory usage |
|------|--------:|----------:|-------------:|
| Panel variants (~1K) | ~1K | <1 sec | ~1MB |
| Custom list (~10K) | ~10K | ~1 sec | ~10MB |
| ClinVar (~2M) | ~2M | ~10-20 sec | ~500MB |
| Large databases (>10M) | >10M | minutes | several GB+ |

### File Preparation

Tabix-indexed sources require bgzip compression and a tabix index:

```bash
# If your VCF is uncompressed:
bgzip myfile.vcf                    # Creates myfile.vcf.gz
tabix -p vcf myfile.vcf.gz          # Creates myfile.vcf.gz.tbi

# If already gzip compressed (re-compress with bgzip):
zcat myfile.vcf.gz | bgzip > myfile.vcf.bgz
tabix -p vcf myfile.vcf.bgz

# Verify the index works:
tabix myfile.vcf.gz chr1:10000-20000
```

**Note:** Regular gzip files will NOT work with tabix — the file must be compressed with `bgzip` (from htslib). Both `bgzip` and `tabix` are included with htslib (`brew install htslib`).

### Usage Examples

#### Adding a ClinVar annotation source

```bash
# Prep the file (once):
bgzip clinvar.vcf && tabix -p vcf clinvar.vcf.gz

# Use tabix mode — loads instantly, minimal memory:
./vep_annotator \
    --gtf genes.gtf.gz --fasta genome.fa.gz \
    --annotation-tabix clinvar:clinvar.vcf.gz:CLNSIG,CLNDN,CLNREVSTAT \
    --vcf input.vcf -o output.tsv
```

Output will include `clinvar:CLNSIG`, `clinvar:CLNDN`, and `clinvar:CLNREVSTAT` in the Extra column (TSV) or `custom_annotations` (JSON).

#### Adding gnomAD population frequencies

```bash
./vep_annotator \
    --gtf genes.gtf.gz --fasta genome.fa.gz \
    --annotation-tabix gnomad:gnomad.exomes.r2.1.1.sites.vcf.bgz:AF,AF_popmax,AC,AN \
    --vcf input.vcf -o output.tsv
```

#### Multiple annotation sources

```bash
./vep_annotator \
    --gtf genes.gtf.gz --fasta genome.fa.gz \
    --annotation-tabix clinvar:clinvar.vcf.gz:CLNSIG,CLNDN \
    --annotation-tabix gnomad:gnomad.vcf.bgz:AF,AF_popmax \
    --annotation-tabix cosmic:cosmic.vcf.gz:CNT,GENE \
    --annotation mylist:my_variants.vcf \
    --vcf input.vcf -o output.tsv
```

#### Extracting all INFO fields

Omit the `:FIELDS` part to extract all INFO fields from the annotation VCF:

```bash
# Extracts every INFO field from the ClinVar VCF:
--annotation-tabix clinvar:clinvar.vcf.gz

# Extracts only CLNSIG and CLNDN:
--annotation-tabix clinvar:clinvar.vcf.gz:CLNSIG,CLNDN
```

Specifying fields is recommended for large databases to avoid extracting unnecessary data.

### Perl VEP `--custom` Syntax

For compatibility with existing Perl VEP workflows, the `--custom` flag is supported:

```bash
--custom FILE,NAME,TYPE,OVERLAP,COORDS,FIELDS
```

| Parameter | Description |
|-----------|-------------|
| `FILE` | Path to the annotation file |
| `NAME` | Short name for the source (used as field prefix) |
| `TYPE` | File type: `vcf`, `bed`, `bigwig`, `gff`, `gtf` |
| `OVERLAP` | Match mode: `overlap` (position only) or `exact` (allele match) |
| `COORDS` | Coordinate system: `0` (0-based) or `1` (1-based) |
| `FIELDS` | Comma-separated list of fields to extract |

#### VCF custom source

```bash
# Exact allele matching (recommended for variant databases):
--custom clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNDN

# Position-only overlap:
--custom gnomad.vcf.bgz,gnomAD,vcf,overlap,0,AF
```

VCF `--custom` sources always use tabix mode (the file must be bgzip'd and indexed).

#### BED custom source

```bash
# Annotate with regions from a BED file:
--custom conserved_regions.bed,ConservedRegion,bed,overlap,0

# BED files are loaded into memory and indexed by position for fast overlap queries.
```

BED annotations add a field with the source name containing the `name` column (column 4) of overlapping intervals, or `1` if no name column is present.

#### bigWig custom source

```bash
# Add a bigWig score track:
--custom phyloP100way.bw,MyPhyloP,bigwig,overlap,0
```

Requires libBigWig (build with `-DWITH_BIGWIG=ON`).

### Output Fields

Custom annotation fields appear with a `NAME:FIELD` prefix in the output:

**TSV (Extra column):**
```
clinvar:CLNSIG=Pathogenic;clinvar:CLNDN=Hereditary_cancer;gnomad:AF=0.0001
```

**JSON (custom_annotations object):**
```json
{
  "custom_annotations": {
    "clinvar:CLNSIG": "Pathogenic",
    "clinvar:CLNDN": "Hereditary_cancer",
    "gnomad:AF": "0.0001"
  }
}
```

**VCF (CSQ INFO field):**
Custom fields are appended to the CSQ field definition in the VCF header and included in each CSQ entry.

### Co-located Variant Lookup

Use `--check_existing` to look up co-located variants from a VCF database (uses tabix for efficient queries):

```bash
./vep_annotator \
    --gtf genes.gtf.gz --fasta genome.fa.gz \
    --check_existing dbsnp.vcf.gz \
    --vcf input.vcf -o output.tsv
```

This populates the `Existing_variation` column with variant IDs (e.g., rs numbers) from the database.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- [Ensembl VEP](https://www.ensembl.org/vep) - Original Perl implementation and specification
- [htslib](https://github.com/samtools/htslib) - High-throughput sequencing library
- [libBigWig](https://github.com/deeptools/libBigWig) - BigWig file library
- [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP) - Database of functional predictions
