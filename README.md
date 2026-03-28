# VEP Annotator

A high-performance C++ implementation of Ensembl's [Variant Effect Predictor (VEP)](https://www.ensembl.org/vep). Achieves ~99.9% feature parity with the Perl VEP while running **75-115x faster**. All annotation is performed locally using standard data files — no API calls required.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://isocpp.org/std/the-standard)
[![Tests](https://img.shields.io/badge/tests-1014%20passing-brightgreen.svg)]()

## Quick Start

```bash
# Install dependencies (macOS)
brew install cmake htslib

# Build
git clone https://github.com/dcho108803/vep-annotator.git
cd vep-annotator && mkdir build && cd build
cmake .. && make -j4

# Run tests (1014 tests)
./vep_tests

# Annotate a single variant
./vep_annotator --gtf genes.gtf.gz --fasta genome.fa.gz -v 7:140753336:A:T

# Annotate a VCF file
./vep_annotator --gtf genes.gtf.gz --fasta genome.fa.gz --vcf input.vcf.gz -o output.tsv
```

<details>
<summary><strong>Ubuntu/Debian dependencies</strong></summary>

```bash
sudo apt-get install cmake libhts-dev zlib1g-dev
```
</details>

### Build Dependencies

| Dependency | Required | Purpose |
|------------|----------|---------|
| C++17 compiler (GCC 7+, Clang 5+) | Yes | Core build |
| CMake 3.16+ | Yes | Build system |
| zlib | Yes | Compressed file support |
| htslib | Optional | Tabix-indexed files (`--dbnsfp`, `--spliceai`, `--annotation-tabix`) |
| libBigWig | Optional | Conservation bigWig files (`--phylop`, `--phastcons`, `--gerp`) |

### Build Options

| CMake Option | Default | Description |
|--------------|---------|-------------|
| `-DWITH_HTSLIB=ON/OFF` | ON | Enable tabix support |
| `-DWITH_BIGWIG=ON/OFF` | ON | Enable bigWig support |
| `-DCMAKE_BUILD_TYPE` | Release | `Release` or `Debug` |

## Features

**Core annotation** — All 39 SO consequence types with impact classification, HGVS notation (HGVSc/HGVSp/HGVSg), CDS/protein position calculation, codon changes, 3' shifting, and structural variant annotation (DEL, DUP, INS, INV, CNV, BND).

**Transcript selection** — `--pick`, `--per_gene`, `--pick_allele`, `--flag_pick`, `--most_severe` with MANE Select/Plus Clinical, canonical, TSL, APPRIS, and CCDS ranking.

**Output formats** — TSV (default), JSON (Perl VEP-compatible schema), and VCF (CSQ INFO field). Use `--everything` for comprehensive output.

**Annotation sources:**

| Source | Flag | Description |
|--------|------|-------------|
| dbNSFP | `--dbnsfp FILE` | 35+ pathogenicity scores (SIFT, PolyPhen-2, CADD, REVEL, AlphaMissense, etc.), conservation, population frequencies. Use `--dbnsfp-fields` to select fields or presets (`essential`/`pathogenicity`/`conservation`/`all`). |
| SpliceAI | `--spliceai FILE` | Deep learning splice predictions |
| MaxEntScan | `--maxentscan` | Splice site scoring (algorithmic, no data file) |
| dbscSNV | `--dbscsnv FILE` | Splice variant predictions (ada_score, rf_score) |
| PhyloP | `--phylop FILE` | Conservation (100-way vertebrate alignment) |
| PhastCons | `--phastcons FILE` | Conservation (probability of negative selection) |
| GERP++ | `--gerp FILE` | Genomic Evolutionary Rate Profiling |
| Regulatory | `--regulatory FILE` | Ensembl Regulatory Build (promoters, enhancers, TFBS). Use `--cell_type LIST` to filter. |
| Pfam | `--pfam FILE` | Protein family domain annotations |
| InterPro | `--interpro FILE` | Integrated protein domain database |
| LOFTEE | `--loftee` | Loss-of-function HC/LC classification |
| NMD | `--nmd` | Nonsense-mediated decay prediction |
| LoFtool | `--loftool FILE` | Gene-level constraint scores |
| Custom VCF | `--annotation NAME:FILE[:FIELDS]` | In-memory VCF annotation |
| Custom VCF (tabix) | `--annotation-tabix NAME:FILE[:FIELDS]` | On-disk tabix VCF queries |
| Custom (Perl-compat) | `--custom FILE,NAME,TYPE,...` | VCF, BED, bigWig, GFF/GTF |

**Plugin system** — Extend with custom C++ annotation sources via shared libraries (`dlopen`).

## Usage Examples

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
    --loftee --nmd \
    --annotation clinvar:clinvar.vcf.gz:CLNSIG,CLNDN \
    --annotation-tabix gnomad:gnomad.exomes.vcf.bgz:AF,AF_popmax \
    --vcf input.vcf.gz \
    -o annotated.tsv
```

### Output Formats

```bash
# TSV (default)
./vep_annotator --gtf genes.gtf.gz --fasta genome.fa.gz --vcf input.vcf.gz -o output.tsv

# JSON
./vep_annotator --gtf genes.gtf.gz --fasta genome.fa.gz --vcf input.vcf.gz --json -o output.json

# VCF with CSQ annotation
./vep_annotator --gtf genes.gtf.gz --fasta genome.fa.gz --vcf input.vcf.gz --vcf-output -o output.vcf

# Everything mode
./vep_annotator --gtf genes.gtf.gz --fasta genome.fa.gz --vcf input.vcf.gz --everything -o output.tsv
```

### Post-Processing with filter_vep

A standalone `filter_vep` binary is built alongside `vep_annotator`:

```bash
# HIGH impact variants
./filter_vep -i output.tsv -o filtered.tsv --impact HIGH

# Rare missense variants in a gene list
./filter_vep -i output.tsv -o filtered.tsv \
    --consequence missense_variant --max-af 0.01 --gene-list genes.txt

# Custom filter expressions with AND/OR logic
./filter_vep -i output.tsv -o filtered.tsv \
    -f 'IMPACT in HIGH,MODERATE' -f 'AF < 0.01' --canonical-only

# Count matching records
./filter_vep -i output.tsv --count -f 'CADD_phred > 20'
```

Supports consequence/impact/gene/biotype filters, numeric thresholds (`--min-af`, `--max-af`, `--min-cadd`, `--min-revel`), boolean flags (`--coding-only`, `--canonical-only`, `--mane-only`, `--pick`), and arbitrary filter expressions. Run `./filter_vep --help` for full usage.

## Performance

Benchmarked on 100,000 chr22 variants (Release build, Apple Silicon M-series):

| Mode | Wall time | Throughput | vs Perl VEP |
|------|----------:|-----------:|------------:|
| Single thread, TSV | 4.0s | **25,000 var/sec** | ~80x faster |
| Single thread, JSON | 4.0s | **25,000 var/sec** | ~80x faster |
| Single thread, --everything | 4.4s | **22,700 var/sec** | ~75x faster |
| --fork 4, --everything | 2.9s | **34,500 var/sec** | ~115x faster |

Use `--fork N` for parallel annotation with N threads. Input is pre-buffered and distributed across threads using atomic work-stealing.

### Memory Usage

| Component | Memory |
|-----------|--------|
| GTF annotations | ~500MB - 2GB |
| Reference FASTA | ~3GB |
| dbNSFP (tabix) | ~100MB |
| Conservation (bigWig) | ~200MB |
| ClinVar (in-memory) | ~500MB |
| gnomAD (tabix) | ~100MB |

Use `--annotation-tabix` instead of `--annotation` for large VCF files to minimize memory.

## Testing

1014 GoogleTest unit tests across 11 test suites:

| Test Suite | Tests | Coverage |
|------------|------:|----------|
| Filter VEP | 206 | Operators, expressions, conditions, pipeline, edge cases |
| Output Writers | 197 | TSV/JSON/VCF formatting, escaping, stats, position formatting |
| CLI Utilities | 144 | Variant parsing, allele trimming, config parsing, format detection |
| Transcript Filter | 121 | Pick modes, filtering, ranking, APPRIS/TSL, complex scenarios |
| HGVS | 90 | Parsing, HGVSg generation, SPDI, RefSeq mapping, edge cases |
| Consequences | 84 | SO terms, impact levels, ranking, variant class, display terms |
| Structural Variants | 60 | SV types, BND parsing, consequences, overlap, properties |
| Annotation Sources | 40 | LOFTEE, NMD, MaxEntScan, dbNSFP, source manager, domains |
| Codon Table | 35 | Translation, MT codons, completeness, case handling, edge cases |
| Exon/Intron Numbers | 25 | Position calculation, formatting |
| SpliceAI | 12 | Score parsing, cutoffs, thread safety |
| **Total** | **1014** | |

```bash
cd build && ./vep_tests
```

## Data Files

### Required

| File | Source |
|------|--------|
| GTF gene annotations | [Ensembl FTP](https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/) |
| Reference genome FASTA | [Ensembl FTP](https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/) |

### Optional

| File | Size | Source |
|------|------|--------|
| dbNSFP (pathogenicity scores) | ~35GB | [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP) |
| SpliceAI (splice predictions) | ~20GB | [Illumina](https://basespace.illumina.com/s/otSPW8hnhaZR) |
| PhyloP (conservation) | ~10GB | [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/) |
| PhastCons (conservation) | ~10GB | [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/) |
| Regulatory Build | ~500MB | [Ensembl FTP](https://ftp.ensembl.org/pub/release-110/regulation/homo_sapiens/) |
| ClinVar | ~100MB | [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) |
| gnomAD | ~30GB+ | [gnomAD](https://gnomad.broadinstitute.org/downloads) |

## Command Line Reference

<details>
<summary><strong>Input & Output</strong></summary>

| Option | Description |
|--------|-------------|
| `--gtf FILE` | GTF/GFF annotation file (required) |
| `--fasta FILE` | Reference genome FASTA (required) |
| `-v, --variant CHR:POS:REF:ALT` | Single variant to annotate |
| `--vcf FILE` | VCF file for batch annotation |
| `-` or `STDIN` | Read VCF from standard input |
| `-o, --output FILE` | Output file (default: stdout) |
| `--json` | JSON output format |
| `--vcf-output` | VCF output format (CSQ INFO field) |
| `--tab` | TSV output format (default) |
| `--compress_output bgzip` | Compress output with bgzip |
| `--no_headers` | Suppress header lines |
| `--fields FIELDS` | Custom field order (comma-separated) |

</details>

<details>
<summary><strong>Transcript Selection</strong></summary>

| Option | Description |
|--------|-------------|
| `--pick` | Pick one annotation per variant (most severe) |
| `--pick_allele` | Pick one annotation per allele |
| `--per_gene` | Pick one annotation per gene |
| `--pick_allele_gene` | Pick one per allele per gene |
| `--flag_pick` | Flag picked annotation but output all |
| `--pick_order LIST` | Custom pick priority order |
| `--most_severe` | Only report most severe consequence |

</details>

<details>
<summary><strong>Display Options</strong></summary>

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
| `--sift [DISPLAY]` | SIFT output (p = prediction, s = score, b = both; default: p) |
| `--polyphen [DISPLAY]` | PolyPhen output (p/s/b; default: p) |
| `--af` | Show allele frequency |
| `--max_af` | Show maximum allele frequency |
| `--variant_class` | Show variant class (SNV/insertion/deletion) |
| `--allele_number` | Show allele number |
| `--total_length` | Show transcript/CDS lengths |
| `--nearest` | Show nearest gene for intergenic variants |

</details>

<details>
<summary><strong>Variant Processing</strong></summary>

| Option | Description |
|--------|-------------|
| `--minimal` | Allele-level minimization for multi-allele VCFs |
| `--check_ref` | Verify reference allele against FASTA |
| `--check_existing FILE` | Co-located variant lookup VCF |
| `--keep_csq` | Preserve existing VCF CSQ annotations |
| `--no_intergenic` | Skip intergenic variants |
| `--coding_only` | Only annotate coding transcripts |
| `--exclude_predicted` | Exclude predicted (XM_/XR_) transcripts |
| `--chr LIST` | Only annotate these chromosomes |
| `--distance N` | Upstream/downstream distance (default: 5000) |
| `--assembly NAME` | Assembly name for output (default: GRCh38) |
| `--synonyms FILE` | Chromosome synonym mapping file |
| `--config FILE` | Load options from config file (key=value) |
| `--stats_file FILE` | Write run statistics to file |
| `--fork N` | Use N parallel annotation threads |
| `--buffer_size N` | Variants to buffer (default: 5000) |
| `--quiet` | Suppress progress messages |
| `--filter_common` | Remove common variants |
| `--freq_pop POP` | Population for frequency filter |
| `--freq_freq FREQ` | Frequency threshold |
| `--freq_gt_lt gt/lt` | Greater than or less than threshold |

</details>

<details>
<summary><strong>Annotation Sources</strong></summary>

| Option | Description |
|--------|-------------|
| `--dbnsfp FILE` | dbNSFP database (tabix-indexed .txt.gz) |
| `--dbnsfp-fields FIELDS` | Fields to extract (or preset: `essential`/`pathogenicity`/`conservation`/`all`) |
| `--spliceai FILE` | SpliceAI VCF file (tabix-indexed) |
| `--maxentscan` | Enable MaxEntScan splice site scoring (algorithmic) |
| `--dbscsnv FILE` | dbscSNV file (tabix-indexed) |
| `--phylop FILE` | PhyloP bigWig file |
| `--phastcons FILE` | PhastCons bigWig file |
| `--gerp FILE` | GERP++ bigWig file |
| `--regulatory FILE` | Ensembl Regulatory Build GFF3 file |
| `--cell_type LIST` | Cell type filter for regulatory (comma-separated) |
| `--pfam FILE` | Pfam domain annotations TSV |
| `--interpro FILE` | InterPro domain annotations TSV |
| `--loftee` | Enable LOFTEE-style LoF classification |
| `--nmd` | Enable NMD prediction |
| `--loftool FILE` | LoFtool gene constraint scores |
| `--annotation NAME:VCF[:FIELDS]` | Add in-memory VCF annotation source |
| `--annotation-tabix NAME:VCF[:FIELDS]` | Add tabix-indexed VCF source |
| `--custom FILE,NAME,TYPE,OVERLAP,COORDS,FIELDS` | Perl VEP-compatible custom source |

</details>

## Output Format

### TSV Columns

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
| Extra | Additional fields (key=value pairs separated by `;`) |

### Consequence Types

Ranked by severity using [Sequence Ontology](http://www.sequenceontology.org/) terms:

| Impact | Consequences |
|--------|-------------|
| **HIGH** | transcript_ablation, splice_acceptor_variant, splice_donor_variant, stop_gained, frameshift_variant, stop_lost, start_lost, transcript_amplification |
| **MODERATE** | inframe_insertion, inframe_deletion, missense_variant, protein_altering_variant |
| **LOW** | splice_donor_5th_base_variant, splice_region_variant, splice_donor_region_variant, splice_polypyrimidine_tract_variant, incomplete_terminal_codon_variant, start_retained_variant, stop_retained_variant, synonymous_variant |
| **MODIFIER** | coding_sequence_variant, mature_miRNA_variant, 5_prime_UTR_variant, 3_prime_UTR_variant, non_coding_transcript_exon_variant, intron_variant, NMD_transcript_variant, non_coding_transcript_variant, upstream_gene_variant, downstream_gene_variant, TFBS_ablation, TFBS_amplification, TF_binding_site_variant, regulatory_region_ablation, regulatory_region_amplification, feature_elongation, regulatory_region_variant, feature_truncation, intergenic_variant |

## Custom Annotation Files

Two loading modes for VCF annotation sources:

**Tabix-indexed** (recommended for large files) — on-disk queries, instant startup:
```bash
--annotation-tabix clinvar:clinvar.vcf.gz:CLNSIG,CLNDN,CLNREVSTAT
```

**In-memory** (for small files < 100K records) — zero per-query I/O:
```bash
--annotation mylist:my_variants.vcf:AF
```

Custom fields appear as `NAME:FIELD` in output (e.g., `clinvar:CLNSIG=Pathogenic` in TSV Extra column).

<details>
<summary><strong>File preparation, examples, and --custom syntax</strong></summary>

### Tabix File Preparation

Tabix sources require bgzip compression + tabix index:

```bash
bgzip myfile.vcf                    # Creates myfile.vcf.gz
tabix -p vcf myfile.vcf.gz          # Creates myfile.vcf.gz.tbi

# Re-compress gzip to bgzip:
zcat myfile.vcf.gz | bgzip > myfile.vcf.bgz
tabix -p vcf myfile.vcf.bgz
```

Regular gzip will NOT work — use `bgzip` from htslib.

### Scaling

| Source | Records | Mode | Startup | Memory |
|--------|--------:|------|---------|--------|
| Panel variants | ~1K | in-memory | <1s | ~1MB |
| ClinVar | ~2M | tabix | instant | ~10MB |
| gnomAD exomes | ~16M | tabix | instant | ~50MB |
| gnomAD genomes | ~750M | tabix | instant | ~200MB |

### Multiple Sources

```bash
./vep_annotator \
    --gtf genes.gtf.gz --fasta genome.fa.gz \
    --annotation-tabix clinvar:clinvar.vcf.gz:CLNSIG,CLNDN \
    --annotation-tabix gnomad:gnomad.vcf.bgz:AF,AF_popmax \
    --annotation-tabix cosmic:cosmic.vcf.gz:CNT,GENE \
    --annotation mylist:my_variants.vcf \
    --vcf input.vcf -o output.tsv
```

Omit `:FIELDS` to extract all INFO fields (not recommended for large databases).

### Perl VEP `--custom` Syntax

```bash
--custom FILE,NAME,TYPE,OVERLAP,COORDS,FIELDS
```

| Parameter | Values |
|-----------|--------|
| TYPE | `vcf`, `bed`, `bigwig`, `gff`, `gtf` |
| OVERLAP | `overlap` (position only) or `exact` (allele match) |
| COORDS | `0` (0-based) or `1` (1-based) |

```bash
# VCF with exact allele matching:
--custom clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNDN

# BED regions:
--custom conserved_regions.bed,ConservedRegion,bed,overlap,0

# bigWig scores (requires libBigWig):
--custom phyloP100way.bw,MyPhyloP,bigwig,overlap,0
```

### Co-located Variant Lookup

```bash
./vep_annotator --gtf genes.gtf.gz --fasta genome.fa.gz \
    --check_existing dbsnp.vcf.gz --vcf input.vcf -o output.tsv
```

Populates `Existing_variation` with variant IDs (e.g., rs numbers) via tabix queries.

### Output Format

**TSV Extra:** `clinvar:CLNSIG=Pathogenic;gnomad:AF=0.0001`

**JSON:**
```json
{ "custom_annotations": { "clinvar:CLNSIG": "Pathogenic", "gnomad:AF": "0.0001" } }
```

**VCF:** Custom fields appended to the CSQ INFO field definition.

</details>

## Using as a Library

```cpp
#include "vep_annotator.hpp"
#include "annotation_sources.hpp"

int main() {
    vep::VEPAnnotator annotator("genes.gtf.gz", "genome.fa.gz");

    auto dbnsfp = vep::create_dbnsfp_source("dbNSFP.txt.gz", "essential");
    annotator.add_source(dbnsfp);
    annotator.initialize_sources();

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

## Differences from Perl VEP

This C++ implementation uses **local files exclusively** (GTF + FASTA) rather than Perl VEP's cache/database/API model. This means faster startup, simpler deployment, and explicit data sources — but no built-in population frequencies or cache-based lookups.

<details>
<summary><strong>Missing features</strong></summary>

| Feature | Notes |
|---------|-------|
| Ensembl cache/database | Use `--gtf` + `--fasta` instead |
| REST API queries | Runs purely locally |
| Built-in population frequencies | Add via `--annotation-tabix gnomad:file.vcf.gz:AF` |
| HTML statistics report | Text stats via `--stats-file` |
| Individual/sample genotypes | Parsed but not functional |
| Motif feature consequences | Not implemented |
| PubMed / gene-phenotype / variant synonyms | Require variant database |
| miRNA / LRG / GA4GH VRS / BAM correction | Not implemented |
| Transcript filter expressions | Use individual filter flags |
| Multi-species support | Single species only |

</details>

<details>
<summary><strong>Behavioral differences</strong></summary>

| Feature | Perl VEP | C++ VEP |
|---------|----------|---------|
| `--everything` + regulatory | Auto-includes from cache | Requires explicit `--regulatory FILE` |
| `--check_existing` | Uses cache | Requires explicit VCF path |
| `--sift` / `--polyphen` | Scores from cache | Scores from `--dbnsfp`; flags control display only |
| `--domains` | From cache | Requires `--pfam` / `--interpro` files |
| `--fork N` | Perl `fork()` processes | C++ `std::thread` + work-stealing |
| `--minimal` VCF output | Preserves original alleles | Writes minimized alleles |
| Plugin system | Perl plugin ecosystem | C++ shared libraries (`.so`/`.dylib`) |
| `filter_vep` regex | Full Perl regex | Substring matching (`string::find`) |
| HGVS multi-base duplications | Full `dup` notation | May use `ins` instead of `dup` |

</details>

<details>
<summary><strong>Known edge cases</strong></summary>

- **MNV multi-codon extraction**: Multi-nucleotide variants spanning multiple codons may have limited codon display (consequence determination is unaffected)
- **UTR-intronic HGVSc**: Variants in UTR introns may have imprecise HGVS notation (very rare)
- **Complex indel codon boundaries**: Complex indels spanning codon boundaries may show approximate codon changes (consequence type is correct)
- **`-v` mode dash alleles**: The single-variant parser's `-` to `:` replacement can misparse dash-containing deletion alleles (use VCF input instead)

</details>

<details>
<summary><strong>Perl VEP compatibility flags (accepted, no-op)</strong></summary>

`--cache`, `--offline`, `--database`, `--merged`, `--refseq`, `--dir-cache`, `--cache-version`, `--host`, `--port`, `--user`, `--password`, `--registry-file`, `--db-version`, `--safe`, `--tmpdir`, `--dir`, `--force`, `--dont-skip`, `--lookup-ref`, `--failed`, `--no-whole-genome`, `--verbose`, `--show-cache-info`, `--shift-length`, `--shift-hgvs`, `--no-check-alleles`, `--exclude-null-alleles`, `--check-svs`, `--custom-multi-allelic`, `--max-sv-size`, `--no-check-variants-order`, `--overlap-cutoff`, `--af-exac`, `--hgvsg-use-accession`

These allow existing Perl VEP command lines to run without modification.

</details>

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for a detailed history of all releases and changes.

## License

MIT License — see [LICENSE](LICENSE) for details.

## Acknowledgments

- [Ensembl VEP](https://www.ensembl.org/vep) — Original Perl implementation and specification
- [htslib](https://github.com/samtools/htslib) — High-throughput sequencing library
- [libBigWig](https://github.com/deeptools/libBigWig) — BigWig file library
- [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP) — Database of functional predictions
