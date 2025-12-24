# VEP Annotator

A high-performance C++ implementation of the Variant Effect Predictor (VEP) for annotating genetic variants. This is a pure local implementation that requires no external API calls - all annotation is performed using local data files.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://isocpp.org/std/the-standard)

## Features

### Core Annotation
- Gene and transcript mapping with spatial indexing
- Consequence prediction using Sequence Ontology (SO) terms
- Impact classification (HIGH, MODERATE, LOW, MODIFIER)
- HGVS notation (coding and protein)
- CDS and protein position calculation
- Codon and amino acid change detection

### Pathogenicity Predictions (via dbNSFP)
- **35+ pathogenicity scores**: SIFT, PolyPhen-2, CADD, REVEL, AlphaMissense, MetaSVM, MetaLR, VEST4, PROVEAN, FATHMM, MutationTaster, MutationAssessor, DANN, Eigen, M-CAP, MPC, PrimateAI, BayesDel, ClinPred, and more
- Conservation scores: PhyloP, PhastCons, GERP++, SiPhy
- Population frequencies: gnomAD, ExAC, 1000 Genomes
- Clinical data: ClinVar annotations

### Splice Predictions
- **SpliceAI**: Deep learning splice predictions (delta scores for acceptor/donor gain/loss)
- **MaxEntScan**: Position weight matrix-based splice site scoring
- **dbscSNV**: Splice site variant predictions (ada_score, rf_score)

### Conservation Scores
- **PhyloP**: Phylogenetic p-values from 100-way vertebrate alignment
- **PhastCons**: Probability of negative selection
- **GERP++**: Genomic Evolutionary Rate Profiling

### Regulatory Annotations
- Ensembl Regulatory Build features
- Promoter overlap detection
- Enhancer overlap detection
- Transcription factor binding sites (TFBS)
- Open chromatin regions
- CTCF binding sites

### Protein Domains
- **Pfam**: Protein family annotations
- **InterPro**: Integrated protein domain database

### Loss-of-Function (LoF) Annotations
- **LOFTEE**: High-confidence (HC) / Low-confidence (LC) classification
- **NMD**: Nonsense-mediated decay prediction (50bp rule)
- **LoFtool**: Gene-level constraint scores

### Custom Annotations
- Add any VCF as an annotation source
- ClinVar clinical significance
- gnomAD population frequencies
- Custom in-house annotations
- Support for tabix-indexed files (memory-efficient for large files)

### Plugin System
- Extend with custom annotation sources
- Dynamic loading of shared libraries
- Simple plugin API

### Performance
- Multi-threaded annotation
- Efficient memory usage with on-disk queries
- Tabix support for large VCF files
- Support for gzipped input files

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

# Optional: Install system-wide
sudo make install
```

### Download Data Files

```bash
# Use the helper script to download required data
./scripts/download_data.sh ./data
```

### Basic Usage

```bash
# Annotate a single variant
./vep_annotator --gtf genes.gtf.gz --fasta genome.fa.gz -v 7:140753336:A:T

# Annotate a VCF file
./vep_annotator --gtf genes.gtf.gz --fasta genome.fa.gz --vcf input.vcf.gz -o output.tsv
```

### Full Annotation Pipeline

```bash
./vep_annotator \
    --gtf Homo_sapiens.GRCh38.110.gtf.gz \
    --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    --dbnsfp dbNSFP4.4a.txt.gz \
    --dbnsfp-fields "SIFT_score,Polyphen2_HDIV_score,CADD_phred,REVEL_score,AlphaMissense_score" \
    --spliceai spliceai_scores.vcf.gz \
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

### Pathogenicity Predictions

| Option | Description |
|--------|-------------|
| `--dbnsfp FILE` | dbNSFP database (tabix-indexed .txt.gz) |
| `--dbnsfp-fields FIELDS` | Fields to extract (default: essential) |

**Field presets:**
- `essential`: SIFT, PolyPhen2, CADD, REVEL, AlphaMissense
- `pathogenicity`: All pathogenicity scores
- `conservation`: Conservation scores only
- `all`: All available fields

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

*Note: Requires libBigWig*

### Regulatory Annotations

| Option | Description |
|--------|-------------|
| `--regulatory FILE` | Ensembl Regulatory Build GFF3 file |

### Protein Domains

| Option | Description |
|--------|-------------|
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
| `--annotation NAME:VCF[:FIELDS]` | Add VCF annotation source (in-memory) |
| `--annotation-tabix NAME:VCF[:FIELDS]` | Add tabix-indexed VCF source (on-disk) |

### Plugins

| Option | Description |
|--------|-------------|
| `--plugin PATH[:CONFIG]` | Load plugin from shared library |
| `--plugin-dir DIR` | Directory to search for plugins |

### Output Options

| Option | Description |
|--------|-------------|
| `-o, --output FILE` | Output file path |
| `--all-transcripts` | Output all transcript annotations |

### Other Options

| Option | Description |
|--------|-------------|
| `-h, --help` | Show help message |
| `--debug` | Enable debug logging |

## Output Format

### TSV Output Columns

| Column | Description |
|--------|-------------|
| CHROM | Chromosome |
| POS | Position (1-based) |
| REF | Reference allele |
| ALT | Alternate allele |
| GENE | Gene symbol |
| TRANSCRIPT | Transcript ID |
| CONSEQUENCE | VEP consequence term(s) |
| IMPACT | Impact level |
| CDS_POS | Position in CDS |
| PROTEIN_POS | Position in protein |
| AMINO_ACIDS | Amino acid change (ref/alt) |
| CODONS | Codon change (ref/alt) |
| HGVSc | HGVS coding notation |
| HGVSp | HGVS protein notation |
| *source:field* | Custom annotation columns |

### Example Output

```
CHROM   POS         REF ALT GENE  TRANSCRIPT       CONSEQUENCE      IMPACT   ...
7       140753336   A   T   BRAF  ENST00000288602  missense_variant MODERATE ...
```

## Consequence Types

Consequences are reported using Sequence Ontology terms:

| Impact | Consequences |
|--------|--------------|
| **HIGH** | transcript_ablation, splice_acceptor_variant, splice_donor_variant, stop_gained, frameshift_variant, stop_lost, start_lost |
| **MODERATE** | inframe_insertion, inframe_deletion, missense_variant, protein_altering_variant |
| **LOW** | splice_region_variant, incomplete_terminal_codon_variant, start_retained_variant, stop_retained_variant, synonymous_variant |
| **MODIFIER** | coding_sequence_variant, 5_prime_UTR_variant, 3_prime_UTR_variant, non_coding_transcript_exon_variant, intron_variant, upstream_gene_variant, downstream_gene_variant, intergenic_variant |

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

### Download Script

```bash
./scripts/download_data.sh ./data
```

This downloads:
- Ensembl GTF annotations
- Reference genome FASTA
- Conservation scores (PhyloP, PhastCons)
- Regulatory Build
- ClinVar

*Note: dbNSFP, SpliceAI, and gnomAD require manual download due to licensing.*

## Building

### Dependencies

**Required:**
- C++17 compiler (GCC 7+, Clang 5+, MSVC 2017+)
- CMake 3.16+
- zlib

**Optional:**
- **htslib**: Tabix support for large VCF files
- **libBigWig**: Conservation scores (bigWig format)

### Build Commands

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `-DWITH_HTSLIB=ON/OFF` | ON | Enable tabix support |
| `-DWITH_BIGWIG=ON/OFF` | ON | Enable bigWig support |
| `-DCMAKE_BUILD_TYPE=Release/Debug` | Release | Build type |

### Installing Dependencies

**macOS:**
```bash
brew install htslib
# libBigWig: build from source
```

**Ubuntu/Debian:**
```bash
sudo apt-get install libhts-dev zlib1g-dev
```

**From Source (htslib):**
```bash
git clone https://github.com/samtools/htslib.git
cd htslib && make && sudo make install
```

## Performance

### Benchmarks

| Operation | Speed |
|-----------|-------|
| Single variant annotation | < 1ms |
| VCF (1,000 variants) | ~2 seconds |
| VCF (100,000 variants) | ~3 minutes |
| VCF with all annotations | ~5x slower |

### Memory Usage

| Component | Memory |
|-----------|--------|
| GTF annotations | ~500MB - 2GB |
| Reference FASTA | ~3GB |
| dbNSFP (tabix) | ~100MB |
| Conservation (bigWig) | ~200MB |
| ClinVar (in-memory) | ~500MB |
| gnomAD (tabix) | ~100MB |
| gnomAD (in-memory) | 30-50GB |

**Recommendation:** Use `--annotation-tabix` for large VCF files to minimize memory usage.

## Using as a Library

```cpp
#include "vep_annotator.hpp"
#include "annotation_sources.hpp"

int main() {
    // Initialize annotator
    vep::VEPAnnotator annotator("genes.gtf.gz", "genome.fa.gz");

    // Add annotation sources
    auto dbnsfp = vep::create_dbnsfp_source("dbNSFP.txt.gz", "essential");
    annotator.add_source(dbnsfp);

    auto phylop = vep::create_phylop_source("phyloP100way.bw");
    annotator.add_source(phylop);

    // Initialize all sources
    annotator.initialize_sources();

    // Annotate a variant
    auto results = annotator.annotate("7", 140753336, "A", "T");

    for (const auto& ann : results) {
        std::cout << "Gene: " << ann.gene_symbol << std::endl;
        std::cout << "Consequence: " << ann.get_consequence_string() << std::endl;
        std::cout << "Impact: " << vep::impact_to_string(ann.impact) << std::endl;

        // Access custom annotations
        for (const auto& [key, value] : ann.custom_annotations) {
            std::cout << key << ": " << value << std::endl;
        }
    }

    return 0;
}
```

## Plugin Development

Create custom annotation plugins by implementing the `VEPPlugin` interface:

```cpp
#include "plugin.hpp"

class MyPlugin : public vep::VEPPlugin {
public:
    vep::PluginInfo get_info() const override {
        return {"myplugin", "1.0.0", "My custom plugin", "Author", "MIT"};
    }

    bool initialize(const vep::PluginConfig& config) override {
        // Initialize with config options
        return true;
    }

    std::vector<std::shared_ptr<vep::AnnotationSource>> create_sources() override {
        return { std::make_shared<MySource>() };
    }
};

VEP_PLUGIN_EXPORT(MyPlugin)
```

Compile as a shared library:
```bash
g++ -std=c++17 -shared -fPIC -I/path/to/include myplugin.cpp -o myplugin.so
```

Load in VEP:
```bash
./vep_annotator --plugin ./myplugin.so:key=value --gtf genes.gtf --fasta genome.fa -v 7:140753336:A:T
```

## Troubleshooting

### Common Issues

**"Cannot open GTF file"**
- Check file path and permissions
- Ensure file is not corrupted

**"Tabix index not found"**
- Create index: `tabix -p vcf file.vcf.gz`
- Ensure `.tbi` file is in same directory

**"libBigWig not found"**
- Build with `-DWITH_BIGWIG=OFF` or install libBigWig

**Out of memory with gnomAD**
- Use `--annotation-tabix` instead of `--annotation`

### Debug Mode

```bash
./vep_annotator --debug --gtf genes.gtf --fasta genome.fa -v 7:140753336:A:T
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Citation

If you use this software in your research, please cite:

```
VEP Annotator: A high-performance C++ variant effect predictor
https://github.com/dcho108803/vep-annotator
```

## Acknowledgments

- [Ensembl VEP](https://www.ensembl.org/vep) - Original Perl implementation
- [htslib](https://github.com/samtools/htslib) - High-throughput sequencing library
- [libBigWig](https://github.com/deeptools/libBigWig) - BigWig file library
- [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP) - Database of functional predictions

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Support

- **Issues**: [GitHub Issues](https://github.com/dcho108803/vep-annotator/issues)
- **Discussions**: [GitHub Discussions](https://github.com/dcho108803/vep-annotator/discussions)
