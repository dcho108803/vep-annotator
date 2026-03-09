# Changelog

All notable changes to the VEP Annotator project are documented here.

## [1.4.0] - 2026-03-09

### Fixed
- **--minimal --fork cache key mismatch**: Indel cache lookups used "-" but cache keys used empty strings, causing all indels to miss the pre-annotation cache and fall back to single-threaded annotation under `--fork N`
- **JSON invalid numbers**: SIFT/PolyPhen scores, TSL, and custom annotation values written as raw unvalidated JSON numbers could produce invalid JSON when values were empty, "NA", or malformed. Added `is_valid_json_number()` validator
- **--config misleading message**: Previously printed "Loaded N arguments" without applying them. Now warns that config file support is not yet implemented
- **Thread-unsafe logging**: Replaced `std::localtime` (shared static buffer) with `localtime_r` for safe multi-threaded logging under `--fork`

## [1.3.0] - 2026-02-21

### Added
- Comprehensive "Differences from Perl VEP" section in README
- Detailed "Custom Annotation Files" documentation with loading modes, file preparation, usage examples, and output format details

### Changed
- Updated performance benchmarks with Release build numbers (25,000 var/sec single-thread, 34,500 var/sec with --fork 4)

## [1.2.0] - 2026-02-20

### Added
- **--fork N multi-threading**: Real parallel annotation using `std::thread` with atomic work-stealing. Pre-buffers all input, deduplicates variants, annotates across N threads, replays from cache
- **--chr LIST**: Filter variants by chromosome
- **--exclude-predicted**: Exclude predicted (XM_/XR_) transcripts
- **--synonyms FILE**: Chromosome synonym mapping
- **--stats-file FILE**: Write run statistics to file
- **--config FILE**: Config file parsing (currently read-only)
- **548 new unit tests** (total: 640): filter_vep (180), output writers (153), transcript filter (104), CLI utilities (76)

### Fixed
- **CDS phase for cds_start_NF transcripts**: `calculate_cds_position()` and `build_cds_sequence()` now account for first CDS phase in incomplete transcripts
- **64KB line buffer**: `gz_read_line()` now reads in chunks and concatenates for arbitrary-length VCF lines
- **parse_variant dash handling**: Fixed dash-to-colon replacement corrupting deletion alleles

## [1.1.2] - 2026-02-20

### Fixed
- **APPRIS format**: Corrected APPRIS annotation output format
- **Consequence ranking**: Fixed most_severe sorting to use ConsequenceType enum (was only sorting by 4-level Impact)
- **Minus-strand codons**: Fixed codon extraction for minus-strand indels
- **BED chr-prefix**: Fixed chromosome name matching for BED custom annotation sources
- **--minimal mode**: Fixed suffix trimming (`size() > 1` to `!empty()`) and display allele computation
- **Filter routing**: Fixed consequence filter case sensitivity

## [1.1.1] - 2026-02-15

### Fixed
- 12 bugs from comprehensive code audit including:
  - CDS-space 3' shift formula for multi-base insertions
  - Genomic right_normalize for insertions
  - Intergenic distance calculation for multi-base variants
  - MaxEntScan diff guard for invalid reference scores
  - VCF header filter substring match precision
  - Multi-allele VCF flush ordering (duplicate output lines)

## [1.1.0] - 2026-02-14

### Fixed
- 8 code review issues: crash safety (HGVS parser stoi, exon_intron_numbers underflow), logic bugs (flag_pick gene_id matching, regulatory cell_type filter), performance improvements

## [1.0.2] - 2026-02-13

### Changed
- Upgraded SpliceAI source implementation for Perl VEP plugin parity

## [1.0.1] - 2026-02-12

### Changed
- Performance optimizations and code structure improvements
- Updated README with accurate performance data and full CLI reference

## [1.0.0] - 2026-01-08

### Added
- Comprehensive VEP feature parity with Perl implementation (~99.9%)
- All 38 Sequence Ontology consequence types with correct ranking
- HGVSc, HGVSp, HGVSg notation with 3' shifting
- Three output formats: TSV, JSON, VCF (CSQ INFO field)
- 8 native annotation sources: dbNSFP, SpliceAI, MaxEntScan, dbscSNV, PhyloP/PhastCons/GERP, regulatory, Pfam/InterPro, LOFTEE/NMD
- Custom annotation sources: VCF (in-memory and tabix), BED, bigWig
- Transcript filtering: --pick, --per_gene, --pick_allele, --flag_pick, --most_severe
- Structural variant annotation (DEL, DUP, INS, INV, CNV, BND)
- --everything, --minimal, --check_existing, --keep_csq, --filter_common
- Perl VEP CLI compatibility (unsupported flags accepted as no-ops)
- filter_vep equivalent for post-filtering annotations
- C++ plugin interface for shared library extensions

## [0.1.0] - 2025-12-23

### Added
- Initial commit: core VEP annotation engine
- GTF/FASTA parsing, consequence determination, basic output
