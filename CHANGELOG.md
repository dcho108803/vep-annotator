# Changelog

All notable changes to the VEP Annotator project are documented here.

## [1.6.1] - 2026-06-23

### Fixed
- **Multi-codon MNV HGVSp** (was a deferred limitation in 1.6.0): a same-length substitution spanning more than one codon reported only the first changed residue (e.g. `p.Val2Lys` when two residues changed). It now emits a protein `delins` over the changed residue range — `p.Val2_Arg3delinsLysPro` — trimming unchanged flanking residues so a single net change still renders as a plain missense. A stop introduced anywhere in the MNV (incl. a non-first codon) now renders as `stop_gained` (`p.Arg3Ter`) instead of a wrong missense/synonymous. Matches Perl VEP and the `AMINO_ACIDS` field. +3 end-to-end tests (1199 → 1202).

### Still deferred
- HGVSg 3'-shift; ≥2-codon in-frame indel HGVSp; deletion spanning an exon/intron boundary (codon math); UTR/intronic insertion duplication detection. These need dedicated Perl-VEP truth vectors.

## [1.6.0] - 2026-05-29

Comprehensive fix sweep addressing the remaining findings from the 2026-05-29 deep
review (60 of 64 applied; 1199 tests passing). Highlights:

### Fixed — correctness
- **SpliceAI gene mismatch (HIGH)**: predictions were selected by VCF ALT-column index into the flat (allele×gene) list, returning the wrong gene's scores when multiple genes overlap. Now matched by allele AND the transcript's gene symbol; INFO key matched at field boundaries.
- **LOFTEE not strand-aware (HIGH)**: donor/acceptor and the NON_CAN_SPLICE dinucleotide check used plus-strand genomic geometry, producing spurious NON_CAN_SPLICE (HC→LC) on nearly every minus-strand splice variant. Now strand-aware. INCOMPLETE_CDS now keys on `cds_start_NF`/`cds_end_NF` (it was dead code).
- **Plugin use-after-unload (HIGH)**: the `PluginLoader` was destroyed (dlclose) while the annotator still held plugin-backed sources. It now outlives all annotators and loads once.
- **SV strand/frame**: strand-aware splice donor/acceptor for SVs; frameshift/inframe now uses CDS-overlap bases (not the full genomic span); symbolic `<INS>` length no longer inflated; DUP/TDUP split from INS; CNV CN=1/unknown handled; symbolic SV start uses POS+1 (VCF spec); BND orientation decoded from bracket character.
- **HGVS parser**: protein `delins`/`dup`/`ext` notation parsed; `inv` (c. and g.); insertion parsing uses `regex_match` (rejects malformed input); 3'UTR `*` positions flagged; end intronic offset stored for ins/dup.
- **Minus-strand non-coding `n.` deletions** now produce correctly-ordered cDNA ranges.
- **CDS 3'-shift off-by-one** for deletions in period>1 tandem repeats.
- **Inframe delins** classified by net length (inframe_insertion/deletion) instead of always protein_altering_variant.
- **--fork** worker exceptions now propagate (non-zero exit) instead of silent partial output.
- **filter_vep**: `match`/`regex` are real regular-expression operators; `regex` reachable from `-f`; `--exclude-intronic` drops compound intronic terms; trailing-`not` negation; numeric/space handling hardened.
- **Transcript pick**: deterministic tiebreaker by transcript_id; per-gene/per-allele output preserves input order; empty `--pick-order` errors.
- **Output**: JSON rejects leading-zero bare numbers; deterministic JSON key order; case-collision guard; `escape_vcf` encodes newlines/control bytes; JSON `close()` flushes; colocated allele uses `-`.
- Consequence terms now emitted most-severe-first.

### Fixed — safety / portability / performance
- UB: `std::toupper`/`isspace` on signed char routed through safe helpers (incl. FASTA load).
- `-march=native` is now opt-in (`-DENABLE_NATIVE=ON`); default builds are portable.
- Uninitialized struct members; `<climits>` include; GFF3 attribute trimming.
- Interval-style transcript lookup (prefix-max-end + binary search) replacing the from-zero linear scan; `build_cds_sequence` returns a reference (no multi-kb copies); JSON per-transcript O(n²) buffer copy removed; gz output buffered (128KB); `sync_with_stdio(false)`; tabix chrom-naming detected once.
- `--output-format` / `--terms` reject unknown values.

### Known limitations (deferred, documented in code)
- HGVSp for multi-codon MNVs and ≥2-codon in-frame indels reports only the first changed residue; deletions spanning an exon/intron boundary use a genomic-span splice. HGVSg is not 3'-shifted (HGVSc is). UTR/intronic insertions are not reported as duplications. These HGVS/codon parity refinements need dedicated test vectors and are deferred to avoid unvalidated changes to clinical notation.

## [1.5.0] - 2026-05-29

### Fixed
- **Minus-strand multi-base CDS position (HIGH)**: a multi-base variant whose genomic-leftmost base lies outside the CDS (e.g. a minus-strand variant spanning the 3'UTR/CDS boundary) was mis-clamped to CDS position 1 and falsely reported as `start_lost`. CDS positions are now computed per affected base (`calculate_cds_position_range`), skipping intron/UTR bases, so the first/last affected CDS positions are correct on both strands and across exon/intron boundaries
- **HGVSg delins dropped bases (HIGH)**: `generate_hgvsg` routed single-base-side delins (e.g. `ATG>C`, `A>GT`) into the deletion/insertion branch and silently dropped the inserted base(s). These now correctly emit `delins`
- **Transcript filters dropping variants (HIGH)**: `--biotype`, `--coding-only`, `--canonical-only`, `--mane-only`, `--no-intergenic`, `--gencode-basic`, `--exclude-predicted`, and include/exclude-consequence filters used without a pick flag ran against only the single most-severe transcript, silently dropping variants whose matching transcript was not the most severe. Both the single-variant and batch paths now annotate all transcripts when any transcript-selecting filter is active (`TranscriptFilterConfig::requires_all_transcripts()`)
- **Gene constraint never emitted**: `--pli` / `--loeuf` / `--constraint` loaded data but it was never queried; `pLI` and `LOEUF` are now looked up per gene (by Ensembl ID then symbol) and emitted in the output

### Tests
- 12 new regression tests (1174 → 1186 total), including end-to-end minus-strand 3'UTR/CDS boundary annotation

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
