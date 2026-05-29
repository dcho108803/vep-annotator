# VEP_ensemble_C — Deep Code Review (2026-05-29)

Multi-agent review: 14 finders across all subsystems + adversarial verification. **72 findings confirmed** (4 rejected by skeptics) plus re-check of 13 documented bugs.

## Known-bug re-check

- **FIXED** — 1. [HIGH] --fork OOM: main.cpp buffered entire VCF in buffered_lines vector, no size limit.
  - /Users/davidcho/Projects/VEP_ensemble_C/src/main.cpp:2211-2218, 2253-2278, 2419-2426
  - Comment: "Pre-buffer input and parallel pre-annotation for --fork N (chunked to bound memory: full-file buffering OOMs on large WGS VCFs)." const size_t FORK_CHUNK_SIZE = 10000;  // lines per chunk. refill_fork_chunk() clears and refills at most FORK_CHUNK_SIZE lines: while (buffered_lines.size() < FORK_CHUNK_SIZE && !fork_input_eof) {...}. Consumer loop calls refill_fork_chunk() when buffered_lin
- **FIXED** — 2. [HIGH] --fork cache race: multiple threads write to annotation_cache unordered_map without synchronization.
  - /Users/davidcho/Projects/VEP_ensemble_C/src/main.cpp:2360-2401, 2736-2744
  - Worker threads write ONLY to distinct indices of a pre-sized vector, never the map: std::vector<std::vector<vep::VariantAnnotation>> results(queries.size()); ... while ((i = work_idx.fetch_add(1)) < queries.size()) { ... results[i] = ann.annotate(...); }. The map is populated only on the main thread AFTER join: for (auto& w : workers) w.join(); ... for (size_t i = 0; i < queries.size(); i++) { ann
- **FIXED** — 3. [HIGH] Multi-allelic VCF ALT: stored single vcf_alt per annotation; VCF output writes wrong ALT for non-first alleles.
  - /Users/davidcho/Projects/VEP_ensemble_C/src/main.cpp:2550, 2891; /Users/davidcho/Projects/VEP_ensemble_C/include/output_writer.hpp:1246-1256, 1270-1274
  - alt is the FULL original ALT column (alt = vcf_fields[4]), and a.vcf_alt = alt stores the complete comma-separated column. Comment at main.cpp: '// Preserve original VCF REF/ALT for multi-allele output'. flush_variant() writes one line per position with the full column: std::string vcf_alt = !first.vcf_alt.empty() ? first.vcf_alt : ...; line << vcf_ref << '\t' << vcf_alt << '\t' ...; and emits one
- **PRESENT** — 4. [MED] CDS 3' shift off-by-one: deletion shifting checks cached_cds[shift_pos + indel_len - 1] (last deleted base) instead of base after deletion.
  - /Users/davidcho/Projects/VEP_ensemble_C/src/vep_annotator.cpp:2078-2092
  - int shift_pos = cds_pos; ... bool is_cds_insertion = (ref.length() < alt.length()); for (int i = 0; i < 1000; ++i) { int check_idx = is_cds_insertion ? shift_pos : (shift_pos + indel_len - 1); ... }. For deletions check_idx is still shift_pos + indel_len - 1 (the last deleted base, since cds_pos is 1-based the deleted region is 0-based [cds_pos, cds_pos+indel_len-1]); shifting right should compare
- **PRESENT** — 5. [MED] Minus-strand multi-base CDS position: subtracts ref.length()-1 but overestimates for variants spanning CDS/intron boundaries.
  - /Users/davidcho/Projects/VEP_ensemble_C/src/vep_annotator.cpp:2917-2920 (and identically at 1999-2002)
  - if (transcript.strand == '-' && ref.length() > 1) { cds_pos = cds_pos - static_cast<int>(ref.length()) + 1; if (cds_pos < 1) cds_pos = 1; }. Still subtracts the full ref.length() naively, assuming all ref bases are contiguous in CDS space; for variants spanning a CDS/intron boundary the intronic ref bases should not count, so the shift overestimates. git -L (commit d0a4bd5) only broadened the guar
- **FIXED** — 6. [MED] Incomplete terminal codon early return: returns INCOMPLETE_TERMINAL_CODON skipping co-occurring HIGH consequences (stop_gained, frameshift).
  - /Users/davidcho/Projects/VEP_ensemble_C/src/vep_annotator.cpp:2925-2937, 2985-2989
  - consequences.push_back(ConsequenceType::INCOMPLETE_TERMINAL_CODON_VARIANT); // ... Indels fall through so co-occurring frameshift / start_lost / stop_lost are also attached. if (ref.length() == alt.length()) { return consequences; }. The early return now only fires for length-preserving variants; indels fall through to the frameshift/start_lost/stop_lost checks (e.g. if (length_diff != 0 && length
- **FIXED** — 7. [MED] Frameshift Ter-distance phase: scan_start uses unadjusted protein_position for cds_start_NF transcripts.
  - /Users/davidcho/Projects/VEP_ensemble_C/src/vep_annotator.cpp:2151 (scan), 3144-3153 (calculate_cds_position phase adjust), 3185-3193 (build_cds_sequence phase trim)
  - scan_start = (ann.protein_position - 1) * 3 where protein_position derives from cds_pos via calculate_cds_position(), which now applies the cds_start_NF phase: 'if (in_cds && transcript.cds_start_NF && !transcript.cds_regions.empty()) { int first_phase = ...; if (first_phase > 0 && first_phase <= 2) { cds_pos -= first_phase; if (cds_pos <= 0) return 0; } }'. build_cds_sequence() trims the SAME fir
- **FIXED** — 8. [MED] filter_vep REGEX uses string::find() not actual regex.
  - /Users/davidcho/Projects/VEP_ensemble_C/include/filter_vep.hpp:23, 246-252
  - #include <regex>; and: } else if (cond.op == FilterOperator::REGEX) { try { std::regex re(cond.value); result = std::regex_search(value, re); } catch (const std::regex_error&) { // Fallback to substring match if regex is invalid; result = (value.find(cond.value) != std::string::npos); } }. Now uses std::regex_search; substring find() is only a fallback for invalid patterns.
- **FIXED** — 9. [MED] --config silently ignored: parsed config file but warned 'not yet implemented'.
  - /Users/davidcho/Projects/VEP_ensemble_C/src/main.cpp:705-720 (load_config_tokens), 904-934 (apply)
  - Config tokens are loaded and PREPENDED to argv so they are actually parsed/applied: std::vector<std::string> cfg_tokens = load_config_tokens(pre_config_path, std::cerr); ... for (auto& t : cfg_tokens) merged_arg_storage.push_back(std::move(t)); for (int i = 1; i < argc; ++i) merged_arg_storage.emplace_back(argv[i]); ... argc = ...; argv = merged_argv_vec.data(); std::cerr << "Loaded config: " << p
- **FIXED** — 10. [LOW] Dead code: norm == "chrM" can never match after chr-prefix stripping.
  - /Users/davidcho/Projects/VEP_ensemble_C/src/vep_annotator.cpp:242-247, 261-268, 289-296
  - Prefix is stripped first, then compared against 'MT'/'M' (no 'chrM' literal): 'std::string norm = chromosome; if (norm.size() > 3 && norm.compare(0, 3, "chr") == 0) norm = norm.substr(3); if (norm == "MT" || norm == "M") {...}'. grep for '"chrM"' across src/vep_annotator.cpp returns no matches (exit 1), so the dead comparison was removed.
- **PARTIAL** — 11. [LOW] std::toupper on signed char without unsigned cast (technically UB). Multiple sites.
  - /Users/davidcho/Projects/VEP_ensemble_C/src/vep_annotator.cpp:2085, 210/232/238/252/258/281/327/1044/1090; /Users/davidcho/Projects/VEP_ensemble_C/src/sources/maxentscan_source.cpp:342
  - Many sites were fixed with static_cast<unsigned char> (e.g. hgvs_parser.hpp, structural_variant.hpp, main.cpp:506). But several still pass raw signed char: vep_annotator.cpp:2085 'if (std::toupper(cached_cds[check_idx]) == std::toupper(indel_bases[0]))'; maxentscan_source.cpp:342 'switch (std::toupper(base)) {' (base is char); and numerous 'std::transform(upper.begin(), upper.end(), upper.begin(),
- **PRESENT** — 12. [LOW] Intronic position equidistant tiebreaker: arbitrary <= when equidistant from both exon boundaries.
  - /Users/davidcho/Projects/VEP_ensemble_C/src/vep_annotator.cpp:3435-3454 (calculate_intronic_hgvsc_position)
  - int dist_to_upstream_exon_end = genomic_pos - transcript.exons[i].end; int dist_to_downstream_exon_start = transcript.exons[i + 1].start - genomic_pos; if (transcript.strand == '+') { if (dist_to_upstream_exon_end <= dist_to_downstream_exon_start) {...+notation...} else {...-notation...} } else { if (dist_to_downstream_exon_start <= dist_to_upstream_exon_end) {...} }. The bare <= still arbitrarily
- **PRESENT** — 13. [LOW] std::move into ostringstream::str(): C++17 str() takes const ref so move is silently a copy.
  - /Users/davidcho/Projects/VEP_ensemble_C/include/output_writer.hpp:1032; /Users/davidcho/Projects/VEP_ensemble_C/CMakeLists.txt:5
  - json.str(std::move(s)); json.seekp(0, std::ios_base::end); — and CMakeLists.txt sets CMAKE_CXX_STANDARD 17. In C++17 ostringstream::str() has only the const-ref setter (the str(string&&) rvalue overload is C++20), so the move binds to const std::string& and copies. Empirically verified under -std=c++17: after oss.str(std::move(s)) the source string s still contains all 66 chars ('INTACT(copy - mov

## Confirmed findings


### HIGH

**[core-consequence/correctness] src/vep_annotator.cpp:2914-2962** — Minus-strand indel spanning the 3'UTR/CDS boundary is mis-called start_lost (HIGH-impact wrong call)
- What: For a minus-strand transcript, cds_start is the genomic minimum (the 3' end of the transcript, near the stop codon) and cds_end is the genomic maximum (5' end, near the start codon). When a length-changing variant (deletion/MNV with multi-base ref) overlaps the CDS but its lowest genomic coordinate pos is below cds_start (i.e. in the genomically-left 3'UTR), calculate_cds_position(pos) correctly returns 0 (line 3110). The minus-strand block then computes cds_pos = 0 - ref_len + 1 (negative) and force-clamps it to 1 (line 2919). With cds_pos==1 and length_diff<0, the start-codon block fires affects_start=true and pushes START_LOST, then returns. The variant is actually at the 3' end of the protein (near the stop), so it should be stop_lost / inframe_deletion / frameshift / 3_prime_UTR_variant — not start_lost. start_lost is a HIGH-impact term, so this produces a clinically wrong, over-severe consequence and short-circuits all correct downstream logic.
- Evidence: if (transcript.strand == '-' && ref.length() > 1) {
    cds_pos = cds_pos - static_cast<int>(ref.length()) + 1;
    if (cds_pos < 1) cds_pos = 1;
}
...
if (length_diff != 0 && cds_pos >= 1 && cds_pos <= 3) {
    bool affects_start = false;
    if (length_diff < 0) { affects_start = true; }
    ...
    if (affects_start) {
        consequences.push_back(ConsequenceType::START_LOST);
        return consequences;
    }
}
- Fix: Do not blindly clamp `cds_pos` to 1 when `calculate_cds_position(pos)` returned 0 due to `pos < cds_start`. Distinguish the two reasons cds_pos can be non-positive after the minus-strand `- ref_len + 1` adjustment:

1) Real start-codon overlap: `pos` is inside the CDS and the variant extends into the genomically-high 5'UTR (var_end > cds_end). Here calculate_cds_position(pos) is a small positive number; clamping to 1 and entering the START_LOST branch is correct.

2) 3'UTR-left overlap (the bug): `pos < cds_start` so calculate_cds_position(pos) == 0. Here the first in-CDS base is actually near the STOP codon. Compute the true first affected CDS base from the clamped overlapping coordinate, e.g. map `min(var_end, cds_end)` and `max(pos, cds_start)` through calculate_cds_position and take the smaller resulting CDS position, instead of `cds_pos - ref_len + 1` with a clamp to 1.

Concretely: only enter the START_LOST branch when the genomic span [pos, var_end] actually overlaps the start codon. On minus strand the start codon is the three highest genomic CDS bases (the cds_end end), i.e. require `var_end >= cds_end - 2` (start codon overlap) before allowing affects_start. When `calculate_cds_position(pos)` returns 0 because `pos < cds_start`, derive cds_pos from `calculate_cds_position(max(pos, cds_start))` (which yields the high CDS position near the stop) so the variant correctly falls into the stop_lost / frameshift / inframe_deletion / 3_prime_UTR logic. Add a regression test for a minus-strand multi-base deletion spanning the genomically-left 3'UTR/CDS boundary asserting it is NOT start_lost.
- Verifier: Verified against current code in /Users/davidcho/Projects/VEP_ensemble_C/src/vep_annotator.cpp.

Entry condition (line 2914): `if (in_exon && pos <= transcript.cds_end && var_end >= transcript.cds_start)`. Note `var_end = pos + ref.length() - 1` (line 2710), so `pos` is the genomically LOWEST coordinate and `var_end` the highest.

Minus-strand adjustment (lines 2915-2920):
```
int cds_pos = calculate_cds_position(pos, transcript);
if (transcript.strand == '-' && ref.length() > 1) {
    cds_pos =

**[core-coord/correctness] src/vep_annotator.cpp:1996-2002 and 2915-2920** — Minus-strand multi-base variant: CDS position derived by subtracting ref length is wrong when the genomic-leftmost base is outside the CDS or the variant spans an intron, producing spurious START_LOST and wrong cds_position/protein_position
- What: For minus-strand transcripts the code computes the first (lowest) affected CDS coordinate as `calculate_cds_position(pos) - ref.length() + 1` where `pos` is the genomic-leftmost base. Two failure modes: (1) On a minus-strand gene the 3'UTR lies at genomic positions BELOW cds_start. A deletion whose genomic-leftmost base sits in the 3'UTR but which extends into the CDS (e.g. spanning the stop-codon/3'UTR boundary) makes `calculate_cds_position(pos)` return 0 (pos < cds_start, see line 3110), then the subtraction yields a negative number that is clamped to 1. The variant is therefore treated as if it hit CDS position 1 (the start codon). At lines 2946-2962 a deletion with `cds_pos>=1 && cds_pos<=3` is unconditionally flagged START_LOST and the function returns, so a variant that should be STOP_LOST/3'UTR is reported as START_LOST, and cds_position/cds_end/protein_position (lines 2031-2035) are all reported as ~1 instead of near the 3' end. (2) Even when `pos` is inside the CDS, if the deletion spans an exon/intron boundary, `ref.length()` counts intronic bases that do not exist in CDS space, so `cds_pos - ref.length() + 1` over-subtracts and lands on a wrong (too-small) CDS coordinate. The correct first CDS position for minus strand is obtained by mapping the genomic-RIGHTMOST in-CDS base, e.g. `calculate_cds_position(std::min(var_end, transcript.cds_end))`, not by arithmetic on ref.length().
- Evidence: if (transcript.strand == '-' && ref.length() > 1) {
            cds_pos = cds_pos - static_cast<int>(ref.length()) + 1;
            if (cds_pos < 1) cds_pos = 1;
        }
        if (cds_pos > 0) {
            annotate_coding_region(chrom, pos, ref, alt, transcript, cached_cds, cds_pos, ann);
- Fix: In both annotate_transcript (vep_annotator.cpp:1996-2002) and determine_consequences (2915-2920), replace the ref.length()-subtraction-and-clamp with a direct mapping of the lowest CDS coordinate. For minus strand, the lowest affected CDS position corresponds to the genomic-rightmost in-CDS base:

  int var_end = pos + static_cast<int>(ref.length()) - 1; // already available in determine_consequences
  int cds_pos;
  if (transcript.strand == '-' && ref.length() > 1) {
      int rightmost = std::min(var_end, transcript.cds_end);
      rightmost = std::max(rightmost, transcript.cds_start);
      cds_pos = calculate_cds_position(rightmost, transcript);
  } else {
      cds_pos = calculate_cds_position(pos, transcript);
  }
  // Do NOT clamp a non-positive result to 1. If cds_pos <= 0, fall through to
  // UTR/intron handling (annotate_noncds_hgvsc / non-coding paths) rather than
  // entering the coding/start-codon block.

This uses calculate_cds_position's intron-aware region walk (fixing Mode 2's over-subtraction) and never fabricates a cds_pos==1 start-codon hit when the leftmost base is in the 3'UTR (fixing Mode 1). Apply identically in both functions. Add an end-to-end test: a minus-strand transcript with a multi-base deletion spanning the stop-codon/3'UTR boundary should yield STOP_LOST or 3'_prime_UTR_variant (not START_LOST) with cds_position near the 3' end, and a deletion spanning an exon/intron boundary should produce a CDS position consistent with calculate_cds_position(var_end).
- Verifier: Both cited locations match the quoted code exactly. In annotate_transcript (vep_annotator.cpp:1996-2002) and in determine_consequences (2915-2920) the identical pattern appears:

  int cds_pos = calculate_cds_position(pos, transcript);
  if (transcript.strand == '-' && ref.length() > 1) {
      cds_pos = cds_pos - static_cast<int>(ref.length()) + 1;
      if (cds_pos < 1) cds_pos = 1;
  }

Premise verified: calculate_cds_position (line 3108-3155) returns 0 when genomic_pos < cds_start || > cds_e

**[gene-exon-plugin/memory-safety] src/main.cpp:1939-1956, 1970** — Plugin shared library is dlclose()'d while annotator still holds plugin-backed sources (use-after-unload)
- What: PluginLoader 'plugin_loader' is a local variable inside the 'setup_annotation_sources' lambda. The lambda receives 'annotator' by reference from the outer scope (created at main.cpp:2096, used for the whole annotation run). Plugin sources are obtained via get_all_sources() and handed to annotator.add_source() (line 1955), where they are stored as shared_ptr<AnnotationSource> with no tie to the dlopen handle. When the lambda returns at line 1970, ~PluginLoader runs and calls unload_plugin() for every plugin, which invokes destroy_func and then dlclose(p.handle), unmapping the plugin's code and vtables. The annotator then continues annotating variants using these now-dangling sources, calling virtual methods whose code/vtable has been unmapped, and eventually destroying the shared_ptr (running the plugin-defined object dtor/allocator deleter that no longer exists). This is a classic dlclose use-after-unload and will crash or corrupt memory whenever a plugin is actually loaded. The DlHandleGuard/RAII work inside load_plugin is correct, but the loader object itself outlives its handles' consumers in the wrong order.
- Evidence: src/plugin.hpp:99 '~PluginLoader();' and src/plugin_loader.cpp:27-33 destructor calls unload_plugin which does dlclose(p.handle) (line 198). In src/main.cpp the loader is local: 'vep::PluginLoader plugin_loader;' (1939) ... 'auto plugin_sources = plugin_loader.get_all_sources();' (1953) ... 'annotator.add_source(plugin_sources[i]);' (1955), inside lambda 'setup_annotation_sources = [&](vep::VEPAnnotator& annotator) -> std::vector<std::string> { ... return cols; };' that returns at line 1970, destroying plugin_loader while 'annotator' (outer scope, line 2096) keeps the sources.
- Fix: Make the PluginLoader outlive the annotator and all plugin-backed sources. Concretely, hoist `vep::PluginLoader` out of the `setup_annotation_sources` lambda into the surrounding scope (the same scope as, or outer to, the `annotator` at line 2096 and the per-thread `thread_annotators`), declared BEFORE the annotator so it is destroyed AFTER it. The lambda should reference this outer loader (captured by `&`) instead of declaring a local. Note the lambda is invoked once per annotator (main at 2106 and per-thread at 2195), so the loader must accumulate/serve sources for all of them without being torn down between calls; consider loading plugins once up front and only calling get_all_sources()/add_source per annotator. Alternatively (option b): attach each dlopen handle's lifetime to the returned sources via an aliasing/custom-deleter shared_ptr so dlclose happens only after the last source shared_ptr is released. Do not destroy PluginLoader (and thus dlclose) until after the annotator and every plugin-created source has been destroyed.
- Verifier: Verified against the actual current code; the claim's structure and line citations are accurate.

1) The PluginLoader is a lambda-local. At src/main.cpp:1939 inside the `setup_annotation_sources` lambda (defined `auto setup_annotation_sources = [&](vep::VEPAnnotator& annotator) -> std::vector<std::string> {` at line 1777): `vep::PluginLoader plugin_loader;`. Plugins are loaded (1943, 1948), then `auto plugin_sources = plugin_loader.get_all_sources();` (1953) and `annotator.add_source(plugin_sour

**[hgvs-parser/correctness] include/hgvs_parser.hpp:794-806** — generate_hgvsg silently drops inserted bases for single-base-alt delins, emitting a plain deletion
- What: The deletion branch is selected whenever `alt.empty() || (ref.size() > alt.size() && alt.size() <= 1)`. This captures genuine delins variants where the alt is a single base that is NOT the VCF anchor of ref. For example ref=AAT alt=C at pos 100 produces `g.100_102del`, and ref=AT alt=C produces `g.100_101del` — both silently delete the inserted base C. The anchor-strip guard `ref[0] == alt[0]` is only used to adjust del_start, never to fall through to the delins branch, so the substituted base is lost entirely. This produces clinically wrong HGVSg (a deletion instead of a substitution/delins) for a common variant class.
- Evidence: } else if (alt.empty() || (ref.size() > alt.size() && alt.size() <= 1)) {
        // Deletion
        int del_start = pos;
        int del_end = pos + static_cast<int>(ref.size()) - 1;
        // Strip VCF anchor base if present
        if (!alt.empty() && alt.size() == 1 && ref.size() > 1 && ref[0] == alt[0]) {
            del_start = pos + 1;
        }
- Fix: Restrict the deletion branch to true deletions. Change the condition at line 794 so a single-base alt only enters the deletion branch when it equals the ref anchor. For example:

} else if (alt.empty() ||
           (ref.size() > 1 && alt.size() == 1 && ref[0] == alt[0])) {
    // Deletion (alt absent, or single-base alt equal to the VCF anchor)
    ...
}

This makes the anchored-deletion case (e.g. ref=AAT alt=A -> g.101_102del) and the pure-deletion case (alt empty) take the deletion branch, while a single-base alt that differs from the anchor (e.g. ref=AT alt=C, ref=AAT alt=C) falls through to the existing delins branch (lines 837-851), which already produces g.100_101delinsC / g.100_102delinsC correctly (its anchor-strip at line 843 requires alt.size() > 1, so single-base alts pass through unstripped). Add unit tests in test/test_hgvs.cpp covering generate_hgvsg for single-base-alt delins to lock in the behavior.
- Verifier: Confirmed against the actual code at include/hgvs_parser.hpp:794-806 (current line numbers match the claim exactly):

    } else if (alt.empty() || (ref.size() > alt.size() && alt.size() <= 1)) {
        // Deletion
        int del_start = pos;
        int del_end = pos + static_cast<int>(ref.size()) - 1;
        // Strip VCF anchor base if present
        if (!alt.empty() && alt.size() == 1 && ref.size() > 1 && ref[0] == alt[0]) {
            del_start = pos + 1;
        }
        ...

Tracing 

**[main-cli/correctness] src/main.cpp:3270-3307** — Single-variant mode: --coding-only / --canonical-only / --mane-only / --biotype / --no-intergenic / --exclude-predicted / --gencode-basic / --include-consequence / --exclude-consequence filter only the single most-severe transcript, silently dropping variants that have matching transcripts
- What: In single-variant mode the decision to fetch ALL transcripts (annotator.annotate) is gated only on all_transcripts and the pick/flag-pick flags. For every other selection flag the code falls through to annotate_most_severe(), which returns exactly ONE annotation (the most severe transcript). The transcript filter at lines 3295-3307 then applies coding_only / canonical_only / mane_only / gencode_basic / biotypes / no_intergenic / exclude_predicted / include_consequences / exclude_consequences via passes_basic_filter(), which DROPS a non-matching annotation. So e.g. `--coding-only` on a variant whose most-severe transcript is non-coding produces ZERO output even though protein_coding transcripts overlap the variant. This is a parity break with Perl VEP (which evaluates all transcripts then filters) and causes silent loss of clinically relevant annotations.
- Evidence: } else if (all_transcripts ||
           filter_config.pick || filter_config.pick_allele ||
           filter_config.pick_allele_gene || filter_config.per_gene ||
           filter_config.flag_pick || filter_config.flag_pick_allele ||
           filter_config.flag_pick_allele_gene) {
    annotations = annotator.annotate(chrom, pos, ref, alt);
} else {
    auto ann = annotator.annotate_most_severe(chrom, pos, ref, alt);
    annotations.push_back(ann);
}
...
if (filter_config.pick || ... filter_config.canonical_only || filter_config.mane_only ||
    filter_config.coding_only || filter_config.gencode_basic || filter_config.exclude_predicted ||
    !filter_config.biotypes.empty() || filter_config.no_intergenic || filter_config.check_frequency ||
    !filter_config.include_consequences.empty() || !filter_config.exclude_consequences.empty()) {
    annotations = transcript_filter.filter(annotations);
}
- Fix: Make the single-variant "fetch all transcripts" condition at src/main.cpp:3270-3274 identical to the filter-application condition at 3295-3307 (and to the existing `need_all`/`need_all_ann` conditions at lines 2222-2233 and 2723-2734). Concretely, add filter_config.most_severe, coding_only, canonical_only, mane_only, gencode_basic, no_intergenic, exclude_predicted, !biotypes.empty(), !include_consequences.empty(), and !exclude_consequences.empty() to the else-if at line 3270 so annotator.annotate() (all transcripts) is called whenever any selection/filter flag is active. Best practice: extract a single shared helper (e.g. bool need_all_transcripts(const TranscriptFilterConfig&, bool all_transcripts)) and reuse it at all four sites (2222, 2723, 3270, and as the basis for 3295) so they cannot drift out of sync again.
- Verifier: The cited code matches the claim exactly. At src/main.cpp:3270-3279 the single-variant (-v/--variant) path gates the "fetch all transcripts" branch only on `all_transcripts` plus the pick/per_gene/flag_pick flags:

  } else if (all_transcripts ||
             filter_config.pick || filter_config.pick_allele ||
             filter_config.pick_allele_gene || filter_config.per_gene ||
             filter_config.flag_pick || filter_config.flag_pick_allele ||
             filter_config.flag_pick_allel

**[main-cli/correctness] src/main.cpp:2723-2734, 2784, 2788** — Batch mode: --biotype LIST (without a pick flag) filters only the single most-severe transcript, dropping variants whose matching-biotype transcripts are not the most severe
- What: In batch mode the need_all decision (lines 2723-2734, mirrored by need_all_ann at 2222-2233) lists coding_only/canonical_only/mane_only/gencode_basic/no_intergenic/exclude_predicted/most_severe/include_consequences/exclude_consequences but OMITS !filter_config.biotypes.empty(). The filter-application condition at line 2784 DOES include !filter_config.biotypes.empty(). Consequence: running with only `--biotype protein_coding` (no --pick/--coding-only) causes annotate_most_severe() to return a single transcript; the biotype filter (passes_basic_filter line 379-382) then drops it if its biotype is not in the list, even though other matching-biotype transcripts overlap. Variants are silently lost / mis-annotated, a parity break vs Perl VEP.
- Evidence: bool need_all = all_transcripts || ... filter_config.coding_only || filter_config.canonical_only ||
    filter_config.mane_only || filter_config.gencode_basic || filter_config.no_intergenic ||
    filter_config.exclude_predicted || !filter_config.include_consequences.empty() ||
    !filter_config.exclude_consequences.empty();   // <-- no biotypes / check_frequency
...
if (... !filter_config.biotypes.empty() || filter_config.no_intergenic || filter_config.check_frequency ...) {
    annotations = transcript_filter.filter(annotations);
}
- Fix: Add !filter_config.biotypes.empty() to all three need-all decision points so biotype filtering always operates over the complete transcript set: src/main.cpp:2723 (need_all, batch path), src/main.cpp:2222 (need_all_ann, fork pre-annotation), and the branch at src/main.cpp:3270-3274 (which is missing many flags, not just biotypes). For check_frequency, add it for output-parity consistency (lower priority). The robust fix is to factor the need-all predicate into a single shared helper (e.g. bool needs_all_transcripts(const TranscriptFilterConfig&)) that exactly mirrors the filter-application condition (main.cpp:2784-2787, 3302-3305), and call it from the batch path, the fork pre-annotation path, and the third single-variant path, eliminating the drift. Add a regression test: a locus where the most-severe consequence falls on a non-protein_coding transcript while an overlapping protein_coding transcript exists, run with only --biotype protein_coding, and assert the protein_coding transcript is annotated.
- Verifier: Confirmed against current code. The two cited need-all predicates omit biotypes (and check_frequency), while the filter-application conditions include them, producing an inconsistency that corrupts biotype filtering.

src/main.cpp:2723-2734 (batch path):
  bool need_all = all_transcripts || filter_config.pick || ... || filter_config.most_severe ||
      filter_config.coding_only || filter_config.canonical_only || filter_config.mane_only ||
      filter_config.gencode_basic || filter_config.no_in

**[perf/performance] src/vep_annotator.cpp:1501-1529** — get_transcripts_in_region linear-scans the whole chromosome transcript vector from index 0 (no binary search)
- What: chrom_index is a vector sorted by transcript start (build_index sorts at line 1184-1186), but get_transcripts_in_region iterates from begin() of the whole per-chromosome vector for every variant. The only optimization is an early break when tr_start > end. For a variant near the end of a large chromosome (e.g. chr1 has tens of thousands of transcripts), this scans every transcript whose start is <= the variant position on EVERY variant lookup. This is the single most-called function in the per-variant hot path. The sibling gene index get_nearby_genes (line 1581) already uses std::lower_bound, proving the codebase knows the right approach; transcripts were left as a from-zero scan. Cost: O(transcripts-before-pos) per variant instead of O(log n + hits). Compounded because annotate() calls it a SECOND time per variant (the wider_transcripts scan at line 1756) with an even wider window.
- Evidence: auto it = pimpl_->chrom_index.find(normalized);
    if (it == pimpl_->chrom_index.end()) {
        return results;
    }

    for (const auto& [tr_start, tr_end, tid] : it->second) {
        // Check overlap
        if (tr_start <= end && tr_end >= start) {
            ...
        }
        // Since sorted by start, can break early
        if (tr_start > end) break;
    }
- Fix: Replace the from-zero linear scan in get_transcripts_in_region. Because transcripts can start before the query window yet still overlap it (long introns/large genes), a plain std::lower_bound on start (as get_nearby_genes uses) is not sufficient on its own. Two sound options: (1) Use the existing IntervalTree<size_t> from file_parsers.hpp — build one per chromosome in build_index() (storing transcript indices/ids), then call query(start, end), bounding cost to O(log n + hits). (2) Keep the sorted-by-start vector but maintain a parallel max-end prefix (running maximum of tr_end up to each index); std::lower_bound to the first start <= end, then scan backward only while the max-end prefix indicates a possible overlap. Option (1) reuses existing, tested code and is preferred. Separately, in annotate() the second wider-window call (lines 1756-1757) doubles the cost; once the index is logarithmic this is acceptable, but consider skipping the wider scan entirely when upstream_distance_ and downstream_distance_ are both 0 / flanks are not requested.
- Verifier: The cited code matches the claim exactly. At src/vep_annotator.cpp:1516, `for (const auto& [tr_start, tr_end, tid] : it->second)` iterates from the beginning of the entire per-chromosome vector, with the only optimization being `if (tr_start > end) break;` at line 1525. The vector IS sorted by start: build_index() (lines 1176-1186) emplaces (transcript.start, transcript.end, tid) tuples and runs std::sort (line 1185), which orders primarily by start. This means the loop visits every transcript w

**[sources-a/parity-gap] src/sources/spliceai_source.cpp:88, 146, 221-266** — SpliceAI ignores transcript/gene; allele index used to select wrong gene's prediction
- What: The Perl VEP SpliceAI plugin emits one prediction per (allele, gene) and reports the prediction for the gene of the transcript being annotated. Here transcript is explicitly discarded ((void)transcript) and the INFO list is indexed purely by alt_index (the index of the matched ALT in the VCF ALT column). SpliceAI INFO is comma-separated by (allele x gene) entries, e.g. 'SpliceAI=A|GENE1|...,A|GENE2|...,C|GENE1|...'. find_alt_index returns the ALT-column allele index (0 for A, 1 for C), but allele_values[alt_index] indexes the flat (allele,gene) entry list. When more than one gene overlaps the position, allele_values[alt_index] selects the wrong entry (e.g. alt_index=1 picks the second gene of allele A rather than allele C), so the reported gene symbol and delta scores can belong to a different allele/gene than the variant. This yields wrong SpliceAI_pred output relative to Perl VEP.
- Evidence: void annotate(... const Transcript* transcript, ...) {
    ensure_initialized();
    (void)transcript;
...
int alt_index = find_alt_index(alt_it->second, alt);  // index within VCF ALT column
...
std::string value = allele_values[alt_index];  // indexes flat (allele,gene) list
- Fix: In parse_spliceai_info, do not index by the VCF ALT-column index. Instead iterate every comma-separated SpliceAI entry, split each on '|', and (1) match field[0] (ALLELE) against the variant alt, and (2) among allele-matching entries, select the one whose field[1] (SYMBOL) equals the transcript's gene symbol (Transcript::gene_name; fall back to gene_id if needed). Pass the transcript through to parse_spliceai_info (stop discarding it via (void)transcript) so its gene is available. If no entry matches both allele and gene, follow Perl VEP behavior (emit nothing / or all matching-allele entries when no transcript context). Update query_reader/query to thread the transcript/gene through. Add a test using a multi-gene + multi-allele INFO string (e.g. SpliceAI=A|G1|...,A|G2|...,C|G1|...,C|G2|...) asserting the correct allele+gene entry is chosen.
- Verifier: Verified against current code in /Users/davidcho/Projects/VEP_ensemble_C/src/sources/spliceai_source.cpp. All three cited facts hold:

1. Line 88: `(void)transcript;` inside `annotate(... const Transcript* transcript, ...)` — the transcript context (which carries `gene_name`/`gene_id`, confirmed available on the Transcript struct in include/vep_annotator.hpp:210-212) is explicitly discarded. The transcript IS passed in: src/file_parsers.cpp:919 calls `source->annotate(chrom, pos, ref, alt, trans

**[sources-b/correctness] src/sources/lof_source.cpp:84-217** — LOFTEE splice donor/acceptor detection is not strand-aware: NON_CAN_SPLICE reads the wrong dinucleotide on minus strand
- What: The LOFTEE source classifies splice sites purely by genomic geometry: a position 2bp before an exon's genomic start is always treated as an acceptor, and 2bp after an exon's genomic end is always treated as a donor (lines 84-97). This is only correct for plus-strand transcripts. The main consequence engine (vep_annotator.cpp ~2726-2789) is correctly strand-aware: for minus-strand transcripts, exon.start-2..-1 is the DONOR side and exon.end+1..+2 is the ACCEPTOR side. Because the LOFTEE source ignores strand, on a minus-strand transcript: (a) is_splice_donor/is_splice_acceptor are swapped relative to the true consequence; and (b) the NON_CAN_SPLICE dinucleotide check reads the wrong genomic side. The donor branch reads reference at exon.end+1..exon.end+2 (line 192) and, after reverse-complement, compares to 'GT' — but on minus strand that genomic location is the acceptor and reverse-complements to 'AG', so the GT comparison fails and a spurious NON_CAN_SPLICE flag is added; symmetrically the acceptor branch reads exon.start-2..-1 (line 207), RCs it, and compares to 'AG' when that location is actually the donor ('GT'). Net effect: nearly every minus-strand splice variant gets a bogus NON_CAN_SPLICE flag and is downgraded from HC to LC, which is wrong and diverges from Perl VEP LOFTEE. Reachable via --loftee with a reference genome (main.cpp:1925-1927).
- Evidence: // Splice acceptor: 2bp before exon start
if (i > 0 && pos >= exon.start - 2 && pos <= exon.start - 1) { is_splice_site = true; is_splice_acceptor = true; ... }
// Splice donor: 2bp after exon end
if (i < transcript->exons.size() - 1 && pos >= exon.end + 1 && pos <= exon.end + 2) { is_splice_site = true; is_splice_donor = true; ... }
...
if (is_splice_donor && splice_donor_exon_idx >= 0) {
    std::string dinuc = reference_->get_sequence(chrom, exon.end + 1, exon.end + 2);
    if (is_minus) { dinuc = loftee_reverse_complement(dinuc); }
    if (dinuc.length() == 2 && dinuc != "GT") { flags.push_back("NON_CAN_SPLICE"); is_hc = false; }
}
if (is_splice_acceptor && splice_acceptor_exon_idx >= 0) {
    std::string dinuc = reference_->get_sequence(chrom, exon.start - 2, exon.start - 1);
    if (is_minus) { dinuc = loftee_reverse_complement(dinuc); }
    if (dinuc.length() == 2 && dinuc != "AG") { flags.push_back("NON_CAN_SPLICE"); is_hc = false; }
}
- Fix: Make donor/acceptor identification strand-aware in the exon loop (lines 84-98) to mirror vep_annotator.cpp. For plus strand: exon.start-2..-1 = acceptor, exon.end+1..+2 = donor (current behavior). For minus strand: swap them — exon.start-2..-1 = donor, exon.end+1..+2 = acceptor. Then in the NON_CAN_SPLICE check, read the donor dinucleotide from the transcript-5' intronic bases and the acceptor from the transcript-3' bases, applying reverse_complement on minus strand so the canonical comparison is GT for donor and AG for acceptor in transcript orientation. Concretely, base which genomic side each branch reads on the strand, e.g. for the donor branch on minus strand read get_sequence(exon.start-2, exon.start-1) (the transcript-5' side) rather than exon.end+1..+2. Add a regression test using a minus-strand transcript with a canonical GT-AG intron (forward-strand donor bases "AC", acceptor bases "CT") plus a reference genome, asserting that a splice_donor/acceptor variant produces NO NON_CAN_SPLICE flag and stays HC; add a complementary non-canonical case to confirm the flag still fires when appropriate.
- Verifier: Confirmed against current code in src/sources/lof_source.cpp.

Splice-site labeling (lines 84-98) is NOT strand-aware:
  if (i > 0 && pos >= exon.start - 2 && pos <= exon.start - 1) { is_splice_site = true; is_splice_acceptor = true; splice_acceptor_exon_idx = i; }
  if (i < size-1 && pos >= exon.end + 1 && pos <= exon.end + 2) { is_splice_site = true; is_splice_donor = true; splice_donor_exon_idx = i; }
This always treats the exon-start side as ACCEPTOR and exon-end side as DONOR. The consequen


### MEDIUM

**[build-quality/portability] src/vep_annotator.cpp:210,232,238,252,258,281,327,1044,1090,2085,3871-3872** — UB: signed char passed to std::toupper/::toupper without unsigned cast (multiple sites, including FASTA load)
- What: std::toupper/std::tolower (and the C ::toupper) have defined behavior only for arguments representable as unsigned char or EOF. Passing a plain (signed) char with a value > 127 (i.e. negative) is undefined behavior. The codebase is INCONSISTENT: many sites correctly do toupper(static_cast<unsigned char>(c)) (e.g. lines 2249-2250, 2534-2535), but several do not. The most concrete is FASTA sequence loading at lines 1044 and 1090, which run std::transform(line.begin(), line.end(), line.begin(), ::toupper) over arbitrary file bytes; any non-ASCII byte (corrupted FASTA, stray high-bit char) triggers UB. Line 2085 calls std::toupper(cached_cds[idx]) / std::toupper(indel_bases[0]) with no cast at all. Lines 3871-3872 mis-apply the cast to the RESULT (static_cast<unsigned char>(std::toupper(next_base))) while still feeding the raw signed char into toupper, so the UB remains. maxentscan_source.cpp:342 has the same raw-char bug. On platforms where char is signed (x86/arm default) this is real UB; on glibc it can index a table out of bounds.
- Evidence: src/vep_annotator.cpp:1044: std::transform(line.begin(), line.end(), line.begin(), ::toupper);
src/vep_annotator.cpp:2085: if (std::toupper(cached_cds[check_idx]) == std::toupper(indel_bases[0])) {
src/vep_annotator.cpp:3871-3872: if (static_cast<unsigned char>(std::toupper(next_base)) ==
            static_cast<unsigned char>(std::toupper(shifted_indel[0]))) {
- Fix: Add a single helper in a shared header, e.g. `inline char to_upper(char c){ return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); }`, and route all sites through it. For the `std::transform(..., ::toupper)` sites (1044, 1090, 210, 232, 238, 252, 258, 281, 327) replace the `::toupper` function pointer with `[](unsigned char c){ return std::toupper(c); }` (or the helper). For line 2085 wrap each argument: `std::toupper(static_cast<unsigned char>(cached_cds[check_idx]))` and likewise for `indel_bases[0]`. For lines 3871-3872 move the cast to the INPUT: `std::toupper(static_cast<unsigned char>(next_base))`. Fix maxentscan_source.cpp:342 the same way: `switch (std::toupper(static_cast<unsigned char>(base)))`. Prioritize the FASTA-load sites (1044/1090) as the only realistically reachable trigger.
- Verifier: Verified all cited locations against current code; every claim matches exactly.

UNSAFE sites confirmed (raw signed char fed to toupper):
- src/vep_annotator.cpp:1044 and :1090 (the two FASTA-load branches, gzipped and plain): `std::transform(line.begin(), line.end(), line.begin(), ::toupper);` operating on raw file bytes. This is the most concrete exposure since `line` holds arbitrary bytes from the input file; a high-bit byte yields a negative `char` index into the ctype table = UB.
- src/vep_

**[build-quality/portability] CMakeLists.txt:18,26-28,181** — -march=native enabled by default for Release combined with install() ships non-portable binaries
- What: ENABLE_NATIVE defaults to ON (line 18) and unconditionally appends -march=native to CMAKE_CXX_FLAGS_RELEASE for Release builds (the default build type, line 13-15). The same CMakeLists then install()s vep_annotator/filter_vep to bin (line 181). A binary built with -march=native on a modern build host (e.g. with AVX-512/AVX2) will execute illegal instructions (SIGILL) when run or deployed on an older/different CPU. For a clinical tool that is frequently built on one machine and deployed/distributed to others (clusters, containers, colleagues), this is a real portability hazard that silently produces crashes rather than wrong output. -march=native should be opt-in, not the default, especially for installed artifacts.
- Evidence: option(ENABLE_NATIVE "Enable -march=native for optimized builds" ON)
...
if(ENABLE_NATIVE AND (CMAKE_BUILD_TYPE STREQUAL "Release"))
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
endif()
...
install(TARGETS vep_annotator filter_vep DESTINATION bin)
- Fix: Change line 18 to `option(ENABLE_NATIVE "Enable -march=native (host-CPU-only, non-portable binaries)" OFF)` so portable binaries are the default and native is explicit opt-in. Alternatively, keep an optimized default that stays portable by emitting `-mtune=native` (tuning only, no host-specific ISA) for Release, and only add `-march=native` when ENABLE_NATIVE is explicitly set. Either way, document in README that -march=native binaries must be run only on the build CPU's microarchitecture or newer, and consider not auto-applying native flags to the install() targets.
- Verifier: All cited lines match the actual file verbatim. CMakeLists.txt line 18: `option(ENABLE_NATIVE "Enable -march=native for optimized builds" ON)` — defaults ON. Lines 13-15: `if(NOT CMAKE_BUILD_TYPE) set(CMAKE_BUILD_TYPE Release) endif()` — Release is the default build type. Lines 26-28: `if(ENABLE_NATIVE AND (CMAKE_BUILD_TYPE STREQUAL "Release")) set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native") endif()` — unconditionally appends -march=native for the default Release build. L

**[build-quality/performance] src/vep_annotator.cpp:1041** — reserve(300000000) per FASTA contig wastes ~300MB for multi-contig assemblies
- What: On every '>' header line the FASTA loader does current_seq.reserve(300000000) (~300MB). This is fine for a handful of human chromosomes, but for a draft/scaffold assembly or transcriptome FASTA with thousands of short contigs it forces a 300MB allocation per contig (and current_seq was just std::move'd out, so it reallocates fresh). This can blow up peak memory and cause heavy allocator pressure / OOM on non-human inputs, while providing no benefit for short sequences.
- Evidence: current_seq.clear();
current_seq.reserve(300000000);  // Reserve ~300MB for human chromosomes
- Fix: Remove the eager reserve(300000000) at both line 1041 and line 1088. std::string growth is geometric, so large chromosomes incur only O(log n) reallocations and acceptable cost. If a hint is desired for the human common case, reserve a modest amount (e.g. a few MB) instead. Critically, to avoid retaining over-allocated buffers in the persistent map, either avoid reserving more than needed, or call current_seq.shrink_to_fit() immediately before each std::move into pimpl_->sequences (lines 1028, 1052, 1077, 1097). The minimal safe change is to delete both reserve(300000000) lines.
- Verifier: The cited code matches the claim. In src/vep_annotator.cpp the FASTA loader reserves 300MB per contig in BOTH code paths:

gzipped path (line 1040-1041):
  current_seq.clear();
  current_seq.reserve(300000000);  // Reserve ~300MB for human chromosomes

plain-text path (line 1087-1088):
  current_seq.clear();
  current_seq.reserve(300000000);

This runs on every '>' header line. The reviewer is correct that current_seq was just std::move'd out (lines 1028/1077: pimpl_->sequences[current_chrom] = 

**[core-consequence/parity-gap] src/vep_annotator.cpp:360-367** — Consequence terms emitted in collection order, not severity order (parity gap with Perl VEP)
- What: get_consequence_string() concatenates ann.consequences in the order they were pushed during determine_consequences(): splice terms first (pushed in the exon loop, then enum-sorted at lines 2825-2828), then UTR terms, then the specific coding term (missense/synonymous/etc.) appended last. The final vector is therefore NOT sorted by severity. Perl VEP emits the Consequence column with terms ordered most-severe first. Example: a missense variant that is also in a splice region yields vector [SPLICE_REGION_VARIANT, MISSENSE_VARIANT] and prints 'splice_region_variant&missense_variant', whereas Perl VEP prints 'missense_variant&splice_region_variant'. This affects the default TSV output (output_writer.hpp format_consequence -> get_consequence_string) and the JSON transcript_consequences array, causing string-level mismatches against Perl VEP for any transcript with more than one consequence term.
- Evidence: std::string VariantAnnotation::get_consequence_string() const {
    std::string result;
    for (const auto& c : consequences) {
        if (!result.empty()) result += '&';
        result += consequence_to_string(c);
    }
    return result;
}
- Fix: Sort `ann.consequences` by get_consequence_rank() ascending (most-severe first) ONCE after the vector is fully assembled — e.g. immediately after the NMD push in annotate_transcript (src/vep_annotator.cpp, right after line 1972), using a stable sort to keep deterministic tie ordering. Doing it there automatically corrects TSV (format_consequence/get_consequence_string), JSON consequence_terms, and VCF CSQ without duplicating sort logic in output_writer.hpp. Caution: get_consequence_rank ranks SPLICE_REGION_VARIANT(16) below MISSENSE_VARIANT(13) which is what we want, but note the enum value order is NOT identical to rank order for splice sub-types (enum places SPLICE_DONOR_5TH_BASE/REGION/POLYPYRIMIDINE before SPLICE_REGION while ranks interleave 15,17,18,16), so the existing enum-value std::sort at lines 2825-2828 must NOT be relied upon for severity ordering — keep it only for dedup, and add the explicit rank-based sort. After fixing, the production tests that build real annotations will match, and the pre-sorted writer unit tests (test_output.cpp:2095, :2743) remain valid.
- Verifier: CONFIRMED. The cited code matches exactly. `VariantAnnotation::get_consequence_string()` (src/vep_annotator.cpp:360-367) joins `consequences` purely in stored-vector order with no severity sort:
```
for (const auto& c : consequences) { if (!result.empty()) result += '&'; result += consequence_to_string(c); }
```
The stored order is set by push order in `determine_consequences()`. Tracing a missense-in-splice-region SNV: the exon loop pushes SPLICE_REGION_VARIANT first (line 2811); the only sort 

**[core-coord/correctness] src/vep_annotator.cpp:3275-3291 and 2137-2162** — Indel codon/frameshift extraction assumes ref bases are contiguous in CDS space; deletions spanning an exon/intron boundary splice the CDS at the wrong offset
- What: get_affected_codons builds the mutated CDS as `cds_seq.substr(0, cds_offset) + alt + cds_seq.substr(cds_offset + ref_len)` using `ref_len = ref.length()` directly as a CDS span. The frameshift fsTer scan in annotate_coding_region (lines 2137-2148) does the same. ref.length() is a GENOMIC span; for a deletion that spans an exon/intron junction it includes intronic bases that are absent from cds_seq, so `cds_offset + ref_len` removes too many CDS bases (and on minus strand combines with the C1 over-subtraction). The resulting alt codon, amino acids, HGVSc deletion length, and the fsTer distance are all computed from a corrupted mutated CDS, yielding wrong codons/p. notation for boundary-spanning indels.
- Evidence: mut_cds = cds_seq.substr(0, cds_offset) + alt + cds_seq.substr(cds_offset + ref_len);
...
            std::string rc_alt = reverse_complement_sequence(alt);
            int cds_offset = cds_pos - 1; // 0-based start of affected region
            if (cds_offset < 0 || cds_offset + ref_len > static_cast<int>(cds_seq.length())) {
                return {ref_codon, ref_codon};
            }
            mut_cds = cds_seq.substr(0, cds_offset) + rc_alt + cds_seq.substr(cds_offset + ref_len);
- Fix: Compute the actual number of deleted CDS bases by intersecting the genomic deletion span [pos, var_end] with each transcript.cds_region (summing overlap lengths), and use that CDS-deletion width as the splice argument (replacing ref_len = ref.length()) in both get_affected_codons (lines 3280/3290) and the fsTer mutated-CDS construction (lines 2142/2147). Likewise base the frameshift/inframe length_diff classification (lines 2939-2986) and the minus-strand cds_pos adjustment (lines 2000/2918) on the CDS-overlap length rather than ref.length(). For deletions whose deleted CDS bases are non-contiguous (span more than one CDS region), do not attempt codon/HGVSp construction from a single contiguous splice — fall back to the coding_sequence_variant / splice-spanning path and emit the splice consequence without precise codon/p. notation, matching Perl VEP behavior. Add tests covering a plus-strand and a minus-strand coding deletion that spans an exon donor/acceptor boundary.
- Verifier: Confirmed against the actual current code. In get_affected_codons (src/vep_annotator.cpp), the indel branch uses the raw genomic ref length as a CDS splice width:

  int ref_len = static_cast<int>(ref.length());                                   // line 3273
  ... mut_cds = cds_seq.substr(0, cds_offset) + alt + cds_seq.substr(cds_offset + ref_len);   // line 3280 (plus)
  ... mut_cds = cds_seq.substr(0, cds_offset) + rc_alt + cds_seq.substr(cds_offset + ref_len); // line 3290 (minus)

The frames

**[core-coord/parity-gap] src/vep_annotator.cpp:2044-2063 and 2166-2174** — Multi-codon MNV produces HGVSp using only the first changed amino acid, dropping the remaining changed residues
- What: For a same-length substitution spanning more than one codon (e.g. a 4-6 bp MNV), the code correctly builds full ref/alt amino-acid strings (ref_aa_str/alt_aa_str, lines 2052-2062) and sets ann.amino_acids to the full multi-residue string, but it then passes only `ref_aa_str[0]` and `alt_aa_str[0]` to generate_hgvsp (lines 2166-2168). The resulting HGVSp reports a single missense substitution (e.g. p.Ala123Thr) even when two adjacent residues change, contradicting the AMINO_ACIDS field (e.g. AG/TV) and Perl VEP, which would emit a delins/range such as p.Ala123_Gly124delinsThrVal. The protein change is therefore silently truncated.
- Evidence: ann.hgvsp = generate_hgvsp(
        CodonTable::get_three_letter(ref_aa_str[0]),
        CodonTable::get_three_letter(alt_aa_str[0]),
        ann.protein_position,
        transcript,
        ann.consequences,
        hgvs_protein_end,
        hgvs_fs_ter_dist,
        hgvs_end_ref_aa
    );
- Fix: In annotate_coding_region, when ref_aa_str.length() > 1 (multi-codon same-length MNV), do not emit a single-residue missense. Instead compute the protein range and emit a protein delins covering all affected residues, mirroring the existing PROTEIN_ALTERING delins branch in generate_hgvsp. Concretely: set hgvs_protein_end = ann.protein_end (already computed at line 2035) and hgvs_end_ref_aa = three-letter name of ref_aa_str.back(), and pass the full alt residues (three-letter concatenation of all alt_aa_str chars) so generate_hgvsp can produce p.<refAA0><protein_position>_<refAAlast><protein_end>delins<allAltAAs> (e.g. p.Ala123_Gly124delinsThrVal). This requires either tagging the consequence as PROTEIN_ALTERING_VARIANT for multi-codon MNVs or adding an explicit multi-residue-substitution branch in generate_hgvsp that accepts multi-char alt AA strings. Also handle the edge case where only some residues in the range change (still emit the full delins range, matching Perl VEP) and where a residue is unchanged so the synonymous '=' branch is not triggered incorrectly. Keep the single-residue missense form only when exactly one codon is affected (ref_aa_str.length() == 1). Add a unit test for a 4-6 bp same-length MNV spanning two codons asserting the delins HGVSp and consistency with AMINO_ACIDS.
- Verifier: Verified against actual code in /Users/davidcho/Projects/VEP_ensemble_C/src/vep_annotator.cpp.

The cited code matches. In annotate_coding_region, the multi-codon MNV branch correctly builds full AA strings:
  } else {
      // Multi-codon MNV: translate each codon in the range
      for (size_t i = 0; i + 2 < ref_codon.length(); i += 3) { ref_aa_str += translate_with_sec(ref_codon.substr(i,3), ...); }
      for (size_t i = 0; i + 2 < alt_codon.length(); i += 3) { alt_aa_str += translate_with_se

**[core-hgvs-gen/correctness] src/vep_annotator.cpp:2083-2085** — CDS-space 3' shift for deletions checks the wrong base (off-by-one), mis-shifting HGVSc in tandem repeats
- What: In annotate_coding_region's CDS-space 3' shifting loop, the deletion branch compares cached_cds[shift_pos + indel_len - 1] against indel_bases[0]. cds_pos/shift_pos is the 1-based CDS position of the VCF anchor base, so the deleted region is 0-based [shift_pos, shift_pos+indel_len-1] and the base immediately 3' of the deletion is 0-based shift_pos+indel_len. The code instead reads shift_pos+indel_len-1 = the LAST deleted base. It therefore compares the last-deleted-base position against the first-deleted-base value. For homopolymer runs this coincidentally yields the right shift, but for tandem repeats of period >1 it shifts incorrectly. Example: CDS '...CATATATG...', VCF ref=CAT/alt=C (delete 'AT'), indel_bases='AT', shift_pos points at 'C'. The base after the deletion is 'A' (==first deleted base) so HGVS should slide right through the AT-repeat, but check_idx lands on the last deleted 'T' (!= 'A'), so the loop stops and HGVSc is emitted at the un-shifted (non-3'-most) position, disagreeing with Perl VEP which 3'-shifts. This produces a clinically incorrect HGVSc for repeat-region deletions. (Matches the long-standing 'CDS 3' shift off-by-one' note in the project bug list.) The insertion branch (check_idx = shift_pos) is correct.
- Evidence: int check_idx = is_cds_insertion ? shift_pos : (shift_pos + indel_len - 1);
                if (check_idx >= cds_len) break;
                if (std::toupper(cached_cds[check_idx]) == std::toupper(indel_bases[0])) {
- Fix: In src/vep_annotator.cpp line 2083, change the deletion branch index to point at the base immediately 3' of the deleted region: `int check_idx = is_cds_insertion ? shift_pos : (shift_pos + indel_len);` keeping the comparison against indel_bases[0]. This makes the loop slide right only when the base after the deletion equals the first deleted base, matching Perl VEP 3'-shifting for tandem repeats. Add an end-to-end test with a CDS containing a period>1 tandem repeat (e.g. ...CATATATG..., ref=CAT/alt=C) asserting the HGVSc is shifted to the 3'-most position, plus a homopolymer regression case to confirm no behavior change there. Also correct the misleading "intentional" note at parity_details.md:21.
- Verifier: Cited code matches exactly (src/vep_annotator.cpp lines 2082-2090):
  for (int i = 0; i < 1000; ++i) {
      int check_idx = is_cds_insertion ? shift_pos : (shift_pos + indel_len - 1);
      if (check_idx >= cds_len) break;
      if (std::toupper(cached_cds[check_idx]) == std::toupper(indel_bases[0])) {
          std::rotate(indel_bases.begin(), indel_bases.begin() + 1, indel_bases.end());
          shift_pos++; } else { break; } }

Index/anchor conventions verified from the actual code:
- calcu

**[core-hgvs-gen/parity-gap] src/vep_annotator.cpp:2014, 2597** — HGVSg is never 3'-shifted, while HGVSc always is, producing inconsistent genomic vs coding HGVS for shiftable indels
- What: populate_transcript_metadata is called with the original, un-shifted pos/ref/alt (line 2014), and generate_hgvsg uses those raw VCF coordinates (line 2597). There is no 3' (right) shift applied to the HGVSg description. Meanwhile annotate_coding_region ALWAYS 3'-shifts the HGVSc (the CDS-space shift at lines 2065-2092 runs unconditionally for any indel, independent of the shift_3prime_/shift_genomic_ flags). Ensembl Perl VEP 3'-shifts all HGVS notations (g., c., n., p.) by default. The result is that for any shiftable indel the C++ tool emits an HGVSc placed at the 3'-most position but an HGVSg placed at the original VCF position, so the two notations describe different loci of the same variant. This is a parity gap that yields clinically misleading HGVSg strings and internally inconsistent records.
- Evidence: populate_transcript_metadata(chrom, pos, ref, alt, transcript, hgvs_offset, ann);
...
    ann.hgvsg = generate_hgvsg(chrom, pos, ref, alt, ref_context);
- Fix: Make HGVSg use the same 3'-shifted coordinates that the rest of HGVS already uses by default. annotate_transcript already computes shifted_pos/shifted_ref/shifted_alt via right_normalize, but currently only when (shift_3prime_ || shift_genomic_) is set (line 1935). Since the codebase intends to "always shift" HGVS (see main.cpp:1588-1589 no-op comment), either (a) compute the right_normalized coordinates unconditionally for indels and pass those to populate_transcript_metadata at line 2014 instead of raw pos/ref/alt, or (b) add a dedicated shift_hgvs setting (default true, matching Perl VEP shift_hgvs=1) that gates the HGVSg shift, decoupled from shift_genomic_ (which should continue to only affect the reported genomic Location). Concretely: change line 2014 to pass the HGVS-shifted coordinates, and have populate_transcript_metadata feed those into generate_hgvsg at line 2597 (recomputing ref_context from the shifted pos). Verify minus-strand handling, since right_normalize/the CDS shift already account for strand. Add an integration test that runs a shiftable indel (e.g. an insertion/deletion in a homopolymer or tandem repeat) through annotate_transcript and asserts HGVSc and HGVSg describe the same 3'-most locus.
- Verifier: I verified every element of the claim against the actual current code.

1. HGVSc is ALWAYS 3'-shifted regardless of flags. In annotate_coding_region (src/vep_annotator.cpp:2065-2094) the shift block is gated only on `if (ref.length() != alt.length())` -- any indel -- with no reference to shift_3prime_/shift_genomic_. It rotates indel_bases and advances shift_pos in CDS space, then `hgvs_cds_pos = shift_pos` feeds generate_hgvsc at line 2096. The non-coding path annotate_noncoding_hgvsc (lines 25

**[core-hgvs-gen/correctness] src/vep_annotator.cpp:2569-2573** — Minus-strand multi-base deletions in non-coding (n.) transcripts produce wrong cDNA coordinates and reversed range
- What: In annotate_noncoding_hgvsc, the exonic n. notation uses pos_str = cdna(hgvs_pos), where hgvs_pos is the genomically-first changed base (pos+1 for anchored indels, possibly 3'-shifted). For a multi-base deletion on a MINUS-strand transcript, the genomically-first deleted base is the transcript/cDNA-LAST (highest cDNA coordinate) base of the deletion, so pos_str is the higher cDNA coordinate, and the offset lambda then ADDS delta (cdna_pos + delta), making the end coordinate even higher. HGVS requires start < end (start_endel), so a 2+ base minus-strand non-coding deletion is emitted as e.g. n.101_102del instead of the correct n.98_99del — wrong start position and inverted range direction. The coding path avoids this by adjusting cds_pos to the FIRST (lowest) CDS position for minus-strand multi-base variants (lines 1999-2002), but the non-coding path performs no such adjustment.
- Evidence: int cdna_pos = calculate_cdna_position(hgvs_pos, transcript);
    if (cdna_pos > 0) {
        std::string pos_str = std::to_string(cdna_pos);
        ann.hgvsc = build_hgvsc(pos_str, ref, alt, transcript, transcript_version_, "n",
            [cdna_pos](int delta) { return std::to_string(cdna_pos + delta); });
- Fix: Mirror the coding-path adjustment (lines 1999-2002) in annotate_noncoding_hgvsc. For minus-strand transcripts with a multi-base deletion (and same logic for delins/MNV where the genomic-first base is the cDNA-last), compute the start from the transcript-first (lowest cDNA) deleted base. Concretely, after computing hgvs_pos, for transcript.strand == '-' and a deletion with ref.length() > alt.length(), use the cDNA position of the genomically-LAST deleted base as the start: e.g. int start_cdna = calculate_cdna_position(hgvs_pos + del_len - 1, transcript); then pass start_cdna as pos_str with the existing +delta lambda. A more robust approach: compute cDNA coords of both deletion endpoints, take std::min as start and std::max as end, and build the range from those so the output is correctly anchored regardless of strand. Add a regression test for a minus-strand non-coding (e.g. lincRNA/NR_) transcript with a 2+ base deletion verifying the expected n.{lower}_{higher}del.
- Verifier: The cited code at src/vep_annotator.cpp:2569-2573 matches the claim exactly:

  int cdna_pos = calculate_cdna_position(hgvs_pos, transcript);
  if (cdna_pos > 0) {
      std::string pos_str = std::to_string(cdna_pos);
      ann.hgvsc = build_hgvsc(pos_str, ref, alt, transcript, transcript_version_, "n",
          [cdna_pos](int delta) { return std::to_string(cdna_pos + delta); });
  }

I traced the full path and confirm the bug for multi-base deletions on minus-strand non-coding transcripts:

1.

**[core-hgvs-gen/correctness] src/vep_annotator.cpp:3268, 3294-3295, 3713** — Multi-codon inframe insertions/delins truncated to a single residue in HGVSp, amino_acids and codons
- What: get_affected_codons (indel branch) always extracts exactly 3 bases for both ref_codon and alt_codon (cds_seq.substr(codon_start,3) and mut_cds.substr(codon_start,3)), so for an inframe insertion or delins that spans more than one codon (e.g. a 6-base in-frame insertion) only the first inserted/changed codon is captured. annotate_coding_region then passes only CodonTable::get_three_letter(alt_aa_str[0]) as alt_aa, and generate_hgvsp's inframe-insertion branch emits 'ins' + that single residue. Consequently a 2-residue in-frame insertion is reported as inserting one amino acid (and the codons/amino_acids columns are likewise truncated), diverging from Perl VEP which lists all inserted residues. Output is wrong for any in-frame insertion/delins of two or more codons.
- Evidence: std::string ref_codon = cds_seq.substr(codon_start, 3);
...
            alt_codon = mut_cds.substr(codon_start, 3);
...
        oss << ref_aa << protein_pos << "_" << next_aa << (protein_pos + 1) << "ins" << alt_aa;
- Fix: In get_affected_codons's indel branch (src/vep_annotator.cpp ~3261-3300), extract the full affected codon range rather than a fixed 3 bases. Compute the last affected codon from the mutated CDS span (codon_start through the codon containing the end of the inserted/changed region), set ref_codon = cds_seq.substr(codon_start, ref_span) and alt_codon = mut_cds.substr(codon_start, alt_span) where the spans cover whole codons. Then annotate_coding_region's existing multi-codon translation branch (line 2052) will translate each codon, producing full ref_aa_str/alt_aa_str. Finally, pass the complete residue strings to generate_hgvsp instead of only alt_aa_str[0]/ref_aa_str[0] (lines 2167-2168), and update generate_hgvsp's inframe-insertion (line 3713) and delins (line 3721) branches to render the full residue string. Optionally add tandem-duplication ('dup') detection as a follow-up, but the core fix is propagating all inserted/changed residues. Add a regression test for a >=2-codon in-frame insertion/delins.
- Verifier: Verified against actual code in src/vep_annotator.cpp. The claim is fully substantiated.

In get_affected_codons, the indel branch (ref.length() != alt.length(), line 3261) always extracts exactly 3 bases:
- line 3268: `std::string ref_codon = cds_seq.substr(codon_start, 3);`
- line 3295: `alt_codon = mut_cds.substr(codon_start, 3);`
The mutated CDS (mut_cds) is correctly built with the full insertion (line 3280: `cds_seq.substr(0, cds_offset) + alt + cds_seq.substr(cds_offset + ref_len)`), but 

**[filter-transcript/parity-gap] include/filter_vep.hpp:422-450** — regex/re operator is unreachable through the -f filter expression parser
- What: parse_filter_operator() recognizes "regex"/"re" and apply_condition() implements a real std::regex_search path, but parse_filter_expression() (the only thing that processes a user '-f EXPR' string) never lists " regex " or " re " in its operators vector. So an expression like 'Consequence regex splice' finds NO operator and falls back to treating the entire string as a field-EXISTS check (field = "Consequence regex splice"). The REGEX code path is therefore dead from the CLI: users cannot ever run a regex match. This is both a parity gap with Perl VEP filter_vep (whose 'match'/regex operator is a true Perl regex) and a silent wrong-result (the expression is misinterpreted as an EXISTS check on a bogus field, which is almost always false, dropping records the user intended to keep). The project's own test (tests/test_filter_vep.cpp:208-223) documents this exact misbehavior.
- Evidence: operators.push_back(" contains "); operators.push_back(" in "); operators.push_back(" match "); operators.push_back(" exists"); ... // no " regex " / " re " ;  and test_filter_vep.cpp: "regex is not in the list, so no operator would be found; it falls back to EXISTS." EXPECT_EQ(cond.field, "Consequence regex splice");
- Fix: Add " regex " (and optionally " re ") to the operators vector in parse_filter_expression (include/filter_vep.hpp ~line 432, near " match "). Caution: the operator-selection loop (lines 444-450) picks the token with the smallest position, not the longest match — so adding bare " re " risks false matches against a field/value containing " re ". Prefer adding only " regex " (low collision risk), or change the selection logic to prefer the longest token at a given position before adding short tokens. After adding it, update tests/test_filter_vep.cpp:207-224 (RegexOperator) to assert cond.op == FilterOperator::REGEX, cond.field == "Consequence", cond.value == "splice", and add a positive apply_condition test exercising std::regex_search. Consider also adding " not contains "/" not in " for fuller parity, but that is separate from this fix.
- Verifier: Verified against current code. CLI path: src/filter_vep.cpp:130-133 handles `-f/--filter EXPR` by calling `vep::parse_filter_expression(expr)` and pushing the result into config.conditions — this is the only entry point for user filter expressions.

parse_filter_expression (include/filter_vep.hpp:418-457) builds its operators vector at lines 422-439: " is ", " eq ", " ne ", " gt ", " ge ", " lt ", " le ", " contains ", " in ", " match ", " exists", ">=", "<=", "!=", ">", "<", "=". There is NO " 

**[filter-transcript/parity-gap] include/filter_vep.hpp:61** — 'match' operator is substring, not regex, contrary to Perl VEP filter_vep semantics
- What: Perl VEP filter_vep's 'match' operator performs a Perl regular-expression match (=~). Here parse_filter_operator maps both 'contains' and 'match' to FilterOperator::CONTAINS, which apply_condition implements as a plain std::string::find substring test (line 228). A user filter such as 'Consequence match ^missense' or 'HGVSc match c\.[0-9]+A>G' will be treated as a literal substring search and silently return wrong matches relative to Perl VEP. Combined with C1 (regex unreachable), there is currently no way to get true regex matching from the CLI.
- Evidence: if (lower == "contains" || lower == "match") return FilterOperator::CONTAINS;  ... } else if (cond.op == FilterOperator::CONTAINS) { result = (value.find(cond.value) != std::string::npos); }
- Fix: Two coordinated changes are needed to restore Perl VEP parity:

1. In parse_filter_operator (include/filter_vep.hpp line 61), separate `match` from `contains`: keep `contains` -> FilterOperator::CONTAINS (substring) and map `match` -> FilterOperator::REGEX. This requires updating the pinning test at tests/test_filter_vep.cpp lines 75-78 (which currently asserts `match` -> CONTAINS) and line 178.

2. Wire regex into the CLI: in parse_filter_expression (lines 422-439), add `" regex "` and `" re "` to the operators list so the REGEX op is reachable (this fixes C1). With change #1, `" match "` will then resolve to REGEX as well. Update the test at lines 207-224 accordingly.

3. Update CLI help text (src/filter_vep.cpp line 27) to advertise `match`/`regex` and README.md line 532 to accurately describe behavior. The existing REGEX implementation (lines 246-253) already falls back to substring on std::regex_error, which is a reasonable safety net for invalid patterns.

If full parity is not desired, the minimum acceptable fix is to remove `match` from parse_filter_operator/parse_filter_expression entirely (so it is not silently accepted with wrong semantics) and correct the misleading README line 532. But mapping `match` -> regex is the parity-correct option.
- Verifier: Confirmed against the actual code in include/filter_vep.hpp.

Line 61 matches the quote exactly: `if (lower == "contains" || lower == "match") return FilterOperator::CONTAINS;`. The CONTAINS branch in apply_condition (line 228) is a plain substring test: `result = (value.find(cond.value) != std::string::npos);`. So `match` is implemented as substring containment, not regex.

In Ensembl Perl VEP filter_vep, the `match` operator performs a Perl regex match (`$value =~ /$pattern/`). Mapping it to s

**[filter-transcript/correctness] include/transcript_filter.hpp:451-456** — Custom --pick-order with no recognized tokens silently disables all ranking
- What: parse_pick_order() simply skips any unrecognized token and returns whatever it accumulated; an all-invalid string (e.g. a typo like '--pick-order canoncial,rnk') returns an EMPTY vector, which main.cpp assigns directly to filter_config.pick_order with no validation (src/main.cpp:1030). is_better() then iterates an empty pick_order and always returns false, so std::sort in pick_one becomes a no-op and the tool just keeps annotations[0] (arbitrary input order). The user gets a wrong/arbitrary 'picked' transcript with no error or warning. Perl VEP validates pick_order tokens and dies on unknown values.
- Evidence: bool is_better(...) const { for (const auto& criteria : config_.pick_order) { int cmp = compare_by_criteria(a,b,criteria); if (cmp!=0) return cmp<0; } return false; }  // empty pick_order -> always false ; main.cpp: filter_config.pick_order = vep::parse_pick_order(argv[++i]); // no emptiness/validity check
- Fix: In parse_pick_order, detect unrecognized tokens (e.g. collect them or return a status). In main.cpp's --pick-order handler, if the resulting pick_order is empty or contained invalid tokens, either (a) exit with a clear error message listing the valid criteria, or (b) emit a warning to stderr and fall back to the default order. The simplest robust fix: after assignment, `if (filter_config.pick_order.empty()) { fprintf(stderr, "Error: --pick-order contained no valid criteria\n"); return 1; }`, and ideally also warn when some-but-not-all tokens were dropped. This mirrors Perl VEP, which validates pick_order tokens.
- Verifier: The cited code matches the claim exactly. `parse_pick_order` (include/transcript_filter.hpp:723-756) tokenizes on commas and only pushes a PickCriteria for recognized tokens; the if/else chain has NO terminal `else`, so any unrecognized token is silently dropped: `if (token=="canonical")...else if(...)...else if (token=="refseq")...`. An all-invalid input (e.g. "canoncial,rnk") therefore returns an empty vector.

main.cpp:1029-1030 assigns the result directly with no validation: `} else if (arg 

**[filter-transcript/correctness] include/transcript_filter.hpp:525-528** — Picking uses non-stable std::sort, so ties are resolved non-deterministically (parity/reproducibility)
- What: pick_one, flag_picked, flag_picked_per_allele, and flag_picked_per_allele_gene all use std::sort with is_better as the comparator. is_better returns false for both (a,b) and (b,a) whenever two transcripts tie on every criterion in pick_order (a legitimately common case, e.g. two non-canonical protein_coding transcripts with equal rank and no length info). std::sort is NOT stable, so the surviving 'picked' transcript among tied candidates is implementation-defined and can differ from input order and from run to run / library to library. Perl VEP's pick is deterministic (it sorts stably and falls back to a defined order, ultimately transcript ID/stable_id). This produces different picked transcripts than Perl VEP and non-reproducible output across platforms.
- Evidence: std::sort(annotations.begin(), annotations.end(), [this](const AnnotationWithMeta& a, const AnnotationWithMeta& b){ return is_better(a, b); });  // is_better returns false on full tie => non-stable reorder
- Fix: Make tie resolution deterministic. Either (a) replace std::sort with std::stable_sort in pick_one (line 525), flag_picked (635), flag_picked_per_allele (683), and flag_picked_per_allele_gene (710) — this preserves input order on ties, which is deterministic given a deterministic input order; or, preferably, (b) add a final deterministic tiebreaker inside is_better/compare_by_criteria after the existing criteria loop, e.g. compare a.annotation.transcript_id vs b.annotation.transcript_id lexicographically (return a.annotation.transcript_id < b.annotation.transcript_id when all pick_order criteria tie). Option (b) most closely matches Perl VEP's deterministic fallback to transcript stable_id and is robust regardless of input ordering or sort algorithm. Add a unit test that asserts the specific surviving transcript_id for a full-tie pair to lock the behavior in.
- Verifier: Verified against the actual current code in include/transcript_filter.hpp.

1) is_better (lines 450-456) returns false when every criterion in pick_order ties:
   for (const auto& criteria : config_.pick_order) { int cmp = compare_by_criteria(a, b, criteria); if (cmp != 0) return cmp < 0; }  return false;
   So a full tie yields is_better(a,b)==false AND is_better(b,a)==false, i.e. the two elements are "equivalent" under this strict-weak-ordering.

2) The sort is std::sort, not std::stable_sort,

**[filter-transcript/correctness] include/filter_vep.hpp:359-368** — --exclude-intronic only drops records whose Consequence is exactly 'intron_variant'
- What: The exclude_intronic filter checks consequence.find("intron_variant") but then only returns false (drops) when consequence == "intron_variant" EXACTLY. VEP commonly emits compound consequence strings like 'intron_variant,non_coding_transcript_variant' or 'splice_region_variant,intron_variant'. For these the inner exact-match check is false, so the record is NOT excluded even though it is an intronic annotation. The result is that --exclude-intronic fails to exclude a large fraction of intronic records, contradicting the option's documented purpose ('Exclude intron variants').
- Evidence: if (consequence.find("intron_variant") != std::string::npos) { // Only filter if it's ONLY an intron variant if (consequence == "intron_variant") { return false; } }
- Fix: Decide and align intent with documentation. Two acceptable options:

1) If the goal is "exclude any annotation that is intronic", change the most-severe-term semantics: split Consequence on ',' and drop the record when the highest-ranked (most severe) term is intron_variant, OR when all terms are intronic/MODIFIER region terms. This correctly excludes "intron_variant,non_coding_transcript_variant" while keeping "splice_region_variant,intron_variant" (splice is more severe). This best matches the documented "Exclude intron variants" while preserving more-severe splice annotations.

2) If the current narrow behavior ("exclude only when the sole consequence is intron_variant") is intended, update the help text at src/filter_vep.cpp:46 to say "Exclude records whose only consequence is intron_variant" so the option name/help is no longer misleading, and update the comment.

Option 1 is preferable given the documented purpose. Either way, add a test for the "intron_variant,non_coding_transcript_variant" case (currently untested) so the chosen semantics are pinned down.
- Verifier: The cited code at include/filter_vep.hpp:359-368 matches the claim verbatim:

  if (config.exclude_intronic) {
      std::string consequence = record.get("CONSEQUENCE");
      if (consequence.empty()) consequence = record.get("Consequence");
      if (consequence.find("intron_variant") != std::string::npos) {
          // Only filter if it's ONLY an intron variant
          if (consequence == "intron_variant") {
              return false;
          }
      }
  }

I verified the semantics of app

**[gene-exon-plugin/parity-gap] include/gene_constraint.hpp:338-355** — Loaded gene-constraint data (pLI/LOEUF/misZ) is never queried or emitted — --constraint/--pli/--loeuf produce no output
- What: main.cpp accepts --constraint, --pli and --loeuf and loads them into the GeneConstraintDB singleton (main.cpp:1728-1745), but get_by_symbol()/get_by_gene_id() (and every GeneConstraint accessor such as pLI, oe_lof_upper, mis_z, get_constraint_level, is_constrained) have zero callers anywhere in src/. The annotator and output writers never look the data back up, and no output column (pLI/LOEUF/mis_z/etc.) is ever produced. The user supplies a constraint file, sees 'Loaded N genes' logged, and gets no constraint annotation in TSV/JSON/VCF output. Versus Perl VEP plugins (LOEUF/pLI/Constraint), this is missing output for a feature the CLI advertises.
- Evidence: Grep across src/ for get_by_symbol|get_by_gene_id|.pLI|oe_lof_upper|.mis_z|get_constraint_level|is_constrained returns no matches outside gene_constraint.hpp. The accessors are defined: 'GeneConstraint get_by_symbol(const std::string& gene_symbol) const { auto it = gene_data_.find(gene_symbol); ... }' (lines 338-344) but never called. main.cpp:1730 'constraint_db.load_gnomad_constraint(constraint_path)' is the only interaction with the DB.
- Fix: Wire the lookup into the annotation path: after per-transcript annotation, look up GeneConstraintDB by gene symbol with fallback to gene_id (get_by_symbol then get_by_gene_id), and when the DB is_loaded(), emit the requested metrics (pLI, oe_lof_upper/LOEUF, mis_z, etc.) into the annotation's custom_annotations/Extra map so they surface in TSV Extra, JSON, and VCF CSQ. Register the corresponding output column names. This requires giving the annotator or output stage access to get_gene_constraint_db() (currently only main.cpp includes the header). Alternatively, if the feature is not intended to ship, remove the --constraint/--pli/--loeuf options and the 'Loaded N genes' log to avoid silently misleading users.
- Verifier: Verified against current code. The CLI advertises and wires the feature but never reads the data back.

1. Usage help advertises the options (src/main.cpp:113-115):
   "  --constraint FILE       gnomAD gene constraint file (TSV)\n"
   "  --pli FILE              pLI scores file (GENE\\tpLI)\n"
   "  --loeuf FILE            LOEUF scores file (GENE\\tLOEUF)\n\n"

2. Options are parsed (src/main.cpp:1203-1207) and data is loaded into the singleton (src/main.cpp:1728-1745), including an informational

**[hgvs-parser/correctness] include/hgvs_parser.hpp:807-836** — generate_hgvsg mis-classifies anchor-mismatched insertions, dropping a base and shifting position
- What: The insertion branch fires whenever `alt.size() > ref.size() && ref.size() <= 1`. When ref is a single base, the code unconditionally strips the first alt char as if it were the VCF anchor (`inserted = alt.substr(1)`), regardless of whether alt[0]==ref[0]. For ref=A alt=GG at pos 100 (which is a delins: replace A with GG -> g.100delinsGG) it emits `g.100_101insG`, dropping a base and producing wrong position and wrong type. The anchor assumption is invalid whenever the ref base is actually replaced.
- Evidence: std::string inserted = alt;
        if (!ref.empty() && ref.size() == 1) {
            inserted = alt.substr(1);
        }
- Fix: Mirror the anchor-strip guard used by the deletion/delins branches. Replace lines 810-812 so the first alt base is only stripped when it actually matches the ref anchor: `if (!ref.empty() && ref.size() == 1 && alt[0] == ref[0]) { inserted = alt.substr(1); }`. Additionally, the branch should only classify as an insertion when that anchor matches; when ref.size()==1 and alt[0] != ref[0], the variant is a replacement and should be routed to the delins branch emitting g.<pos>delins<alt> (e.g. g.100delinsGG for ref=A alt=GG). The simplest correct approach: in the insertion-branch condition, add `&& (ref.empty() || alt[0] == ref[0])`, and let the mismatched case fall through to the complex/delins branch (which would then also need to handle the single-base-ref delins, currently its anchor-strip at line 843 requires ref.size()>1 so a single-base ref would correctly emit g.100delinsGG with no strip). Add a regression test asserting generate_hgvsg("7",100,"A","GG","") == "NC_000007.14:g.100delinsGG".
- Verifier: Verified against current code in /Users/davidcho/Projects/VEP_ensemble_C/include/hgvs_parser.hpp. The insertion branch (lines 807-836) fires whenever `ref.empty() || (alt.size() > ref.size() && ref.size() <= 1)`. Inside it, the anchor strip at lines 810-812 is:

    std::string inserted = alt;
    if (!ref.empty() && ref.size() == 1) {
        inserted = alt.substr(1);
    }

This strips alt[0] based solely on `ref.size() == 1`, with NO check that `alt[0] == ref[0]`. This is inconsistent with th

**[hgvs-parser/correctness] include/hgvs_parser.hpp:507-525** — parse_protein_hgvs uses regex_search (not regex_match) for insertions, accepting malformed input as valid
- What: The protein insertion parser uses `std::regex_search` while every other branch in this function uses `std::regex_match`. regex_search matches any substring, so leading/trailing garbage is accepted and reported as a valid parse. e.g. `garbage_Gly12_Ala13insValTrp` parses as valid=true with ref_aa=Gly, pos=12, alt=ValTrp. This silently accepts corrupt HGVS strings, defeating the validation contract callers rely on (result.valid).
- Evidence: static const std::regex ins_regex("([A-Z][a-z]{2}|[A-Z*])(\\d+)_([A-Z][a-z]{2}|[A-Z*])(\\d+)ins([A-Z][a-z]{2}|[A-Z*])+");
    if (std::regex_search(change, match, ins_regex)) {
- Fix: Replace std::regex_search with std::regex_match at hgvs_parser.hpp:508, and capture the inserted AAs in a group instead of using substr. Change the regex to: "([A-Z][a-z]{2}|[A-Z*])(\\d+)_([A-Z][a-z]{2}|[A-Z*])(\\d+)ins((?:[A-Z][a-z]{2}|[A-Z*])+)" and set result.alt_aa = match[5].str() instead of change.substr(ins_pos + 3). I verified this variant: it accepts Gly12_Ala13insVal (alt=Val) and Gly12_Ala13insValTrp (alt=ValTrp) while rejecting garbage_Gly12_Ala13insValTrp, Gly12_Ala13insValTrp_TRAILING, and Gly12_Ala13insValTrpXYZ123 (all match=0). This both restores the validation contract and eliminates the fragile substr-based alt_aa extraction. Add the three malformed inputs as negative regression tests.
- Verifier: CONFIRMED by direct code inspection and runtime testing. At include/hgvs_parser.hpp:507-525, the insertion branch reads:

  static const std::regex ins_regex("([A-Z][a-z]{2}|[A-Z*])(\\d+)_([A-Z][a-z]{2}|[A-Z*])(\\d+)ins([A-Z][a-z]{2}|[A-Z*])+");
  if (std::regex_search(change, match, ins_regex)) { ... result.valid = true; return result; }

This is the ONLY branch in parse_protein_hgvs using regex_search; every other branch (missense_3letter L451, missense_1letter L462, fs_regex L473, del_regex L

**[hgvs-parser/correctness] include/hgvs_parser.hpp:302-310** — 3'UTR star (*) position sign is discarded in parse_coding_hgvs, conflating *N with N
- What: For c. notation, a leading '*' marks a 3'UTR position (e.g. c.*51 is 51 bases past the stop codon). The parser strips the '*' and stores only the magnitude (start_pos=51) with no flag to mark it as 3'UTR. The HGVSParseResult struct has no field to record the star, so c.*51A>G is indistinguishable from c.51A>G after parsing. The inline comment even says 'Mark as 3'UTR position (positive value with flag)' but no flag exists. This maps 3'UTR variants to the wrong coding coordinate. Verified: `parse_coding_hgvs("X","*51+2A>G")` returns start_pos=51 with no UTR indicator.
- Evidence: std::string pos_str = match[1].str();
        // Handle UTR positions (*123 for 3'UTR, -123 for 5'UTR)
        if (!pos_str.empty() && pos_str[0] == '*') {
            result.start_pos = std::stoi(pos_str.substr(1));
            // Mark as 3'UTR position (positive value with flag)
        } else {
- Fix: Add a boolean field (e.g. bool is_utr3 = false;) to HGVSParseResult and set it to true in each branch when pos_str[0] == '*' (and a separate end_is_utr3 for the second coordinate in del/ins/dup). Then in src/main.cpp, before calling map_cds_to_genomic, translate a 3'UTR position into the correct CDS offset by adding the CDS length (i.e. genomic position = stop-codon-end offset + star value), or add a dedicated map_utr3_to_genomic helper. Until distinguished, the HGVS-input path should at minimum reject or error on c.* positions rather than silently mapping them into the CDS. Apply the same flag-setting at the del/ins/dup star-stripping sites (lines ~338, 350, 374, 385, 403, 415).
- Verifier: The cited code matches the claim exactly. In include/hgvs_parser.hpp parse_coding_hgvs (substitution branch, lines 302-310):

    std::string pos_str = match[1].str();
    // Handle UTR positions (*123 for 3'UTR, -123 for 5'UTR)
    if (!pos_str.empty() && pos_str[0] == '*') {
        result.start_pos = std::stoi(pos_str.substr(1));
        // Mark as 3'UTR position (positive value with flag)
    } else {
        result.start_pos = std::stoi(pos_str);
    }

The HGVSParseResult struct (lines 56-

**[hgvs-parser/parity-gap] include/hgvs_parser.hpp:441-543** — No support for protein delins or extension (ext) notation — parse fails on valid HGVS
- What: parse_protein_hgvs has no branch for protein delins (p.Cys28delinsTrpVal) or protein extension/stop-loss (p.Ter110GlnextTer17, p.*110Qext*17), both of which are standard HGVS protein notations produced by Ensembl VEP. These inputs return valid=false with 'Unrecognized protein HGVS pattern'. Frameshift is recognized but degenerately collapsed (see C7).
- Evidence: result.error_message = "Unrecognized protein HGVS pattern: " + change;
    return result;
- Fix: Add two regex branches to parse_protein_hgvs before the final error return (after the deletion branch). (1) Protein delins: match `([A-Z][a-z]{2}|[A-Z*])(\d+)(?:_([A-Z][a-z]{2}|[A-Z*])(\d+))?delins((?:[A-Z][a-z]{2}|[A-Z*])+)` -> set variant_type=DELINS, populate ref_aa/protein_pos/start_pos/end_pos and alt_aa from the inserted residues. (2) Extension/stop-loss: match `(Ter|\*)(\d+)((?:[A-Z][a-z]{2})|[A-Z*])ext(Ter|\*)(\d+|\?)` -> set ref_aa=Ter, protein_pos, alt_aa to the substituted residue, and store the extension length (e.g. in a dedicated field or note). Add round-trip tests asserting parse_protein_hgvs succeeds on generate_hgvsp's own outputs (`p.Ala123_Gly125delinsVal`, `p.Ter123AlaextTer?`).
- Verifier: Verified directly against current code and by empirical execution. `parse_protein_hgvs` (include/hgvs_parser.hpp:441-543) contains branches only for: 3-letter missense (line 448), 1-letter missense (461), frameshift (472), deletion (487), insertion (507), and synonymous (528). There is NO branch for protein delins (e.g. `Cys28delinsTrpVal`, `Ala123_Gly125delinsVal`) nor for extension/stop-loss (`Ter110GlnextTer17`, `Ter123AlaextTer?`). Any unmatched input falls through to line 542: `result.error

**[main-cli/correctness] src/main.cpp:2379-2402** — --fork worker exception is swallowed: failed query yields empty annotation and program exits 0 (silent data loss)
- What: In the parallel pre-annotation, each worker catches all exceptions into thread_exceptions[t]. After join, the exceptions are merely logged at ERROR level; processing then continues, results[i] for the query that threw stays default-constructed (empty vector), it is moved into annotation_cache, and the consumer emits zero annotations for that variant. main() still returns 0. Unlike the single-thread path (where annotate() throwing propagates to the outer try/catch and returns exit code 1), a fork run that hits the same fault produces a partial output file with no error exit code, so downstream pipelines cannot detect the failure. This is a divergence in error handling/exit codes between --fork 1 and --fork N.
- Evidence: } catch (...) {
    thread_exceptions[t] = std::current_exception();
}
...
for (int t = 0; t < fork_count; t++) {
    if (thread_exceptions[t]) {
        try { std::rethrow_exception(thread_exceptions[t]); }
        catch (const std::exception& e) {
            vep::log(vep::LogLevel::ERROR, "Thread " + std::to_string(t) + " failed: " + e.what());
        }
    }
}
... return true;
- Fix: After the post-join loop that inspects `thread_exceptions`, propagate the failure instead of only logging it. Concretely, in the loop at src/main.cpp:2387-2396 either (a) rethrow the first non-null `thread_exceptions[t]` so it reaches the outer try/catch in main() and returns 1, or (b) set a `bool fork_had_error = true` flag, and after the consumer loop completes have main() `return 1` when the flag is set. Rethrowing is simplest and matches single-thread semantics exactly: e.g. keep the ERROR log for each failed thread, then after the loop `for (auto& ep : thread_exceptions) if (ep) std::rethrow_exception(ep);`. Optionally also avoid emitting the empty cached result for the specific throwing query, but the primary fix is restoring the non-zero exit code so downstream pipelines can detect the failure.
- Verifier: The cited code matches the claim and the data-loss/exit-code divergence is real and verifiable end-to-end.

Worker (src/main.cpp:2365-2383): each thread runs `while ((i = work_idx.fetch_add(1)) < queries.size()) { ... results[i] = ann.annotate(...); }` wrapped in `try { ... } catch (...) { thread_exceptions[t] = std::current_exception(); }`. If `ann.annotate()` throws on query `i`, the catch terminates that thread's loop. `results[i]` for the throwing query stays default-constructed (empty `std:

**[output-writers/correctness] include/output_writer.hpp:1053-1071, 1017-1018, 979-980, 907-908, 928-929** — JSON numeric type-coercion emits invalid JSON for leading-zero values
- What: is_valid_json_number() returns true for strings with leading zeros such as "01", "0123", "007", "00", and "-0123". When a custom-annotation value (or SIFT/PolyPhen/TSL score) matches this, it is written as a BARE token (json << pair.second) instead of a quoted string. RFC 8259 forbids leading zeros in numbers, so a value like "0123" produces `"key": 0123` which a strict JSON parser (jq, Python json, JS JSON.parse) rejects, corrupting the entire NDJSON line. Any all-digit, leading-zero custom field (zero-padded IDs, certain allele/count strings, dbNSFP fields) triggers this. Confirmed by compiling the function in isolation: all of "0123","01","-0123","00","007" returned NUMBER(bare).
- Evidence: if (i >= s.size() || !std::isdigit(static_cast<unsigned char>(s[i]))) return false;
    while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) ++i;  // accepts leading zeros
...
if (is_valid_json_number(pair.second)) {
    json << "        \"" << escape_json(key) << "\": " << pair.second << ",\n";
- Fix: Add a leading-zero guard to is_valid_json_number() right after the integer-part's first digit is confirmed (line ~1057, before the digit-consuming while loop): `if (s[i] == '0' && i + 1 < s.size() && std::isdigit(static_cast<unsigned char>(s[i+1]))) return false;`. Verified in isolation that this rejects "0123","01","-0123","00","007","09.5" while still accepting all valid JSON numbers ("0","0.5","-0","-0.0","0.0","10","1.5e3","123"). No call-site changes needed since they all already fall through to the quoted-string branch when the function returns false. Consider adding a GoogleTest case covering leading-zero inputs.
- Verifier: The cited code matches exactly. is_valid_json_number() at include/output_writer.hpp:1053-1071 scans an optional '-' then `while (i < s.size() && std::isdigit(...)) ++i;` with no leading-zero guard. I compiled the function verbatim in isolation: "0123","01","-0123","00","007" ALL return true (NUMBER/bare), exactly as claimed. RFC 8259's number grammar is `int = zero / (digit1-9 *DIGIT)`, which forbids leading zeros, so these are invalid JSON numbers.

The four call sites confirmed:
- Line 1017-10

**[output-writers/correctness] include/output_writer.hpp:1010-1023** — JSON per-transcript custom fields iterate unordered_map -> non-deterministic key order
- What: The 'remaining custom annotations' loop iterates ann.custom_annotations directly, which is a std::unordered_map<std::string,std::string> (vep_annotator.hpp:325). Iteration order over unordered_map is unspecified and varies across runs, libstdc++/libc++ versions, and platforms. The TSV and VCF writers iterate the ordered custom_columns vector and are deterministic, but the JSON writer is not. This breaks output reproducibility and makes byte-for-byte parity testing against Perl VEP (and regression diffing) impossible for any transcript carrying more than one extra custom field.
- Evidence: // include/vep_annotator.hpp:325
    std::unordered_map<std::string, std::string> custom_annotations;
// include/output_writer.hpp:1010
        for (const auto& pair : ann.custom_annotations) {
            if (skip_custom.count(pair.first) || pair.first.substr(0, 11) == "_colocated:") continue;
- Fix: In the JSON writer at output_writer.hpp:1010-1023, do not iterate the unordered_map directly. Collect surviving (lowercased-key, value) pairs into a std::vector<std::pair<std::string,std::string>>, std::sort by key, then emit. Alternatively, iterate the ordered custom_columns list (as the TSV/VCF writers do) and look up each column in ann.custom_annotations, which would also align JSON field ordering with the other two formats. Either approach makes JSON output stable and reproducible across platforms/stdlib versions and across runs with differing source-firing order. Add a regression test that annotates a transcript carrying several custom fields (e.g., MAX_AF + MAX_AF_POP + CCDS + a custom-VCF field) and asserts the emitted JSON key order is fixed.
- Verifier: Verified against current code. include/vep_annotator.hpp:325 declares `std::unordered_map<std::string, std::string> custom_annotations;` exactly as quoted. The JSON writer's "remaining custom annotations" loop at include/output_writer.hpp:1010 iterates this map directly:

    for (const auto& pair : ann.custom_annotations) {
        if (skip_custom.count(pair.first) || pair.first.substr(0, 11) == "_colocated:") continue;
        if (!pair.second.empty()) {
            std::string key = pair.firs

**[perf/performance] src/vep_annotator.cpp:3158-3202** — build_cds_sequence cache returns the full CDS string by value on every call (per-transcript, per-variant copy)
- What: build_cds_sequence is the hot accessor for coding annotation. On a cache hit it returns it->second BY VALUE (a full heap copy of the CDS string). For a coding SNV it is invoked at least twice per transcript per variant: once inside determine_consequences (line 3007) and once in annotate_transcript (line 1995, stored into 'const std::string cached_cds'). CDS sequences can be many kilobases (TTN > 100kb), so each variant pays multiple full-length string allocations+copies. The mutex (cds_cache_mutex_) guards nothing in practice: under --fork each thread uses a separate VEPAnnotator (main.cpp get_annotator returns thread_annotators[t-1]), so cds_cache_ is per-thread and never shared. The lock + copy are pure overhead.
- Evidence: std::lock_guard<std::mutex> lock(cds_cache_mutex_);
        auto it = cds_cache_.find(transcript.id);
        if (it != cds_cache_.end()) {
            return it->second;   // full string copy on every hit
        }
    ...
    const std::string cached_cds = build_cds_sequence(chrom, transcript);  // caller copies again
- Fix: Change build_cds_sequence to return `const std::string&` (a reference into cds_cache_) instead of `std::string` by value. The cache entry outlives every per-variant call within a single annotate() invocation, so the reference is safe. Update the three hot callers to bind references: line 1995 `const std::string& cached_cds = ...`, line 3007 `const std::string& cds_seq = ...`, and the wrapper calls at 3211/3310 (the get_affected_codons overload already takes `const std::string&`, so it can forward the reference directly). On the build/miss path, build into a local, insert into the map, then return `cds_cache_[transcript.id]` (or the iterator from insert) by reference. Since each annotator's cache is per-thread, the per-call mutex can be dropped, or replaced with a one-time build guard if any future shared-annotator path is anticipated. This eliminates 2-3 full multi-kb string copies per coding variant.
- Verifier: All four sub-claims are verified against the actual current code.

1) Value-return on cache hit (src/vep_annotator.cpp:3158-3168): the signature is `std::string VEPAnnotator::build_cds_sequence(...) const`, and the cache hit path is:
    std::lock_guard<std::mutex> lock(cds_cache_mutex_);
    auto it = cds_cache_.find(transcript.id);
    if (it != cds_cache_.end()) { return it->second; }
Returning `it->second` (a member of the map) is an unavoidable full-string copy — copy elision/RVO cannot app

**[perf/performance] src/file_parsers.cpp:235-244** — TabixTSVReader::query builds a std::map<string,string> per record with all columns, allocated per variant
- What: For every overlapping record returned by a tabix query, the reader builds a std::map<std::string,std::string> keyed by column NAME, copying every column name and value into red-black-tree nodes. dbNSFP files have 300+ columns, yet callers (src/sources/dbnsfp_source.cpp line 96-137) only read a handful of fields ('ref','alt','Ensembl_transcriptid', a few scores) via .find(). query() is called once per variant per enabled tabix source (dbnsfp_source.cpp:93, dbscsnv_source.cpp:66, spliceai_source.cpp:139), so each coding variant pays ~300 tree-node allocations + ~300 string copies of column names that are then mostly discarded. The column names are identical across all records and could be interned once.
- Evidence: auto fields = split_line(std::string(str.s, str.l), '\t');

            std::map<std::string, std::string> row;
            for (size_t i = 0; i < fields.size() && i < pimpl_->columns.size(); ++i) {
                row[pimpl_->columns[i]] = fields[i];
            }

            results.push_back(std::move(row));
- Fix: Avoid building a name-keyed std::map per record. Lowest-churn option: change the per-record representation to positional std::vector<std::string> and expose the already-existing column_index (or add a TabixTSVReader::index_of(name)/get(record, name) helper), so callers resolve a name to an index once (at initialize time) and then index into the vector. dbNSFP already computes field_to_column_ indices in initialize(), so it can switch directly to positional access. This reduces per-variant cost from O(columns) tree-node allocations + name copies down to one split_line vector. If changing the public return type across the 6 call sites is undesirable, an interim improvement is to pass the restricted requested-column list into the TabixTSVReader constructor so pimpl_->columns (and thus the per-record map) only contains the few needed fields rather than all 300+.
- Verifier: The cited code at src/file_parsers.cpp:235-244 matches the claim exactly:

    while (tbx_itr_next(pimpl_->fp, pimpl_->tbx, itr, &str) >= 0) {
        auto fields = split_line(std::string(str.s, str.l), '\t');
        std::map<std::string, std::string> row;
        for (size_t i = 0; i < fields.size() && i < pimpl_->columns.size(); ++i) {
            row[pimpl_->columns[i]] = fields[i];
        }
        results.push_back(std::move(row));
    }

For every overlapping record, a std::map<string,st

**[perf/performance] include/output_writer.hpp:1024-1035** — JSON write_transcript_consequence does json.str() copy of the entire growing buffer per transcript (O(n^2))
- What: write_transcript_consequence receives the shared variant-level ostringstream 'json' by reference and is called once per overlapping transcript (loop at line 683-688). To strip a trailing comma it copies the ENTIRE accumulated buffer out with json.str() (full string copy), resizes by 2 bytes, then moves it back with json.str(std::move(s)) and re-seeks to end. Because the buffer grows with each transcript, total work is O(n^2) in emitted bytes for a variant overlapping n transcripts. Genes/regions with many isoforms (dozens of transcripts) make this quadratic. This is on the JSON output hot path (default per-transcript emission).
- Evidence: std::string s = json.str();
            if (s.size() >= 2 && s.back() == '\n' && s[s.size()-2] == ',') {
                s.resize(s.size() - 2);
                s += '\n';
            }
            json.str(std::move(s));
            json.seekp(0, std::ios_base::end);
- Fix: Avoid copying the whole accumulated buffer per transcript. Two clean options: (a) Build each transcript_consequence object into a LOCAL std::ostringstream/std::string inside write_transcript_consequence, strip its own trailing comma locally (O(object size)), then append the finished object once to the parent `json` buffer — drop the json.str()/str(std::move)/seekp dance entirely. (b) Use comma-BEFORE-field emission: track whether any field has been written and prepend ",\n" before each field after the first, so no trailing comma is ever produced and no fixup is needed. Either removes the per-transcript full-buffer copy and the O(n²) blowup, reducing it to O(total bytes). Option (a) is the smallest, most localized change since the function already takes `json` by reference.
- Verifier: Verified against the actual current code in /Users/davidcho/Projects/VEP_ensemble_C/include/output_writer.hpp.

(1) The buffer IS shared and growing. At line 572 a single `std::ostringstream json;` is declared and accumulates the ENTIRE variant-level JSON object: input/id/allele_string/colocated_variants header (lines 573-656), then the transcript_consequences array (lines 680-692), then regulatory/intergenic arrays. There is no per-transcript reset.

(2) `write_transcript_consequence(std::ostri

**[sources-a/parity-gap] src/sources/maxentscan_source.cpp:30-87, 295-328** — MaxEntScan uses a position-independent product model, not the real Yeo & Burge MaxEnt model
- What: The real MaxEntScan (Yeo & Burge) 5'ss score is a maximum-entropy model with first-order dependencies, and the 3'ss is a more complex model over 23 positions; scores are not a sum of per-position log-odds from a position-independent PWM. This implementation computes score = sum_i log2(P[i][base]/0.25) using a single position-independent matrix, and the 3' matrix (MES_3SS_SCORES) uses hand-written generic consensus probabilities with positions -20..-4 all set to a uniform polypyrimidine guess (0.15/0.35/0.10/0.40). The emitted maxentscan:5ss_ref/alt/diff and 3ss_* values therefore do not match Perl VEP's MaxEntScan plugin numerically, so any downstream threshold (e.g. interpret_diff -3/-1.5) is applied to non-comparable scores. This is acknowledged in the project's own gap notes ('3' PWM is ~80% simplified (not real Yeo & Burge model)', '5' PWM is ~50% approximate').
- Evidence: // Positions -20 to -4 (polypyrimidine tract region)
{{ 0.15, 0.30, 0.15, 0.40 }},  // -20
... (all -16..-4 identical {{0.15,0.35,0.10,0.40}})
...
double prob = MES_5SS_SCORES[i][base_idx];
if (prob <= 0) return -999.0;
score += std::log2(prob / 0.25);  // Log-odds vs uniform
- Fix: Port the real Yeo & Burge MaxEntScan score5/score3 algorithm using the published probability/score tables (me2x5, me2x3acc and the splice site consensus matrices), which implement the maximum-entropy models with dependency structure, instead of the position-independent PWM product. The original Perl matrices can be loaded at startup the same way Perl VEP's plugin does. If exact parity is not feasible, stop labeling these outputs with the maxentscan: prefix as if comparable to Perl VEP — rename them or clearly document in --help/README that they are an internal approximation, and remove the unused interpret_diff thresholds (or recalibrate them to the approximate scale) so future code does not apply Perl-VEP-calibrated cutoffs to these non-comparable values. Additionally, add at least one test that scores a known splice-site sequence and compares against the published MaxEntScan reference value to lock in parity once implemented.
- Verifier: The cited code matches the claim exactly. In /Users/davidcho/Projects/VEP_ensemble_C/src/sources/maxentscan_source.cpp:

5'ss: MES_5SS_SCORES (lines 30-49) is a single 9-position position-weight matrix, and score_5ss (lines 295-308) scores it as a position-INDEPENDENT product of per-position log-odds: `score += std::log2(prob / 0.25);` over each of 9 positions independently. The real Yeo & Burge MaxEnt 5'ss model is a maximum-entropy distribution capturing adjacent and non-adjacent first-order d

**[sources-a/parity-gap] src/sources/dbscsnv_source.cpp:115-123, 195-198** — dbscSNV adds a 'possible'/0.4 tier not present in Perl VEP and drops scores <= 0
- What: Perl VEP's dbscSNV usage reports the raw ada_score and rf_score and treats >=0.6 as the splice-altering cutoff; there is no standard 'possible' (>=0.4) interpretation tier. interpret_score invents a 'possible' category that has no Perl counterpart, so the dbscsnv:prediction field is non-parity output. Separately, max_score is only emitted when max_score > 0; ada_score/rf_score are probabilities in [0,1] and a legitimate near-zero score (e.g. 0.0, or values just above 0) is silently suppressed, and an exactly-0 valid score is treated identically to missing. The dbscsnv:max_score / prediction fields can therefore differ from a reference and omit valid low scores.
- Evidence: double max_score = std::max(ada, rf);
if (max_score > 0) { ... interpret_score(max_score) ... }
...
static std::string interpret_score(double score) {
    if (score >= 0.6) return "splice_altering";
    if (score >= 0.4) return "possible";   // no Perl VEP equivalent
    return "";
}
- Fix: For strict parity, remove the dbscsnv:max_score and dbscsnv:prediction fields entirely (and the interpret_score helper), emitting only the raw ada_score and rf_score as the Perl VEP plugin does; update get_fields() and the doc comment accordingly. If a prediction tier is desired as an explicit extension, drop the non-standard 0.4 "possible" tier and use only the documented 0.6 cutoff reporting a single PASS/FAIL (or splice_altering/"") classification, and clearly label it as a non-Perl extension. Separately, if max_score is kept, guard its emission on the presence of at least one parsed score (e.g. a parsed-success flag) rather than on `max_score > 0`, so a legitimate 0.0 is not conflated with missing data.
- Verifier: The cited code matches the claim exactly. At src/sources/dbscsnv_source.cpp:195-199:

  static std::string interpret_score(double score) {
      if (score >= 0.6) return "splice_altering";
      if (score >= 0.4) return "possible";
      return "";
  }

And at lines 115-123:

  double max_score = std::max(ada, rf);
  if (max_score > 0) {
      char buf[32]; std::snprintf(buf, sizeof(buf), "%.4f", max_score);
      annotations["dbscsnv:max_score"] = buf;
      std::string pred = interpret_score(m

**[sources-b/parity-gap] src/sources/lof_source.cpp:153-157** — LOFTEE INCOMPLETE_CDS flag is effectively dead code; uses genomic CDS coordinates instead of cds_start_NF/cds_end_NF
- What: The INCOMPLETE_CDS flag tests transcript->cds_start <= 0 || transcript->cds_end <= 0. But cds_start/cds_end hold genomic coordinates of the CDS boundaries (set during GTF load in vep_annotator.cpp:1448-1463), not a 'start found' boolean. The LOFTEE source already returns early unless transcript->is_coding() (which is !cds_regions.empty()), so any transcript that reaches this code has CDS regions and therefore cds_start/cds_end set to real (positive) genomic positions. As a result this branch can essentially never fire, and genuinely incomplete-CDS transcripts are never flagged. The actual incomplete-CDS indicators are the cds_start_NF / cds_end_NF fields (Transcript struct, vep_annotator.hpp:227-228; set at vep_annotator.cpp:1369-1372), which this code never consults. Perl VEP LOFTEE downgrades transcripts with incomplete CDS, so this is a missing-output parity gap.
- Evidence: // Flag: Incomplete CDS
if (transcript->cds_start <= 0 || transcript->cds_end <= 0) {
    flags.push_back("INCOMPLETE_CDS");
    is_hc = false;
}
- Fix: Replace the genomic-coordinate test at lof_source.cpp:154 with the actual incomplete-CDS booleans: `if (transcript->cds_start_NF || transcript->cds_end_NF) { flags.push_back("INCOMPLETE_CDS"); is_hc = false; }`. The existing `is_coding()` guard at line 69 already handles the empty-cds_regions case. Add a regression test that loads a transcript with cds_start_NF (or cds_end_NF) set and a LoF variant, then asserts the INCOMPLETE_CDS flag is emitted and is_hc is downgraded to LC.
- Verifier: Confirmed against actual code. At src/sources/lof_source.cpp:153-157 the flag reads exactly: `// Flag: Incomplete CDS\nif (transcript->cds_start <= 0 || transcript->cds_end <= 0) {\n    flags.push_back("INCOMPLETE_CDS");\n    is_hc = false;\n}`.

The branch is effectively dead because:
1. lof_source.cpp:69 returns early unless `transcript->is_coding()`, and `is_coding()` is defined at vep_annotator.hpp:238 as `return !cds_regions.empty();`. So any transcript reaching line 154 has at least one CD

**[structural-variants/correctness] include/structural_variant.hpp:228-244, 258-273** — Splice donor/acceptor assignment ignores transcript strand (always plus-strand semantics)
- What: For DEL/INS/DUP SVs, the code labels the 2bp before an exon start (exon.start-2..exon.start-1) as the splice ACCEPTOR and the 2bp after an exon end (exon.end+1..exon.end+2) as the splice DONOR, with no strand awareness. This is only correct for plus-strand transcripts. On a minus-strand transcript the genomic upstream boundary (exon.start side) is the DONOR and the genomic downstream boundary (exon.end side) is the ACCEPTOR. The main small-variant annotator does this correctly with `bool is_minus = (transcript.strand == '-')` (src/vep_annotator.cpp:2727-2789), so SVs on minus-strand genes will report splice_donor where Perl VEP reports splice_acceptor and vice versa. Roughly half of all genes are on the minus strand, so this systematically swaps donor/acceptor consequences for SVs hitting minus-strand splice sites — a clinically meaningful mislabeling.
- Evidence: for (size_t i = 0; i < transcript.exons.size(); ++i) {
    const auto& exon = transcript.exons[i];
    // Acceptor site: 2bp before exon start (exon_start-2 to exon_start-1)
    if (i > 0 && sv.overlaps(exon.start - 2, exon.start - 1)) {
        hits_acceptor = true;
    }
    // Donor site: 2bp after exon end (exon_end+1 to exon_end+2)
    if (i + 1 < transcript.exons.size() && sv.overlaps(exon.end + 1, exon.end + 2)) {
        hits_donor = true;
    }
}
if (hits_donor) {
    consequences.push_back(ConsequenceType::SPLICE_DONOR_VARIANT);
}
if (hits_acceptor) {
    consequences.push_back(ConsequenceType::SPLICE_ACCEPTOR_VARIANT);
}
- Fix: In include/structural_variant.hpp, make donor/acceptor assignment strand-aware in both the DEL block (228-244) and the INS/DUP/TDUP block (258-273), mirroring src/vep_annotator.cpp:2727-2789. Compute `bool is_minus = (transcript.strand == '-');` once. The exon.start-side dinucleotide (exon.start-2..exon.start-1) is the acceptor on plus strand but the donor on minus strand; the exon.end-side dinucleotide (exon.end+1..exon.end+2) is the donor on plus strand but the acceptor on minus strand. Simplest implementation: keep the overlap tests for the two genomic boundaries but assign the consequence based on is_minus, e.g. set hits at exon.start-side -> (is_minus ? donor : acceptor) and exon.end-side -> (is_minus ? acceptor : donor). Add a regression test with a minus-strand transcript and an SV deleting/inserting at each canonical splice dinucleotide, asserting the swapped labels match the small-variant annotator and Perl VEP.
- Verifier: Confirmed against the actual code. In include/structural_variant.hpp the DEL block (lines 228-244) and the INS/DUP/TDUP block (lines 258-273) assign splice labels with no strand awareness:

  for (size_t i = 0; i < transcript.exons.size(); ++i) {
      const auto& exon = transcript.exons[i];
      if (i > 0 && sv.overlaps(exon.start - 2, exon.start - 1)) { hits_acceptor = true; }   // labeled ACCEPTOR
      if (i + 1 < transcript.exons.size() && sv.overlaps(exon.end + 1, exon.end + 2)) { hits_do

**[structural-variants/correctness] include/structural_variant.hpp:246-254, 278-289** — Frameshift/inframe determination for DEL/INS/DUP uses total SV span, not bases inside CDS
- What: When a DEL (or INS/DUP) affects the CDS, frameshift vs inframe is decided by `sv.length() % 3`, where length() is the entire SV span (end-start+1 or |sv_len|). But an SV usually overlaps only part of the CDS (it can span introns, UTRs, and intergenic flanks). The number of *coding* bases removed/added is what determines reading-frame disruption, not the total SV length. E.g. a deletion spanning 1000bp of which only 4 fall in a coding exon: 1000%3==1 -> frameshift_variant, but the true effect on the protein depends on the 4 coding bases (and on intron content). Conversely a 1002bp deletion (1002%3==0 -> inframe_deletion) that actually removes a non-multiple-of-3 number of coding bases is mislabeled inframe. For intron-spanning SVs the intronic bases must not be counted at all. Perl VEP computes the inframe/frameshift status from the CDS-coordinate effect, so this diverges and can flip the IMPACT (HIGH frameshift vs MODERATE inframe).
- Evidence: if (affects_cds) {
    // Check if frameshift
    int deleted_bp = sv.length();
    if (deleted_bp % 3 != 0) {
        consequences.push_back(ConsequenceType::FRAMESHIFT_VARIANT);
    } else {
        consequences.push_back(ConsequenceType::INFRAME_DELETION);
    }
}
- Fix: In get_sv_consequences, replace `sv.length()` in the affects_cds blocks (lines 248 and 283) with the count of coding bases the SV overlaps: int coding_bp = 0; for (const auto& cds : transcript.cds_regions) coding_bp += std::max(0, std::min(sv.end, cds.end) - std::max(sv.start, cds.start) + 1); then use `coding_bp % 3` for the DEL frameshift/inframe decision. For DUP/TDUP/INS, base the inserted-coding-base count on the duplicated/inserted copy restricted to CDS (for a DUP the inserted copy length within CDS equals coding_bp; for INS with known sv_len it is the inserted length, which only disrupts frame if it lands inside a coding exon). Guard against coding_bp == 0 (SV touches CDS region edge but adds/removes no coding bases) by falling back to CODING_SEQUENCE_VARIANT rather than emitting a spurious inframe call. Add tests with multi-region cds_regions where the SV spans an intron so the full span and coding-base count differ in their mod-3 residue.
- Verifier: Verified the cited code in include/structural_variant.hpp. The DEL path (lines 246-254) does `int deleted_bp = sv.length(); if (deleted_bp % 3 != 0) FRAMESHIFT_VARIANT else INFRAME_DELETION;` and the INS/DUP/TDUP path (lines 278-289) does the same with `int inserted_bp = sv.length();`. The `length()` method (lines 99-102) returns `abs(sv_len)` else `end - start + 1`, i.e. the ENTIRE genomic span of the SV — and parse_sv_from_vcf (lines 384-414) sets sv_len/end to the full span including introns/

**[structural-variants/correctness] include/structural_variant.hpp:411-414** — Insertion length for symbolic <INS> with explicit END is computed as a span, inflating inserted length
- What: For an INS with no SVLEN but an explicit END, sv_len is set to `end - start + 1`. For a true insertion the END normally equals POS (insertions do not consume reference bases), so this branch only triggers when END>POS, but it then treats the reference span as the inserted length. The inserted length of an INS is the length of novel sequence (SVLEN / ALT length), not the reference interval. Downstream (lines 283-288) this wrong inserted_bp feeds the frameshift/inframe %3 test, so the inframe-vs-frameshift call for symbolic insertions with END can be wrong. Insertions should derive length from SVLEN or the explicit inserted sequence, not from end-start.
- Evidence: } else if (sv.sv_type == SVType::INS && sv.end != sv.start) {
    // Only calculate INS length when END was explicitly provided
    sv.sv_len = sv.end - sv.start + 1;
}
- Fix: Drop the `else if (sv.sv_type == SVType::INS && sv.end != sv.start)` branch (lines 411-413) entirely. For symbolic INS, derive length only from SVLEN (already handled at lines 392-404) or, if available, from a literal inserted ALT sequence. When neither is present, leave `sv_len = 0` (unknown), which already routes to CODING_SEQUENCE_VARIANT at line 279 instead of guessing frameshift/inframe from the reference span. Optionally add a regression test: parse_sv_from_vcf("chr1", 1000, "N", "<INS>", {{"SVTYPE","INS"},{"END","1100"}}) and assert sv.sv_len == 0.
- Verifier: The cited code matches exactly. In include/structural_variant.hpp, within the `else` branch (no SVLEN in INFO, lines 405-416):

```cpp
} else if (sv.sv_type == SVType::INS && sv.end != sv.start) {
    // Only calculate INS length when END was explicitly provided
    sv.sv_len = sv.end - sv.start + 1;
}
```

Tracing how `sv.end` is set: for an INS without explicit END, line 387 forces `sv.end = pos` (== start). END is only != start when it was explicitly provided in INFO (lines 375-377). So this 


### LOW

**[build-quality/correctness] include/vep_annotator.hpp:191-193,201-202,214-216,250-252** — Uninitialized scalar members in Exon/CDS/Transcript/Gene structs
- What: Several core structs leave scalar members without in-class default initializers: Exon::start, Exon::end, Exon::exon_number (lines 191-193); CDS::start, CDS::end (201-202); Transcript::start, Transcript::end, Transcript::strand (214-216); Gene::start, Gene::end, Gene::strand (250-252). The current GTF loader happens to set Exon/CDS/Transcript fields before use, so the bug is latent, but any default-constructed instance (a test, a future code path, a transcript built without a full GTF record, or strand left unset when the GTF strand column is malformed) reads indeterminate values -> UB and nondeterministic coordinate math. strand in particular is only assigned conditionally and a bad GTF strand field can leave it garbage. Apple clang (the dev compiler here) does not emit -Wmaybe-uninitialized, so this will not surface during normal builds.
- Evidence: struct Exon {
    int start;          // 1-based, inclusive
    int end;            // 1-based, inclusive
    int exon_number;
    int phase = -1;
};
struct CDS {
    int start;
    int end;
    int phase = 0;
};
struct Transcript {
    ...
    int start;
    int end;
    char strand;
- Fix: Add in-class default initializers to the scalar members in include/vep_annotator.hpp to remove the latent UB and make default construction deterministic, matching the style already used for the other members (phase, cds_start, tsl): Exon { int start = 0; int end = 0; int exon_number = 0; int phase = -1; }; CDS { int start = 0; int end = 0; int phase = 0; }; Transcript { ... int start = 0; int end = 0; char strand = '+'; ... }; Gene { ... int start = 0; int end = 0; char strand = '+'; ... }. This is a no-op for all current code paths (which already assign these fields before use) but hardens future/test code. Optionally, in the GTF loader, validate that strand_str[0] is '+' or '-' (e.g., default unknown '.' to '+') to avoid a defined-but-nonsensical strand value propagating into coordinate math.
- Verifier: The cited code matches the claim exactly in /Users/davidcho/Projects/VEP_ensemble_C/include/vep_annotator.hpp:

  struct Exon { int start; int end; int exon_number; int phase = -1; };  // lines 190-195
  struct CDS { int start; int end; int phase = 0; };                     // lines 200-204
  struct Transcript { ... int start; int end; char strand; ... };        // lines 209-216 (start/end/strand uninitialized; cds_start/cds_end/tsl etc. DO have initializers)
  struct Gene { ... int start; int e

**[core-consequence/parity-gap] src/vep_annotator.cpp:2992-3002** — Inframe delins with multi-base ref and alt classified as protein_altering_variant instead of inframe_insertion/deletion
- What: For an in-frame length-changing variant (length_diff % 3 == 0), the code emits PROTEIN_ALTERING_VARIANT whenever both ref and alt have length > 1 and ref_len != alt_len (a delins). For the common case of a net inframe deletion delins (e.g. ref 6bp, alt 3bp) or net inframe insertion delins, Perl VEP classifies the predominant effect as inframe_deletion / inframe_insertion respectively; it reserves protein_altering_variant for cases where the precise inframe insertion/deletion cannot be resolved. This produces a less-specific (sometimes incorrect) term for clean net inframe delins. Note: Perl does emit protein_altering_variant in some genuinely indeterminate inframe cases, so this is a partial parity gap rather than a universal error.
- Evidence: if (length_diff != 0 && length_diff % 3 == 0) {
    if (ref_len > 1 && alt_len > 1 && ref_len != alt_len) {
        // Complex inframe change (delins)
        consequences.push_back(ConsequenceType::PROTEIN_ALTERING_VARIANT);
    } else if (length_diff > 0) {
        consequences.push_back(ConsequenceType::INFRAME_INSERTION);
    } else {
        consequences.push_back(ConsequenceType::INFRAME_DELETION);
    }
    return consequences;
}
- Fix: Replace the unconditional PROTEIN_ALTERING_VARIANT for multi-base in-frame delins with classification by net length change, matching Perl VEP. Concretely, in the inframe branch (lines 2992-3002), drop the `ref_len > 1 && alt_len > 1 && ref_len != alt_len` special case and classify by sign of length_diff: length_diff > 0 -> INFRAME_INSERTION, length_diff < 0 -> INFRAME_DELETION. Reserve PROTEIN_ALTERING_VARIANT for cases where the inframe change cannot be cleanly resolved to insertion or deletion (e.g. the affected codons cannot be reconstructed from the cached CDS sequence, or the variant spans an exon/intron boundary), mirroring Perl's BaseVariationFeatureOverlapAllele fallback. Add tests in tests/test_consequences.cpp for a net-loss delins (ref 6bp/alt 3bp -> inframe_deletion) and net-gain delins (ref 3bp/alt 6bp -> inframe_insertion) to lock in the behavior.
- Verifier: The cited code at src/vep_annotator.cpp:2991-3002 matches the reviewer's quote exactly:

    // Inframe indel
    if (length_diff != 0 && length_diff % 3 == 0) {
        if (ref_len > 1 && alt_len > 1 && ref_len != alt_len) {
            // Complex inframe change (delins)
            consequences.push_back(ConsequenceType::PROTEIN_ALTERING_VARIANT);
        } else if (length_diff > 0) {
            consequences.push_back(ConsequenceType::INFRAME_INSERTION);
        } else {
            consequen

**[core-hgvs-gen/parity-gap] src/vep_annotator.cpp:3507-3514** — UTR/intronic/non-coding insertions never reported as duplications (no dup detection in append_hgvs_change)
- What: append_hgvs_change (used by build_hgvsc for all UTR c.-/c.*, intronic +/-, and non-coding n. notations) always emits insertions as 'ins' and has no duplication detection, unlike the coding CDS generate_hgvsc which checks the preceding sequence and emits 'dup'. Per HGVS and Perl VEP, an inserted sequence identical to the immediately preceding reference bases must be described as a duplication (e.g. c.-12dup, n.345_346dup) rather than an insertion. As a result, duplications occurring in 5'/3' UTRs, introns, and non-coding transcripts are emitted as 'ins' notation that does not match Perl VEP output.
- Evidence: } else if (hgvs_alt.length() > hgvs_ref.length()) {
        std::string inserted = a;
        std::string end_pos = offset_fn(1);
        if (!end_pos.empty()) {
            oss << start_pos << "_" << end_pos << "ins" << inserted;
        } else {
            oss << start_pos << "ins" << inserted;
        }
- Fix: Add duplication detection to the insertion branch used by UTR/intronic/non-coding HGVSc. The cleanest approach mirrors the CDS logic but requires the relevant transcript sequence for the comparison: (1) For the exonic non-coding `n.` path (line 2572), build a full cDNA sequence (analogous to build_cds_sequence) and check whether the inserted bases equal the `ins_len` cDNA bases immediately preceding the cDNA position; if so emit `dup` (single base `Ndup`, multi-base `start_enddup`) instead of `ins`. (2) For UTR (`c.-N`/`c.*N`) and intronic (`c.N+/-M`) paths, dup detection cannot reuse a simple cDNA index — it needs the genomic reference and position-aware arithmetic (the offset_fn already encodes how to step the position). Pass the inserted bases and a preceding-base lookup (genomic reference window, strand-complemented for '-' transcripts) into append_hgvs_change (or do the check in build_hgvsc before delegating) so that when the inserted sequence equals the immediately preceding bases of equal length, the formatter emits the duplicated range with `dup`. Add unit tests covering c.-12dup, c.*5_*6dup, intronic c.N+Mdup, and n.345_346dup against expected Perl VEP output. Given the narrow impact, this is reasonable to schedule as a parity-polish item rather than urgent.
- Verifier: The claim accurately describes the current code. The shared formatter `append_hgvs_change` (src/vep_annotator.cpp:3478-3524) handles the insertion branch at lines 3507-3514 exactly as quoted, with no duplication detection — it unconditionally emits `ins`:
```cpp
} else if (hgvs_alt.length() > hgvs_ref.length()) {
    std::string inserted = a;
    std::string end_pos = offset_fn(1);
    if (!end_pos.empty()) {
        oss << start_pos << "_" << end_pos << "ins" << inserted;
    } else {
        o

**[file-parsers/portability] src/file_parsers.cpp:755** — INT_MAX/INT_MIN used without including <climits>
- What: IntervalTree::build_tree uses the C macros INT_MAX and INT_MIN, but neither file_parsers.cpp nor any header it includes pulls in <climits>/<limits.h>. I verified no project header includes <climits>, <limits.h>, or <limits>. It currently compiles only because some other standard-library header transitively defines these macros, which is not guaranteed by the C++ standard. A change in include order, a libstdc++/libc++/MSVC upgrade, or moving this code could turn it into a hard compile error. This is a fragile portability dependency on transitive includes in a hot-path data structure.
- Evidence: int min_start = INT_MAX, max_end = INT_MIN;
    for (size_t idx : indices) {
        min_start = std::min(min_start, intervals_[idx].start);
        max_end = std::max(max_end, intervals_[idx].end);
- Fix: Add `#include <climits>` to the include block at the top of src/file_parsers.cpp (after <cmath>). Alternatively, for full standard-library style consistency with filter_vep.hpp, replace INT_MAX/INT_MIN with std::numeric_limits<int>::max()/std::numeric_limits<int>::min() and add `#include <limits>`. Either change is one line, zero runtime impact, and removes the reliance on transitive includes.
- Verifier: The cited code matches exactly. Line 755 of src/file_parsers.cpp reads: `int min_start = INT_MAX, max_end = INT_MIN;` inside IntervalTree<T>::build_tree, followed by the loop computing std::min/std::max over intervals_[idx].start/.end (which are `int`, per include/file_parsers.hpp lines 335-336). A repo-wide grep (`grep -rn "climits\|limits.h\|<limits>" include/ src/`) returns ZERO matches, confirming neither file_parsers.cpp nor any header it includes pulls in <climits>/<limits.h>/<limits>. The

**[file-parsers/performance] src/file_parsers.cpp:215-247** — Tabix chromosome-format fallback issues two tabix queries per lookup on the hot path
- What: query_range always builds a 2-element chrom_variants list and tries each region in turn, breaking only when results are non-empty. For the common dbNSFP/SpliceAI/dbscSNV layout that uses bare contig names (1, 2, X), the first attempt is 'chr'+chrom which never matches, so every single per-variant query performs a wasted tbx_itr_querys/iterator-create/destroy before the second attempt succeeds (and for positions with no record at all, both attempts run). This doubles tabix index lookups on a per-variant hot path used by multiple annotation sources. Detecting the file's naming convention once (e.g. from tbx seq names at open time) would avoid the redundant query.
- Evidence: std::vector<std::string> chrom_variants;
    if (chrom.substr(0, 3) == "chr") {
        chrom_variants.push_back(chrom);
        chrom_variants.push_back(chrom.substr(3));
    } else {
        chrom_variants.push_back("chr" + chrom);
        chrom_variants.push_back(chrom);
    }

    for (const auto& try_chrom : chrom_variants) {
        std::string region = try_chrom + ":" + ...
- Fix: At TabixTSVReader construction, call tbx_seqnames(pimpl_->tbx, &n) to record whether the index uses 'chr'-prefixed names (store a bool/enum and the seqname set in Impl). In query_range, build only the single region string matching the index's convention (translate the input chrom accordingly), eliminating the speculative first attempt. Alternatively, keep the fallback but cache the first successfully-resolved chrom style per reader and try that format first on subsequent queries, falling back only if it ever misses. Either approach removes the redundant per-variant tbx_itr_querys on the common bare-contig layout.
- Verifier: The cited code matches the claim exactly. In src/file_parsers.cpp query_range (lines 215-247), the reader unconditionally builds a 2-element chrom_variants list and loops over it, breaking only on non-empty results:

```
std::vector<std::string> chrom_variants;
if (chrom.substr(0, 3) == "chr") {
    chrom_variants.push_back(chrom);
    chrom_variants.push_back(chrom.substr(3));
} else {
    chrom_variants.push_back("chr" + chrom);   // tried FIRST for bare input
    chrom_variants.push_back(chro

**[file-parsers/parity-gap] src/file_parsers.cpp:48-62** — GFF3 attribute parser does not trim whitespace around key/value, silently dropping ID/Name/Parent on lenient files
- What: parse_gff3_attributes splits the column-9 string on ';' and '=' and uses the raw substrings as key/value with no trimming. Strictly spec-compliant Ensembl GFF3 has no spaces around '=' or after ';', so this works for the intended input. However, GFF3 files produced by some tools include spaces (e.g. 'ID=foo; Name=bar' or 'ID = foo'), in which case the key becomes ' Name' / 'ID ' and the reserved-attribute lookups feat.attributes.count("ID"/"Name"/"Parent") fail, leaving feat.id/name/parent empty. For the regulatory source that means regulatory:feature_id and feature names are silently omitted from output for such inputs rather than raising an error. Low impact for canonical Ensembl files; a real output/parity gap for non-canonical GFF3.
- Evidence: auto fields = split_line(attrs, ';');
    for (const auto& field : fields) {
        size_t eq = field.find('=');
        if (eq != std::string::npos) {
            std::string key = field.substr(0, eq);
            std::string value = url_decode(field.substr(eq + 1));
            result[key] = value;
        }
    }
- Fix: In parse_gff3_attributes (src/file_parsers.cpp:48-62), trim leading/trailing whitespace (' ', '\t', and also '\r' to be safe against CRLF line endings) from each field before splitting on '=', and from the key and value substrings before insertion. A small helper, e.g. auto trim = [](std::string s){ size_t b = s.find_first_not_of(" \t\r"); if (b==std::string::npos) return std::string(); size_t e = s.find_last_not_of(" \t\r"); return s.substr(b, e-b+1); }; then key = trim(field.substr(0, eq)) and value = url_decode(trim(field.substr(eq+1))). This is low-risk, leaves canonical Ensembl GFF3 unchanged, and restores regulatory:feature_id/name output for lenient GFF3. Optionally add a unit test feeding "ID = foo; Name = bar" and asserting id/name are populated.
- Verifier: The cited code matches the claim exactly. src/file_parsers.cpp:48-62, parse_gff3_attributes:

    auto fields = split_line(attrs, ';');
    for (const auto& field : fields) {
        size_t eq = field.find('=');
        if (eq != std::string::npos) {
            std::string key = field.substr(0, eq);
            std::string value = url_decode(field.substr(eq + 1));
            result[key] = value;
        }
    }

No trimming is applied to either key or value. The consumption chain is confirmed 

**[filter-transcript/correctness] include/filter_vep.hpp:203-226** — Numeric gt/ge/lt/le comparisons 'pass through' (return true) for non-numeric values
- What: When a field value is non-numeric (e.g. '.', empty, or text), the GREATER/GREATER_EQ/LESS/LESS_EQ branches unconditionally set result = true (commented 'Pass through ... match Perl VEP'). This means a filter like 'CADD_phred > 20' KEEPS every record whose CADD_phred is '.' or missing, rather than dropping it. For a clinical 'show me variants with CADD > 20' filter this silently includes unscored variants, the opposite of the user's intent. Perl VEP's filter_vep instead treats a value that is undef/non-numeric under a numeric operator as a non-match (it skips/excludes), so this is also a parity gap. Note the inconsistency: EQUALS on a non-numeric value correctly does a string compare, but the ordering operators do not.
- Evidence: } else if (cond.op == FilterOperator::GREATER) { if (is_numeric) { result = (num_value > num_target); } else { result = true; // Pass through when value is non-numeric (match Perl VEP) } }
- Fix: Do not change behavior blindly. First verify against real Perl VEP filter_vep: run `filter_vep -filter "CADD_phred > 20"` on a record whose CADD_phred is "." and on one where the field is absent, and observe whether the record is kept or dropped. If Perl drops it (as the reviewer believes), change the four non-numeric branches (lines 207, 213, 219, 225) to `result = false;` so that missing/non-numeric values fail an ordering threshold. If Perl keeps it (as the current comment claims), leave the code but document the behavior explicitly in the README filter_vep section so clinical users know that `> threshold` retains unscored variants and must pair it with an EXISTS check (e.g. `--filter "CADD_phred and CADD_phred > 20"`). Either way, replace the bare "(match Perl VEP)" comment with the actual verified evidence (a link/citation to the observed Perl output) so the rationale is auditable.
- Verifier: The cited code matches the reviewer's quote exactly. include/filter_vep.hpp lines 203-226 contain four ordering branches that each do, for a non-numeric value: `result = true; // Pass through when value is non-numeric (match Perl VEP)`. The preceding numeric-detection block (lines 176-189) marks a value as non-numeric when it is empty, ".", "NA", or any text that fails std::stod, and the EXISTS/NOT_EXISTS operators are handled separately (lines 166-174) so they do not short-circuit the ordering 

**[filter-transcript/correctness] include/filter_vep.hpp:307-340** — --max-af / --min-af / --min-cadd / --min-revel keep records when the score field is absent
- What: Each numeric quick-filter only applies its threshold when the looked-up value is non-NaN. If the AF/CADD/REVEL field is missing or '.', the record passes the filter unchanged. So '--max-af 0.001' retains all variants with no recorded frequency, and '--min-cadd 20' retains all variants with no CADD score. Whether that matches Perl VEP depends on the field, but for --max-af in particular this is dangerous: a rare-variant filter silently keeps every variant lacking a frequency annotation, which may be the majority. At minimum this divergence from a strict threshold is undocumented and clinically risky.
- Evidence: if (!std::isnan(af)) { if (config.min_af >= 0 && af < config.min_af) return false; if (config.max_af >= 0 && af > config.max_af) return false; }  // af==NaN => no check, record passes
- Fix: Make missing-value semantics explicit and document them. Concretely: (1) Document in `filter_vep --help` (src/filter_vep.cpp around lines 39-42) and README that records lacking AF/CADD/REVEL annotation pass the numeric threshold filters unchanged. (2) Reconsider per-filter semantics: keeping "missing AF = novel = pass" for --max-af is reasonable, but for score-threshold filters (--min-cadd, --min-revel, and arguably --min-af) consider an opt-in flag (e.g. --exclude-missing) to reject records lacking the score, since silently retaining unscored variants under a pathogenicity/frequency floor is clinically surprising and diverges from Perl VEP's FilterSet behavior, which drops records whose comparison field is undefined. (3) If strict Perl-VEP parity for comparison operators is the goal, correct the test comments at tests/test_filter_vep.cpp:674/885 and apply_condition lines 203-226, since Perl VEP's FilterSet does not pass undefined fields through `<`/`>` comparisons.
- Verifier: The cited code in include/filter_vep.hpp matches the quote exactly. At lines 307-340:

  if (config.min_af >= 0 || config.max_af >= 0) {
      double af = record.get_numeric("AF"); ... fallbacks to gnomAD_AF, gnomAD:AF ...
      if (!std::isnan(af)) {
          if (config.min_af >= 0 && af < config.min_af) return false;
          if (config.max_af >= 0 && af > config.max_af) return false;
      }
  }
  if (config.min_cadd >= 0) { double cadd = ...; if (!std::isnan(cadd) && cadd < config.min_cadd

**[filter-transcript/parity-gap] include/filter_vep.hpp:469-473** — NOT_CONTAINS and NOT_IN operators cannot be produced by the expression parser; 'not' prefix only negates whole conditions
- What: FilterOperator::NOT_CONTAINS and NOT_IN are defined and implemented but parse_filter_operator never returns them and parse_filter_expression never produces them; they are only constructible by manual struct manipulation (as the tests do). The only negation route is a leading 'not ' on the field, which is stripped from cond.field via substr(0,4)=="not ". This means a natural Perl-VEP-style expression like 'Feature not in ENST1,ENST2' parses with the " in " operator first, leaving field = "Feature not" (the 'not ' prefix check looks at the START of the field, sees 'Feat', and does NOT negate), so the negation is silently dropped and the dangling 'not' becomes part of the field name (which then never matches a real column). Result: 'X not in ...' / 'X not contains ...' filters behave as a no-op or wrong match rather than a negated test.
- Evidence: if (cond.field.size() > 4 && cond.field.substr(0, 4) == "not ") { cond.field = cond.field.substr(4); cond.negated = true; }  // only matches 'not ' at field START, not trailing 'X not'
- Fix: Two reasonable options. Minimal/cleanest: handle "not " adjacent to the operator in parse_filter_expression. After extracting and trimming cond.field, also check for a trailing " not" token (e.g. if field ends with " not", strip it and set cond.negated = true), so "Feature not in [...]" and "X not contains Y" negate correctly. Alternatively, since NOT_CONTAINS/NOT_IN are already fully implemented in apply_condition, map a stripped-not + CONTAINS/IN to those operators (or just rely on the negated flag, which already inverts the result at line 256). Either way, also add a regex-free unit test for trailing-not parsing. If parity with Perl VEP is the only concern, note that leading-not already works, so this is optional; at minimum, document that only leading "not" is supported and consider removing the unreachable NOT_CONTAINS/NOT_IN enum values to avoid dead-code confusion.
- Verifier: The claim is factually verified against the actual code in include/filter_vep.hpp.

1) The enum defines NOT_CONTAINS (line 38) and NOT_IN (line 40), and apply_condition implements them at lines 229-230 (NOT_CONTAINS) and 238-245 (NOT_IN). But parse_filter_operator (lines 49-67) never returns either: it only maps "contains"/"match" -> CONTAINS and "in" -> IN. parse_filter_expression (lines 418-518) calls parse_filter_operator at line 485, so it can never produce NOT_CONTAINS or NOT_IN. Confirmed 

**[filter-transcript/parity-gap] include/transcript_filter.hpp:540-577** — per_gene / pick_allele output order is alphabetized by std::map key, not input/Perl order
- What: pick_per_allele, pick_per_gene, and pick_per_allele_gene group annotations into a std::map<std::string, ...> keyed by allele/gene string and then emit results by iterating the map. std::map iterates in sorted key order, so the surviving annotations are emitted sorted alphabetically by allele or gene_id/symbol rather than in the order Perl VEP emits them (order of first appearance / transcript processing order). This changes per-variant line ordering relative to Perl VEP, which breaks line-by-line output diffing/parity even when the selected set is identical.
- Evidence: std::map<std::string, std::vector<AnnotationWithMeta>> by_gene; for (auto& ann : annotations) { ... by_gene[gene_key].push_back(std::move(ann)); } ... for (auto& pair : by_gene) { auto picked = pick_one(std::move(pair.second)); ... }
- Fix: Preserve first-appearance (input/transcript-processing) order instead of std::map's sorted key order. Replace the std::map grouping with insertion-order tracking: keep a std::unordered_map<std::string, size_t> mapping key to a slot index, and a parallel std::vector<std::vector<AnnotationWithMeta>> groups that preserves first-seen order; push new keys to the back. Then iterate groups in vector order, calling pick_one on each and appending the winner to result. Apply this identically to pick_per_allele (lines 535-555), pick_per_gene (557-577), and pick_per_allele_gene (643-666). The flag_picked_per_* variants already preserve order and need no change. Add a test asserting that, given input annotations whose gene_ids/alleles are not in alphabetical order, the filtered output preserves first-appearance order.
- Verifier: The cited code matches the claim exactly. In include/transcript_filter.hpp:

pick_per_allele (lines 535-555):
  std::map<std::string, std::vector<AnnotationWithMeta>> by_allele;
  for (auto& ann : annotations) { std::string key = ann.annotation.ref_allele + ">" + ann.annotation.alt_allele; by_allele[key].push_back(std::move(ann)); }
  std::vector<AnnotationWithMeta> result;
  for (auto& pair : by_allele) { auto picked = pick_one(std::move(pair.second)); if (!picked.empty()) result.push_back(std:

**[filter-transcript/portability] include/filter_vep.hpp:462-540** — std::isspace called on raw char (possible negative value) — undefined behavior on high-bit bytes
- What: parse_filter_expression and load_gene_list call std::isspace(cond.field[0]) etc. with a plain char argument. If the input contains bytes >= 0x80 (e.g. UTF-8 field names/values), char may be negative and passing a negative value other than EOF to std::isspace is undefined behavior (the standard requires the argument be representable as unsigned char or EOF). This can crash or misbehave depending on the C library. The same file already shows the correct pattern elsewhere (transcript_filter.hpp line 26 casts to unsigned char), so this is an internal inconsistency.
- Evidence: while (!cond.field.empty() && std::isspace(cond.field[0])) { cond.field.erase(0, 1); }  // cond.field[0] is char, may be negative; cf. correct cast at transcript_filter.hpp:26
- Fix: Wrap each operator[] argument in a static_cast<unsigned char> at all 10 sites, matching the existing pattern at filter_vep.hpp:52 and transcript_filter.hpp:26. For example: `while (!cond.field.empty() && std::isspace(static_cast<unsigned char>(cond.field[cond.field.size() - 1])))` and `std::isspace(static_cast<unsigned char>(cond.field[0]))`. Apply identically to lines 462, 465, 478, 481, 491, 494, 504, 507, 537, 540. Optionally factor the repeated front/back-trim loops into a small helper (e.g. trim(std::string&)) that performs the cast once, eliminating the duplication and the chance of future regressions.
- Verifier: Verified against the actual current code in include/filter_vep.hpp. There are exactly 10 std::isspace calls at lines 462, 465, 478, 481, 491, 494, 504, 507, 537, 540, each passing a value obtained via std::string::operator[]. Examples:
- Line 462/465: `while (!cond.field.empty() && std::isspace(cond.field[cond.field.size() - 1]))` and `std::isspace(cond.field[0])`
- Line 491/494: `std::isspace(cond.value[...])`
- Line 504/507: `std::isspace(item[...])`
- Line 537/540: `std::isspace(line[...])`



**[gene-exon-plugin/correctness] include/gene_constraint.hpp:400-403** — format_constraint_score() suppresses legitimate negative Z-scores (mis_z/syn_z), printing '.'
- What: format_constraint_score uses value<0 as a 'missing' sentinel and returns '.'. That sentinel is valid for ratios/probabilities (pLI, oe_*), but mis_z and syn_z are Z-scores that are routinely negative for tolerant genes (e.g. mis_z = -1.5). The parser stores them faithfully (parse_double can return negative values), but any attempt to format mis_z/syn_z via this helper would turn every valid negative Z-score into '.', i.e. data loss / wrong output. The -1.0 'unset' default for mis_z/syn_z is also indistinguishable from a real Z-score of -1.0. This becomes an active correctness bug the moment H2 is fixed and these fields are emitted through this formatter.
- Evidence: Struct defaults: 'double mis_z = -1.0;' (line 42), 'double syn_z = -1.0;' (line 48). Formatter: 'inline std::string format_constraint_score(double value, int precision = 4) { if (value < 0 || std::isnan(value)) { return "."; } ...' (lines 400-403). parse_double returns the parsed value including negatives (line 176 'return std::stod(val);').
- Fix: Stop overloading negative values as a "missing" sentinel. Two complementary fixes: (1) For the Z-score fields (mis_z, syn_z), do not treat negatives as missing — they are legitimately negative. Initialize them to NaN instead of -1.0 (e.g. `double mis_z = std::numeric_limits<double>::quiet_NaN();`) and rely on the existing parse_double NA/./empty handling to mark them missing. (2) Make format_constraint_score's "missing" test use NaN (or an explicit has_value flag) rather than `value < 0`: `if (std::isnan(value)) return ".";`. Keep the `< 0` test only for ratio/probability fields if a separate formatter is desired, or pass a flag indicating whether negatives are valid. Since the formatter is currently dead code and the constraint DB is never queried for output, this can be fixed together with the downstream emission work (the reviewer's "H2") rather than urgently.
- Verifier: The cited code matches the claim exactly. include/gene_constraint.hpp line 401: `if (value < 0 || std::isnan(value)) { return "."; }`. Struct defaults at lines 42 and 48: `double mis_z = -1.0;` and `double syn_z = -1.0;`. parse_double at line 176 returns `std::stod(val)`, faithfully preserving negative values (it only maps empty/NA/./NaN to -1.0 at line 174). The reviewer's domain reasoning is correct: mis_z and syn_z are gnomAD missense/synonymous Z-scores that are routinely negative for constr

**[gene-exon-plugin/correctness] include/exon_intron_numbers.hpp:104, 126** — ExonIntronInfo.position_in_feature is computed from genomic start regardless of strand (wrong on minus strand)
- What: In get_exon_intron_number, position_in_feature is always 'position - start + 1' (and 'position - intron_start + 1' for introns), measuring offset from the low genomic coordinate. The exon/intron *number* is correctly strand-flipped (lines 99, 122), but the within-feature position is not: on the minus strand the transcript-relative offset should be measured from the high coordinate (end - position + 1). So for a minus-strand transcript the reported position_in_feature is the reverse-complement distance, inconsistent with the strand-aware numbering it sits next to. get_cds_exon_number has the same issue (lines 222, 242). Impact is currently limited because this module is only exercised by tests (production exon/intron numbering lives in vep_annotator.cpp), and the tests only cover '+' strand for position_in_feature; but if this helper is ever used for output it yields wrong coordinates on minus-strand features.
- Evidence: Lines 98-106: 'if (strand == '-') { info.number = num_exons - i; } else { info.number = i + 1; } info.position_in_feature = position - start + 1;' — the position_in_feature line is outside the strand branch and never uses 'end'. Same pattern at lines 216-222 of get_cds_exon_number.
- Fix: In get_exon_intron_number and get_cds_exon_number, make position_in_feature strand-aware to match the strand-flipped number sitting next to it. For exons on the minus strand set info.position_in_feature = end - position + 1; for introns set info.position_in_feature = intron_end - position + 1 (using the appropriate intron_end / coding-exon end). Keep the '+' branch as position - start + 1. feature_length is already correct and needs no change. Also add minus-strand assertions on position_in_feature in tests/test_exon_intron.cpp so the behavior is locked in. Optionally, drop the unused #include "exon_intron_numbers.hpp" in src/main.cpp, or wire this helper into production if it is intended to be the source of truth (currently vep_annotator.cpp has its own implementation).
- Verifier: The cited code matches the claim exactly. In include/exon_intron_numbers.hpp:

get_exon_intron_number, exon branch (lines 98-105):
  if (strand == '-') { info.number = num_exons - i; } else { info.number = i + 1; }
  info.position_in_feature = position - start + 1;
  info.feature_length = end - start + 1;

intron branch (lines 120-127):
  if (strand == '-') { info.number = num_exons - i - 1; } else { info.number = i + 1; }
  info.position_in_feature = position - intron_start + 1;

get_cds_exon_n

**[hgvs-parser/correctness] include/hgvs_parser.hpp:368-426** — end_intron_offset not stored for c. insertions and duplications, losing the end intronic offset
- What: The coding insertion and duplication regexes both capture the end-position intronic offset (group 4, `([+-]\d+)?`), but neither branch assigns it to result.end_intron_offset — only the deletion branch (line 354-356) does. So c.100+1_100+5dup parses with end_pos=100 but end_intron_offset=0, dropping the +5, and c.100+1_100+2insA loses its end offset. This corrupts the range for intronic insertions/duplications.
- Evidence: // ins branch:
        std::string end_str = match[3].str();
        if (!end_str.empty() && end_str[0] == '*') {
            result.end_pos = std::stoi(end_str.substr(1));
        } else {
            result.end_pos = std::stoi(end_str);
        }

        result.ref_allele = "";
        result.alt_allele = match[5].str();
- Fix: Apply the suggested fix for consistency: in the ins branch after parsing match[3] add `if (match[4].matched) result.end_intron_offset = std::stoi(match[4].str());`, and likewise inside the dup branch's `if (match[3].matched)` block. This mirrors the deletion branch and keeps the parser internally consistent. Note this is currently a latent fix (end_intron_offset is not yet consumed); the higher-value follow-up is to actually use end_pos/end_intron_offset when mapping ranged coding HGVS (del/dup/ins) to genomic coordinates in main.cpp:3045-3092, which today only maps the start position.
- Verifier: The cited code matches the claim. In include/hgvs_parser.hpp:

ins_regex (line 368): "(-?\\*?\\d+)([+-]\\d+)?_(-?\\*?\\d+)([+-]\\d+)?ins([ACGTacgt]+)" — group 4 = ([+-]\d+)? is the end-position intronic offset. The ins branch (lines 383-391) reads match[3] (end_pos) then match[5] (alt_allele), and never reads match[4]:
    std::string end_str = match[3].str();
    ... result.end_pos = std::stoi(end_str); ...
    result.ref_allele = "";
    result.alt_allele = match[5].str();

dup_regex (line 397

**[hgvs-parser/parity-gap] include/hgvs_parser.hpp:441-543** — Protein dup notation (p.Gly4_Gln6dup) not recognized
- What: parse_protein_hgvs handles del, ins, missense, frameshift and synonymous, but has no branch for protein duplication (p.Gly4dup, p.Gly4_Gln6dup), a standard HGVS protein form. These return valid=false. While protein-level parsing is less load-bearing than c./g., this is a parity gap with Perl VEP output that round-trips its own HGVSp.
- Evidence: // Parse insertion: Gly12_Ala13insVal
    static const std::regex ins_regex("([A-Z][a-z]{2}|[A-Z*])(\\d+)_([A-Z][a-z]{2}|[A-Z*])(\\d+)ins([A-Z][a-z]{2}|[A-Z*])+");
    ... (no dup branch follows)
- Fix: Add a protein dup branch in parse_protein_hgvs, mirroring the existing del_regex branch (lines 487-504). Suggested: static const std::regex dup_regex("([A-Z][a-z]{2}|[A-Z*])(\\d+)(?:_([A-Z][a-z]{2}|[A-Z*])(\\d+))?dup"); on match, set result.variant_type = HGVSVariantType::DUPLICATION (already in the enum), populate ref_aa (one->three conversion when single-char, as the del branch does), protein_pos/start_pos from group 2, and end_pos from group 4 if matched else equal to start_pos, then result.valid = true. Place it adjacent to the del branch. Add unit tests asserting valid==true and correct start/end positions for both "Gly4dup" and "Gly4_Gln6dup".
- Verifier: Confirmed against actual current code in /Users/davidcho/Projects/VEP_ensemble_C/include/hgvs_parser.hpp. The function parse_protein_hgvs (lines 441-544) has branches for: missense 3-letter (line 448), missense 1-letter (line 461), frameshift (line 472), deletion ("del_regex" line 487), insertion ("ins_regex" line 507), and synonymous ("syn_regex" line 528). There is NO branch matching "dup". Any input not matching these falls through to line 542: result.error_message = "Unrecognized protein HGV

**[hgvs-parser/parity-gap] include/hgvs_parser.hpp:282-430** — No support for inversion (inv) notation in c. or g. parsers despite INVERSION enum
- What: HGVSVariantType::INVERSION is defined but parse_coding_hgvs and parse_genomic_hgvs have no branch to parse `inv` notation (e.g. g.100_110inv, c.100_110inv). Such inputs return valid=false 'Unrecognized ... HGVS pattern'. This is a parity gap for a supported variant class.
- Evidence: result.error_message = "Unrecognized coding HGVS pattern: " + change;
    return result;
- Fix: Add an `inv` regex branch to both parse_genomic_hgvs and parse_coding_hgvs that sets variant_type = HGVSVariantType::INVERSION and captures start/end coordinates (e.g. genomic: `(\d+)_(\d+)inv`; coding with optional intron offsets mirroring the existing del/dup patterns). Note this only captures coordinates — to emit a correct alt allele the downstream code must reverse-complement the reference sequence spanning start..end, so the conversion-to-genomic stage (and CDS/transcript sequence lookup) must also be extended; otherwise the parsed result will have an empty/incorrect alt. If full inversion support is out of scope, document this as a known parity gap rather than leaving the unused INVERSION enum value implying support. Add a unit test asserting parse_genomic_hgvs and parse_coding_hgvs return valid=true with variant_type==INVERSION for inv inputs once implemented.
- Verifier: Verified against actual current code in include/hgvs_parser.hpp.

(1) HGVSVariantType::INVERSION is declared at line 47 ("INVERSION, // inv") but a codebase-wide grep shows it is NEVER assigned anywhere — its only occurrence is the enum declaration itself.

(2) parse_genomic_hgvs (lines 192-280) has branches only for SUBSTITUTION (sub_regex), DELETION (del_regex), INSERTION (ins_regex), DUPLICATION (dup_regex), and DELINS (delins_regex). There is no `inv` branch. Unmatched input reaches line 278

**[main-cli/quality] src/main.cpp:973-974** — --output-format and --terms accept invalid values silently, falling back to TSV/SO instead of erroring
- What: parse_output_format() (include/output_writer.hpp:38-47) returns OutputFormat::TSV for any unrecognized string, so a typo like `--output-format jsonn` or `--output-format gtf` silently produces TSV output with exit code 0 rather than reporting the error. Likewise --terms accepts any arbitrary style string. For a clinical tool this can cause the wrong output format to be generated without any warning. Perl VEP rejects unknown --format/--output_format values.
- Evidence: else if (arg == "--output-format" && i + 1 < argc) {
    output_format = vep::parse_output_format(argv[++i]);
}
- Fix: In src/main.cpp at the --output-format handler (line 973), validate the raw argument before calling parse_output_format. Either (a) compare argv[++i] (lowercased) against the known set {tsv, json, vcf} and on mismatch emit `std::cerr << "Error: unknown output format '" << val << "'" << std::endl;` and `return 1;`, or (b) change parse_output_format to signal failure (e.g. return std::optional<OutputFormat> or take a bool& ok out-param) and have the caller error out. Apply the same validation to --terms at line 1378-1379: accept only {SO, display} (the implemented set; drop or implement "NCBI" to match the comment at line 859), erroring on anything else. This matches Perl VEP, which rejects unknown --format/--output_format values.
- Verifier: Confirmed against the actual code. At src/main.cpp:973-974 the handler is `else if (arg == "--output-format" && i + 1 < argc) { output_format = vep::parse_output_format(argv[++i]); }`. The parser at include/output_writer.hpp:38-47 lowercases the string and returns OutputFormat::JSON only for "json", OutputFormat::VCF for "vcf", and `return OutputFormat::TSV;` for EVERYTHING else, with no error signal. So `--output-format jsonn` or `--output-format gtf` silently yields TSV.

Crucially, the unknow

**[output-writers/correctness] include/output_writer.hpp:1013-1021** — JSON key lowercasing can collapse distinct custom keys into duplicate JSON keys
- What: Each remaining custom-annotation key is force-lowercased before emission. If two distinct custom_annotations keys differ only in case (e.g. an upstream 'AF' field and a 'af' field, or 'gnomAD' vs 'gnomad' style keys from different sources), both lowercase to the same JSON key, producing two members with the same name in one object. JSON parsers keep only one, silently dropping the other; combined with the unordered_map iteration order (C2) which value survives is non-deterministic. This is silent annotation data loss.
- Evidence: std::string key = pair.first;
                for (auto& c : key) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
...
                    json << "        \"" << escape_json(key) << "\": \"" << escape_json(pair.second) << "\",\n";
- Fix: In the loop at include/output_writer.hpp:1010-1023, track already-emitted lowercased keys in a local std::set<std::string>. After computing the lowercased `key`, skip emission (or merge) if it is already in the set, and warn so data loss is not silent. For example:

  std::set<std::string> emitted_lc;
  for (const auto& pair : ann.custom_annotations) {
      if (skip_custom.count(pair.first) || pair.first.substr(0, 11) == "_colocated:") continue;
      if (pair.second.empty()) continue;
      std::string key = pair.first;
      for (auto& c : key) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
      if (!emitted_lc.insert(key).second) continue; // duplicate after lowercasing — skip (optionally log a warning)
      ... emit ...
  }

This guarantees a single transcript object never contains two members with the same name. Note that std::unordered_map iteration order is unspecified, so to make WHICH entry wins deterministic, consider iterating over an ordered copy (e.g. collect into a std::map) before emission. Optionally preserve original casing when a collision is detected to avoid losing the distinction entirely.
- Verifier: The cited code in include/output_writer.hpp (lines 1010-1023) matches the claim exactly:

  for (const auto& pair : ann.custom_annotations) {
      if (skip_custom.count(pair.first) || pair.first.substr(0, 11) == "_colocated:") continue;
      if (!pair.second.empty()) {
          // Perl VEP lowercases all JSON keys
          std::string key = pair.first;
          for (auto& c : key) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
          if (is_valid_json_number(pair.sec

**[output-writers/correctness] include/output_writer.hpp:1423-1441** — escape_vcf does not escape newline/CR, allowing a CSQ value to split the VCF data line
- What: escape_vcf percent-encodes %,|,;,=,&,space,tab but NOT '\n' or '\r'. A control character in any annotation value carried into the CSQ field (custom BED/text annotation 'name', or any free-text source value) passes through verbatim and breaks the single VCF data line into two physical lines, corrupting the record and every downstream tool. VCF 4.x requires INFO values be free of embedded line breaks. Verified by compiling escape_vcf in isolation: input containing '\n' produced output still containing a literal newline.
- Evidence: switch (c) {
                case '%': result += "%25"; break;
                case '|': result += "%7C"; break;
...
                case '\t': result += "%09"; break;
                default: result += c; break;  // '\n' and '\r' fall through unescaped
            }
- Fix: In escape_vcf (include/output_writer.hpp ~1428), add cases so embedded control characters cannot break VCF line structure:
  case '\n': result += "%0A"; break;
  case '\r': result += "%0D"; break;
Ideally also encode any other byte < 0x20 (e.g. a default branch: `if ((unsigned char)c < 0x20) { append %XX } else result += c;`). Separately, harden domain_source.cpp's load loop (and split_line callers generally) to strip a trailing '\r' after std::getline, matching the VCF reader at vep_annotator.cpp:747-749, so CRLF files do not leak '\r' into field values in the first place.
- Verifier: The cited code at include/output_writer.hpp:1423-1441 matches the claim exactly. escape_vcf() (default no_escape_=false) percent-encodes %, |, ',', ;, =, &, space, and tab, but '\n' and '\r' hit `default: result += c;` and are emitted verbatim:

  case '\t': result += "%09"; break;
  default: result += c; break;

So the literal defect (no escaping of CR/LF) is real, and the reviewer's isolated-compile observation is accurate. escape_vcf is precisely the centralized chokepoint that should guarant

**[output-writers/correctness] include/output_writer.hpp:473-475, 507-523, 539-540** — JSONWriter::close() does not flush the buffered final variant; data loss if write_footer() is skipped
- What: JSONWriter buffers the last variant's annotations in buffered_annotations_ and only emits them in flush_current_variant(), which is called solely from write_footer(). Unlike VCFWriter::close() (which calls flush_variant()), JSONWriter::close() and the destructor do NOT flush. If the program closes the writer without calling write_footer() (e.g. during stack unwinding from an exception, or any future caller that forgets the footer), the entire last variant's JSON object is silently dropped from the output. The normal main.cpp path calls write_footer() first so production output is currently complete, but the class invariant is fragile and inconsistent with VCFWriter.
- Evidence: void write_footer() override {
        flush_current_variant();
        // NDJSON format: no closing bracket
    }
    void close() override {
        if (gz_file_) { ... }   // no flush_current_variant() here
        if (output_.is_open()) { output_.close(); }
    }
- Fix: Add flush_current_variant(); as the first statement of JSONWriter::close() (include/output_writer.hpp line 512), mirroring VCFWriter::close() at line 1292. It is idempotent (empty-guard at line 540 and buffer clear at line 780), so it remains correct even though write_footer() also calls it on the normal path. This restores the class invariant that destruction always emits all buffered output, consistent with VCFWriter, with zero risk of duplicate variants.
- Verifier: Verified against actual code in /Users/davidcho/Projects/VEP_ensemble_C/include/output_writer.hpp. The claim is accurate.

JSONWriter::close() (lines 512-523) does NOT flush the buffered final variant:
    void close() override {
        if (gz_file_) { int ret = gzclose(gz_file_); gz_file_ = nullptr; if (ret != Z_OK) {...} }
        if (output_.is_open()) { output_.close(); }
    }
The only flush of the last variant happens in write_footer() (lines 507-510):
    void write_footer() override { f

**[output-writers/parity-gap] include/output_writer.hpp:629-630** — JSON colocated_variants allele_string omits Ensembl '-' substitution used at variant top level
- What: The top-level allele_string substitutes '-' for empty ref/alt (lines 590-593) to match Ensembl format, but the colocated_variants object emits raw first.ref_allele/first.alt_allele with no substitution. For an insertion (empty ref) this yields allele_string ":/A" style values like "/A" instead of "-/A", diverging from the top-level allele_string and from Perl VEP's Ensembl-format colocated allele_string.
- Evidence: json << "      \"allele_string\": \"" << escape_json(first.ref_allele) << "/"
                     << escape_json(first.alt_allele) << "\",\n";
- Fix: Mirror the top-level substitution when building the colocated_variants allele_string. Replace lines 629-630 with:
  std::string cv_ref = first.ref_allele.empty() ? "-" : first.ref_allele;
  std::string cv_alt = first.alt_allele.empty() ? "-" : first.alt_allele;
  json << "      \"allele_string\": \"" << escape_json(cv_ref) << "/"
       << escape_json(cv_alt) << "\",\n";
This makes colocated allele_string consistent with the top-level allele_string (lines 590-593) and the id field (lines 579-583), matching Ensembl format. Add a JSON-output regression test for an insertion (empty ref) with a populated existing_variation to lock in "-/A".
- Verifier: Confirmed against the actual current code in include/output_writer.hpp.

Top-level allele_string (lines 590-593) applies the empty->"-" substitution to match Ensembl format:
  std::string as_ref = first.ref_allele.empty() ? "-" : first.ref_allele;
  std::string as_alt = first.alt_allele.empty() ? "-" : first.alt_allele;
  json << "    \"allele_string\": \"" << escape_json(as_ref) << "/" << escape_json(as_alt) << "\",\n";

The colocated_variants allele_string (lines 629-630) emits raw values with

**[perf/performance] include/output_writer.hpp:432-441** — Output to stdout uses std::cout without sync_with_stdio(false); every record string goes through synchronized C/C++ stdio
- What: All three writers' write_string() send output to stdout via 'std::cout << s' when use_stdout_ is true (the default). main.cpp never calls std::ios::sync_with_stdio(false) (confirmed: no occurrence anywhere in src/ or include/). With the default sync flag, libstdc++/libc++ route every std::cout operation through the C stdio layer with per-operation synchronization, which is materially slower for high-volume output (one string per annotation line, potentially millions of lines for WGS). This is the common case since the default output target is stdout.
- Evidence: void write_string(const std::string& s) {
        if (use_stdout_) {
            std::cout << s;
        } else if (compress_ && gz_file_) {
            ...
        } else {
            output_ << s;
        }
    }
- Fix: Add std::ios::sync_with_stdio(false); once at the top of main() in src/main.cpp. It is safe because no printf/C-stdio output is interleaved with std::cout anywhere in the codebase. This benefits the explicit stdout/pipe workflow (-o - or -o STDOUT) on large batch runs. Note that the default no-flag batch path writes to a derived file via std::ofstream and is already unaffected, so do not expect a speedup there. Optionally, for the file path, the existing std::ofstream could be given a larger buffer, but that is separate from this stdio-sync issue.
- Verifier: The cited code is accurate. include/output_writer.hpp:432-441 contains exactly the quoted TSVWriter::write_string, and the same std::cout << s pattern appears in JSONWriter (line 1040) and VCFWriter (line 1322). use_stdout_ is true when output_path is empty, "-", or "STDOUT" (lines 205, 455, 1108). I confirmed via grep that sync_with_stdio is never called anywhere in src/ or include/, and that NO C stdio output (printf/fprintf/fputs/fwrite/putchar) is mixed on stdout — the program uses std::cout

**[perf/performance] include/output_writer.hpp:435-438** — gzwrite called once per record with no user-space buffering (unbuffered compressed output)
- What: When writing compressed output (.gz), write_string() issues one gzwrite() call per record string (per annotation line). gzwrite incurs the deflate state machine entry and a function-call boundary per call; for millions of small per-line writes this is slower than accumulating into a buffer and flushing in large chunks. The same pattern repeats in JSONWriter (line 1042) and VCFWriter (line 1324). gzbuffer() is never called to enlarge the internal buffer either.
- Evidence: } else if (compress_ && gz_file_) {
            if (gzwrite(gz_file_, s.c_str(), static_cast<unsigned int>(s.size())) != static_cast<int>(s.size())) {
                throw std::runtime_error("Failed to write to compressed output");
            }
        }
- Fix: Two complementary, low-risk changes. (1) After each gzopen("wb") call (lines ~211, ~461, ~1114), call gzbuffer(gz_file_, 1 << 17) to enlarge zlib's internal buffer from the 8 KB default to 128 KB, reducing per-call processing. (2) Optionally add a std::string accumulation buffer in each write_string: append s to the buffer, and only issue a single gzwrite when the buffer exceeds ~256 KB, with a final flush in close() (the close() at line 411 already exists and would need a flush-before-gzclose). Since the three write_string implementations are identical triplicated code, consider factoring the gz buffering into a shared helper to avoid divergence. Keep the existing return-value check on gzwrite for error handling.
- Verifier: The cited code matches exactly. In include/output_writer.hpp the TSVWriter::write_string (lines 432-442) does:

  } else if (compress_ && gz_file_) {
      if (gzwrite(gz_file_, s.c_str(), static_cast<unsigned int>(s.size())) != static_cast<int>(s.size())) {
          throw std::runtime_error("Failed to write to compressed output");
      }
  }

The identical pattern is duplicated in JSONWriter::write_string (lines 1038-1048) and VCFWriter::write_string (lines 1320-1330). The per-record call pat

**[sources-a/correctness] src/sources/spliceai_source.cpp:226-229** — SpliceAI INFO key matched by substring, can match a different key ending in 'SpliceAI='
- What: The INFO field for SpliceAI is located with info.find("SpliceAI="), an unanchored substring search. A different INFO key whose name ends in 'SpliceAI' (e.g. a custom 'AnnotSpliceAI=...' or a key appearing before the real one) could be matched, or if the real SpliceAI value contains the substring within an earlier field it could anchor at the wrong offset. INFO keys are ';'-delimited and key boundaries should be respected (start of INFO or after a ';').
- Evidence: size_t start = info.find("SpliceAI=");
if (start == std::string::npos) return;
std::string spliceai_value = info.substr(start + 9);
- Fix: Match the SpliceAI key only at an INFO field boundary. Minimal fix: after each find, require start==0 || info[start-1]==';'; if not, search forward from start+1 in a loop until a boundary match or npos. Robust fix: split info on ';' into KEY=VALUE tokens, then select the token whose key (substring before the first '=') equals exactly \"SpliceAI\", and parse its value. This also naturally handles the existing ';' termination. Add a unit test feeding parse_spliceai_info an INFO string like \"AnnotSpliceAI=foo;SpliceAI=A|GENE|0.01|...\" to assert the correct key is selected.
- Verifier: Verified against the actual current code in src/sources/spliceai_source.cpp. The cited lines are present (now at 226-229) inside parse_spliceai_info():

    size_t start = info.find("SpliceAI=");
    if (start == std::string::npos) return;
    std::string spliceai_value = info.substr(start + 9);
    size_t end = spliceai_value.find(';');
    if (end != std::string::npos) spliceai_value = spliceai_value.substr(0, end);

I traced the data flow: `info` is the raw VCF INFO column (8th field). In fil

**[sources-a/correctness] src/sources/dbnsfp_source.cpp:226-247** — dbNSFP query() path returns first value of multi-value fields regardless of allele/transcript
- What: The const query() entry point has no transcript context and explicitly selects the first valid value of any ';'-separated field via select_value_for_transcript(values, -1, alt_index). For per-transcript score columns (SIFT_score, Polyphen2_*), the first value corresponds to whatever transcript dbNSFP listed first, not the transcript/allele relevant to the caller. Any consumer using query() (as opposed to annotate()) thus receives an arbitrary transcript's score, which can be the wrong prediction for the variant. Combined with C1, the alt_index passed here is also ignored.
- Evidence: // No transcript context available, use first valid value
value = select_value_for_transcript(values, -1, alt_index);
- Fix: Document VariantAnnotationSource::query() (and the dbNSFP override) as transcript-agnostic in its doc comment, and make the behavior non-misleading: for known per-transcript columns either (a) return all per-transcript values joined with ';' to mirror the raw dbNSFP field, or (b) return "." / require the caller to disambiguate, rather than silently returning the first transcript's value. Separately, since the alt_index parameter is entirely unused in select_value_for_transcript, either implement per-alt selection or drop the parameter to remove the dead argument at lines 161-162 and 240-241. Because no caller currently exercises this query() path, this can be treated as low-priority cleanup/documentation rather than an urgent fix.
- Verifier: The cited code matches the claim. In src/sources/dbnsfp_source.cpp the VariantAnnotationSource::query() override (lines 188-253) has no transcript context by interface (declared at include/annotation_source.hpp:202-207 with no Transcript* param). At lines 239-241 it does: `// No transcript context available, use first valid value` / `value = select_value_for_transcript(values, -1, alt_index);`. The helper select_value_for_transcript (lines 314-337) with transcript_index == -1 skips the matched b

**[sources-b/correctness] src/sources/regulatory_source.cpp:105-135** — regulatory:count under-reports overlapping features (counts only features with a non-empty ID)
- What: regulatory:count is set to ids.size(), but ids is only appended when feature->id is non-empty (line 105-107). GFF3 features without an ID attribute get feature->id == "" (file_parsers.cpp:634), so they are counted in feature_type (via the types set) but excluded from the count. This produces inconsistent output: a variant can have regulatory:feature_type populated while regulatory:count is 0 or lower than the true number of overlapping features. Ensembl Regulatory Build TFBS/motif features in particular may lack an ID attribute.
- Evidence: types.insert(feature->type);
if (!feature->id.empty()) {
    ids.push_back(feature->id);
}
...
annotations["regulatory:count"] = std::to_string(ids.size());
- Fix: Introduce a separate counter (e.g., `int matched_count = 0;`) incremented once per feature that survives the cell-type filter, right alongside `types.insert(feature->type);` at line 104. Then report `annotations["regulatory:count"] = std::to_string(matched_count);` at line 135 instead of `ids.size()`. Keep `ids` solely for building the truncated `regulatory:feature_id` list. Optionally, to avoid the inverse double-count problem, note that `ids` is a vector so duplicate IDs across features are listed twice; using matched_count for the count avoids relying on ids at all, which is the cleaner fix.
- Verifier: The cited code matches the claim exactly. In src/sources/regulatory_source.cpp the per-feature loop (line 76) unconditionally does `types.insert(feature->type);` (line 104) for every feature that passes the cell-type filter, but only does `ids.push_back(feature->id);` when `if (!feature->id.empty())` (lines 105-107). Then line 135 sets `annotations["regulatory:count"] = std::to_string(ids.size());`. So the reported count is the number of features that carry a non-empty ID, NOT the number of over

**[sources-b/correctness] src/sources/conservation_source.cpp:78-88** — Conservation indel mean includes the unchanged VCF anchor base
- What: For indels the source computes the conservation mean over start=pos .. end=pos+ref.length()-1. The variant ref/alt reaching annotate are raw VCF-style alleles (passed straight through from annotate_all without anchor stripping, vep_annotator.cpp:1844), so for a typical VCF deletion (e.g. REF=ATG, ALT=A) pos refers to the unchanged anchor base 'A' and the window includes it, slightly biasing the reported phyloP/phastCons/GERP mean toward the wrong (unaffected) base. For pure insertions (REF=A, ALT=ATG) the window collapses to the single anchor base only. This affects the numeric conservation values reported for indels.
- Evidence: } else {
    // For indels, get mean score across affected region
    int start = pos;
    int end = pos + static_cast<int>(ref.length()) - 1;
    auto mean = reader_->get_mean(chrom, start, end);
    if (mean.has_value()) { ... }
}
- Fix: In BigWigScoreSource::annotate (conservation_source.cpp lines 78-88), strip the shared leading anchor base before computing the scoring window, mirroring the anchor convention already used in vep_annotator.cpp (e.g. line 2190). Concretely: compute a shared-prefix length when `!ref.empty() && !alt.empty() && ref[0] == alt[0] && ref.length() != alt.length()`, then for a deletion score only the deleted bases (start = pos + shared_prefix_len, end = pos + ref.length() - 1). For a pure insertion (deleted region has zero genomic width after stripping), fall back to scoring the two genomic bases flanking the insertion point rather than the single unchanged anchor. Keeping the implementation consistent with the existing --minimal trimming logic (main.cpp lines 2648-2664) would avoid duplicating two divergent conventions.
- Verifier: The cited code in src/sources/conservation_source.cpp (lines 78-88) matches the claim exactly. For non-SNP variants: `int start = pos; int end = pos + static_cast<int>(ref.length()) - 1; auto mean = reader_->get_mean(chrom, start, end);`. I verified BigWigReader::get_mean (src/file_parsers.cpp:418) calls `bwStats(..., start - 1, end, 1, mean)`, a 0-based half-open [start-1, end) interval == 1-based inclusive [start, end]. So the window is inclusive of `pos`.

Premise verified: VEPAnnotator::anno

**[structural-variants/correctness] include/structural_variant.hpp:453-456** — BND mate orientation derived solely from bracket position, ignoring bracket type, producing wrong strand for half of breakends
- What: bnd_mate_forward is set to `(bracket_pos > 0)`, i.e. true whenever the breakend sequence appears before the bracket (forms t[p[ or t]p]) and false otherwise. But VCF breakend strand/join orientation is encoded by the bracket *character* ('[' vs ']'), not by whether the base is before or after the brackets. The four canonical forms (t[p[, t]p], ]p]t, [p[t) encode distinct join orientations; `(bracket_pos>0)` collapses them to just 'base-before' vs 'base-after' and loses the bracket-type information. Two of the four orientations are therefore labeled incorrectly, which matters for any fusion/strand-aware downstream interpretation of the mate.
- Evidence: // VCF BND orientation: sequence before brackets = forward on local end
// t[p[ or t]p] = forward (bases before first bracket)
// ]p]t or [p[t = reverse complement (bases after last bracket)
sv.bnd_mate_forward = (bracket_pos > 0);
- Fix: Decode the mate orientation from the bracket CHARACTER, not just position. Record the bracket char actually found (e.g. char br = (alt.find('[') != npos) ? '[' : ']') and set sv.bnd_mate_forward = (br == '['). If both the local-end and mate orientations are needed for fusion/strand work, add a second field (e.g. bnd_local_forward = (bracket_pos > 0)) so the four canonical forms (t[p[, t]p], ]p]t, [p[t) are each represented distinctly per VCF 4.x. Also update/extend the tests in tests/test_structural_variant.cpp to cover all four forms, since the current two cases mask the bug.
- Verifier: The cited code matches the claim exactly. At include/structural_variant.hpp:456:

    // VCF BND orientation: sequence before brackets = forward on local end
    // t[p[ or t]p] = forward (bases before first bracket)
    // ]p]t or [p[t = reverse complement (bases after last bracket)
    sv.bnd_mate_forward = (bracket_pos > 0);

bracket_pos is found at lines 437-440 as the first '[' (or, if absent, the first ']'). So bracket_pos > 0 is true exactly when the base 't' precedes the bracket, i.e. it

**[structural-variants/correctness] include/structural_variant.hpp:360, 384-389, 108-115** — Symbolic SV span overlaps the padding base, causing spurious 5' UTR / upstream / exon overlaps off by one
- What: parse_sv_from_vcf sets sv.start = pos directly. For VCF symbolic alleles (<DEL>, <DUP>, <INV>, <INS>, <CNV>) POS is the base *immediately before* the event; the affected reference span is POS+1..END (per VCF spec). Because start is left at POS, overlaps()/contains() include one extra base to the left of the true event. This produces off-by-one false positives at feature boundaries: e.g. a deletion whose true span begins exactly at an exon/CDS/UTR start, but whose padding base sits in the preceding intron/UTR, will be reported as also touching that preceding region (extra splice_acceptor / five_prime_UTR / upstream calls), and contains() can fail to register full transcript ablation when the transcript starts one base after POS. Perl VEP works on the POS+1..END span for symbolic SVs.
- Evidence: sv.start = pos;
...
} else {
    // Estimate end from ref/alt lengths
    if (sv.sv_type == SVType::DEL) {
        sv.end = pos + static_cast<int>(ref.size()) - 1;
    } else {
        sv.end = pos;
    }
}
- Fix: In parse_sv_from_vcf, detect symbolic alleles (alt non-empty and alt[0]=='<') and, for span-based SV types where POS is the padding base (DEL, DUP, INV, CNV), set sv.start = pos + 1 so the event span is POS+1..END, matching the VCF spec and Perl VEP. Leave non-symbolic small indels untouched (they already use ref-length-derived spans). Update ParseSVFromVCF.DeletionWithEnd and the CNV tests to expect start==pos+1. Do NOT shift start for BND (POS is the actual breakend base) or for symbolic INS (zero-length on the reference between POS and POS+1, handled separately). Optionally verify get_transcripts_in_region callers still pass the corrected sv.start.
- Verifier: VERIFIED AGAINST CURRENT CODE. The cited code matches the claim:

include/structural_variant.hpp:360 `sv.start = pos;` sets start to the raw VCF POS unconditionally, with no special handling for symbolic alleles. Lines 384-389 estimate END only from ref length for non-symbolic DEL. overlaps()/contains() (lines 108-115) use start inclusively: `start <= e && end >= s` and `start <= s && end >= e`.

Both real call sites pass the raw POS with no padding adjustment: src/main.cpp:2673 `parse_sv_from_v

**[structural-variants/parity-gap] include/structural_variant.hpp:295-309** — CNV consequence logic omits copy_number==1 and unknown copy number (-1)
- What: The CNV branch only handles copy_number==0 (ablation/frameshift) and copy_number>2 (amplification). copy_number==1 (a heterozygous loss / single-copy state, a genuine deletion relative to diploid) and copy_number==-1 (unknown, the default when neither CN INFO nor <CNn> ALT is present) fall through with no CDS/splice consequence assigned in this block. A CNV overlapping CDS with CN=1 or unknown CN therefore gets no coding consequence here and relies on the generic exon/intron fallback (only emitting coding_sequence_variant if it happens to land in an exon and consequences are otherwise empty), losing the deletion semantics Perl VEP would assign.
- Evidence: } else if (sv.sv_type == SVType::CNV) {
    if (sv.copy_number == 0) {
        ...
    } else if (sv.copy_number > 2) {
        // Duplication/amplification
        if (affects_cds) {
            consequences.push_back(ConsequenceType::CODING_SEQUENCE_VARIANT);
        }
    }
}
- Fix: In the CNV branch (lines 295-308), add handling for copy_number==1 and copy_number==-1 (unknown). Treat copy_number==1 as a single-copy loss and copy_number==-1 as an indeterminate copy-number change: when affects_cds, push CODING_SEQUENCE_VARIANT (and FEATURE_TRUNCATION when the CNV partially overlaps the transcript), so the coding effect is assigned directly in the block rather than relying on the empty()-gated exon/intron fallback that is skipped whenever a UTR is also overlapped. Optionally also add splice-dinucleotide detection to the CNV branch to mirror the DEL/INS branches. A simple robust fix is to make the CNV branch default (else, covering CN==1 and CN==-1) emit CODING_SEQUENCE_VARIANT when affects_cds, ensuring parity regardless of concurrent UTR overlap.
- Verifier: Verified against current code in include/structural_variant.hpp. The CNV consequence branch (lines 295-308) is exactly:

  } else if (sv.sv_type == SVType::CNV) {
      if (sv.copy_number == 0) {
          if (sv.contains(transcript.start, transcript.end)) {
              consequences.push_back(ConsequenceType::TRANSCRIPT_ABLATION);
          } else if (affects_cds) {
              consequences.push_back(ConsequenceType::FRAMESHIFT_VARIANT);
          }
      } else if (sv.copy_number > 2) {
   

**[structural-variants/correctness] include/structural_variant.hpp:255, 278-289** — DUP/TDUP affecting CDS reuses INS frameshift logic, mislabeling tandem duplications
- What: DUP and TDUP are grouped with INS in the same branch, so a duplication affecting the CDS is classified as inframe_insertion or frameshift_variant based on sv.length()%3. A duplication is not an insertion of arbitrary novel sequence; its effect on reading frame depends on the duplicated CDS bases and the junction, and Perl VEP classifies whole-feature DUPs as transcript_amplification (handled above) and partial DUPs differently. Combined with C2, partial coding duplications get an inframe/frameshift label computed from the full genomic span, which can disagree with Perl VEP for tandem duplications.
- Evidence: } else if (sv.sv_type == SVType::INS || sv.sv_type == SVType::DUP || sv.sv_type == SVType::TDUP) {
    ...
    int inserted_bp = sv.length();
    if (inserted_bp % 3 != 0) {
        consequences.push_back(ConsequenceType::FRAMESHIFT_VARIANT);
    } else {
        consequences.push_back(ConsequenceType::INFRAME_INSERTION);
    }
- Fix: Split DUP/TDUP out of the INS branch at line 255. For a partial DUP/TDUP that affects the CDS, emit feature-level terms consistent with Perl VEP — e.g. CODING_SEQUENCE_VARIANT (and FEATURE_ELONGATION, which is currently only added for INS at line 170) — instead of computing inframe_insertion/frameshift_variant from `sv.length() % 3`. If a frame determination is desired, base it on the count of duplicated CDS bases overlapping the transcript's cds_regions rather than the full genomic span. This also aligns the DUP path with the existing DEL partial-overlap convention documented in the DELPartialOverlapExon test.
- Verifier: The cited code matches the claim exactly. At include/structural_variant.hpp line 255: `} else if (sv.sv_type == SVType::INS || sv.sv_type == SVType::DUP || sv.sv_type == SVType::TDUP) {`, and at lines 282-289 the CDS-affecting case computes `int inserted_bp = sv.length();` then `if (inserted_bp % 3 != 0) { FRAMESHIFT_VARIANT } else { INFRAME_INSERTION }`. So a partial coding DUP/TDUP is classified using the same arbitrary-novel-insertion logic as INS.

Two supporting facts check out:
1. Whole-fe


## Rejected by adversarial verification

- src/sources/dbnsfp_source.cpp:161-163, 314-337 — dbNSFP alt-index is computed but discarded; per-alt fields use wrong value: The cited code matches the textual claim but the HIGH "wrong alternate allele" framing is refuted by the actual data model.

Verified code facts: at src/sources/dbnsfp_source.cpp:161-162, `value = select_value_for_transcript(values, transcript_index, alt_index);` does pass alt_index. At line 314-317
- src/sources/dbnsfp_source.cpp:104-123, 204-223 — dbNSFP alt-column matching splits on ';' but real dbNSFP alt is a single base per row: The cited code matches the quote. In annotate() (lines 104-123) and query() (lines 204-223), the 'alt' column is split on ';' and compared base-by-base to find found_alt and an alt_index. The reviewer is factually correct that standard dbNSFP rows hold a single-base alt per row (the ';' separation a
- src/main.cpp:638-659 — gz_read_line can split a line at a 65535-byte buffer boundary, corrupting very long gzipped VCF lines: The cited code at src/main.cpp:638-659 matches the quote exactly:

  if (gzgets(gz, buf, CHUNK) == nullptr) { eof = result.empty(); break; }
  result += buf;
  if (!result.empty() && result.back() == '\n') break;
  size_t len = std::strlen(buf);
  if (static_cast<int>(len) < CHUNK - 1) break;

Howev
- include/structural_variant.hpp:351-477 — CIPOS/CIEND confidence intervals are never parsed or applied: The factual code observation is accurate but the parity-gap impact claim is not substantiated. Verified facts from the actual code:

1. `StructuralVariant` struct (include/structural_variant.hpp:83-116) has no CIPOS/CIEND fields — only chromosome, start, end, ref_allele, alt_allele, sv_type, sv_len,