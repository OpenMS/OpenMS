# ProSE: fixes for issue #9197 items 1 and 3

**Date:** 2026-04-21
**Scope:** `src/topp/ProSE.cpp` (per-file output sequencing, merged-write isolation,
deferred FDR:PSM application, dangling protein reference cleanup).
**Source issue:** [OpenMS#9197](https://github.com/OpenMS/OpenMS/issues/9197) items 1 and 3.
**Related follow-up:** [OpenMS#9198](https://github.com/OpenMS/OpenMS/issues/9198) (internal PEP/FDR fallback when Percolator fails).
**Note on line numbers:** all `ProSE.cpp:NNN` references in this spec point to the pre-patch state at branch `develop` HEAD; they will shift during implementation.

## Background

Two distinct ProSE bugs were observed in a multi-file MHC-I run:

1. **Data loss on merged-write failure (item 1a).** `ProSE.cpp` writes the merged
   idXML at line 397 BEFORE the per-file `out_idxml` writes at line 441. When the
   merged write threw, the exception aborted `main_` and all per-file outputs were
   discarded — even though the per-file `protein_ids`/`peptide_ids` were already
   complete in memory.

2. **The merged write itself crashed on `target+decoy` PSMs (item 1b).** Root
   cause: in the merged path, `applyPickedProteinFDR` (with default
   `add_decoy_proteins=false`, routing through `IDScoreGetterSetter::
   setScoresAndRemoveDecoys_`) drops decoy `ProteinHit` entries while
   reassigning q-value scores; `filterHitsByScore` immediately afterward drops
   the remaining target proteins above the q-value cutoff. Net: decoy proteins
   are gone from the merged protein list. But `IDFilter::removeDecoyHits` tests
   `target_decoy == "decoy"` exactly via `HasDecoyAnnotation`, so it **keeps**
   PSMs annotated as `target+decoy` (peptide sequences that map to both a
   target and a decoy protein, common in short non-tryptic immunopeptidomics).
   Those kept PSMs retain `PeptideEvidence` entries pointing to the
   now-deleted decoy proteins. `IdXMLFile::store` validates every evidence
   accession against the protein list and throws "Invalid protein reference
   'DECOY_...'", causing the merged write to fail at the same time as item
   1a's data loss. Bug only triggers when `FDR:protein > 0`; with
   `FDR:protein == 0` no decoy proteins are removed and there's nothing to
   dangle against.

3. **`FDR:PSM` strips decoys before Percolator (item 3).** The wrapper at
   `ProSE.cpp:191-195` zeros `FDR:protein` to 0 when `-percolator_executable`
   is set but does NOT zero `FDR:PSM`. Inside `ProSEAlgorithm::run()`,
   `removeDecoyHits` is called immediately after PSM-level FDR filtering,
   conditionally on `fdr_protein_ == 0.0` (`ProSEAlgorithm.cpp:1264, 1538`).
   The wrapper's zeroing of `FDR:protein` makes the algorithm see
   `fdr_protein_ == 0.0` → it triggers `removeDecoyHits` after PSM FDR. With
   `FDR:PSM > 0`, decoys are stripped from the algorithm's output before
   Percolator runs. `PercolatorAdapter` then aborts with "No decoys found" and
   ProSE silently falls back to raw HyperScore. The user only sees a generic
   "rescoring failed" message and has to dig into stderr to find the cause.

## Goal

Fix all three behaviors in a single PR. Output writes become independent and
idempotent: a failure in one mode (per-file or merged) does not discard the
others. PSM FDR is correctly deferred to post-Percolator and re-applied using
Percolator q-values. The merged-write protein-reference invariant is restored.

## Non-goals

- Item 2 from issue #9197 (sibling executable PATH resolution) — re-scoped as a
  separate cross-cutting OpenMS-CLI issue (`OpenNuXL` uses the same pattern).
- Item 4 (`.pin` missing standard columns) — independent fix; out of scope.
- Reworking `ProSEAlgorithm` internals — all changes here are confined to
  `ProSE.cpp` (the TOPP wrapper).
- Falling back to internal PEP/FDR when Percolator fails (or is skipped due to
  no decoys for that file). The score-scale heterogeneity that results — merged
  output mixing Percolator q-values and raw HyperScore PSMs across input files
  when some Percolator runs fail — is a direct consequence of this PR's deferral
  design, not a pre-existing concern. Tracked separately as #9198.

## Architecture

### Per-file output sequencing

```
for i in 0..N-1:
  result = mfres.per_file[i]                  (raw HyperScore PSMs)

  # 1. pin first — preserves the "raw input for external rescoring" semantic
  if (-out_pin set):
    try:    write -out_pin[i] from raw result
    catch:  log error, set input_failed

  # 2. Percolator rescoring (if requested) — mutates result.peptide_ids in place
  if (percolator_executable set AND has_decoys AND PSMs >= 100):
    try:    run PercolatorAdapter, load rescored idXML back into result
    catch:  log warn (mention deferred FDR:PSM not applied), keep raw result

  # 3. Re-apply deferred PSM FDR using Percolator q-values
  if (user_psm_fdr > 0.0 AND percolator_executable set AND rescored OK):
    IDFilter::filterHitsByScore(result.peptide_ids, user_psm_fdr)
    IDFilter::removeEmptyIdentifications(result.peptide_ids)

  # 4. Other outputs from the rescored, FDR-filtered result
  if (-out_idxml set):    try {...} catch {log}
  if (-out_qpx set):      try {...} catch {log}
  if (-out_parquet set):  try {...} catch {log}
```

Pin moves above Percolator. The current code writes pin AFTER Percolator runs;
the change makes pin contain raw HyperScore PSMs, which matches the parameter
description (`"for external rescoring"`).

### Merged output isolation

```
if (-out_merged set AND N > 1):
  try:
    merge per-file results via IDMergerAlgorithm
    BasicProteinInferenceAlgorithm.run(merged)
    if (user_protein_fdr > 0.0):
      applyPickedProteinFDR + filterHitsByScore     (deletes decoy ProteinHits)
    IDFilter::removeDecoyHits(merged_peptides)      (keeps target+decoy PSMs)
    IDFilter::removeEmptyIdentifications(merged_peptides)
    IDFilter::removeUnreferencedProteins(merged_protein_ids, merged_peptides)
    IDFilter::removeDanglingProteinReferences(merged_peptides, merged_protein_ids)  [NEW]
    storeIdentifications(out_merged, ...)
  catch:
    log error, do NOT propagate — per-file outputs are already on disk
```

The whole merged block is wrapped in try/catch matching the per-mode try/catch
already used in the per-file loop. The new `removeDanglingProteinReferences`
call closes the `target+decoy` evidence-pointing-to-deleted-decoy gap that
caused the original crash.

The merged path moves AFTER the per-file loop. Per-file writes complete first;
merged is the last output produced. Order of operations within the per-file
loop is unchanged from current except for pin moving above Percolator and the
new FDR-re-apply step.

### FDR deferral

Extend the existing `FDR:protein` deferral (`ProSE.cpp:191-195`) to also handle
`FDR:PSM`:

```cpp
const String percolator_executable = getStringOption_("percolator_executable");
const double user_protein_fdr = static_cast<double>(search_params.getValue("FDR:protein"));
const double user_psm_fdr     = static_cast<double>(search_params.getValue("FDR:PSM"));

if (!percolator_executable.empty())
{
  if (user_protein_fdr > 0.0)
  {
    OPENMS_LOG_INFO << "[ProSE] Percolator rescoring enabled: deferring protein FDR to post-rescoring." << endl;
    search_params.setValue("FDR:protein", 0.0);
  }
  if (user_psm_fdr > 0.0)
  {
    OPENMS_LOG_INFO << "[ProSE] Percolator rescoring enabled: deferring PSM FDR to post-rescoring." << endl;
    search_params.setValue("FDR:PSM", 0.0);
  }
}
```

Re-apply per-file after Percolator returns success (see per-file step 3 above).
Merged-path FDR:PSM is inherent — per-file results are already filtered before
`IDMergerAlgorithm::insertRuns()`, so no additional re-application is needed in
the merged block.

When Percolator fails for a given file, the warning log explicitly notes:
`"FDR:PSM=<X> was not applied for <file>; falling back to raw HyperScore."`
so users can detect the silent-fallback case from logs.

## Files touched

- `src/topp/ProSE.cpp` — main change (~150 lines reorganized; ~30 net additions
  for the new try/catch wrappers, FDR deferral block, and dangling-ref call).
- `CHANGELOG` — entry under `OpenMS 3.6.0 (under development) → TOPP tools →
  Changes → ProSE`.

No header changes. No new public API. `IDFilter::removeDanglingProteinReferences`
already exists at `IDFilter.h:787`.

## Testing

### What we add

- Topp-style test fixture for the merged-write isolation case if feasible:
  multi-input PQP+mzML where the merged path triggers `removeDanglingProteinReferences`
  (i.e., contains a peptide mapping to both target and decoy proteins with
  protein FDR enabled). Verify both per-file outputs and merged output exist
  on disk after the run. If constructing this fixture is too much yak-shaving,
  document a manual repro in the PR description.

- No new unit test for `removeDanglingProteinReferences` itself — already
  covered in `IDFilter_test.cpp`.

- FDR deferral re-application is hard to test without a real Percolator binary.
  Existing `TOPP_ProSE_*` tests skip Percolator; we don't add new automation
  here but document the deferral change in the PR description.

### What we keep

All existing `TOPP_ProSE_*` tests must still pass. The per-file write block
already uses per-mode try/catch; merged-write try/catch is additive. Pin
moving above Percolator changes pin output content only when `-percolator_executable`
is also set — verify any existing test that sets both flags simultaneously.

## Risks

- **Reordering output writes** could surface latent crashes in pin/qpx/parquet
  writers that previously failed AFTER idXML succeeded. Mitigation: per-mode
  try/catch isolates each writer's failure; the worst case is `input_failed=true`
  for that input, not aborting the whole run.

- **`removeDanglingProteinReferences` changes merged idXML output** when
  target+decoy PSMs were present. Existing test reference files may need
  refresh if any fixture had the bug latent — risk is low (test fixtures
  rarely include shared target+decoy peptides).

- **`FDR:PSM` deferral is a user-visible behavioral change.** With
  `-percolator_executable` set and `FDR:PSM > 0`, the cutoff now correctly
  applies post-Percolator; before this PR it stripped decoys pre-Percolator
  and silently fell back to raw HyperScore. CHANGELOG entry must mention this
  as a fix — users with `FDR:PSM > 0 + Percolator` who saw "rescoring failed"
  warnings will now see successful Percolator runs and a different (correct)
  output.

- **Pin moving above Percolator changes pin content** when both
  `-out_pin` and `-percolator_executable` are set. Today: pin contains
  Percolator-rescored PSMs. After: pin contains raw HyperScore PSMs (the
  intended "external rescoring input" semantic). Document in CHANGELOG.

- **`removeDanglingProteinReferences` default** uses
  `remove_peptides_without_reference = false`, so a `target+decoy` PSM whose
  only remaining evidence was the decoy gets an empty `getPeptideEvidences()`
  list rather than being deleted entirely. OpenMS handles this safely (no
  downstream crash). Choosing the default avoids silent PSM loss.

## Edge cases

- **All Percolator runs fail**: per-file results retain raw HyperScore;
  `FDR:PSM` filter is not applied (deferred and never re-applied because
  Percolator fallback). Logged in the per-file fallback warning so users
  notice. Merged path runs on the raw per-file results.

- **Some Percolator runs fail**: merged output mixes Percolator-q-value and
  raw HyperScore PSMs across input files. **Direct consequence of the
  deferral design** (raw-HyperScore files were never FDR-filtered). Tracked
  for follow-up in #9198 (internal PEP/FDR fallback when Percolator fails).
  This PR knowingly accepts this trade-off as part of fixing the silent
  decoy-stripping bug.

- **`has_decoys_for_percolator == false`** (Percolator skipped entirely
  before per-file loop): with `FDR:PSM > 0`, the deferred PSM FDR is
  silently never applied to ANY file. Raw HyperScore PSMs are written.
  Surfaced via a new warning at the decoy-check site:
  `"Percolator skipped (no decoys); deferred FDR:PSM=<X> was not applied
  to any file."` Long-term fix tracked in #9198.

- **`-out_pin` only, no `-percolator_executable`**: pin is written from raw
  result; no Percolator runs. Same as today minus the now-meaningless work of
  "writing pin AFTER Percolator that didn't run."

- **`-out_pin` + `-percolator_executable`**: pin from raw HyperScore (NEW —
  was Percolator-rescored before); idXML/qpx/parquet from rescored, FDR-filtered
  result. Both outputs derived from the same input but at different stages.

- **`-out_pin` only + `-percolator_executable` set + Percolator fails**:
  pin file is written from raw result (step 1, before the Percolator attempt).
  No idXML/qpx/parquet outputs are requested. User gets pin file as expected;
  the per-file Percolator fallback warning is logged but no actionable
  output is missing.

## Rollout

Single PR. CHANGELOG entry under `OpenMS 3.6.0 → TOPP tools → ProSE`:

```
- Fix: -out_merged write failure no longer discards per-file outputs.
  Per-file -out_idxml/-out_pin/-out_qpx/-out_parquet writes now complete before
  the merged write is attempted, and the merged write is wrapped in try/catch.
- Fix: merged idXML write no longer crashes with "Invalid protein reference
  'DECOY_...'" when target+decoy shared peptides are present and protein FDR
  is enabled. Added IDFilter::removeDanglingProteinReferences pass after
  decoy/protein cleanup.
- Fix: FDR:PSM is now correctly deferred when -percolator_executable is set
  (mirroring the existing FDR:protein deferral), then re-applied per-file
  using Percolator q-values. Previously FDR:PSM stripped decoys pre-Percolator
  and the rescoring silently failed with "No decoys found".
- BREAKING (minor): -out_pin now contains raw HyperScore PSMs even when
  -percolator_executable is set, matching the parameter's "external rescoring
  input" semantic. Use -out_idxml for Percolator-rescored output.
```
