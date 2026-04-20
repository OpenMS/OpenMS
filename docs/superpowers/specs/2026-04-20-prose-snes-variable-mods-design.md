# SNES Variable-Modification Support (v1.1) — Design Spec

**Date**: 2026-04-20
**Author**: Timo Sachsenberg
**Status**: Design approved, pending implementation plan
**Target**: ProSE SNES v1.1 (follow-up to PR #9189)

## 1. Context & Motivation

PR #9189 shipped SNES v1: the ProSE FragmentIndex can be built in a mode where
mother peptides (the longest peptide anchored at one terminus of each protein)
are indexed instead of the full O(L²) set of non-specific sub-peptides. This
preserves index compactness and query speed for immunopeptidomics-style
searches with `peptide:enzyme_specificity=none`.

SNES v1 has a hard restriction: **no variable modifications are supported**.
Users running non-specific searches must forgo `modifications:variable` when
opting into SNES. This blocks common immunopeptidomics workflows that rely on
`Oxidation (M)`, `Deamidation (NQ)`, or N-term PTM labels.

v1.1 adds variable-mod support at the query path without sacrificing the v1
index-size / build-time advantages.

## 2. Scope

### In scope
- Query-time enumeration of variable modifications on realized SNES
  sub-peptides (not on mothers at build time).
- Full parity with non-SNES variable-mod semantics: `modifications:variable`,
  `modifications:variable_max_per_peptide`, ANYWHERE / N-term / C-term /
  protein-N-term / protein-C-term specificity, position conflicts, per-AA
  applicability.
- Σ_delta enumeration (MetaMorpheus-style delta-bag precursor filter).
- Emitting multiple SpectrumMatches per (mother, realized-k, Σ_delta) when
  multiple subsets satisfy the Σ constraint.

### Out of scope (roadmap)
- Shift-aware fragment-match filter (Approach B′) — v1.2 if benchmarks show
  A's rescore load dominates.
- X/B/Z ambiguous-residue span splitting (CodeRabbit #1) — tracked
  separately on the SNES roadmap.
- Variable modifications stored in the fragment index (Approach B from
  brainstorm) — rejected due to index-size explosion.
- Open-search semantics (`search:open_search_mode`) — v1 already disables
  this; unchanged in v1.1.

## 3. Architecture Overview

Build time: **unchanged**. Mothers are indexed once, unmodified. Σ_delta
enumeration is not materialized in the index.

Query time: two new enumeration passes layered on top of v1:

1. **Precursor-filter expansion** (`querySpectrumSNES_`)
   - Precompute `snes_sigma_delta_set_` at `updateMembers_` time: distinct
     Σ values reachable by any valid subset of variable mods capped by
     `max_variable_mods_per_peptide`.
   - The existing two per-(charge, iso_err) bin walks (Single-N target,
     Single-C target) expand into `2 × |snes_sigma_delta_set_|` walks, each
     at a Σ-shifted target.
   - Cost multiplier: linear in D = |Σ_delta_set|. Typical D = 5–15.

2. **Per-candidate subset enumeration** (`querySpectrumSNES_`)
   - Each bin-walk hit returns `(mother_idx, Σ_delta)` — Σ_delta is known
     from which outer-loop iteration produced the hit.
   - Realize k via `realizeSNESLength(..., shifted_mh − Σ_delta, ...)`.
   - Call `buildModSlots_` scoped to the realized sub-peptide's position
     range (Single-N: mother positions `[0, k)`; Single-C: `[L−k, L)`).
   - Enumerate bitmask subsets satisfying: (popcount ≤ max_per_peptide)
     ∧ (no two bits share a residue position) ∧ (Σ over the subset ≈
     Σ_delta within tolerance).
   - Each valid subset → one `SpectrumMatch` with `subset_bitmask_` set.

Storage: `mod_bitmask_` in `Peptide` entries remains kind-bit-only in SNES
mode (bit 31). Variable mod subsets are transient per-query state, not
materialized in the index.

## 4. Components

### 4.1 `src/openms/include/OpenMS/ANALYSIS/ID/FragmentIndex.h`

- **`SpectrumMatch` struct** (line 71): add two new fields. The bitmask
  encodes which slots of `buildModSlots_`'s slot list are active for this
  match — it indexes into the slot list, **not** raw peptide positions.
  Example: peptide with one eligible mod site → slot list has one entry →
  activated bitmask is `1u << 0 = 1`, regardless of where that residue
  sits in the peptide. `sigma_delta_` is the Σ of variable-mod deltas
  that led to this match; ProSE needs it to subtract from the realization
  target. Layout (both new fields placed to avoid padding growth):
  ```cpp
  uint32_t num_matched_{};       // 4
  uint32_t subset_bitmask_{};    // 4 (NEW)
  float    sigma_delta_{};       // 4 (NEW) — Σ of variable-mod deltas for this match
  uint16_t precursor_charge_{};  // 2
  int16_t  isotope_error_{};     // 2
  size_t   peptide_idx_{};       // 8
  ```
  Total size: 4+4+4+2+2+8 = 24 B (natural 8-byte alignment; no padding).
  Grows struct from 16 B to 24 B; per-spectrum overhead ≤ 400 B at default
  `max_candidates_per_spectrum=50`. No external ABI consumers.

- **New private members**:
  ```cpp
  std::vector<double> snes_sigma_delta_set_;                   // ANYWHERE + peptide-term mods only
  std::vector<double> snes_sigma_delta_set_with_prot_nterm_;   // above + PROTEIN_N_TERM mods
  std::vector<double> snes_sigma_delta_set_with_prot_cterm_;   // above + PROTEIN_C_TERM mods
  ```
  Populated in `updateMembers_`. Always contain `0.0`. Always sorted ascending,
  distinct. The two protein-term sets apply only to mothers anchored at the
  corresponding protein terminus (Single-N @ protein pos 0, Single-C @
  protein end). The sets are configuration-global — they enumerate every
  Σ reachable from the variable-mod configuration capped by
  `max_variable_mods_per_peptide`, without regard to any specific peptide's
  residue inventory. Per-peptide applicability (does this subset actually
  fit on this realized k-mer?) is checked later in the subset-enumeration
  step. Many bin-walk hits will produce zero valid subsets on their
  realized k-mer — this is expected, not an error.

- **New helper**:
  ```cpp
  std::vector<double> computeSnesSigmaDeltaSet_(bool include_prot_nterm_mods,
                                                bool include_prot_cterm_mods) const;
  ```
  Returns sorted distinct Σ_delta values. Uses `modifications_variable_`,
  `max_variable_mods_per_peptide_`, per-mod ResidueModification specificity.
  Dedup tolerance: 1e-6 Da (chosen to absorb compounded FP error across up
  to ~16 summed deltas in double precision, not just database duplicates).

- **Extend `reconstructRealizedSubSequence`**: add a `uint32_t subset_bitmask`
  parameter. Applies the corresponding mods via `setModification` /
  `setNTerminalModification` / `setCTerminalModification` during AASequence
  construction. Old signature remains as a thin wrapper delegating with
  `subset_bitmask=0`.

- **Defensive masking in `reconstructModifiedSequence`** (line 416): replace
  `if (peptide.mod_bitmask_ != 0)` with a loop over
  `peptide.mod_bitmask_ & SNES_SLOT_MASK` to ensure the kind bit never
  accidentally drives slot iteration. No-op for non-SNES callers.

### 4.2 `src/openms/source/ANALYSIS/ID/FragmentIndex.cpp`

- **`updateMembers_`**: after SNES mode derivation (line 1836), populate the
  three Σ_delta sets via three calls to `computeSnesSigmaDeltaSet_`. If the
  largest set has > 64 entries, emit `OPENMS_LOG_WARN` recommending the user
  reduce variable mods or `max_per_peptide`. No hard reject.

- **`computeSnesSigmaDeltaSet_`** (new): iterate over
  `modifications:variable` entries, filter by specificity (include or
  exclude PROTEIN_N_TERM-only and PROTEIN_C_TERM-only mods based on the
  parameters), enumerate all `m_count`-combinations where
  `m_count ≤ max_per_peptide`, compute Σ for each, deduplicate within
  1e-6 Da, sort ascending. Always includes `0.0`.

- **`querySpectrumSNES_`**: restructure the per-(charge, iso_err) block into:
  ```
  for charge in charges:
    for iso_err in [min, max]:
      compute shifted_mh = mh_plus + iso_err * C13C12
      // Baseline walks: hit all mothers regardless of protein anchor.
      for sigma in snes_sigma_delta_set_:
        collect_candidates(shifted_mh - water - fixed_cterm_delta_ - sigma,
                           ..., expect_single_c=false,
                           require_prot_anchor=NONE, sigma_tag=sigma)
        collect_candidates(shifted_mh - fixed_nterm_delta_ - sigma,
                           ..., expect_single_c=true,
                           require_prot_anchor=NONE, sigma_tag=sigma)
      // Additional walks for protein-anchored mothers; gated by require_prot_anchor.
      for sigma in (snes_sigma_delta_set_with_prot_nterm_ \ snes_sigma_delta_set_):
        collect_candidates(shifted_mh - water - fixed_cterm_delta_ - sigma,
                           ..., expect_single_c=false,
                           require_prot_anchor=PROT_NTERM, sigma_tag=sigma)
      for sigma in (snes_sigma_delta_set_with_prot_cterm_ \ snes_sigma_delta_set_):
        collect_candidates(shifted_mh - fixed_nterm_delta_ - sigma,
                           ..., expect_single_c=true,
                           require_prot_anchor=PROT_CTERM, sigma_tag=sigma)
  ```
  **Extension to `collect_candidates`**: gains two new parameters:
  `require_prot_anchor` (enum: NONE / PROT_NTERM / PROT_CTERM) and
  `sigma_tag` (double). The former adds a post-kind-check filter:
  `if (require_prot_anchor == PROT_NTERM && mother.sequence_.first != 0) continue;`
  (symmetric check for PROT_CTERM). The latter stores σ on each emitted
  `SpectrumMatch` so downstream subset enumeration knows the target Σ. The
  set-difference notation `A \ B` means "only values in A but not B" to
  avoid double-counting the shared baseline Σ values.

- **`querySpectrumSNES_`**: candidate collection is keyed by the pair
  `(mother_idx, Σ_delta)` — the same mother may appear once per Σ value,
  each treated as a separate rescore target. After collection, for each
  unique (mother_idx, Σ_delta):
  - Realize k via `realizeSNESLength(mother, entries, shifted_mh − Σ_delta,
    tol, ppm)`. Important: `shifted_mh` here is the same iso-corrected
    observed (M+H)+ computed at the top of the (charge, iso_err) block
    (line 1505-1506 of the v1 code), **not** the fragment-bin target m/z.
    This is a distinct arithmetic step from the bin-walk subtraction,
    even though both subtract Σ_delta. `realizeSNESLength` compares
    `target_mh_plus` against the *unmodified* realized mass (water + proton
    + fixed_nterm + fixed_cterm + Σ internal residue masses + Σ internal
    fixed mod deltas — see v1 code lines 461-463 and 482-483), so
    pre-subtracting Σ_delta before the call tells it "find the k that
    would yield this mass after removing the variable-mod contribution".
  - If k < 0 (no valid realization): skip.
  - Call `buildModSlots_(seq_ptr, k, slots, is_prot_nterm, is_prot_cterm)`
    where `seq_ptr` points at the realized sub-peptide and `is_prot_nterm`
    / `is_prot_cterm` reflect the mother's anchor in the protein.
  - Enumerate bitmask subsets over the returned slot list (bit `s` = slot
    `s` active). Accept subsets satisfying: (popcount ≤ max_per_peptide) ∧
    (no two active bits share `slots[s].position`) ∧
    (|Σ_subset − Σ_delta| < 1e-6 Da).
  - Apply **per-mother cap**: ≤ 16 subsets total emitted per mother across
    all (k, Σ_delta) tuples for that mother in a single query. Rationale:
    capping per-tuple would interact poorly with `trimHits`, which ranks
    by `num_matched_` — subsets of the same mother share that value, so
    tuple-scoped caps deterministically drop entire mother-subset groups
    rather than diluting across mothers. Per-mother capping distributes
    pressure evenly. If more than 16 valid subsets exist across the
    mother's (k, Σ) tuples, emit the first 16 in visitation order
    (deterministic; log at debug level when triggered).
  - Each accepted subset → one `SpectrumMatch` with `subset_bitmask_` set
    to the slot-list bitmask. Other fields: `peptide_idx_ = mother_idx`,
    `isotope_error_`, `precursor_charge_` from the outer loop;
    `num_matched_` from the score_table entry at `mother_idx`.

- **min_matched_ions pre-filter**: applied uniformly to all candidates
  regardless of Σ_delta. The score_table counts matches against the
  mother's *unmodified* indexed b/y fragments, so for Σ>0 candidates the
  count systematically under-reports true matches. This is the accepted
  tradeoff for v1.1: modified candidates with weak unmodified-fragment
  support will be filtered before rescore. The roadmap item B′
  (shift-aware fragment filter) is the lever for recovering those
  candidates if benchmarks show the filter is dropping true positives.
  Add an in-code comment documenting this tradeoff and the roadmap link.

- **`trimHits(sms)`** (line 1583): unchanged; applies to the pooled hit list
  including multi-subset emissions.

### 4.3 `src/openms/source/ANALYSIS/ID/ProSEAlgorithm.cpp`

- At the SNES-hit call site for `reconstructRealizedSubSequence`
  (ProSEAlgorithm.cpp lines 917-920, branch gated by `snes_mode`),
  thread `sm.subset_bitmask_` and `sm.sigma_delta_` through. The existing
  site computes
  ```cpp
  iso_shifted_target = observed_mh_plus + isotope_error * c13c12;
  realizeSNESLength(sms_pep, db, iso_shifted_target, ...);
  ```
  v1.1 subtracts `sm.sigma_delta_` from the target:
  ```cpp
  iso_shifted_target = observed_mh_plus + isotope_error * c13c12
                       - sm.sigma_delta_;
  realizeSNESLength(sms_pep, db, iso_shifted_target, ...);
  reconstructRealizedSubSequence(sms_pep, db, realized_len, sm.subset_bitmask_);
  ```
  The realization target uses the subtracted value so it lands on the
  unmodified realized mass, consistent with `realizeSNESLength`'s internal
  algebra. The subset_bitmask_ parameter is forwarded to the reconstruction
  so variable mods are applied during AASequence construction.
- **No changes to HyperScore / scoring logic**: it operates on the
  reconstructed AASequence. The modified residues are applied before
  scoring via the same code path as non-SNES peptides.
- Calibration disable warning (line 2301): unchanged.

### 4.4 `src/tests/class_tests/openms/source/FragmentIndex_test.cpp`

New sections after the current SNES block:

- **Σ_delta set computation (values)**: given `modifications:variable =
  [Oxidation (M), Deamidation (NQ)]`, `max_per_peptide = 2` → verify set
  contains, in ascending order, the six distinct values
  `{0.0, 0.984016 (1 deamid), 1.968032 (2 deamid), 15.994915 (1 ox),
  16.978931 (1 ox + 1 deamid), 31.989830 (2 ox)}`.
  Use `TEST_REAL_SIMILAR` with 1e-4 Da tolerance (absolute values from
  Unimod; delta precision varies by DB version).
- **Query returns modified candidate (single-subset)**: synthesize a
  spectrum from `"ACDEFMGR"` with `Oxidation (M)` applied at the M residue
  (position index 5, 0-based), build SNES index with matching config,
  assert at least one candidate is returned with `subset_bitmask_ != 0`
  and reconstructs correctly. (Canonical example peptide — same used in §5.)
- **Multiple subsets per Σ (emit-both)**: peptide with two M residues,
  `max_per_peptide = 1`, Σ=15.995 → assert two distinct SpectrumMatches
  with different `subset_bitmask_` values.
- **Position conflict rejection**: two variable mods claiming the same
  residue → activating both produces no SpectrumMatch.
- **max_mods cap**: three eligible mod sites, `max_per_peptide = 1` →
  only single-mod subsets emitted.
- **Protein-N-term Σ-set computation (unit)**: call
  `computeSnesSigmaDeltaSet_` with `include_prot_nterm_mods=false` →
  assert Acetyl delta absent; call again with `true` → assert Acetyl delta
  present. Symmetric test for `include_prot_cterm_mods` with e.g. Amidation
  (Protein C-term).
- **Protein-N-term specificity (query-path)**: configure
  `Acetyl (Protein N-term)` and build an index with two proteins: one
  where a Single-N mother anchors at protein position 0, one where the
  Single-N mother is mid-protein. Query a spectrum matching the Σ that
  requires the Acetyl: only the protein-anchored mother admits the
  candidate; the mid-protein mother emits no SpectrumMatch for that Σ.
- **Protein-C-term specificity (query-path)**: symmetric test for
  `Amidated (Protein C-term)` with Single-C mothers.
- **Identical-delta mods**: two variable mods with identical delta (e.g.,
  two synthetic mods both at +15.995 on different residues) produce two
  distinct `subset_bitmask_` values at the same Σ (regression guard
  against dedup-tolerance swallowing legitimate subsets).
- **Regression**: existing PR #9189 review tests (Acetyl fixed N-term
  fixed-mod, b/y rejection) still pass unchanged.

## 5. Data Flow Example

Input: spectrum from peptide `"ACDEFMGR"` (an 8-mer) with `Oxidation (M)`
at the M residue (position index 5, 0-based). Observed (M+H)+ equals the
unmodified (M+H)+ of `"ACDEFMGR"` shifted by +15.995 Da (the Oxidation delta).
Configuration: `modifications:variable = [Oxidation (M)]`,
`max_per_peptide = 1`, `precursor:isotope_error_max = 0`.

Notation convention for this example: `shifted_mh` is
`mh_plus + iso_err × C13C12`, matching the v1 code (line 1505-1506).
At iso_err = 0, `shifted_mh = mh_plus = observed (M+H)+`.

1. `updateMembers_` populates `snes_sigma_delta_set_ = {0, 15.995}`.
2. `querySpectrumSNES_` at (charge=1, iso_err=0):
   - Outer σ=0: bin walks at `shifted_mh − water − fixed_cterm` (Single-N)
     and `shifted_mh − fixed_nterm` (Single-C). No hit — the mother's
     indexed b_k/y_k are unmodified, but the observed shifted_mh carries
     the +15.995 from the Oxidation, so these targets land ~16 Da off the
     indexed values.
   - Outer σ=15.995: bin walks at `shifted_mh − water − fixed_cterm −
     15.995` and `shifted_mh − fixed_nterm − 15.995`. The Single-N target
     now equals the indexed b_8 of a mother containing "ACDEFMGR" as a
     prefix. Hit returned: `(mother_idx, Σ=15.995)`.
3. Per-candidate enumeration for `(mother, Σ=15.995)`:
   - Call `realizeSNESLength(mother, entries, shifted_mh − 15.995, tol, ppm)`.
     The function compares the argument against the unmodified realized
     mass and returns `realized_len = 8` (full sub-peptide `"ACDEFMGR"`).
   - `buildModSlots_` on the 8-mer at positions `[0, 8)` returns one slot:
     slot 0 = Oxidation on the M residue at sub-peptide position 5.
   - Subset enumeration: bitmasks 0b0 (Σ=0, skip — doesn't match Σ=15.995)
     and 0b1 (Σ=15.995, accept).
   - Emit one `SpectrumMatch` with `subset_bitmask_ = 1`, `sigma_delta_ =
     15.995`, `peptide_idx_ = mother_idx`.
4. ProSE picks up the match and computes
   `iso_shifted_target = observed_mh_plus + 0 × c13c12 − 15.995`.
   Calls `realizeSNESLength(...)` → 8 (same result as the FI-side call).
   Calls `reconstructRealizedSubSequence(mother, db, 8, subset_bitmask=1)`.
   The function reconstructs `"ACDEFMGR"` and applies
   `setModification(5, "Oxidation")` (slot 0 → sub-peptide position 5).
5. HyperScore scores the modified theoretical spectrum against the observed.

Note: the Σ subtraction happens in two separate arithmetic contexts — the
bin-walk target in step 2 and the realization target in steps 3 and 4.
Both subtract the same Σ value from `shifted_mh` but the results flow into
different code paths. Do not conflate them.

## 6. Performance Expectations

Rough estimates for typical immunopeptidomics config (8–12 AA,
`Oxidation (M)` + `Deamidation (NQ)`, `max_per_peptide = 2`, 10 ppm prec,
human proteome):

| Metric | Non-SNES v1 | SNES v1 (no mods) | SNES v1.1 M1 |
|---|---|---|---|
| Build time | 1–5 min | 1–2 min | 1–2 min |
| Index size | ×5–10 | baseline | baseline |
| Per-spectrum query | ~1–3 ms | ~0.2 ms | ~2–5 ms |
| 50k-spec run (16 thr) | 3–10 min | 30 sec | 5–10 min |

**vs SNES v1 no-mods**: multiplicative slowdown D × K. For common configs:
- D = 1, K = 1 (no variable mods): collapses exactly to v1. Zero overhead.
- D = 6, K = 2 (Oxidation M + Deamidation NQ, max=2): roughly 12× per query.
- D = 15, K = 5 (richer config, many eligible residues on long peptides):
  up to 30–75× per query. This is the pathological end — the 64-entry
  Σ-set warning (§4.2) should fire well before this.

**vs non-SNES v1 (same mod config)**: same *order of magnitude* per-query
time — both paths pay the mod cost somewhere. Non-SNES pays at build
(pre-expanded peptide entries) and enjoys direct fragment-index lookups
at query. SNES v1.1 pays D bin walks at query but builds faster and keeps
a smaller index. The relative performance depends on workload — narrow
precursor tolerance favors SNES (D × narrow ≪ non-SNES's wider candidate
admission); wide tolerance / high iso_err favors non-SNES. **No
trustworthy general claim without benchmarks.**

These are algorithmic extrapolations. Benchmarking as part of the
implementation plan must validate the D and K estimates on a realistic
immunopeptidomics workload before the code path ships as default-safe.

## 7. Backward Compatibility & Non-SNES Impact

The M1 implementation is **additive** to SNES-only code paths. Non-SNES
tryptic / semi-tryptic searches are unaffected:

- `Peptide` struct: same 32-bit `mod_bitmask_`, same meaning.
- `buildModSlots_`: signature + behavior unchanged.
- `reconstructModifiedSequence`: gains a defensive `& SNES_SLOT_MASK`, which
  is a no-op for non-SNES callers (bit 31 is never set there).
- `MAX_MOD_SLOTS = 32`: non-SNES still uses all 32 bits.
- `generatePeptides`, `querySpectrum` (non-SNES): untouched.
- `modifications:variable`, `variable_max_per_peptide`: same parameters,
  same semantics, same output for non-SNES.

User-visible changes in non-SNES: **zero**.

## 8. Failure Modes

| Condition | Behavior |
|---|---|
| `modifications:variable=[]` + `snes_enabled=true` | Σ_set = {0}, D=1, identical to SNES v1 path. No overhead. |
| `|Σ_set| > 64` | `OPENMS_LOG_WARN` recommending config reduction. No reject; performance degrades. |
| Variable mod position conflicts with fixed terminal mod | Handled by existing `buildModSlots_` conflict path. |
| Σ from bin walk has no matching subset on realized k-mer | Candidate silently discarded (no SpectrumMatch emitted). Common; not an error. |
| Degenerate subset explosion (many eligible residues + high max) | Per-mother cap ≤ 16 subsets across all (k, Σ) tuples; global `trimHits(sms)` cap unchanged. |
| Σ from a protein-term-only mod on a mid-protein mother | Σ not in the baseline set; additional protein-term walks gated by `require_prot_anchor` skip this mother. No ghost candidate emitted. |

## 9. Open Items

- **Benchmark validation**: the performance estimates in §6 are algorithmic
  extrapolations. Initial implementation should include a micro-benchmark
  (single spectrum, varying D and candidate count) and a small-proteome
  end-to-end benchmark before tuning the per-tuple subset cap.
- **Per-mother cap value**: 16 is a conservative guess. Actual value may
  need tuning based on benchmark observations.
- **Implementation split**: a reviewer raised the option of landing this
  in two PRs (infrastructure + Σ-set helper + defensive masking first,
  then query path + ProSE integration). The combined spec reads as one
  PR; the decision on split can happen at plan-writing time.
- **Σ_delta dedup tolerance**: 1e-6 Da is tight. For mods with very close
  but distinct deltas (rare), may collapse to one. Document the tolerance
  in the parameter docstring.

## 10. Roadmap (v1.2+)

- **Approach B′** (shift-aware fragment filter): if A's rescore load
  dominates for realistic workloads, B′ adds a shift-aware score_table
  pre-filter. Estimated 2–3× query speedup for mod-heavy searches at the
  cost of D_prefix-multiplied per-peak cost.
- **X/B/Z span splitting** (CodeRabbit #1): separate from variable-mod work;
  reduces false negatives in proteins with ambiguous residues.
- **Precursor-delta open searches**: MSFragger-style open search via a
  continuous Σ window instead of a discrete set. Requires different query
  architecture; probably v2.
