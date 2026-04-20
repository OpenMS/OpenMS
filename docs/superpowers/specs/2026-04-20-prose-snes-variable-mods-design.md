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

- **`SpectrumMatch` struct** (line 71): add
  ```cpp
  uint32_t subset_bitmask_{};  ///< variable-mod subset for SNES v1.1+ (SNES hits) / 0 otherwise
  ```
  Grows struct from 16 B to 24 B with 8-byte padding. Per-spectrum overhead
  ≤ 400 B (at default `max_candidates_per_spectrum=50`). No external ABI
  consumers.

- **New private members**:
  ```cpp
  std::vector<double> snes_sigma_delta_set_;
  std::vector<double> snes_sigma_delta_set_with_prot_nterm_;
  ```
  Populated in `updateMembers_`. Always contain `0.0`. Always sorted ascending,
  distinct.

- **New helper**:
  ```cpp
  std::vector<double> computeSnesSigmaDeltaSet_(bool include_prot_nterm_mods) const;
  ```
  Returns sorted distinct Σ_delta values. Uses `modifications_variable_`,
  `max_variable_mods_per_peptide_`, per-mod ResidueModification specificity.
  Dedup tolerance: 1e-6 Da.

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

- **`updateMembers_`**: after SNES mode derivation (line 1836), populate both
  Σ_delta sets by calling `computeSnesSigmaDeltaSet_`. If `|set| > 64`, emit
  `OPENMS_LOG_WARN` recommending the user reduce variable mods or
  `max_per_peptide`. No hard reject.

- **`computeSnesSigmaDeltaSet_`** (new): iterate over
  `modifications:variable` entries, filter by specificity (include or
  exclude protein-N-term-only mods based on the parameter), enumerate all
  `m_count`-combinations where `m_count ≤ max_per_peptide`, compute Σ for
  each, deduplicate within 1e-6 Da, sort ascending. Always includes `0.0`.

- **`querySpectrumSNES_`**: restructure the per-(charge, iso_err) block into:
  ```
  for charge in charges:
    for iso_err in [min, max]:
      for sigma in snes_sigma_delta_set_:
        collect_candidates(shifted_mh - water - fixed_cterm_delta_ - sigma,
                           ..., expect_single_c=false, ...)
        collect_candidates(shifted_mh - fixed_nterm_delta_ - sigma,
                           ..., expect_single_c=true, ...)
      for sigma in snes_sigma_delta_set_with_prot_nterm_:
        # additional Single-N walks for protein-N-term mothers only
        ...
  ```
  `collect_candidates` annotates hits with the originating `sigma` so
  downstream subset enumeration knows the target Σ.

- **`querySpectrumSNES_`**: candidate collection is keyed by the pair
  `(mother_idx, Σ_delta)` — the same mother may appear once per Σ value,
  each treated as a separate rescore target. After collection, for each
  unique (mother_idx, Σ_delta):
  - Realize k via existing `realizeSNESLength(..., shifted_mh - Σ_delta, ...)`.
  - If k < 0 (no valid realization): skip.
  - Call `buildModSlots_(seq_ptr, k, slots, is_prot_nterm, is_prot_cterm)`
    where `seq_ptr` points at the realized sub-peptide and `is_prot_nterm`
    / `is_prot_cterm` reflect the mother's anchor in the protein.
  - Enumerate bitmask subsets over the returned slot list (bit `s` = slot
    `s` active). Accept subsets satisfying: (popcount ≤ max_per_peptide) ∧
    (no two active bits share `slots[s].position`) ∧
    (|Σ_subset − Σ_delta| < 1e-6 Da).
  - Apply per-tuple cap: ≤ 16 subsets emitted per (mother_idx, k, Σ_delta).
    If more valid subsets exist, emit the first 16 in bitmask numerical
    order (arbitrary but deterministic; log at debug level when triggered).
  - Each accepted subset → one `SpectrumMatch` with `subset_bitmask_` set
    to the slot-list bitmask. Other fields: `peptide_idx_ = mother_idx`,
    `isotope_error_`, `precursor_charge_` from the outer loop;
    `num_matched_` from the score_table entry at `mother_idx`.

- **min_matched_ions interaction**: Σ=0 candidates honor the v1
  `fragment:min_matched_ions` score_table pre-filter. Σ>0 candidates bypass
  it (rescored unconditionally). Comment explaining the tradeoff + pointer
  to B′ roadmap.

- **`trimHits(sms)`** (line 1583): unchanged; applies to the pooled hit list
  including multi-subset emissions.

### 4.3 `src/openms/source/ANALYSIS/ID/ProSEAlgorithm.cpp`

- At the SNES-hit call site for `reconstructRealizedSubSequence`
  (in `scoreSpectraAgainstIndex_`, the SNES branch), thread
  `sm.subset_bitmask_` into the extended call. The exact line number will
  shift with ongoing ProSE work; the call site is uniquely identified by
  the surrounding `is_snes_mode_` branch.
- **No changes to HyperScore / scoring logic**: it operates on the
  reconstructed AASequence. The modified residues are applied before
  scoring via the same code path as non-SNES peptides.
- Calibration disable warning (line 2301): unchanged.

### 4.4 `src/tests/class_tests/openms/source/FragmentIndex_test.cpp`

New sections after the current SNES block:

- **Σ_delta set computation**: given `modifications:variable =
  [Oxidation (M), Deamidation (NQ)]`, `max_per_peptide = 2` → verify set
  equals the hand-computed
  `{0, +0.984 (1 deamid), +1.968 (2 deamid), +15.995 (1 ox),
  +16.979 (1 ox + 1 deamid), +31.990 (2 ox)}`.
- **Query returns modified candidate (single-subset)**: synthesize a spectrum
  from `"AKACDEFMGR"` with `Oxidation (M)` applied at position 7, build SNES
  index with matching config, assert at least one candidate is returned
  with `subset_bitmask_ != 0` and reconstructs correctly.
- **Multiple subsets per Σ (emit-both)**: peptide with two M residues,
  `max_per_peptide = 1`, Σ=15.995 → assert two distinct SpectrumMatches
  with different `subset_bitmask_` values.
- **Position conflict rejection**: two variable mods claiming the same
  residue → activating both produces no SpectrumMatch.
- **max_mods cap**: three eligible mod sites, `max_per_peptide = 1` →
  only single-mod subsets emitted.
- **Protein-N-term specificity**: `Acetyl (Protein N-term)` appears in
  Σ_delta set only for Single-N mothers at protein position 0.
- **Regression**: existing PR #9189 review tests (Acetyl fixed N-term
  fixed-mod, b/y rejection) still pass unchanged.

## 5. Data Flow Example

Input: spectrum from peptide `"ACDEFMGR"` with `Oxidation (M)` at position 6.
Observed precursor (M+H)+ = unmodified mass + 15.995 + proton.
Configuration: `modifications:variable = [Oxidation (M)]`,
`max_per_peptide = 1`.

1. `updateMembers_` populates `snes_sigma_delta_set_ = {0, 15.995}`.
2. `querySpectrumSNES_` at (charge=1, iso_err=0):
   - Outer loop σ=0: bin walk at shifted_mh − water − fixed_cterm (Single-N)
     and shifted_mh − fixed_nterm (Single-C). No hit (observed is shifted
     by 15.995).
   - Outer loop σ=15.995: bin walks at `shifted_mh − water − fixed_cterm −
     15.995` and `shifted_mh − fixed_nterm − 15.995`. Hit returned: Single-N
     mother containing "ACDEFMGR".
3. Per-candidate enumeration for `(mother, Σ=15.995)`:
   - Realize k → k=8 (full sub-peptide "ACDEFMGR").
   - Build slot list on positions `[0, 8)` → one slot: M@position_5
     (relative to sub-peptide start).
   - Enumerate subsets with Σ = 15.995: `{1u << 0}` (activate the M slot).
   - Emit one `SpectrumMatch` with `subset_bitmask_ = 1u << 0`.
4. ProSE's `reconstructRealizedSubSequence(..., subset_bitmask=1)` applies
   `setModification(5, "Oxidation")` to the AASequence.
5. HyperScore scores the modified theoretical spectrum against the observed.

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

**vs SNES v1 no-mods**: ~5–30× slower per query. Two stacked multipliers:
D ≈ 5–15 Σ values × K ≈ 2–5 subsets per candidate at rescore.

**vs non-SNES v1 (same mod config)**: roughly comparable query time,
2–3× faster build, 5–10× smaller index. The SNES index-size win is
preserved.

These are algorithmic extrapolations, not measured. Benchmarking as part of
the implementation plan will validate.

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
| Degenerate subset explosion (many eligible residues + high max) | Per-tuple cap ≤ 16 subsets; global `trimHits(sms)` cap unchanged. |

## 9. Open Items

- **Benchmark validation**: the performance estimates in §6 are algorithmic
  extrapolations. Initial implementation should include a micro-benchmark
  (single spectrum, varying D and candidate count) and a small-proteome
  end-to-end benchmark before tuning the per-tuple subset cap.
- **Per-tuple cap value**: 16 is a conservative guess. Actual value may
  need tuning based on benchmark observations.
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
