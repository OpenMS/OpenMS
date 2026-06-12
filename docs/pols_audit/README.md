# OpenMS POLS (Principle of Least Surprise) API Audit — Master Report

A module-by-module audit of the OpenMS public API for places where a declaration's name/signature surprises a competent caller. 
Every finding was produced by a header-cluster finder agent, then **independently re-verified against the actual source by a skeptical agent** (default = reject), which also re-graded severity/category/ABI. False positives and domain conventions were filtered out.

## Headline

- **827 verified findings** across all 18 OpenMS modules (~787 public headers, 102 header clusters).
- Severity: **62 high · 359 medium · 406 low**.
- ABI impact of the recommended fix: **528 need no ABI change**, **227 source-compatible** (deprecate-and-alias / add overload / mark const), only **72 truly breaking**. → ~91% are fixable without breaking ABI.
- **257 high/medium findings are ABI-neutral** — the highest-value quick wins (see `all_findings.csv`, filter abi=none).

## Per-module breakdown

| Batch | Scope | Findings | High | Med | Low | Report |
|---|---|--:|--:|--:|--:|---|
| B1 | CONCEPT/DATASTRUCTURES/KERNEL | 133 | 3 | 49 | 81 | [batch1_concept_datastructures_kernel.md](batch1_concept_datastructures_kernel.md) |
| B2 | METADATA/CHEMISTRY | 124 | 4 | 69 | 51 | [batch2_metadata_chemistry.md](batch2_metadata_chemistry.md) |
| B3 | FORMAT | 150 | 14 | 60 | 76 | [batch3_format.md](batch3_format.md) |
| B4a | ANALYSIS/ID | 73 | 5 | 30 | 38 | [batch4a_analysis_id.md](batch4a_analysis_id.md) |
| B4b | ANALYSIS/OPENSWATH+TARGETED | 82 | 7 | 39 | 36 | [batch4b_openswath.md](batch4b_openswath.md) |
| B4c | ANALYSIS/MAPMATCHING/QUANT/NUXL/TOPDOWN | 69 | 8 | 29 | 32 | [batch4c_analysis_rest.md](batch4c_analysis_rest.md) |
| B5 | FEATUREFINDER/PROCESSING | 89 | 9 | 35 | 45 | [batch5_featurefinder_processing.md](batch5_featurefinder_processing.md) |
| B6 | MATH/ML/COMPARISON | 64 | 8 | 33 | 23 | [batch6_math_ml_comparison.md](batch6_math_ml_comparison.md) |
| B7 | QC/SYSTEM/APPLICATIONS/IM/IF/IMG | 43 | 4 | 15 | 24 | [batch7_qc_system_apps.md](batch7_qc_system_apps.md) |
| **Total** | **18 modules** | **827** | **62** | **359** | **406** | [ALL_HIGH_SEVERITY.md](ALL_HIGH_SEVERITY.md) · [all_findings.csv](all_findings.csv) |

## Finding categories (whole codebase)

| Count | Category |
|--:|---|
| 150 | silent-failure |
| 101 | misleading-name |
| 95 | hidden-side-effect |
| 64 | inconsistent-convention |
| 61 | const-correctness |
| 50 | param-order-or-bool |
| 43 | return-value |
| 41 | surprising-throw |
| 34 | asymmetric-api |
| 28 | surprising-default |
| 16 | unit-or-index |
| 14 | ownership-lifetime |
| 12 | other |
| 12 | misleading-doc |
| 11 | implicit-conversion |
| 9 | documentation |

## The systemic patterns (fix the pattern, not just the instance)

These recur across many modules; each is worth a focused, codebase-wide pass:

1. **Mutating methods named like queries/getters.** Non-const `get*`/`is*` that insert-on-miss via `map_[key]`, lazily recompute, re-sort, or rewrite scores in place. _e.g._ `MRMFeature::getFeature`, `FalseDiscoveryRate::apply`, `QTCluster::getAnnotations`, `TargetedExperiment::get*ByRef`, `MultiplexDeltaMassesGenerator::getLabelShort`, `SpectrumCheapDPCorr::operator()` (const but mutates), `Math::quantile*`/`median` (sort the caller's range).
2. **`load()`/import: clear-vs-append is unspecified and inconsistent.** Reusing an output object across two loads silently accumulates. _e.g._ `MzIdentMLFile`, `MSPFile`, `OMSFile`, `ParamXMLFile`, `InspectOutfile`, `ConsensusMapArrowIO`, `MzMLSqliteHandler`, `Ms2IdentificationRate::compute`. → define the contract on the base class and audit every loader.
3. **Silent-failure sentinels instead of the codebase's own loud-throw idiom.** empty/`-1`/`end()`/`NaN`/`0.0`/perfect-`1.0`/`0`-score on error. _e.g._ `LPWrapper::getRowIndex`, `getDataArrayByName`, `HyperScore`/`MorpheusScore`, `PrecursorPurity`, store() on unwritable path (GNPS*), `SwathWindowLoader`, `MRMScoring` returns 0='perfect coelution' on degenerate input, `File::fileSize` (-1 through unsigned).
4. **`==` by pointer identity while `<`/hash use value; comparators that aren't strict-weak-orderings.** → silent corruption / UB in std::map/set/sort. _e.g._ `ModificationDefinition`, `MetaInfoDescription`, `GridBasedCluster::operator<` (y only), `FuzzyDoubleComparator`, `OPXLHelper::PeptideIDScoreComparator`.
5. **`@param[in]` that are really outputs, `[out]` that append, and out-params the callee doesn't size.** Plus **bare same-typed bool/numeric params** trivially swapped at call sites (XICParquetFile's six `Int64=-1`, NuXLDeisotoper's 16 params).
6. **Documented `@throw`/contracts that never fire** (assert-only guards → release-build UB or garbage). _e.g._ the Binned comparators' `IncompatibleBinning`, `Matrix`/`DistanceMatrix::operator==`, `RangeManager::getRangeForDim`, `QCBase::names_of_requires` OOB.
7. **Hidden global-singleton mutation as a side effect** of a read/generate call. _e.g._ `NuXLPresets::getPresets` (mutates ResidueDB), `EmpiricalFormula` proton folding, exception construction writing process-global state.
8. **Unit/dimension erosion** — m/z stored under 'mass' names, ppm under generic names, tolerances silently halved/doubled, `RangeBase` implicitly convertible across dimensions, IM 1/K0 vs ms.

## A few genuine bugs (not just smells) surfaced

- `MzTabBoolean::setNull(bool)` — **inverted polarity** (`setNull(true)` makes it not-null).
- `Math::absdev` — sums **signed** deviations; a function named 'absolute deviation' returns ≈0.
- `QCBase::names_of_requires[]` — missing the `Requires::ID` entry → **out-of-bounds read**.
- `IDFilter::filterHitsByScore` — doc says removes hits when score type missing; silently **keeps** them.
- `NuXLFDR::splitIntoPeptidesAndXLs` — keeps the single first hit overall, not first-of-each-class as documented (wrong FDR split).
- `TIC::getResults` — declared, **never defined**.

## How to use this

- `ALL_HIGH_SEVERITY.md` — the 62 highs, full detail, grouped by module. Start here.
- `all_findings.csv` — every finding (id, module, file, symbol, severity, category, ABI, fix). Sort by `abi=none` + `severity` for the quick-win backlog.
- `batch*.md` — full per-module reports with evidence + the verifier's reasoning for every finding.

## Method notes / caveats

- 102 finder clusters → adversarial per-finding verification (skeptical, default-reject) → per-module synthesis. ~1500+ agents total.
- Severity/category/ABI shown are **post-verification** (the verifier frequently down-graded finder over-claims; those corrections are preserved in the per-batch reports).
- Running two workflows concurrently tripped server-side rate limiting; all batches were ultimately (re-)run **solo** to ensure complete coverage.
- This is an API-surprise audit, not a full bug hunt or security review; low-severity items are largely doc/naming papercuts. Recommendations favour ABI-safe fixes (deprecate-and-alias, add overload, mark const/explicit, fix docs).