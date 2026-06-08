# OpenMS 3.6.0 — Pre-release Testing Plan

> Derived from the `CHANGELOG` "OpenMS 3.6.0 (under development)" section (lines 16–531) and the
> in-repo test/CI infrastructure. Prioritized by **blast radius × risk**, with a dedicated
> ProteoBench end-to-end strategy.

## Summary

OpenMS 3.6.0 is a large, infrastructure-heavy release. The dominant themes are:

1. Removing Qt from the core library (std::filesystem / std::chrono / boost::process / libcurl / nlohmann-json / libzip).
2. A complete pyOpenMS rewrite (Autowrap/Cython → nanobind, ~544 classes).
3. Native Parquet I/O (`.idparquet` / `.featureparquet` / `.consensusparquet`; Arrow now mandatory).
4. Native vendor-format readers (Bruker `.d` / Thermo `.raw`).
5. Heavy churn in the search and quantification workflows (ProSE, in-process Percolator, ProteomicsLFQ, OpenSwathWorkflow).

Most of the risk is in **behavioral regressions that unit tests cannot catch** — cross-platform path
handling, round-trip fidelity, and quantification accuracy.

**Already covered (verified in-repo):** every major new class has a `*_test.cpp`
(BrukerTimsFile, ThermoRawFile, Percolator, FragmentIndex, IsoelectricPoint, the PyProphet stats
classes, ProForma/USI/mzPAF/PEFF, WNet, …). TOPP round-trip tests exist for `.idparquet`
(IDMerger/IDRipper) and FileConverter parquet cells, plus `OpenSwathWorkflow_1..23`. CI matrix =
Windows 2025 / macOS-15 (ARM) / Ubuntu 24.04 (g++ & clang) / Ubuntu 24.04-ARM. Unit coverage is
good; the gaps are **integration, cross-platform, and scientific accuracy**.

---

## Tier 1 — Release-blocking, broad blast radius

### 1. Qt removal / core infrastructure modernization (BREAKING, foundational)

**Changed:** `QDate/QDateTime → std::chrono`; `QDir/QFile/QFileInfo → std::filesystem` (+ new
`PathUtils::to_path()`); `QProcess → boost::process` (ExternalProcess, PythonInfo, JavaInfo,
RWrapper); Qt string interop removed from String/DataValue; `QJson → nlohmann/json`;
`Qt6::Network → libcurl`; `minizip-ng → libzip`; `Qt6::Core` removed from PUBLIC link deps. Touches
**every file op, every timestamp, every subprocess launch, every HTTP call**.

**Covered:** compiles on the CI matrix. The recent `Fix non-ascii paths (#9452)` confirms this is a
live regression surface.

**Strategy:**
- **Unicode / edge-case path matrix** (the #1 risk): inputs *and* outputs with CJK/accented chars,
  spaces, very long paths, UNC paths on Windows — run FileConverter + an adapter on each of the 3
  OSes. Turn #9452 into a permanent regression test.
- **Subprocess launching** via boost::process for *every* adapter (Comet, Sage, MSGF+, Percolator,
  R/Java tools): executable resolution (PATH + sibling-binary, cf. #9204), arg quoting with
  spaces/unicode, stdout/stderr capture, exit codes, and the "exit code 14 when tool missing"
  contract.
- **libcurl UpdateCheck**: offline, proxy, and the `$HOME`-not-writable fix (#9075).
- **libzip**: ZIP64, transparent `.zip`/`.d.zip` decompression (#9139).
- **std::chrono**: locale/timezone-independent mzML run-start timestamps round-trip.
- **Downstream CMake consumer test**: a minimal external project linking libOpenMS must now link
  `Qt6::Core` itself (no longer PUBLIC) — verify it builds.

### 2. pyOpenMS nanobind rewrite (BREAKING, ~544 classes)

**Changed:** entire Autowrap/Cython backend replaced by hand-written nanobind; ~50 mutable-ref
output-param bugs fixed, ~100 missing default args restored, copy-constructor dispatch fixed for 12
classes; many BREAKING renames (`get_df→to_df`, `_mv→_view`, enum classes).

**Covered:** `src/pyOpenMS/tests/` pytest suite; docstring-format tests (#8674).

**Strategy:**
- Run the **full pytest suite on every OS × every supported Python (3.11–3.14)** via cibuildwheel;
  smoke-import + run one real workflow inside each built wheel.
- **3.5→3.6 API regression corpus**: run the published pyOpenMS tutorials/notebook examples against
  both and diff outputs — catches silent semantic drift in type casters, default args, and
  reference/ownership semantics.
- Targeted tests for each *fixed* class (AASequence copy, the 11 copy-ctor classes, the 50
  output-param fixes) so they don't regress again.
- **Zero-copy lifetime**: `get_peaks_struct`/`to_df`/pyarrow zero-copy must not use-after-free when
  the owning C++ object is GC'd.
- Verify deprecated aliases still work, and `.pyi` stubs match the runtime API. Wheels are built
  `WITH_THERMO_RAW=OFF` → test the *graceful absence* path (`hasattr` feature detection).

### 3. Native Parquet I/O + Arrow as a required dependency

**Changed:** new directory-based `.idparquet` / `.featureparquet` / `.consensusparquet` bundles
wired into ~20 TOPP tools; Arrow/Parquet now mandatory (`WITH_PARQUET` removed). Changelog already
shows several *found* fidelity bugs (dropped `fmap_metavalues`, run-id collisions, duplicate-id
handling).

**Covered:** IDMerger/IDRipper idparquet round-trip + FileConverter_32–36 cells.

**Strategy:**
- **Lossless round-trip fidelity matrix** for all three entity types: load→store→load and assert
  *semantic* equality including all metavalues, data-processing, run identifiers, search params,
  protein groups — not byte equality.
- Lock in the *already-fixed* bugs as regression tests: `setPrimaryMSRunPath` survives
  `.featureparquet`; IDRipper→IDMerger produces no run-id collision; duplicate-identifier rejection
  on store.
- **Cross-serialization** (parquet-in → XML-out and vice-versa) through each claiming tool, including
  the `.unknown` size-1 fallback contract in TOPPAS pipelines (directory-type formats are unusual
  here).
- **Schema-drift guard** (ArrowSchemaRegistry strict mode) + **quantms QPX interop** (read OpenMS QPX
  parquet in quantms). Build against **Arrow 23 and 24+** on all platforms (#9196 / #9221).

### 4. Quantification workflows (ProteomicsLFQ / OpenSwathWorkflow / IsobaricWorkflow)

**Changed:** ProteomicsLFQ — `tree_guided` re-enabled after the MapAlignerTreeGuided trafo fix
(which explicitly *shifts feature RTs slightly*), map-alignment reference-leak fix, BREAKING
`-out_qpx`, Bruker/Thermo input, Biosaur2 seeding. OpenSwathWorkflow — new **auto SWATH wave
scheduler** (default behavior change), Parquet outputs, PyProphet-ported inference (new
OpenSwathInfer/OpenSwathExport), OSW `RUN.ID` BLOB→INTEGER fix. IsobaricWorkflow — TMT 32/35-plex,
SPS-MS3 handling.

**Covered:** small deterministic `OpenSwathWorkflow_*` TOPP tests (calibration/estimation disabled).
**No end-to-end DDA-LFQ accuracy gate in CI.**

**Strategy → see the ProteoBench section.** Plus:
- Scheduler **determinism** test (auto wave scheduler at varying thread counts / memory thresholds
  must yield *identical* features — it changes scheduling, not results).
- Confirm the `RUN.ID INTEGER` fix yields non-empty FEATURE↔RUN joins downstream.
- Cross-validate the ported inference numerically against **PyProphet** on a shared `.osw`.

---

## Tier 2 — High risk: scientific correctness + can't run in public CI

### 5. Search & rescoring: ProSE + in-process PercolatorAdapter

**Changed:** PercolatorAdapter now defaults to **in-process** vendored Percolator 3.08 (no external
binary) for idXML/mzid PSM-FDR; subprocess still required for OSW/protein/peptide FDR. ProSE:
BREAKING tolerance schema, new outputs (idxml/qpx/parquet/pin), enzyme specificity, SNES, DB
chunking, two-pass calibration, protein FDR, ETD/EThcD ion fixes, isotope-range change.

**Strategy:**
- **In-process vs subprocess parity** (key gate for the new default): identical inputs → matching
  q-values/PEPs/PSM counts *and* post-filters (`-score:fdr`, `-best_per_spectrum_only`, per #9240).
- **ProSE correctness**: on a standard dataset, compare PSM/peptide/protein counts and FDR
  calibration vs Comet/Sage/MSGF+ at matched FDR; verify target/decoy separation and **no decoy
  leakage in merged output** (#9205). Migration tests for the BREAKING param renames (old names must
  error clearly).
- Dedicated **ETD/EThcD c/z dataset** to confirm the ion-forwarding + HyperScore fixes.
- **`.pin` validity**: feed `-out_pin` into the real percolator CLI (#9195).
- **ASan/valgrind** over repeated train/rescore for the in-process leak fixes (#9257).

### 6. Native vendor readers: Bruker `.d` (BrukerTimsFile) + Thermo `.raw` (ThermoRawFile)

**Changed:** BrukerTimsFile (DDA/DIA-PASEF, MS1 aggregation, IM + hill centroiding, **BREAKING
native-ID format** `frame=…scan=…precursor=…` rippling through
SpectrumLookup/USI/PepXML/Percolator/IDMapper); ThermoRawFile (.NET bridge, **license-gated**, OFF on
Linux/aarch64 and in wheels); CometAdapter mzParser segfault workaround; new IMAGING module +
BrukerTimsImagingFile.

**Core problem:** the Thermo bridge is license-gated (CMakeLists warns + links a GH issue) and the
Bruker SDK is **not bundled** → these **cannot run on public GH runners**. Class tests exist but with
tiny/mock data.

**Strategy:**
- **Dedicated internal/self-hosted CI lane** with the Bruker SDK + Thermo bridge installed and small
  real `.d`/`.raw` fixtures.
- **Validate against msconvert/pwiz as reference**: convert the same `.d`/`.raw` via msconvert→mzML
  and via OpenMS native readers; diff spectrum counts, m/z, intensities, IM arrays, precursor
  assignments. Verify **native-ID round-trips** (MS:1002818) through
  SpectrumNativeIDParser/USI/scan-number extraction/IDMapper.
- **CometAdapter native-ID rewrite**: regression on the real Bruker DDA case that previously
  segfaulted Comet's mzParser (mixed `frame=`/`scan=`).
- Exercise every integrated adapter (Comet/Sage/MSGF+/ProSE/ProteomicsLFQ/OpenSwath/FileConverter/
  PeakPickerIM) with `.d` and `.raw` input. **Runtime graceful-degradation** test when the SDK is
  absent (clear error, not crash). Imaging: pixel geometry + ion-image extraction on a small MALDI
  `.d`.
- These vendor files also feed the ProteoBench timsTOF/Astral modules → reuse end-to-end.

---

## Tier 3 — Medium risk

### 7. New scientific algorithms & parsers

- **IsoelectricPoint** — validate pI/charge vs EMBOSS/Sillero references.
- **PyProphet-ported stats** (KernelDensityEstimation/RankData/MultipleTesting) — numeric cross-check
  vs PyProphet/scipy.
- **HydrophobicityProfile**, **ModifiedSincSmoother**.
- **FeatureLinkerWNet** — compare grouping vs existing linkers on a known multi-map set.
- **ProForma v2 / USI / mzPAF / PEFF** — HUPO-PSI spec test vectors + round-trip.
- *Bug-fix regression tests*: TheoreticalSpectrumGenerator b/a neutral-loss m/z (#9078),
  FeatureFinderAlgorithmPicked mz_tolerance swap (#9247), LinearResampler intensity-spike (#9127).

Strategy: extend existing class tests with **authoritative external reference vectors**.

### 8. TOPPView / GUI

Crash fixes (#8997 large-mzML right-click, #9206/#9263 SNAP-zoom divide-by-zero, #9447/#9448 Apply
TOPP tool, #9452 non-ascii). Not covered by ctest. Strategy: **manual GUI smoke-test matrix** on all
3 OSes (large mzML, 1D/2D/3D + IM views, extreme SNAP zoom, Apply TOPP tool, TOPPAS), and verify
`WITH_GUI=OFF` builds all non-GUI tools **without Qt**.

### 9. Build / packaging / dependency permutations

contrib-as-submodule (**fresh-clone build test**), Boost≥1.81 / Eigen≥3.4 / nanobind / libzip /
libcurl / HiGHS / wnetalign / opentims-system-vs-fetch, `ENABLE_TDL` now OFF (CWL off by default).
Strategy: **build-permutation matrix** (`WITH_GUI`, `WITH_THERMO_RAW`, `WITH_OPENTIMS`,
`WITH_WNETALIGN`, `ENABLE_TDL`, each `LP_SOLVER`) per platform; installer/notarization/Gatekeeper
smoke tests on macOS; and **re-validate Windows Release performance** — #8961 reveals Windows Release
binaries were unoptimized (`/Od`) since 2021, so any historical Windows perf baseline is invalid and
must be re-measured.

---

## ProteoBench strategy (the end-to-end accuracy gate)

**Availability: yes, and OpenMS is already integrated.** ProteoBench is public, uses a
human/yeast/E.coli hybrid (DDA + DIA, ProteomeXchange **PXD028735**), and **ProteomicsLFQ results can
already be submitted** to the DDA ion-level module (via the quantms path). The benchmark scores each
precursor's log2 fold-change A vs B against expected **Human = 0, Yeast = +1, E. coli = −2**, plotting
*number of quantified precursors* vs *weighted mean absolute error*. This is exactly the regression
signal unit tests can't provide.

| ProteoBench module        | Data / instrument        | OpenMS 3.6 changes it exercises |
|---------------------------|--------------------------|---------------------------------|
| DDA quant – ion level     | Orbitrap (PXD028735)     | ProteomicsLFQ (tree_guided, alignment fixes, QPX), in-process Percolator default, FeatureFinder/Biosaur2, ProSE/Comet/Sage |
| DDA quant – Astral        | Thermo Astral `.raw`     | **ThermoRawFile** native reader end-to-end |
| DIA quant – ion level     | Q Exactive HF-X          | OpenSwathWorkflow wave scheduler + PyProphet inference/export |
| DIA quant – diaPASEF      | timsTOF SCP `.d`         | **BrukerTimsFile** direct `.d` input + OpenSwath DIA |

**How to use it as a release gate:**

1. **Regression, not just absolute score:** run the OpenMS workflows on the ProteoBench sets under
   **3.6 and 3.5**, and gate on *no regression* in (a) quantified-precursor depth and (b) weighted
   MAE from expected ratios. Any drop gets root-caused — not waved through.
2. **Validate the changed code paths specifically:** run ProteomicsLFQ with both `star` (default) and
   the re-enabled `tree_guided` alignment; confirm the latter now yields informative residual stats
   and doesn't degrade ratios. Run OpenSwath with the auto wave scheduler and confirm depth/accuracy
   match the legacy batch path.
3. **Two-birds reuse:** the Astral `.raw` and diaPASEF `.d` ProteoBench files are ideal real-world
   fixtures for the Tier-2 native-reader validation — one dataset validates both the reader and the
   quantification accuracy.
4. **Sanity-anchor against the community:** compare OpenMS points to the public DIA-NN/Spectronaut
   points already on ProteoBench for each module to catch gross errors.

---

## Suggested release-gate ordering

1. Green CI matrix + **fresh-clone (submodule) build** + build-permutation matrix.
2. **Unicode-path** + **subprocess/adapter** regression suite (§1) — gates everything downstream.
3. **pyOpenMS** full pytest × Python 3.11–3.14 wheels (§2).
4. **Parquet round-trip fidelity** matrix (§3).
5. **Percolator in-process↔subprocess parity** + ProSE correctness (§5).
6. **ProteoBench DDA + DIA** regression vs 3.5 (§4 + ProteoBench).
7. Internal-lane **Bruker/Thermo** validation vs msconvert (§6), reusing ProteoBench raw files.
8. GUI manual smoke + macOS notarization / Windows-perf checks.

**Explicit CI gaps to flag to the team:** Thermo `.raw` (license) and Bruker `.d` (SDK not bundled)
won't run on public runners — they need an internal lane; the GUI isn't in ctest; and there is
currently no automated DDA-LFQ accuracy gate (ProteoBench would fill it).

---

## Sources

- [ProteoBench documentation (home)](https://proteobench.readthedocs.io/en/stable/)
- [DDA quantification – precursor ions module](https://proteobench.readthedocs.io/en/stable/available-modules/2-DDA-Quantification-ion-level/)
- [DIA quantification – precursor ions module](https://proteobench.readthedocs.io/en/stable/available-modules/4-DIA-Quantification/)
- [DIA quantification – diaPASEF module](https://proteobench.readthedocs.io/en/v0.8.19/available-modules/5-quant-lfq-ion-dia-diapasef/)
- [ProteoBench preprint (bioRxiv, 2025)](https://www.biorxiv.org/content/10.64898/2025.12.09.692895v2.full)
- [A comprehensive LFQ benchmark dataset (Nature Sci Data, PXD028735)](https://www.nature.com/articles/s41597-022-01216-6)
