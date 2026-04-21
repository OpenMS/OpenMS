# MSGFPlusAdapter: PIN-direct output support

**Date:** 2026-04-21
**Scope:** `src/topp/MSGFPlusAdapter.cpp`, `src/openms/source/FORMAT/PercolatorInfile.cpp`,
`src/openms/source/ANALYSIS/ID/PercolatorFeatureSetHelper.{h,cpp}`, related tests.

## Background

[bigbio/msgfplus PR #23](https://github.com/bigbio/msgfplus/pull/23) is a fork of MS-GF+
that drops mzIdentML output entirely. The new fork:

- Makes `.pin` (Percolator INput) the default output format and the only one carrying
  full feature information; `.tsv` remains as a secondary text format.
- Hard-rejects `-o <path>.mzid` at the CLI.
- Removes ~2,100 lines of mzid-generation Java (MZIdentMLGen, MzIDParser, MzIDToTsv,
  ScoringParamGen, AnnotatedSpectra) and the `-useFragmentIndex` flag.
- Adds new ion-series run-length features in pin output: `longest_b`, `longest_y`,
  `longest_y_pct`. Also exposes Java-side `enzN`/`enzC`/`enzInt` and
  `lnDeltaSpecEValue`/`matchedIonRatio` in pin.

The current `MSGFPlusAdapter` requests mzid output and parses it (or runs MzIDToTsv as a
legacy fallback). Against the new fork the adapter would fail at MS-GF+ invocation
because of the `.mzid` extension.

## Goal

Make `MSGFPlusAdapter` work transparently with both upstream MS-GF+ (mzid) and the
fork (pin). Single PR. No flag-day for users.

## Non-goals

- Removing the existing mzid code path. Kept for upstream MS-GF+ users; can be removed
  in a follow-up once the fork is mainstream.
- Changing the public idXML schema or score type for downstream consumers
  (`PercolatorAdapter`, `IDFilter`, `IDMapper`, `ConsensusID`, `MzTabExporter`, etc.).
- Building a general-purpose pin reader. `PercolatorInfile::load()` already exists and
  is the right reuse target.

## Architecture

```
main_()
  |- parse params
  |- detectMSGFVariant(jar, java)            [new]   -> {MZID, PIN}
  |- validateParamsForVariant(variant)       [new]   -> reject mzid_out on PIN, pin_out on MZID, legacy_conversion on PIN
  |- createLockedDBIndex                     (unchanged)
  |- build process_params
  |    |- if MZID: -o tmp.mzid               (unchanged)
  |    \- if PIN:  -o tmp.pin                [new]
  |- runExternalProcess_
  |- if MZID:
  |    |- existing FileHandler.loadIdentifications(mzid) path
  |    |- switchScores_ / addMissingRTs / addMissingFAIMS
  |    \- PercolatorFeatureSetHelper::addMSGFFeatures (existing 5-feature set)
  \- if PIN:
       |- loadFromPin_(tmp.pin)              [new]   -> protein_ids, peptide_ids
       |- addMissingRTs / addMissingFAIMS    (reused)
       |- reindex_                           (reused)
       \- PercolatorFeatureSetHelper::addMSGFFeatures(..., from_pin=true)  [extended]
```

### Variant probe

`detectMSGFVariant(jar, java)` runs `java -jar MSGFPlus.jar` (no args). Both upstream
and fork print a usage banner including the `-outputFormat` description. The probe
captures stdout+stderr and:

1. Locates the line containing `-outputFormat` (case-insensitive).
2. Scans that line plus the next ~5 lines for the keywords `pin` and `mzIdentML`.
3. Returns `PIN` if `pin` is found and `mzIdentML` is not; returns `MZID` if
   `mzIdentML` is found and `pin` is not.
4. If both or neither are found, aborts with `EXTERNAL_PROGRAM_ERROR` and includes
   the captured banner in the error message (do not guess).

Implementation note: confirm exact banner wording against both the latest upstream
release and the fork's HEAD before finalizing the keyword set. The fixture-based
unit test (see Testing) pins both banners against the same probe code.

Probe runs once, result cached in a member field. ~1s overhead.

### Output extension

The fork errors out on `.mzid`. We always pass the matching extension
(`tmp.mzid` or `tmp.pin`) based on detected variant. Old MS-GF+ writes mzid regardless
of extension; new fork is happy because we never give it `.mzid`.

## PIN parsing (reuse `PercolatorInfile::load`)

Mirror the SageAdapter approach. In `MSGFPlusAdapter::loadFromPin_()`:

```cpp
StringList filenames;
StringList extra_scores = {
  "RawScore", "DeNovoScore", "lnSpecEValue", "lnEValue",
  "lnDeltaSpecEValue", "matchedIonRatio",
  "longest_b", "longest_y", "longest_y_pct",
  "enzN", "enzC", "enzInt",
  "isotope_error"
};
PeptideIdentificationList peptide_ids = PercolatorInfile::load(
    pin_temp,
    /*higher_score_better=*/false,
    /*score_name=*/"lnSpecEValue",
    extra_scores,
    filenames,
    /*decoy_prefix=*/"",                    // pin's Label column already has target/decoy
    /*threshold=*/std::numeric_limits<double>::infinity(),  // unused outside SageAnnotation
    /*SageAnnotation=*/false);
```

`PercolatorInfile::load()` already handles header parsing, charge one-hot decoding,
peptide+brackets -> `AASequence`, target/decoy labels via `Label`, ExpMass -> m/z,
extra-score extraction, and per-spectrum grouping.

### Post-processing in MSGFPlusAdapter

Done in a small loop, not in `PercolatorInfile`:

1. **Score-type continuity:** for each `PeptideHit`, `setScore(std::exp(getScore()))`;
   for each `PeptideIdentification`, `setScoreType("SpecEValue")`. Preserves the
   public score contract used by every downstream OpenMS tool.
2. **Mirror MS-GF+ scores into PSI CV terms** (back-compat for code reading these
   meta values, including `addMSGFFeatures`):
   - `MS:1002049` <- `RawScore`
   - `MS:1002050` <- `DeNovoScore`
   - `MS:1002052` <- `clamp(exp(lnSpecEValue), DBL_MIN, +INF)` (clamp prevents
     underflow to 0 from very strong hits with `lnSpecEValue ~ -700`)
   - `MS:1002053` <- `clamp(exp(lnEValue), DBL_MIN, +INF)`
3. **Rename** `isotope_error` meta -> `Constants::UserParam::ISOTOPE_ERROR` (matches
   mzid path).
4. **Construct `ProteinIdentification` + `SearchParameters`:** copy the construction
   block from the existing `legacy_conversion` branch (db, charges, mods, tolerances,
   enzyme, term-specificity). Set `missed_cleavages = 1000` if
   `max_missed_cleavages == -1`.
5. **Fill missing RTs and FAIMS:** call `SpectrumMetaDataLookup::addMissingRTsToPeptideIDs`
   and `addMissingFAIMSToPeptideIDs` on the loaded mzML (same as mzid path).

### Required upstream change

`PercolatorInfile::load()` line 259 unconditionally accesses
`row[to_idx.at("retentiontime")]`. MS-GF+ pin has no such column; the load would
throw. Make `retentiontime` optional - default RT to 0.0 when the column is absent;
RT then gets populated by `addMissingRTsToPeptideIDs` from the mzML. SageAdapter is
unaffected because Sage emits the column.

## Parameter surface

### New / changed

| Param              | Change |
|--------------------|--------|
| `pin_out`          | **NEW** optional output file (`.pin`). Mirrors `mzid_out`. Only valid on PIN variant. |
| `mzid_out`         | Unchanged. Only valid on MZID variant. |
| `legacy_conversion`| Unchanged. Only valid on MZID variant. Hard-error if set on PIN. |
| `add_features`     | On PIN variant: still passed to MS-GF+ (new Java code respects `-addFeatures`); doc note added. |
| All others         | Unchanged. |

### Validation matrix (after probe, before search)

|                       | MZID variant         | PIN variant          |
|-----------------------|----------------------|----------------------|
| `mzid_out` set        | OK                   | `ILLEGAL_PARAMETERS` |
| `pin_out` set         | `ILLEGAL_PARAMETERS` | OK                   |
| `legacy_conversion`   | OK                   | `ILLEGAL_PARAMETERS` |
| `out` (idXML)         | OK (optional)        | OK (optional)        |
| no `out` and no \*_out| error                | error                |

Error messages name both the conflicting param and the detected variant so the user
understands which MS-GF+ jar is in use.

## Percolator feature-set declaration

Extend `PercolatorFeatureSetHelper::addMSGFFeatures` with a `bool from_pin` flag:

```cpp
void addMSGFFeatures(PeptideIdentificationList& peptide_ids,
                    StringList& feature_set,
                    bool from_pin = false);
```

- `from_pin == false` (default; existing callers): unchanged - declares
  `MS:1002049`, `MS:1002050`, `MS:1002052`, `MS:1002053`, `ISOTOPE_ERROR`.
- `from_pin == true` (new MSGFPlusAdapter pin path): additionally declares
  `lnDeltaSpecEValue`, `matchedIonRatio`, `longest_b`, `longest_y`, `longest_y_pct`.

`enzN`/`enzC`/`enzInt`/`dm`/`absdm`/`peplen` are intentionally **not** declared.
`PercolatorAdapter`'s `PercolatorInfile::preparePin_` recomputes them from sequence
unconditionally, so any pre-stored values would be overwritten - extracting them from
MS-GF+ pin would be wasted work.

`PercolatorAdapter` reads `extra_features` from `SearchParameters` and emits these
columns into its own pin -> the richer features land in Percolator without further
changes there.

## What becomes redundant on the OpenMS side

| Item                                                                 | Status with PIN path |
|----------------------------------------------------------------------|----------------------|
| `MzIDToTsv` legacy conversion path (~150 lines, `if(legacy_conversion)`) | Still used on MZID variant; not removed (back-compat). |
| `splitSequence_`, `modifyNTermAASpecificSequence_`, `modifySequence_`, `generateInputfileMapping_` | Used only by the legacy_conversion mzid path. Kept. |
| `switchScores_`                                                      | Not called on PIN path (score set directly during load). Still used on MZID path. |
| `PercolatorAdapter`'s recomputation of `dm`/`absdm`/`peplen`/`enzN`/`enzC`/`enzInt` | No change. Recomputed from sequence regardless of source. |
| `IsotopeError` -> `ISOTOPE_ERROR` rename (mzid postprocess)            | Still used on MZID path. PIN path renames inline during load. |
| `addMSGFFeatures` no-op zero-fill of missing scores                  | Kept - applies to both paths defensively. |

**Net redundancy after this change:** none of the existing OpenMS code is removed.
The new pin path is purely additive. The fork eliminates ~2,100 lines of mzid Java,
but that's on the MS-GF+ side. If MS-GF+ moves entirely to the fork in a future
release, a follow-up PR can delete the mzid path and its helpers.

## Testing

### Topp tests

`src/tests/topp/THIRDPARTY/` already contains `MSGFPlusAdapter_1.ini` (mzid path,
default) and `MSGFPlusAdapter_2.ini` (mzid path with legacy_conversion).
`third_party_tests.cmake` gates them on detection of `MSGFPlus.jar`.

Add `MSGFPlusAdapter_3.ini`:
- Same FASTA/mzML inputs as test 1.
- Gated on detection of the fork jar (probe finds `pin` in `-outputFormat`).
- `out` = idXML, no `mzid_out`.
- Reference idXML built once with the fork; checked into THIRDPARTY/.

Tests 1 and 2 keep their existing reference files unchanged - back-compat regression
net.

### Unit tests

- `PercolatorInfile_test.cpp`: add a case loading a small synthetic pin file lacking
  the `retentiontime` column - verifies optional-column behavior, RT defaults to 0,
  no exception. Self-contained; no MS-GF+ needed.
- `PercolatorFeatureSetHelper_test.cpp`: add a case calling
  `addMSGFFeatures(..., from_pin=true)` and asserting the extended feature list
  contains the new names.
- Variant probe: two checked-in fixture banner files in `src/tests/topp/THIRDPARTY/`
  (`msgfplus_banner_mzid.txt`, `msgfplus_banner_pin.txt`). A small unit test parses
  each through the same regex/string-match the probe uses. No Java needed.

### Documentation

- Update the doxygen block at the top of MSGFPlusAdapter.cpp:
  - Replace the "three-step mzid -> tsv -> idXML" description with a fork-aware
    version: results are read either from `.mzid` (legacy) or `.pin` (fork),
    selected automatically.
  - Note that `legacy_conversion` and `mzid_out` are MZID-only; `pin_out` is
    PIN-only.
- CHANGELOG entry: "MSGFPlusAdapter: support for pin-direct MS-GF+ output
  (bigbio/msgfplus PR #23)".

## Rollout / sequencing

Single PR. No flag-day for users:
- Existing upstream MS-GF+ users: zero behavior change (mzid path identical).
- Fork users: pin path activates automatically via the probe.

No deprecation needed yet. Once the fork is mainstream, follow-up PR removes the
mzid path and its helpers.

## Risks / edge cases

- **Probe parsing fragility:** the banner format could change in future fork
  releases. Mitigation: probe matches on `pin` keyword near `-outputFormat`, not
  exact string; if neither marker found, abort with clear error including captured
  output.
- **Future PercolatorInfile changes:** the optional-RT change touches code SageAdapter
  relies on. Mitigation: existing PercolatorInfile_test cases continue exercising
  present-RT path; new test exercises absent-RT path.
- **Score underflow:** `exp(lnSpecEValue)` can underflow to 0 for very strong hits
  (`lnSpecEValue ~ -700`). Mitigation: clamp at `std::numeric_limits<double>::min()`
  before storing as `MS:1002052` so downstream `log()` doesn't NaN.
- **JAR location for probe:** if `executable` path is invalid the probe fails the
  same way the existing `JavaInfo::canRun` check does. Probe runs *after* that check,
  so user gets the existing clear "Java is needed" error first.
