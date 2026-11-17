# MassTraceDetection drift fixture for issue #10905

This is a characterization of the **known unbounded extraction bug**, not a
chemically validated mass trace. It provides a small shared input for the
association/guard work discussed in
[OpenMS #10905](https://github.com/OpenMS/OpenMS/issues/10905).
A later algorithm fix should replace the wide-span expectations in
`MassTraceDetection_test.cpp` with scientifically validated expectations.

## Source and reduction

- Source: `doc/pyopenms/src/data/Metabolomics_2.mzML` in this repository at
  `d8286837b002763579ec6e036cab0b5c294f3ec1`.
- Source SHA256:
  `b3762268d923944f76863f6b285c5d5783fcf34e944669d9aa6ae7faf5243d2c`.
- Fixture: `MassTraceDetection_drift_10905.mzML.gz` (16,026 bytes).
- Fixture SHA256:
  `a0cd1a4cc8a5f27006ccb99ee9534163c09c693d7eae9e031ec310f96acd764f`.

Keep the contiguous MS1 slice with RT >=380 s and the inclusive m/z range
694.2-695.6. The source ends at RT 599.920501983 s. Preserve scan order, original
native IDs, RTs, m/z values and intensities, including the 21 scans that become
empty after mass filtering. The fixture has 181 MS1 scans and 778 peaks; the
full input has 329 MS1 scans. The crop retains competing peaks, rather than
only the 119 observations ultimately assigned to the target trace.

Reduction was checked in stages: m/z ranges 690-700, 693-697, 694-696 and
694.2-695.6, first with all MS1 scans and then with RT >=380 s. Each crop retained
the same 119 accepted RT/m/z positions. After preparing the final fixture, a
native diagnostic run confirmed that **every accepted RT, m/z and intensity**,
the extraction seed, and the sequence of centroid/SD values before accepting
each peak are identical between the full input and fixture. No scans within
the retained interval were removed. This is a bounded reduction, not a claim
of a globally minimal reproducer.

Unrelated acquisition metadata and auxiliary scan descriptions are omitted.
Only native ID, RT, MS level, spectrum type, instrument settings, and the two
peak arrays are copied into fresh spectra. The file uses the standard OpenMS
mzML writer and gzip compression; reading the test requires no Python packages.

The following recipe reproduces the decoded fixture using pyOpenMS 3.6.0 and
NumPy, starting in the repository root. Writer metadata and compressed bytes
may differ between library versions; compare decoded peaks when regenerating.

```python
import gzip
import tempfile
from pathlib import Path
import pyopenms as oms

source = oms.MSExperiment()
oms.MzMLFile().load("doc/pyopenms/src/data/Metabolomics_2.mzML", source)
source.sortSpectra(True)
spectra = []
for spectrum in source:
    if spectrum.getMSLevel() != 1 or spectrum.getRT() < 380.0:
        continue
    mzs, intensities = spectrum.get_peaks()
    keep = (mzs >= 694.2) & (mzs <= 695.6)
    selected = oms.MSSpectrum()
    selected.setNativeID(spectrum.getNativeID())
    selected.setRT(spectrum.getRT())
    selected.setMSLevel(1)
    selected.setType(spectrum.getType())
    selected.setInstrumentSettings(spectrum.getInstrumentSettings())
    selected.set_peaks((mzs[keep], intensities[keep]))
    spectra.append(selected)
reduced = oms.MSExperiment()
reduced.setSpectra(spectra)
assert len(spectra) == 181
assert sum(s.size() for s in spectra) == 778
assert sum(s.size() == 0 for s in spectra) == 21
target = Path("src/tests/class_tests/openms/data/"
              "MassTraceDetection_drift_10905.mzML.gz")
with tempfile.TemporaryDirectory() as temporary:
    plain = Path(temporary) / "fixture.mzML"
    oms.MzMLFile().store(str(plain), reduced)
    with target.open("wb") as output:
        with gzip.GzipFile(filename="", mode="wb", fileobj=output,
                           mtime=0) as compressed:
            compressed.write(plain.read_bytes())
```

## Effective extraction parameters

Call `MassTraceDetection::run(const PeakMap&, ...)` with no `max_traces` limit.
No experimental parameters from the closed, unmerged PR #10903 are used.

| Parameter | Value |
| --- | --- |
| `mass_error_ppm` | 20.0 |
| `noise_threshold_int` | 10.0 |
| `chrom_peak_snr` | 3.0 |
| `reestimate_mt_sd` | `true` |
| `quant_method` | `area` |
| `trace_termination_criterion` | `outlier` |
| `trace_termination_outliers` | 5 |
| `min_sample_rate` | 0.5 |
| `min_trace_length` | 5.0 s |
| `max_trace_length` | -1.0 (disabled) |
| `ion_mobility_tolerance` | 0.01 (default; this fixture has no IM data) |

## Accepted trajectory and ambiguity

| Measurement | Full input and fixture |
| --- | --- |
| Accepted observations | 119 |
| First/last accepted RT | 388.858178160 / 599.920501983 s |
| RT span | 211.062323823 s |
| Minimum/maximum m/z | 694.24554443359375 / 695.4852294921875 |
| Full m/z span | 1.23968505859375 Th |
| Intensity-weighted centroid m/z | 694.8168204114992 |
| Final trace m/z SD | 0.2535357013084262 Th |
| Integrated area | 681749.4925898777 |
| Extraction seed RT / m/z / intensity | 393.578101856 / 694.822998046875 / 7240.923828125 |
| Intensity maximum RT / m/z / intensity | 440.696651312 / 694.59539794921875 / 12166.330078125 |

The full-input label is T1930; the fixture label is T14. Neither label is a
stable identity. The unit test finds the trace by its intensity maximum and
checks representative original observations, including the seed. The extraction
seed is not the final trace's intensity maximum: failed extraction attempts can
leave stronger observations available to later seeds.

`MassTraceDetection_drift_10905_trajectory.tsv` lists all 119 accepted observations
in RT order. `extension_step` gives the actual processing order (0 is the seed).
The initial extension goes down in RT; subsequent accepted observations extend
upward. `centroid_before_mz` and `sd_before_Th` were recorded immediately before
`processPeak_()` updates them. For the seed row they contain the initial mass
and SD. These values were collected using temporary logging at seed initialization
and in `processPeak_()`; no diagnostic logging is added to production code.
These diagnostic values are documentation, not additional test goldens.

`other_retained_peaks_in_gate_mz:intensity` lists the other retained peaks above
the noise threshold within the moving +/-3 SD mass window, as semicolon-separated
`m/z:intensity` pairs. It ignores ownership by earlier traces, so these peaks are
not necessarily available alternatives. The seed row shows its initial window;
the seed itself is selected by intensity ranking. Peaks outside the cropped mass
interval are not listed.

For example, the seed's initial SD is 0.01389646 Th (a +/-0.04168938 Th window).
Before accepting m/z 695.485229492 at RT 577.907320015 s, the center is
694.829644429 and the SD has reached 0.223486062 Th. Before accepting m/z
694.245544434 at RT 598.670037503 s, the center is 694.820524654 with SD
0.248043824 Th. The growing scale admits large excursions on both sides of a
center that has moved little from the seed. This case therefore records more
than a shift of the final centroid.

The table marks numerical competition and acceptance, not chemical identities.
No individual join is established here as an isotope crossing rather than
unrelated interference. Such assignments require further annotation; the full
span alone cannot identify the ions or prove quantitative accuracy. This tutorial
excerpt also does not establish the separate near-whole-acquisition duration
example in #10905.

## Validation

The full-input and fixture measurements were reproduced with the checkout's
unmodified `MassTraceDetection.cpp` at the source revision above, compiled with
GCC 13.3 and linked to the official Linux x86_64 OpenMS SDK dated 2026.11.14
(library revision `339429f`). The updated `MassTraceDetection_test.cpp`, the
existing `MassTrace_test.cpp`, and the checkout's test framework were compiled
and both class tests passed in that native harness. This is not a claim that a
complete source-tree build or the full OpenMS suite was run.
The decompressed fixture also passes `MzMLFile::isValid()` schema validation.

In a configured source-tree build, run:

```bash
cmake --build <build-dir> --target MassTraceDetection_test MassTrace_test --parallel 2
ctest --test-dir <build-dir> -R '^(MassTraceDetection|MassTrace)_test$' \
  --output-on-failure --no-tests=error
```
