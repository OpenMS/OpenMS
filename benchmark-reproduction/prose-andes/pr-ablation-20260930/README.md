# Isolated ProSE PR benchmarks — 2026-09-30

This supplement to OpenMS issue [#10364](https://github.com/OpenMS/OpenMS/issues/10364)
measures PRs #10365, #10366 and #10368 independently, using the same 20 selected
mzML files, three target/decoy databases and Percolator protocol as the frozen
ProSE / ANDES / Sage comparison. Each file contains 8,000 selected DDA MS2 spectra.
The suite covers six instrument models from two vendors and seven acquisition
groups. Exact URLs, scan selections and input hashes are in `reference/inputs.json`
and the original reproduction bundle. The archive paths called `eclipse_tmtpro`
actually contain Exploris 480 data, as identified by their mzML headers.

## Source isolation

All PR patches are applied separately to the common `codex/native-cid` base
`c1b2f202742270273a045ef4b56e7aa346e9871f`. These are independent changes, not a
cumulative stack. The PR branches originated at
`0af4cf62dc9e0ab944102a0647df1ecd431e32bb`; searching their original heads directly
would also remove intervening base fixes and confound attribution. The original
PR delta was applied with `git apply --3way`; every patch applied cleanly.

| Arm | PR source head | Activation / purpose |
| --- | --- | --- |
| `base` | `c1b2f202742270273a045ef4b56e7aa346e9871f` | Frozen search parameters; calibration=false |
| `pr65_auto` | `160412c13122c1ae2205580eeec5b43018f7a006` | calibration:enabled=auto |
| `pr66_priors` | `5f66ab7ee2b7e3487081838964b9491e7696a4b8` | annotate:self_trained_ion_priors=true |
| `pr68_cutoff` | `9a3d57796d3ef685f66f9962e49e4c3c90e5ac58` | Initial positive-intensity cutoff fix |
| `pr68_current` | `d21f55b2575c523b4f8be13f407a62b276819496` | Updated cutoff including positive subnormal floats |
| `base_mass_cal` | Common base | calibration=true, scoring:method=mass_accuracy |
| `pr65_mass_cal` | #10365 head above | Same parameters; isolate the fitted kernel |

During execution the PRs received follow-ups. #10365 head
`79e08f2948f6000787d4eef992bb6b33572f5419` and #10366 head
`6a29032936bb6e4d49c4fa40e712f4bbefd1cee0` change only `CHANGELOG` relative to
their measured heads; their algorithm and test code is identical. #10368 changes
`float::min()` to `float::denorm_min()` and was built and benchmarked again on all
20 files. Both measured versions remain in the records. `source_provenance.json`
contains the original heads, merge bases, exact applied Git tree IDs and patch
SHA256 hashes. The original and effective patches, including follow-up diffs,
are bundled. Reproduction uses the effective patch on the common base.

## Protocol

There are 128 unscaled searches: common base plus four independently compiled PR
variants on 20 files, and two kernel-control arms on 14 high-resolution files.
There are also four searches of intensity-scaled controls (two files, base and
latest #10368). Every search gets three Percolator seeds: 1, 42 and 137. This is
132 searches and 396 Percolator invocations; seeds are algorithmic repeats,
not independent biological replicates. No combined-PR arm is inferred from the
individual results.

Search parameters are fully recorded for every job in `reference/jobs.json`:
full Trypsin/P, two missed cleavages, peptide length 7–40, mass 100–9000 Da,
fixed CAM(C), variable oxidation(M), at most one variable modification; fixed
TMT6plex or TMTpro on lysines and peptide N-termini for labelled inputs.
Known precursor charges 2–5; precursor window +/-20 ppm except TMTpro +/-50 ppm;
isotope errors -2..1 in ProSE's sign convention. Fragments use 0.5 Da for Velos
CID and Lumos CID-TMT and 20 ppm for the other groups. Fragment m/z 150–2000,
minimum index ion 2, minimum five matches, up to 50 preliminary candidates,
top ten reported hypotheses, peptide deduplication=true, processed-spectrum
retrieval, auto scoring/fragment-charge/deisotoping/window settings. Each arm
changes only the setting listed above. Initial methionine clipping PR #10358 is
not included.

`pr65_auto` tests the new default calibration behaviour. High-resolution auto
scoring still resolves to HyperScore, so this arm does **not** measure the new
mass-accuracy kernel. The paired `*_mass_cal` arms activate mass-accuracy scoring
and calibration on both sides; the calibration pass itself uses the configured
7-ppm zero-centred kernel. The evaluator verifies identical fitted precursor and
fragment tolerances within each pair. Any yield difference between these two
arms therefore measures the fitted main-search kernel, not turning calibration
on. Requested and resolved settings, including fitted width and shift, are saved.

The priors are opt-in. For #10366, the evaluator removes only its three added
`ion_prior_*` PIN columns and compares all remaining bytes with the common base.
It also compares native TSV hashes. This tests that gains originate in the new
Percolator features with unchanged candidates and native scores. Training PSM
counts and fallback status are saved per file.

All arms use the same normalized PIN exporter as the frozen benchmark. Spectrum
metadata comes from the original scan provenance; ExpMass is neutral observed
mass, CalcMass is reconstructed by the pinned pyOpenMS AASequence code, and dm
uses ProSE's isotope sign. All candidates of a spectrum share its observed mass.
Calibration adjusts search tolerances rather than rewriting observed precursor
m/z. RT features are excluded. Native features are retained. Percolator 3.09.0
uses `-Y -U --seed SEED`, independently per file. Count target PSMs with q <= 0.01;
average seeds within each file, then files within an instrument group. Native
counts are separate rank-one target/decoy competition, decoy-favouring ties,
whole score groups and `(D+1)/T <= 0.01`.

`evaluate.py` recounts the actual Percolator outputs, checks accepted scan
uniqueness, complete PIN consumption, unique candidate hypotheses, unchanged
PIN hashes, raw output hashes, native TDC and the PR-isolation checks. Existing
published reference counts and PIN hashes are checked on reproduction. Binary,
mzML and FASTA hashes are checked before execution. `historical_pin_identical`
explicitly checks the bridge to the earlier ProSE reference.

## Intensity-scale control for #10368

The unscaled measurements mostly have large raw intensities, so a zero yield
change cannot exercise the old absolute 0.05 cutoff. `intensity_control.py`
derives one Velos CID and one Exploris TMTpro input by multiplying each spectrum's
intensities by a power of two so its maximum lies in [1,2). It edits only mzML
intensity binary arrays, preserving their stored float precision. It verifies
exact float32 reversibility and unchanged m/z arrays, retention times, precursor
metadata, scan IDs, spectrum counts and peak counts after parsing the result.
The transformation, every exponent and source/derived file hashes are saved.

The base and latest #10368 use the identical derived file and all other original
settings. The expected property is agreement of #10368's unscaled and scaled
normalized PIN, not merely a higher PSM count. The old base can lose low-intensity
peaks before normalization. Derived inputs are reproducible controls, not two
additional independent measurements. An initial pyOpenMS load/store attempt
rounded RT metadata and was rejected before any search; only the validated
binary-array transformation was searched.

## Build and environment

Linux x86_64, Ubuntu 24.04-compatible runtime, Python 3.12.14, GCC 13.3.0, C++23,
CMake 4.4.3, OpenMP. Algorithms are compiled from each isolated tree and linked
against the same pyOpenMS 3.6.0.dev20260928 OpenMS support library. This is the
targeted native harness used by the original comparison; it is not a full TOPP
or pyOpenMS binding build. Harness sources and generated headers are included.
Searches and Percolator use four OpenMP threads; CTest uses two. At most two
search pipelines run concurrently under an 8-GiB memory limit. Recorded times
are diagnostic and are not a controlled speed benchmark.

`HyperScore_test`, `FragmentIndex_test` and `ProSEAlgorithm_test` run for all five
builds; #10366 also runs `FragmentIonLikelihoodModel_test`: 16 CTest executions.
All must pass before their binary can be searched. Fixture data is pinned to
OpenMS commit `e25ac614f495ce1e4933f46b41da18199eeba3ef`. Build logs and CTest logs
are included, with executable hashes and applied source trees.

## Reproduce

Allow at least 8 GiB RAM and 25 GB free disk. Install `g++-13`, `git`, `curl`,
`unzip`, `libboost-filesystem1.83.0` and `libgomp1`. Unpack the supplement as
`pr-ablation`, beside `prose-andes-reproduction`. From their common parent:

```bash
curl -fL https://raw.githubusercontent.com/OpenMS/OpenMS/a5705ec8e1cb19d40ebfa7bb16f6b5651aa1dc65/benchmark-reproduction/prose-andes/prose-andes-reproduction.zip -o reproduction.zip
echo 'e24cbe748001df4e77bd03b98b2bc47d18463948a6ce05e6f1f422806f5075b2  reproduction.zip' | sha256sum -c -
unzip reproduction.zip
(cd prose-andes-reproduction && sha256sum -c MANIFEST.sha256)
(cd pr-ablation && sha256sum -c MANIFEST.sha256)

python3 -m pip install --target prose-andes-reproduction/nightly --no-deps --index-url https://pypi.openms.de pyopenms==3.6.0.dev20260928
python3 -m pip install --target prose-andes-reproduction/deps numpy==2.5.3 pandas==3.0.6 lxml==6.1.3 cmake==4.4.3
python3 -m pip install --target prose-andes-reproduction/header_deps cmeel-boost==1.84.0 cmeel-eigen==3.4.0.2
curl -fL https://github.com/percolator/percolator/releases/download/rel-3-09/percolator-v3-09-linux-amd64.deb -o prose-andes-reproduction/percolator.deb
echo '3488743548d607d468f5b1bdbc06e7d99d03af4f0bf00264a0a086e32d662cf1  prose-andes-reproduction/percolator.deb' | sha256sum -c -
dpkg-deb -x prose-andes-reproduction/percolator.deb prose-andes-reproduction/percolator

python3 prose-andes-reproduction/reproduce.py inputs
python3 prose-andes-reproduction/reproduce.py check-inputs
python3 pr-ablation/prepare_sources.py
python3 pr-ablation/build.py
python3 pr-ablation/run.py
python3 pr-ablation/evaluate.py --require-complete --write
python3 pr-ablation/intensity_control.py prepare
python3 pr-ablation/intensity_control.py run
python3 pr-ablation/intensity_control.py evaluate
python3 pr-ablation/diagnostics/calibration_losses.py
```

`prepare_sources.py` creates a dedicated OpenMS fixture clone when absent; it
refuses to modify an existing checkout with a different HEAD or local changes.
Source worktrees and effective patch tree IDs are verified. Use a fresh directory
if an existing fixture clone has another purpose. `build.py` accepts build names
and `--jobs 1` to reduce peak memory. `run.py velos_125_R1 --arms base pr66_priors`
runs a selected pair. `--output-root` preserves other results, and `--resume`
skips completed output only after checking its PIN. Partial output is never
silently overwritten.

The ZIP contains scripts, patches, parameters, expected counts and hashes,
resolved metadata, build/test logs, search and Percolator command/log records,
and all per-file/per-group tables. Large mzML/FASTA, PINs, binaries and result
tables are regenerated by the commands above. The checksum manifest covers all
bundled files. `reference/jobs.json` retains the original measured binary and
result hashes even when a local rebuild rewrites `build_records.json`.

A conversation/session transition interrupted two in-progress searches after
83 complete jobs. Their partial logs/parameters were archived and excluded;
the jobs were rerun, while completed jobs were resumed only after PIN checks.
`interruption.json` records the excluded attempts. After measured peak RSS was
available, the two-worker scheduler's memory reservation was reduced to allow
two HYE searches concurrently while keeping at least 1.2 GiB headroom.

One completed job (`eclipse_tmtpro_10863/base_mass_cal`) retained a truncated
search log despite complete, hash-validated outputs. It was independently
searched and rescored again: normalized PIN, all three seed counts, native
outputs, parameters and resolved settings agree exactly. The original log is
preserved; `log_recoveries.json` points to the complete replay log used for its
window-bound check. This adds one search and three rescoring runs beyond the
132/396 matrix, separately from the three historical TMTpro baseline replays.
The runner now publishes logs with an atomic rename after the process exits.

To independently repeat that replay without overwriting its bundled metadata:

```bash
python3 pr-ablation/run.py eclipse_tmtpro_10863 --arms base_mass_cal --output-root pr-ablation/replay-check
```

## Sage gap diagnostic and limitations

`diagnostics/analyze_candidates.py` compares the three TMTpro files with the
previous Sage supplement. It requires Sage outputs reproduced from the pinned
Sage bundle next to this directory, then runs with:

```bash
python3 pr-ablation/diagnostics/analyze_candidates.py
```

Comparison is by scan and modified peptide, retaining modification sites and
collapsing indistinguishable I/L. A Sage assignment absent from the accepted
ProSE set can still be present in ProSE's top ten candidates. These additional
accepted assignments include conflicts; they are not independently validated
true positives. The original three TMTpro ProSE runs were replayed exactly
(PIN hashes and all seeds) before the new common-base matrix.

`diagnostics/calibration_losses.py` separately compares base and auto-calibrated
accepted assignments, recording whether lost assignments remain in the calibrated
top ten. It tests precursor exclusion against the rounded window in the log,
considering every isotope offset and requiring a 0.05-ppm margin beyond the bounds.
This is a conservative diagnostic, not an exact reimplementation of FragmentIndex's
single-precision window arithmetic. Isotope-offset histograms are included.

The scope remains sampled DDA, two vendors, six models, no ETD/DIA or independent
true-FDR validation. Per-file calibration uses the selected 8,000 spectra rather
than the full original run. The prior model trains on confident PSMs from the
same file before Percolator cross-validation; this benchmark is not an out-of-fold
assessment of prior training. Nominal q-value gains alone do not establish true
FDR control. The comparison with ANDES and Sage uses their published, pinned
results; those engines are not rerun in this PR supplement.

## Measured results

Mean accepted PSMs per file over seeds 1/42/137 at q <= 0.01. PRs are independent.

| Instrument / acquisition | n | Common base | #10365 auto | #10366 priors | #10368 latest |
| --- | --- | --- | --- | --- | --- |
| Velos CID | 3 | 2876.6 | 2876.6 | 2890.1 | 2876.6 |
| HF-X HCD | 3 | 4383.3 | 4231.7 | 4386.6 | 4383.3 |
| Astral HCD | 3 | 2162.1 | 2262.3 | 2276.3 | 2162.1 |
| Lumos HCD LFQ | 3 | 5479.4 | 5261.7 | 5481.3 | 5479.4 |
| Lumos CID TMT | 3 | 2398.9 | 2398.9 | 2415.6 | 2398.9 |
| Exploris 480 TMTpro | 3 | 2247.4 | 2285.7 | 2285.9 | 2247.4 |
| timsTOF HT | 2 | 1038.0 | 1011.8 | 1030.7 | 1038.0 |

Kernel isolation: both sides use calibration=true and mass-accuracy scoring.

| Group | Fixed kernel PSMs | Fitted kernel PSMs | Delta | Native TDC delta |
| --- | --- | --- | --- | --- |
| HF-X HCD | 4206.7 | 4217.0 | +0.25% | -566.3 |
| Astral HCD | 2296.6 | 2292.2 | -0.19% | -571.7 |
| Lumos HCD LFQ | 5294.6 | 5299.1 | +0.09% | -160.3 |
| Exploris 480 TMTpro | 2267.3 | 2270.4 | +0.14% | -137.0 |
| timsTOF HT | 1003.8 | 1000.3 | -0.35% | +41.0 |

Intensity-scaled controls (latest #10368):

| File | Arm | Unscaled seeds | Scaled seeds | Identical normalized PIN |
| --- | --- | --- | --- | --- |
| velos_125_R1 | base | [2927, 2848, 2946] | [2410, 2340, 2375] | False |
| velos_125_R1 | pr68_current | [2927, 2848, 2946] | [2927, 2848, 2946] | True |
| eclipse_tmtpro_10855 | base | [2342, 2387, 2363] | [2350, 2336, 2328] | False |
| eclipse_tmtpro_10855 | pr68_current | [2342, 2387, 2363] | [2342, 2387, 2363] | True |

TMTpro diagnostic: 3268/4393 (74.4%) of Sage's additional accepted assignments are already in ProSE's top ten, summed over the three files and three seeds. Assignments are not independent observations across seeds and are not validated true positives.

Download [pr-ablation-reproduction.zip](pr-ablation-reproduction.zip), SHA256 `ec8e4813196141bd6394690fe3a3bac96841f9c0c0e29b577f58a43fa919ab09`.
