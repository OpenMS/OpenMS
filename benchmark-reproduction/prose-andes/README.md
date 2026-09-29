# Reproduce the benchmarks in OpenMS PR #10335

This bundle makes the experiments reproducible without access to a chat or private
file store. It contains **182 saved search configurations and 546 per-seed result
counts**, original source patches, the C++ harness, full sampling metadata,
public source URLs, checksums, PIN normalization and evaluation code. It does not
contain the large source mzML files, native result tables or engine executables;
the commands below regenerate them. Expected PIN SHA256 and per-seed counts allow
checking the regenerated results. Retained/rejected experiments are both recorded;
the 182 includes 86 broad, 49 retrieval, 21 ranking and 26 deduplication experiments.

## 1. Environment and dependencies

Recorded platform: Linux x86_64, Ubuntu 24.04-compatible runtime, Python 3.12.14,
GCC 13.3, C++23, OpenMP, CMake 4.4.3, Rust 1.87.0. Search processes used four OpenMP
threads and at most two searches at once within an 8-GiB memory limit. The public
launcher is deliberately serial. Allow at least 15 GB of free disk for a small
subset; substantially more for retaining outputs from all 182 runs. Source data
total 13.16 GB; the largest download is 2.73 GB. The downloader processes files
one at a time and removes verified originals after sampling unless --keep-raw is
set. Time measurements from the original concurrent runs are not speed claims.

From the extracted bundle root, install system prerequisites (`g++-13`, `make`,
`git`, `curl`, `pkg-config`, `libssl-dev`, and the Ubuntu 24.04 runtime package
`libboost-filesystem1.83.0`). Install Rust using the official rustup instructions
if needed; ANDES' pinned rust-toolchain.toml selects 1.87.0. Then:

```bash
python3 -m pip install --target nightly --no-deps --index-url https://pypi.openms.de pyopenms==3.6.0.dev20260928
python3 -m pip install --target deps cmake==4.4.3 numpy==2.5.3 pandas==3.0.6 lxml==6.1.3
python3 -m pip install --target header_deps --no-deps cmeel-boost==1.84.0 cmeel-eigen==3.4.0.2
export PYTHONPATH="$PWD/nightly:$PWD/deps"
export OPENMS_DATA_PATH="$PWD/nightly/pyopenms/share/OpenMS"
export OMP_NUM_THREADS=4

curl -fL https://github.com/percolator/percolator/releases/download/rel-3-09/percolator-v3-09-linux-amd64.deb -o percolator.deb
echo '3488743548d607d468f5b1bdbc06e7d99d03af4f0bf00264a0a086e32d662cf1  percolator.deb' | sha256sum -c -
dpkg-deb -x percolator.deb percolator
echo '1f067b5d438a3a88be8a88f636844baea824e239fd2c5c053462ae56fd0e7c15  percolator/usr/bin/percolator' | sha256sum -c -
```

`native_dedup/provenance/python_packages.json` records the complete original Python
package inventory, including plotting packages not needed for rerunning searches.
The harness compiles ProSEAlgorithm, FragmentIndex, HyperScore and
TheoreticalSpectrumGenerator from source and links the unchanged support library
from the pinned official pyOpenMS wheel. This is **not a full clean OpenMS build**,
and Python bindings were not rebuilt for this deduplication follow-up.

## 2. Recreate the exact inputs

```bash
# Start with one CID measurement, or omit the dataset name for all 20.
python3 reproduce.py inputs velos_125_R1
python3 reproduce.py check-inputs velos_125_R1
```

The manifest is `broad_benchmark/manifest.json`. `inputs.json` freezes original
download sizes/SHA256 and selected mzML hashes. Full original IDs, acquisition
metadata, precursor m/z/charge and selected scan identities are preserved under
`broad_benchmark/provenance/<dataset>.json`.

Selection: take the 8,000 smallest SHA256(dataset ID + NUL + original native ID)
values among nonempty MS2 spectra with selected-ion m/z, known charge 2–5 and no
profile-spectrum flag; sort the selected spectra by RT. Preserve peak arrays and
instrument references. Normalize non-scan native IDs to scan=N for both engines,
keeping originals in provenance. For timsTOF, represent exporter userParam CID
also as CV MS:1000133; analyzer metadata remains TOF. See the exact `prepare.py`.
Any downloaded/prepared hash mismatch aborts; do not update expected checksums to
silence it. Existing mismatched inputs must be investigated before searching.

FASTA source:
https://archive.openms.de/openms/benchmarks/pride-benchmarks/lfq/QExactiveHF/ProteoBench_Module_2/ProteoBenchFASTA_MixedSpecies_HYE.fasta

The complete FASTA is **16,670,994 bytes / 31,889 targets**, SHA256
`d9ac434d88492c10c8e9a587ee7dbc9480fa0995fa07a6ba35a7da8abf39aa25`, Git blob
`d832aa92d3f3869ab17ba24a0a8b190b38e64858` in bigbio/quantms-test-datasets.
An earlier truncated copy invalidated earlier HYE comparisons; those are excluded.
Human/yeast selects accessions ending _HUMAN/_YEAST or containing Cont_; human-only
selects _HUMAN or Cont_. Decoys are one full protein-sequence reversal with DECOY_
prefix. Both engines receive the same concatenated target-decoy FASTA, and ANDES
is explicitly forbidden to generate another set. `prepare_databases.py` is exact.

| Database | Targets | Decoys | Combined FASTA SHA256 |
|---|---:|---:|---|
| Human/yeast | 27488 | 27488 | `35614831eddb913fab678c63461d6b4cb44ca73f5aebe2ead1a00020d791def6` |
| Human/yeast/E. coli | 31889 | 31889 | `efbf5ac224bf51fca252f58cd55df5e2aad7eddd824d0dee9959b086e147bb59` |
| Human | 20767 | 20767 | `6f86369103a9003560661f559cbd5bd1c37e55187670f662d69f91f81359a3ab` |

## 3. Build the exact source variants

`versions.json` defines every build key and patch. `build` uses separate detached
worktrees and preserves existing outputs. Do not edit these source worktrees.

| Build key | Base commit | Additional patch |
|---|---|---|
| baseline | `bcbc2f13052bbfb90b9441bed1e298a67b0981c7` | none |
| window | same baseline | `native_gap/provenance/window.patch` |
| query | same baseline | `native_gap/provenance/query.patch` |
| charge_bound (rejected) | same baseline | `native_gap/provenance/charge_bound.patch` |
| mass | `46e3f6573fc70877e863bb3fa8a87af89a13ac50` | `native_rank/provenance/mass.patch` |
| local (rejected native score) | same retrieval commit | `native_rank/provenance/rejected_local_with_tests.patch` |
| dedup | `48dc1c22adc25ca2656448676f0ba2a663ccb6fe` | `native_dedup/provenance/first.patch` |
| final | `0af4cf62dc9e0ab944102a0647df1ecd431e32bb` | none |
| andes | `b7eaece219ef45198eb233030cf9d69075988f72` | none, shipped model store |

Published retrieval/mass/dedup trees are recorded in the corresponding provenance
folders. Experiments used frozen prototypes; their patches above preserve the
source actually measured. The final compact dedup implementation has tested tree
`3ffb32ea9ac43d3e0a833863d6b6c38e602fcb72`. Do not replace all historical binaries
with final defaults: deduplication changed the default. For final-code validation,
use --binary bin/final explicitly on native_dedup jobs; the original follow-up
verified eight final native/PIN comparisons, including legacy controls.

```bash
python3 reproduce.py build baseline
python3 reproduce.py build dedup
python3 reproduce.py build andes

# To test the published implementation itself:
python3 reproduce.py build final
deps/cmake/data/bin/cmake --build builds/final -j2
deps/cmake/data/bin/ctest --test-dir builds/final --output-on-failure
```

Only final source is expected to pass the final regression-test expectations.
Some frozen exploratory patches contain early test expectations later corrected.
The helper builds only the search target for historical variants.

## 4. Run the saved experiments

```bash
python3 reproduce.py list
python3 reproduce.py list --stage native_dedup
python3 reproduce.py run --job broad_benchmark:velos_125_R1:prose_auto
python3 reproduce.py run --job broad_benchmark:velos_125_R1:andes
python3 reproduce.py run --job native_dedup:velos_125_R1:default

# Final-source byte/result check, in a separate output directory:
python3 reproduce.py run --job native_dedup:velos_125_R1:default --binary bin/final --output-root final-checks

# All 26 dedup experiments (after preparing all inputs and building dedup):
python3 reproduce.py inputs
python3 reproduce.py run --stage native_dedup
```

Do not repeat a job into a nonempty output directory; use a fresh --output-root.
The last command assumes those jobs have not already been run there.
To reproduce **every** recorded experiment from scratch, build baseline, window,
query, charge_bound, mass, local, dedup and andes, prepare all inputs, then run:

```bash
python3 reproduce.py run --output-root full-rerun
python3 evaluate.py --output-root full-rerun
```

Every job in `jobs.json` includes its source key, full parameter dictionary,
original command, expected normalized PIN SHA256, candidate row count, resolved
settings, input hashes, native count and seed counts. Paths in the original logs
describe the original machine; the launcher rewrites input/output paths and runs
the pinned binary. Outputs are retained under
`reruns/<stage>/<dataset>/<arm>/`: native.tsv, search.idXML (ProSE), native.pin,
input.pin, target/decoy tables for s1/s42/s137, logs, parameters and comparison.json.
Any PIN/count discrepancy is recorded and causes failure; floating-point/platform
differences are possible and must be reported rather than silently accepted.

## 5. Shared search and rescoring protocol

Saved per-job parameters are authoritative. The common space is full Trypsin/P
(ANDES Trypsin), two missed cleavages, length 7–40, known charges 2–5, fixed CAM(C),
at most one variable Oxidation(M), and fixed TMT6plex or TMTpro on K and peptide
N-terminus where indicated. Precursors ±20 ppm, except TMTpro ±50 ppm. ProSE isotope
errors −2..1 correspond to ANDES −1..2 because their signs differ. Calibration is
off. Fragments: 0.5 Da for Velos CID/Lumos CID-TMT, 20 ppm for high resolution.
ProSE candidate cap 50 and report top 10, except explicitly recorded diagnostic
variants. Engine-specific peak processing and scoring features remain native.

ANDES auto uses its pinned models, --candidate-index mmap,
--mmap-window-cache-candidates 200000, --fragment-index-slice-da 10,
--decoy-strategy none, --decoy-prefix DECOY, --threads 4, --precursor-cal off,
--protocol standard for unlabelled/TMT for labelled data, --top-n 10. Full commands
are saved in jobs.json. timsTOF auto selects the reported LowRes fallback; the
separate andes_highres arm forces hcd_qexactive_tryp and is diagnostic only.

PIN normalization (`native_dedup/export.py`): assign every candidate of a scan
the same observed neutral ExpMass = (precursor_mz − 1.007276466771) × charge;
CalcMass is AASequence.getMonoWeight(). Use isotope difference 1.0033548378 with
the engine-specific sign for dm/absdm. Normalize sequence/modification notation,
scan IDs, lengths and internal cleavage counts. Remove RT-prediction columns and
verify every numeric feature is finite. No external predicted RT/intensity
features are added. The exact feature columns come from each engine's native PIN.

```bash
for seed in 1 42 137; do
  mkdir -p "s$seed"
  percolator/usr/bin/percolator -Y -U --seed "$seed" \
    --results-psms "s$seed/target.tsv" \
    --decoy-results-psms "s$seed/decoy.tsv" input.pin
done
```

The launcher runs this automatically and checks that Percolator read every PIN
row. Accepted PSMs have q-value ≤0.01; average three seeds within a file, then
average files within an instrument group. Do not pool separate-file FDR estimates.
Native TDC uses one top candidate per scan, favors decoys at exact score ties,
adds complete equal-score groups, and uses (D+1)/T≤0.01. Use score for ProSE;
ANDES strong mode uses RawScore, otherwise RankScore+EdgeScore, as in analyze.py.

```bash
python3 evaluate.py --expected   # published reference group means
python3 evaluate.py             # corresponding newly generated means
python3 gap.py --output-root full-rerun --stage native_dedup --arm mass7_evidence
```

`gap.py` recreates the five Astral sequence-assignment categories at seed 42;
both the requested ProSE arm and the three ANDES controls must have been run.

Native/PCL counts alone do not establish actual FDR. ANDES agreement is a comparison,
not ground truth; sequence agreement ignores modifications and equates I/L.
Seeds are not biological replicates. Astral A2 was exploratory. ANDES has known
benchmark-series overlap, so no blind-test or training-set-independence claim is
made. The archive's Eclipse folder contains Exploris 480 data. No performance
validation is claimed for SCIEX/Waters, ETD/DIA, full raw conversion, quantification,
or broad PTM localization. All these limits remain in the PR description.

`MANIFEST.sha256` hashes the bundle files. Source mzML files and FASTAs have their
own hashes in inputs.json/databases.json; newly generated search results are
checked against jobs.json. Public source URLs and checksums, not private artifact
links or unspecified engine HEAD revisions, define the reproduction contract.
