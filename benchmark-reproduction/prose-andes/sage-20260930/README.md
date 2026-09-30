# Sage supplement to OpenMS benchmark issue #10364

This supplement runs Sage on the same 20 hash-selected mzML inputs and the same
three target/decoy FASTA files as the pinned ProSE/ANDES benchmark. It measures
20 Sage searches and three Percolator seeds per search. ProSE/ANDES comparison
counts are historical references, not new searches in this supplement.

## Exact Sage version

Official release asset:
https://github.com/lazear/sage/releases/download/v0.14.7/sage-v0.14.7-x86_64-unknown-linux-gnu.tar.gz

- Tag: `v0.14.7` (the stable release selected on 2026-09-30).
- Source commit: `99407db6e3754b31a9b88b7316a0aee67293c93f`.
- Archive SHA256: `e3dc6b41015cb167574f6c82525b75e946c094f30bd700271b05c051c30cbe8a`.
- Binary SHA256: `d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625`.
- The official binary prints `sage 0.14.6`; the v0.14.7 source tag's Cargo package
  metadata also retains `0.14.6`. The asset and hashes above identify the binary.
- The official precompiled Linux x86_64 binary was used; no Sage code was modified.

## Matched search space and native engine settings

The complete requested and resolved configurations for each file are recorded
in `reference/jobs.json`. `run_sage.py` regenerates their paths for a new workspace.

Matched with the existing suite: full Trypsin/P (`cleave_at=KR`, `restrict=null`,
`c_terminal=true`, `semi_enzymatic=false`); two missed cleavages; length 7–40;
peptide mass 100–9000 Da; fixed CAM(C); variable oxidation(M), at most one;
fixed TMT6plex or TMTpro on K and peptide N-termini where applicable. Both engines'
FASTA target and decoy sequences are passed unchanged to Sage with
`generate_decoys=false`, `decoy_tag=DECOY_`.

Precursor tolerance: +/-20 ppm, except TMTpro +/-50 ppm. Known charges 2–5.
Sage `isotope_errors=[-1,2]` subtracts the isotope mass from the observed mass,
matching ANDES' sign and ProSE's `[-2,1]` opposite sign. Sage's output
`isotope_error` is a mass in Da, not an integer isotope index. No precursor
calibration is applied by this Sage version.

Fragments: +/-0.5 Da for Velos CID/Lumos CID-TMT; +/-20 ppm for every other group,
including timsTOF. The historical ANDES auto timsTOF fallback is not matched at
0.5 Da; its separate high-resolution reference is also reported.

Sage native settings are explicit: deisotope=true (the actual code default in
this release), min_peaks=15, max_peaks=150, min_matched_peaks=4,
max_fragment_charge=null (up to precursor charge minus one), fragment_min_mz=150,
fragment_max_mz=2000, ion_kinds=[b,y], min_ion_index=2, bucket_size=8192.
`min_ion_index=2` excludes b1/b2/y1/y2 in Sage's preliminary index; full scoring
regenerates all b/y ions. Do not equate this parameter's semantics with ProSE's.
Despite the `fragment_min_mz`/`fragment_max_mz` names, this release applies these
bounds to neutral fragment masses in the index, before charge expansion.
Sage's internal preliminary candidate cap is 50. Report the top 10 PSMs;
chimera=false, wide_window=false. Quantification is disabled.

RT/ion-mobility prediction is disabled (`predict_rt=false`). Four Rayon threads,
four OpenMP threads, one file per invocation. Telemetry disabled by the CLI flag.
These are matched benchmark settings with native Sage scoring/preprocessing,
not an unmodified out-of-the-box default search. Initial methionine clipping is
absent from this Sage release and the historical ProSE builds; it is present in
the pinned ANDES engine. Sage removes target-identical decoy sequences internally
and deduplicates peptide forms; native engine candidate handling is retained.

## Percolator and validation

Use the same Percolator 3.09.0 binary and `-Y -U --seed 1/42/137` commands as the
original suite. Count target PSMs with q <= 0.01, separately within each file and
seed. Average seeds per file, then files per instrument group. Do not pool files'
q-values or call the seeds biological replicates.

`run_sage.py` converts the native PIN with these explicit rules:

- Resolve each ScanNr against the frozen selected-spectrum provenance. Validate
  charge, decoy labels, isotope mass, missed cleavages and modification identities.
- Reconstruct peptide modifications using the pinned pyOpenMS modification names
  (including TMT/TMTpro N-termini), preserve them in the normalized peptide string,
  and use `AASequence.getMonoWeight()` for CalcMass.
- Set ExpMass to `(precursor_mz - 1.007276466771) * charge`, identical for every
  candidate of a spectrum. Verify the native single-precision masses agree within
  0.02 Da; save the maximum observed differences. Do not recompute Sage's native
  mass-error/score features.
- Remove `FileName`, raw/predicted RT and mobility fields, and `posterior_error`.
  The last field is an output of Sage's built-in LDA/FDR; it is not passed into
  Percolator as an extra learned feature. Keep the other native numeric features.
- Sort by scan, native rank, peptide, label and isotope mass before assigning
  stable SpecId values. Sage's LDA output order therefore does not control
  Percolator row order. Check every numeric feature is finite.
- Verify both target/decoy classes, no redundant hypotheses and at most ten
  candidates per spectrum. Check all PIN rows are read by Percolator, input PIN
  hashes remain unchanged, and accepted results contain at most one PSM per scan.

The common native result is a separate raw-HyperScore rank-one competition,
with decoy-favouring ties, whole score groups and (D+1)/T <= 0.01. Sage's own
LDA-based spectrum q-value results are saved only as a secondary diagnostic:
with top-10 reporting they can accept multiple peptide assignments per scan and
are not equivalent to the main Percolator comparison. Native HyperScore TDC can
return zero even when Percolator subsequently identifies many PSMs.

## Reproduce

Linux x86_64 with an Ubuntu 24.04-compatible runtime, Python 3.12, curl, unzip,
tar and dpkg-deb. Install `libboost-filesystem1.83.0` and `libgomp1` for the
pyOpenMS/Percolator runtime if absent. The measured environment has an 8-GiB
memory limit; allow at least 15 GB free disk. No Rust compiler or OpenMS source
build is needed for this supplement.
Unpack this supplement as `sage-benchmark` next to `prose-andes-reproduction`.
Start in their common parent directory:

```bash
curl -fL https://raw.githubusercontent.com/OpenMS/OpenMS/a5705ec8e1cb19d40ebfa7bb16f6b5651aa1dc65/benchmark-reproduction/prose-andes/prose-andes-reproduction.zip -o reproduction.zip
echo 'e24cbe748001df4e77bd03b98b2bc47d18463948a6ce05e6f1f422806f5075b2  reproduction.zip' | sha256sum -c -
unzip reproduction.zip
(cd prose-andes-reproduction && sha256sum -c MANIFEST.sha256)
(cd sage-benchmark && sha256sum -c MANIFEST.sha256)

python3 -m pip install --target prose-andes-reproduction/nightly --no-deps --index-url https://pypi.openms.de pyopenms==3.6.0.dev20260928
python3 -m pip install --target prose-andes-reproduction/deps numpy==2.5.3 pandas==3.0.6 lxml==6.1.3

curl -fL https://github.com/percolator/percolator/releases/download/rel-3-09/percolator-v3-09-linux-amd64.deb -o prose-andes-reproduction/percolator.deb
echo '3488743548d607d468f5b1bdbc06e7d99d03af4f0bf00264a0a086e32d662cf1  prose-andes-reproduction/percolator.deb' | sha256sum -c -
dpkg-deb -x prose-andes-reproduction/percolator.deb prose-andes-reproduction/percolator

mkdir -p sage-benchmark/bin
curl -fL https://github.com/lazear/sage/releases/download/v0.14.7/sage-v0.14.7-x86_64-unknown-linux-gnu.tar.gz -o sage-benchmark/sage-v0.14.7-linux.tar.gz
echo 'e3dc6b41015cb167574f6c82525b75e946c094f30bd700271b05c051c30cbe8a  sage-benchmark/sage-v0.14.7-linux.tar.gz' | sha256sum -c -
tar --no-same-owner -xzf sage-benchmark/sage-v0.14.7-linux.tar.gz -C sage-benchmark/bin --strip-components=1

# Optional: download the two large PRIDE sources using 8 validated HTTP ranges at a time.
# This is the transport used for this run; complete source SHA256 checks still apply.
python3 sage-benchmark/fetch_tims.py
python3 prose-andes-reproduction/reproduce.py inputs
python3 prose-andes-reproduction/reproduce.py check-inputs
python3 sage-benchmark/run_sage.py
python3 sage-benchmark/evaluate_sage.py --require-complete --write
```

Use `run_sage.py velos_125_R1` for a single measurement; input preparation also
accepts a dataset ID. Use `--output-root` for another output directory. Existing
job output is never overwritten. `--resume` skips only completed summaries whose
PIN checksum still matches. For an independent rerun use a fresh output root.
The runner checks executable, mzML and FASTA hashes before searching.

The supplement contains scripts, all parameters, expected per-seed counts,
normalized-PIN checksums, search/Percolator logs and command records. It does not
embed raw mzML, FASTA, search-result tables, PINs or executables. The commands
regenerate those files. Use the comparison helper to verify counts and hashes;
report any discrepancy instead of replacing expectations.

## Independent replay checks

Two files were independently searched and rescored again: `velos_125_R1` (CID)
and `eclipse_tmtpro_10855` (labelled HCD). These add two searches and six
Percolator runs beyond the primary 20/60 matrix. `replay_checks.json` records
canonical PIN, three-seed-count and native-count comparisons; the replay logs
and configurations are under `reference/replays`.

```bash
python3 sage-benchmark/run_sage.py velos_125_R1 eclipse_tmtpro_10855 --output-root sage-benchmark/replay
python3 sage-benchmark/evaluate_sage.py --output-root sage-benchmark/replay
```

## Interpretation

The parent issue's limitations apply: two vendors, six instrument models, seven
acquisition groups; sampled DDA MS2 only; no SCIEX/Waters, ETD/DIA, quantification,
full-run calibration or independently measured true FDR. Results depend on the
exact source, native preprocessing, features and pinned datasets. This is not a
runtime benchmark or evidence that one engine is universally superior. The new
ProSE methionine PR is not included in the historical ProSE reference results.

## Measured results

Mean accepted PSMs per file over Percolator seeds 1/42/137, q <= 0.01.

| Instrument / acquisition | n | ProSE historical default | ANDES auto | Sage |
| --- | --- | --- | --- | --- |
| Velos CID | 3 | 2876.6 | 3000.6 | 2919.3 |
| HF-X HCD | 3 | 4383.3 | 4703.9 | 4480.9 |
| Astral HCD | 3 | 2162.1 | 3073.8 | 2204.4 |
| Lumos HCD LFQ | 3 | 5479.4 | 5616.2 | 5604.3 |
| Lumos CID TMT | 3 | 2398.9 | 2468.9 | 2399.8 |
| Exploris 480 TMTpro | 3 | 2247.4 | 2502.6 | 2596.3 |
| timsTOF HT | 2 | 1038.0 | 964.0 | 1050.0 |

The original ProSE/ANDES sources and results remain pinned at `a5705ec8e1cb19d40ebfa7bb16f6b5651aa1dc65`. Sage results were measured separately with identical mzML/FASTA hashes. Full per-file/seed results are in the supplement.

Download [sage-reproduction.zip](sage-reproduction.zip). SHA256: `d6583638788543ce065c3d06d7cd86a9394f083be66b47f1ae330bb490562bf3`.
