# Local result: memory round 5 and runtime rounds 1 and 2, 25 November 2026

**Verdict: all three arms byte-identical to the `3b1d90ed67` baseline on raw01 seed0 (3/3 files each). 19/19 OpenDIA CTest gates pass.** This is a local laptop diagnostic with the 603k-precursor library. The full-library measurement is frozen separately in [PLAN_RT2](../astral-full-rt2/PLAN_RT2.md).

## Source

| Arm | OpenMS commit | Change |
|---|---|---|
| mem5 | `a11a95ec05` | Extraction collectors and stage-2 blocks released before scoring; features moved into the scorer and back; sibling columns widened in place |
| rt1 | (uncommitted intermediate, diff in evidence) | Scorer forward pass with transposed weights: loops interchanged so independent hidden and embedding units vectorise; every unit sum keeps its input order |
| rt2 | `bdb0159ff4` | rt1, plus: decoy covariance in `readout` split across threads by entry (each entry still sums rows in ascending order); calibration scans a compact sorted (m/z, window, closed) slot array before touching panel members (same skip and close rules); by-ion counting uses the bucketed `SpectrumMzIndex` (returns exactly `MZBegin`) |

New unit tests: transposed forward pass equals the plain one bit for bit (three network shapes including hidden > 256), and `readout` with 1 and 4 threads gives identical `w` and objective above the parallel threshold.

## Profile that drove the changes (mem5 binary, `perf record`, 16 threads, `-no_rt_refine`)

The profiled run took 8:20 wall at 451 % CPU. Self time by stage window:

- Scoring (two passes, about 190 s of 500 s): `WhitenedNetwork::objective` forward and backward 69 %, `embed` 12 %, `expm1`/`tanh` 12 %. The forward loop was scalar (`cvtss2sd`, `mulsd`, `addsd`): the inner loop ran over inputs, a reduction that `-ffp-contract=off` without reassociation cannot vectorise.
- Calibration (144 s): `CalibrationCollector::consumeSpectrum` 74 %, hot instructions the loads of `member.window_low` and `member.rt.back()` across about 600k 300-byte members per isolation window. Most members are outside their window and are skipped, but only after a cache miss.
- Extraction: spread over the trace and grid collectors, `MSSpectrum::MZBegin` (5 %, from `countByIons_`, the one caller not using the m/z index) and malloc/free (about 12 %).

## Local measurements (raw01 seed0, 16 threads, same inputs, desktop not idle)

| Arm | Calibration s | Search s | Extraction s | Features s | Scoring pass 1 / 2 s | Wall | CPU | Peak RSS GB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| mem5 | 116 | 55 | 177 | 22 | 197 / 223 | 13:20 | 475 % | 16.6 |
| rt1 | 141 | 59 | 184 | 23 | 167 / 180 | 12:49 | 402 % | 16.3 |
| rt2 | 100 | 55 | 184 | 22 | 185 / 182 | 12:20 | 521 % | 13.8 |

Calibration pass-2 consume time fell from 54 s (rt1) to 36 s (rt2). Scoring fell from 420 s (mem5) to 347 to 367 s. The machine was in interactive use during all three runs, and the three arms ran one after another, not interleaved, so per-stage differences of 10 to 20 % are within noise and are not attributed. The rt2 peak RSS of 13.8 GB against 16.3 to 16.6 GB was not expected from these changes and is not claimed: the peak is set during calibration, whose transient allocations depend on timing. Byte identity is the result. Speed is measured on Codon.

## Evidence

`evidence/{mem5,rt1,rt2}/` hold the exact run script, binary SHA256, empty `diff -rq` output (`identity.txt`), OpenDIA log and GNU time per arm, plus the perf stage tables for mem5. `container-evidence/` holds the `opendia-rt2-20261125.sif` build script and hashes. Local, Docker and SIF OpenDIA SHA256 are all `1587283d3c6fe1d705d145576b44d60ec2d9803a4a8eca1061419dc7515b6794` (the rt2 binary). SIF SHA256 is `775a7bcd493820b513abee1b4f81b2f7e692e5862649009e583fd0be17092e1f`. Container-native `OpenDIAWorkflow_test`, `OpenDIALibraryParquet_test` and `OpenDIAScorer_test` pass.
