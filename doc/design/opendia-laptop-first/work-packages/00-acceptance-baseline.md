# WP-00: Acceptance baseline and datasets

- **Track:** baseline.
- **Status:** proposed, full workstation A/B output comparisons passed, public evidence release blocked by archive path cleanup, updated 25 October 2026. **Owner / reviewer:** unassigned.
- **Parent epic:** [#10661](https://github.com/OpenMS/OpenMS/issues/10661). **Issue:** [#10662](https://github.com/OpenMS/OpenMS/issues/10662).
- **Shared policies:** [architecture](../architecture.md) and [plan](../README.md).

## Objective

Establish reproducible scientific, runtime and storage comparisons without blocking independent implementation.

## Scope

- Pin an upstream/CI `develop` revision and the applicable #10642 revision/state, commands, dependencies and input checksums.
- Pin the established TimSim and exploratory Astral inputs below. Select separate known-site data for localization validation, which neither current dataset supplies.
- Record cold preparation and warm reuse, stage/total runtime, peak memory, actual peak disk and output row cardinalities.
- Inventory OSW `FEATURE_TRANSITION` fanout and rescoring cost for WP-08, alongside identification, localization and quantitative error.
- The local 600pt fixture is a regression aid, not the baseline. Repairing local `libnethost` is not a prerequisite.

## Established datasets and scripts

| Dataset | Current search and validation scope |
|---|---|
| TimSim fixture `validation_3x3` | Six diaPASEF inputs: three controls and three treatments. Run target-only and target-plus-entrapment transition libraries against packaged realized identification/abundance and entrapment truth. |
| Astral HeLa 50 ng, 15 min | Three technical replicates, currently searched with a randomly sampled 100k-peptide FASTA. Exploratory identification, within-run entrapment diagnostics and repeatability only, not official full ProteoBench global precursor validation or biological protein performance. No abundance truth is supplied. |

Use [the benchmark scripts](https://github.com/singjc/OpenMS/blob/9abc264f4c95cb550d67490030dcac649750060c/tools/scripts/opendia_benchmark/README.md): `build.sh`, `run_timsim.sh`, `run_astral.sh`, `validate_timsim.py`, `metrics_astral.py`, `make_report.py` and `render_report.py`. The runners capture commands, source/binary/input provenance, elapsed time, process peak RSS and sampled output-directory disk use. `monitor_host.py` and optional systemd memory/swap limits add system-resource context, not proof of whole-process/cgroup peaks or stage attribution. Machine, total RAM/logical threads, available RAM and idle CPU are now captured per search. Fresh target-only/mixed six-file TimSim and three-file Astral searches completed successfully on 21 October 2026 at source `d3b00d382efc`, with fresh logs, metrics and a regenerated report/PDF. TimSim target/mixed wall times were 6 min 55 s / 11 min 11 s with 0.84 / 0.85 GiB process peak RSS. Astral took 3 h 8 min 31 s with 10.41 GiB process peak RSS. All used four configured threads on a 30.80 GiB, 14-thread laptop, with other host activity, so these are not isolated speed comparisons. Results and reporting artifacts are tracked in [#10662](https://github.com/OpenMS/OpenMS/issues/10662). Repeat-tolerance, stage/cgroup/peak-disk and scientific acceptance gates below remain open.

The standard view requires 1% run-level precursor, global peptide and global protein q-values. Astral FDP diagnostics instead use run-level precursor `m_score` without the extra peptide/protein filters. Global precursor q-values and the exact modified-peptide pair mapping are unavailable, so global precursor validation and paired FDP remain blocked. Combined FDP also requires an independently verified effective entrapment:target search-space ratio, which the sampled FASTA does not establish. Neither current dataset closes the known-site localization gate.

## Historical reproducibility package (21 October 2026)

Published final package and failure summary: [WP-00 comment](https://github.com/OpenMS/OpenMS/issues/10662#issuecomment-6323628145).

`package_reproducibility.py` builds a local bundle with completed-job commands,
logs, reference TSVs, metrics, optional report/evidence, dataset-relative input
hashes, a script snapshot, replay instructions and a SHA-256 inventory. Its
verifier checks file sizes, hashes and the complete inventory. Raw inputs,
intermediates, binaries and dependencies are not redistributed. Generated bundles
remain local, not committed to this planning branch.

Final local archive: `opendia-reproducibility-20261021-complete.tar.gz`
(10,511,211 bytes), SHA-256
`f31ba2a4e39f9fe255cd47961bcf29e0a5dcad2df2665712869166faf88d1bb0`.
It includes all three completed independent repeats and their failed comparisons.
Extraction into another directory, bundled verification and GNU checksum checks
passed. Replaying all three comparisons from the extracted archive reproduced
the complete saved file diagnostics and expected failed verdicts. All 43 Python
regression tests pass. Package acceptance is complete, not scientific acceptance.

`compare_repeats.py` compares scientific row identities and accepted identities,
aligns matrix rows and run columns by name, and checks scores, quantities and RT
against a policy declared before the fresh repeat. It ignores volatile database
IDs, not scientific changes. This is output comparison, not a substitute for
input/settings/dependency provenance checks.

The first independent target-only TimSim repeat retained the same 6,000 accepted
precursor identities, but failed the frozen numerical/raw-identity checks:
raw rows changed from 11,485 to 11,458 and 28 of 6,000 precursor abundance cells
exceeded tolerance. Captured executable/library hashes, resource limits and
commands match after substituting only the new output-root path. Preserve this
failure, do not widen tolerances or mark the repeat gate complete. The mixed
TimSim repeat also completed, but standard-filter accepted identities changed
from 5,949 to 5,935 (18 missing, 4 added), and raw rows from 6,682 to 6,709.
Its captured commands, executable/library hashes and limits also match. The
three-file Astral repeat completed successfully in 3 h 13 min 10 s, but failed
strict comparison: raw rows 6,236 to 6,206, standard-filter accepted identities
4,337 to 4,291 (61 missing, 15 added). These identity totals include entrapment,
not only the target-only report view. Astral captured command, binary, RAW/FASTA
hashes and limits match, and 57 dependency/model hashes remained unchanged
between the pre-Astral and post-repeat captures. This does not recover missing
historical baseline hashes. All three strict repeat comparisons fail.

Read-only prepared-library inspection proves 998 of 1,000 target-only TimSim
decoy sequences differ, while all 1,000 target precursor identities match.
In that historical build, default transition-library decoy shuffling was
time-seeded, without a wired general shuffle-seed option, so the Percolator seed
alone did not establish deterministic inputs to rescoring. #10698 subsequently
addressed the seed and ordering paths; see the current verification section.

The environment supplement pins cached upstream `develop`
`f9f81826b9bd4c8dcb7166af709f4e39b23108b1` and applicable #10642 revision
`2d5f4a76d5837e146e124949bc420fba802b90cd`, both ancestors of the measured source.
There is no captured CI-run attestation. The original dirty-tree changes and
historical dynamic-library/model hashes cannot be reconstructed from the source
commit alone. Packaging-time hashes are explicitly labelled, not retrospective
proof. These historical provenance gaps cannot be closed by a new successful run.
The current-source merge criteria require separately captured evidence below.

## Current merged-source verification (24 October 2026)

The selected clean source pin is `3e47bc7e7468656d0d53237f11fc4f049d278efe`.
It contains #10642 and the repeat-determinism fixes from #10698, merged on
23 October 2026. The independently frozen harness pin is
`5a804f5e7c92f5c7e541050123a66f2978d23368`. Full fresh A/B TimSim target-only,
mixed/entrapment and three-replicate Astral output comparisons have now passed.
The final baseline task `5858088rjg` exited 0. This is baseline output evidence,
not candidate-engine acceptance or completed provenance review.

The local evidence root is `/home/singjc/.jcode/scratch/opendia-wp00-workstation-20261024`.
Its `timsim-target_only-comparison.json`, `timsim-entrapment-comparison.json` and
`astral-comparison.json` each report four passing files: results plus precursor,
peptide and protein matrices. Counts below are identical in A and B.

| Pair | Results rows | Precursor-only 1% accepted | Standard three-filter 1% accepted | Matrix rows: precursor / peptide / protein | Numeric cells checked across four files |
|---|---:|---:|---:|---:|---:|
| TimSim target-only | 11,484 | 6,000 | 6,000 | 1,000 / 1,000 / 250 | 426,924 |
| TimSim mixed/entrapment | 6,783 | 5,958 | 5,919 | 993 / 993 / 250 | 257,604 |
| Astral, all three replicates | 6,274 | 4,590 | 4,343 | 1,747 / 1,624 / 1,595 | 240,762 |

All 12 files have zero changed numeric cells, zero missing/added scientific
identities and zero maximum numeric differences. Both accepted sets are unchanged.
These counts include entrapment where present, not only a target-only report view.
The comparator checks binary SHA-256, successful exits and explicit non-aborted
monitor summaries. It explicitly does not establish complete input/settings/model/
dependency equivalence. Independent frozen provenance checks and the report
worker's evidence package remain separate and pending review. Do not infer them
from output equality or modify frozen harness/provenance artifacts.

The separate profiling branch's native build and OpenDIA integration tests passed: the requested native target set built successfully, the profiling fixture passed, and all 29 OpenDIA TOPP checks passed. Four compile-time accesses to a private parameter helper were caught and corrected to use resolved parameters before the successful build. The instrumented binary is separate from the clean merged-source baseline. Four full instrumented runs (TimSim target-only cold/warm and Astral predicted cold/warm) exited successfully, produced valid stage summaries, and matched each frozen baseline output using the same table-comparison policy; the two builds have different binary hashes, which were recorded. These are different-binary output-equivalence checks, not same-binary repeats or candidate acceptance.

Instrumented stage/resource alignment is valid at recorded boundaries. Stage resource maxima are sampled observations, not exact peaks. Standalone TimSim `ms1ms2` rescoring completed on a copied workflow and processed 16,441 rows (8,429 target, 8,012 decoy), writing 13,421,129 apparent bytes. Astral standalone `ms1ms2` rescoring also completed on a disposable copy. It ran for about 9m09s and used 10,515 MiB (about 10.27 GiB) peak process RSS. It reported 4,063 test-set PSMs at q < 0.01. It changed only the copy. Checks before and after confirmed that the frozen SQLite source and its WAL/SHM sidecars did not change. These descriptive counts do not establish scientific acceptance. A checksummed baseline bundle and a followup bundle both passed verification after clean extraction. The followup archive SHA-256 is `0ae95842b063ff37003ba7d4a9a100591deb0fa97bf4ea75080ad24f8142381e`; the baseline archive SHA-256 is `4d45c0177bbcdc76124fe4821579137e9e703dc136b125dfb908220f791257ca`. A safety scan found workstation-specific absolute paths in the archives. Do not publish these files as-is. Remove or replace local paths, rebuild the checksums and manifest, and repeat the clean-extraction check before any public upload. Durable publication remains pending.

The declared workstation protocol uses `THREADS=17`; historical four-thread
laptop timings are context, not matched performance evidence. Record actual
affinity, prediction and Percolator worker settings separately. Pin engine source
and harness revision/status independently, including frozen harness hashes:
a clean engine does not establish harness identity or cleanliness. No fresh
workstation performance-target compliance is asserted by this protocol declaration.

### Resource evidence semantics

See the [benchmark measurement contract](https://github.com/singjc/OpenMS/blob/9abc264f4c95cb550d67490030dcac649750060c/tools/scripts/opendia_benchmark/README.md#dedicated-workload-accounting-and-limitations).
`BENCH_CGROUP=1` or either memory/swap limit requests a dedicated verified user
scope. `host_workload.cgroup_final` is the completed-workload kernel peak read,
including file cache and keeper overhead, not process RSS. Sampled counters are
diagnostic. `resources.txt` measures the inner workload; `scope_resources.txt`
also includes outer scope-launcher/keeper handshake overhead. Sampled process-tree
RSS may double-count shared pages and miss children/peaks; individual observed
`VmHWM` values must not be summed. Sampled apparent output/temp bytes are a lower
bound, not an exact filesystem peak, and temp is a subset of output. These
semantics improve the evidence contract, but do not close the actual-peak-disk
or other pending measurements without new reviewed evidence.

The 22 October evidence reported in #10662 passed all 12 output comparisons for
TimSim target, TimSim mixed and Astral at premerge source `08c200dba981` with zero
changed numeric cells or scientific identities. This establishes only that tested
same-binary/host/settings/thread-count repeat check. TimSim used memory-loaded
cache, Astral disk cache with three concurrent SWATH windows. It does not prove
cross-cache, cross-thread-count or cross-platform equivalence, or constitute a
rerun of the selected merged source. Preserve the earlier failed comparisons.
The Parquet correction changes scientific score values, so affected historical
scores are not an equivalent corrected-source reference.

### Frozen repeat-comparison policy

Reuse the existing policy without widening limits after observing results:

| Value class | Relative tolerance | Absolute tolerance |
|---|---:|---:|
| Scores | `1e-5` | `1e-8` |
| Quantities | `1e-5` | `1e-6` |
| Retention time | `1e-6` | `1e-4` |

Require exact scientific row identities and accepted identity sets, matching
missingness and finite numeric values. Compare raw results and precursor,
peptide and protein matrices. Normalize commands only for repeat output paths,
and independently verify frozen input/model/binary/dependency/harness hashes.
Standard accepted sets use the three 1% filters described above. These filters
are not proof of calibrated 1% error control or localization accuracy.

### Requirement-to-check map after full A/B output comparison

| Merge requirement | Evidence still required for the selected source |
|---|---|
| Pinned provenance | Complete clean-source/build/input/environment manifest and explicit CI status. Keep historical dirty-tree/dependency gaps labelled, not retroactively repaired. |
| Fresh producers and reference artifacts | Full fresh A/B target/mixed TimSim and all-three Astral producers completed and all 12 output comparisons passed. The local evidence archives passed clean-extraction checks, but a safety review found workstation paths and ownership metadata. They are not ready for public release. Four separate instrumented cold/warm followups also passed output-equivalence against their frozen outputs; they do not replace same-binary repeats. |
| Repeat baseline | Output reproduction is demonstrated under the unchanged frozen policy with zero differences. Independent input/settings/model/dependency/harness provenance review remains separate before full requirement sign-off. |
| Complete measurement table | Cold/warm and stage/total runtime, stage memory, process/cgroup peaks including children, swap, rescoring share, table cardinalities, final bytes and observed peak disk including sidecars/staging. Sampled directory sizes and configured limits alone do not close this check. Instrumented TimSim/Astral cold/warm stage records and standalone rescoring for both datasets exist. The Astral rescore ran on a copy because the frozen source has WAL/SHM sidecars. Both local evidence bundles passed clean-extraction verification, but contain workstation-specific absolute paths and are not safe to publish as-is. The two full Astral baseline runs used 13,875.97 and 13,858.52 MiB single-process peak RSS, above the proposed 4 GiB process target. Their cgroup peaks were about 20.5 and 22.4 GiB. The 2 GiB low-memory mode was not tested. Comprehensive peak-disk measurement and resource acceptance remain open. |
| Predeclared quality/resource/size gates | Numerical identification, retrieval, localization and quantitative-error margins, runtime and matched-population results/rescore size limits require explicit review and freezing before candidate-engine evaluation. No such undecided margin is accepted by this update. |

### Owner/reviewer decision matrix

| Decision | Recorded policy or proposal | Status / decision owner |
|---|---|---|
| Repeat comparator | Existing frozen tolerances, exact scientific/accepted identities, matching missingness and finite values | Enforced repeat check; do not widen after results |
| Behavior-preserving candidate refactors | Exact accepted identity sets, numeric results within the frozen comparator, no added missing/nonfinite cells | Proposed candidate gate; WP-00 owner/reviewer agreement required before evaluation |
| Runtime | Median wall time no greater than matched same-host baseline, separating cold preparation and warm reuse | Proposed, not accepted or demonstrated; WP-00/WP-10 owners must agree repetitions and matched configuration before evaluation |
| Output footprint | Final results plus rescore bytes no greater than matched-population baseline | Proposed, not accepted or demonstrated; WP-00/WP-08/WP-10 owners must agree artifact capabilities/population |
| Peak memory and disk | Observed complete-process/cgroup memory and peak disk no greater than matched same-host/configuration baseline | Proposed, not accepted or demonstrated; WP-00/WP-10 owners must approve measurement scope |
| Scientific accuracy and FDP | Report observed identification, abundance and entrapment metrics; no new numerical truth-accuracy or FDP cutoff is invented here | Numerical acceptance margins await scientific owner/reviewer agreement; unavailable truth remains blocked |
| Localization truth | Select and pin separate known-site DIA data and truth mapping, then review exhaustive/IPF and null/decoy experiment | Open decision for WP-00/WP-07 owners; neither existing dataset closes it |

These proposals are recorded before future candidate evaluation, not approved by
being documented. A behavior-changing scientific candidate needs separately
reviewed scientific margins, not automatic exact-set acceptance. No observed
speedup, footprint reduction or scientific calibration is claimed by this matrix.

The existing 4 GiB process and 2 GiB low-memory targets on a 16 GiB laptop remain
targets, not measured caps. Baseline characterization may proceed using the
frozen repeat policy, but candidate scientific/performance acceptance remains
blocked until its additional thresholds are approved and recorded before runs.
Proposals must remain labelled as proposals, with no retrospective relaxation.

Known-site localization remains blocked. PXD004573 is the proposed follow-up:
all 579 exact published modified-peptide strings have been authenticated, so
published-list authentication is no longer a blocker. Proposed raw data have not
been downloaded or locally checksum-verified. Site/S3-S4 mapping, synthetic purity,
modification/isotope/charge/transition joins and nominal-zero/blank interpretation
remain unresolved, not approved per-run truth. Coordinate the exhaustive/IPF and localization
null/decoy experiment with WP-07 (#10669) before freezing scientific contracts.
TimSim supplies realized identification/abundance truth, not known-site truth.
Astral supplies repeatability and exploratory within-run entrapment diagnostics,
not abundance truth or biological-protein validation. Missing global precursor
q-values, exact modified-peptide pairing and independently verified effective
entrapment:target search-space ratio keep the respective Astral checks blocked.

The comparison subcheck is passed, not pending. Full merge sign-off still awaits
the requirement-specific evidence and approvals above, rather than a claim that
all five requirements remain wholly untested.
Complete evidence checks and adversarial review before opening a PR. This
workstation-verification update does not approve candidate acceptance or close
any open scientific, resource or provenance gates.

## Dependencies

- **Start / local merge:** none.
- Supply pinned inputs and thresholds to WP-10. WP-01 and the WP-07 early experiment can proceed independently.

## Local merge criteria

- [ ] A manifest records exact source/CI and #10642 provenance, checksums, parameters and executable commands for both representative datasets.
- [ ] Fresh successful producer logs and reference outputs are attached for Astral and diaPASEF, with no stale-output reuse.
- [ ] A repeat baseline run reproduces identities and scores/quantities within recorded numerical tolerances.
- [ ] The measurement table includes stage runtime and memory, process/cgroup peak memory, rescoring share, table row counts, final bytes and observed peak disk including sidecars/staging.
- [ ] Numerical quality, runtime and results/rescore size gates are recorded before candidate-engine evaluation, with unavailable truth explicitly blocked.

**Integration gates:** WP-09 owns public workflow checks, and WP-10 owns representative quality/resource/performance acceptance.
