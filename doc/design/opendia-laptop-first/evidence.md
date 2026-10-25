# Historical evidence and baseline selection

[Plan index](README.md) | [Architecture](architecture.md)

## Current verification status (updated 25 October 2026, searches ran 24 October 2026)

The selected clean merged source is `3e47bc7e7468656d0d53237f11fc4f049d278efe`,
including #10642 and #10698, with harness `5a804f5e7c92f5c7e541050123a66f2978d23368`
pinned separately. Full fresh A/B TimSim target-only/mixed and all-three Astral
output comparisons passed all 12 files with zero numeric or scientific-identity
differences. Results rows per repeat are 11,484 / 6,783 / 6,274 and standard accepted
counts are 6,000 / 5,919 / 4,343. Independent provenance checks and final report/package
review remain separate from output equality. The separate profiler branch
now builds and passes its profiling fixture plus all 29 OpenDIA TOPP checks. Four
instrumented cold/warm followups produced valid stage summaries and passed frozen
output-equivalence comparisons across results and all three matrices. Standalone
TimSim rescoring completed. Astral rescoring completed on a disposable copy of the
frozen database and its WAL/SHM sidecars. It exited successfully in 9m09s with
10.27 GiB maximum RSS. Checks confirmed that the frozen source files did not change.
This single rescore is descriptive, not a baseline, repeatability result or acceptance
evidence. These are separate-binary measurements, not same-binary repeats, exact
resource peaks for the full workflow, or candidate acceptance.
See [WP-00](work-packages/00-acceptance-baseline.md) for exact matrix counts, frozen
comparison tolerances, remaining evidence and proposed gates awaiting owner/reviewer agreement.

The declared workstation protocol uses 17 OpenDIA threads rather than the
historical four. Engine and harness provenance must be independently captured.
Dedicated verified scope `cgroup_final` kernel peaks include file cache and keeper
overhead; sampled RSS/cgroup counters are diagnostics. Inner workload timing is
separate from outer scope/keeper timing. Sampled apparent output/temp sizes are
lower bounds, not exact filesystem peaks. See the [measurement contract](https://github.com/singjc/OpenMS/blob/9abc264f4c95cb550d67490030dcac649750060c/tools/scripts/opendia_benchmark/README.md#dedicated-workload-accounting-and-limitations).
These semantics do not establish resource-target compliance or exact peak disk.
The two full Astral baseline searches used 13,875.97 and 13,858.52 MiB
single-process peak RSS, above the proposed 4 GiB process target. Their dedicated
cgroup peaks were about 20.5 and 22.4 GiB. The 2 GiB low-memory mode was not tested.
Sampled disk values are lower bounds and do not prove peak-disk compliance.

The 21 October dirty-source searches and failed repeats remain historical. The
22 October premerge same-binary comparisons reported in #10662 passed all 12
output comparisons, without proving cross-thread/platform/cache invariance or
rerunning merged develop. #10698 merged on 23 October 2026 and corrects Parquet
scientific score values, so affected earlier scores are not an equivalent
corrected-source baseline. Preserve original logs and failures rather than
retrospectively replacing them with successful corrected results.

## Historical project baseline selection for WP-00

Use pinned upstream/CI builds of develop with [#10642](https://github.com/OpenMS/OpenMS/pull/10642), the ProteoBench Astral entrapment run and one diaPASEF dataset. Timo's [plan review](https://github.com/OpenMS/OpenMS/issues/10661) reports approximately 33 min wall time and 8.0 GB peak RSS for 600k precursors, 300 windows and four threads: 15 min extraction/scoring, 13.5 min rescoring/inference/export, and 5 min loading/calibration. The working OSW contains 2.9M features and 17.3M `FEATURE_TRANSITION` rows. These rounded review figures motivate early WP-08 work, not a new measurement by this documentation update. WP-00 must reproduce stage time/RSS with pinned inputs and build provenance.

## Historical local observations

The observations below are retained for reproducibility, not adopted as the project baseline. The tiny fixture and mixed local binaries cannot establish representative performance or scientific quality.

### Provenance and interpretation

These observations were recorded during the 20 October 2026 design investigation. The inspected source was `6b6c341434`; available OpenMS binaries reported `e85d291`, and the installed extraction prototype reported version `0.2.14`. These are baseline experiments, **not current-source or greenfield-engine validation**. No production code or defaults were changed during that investigation.

The original detailed reports and logs are local scratch artifacts, not dependencies of this portable plan. WP-00 must archive reproducible commands, input checksums, executable/library provenance and outputs in an agreed durable location before using them as release evidence.

## Measured configuration improvement

The actual `DIATraceExtractor` CLI was run nine times, interleaving three repetitions at each buffer size. Input was the checked-in `src/tests/topp/OpenSwathWorkflow_17_input.mzML` and matching `.tsv`. Settings included scheme B, 250 ppm MS2 and 25 ppm MS1 tolerances, identity RT, unrestricted RT window, four threads and one precursor per library part.

| Requested buffer | Median peak RSS (KiB) | Median elapsed (s) | Output check |
|---|---:|---:|---|
| 512 MiB | 580,384 | 0.12 | All 600 rows equal to reference |
| 64 MiB | 121,564 | 0.04 | All 600 rows equal to reference |
| 16 MiB | 72,420 | 0.04 | All 600 rows equal to reference |

All fields and values were compared after row sorting against a one-thread nonzero reference. Each output had 12 traces, 600 points and 294 nonzero intensities. The 16 MiB setting reduced measured median RSS by **87.52%**. This supports demand-sensitive committed memory rather than eagerly materializing maximum capacity.

This tiny fixture is dominated partly by setup/pool initialization. It does not prove a 3x large-data throughput gain, automatic adaptation, or a whole-engine cap. A separate 1 MiB request actually created a 4 MiB pool due to minimum allocation granularity, so user settings must not be mistaken for enforced limits.

## Additional public check: page and worker scheduling

A follow-up on 20 October tested the hypothesis that library-page boundaries and worker scheduling can change without changing the extracted scientific evidence. The installed `DIATraceExtractor` was executed 18 times: `library_part_rows` of 1, 2 and 65,536, one/four workers, and three interleaved repetitions of each combination. The buffer stayed at 16 MiB and all scientific extraction settings and inputs matched the nonzero experiment above.

| Library part rows | Workers | Median peak RSS (KiB) | Median elapsed (s) | Output |
|---|---:|---:|---:|---|
| 1 | 1 | 72,420 | 0.04 | Exact reference equality |
| 1 | 4 | 72,504 | 0.03 | Exact reference equality |
| 2 | 1 | 72,456 | 0.04 | Exact reference equality |
| 2 | 4 | 72,576 | 0.03 | Exact reference equality |
| 65,536 | 1 | 72,460 | 0.03 | Exact reference equality |
| 65,536 | 4 | 72,712 | 0.03 | Exact reference equality |

For every execution, the full TSV header and all 600 rows matched the previous one-thread reference exactly after row sorting, including all metadata fields and 294 nonzero intensity points. Actual CLI logs showed 6/3/1 library parts and read/plan tasks for part sizes 1/2/65,536 respectively. Part size 1 produced two extraction tasks with up to two concurrent workers; the larger settings produced one extraction task. The trace pool peaked at two versus one blocks, with no heap-fallback units.

**Interpretation:** this exercises actual reader/planner/extractor/writer boundaries and supports scheduling independence for this small public workflow. It does not show a meaningful throughput improvement, many-worker scalability, bounded large-library memory, or a mobility-preserving/peptidoform-aware greenfield implementation. The design's whole-engine gates remain open. WP-04/WP-05/WP-09 should extend this check to representative data rather than assume the tiny fixture is sufficient.

### Reproduction details

The originating machine retains `opendia-page-experiment.sh` and `opendia-page-*.{tsv,stats.json,resources,log}` under `/home/sing/.jcode/scratch`. The per-run command below is the tested public interface. Set `EXTRACTOR` to the pinned installed executable, `OUT` to a fresh writable directory, `PAGE` to 1/2/65536, `THREADS` to 1/4 and `REP` to 1/2/3. Run from the repository root, preserving outputs for each combination:

```bash
/usr/bin/time -f 'elapsed_s=%e user_s=%U system_s=%S max_rss_kib=%M' \
  -o "$OUT/page-${PAGE}-t${THREADS}-r${REP}.resources" \
  "$EXTRACTOR" \
  -in src/tests/topp/OpenSwathWorkflow_17_input.mzML \
  -lib src/tests/topp/OpenSwathWorkflow_17_input.tsv \
  -out_tsv "$OUT/page-${PAGE}-t${THREADS}-r${REP}.tsv" \
  -out_stats "$OUT/page-${PAGE}-t${THREADS}-r${REP}.stats.json" \
  -scheme B -mz_tol 250 -mz_tol_ms1 25 -rt_mode identity -rt_window -1 \
  -threads "$THREADS" -library_part_rows "$PAGE" -buffer_mb 16
```

Captured SHA-256 identities:

```text
mzML:      d95a281e78aa440f3081d441d21094a8ba2cd8a724fb55bd0b30737e368c2d00
TSV:       92d726b4e158b0a8c2fa4ccf2fab08230256b031ef1966cedfcd9b7532a5838a
executable:76e14ffa4634480800808a0d5ffb4d50921488c07a4bdbce47c64d7546933033
```

The executable hash alone does not pin its dynamically linked dependencies. A portable benchmark release still needs the complete environment manifest required by WP-00.

## Public scientific and output checks

| Requirement | Observed check | Result and limitation |
|---|---|---|
| OpenDIA SQLite/Parquet outputs | Seven selected success/help CTests | Passed; the small exported tables agreed. Not a large-library or localization benchmark |
| Mobility extraction/scoring | `TOPP_OpenSwathWorkflow_17` and `_17_cache` | Both producers succeeded with two features, and chromatogram comparisons passed |
| PASEF geometry/scoring | `TOPP_OpenSwathWorkflow_23` | Producer succeeded with three features and nonzero IM score fields; chromatogram comparison passed |
| Feature-output compatibility | `_17_out1`, `_17_cache_out1`, `_23_out1` | Failed because emitted featureXML lacked the IdentificationRun block expected by current references. Older binary/newer reference mismatch is plausible, not a proven root cause |
| Baseline IPF and export | `TOPP_OpenSwathInfer_1_prepare`, `_1`, `_1_export`, `_1_out1` | Passed with 333 exported data rows. This fixture disables optional MS1/MS2 precursor scoring |
| Full precursor-evidence IPF | Public `OpenSwathInfer` with both MS1/MS2 switches enabled, one and four threads | 291 exactly equal SCORE_IPF rows, including IDs and PEP/QVALUE fields; public `OpenSwathExport` exported 291 rows |
| Missing IM/localization prerequisites | OpenDIA negative tests | Rejected missing required evidence as expected. Negative tests are not successful identification workflows |
| Storage/localization capability | `TOPP_OpenDIA_parquet_peptidoform_precondition_failure` | Confirmed rejection of the unsupported Parquet/IPF path. New canonical storage must avoid this capability split |

Full-evidence IPF used copies of `src/tests/class_tests/openms/data/PyProphet_inference_test.osw`. It consumed 294 precursor/IPF rows and 4,655 transition rows, expanded to 9,321 hypotheses. One/four-thread runs each took 0.07 s, with 49,088/49,100 KiB peak RSS. This demonstrates determinism on stored evidence, not joint mobility/IPF acquisition processing, large-family scaling or true localization accuracy. Disabling scientific evidence changed output counts and is not an equivalent optimization.

## Extraction correctness caveat

An initial cross-engine comparison produced only zeros and was rejected as accuracy evidence. A widened, matched-policy comparison produced 294 nonzero points. Existing OpenSWATH wave versus legacy batching agreed exactly on all 600 values/times; prototype thread configurations also agreed.

The two extractors did **not** agree with one another: the largest difference was approximately 20.042 intensity units. A raw-spectrum sum of 120.0000002682 became about 119.2002182 in existing OpenSWATH output, consistent with omitting the first centroid of intensity 0.8000000119 plus output encoding effects. Source inspection confirmed a left-walk first-centroid edge case. No fix was made during that investigation. Test first/last centroids, one-element spectra and targets beyond the final centroid against explicit mathematical semantics, not only legacy output.

Timo additionally reported beyond-last-centroid double counting in both tophat overloads. These historical bugs were subsequently fixed in standalone core PR #10673, included in the selected source, with hand-computed IM/non-IM edge tests. Do not add a legacy-bug compatibility mode. This merge status does not retroactively validate the historical measurements above.

## Historical local build limitations

- The standalone picker/analyzer test selection had 24 failures out of 29, including producer aborts with `stack smashing detected`. Passing comparisons against stale outputs were not accepted as success.
- Executables and shared libraries had different build dates. A stale ABI was suspected, not established.
- A coherent rebuild was blocked by a configured missing Thermo/.NET `libnethost.a` under the .NET 8.0.31 host pack. The design exercise did not bypass or modify dependencies to hide this.
- These are local-environment limitations, not prerequisites for WP-00. A coherent rebuild and rerun would be necessary only to promote these particular historical results to current-source evidence. Use upstream/CI provenance for the project baseline instead.

## Local investigation archive

The following paths document where the original investigation is stored on the originating machine. They are not portable links and may be cleaned up:

```text
/home/sing/.jcode/scratch/opendia-laptop-first-greenfield-2026-10-20.md
/home/sing/.jcode/scratch/opendia-evidence-pipeline-design-2026-10-20.md
/home/sing/.jcode/scratch/opendia-pool-experiment.sh
/home/sing/.jcode/scratch/opendia-pool-*.{resources,log,tsv}
/home/sing/.jcode/scratch/opendia-ipf-full-*
/home/sing/.jcode/scratch/opendia-mobility-ipf-baseline-ctest.log
/home/sing/.jcode/scratch/opendia-standalone-stages-ctest-20261020.log
```

## Not yet demonstrated

Whole-engine 2/4 GiB limits, representative large-library speed, candidate-retrieval recall, joint mobility/localization quality, localization error calibration, oversized-family behavior, and release packaging remain unverified for the proposed engine. The [work packages](README.md) turn these gaps into owned implementation and acceptance work rather than treating baseline results as proof.

The proposed `results` and `rescore` profiles and capability-matched output-size reductions remain unmeasured. WP-00/WP-08/WP-10 must agree numerical footprint targets and test retained artifacts before claiming improvement. Detailed disk quotas and diagnostic profiles are deferred beyond v1.
