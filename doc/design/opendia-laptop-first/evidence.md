# Baseline evidence and unresolved limits

[Plan index](README.md) | [Architecture](architecture.md)

## Provenance and interpretation

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

The two extractors did **not** agree with one another: the largest difference was approximately 20.042 intensity units. A raw-spectrum sum of 120.0000002682 became about 119.2002182 in existing OpenSWATH output, consistent with omitting the first centroid of intensity 0.8000000119 plus output encoding effects. Source inspection confirmed a left-walk first-centroid edge case. No fix was made. Test first/last centroids, one-element spectra and targets beyond the final centroid against explicit mathematical semantics, not only legacy output.

## Build and execution blockers

- The standalone picker/analyzer test selection had 24 failures out of 29, including producer aborts with `stack smashing detected`. Passing comparisons against stale outputs were not accepted as success.
- Executables and shared libraries had different build dates. A stale ABI was suspected, not established.
- A coherent rebuild was blocked by a configured missing Thermo/.NET `libnethost.a` under the .NET 8.0.31 host pack. The design exercise did not bypass or modify dependencies to hide this.
- WP-00 must restore a coherent local build and rerun affected public workflows before interpreting these results as current-source evidence.

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
