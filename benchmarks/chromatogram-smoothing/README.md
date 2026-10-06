# Chromatogram smoother benchmark for OpenMS issue #10425

This benchmark compares OpenMS' ModifiedSincSmoother with SavitzkyGolayFilter on the same known synthetic traces, tracked chromatogram fixtures, and extracted traces from one public DIA run. It tunes each smoother independently, holds synthetic noise seeds out of tuning, applies identical downstream peak-picking settings, and records per-trace results.

This is a reproducible benchmark contribution; it does not change OpenMS production code.

## Reproduce the benchmark

From this directory in PowerShell, create a local virtual environment and run the pinned pyOpenMS wheel:

~~~powershell
py -3.14 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe benchmark_smoothers.py
~~~

The wheel supplies OpenMS' C++ implementations, so no OpenMS project build is needed for this benchmark. The output records the Python, NumPy, pyOpenMS, and checkout versions. The harness verifies that the smoother and peak-picker source files in the checkout match the wheel's OpenMS source revision.

The tracked fixtures are hashed in results/summary.json. Synthetic traces use 401 points at 5.26686-second intervals, matching the median interval measured in the DIA extraction. Peak centers and widths below are specified in samples; the generated RT axis is in seconds. All traces have a baseline of 100 and seeded additive Gaussian noise.

The 64 DIA fragment chromatograms used for the real-data comparison are included as `inputs/extracted_fragment_chromatograms.npz` (748,418 bytes). Its adjacent JSON file records the source mzML checksum, extraction settings, and archive SHA-256. This compact input lets you repeat the smoother comparison without downloading the original run.

| Shape | Peak center, height, sigma (samples) | Noise SD |
| --- | --- | ---: |
| Narrow | (200, 1000, 2.6) | 24 |
| Broad | (200, 1000, 16) | 24 |
| Overlapping | (188, 950, 6.5), (212, 720, 8.5) | 22 |
| Low intensity | (200, 75, 6) | 18 |

Training uses seeds 101, 202, 303, and 404 per shape; held-out validation uses seeds 1001 through 1012 per shape.

## Parameter selection and measurements

Each method gets one global setting selected independently by its mean normalized RMSE across all 16 training traces. The validation seeds do not participate in selection.

| Smoother | Search range | Selected setting |
| --- | --- | --- |
| Modified sinc | is_ms1 false/true; degree 4, 6, 8, 10; m 2, 3, 4, 5, 7, 9, 11, 13, 15, 19, 23, 27, omitting degree-invalid pairs | is_ms1=true, degree 4, m=19 |
| Savitzky-Golay | Frame length 5, 7, 9, 11, 15, 21, 31; polynomial order 2, 3, 4 where less than frame length | Frame length 15, polynomial order 2 |

The is_ms1 entry is the smoother's parameter value; it does not describe the MS level of the extracted chromatograms. Every candidate and score is in results/parameter_search.csv.

Synthetic fidelity metrics are computed against the known clean signal. Apex RT and intensity errors compare each output peak with the clean signal's local maximum in its assigned region. For overlapping peaks, regions meet at the midpoint between their known centers. Width is measured at half prominence within each region, using the larger edge intensity as its local baseline; the reported error compares the filtered width with the clean-signal width. Area error compares the integrated signal above the known baseline with the sum of the known Gaussian areas. Noise reduction is the RMS-residual reduction in samples more than four standard deviations from every peak.

Raw and smoothed inputs use the same PeakPickerChromatogram configuration: corrected mode, 1.0 signal-to-noise threshold, use_gauss=false, and a 3-point/order-2 internal Savitzky-Golay stage. The resulting picked counts, integrated intensities, and median picked FWHM are a downstream sensitivity check, not a production-parameter comparison or ground-truth quantification.

## Synthetic validation results

Values below are medians over the 12 held-out seeds per shape. Apex and width metrics are per-trace means across the known component peaks before taking the median across seeds.

| Shape | Smoother | Apex RT error (s) | Apex intensity error | Area error | Width error |
| --- | --- | ---: | ---: | ---: | ---: |
| Narrow | Modified sinc | 0.0 | 24.97% | 1.98% | 48.22% |
| Narrow | Savitzky-Golay | 0.0 | 21.11% | 0.52% | 37.27% |
| Broad | Modified sinc | 0.0 | 0.63% | 0.35% | 1.33% |
| Broad | Savitzky-Golay | 5.27 | 0.66% | 0.35% | 1.28% |
| Overlapping | Modified sinc | 0.0 | 1.64% | 0.49% | 2.36% |
| Overlapping | Savitzky-Golay | 0.0 | 1.46% | 0.48% | 2.16% |
| Low intensity | Modified sinc | 5.27 | 6.02% | 12.51% | 8.68% |
| Low intensity | Savitzky-Golay | 5.27 | 7.21% | 12.31% | 6.55% |

| Shape | Modified sinc normalized RMSE | Savitzky-Golay normalized RMSE | Modified sinc noise reduction | Savitzky-Golay noise reduction |
| --- | ---: | ---: | ---: | ---: |
| Narrow | 0.0284 | 0.0237 | 59.7% | 60.2% |
| Broad | 0.0087 | 0.0094 | 62.8% | 60.2% |
| Overlapping | 0.0091 | 0.0095 | 62.8% | 60.2% |
| Low intensity | 0.0872 | 0.0944 | 62.9% | 60.1% |

There is no universal winner in these synthetic results. Savitzky-Golay is closer on the narrow trace; modified sinc has lower normalized RMSE on broad, overlapping, and low-intensity traces. Both reduce noise, but the narrow-peak intensity and width errors are substantial. With the fixed picker settings, median picked counts on narrow traces were 12 raw, 20 modified sinc, and 31.5 Savitzky-Golay; smoothing can create additional downstream picks under this configuration. During some downstream picker calls, the default SignalToNoiseEstimatorMedian warned that its highest histogram bin was reached. The same picker settings were retained for every method; treat the picked outputs as sensitivity indicators, not detection-accuracy or validated quantification results. Inspect results/synthetic_validation.csv for every seed, including picked intensity and width.

## Tracked chromatogram fixtures

The benchmark processes every eligible chromatogram in these tracked mzML fixtures:

- NoiseFilterSGolay_2_input.chrom.mzML
- OpenSwathWorkflow_17_output.chrom.mzML
- OpenSwathWorkflow_22_output.chrom.mzML
- OpenSwathWorkflow_23_output.chrom.mzML

results/real_chromatograms.csv reports each fixture chromatogram's point count, RT spacing, picked count, picked integrated intensity, median picked FWHM, and trace area for raw and smoothed inputs. Total picked counts were unchanged on the three OpenSwathWorkflow fixtures. On the NoiseFilterSGolay fixture, counts changed from 12 raw picks to 8 with either smoother.

## Full-size DIA validation

The optional raw-data extraction phase checks how the synthetic-trained settings behave on a long real acquisition. It uses **Fig2HeLa-4h_MHRM_R01_T0.mzML** from the [DIA-CLIP Zenodo record](https://zenodo.org/records/18863866), which identifies PXD005573 as its source. This is one 2,368,588,188-byte run file, not a multi-run dataset. Its MD5 19b814e1bcc9b67afbdac6624428eb31 matches the published archive checksum. The original run is kept outside the repository because it is 2.37 GB.

The extracted comparison input is included in `inputs/` so the reported DIA analysis can be rerun directly:

~~~powershell
.\.venv\Scripts\python.exe .\benchmark_full_dia.py --chromatograms .\inputs\extracted_fragment_chromatograms.npz --repeats 3
~~~

This loads the included RT axis, fragment m/z targets, and 64 intensity traces; it verifies the archive checksum against `inputs/extracted_fragment_chromatograms.json`, then reruns both smoothers and the downstream picker. It refreshes `results/full_dia_summary.json` and `results/full_dia_chromatograms.csv`. The included NPZ is a derived, compact input from the cited public run, not the original spectra. The original mzML is 2.37 GB and remains downloadable separately when you want to reproduce extraction too.

Download and verify the file in PowerShell:

~~~powershell
$dataDir = Join-Path $env:TEMP "OpenMS-DIA-benchmark-data\PXD005573"
New-Item -ItemType Directory -Force $dataDir | Out-Null
$file = Join-Path $dataDir "Fig2HeLa-4h_MHRM_R01_T0.mzML"
curl.exe --fail --location --retry 3 --continue-at - --output $file "https://zenodo.org/api/records/18863866/files/Fig2HeLa-4h_MHRM_R01_T0.mzML/content"
Get-FileHash -Algorithm MD5 $file
~~~

Run the external phase from this benchmark directory, reusing the venv created above:

~~~powershell
.\.venv\Scripts\python.exe .\benchmark_full_dia.py --input $file --output-dir $dataDir --targets 64 --repeats 3
~~~

benchmark_full_dia.py uses OnDiscMSExperiment to stream the indexed mzML rather than loading the 2.37-GB file into memory. It selects the most frequent MS2 isolation window (center 364.5 m/z, lower/upper offsets 14.5 m/z), scores up to 128 strong peaks per scan into 0.02 m/z bins, then selects the 64 highest cumulative-intensity m/z bins separated by at least 0.25 Da and extracts those fragment-ion traces from the same 2,892 scans using 10 ppm tolerance with a 0.005 Da minimum. This favors strong recurring fragments rather than an abundance-stratified sample; low-intensity fidelity is covered by the synthetic traces. The resulting traces have 2,892 points and a median RT interval of 5.267 seconds (3.67% interval CV). Both filters process the identical extracted arrays using their synthetic-trained settings. When extraction starts from the raw mzML, the compressed arrays and provenance sidecar are saved under the external output directory. With the included NPZ, the script reads that file directly. In both cases per-trace and aggregate results are written under results/.

| Full-DIA result over 64 traces | Modified sinc | Savitzky-Golay |
| --- | ---: | ---: |
| Picked peaks (raw total: 4,800) | 3,217 | 4,136 |
| Median picked-peak count change per trace | -29 | -17 |
| Median picked integrated-intensity change vs raw | +12.52% | +13.69% |
| Median picked FWHM change vs raw | +208.44% | +169.89% |
| Median filter time per 2,892-point trace | 0.0327 ms | 0.0171 ms |
| Filter-worker working-set increase | 2.44 MB | 2.25 MB |

On these traces, both methods substantially changed peak-picking output and broadened median picked widths; modified sinc took about twice as long per trace. The mzML has no known clean chromatogram, so the real-data area and width figures are changes relative to raw traces, not evidence of improved true quantification. This covers one run and one isolation window, not multiple instruments or a complete OpenSWATHWorkflow execution.

## Runtime, memory, and result files

On the 100,000-point synthetic stress trace, median filter times were 2.18 ms for modified sinc and 0.92 ms for Savitzky-Golay; process working-set increases were 1.54 MB and 1.55 MB, respectively. The external-DIA worker snapshots are reported separately above. All memory readings are process-level Python/OpenMS working-set and private-byte measurements, not allocation attribution. These timings are from one Windows machine and are not workflow-level performance guarantees.

- results/summary.json: environment, source revision, fixture hashes, selected settings, and synthetic RMSE summaries.
- results/parameter_search.csv: all candidate settings and training scores.
- results/synthetic_validation.csv: per-seed held-out fidelity and picker metrics, including raw baselines.
- results/real_chromatograms.csv: per-chromatogram fixture results.
- results/performance_memory.csv: per-process filter timing and memory readings.
- results/full_dia_summary.json: external dataset provenance, extraction details, and aggregate measurements.
- results/full_dia_chromatograms.csv: per-trace external DIA results.
- inputs/extracted_fragment_chromatograms.npz and its JSON sidecar: the exact 64-trace DIA input used for the comparison, with source and extraction provenance.
- benchmark_smoothers.py and benchmark_full_dia.py: reproducible harnesses.

The original 2.37-GB mzML is not committed; the compact extracted chromatogram input is included. The NPZ can be regenerated from the original mzML with the extraction command above.

## How this matches issue #10425

- **Same chromatograms:** both smoothers run on the same synthetic traces, each tracked fixture trace, and the same 64 DIA XICs.
- **Independent parameter tuning:** separate grids are searched on training seeds; chosen parameters and all scores are recorded.
- **Signal fidelity and noise:** held-out synthetic traces have known clean signals and report apex, area, width, RMSE, and noise metrics for all four requested peak shapes.
- **Downstream peak-picking and quantification effects:** the same picker settings are applied to raw and smoothed traces; picked counts, integrated intensities, and widths are saved per trace. These intensities are proxies, not validated peptide quantification.
- **Runtime and memory:** repeated filter timings and process-memory readings cover fixtures, a 100k-point stress trace, and the extracted DIA traces.
- **Inputs and repeatability:** generated synthetic traces are seeded; tracked mzML fixtures are identified by hash; the extracted DIA arrays and source/extraction metadata are included and checksum-verified; the original run has a source link, published checksum, and download command; parameter ranges and per-seed/per-trace results are included.

Limitations: the external-data comparison uses one run and one isolation window, and it has no ground truth; the tracked files are small regression fixtures; synthetic peaks are Gaussian with additive noise. This is a comparative benchmark and does not establish a universally preferable smoother or predict all production workflows.
