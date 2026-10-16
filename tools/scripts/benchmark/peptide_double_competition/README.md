# Peptide-level double-competition benchmark

This benchmark compares the existing OpenMS peptide-level target-decoy calculation with the opt-in `PSM_and_peptide` method on the same Crux Tide PSMs. It uses the UPS1-in-yeast spike-in data recommended in [OpenMS issue #10490](https://github.com/OpenMS/OpenMS/issues/10490) and the method in [Lin et al. (2022)](https://doi.org/10.1021/acs.jproteome.2c00282).

## Dataset and exact inputs

The four mzML files are from the public OpenMS ProteomicsLFQ example data directory:

`https://abibuilder.cs.uni-tuebingen.de/archive/openms/Tutorials/Example_Data/ProteomicsLFQ/`

They cover two UPS1 spike concentrations and two replicate runs. The accompanying target/decoy FASTA is from the same directory. `SHA-256` values below identify the exact downloaded inputs.

| Input | Size (bytes) | SHA-256 |
|---|---:|---|
| `UPS1_2500amol_R1.mzML` | 735183675 | `445A4BE3523EC7033FEFFBDFB7FAC56DC6D2FA6628490A60B38380FE4C5132E3` |
| `UPS1_2500amol_R2.mzML` | 733797294 | `BDF2498BC73AEB46C2ABD3357EF1CCA3AFA045035748A76789E7EA3B1236CC11` |
| `UPS1_50000amol_R1.mzML` | 767043386 | `0A9817E7BDCC54402C5809E54F52A4241F1A52AB9413C25C5873DC4572EB659F` |
| `UPS1_50000amol_R2.mzML` | 767355283 | `947B92F3ABB8CEF22A19055832953843786C19B1EE4170B8DFF7165A1D79F914` |
| `uniprot_yeast_reviewed_isoforms_ups1_crap.fasta_td.fasta` | 9671332 | `516E7CBE43DFEB94E6CE9FE318E6F4B31E0A9076933196C271EE83F2702457CC` |

The mixture contains both yeast background and the 48-protein UPS1 standard. UPS-positive counts below are a known-positive recovery measure; yeast target hits are not treated as false positives.

## Reproduce

Commands below are for PowerShell from the OpenMS repository root. They use Crux Tide `4.3.2-835f0f3-2025-09-02`, CMake-built OpenMS `FalseDiscoveryRate.exe`, and the Python standard library. Download Crux from its [official releases](https://github.com/crux-toolkit/crux-toolkit/releases) and set `$crux` to its `crux.exe` path. Build the Release target with `cmake --build .\build --config Release --target FalseDiscoveryRate --parallel 8` first; adjust the build directory in `$fdr` if needed.

```powershell
$python = 'python'
$crux = 'C:\path\to\crux.exe'
$fdr = '.\build\bin\Release\FalseDiscoveryRate.exe'
$data = 'https://abibuilder.cs.uni-tuebingen.de/archive/openms/Tutorials/Example_Data/ProteomicsLFQ'
$work = '.\peptide-double-competition-benchmark'
New-Item -ItemType Directory -Force "$work\raw", "$work\search", "$work\converted", "$work\index_metadata" | Out-Null

# Download the four mzMLs with checked parallel HTTP byte ranges (sizes are in the table above).
& $python tools/scripts/benchmark/peptide_double_competition/download_ranges.py "$data/UPS1_2500amol_R1.mzML" "$work/raw/UPS1_2500amol_R1.mzML" --size 735183675
& $python tools/scripts/benchmark/peptide_double_competition/download_ranges.py "$data/UPS1_2500amol_R2.mzML" "$work/raw/UPS1_2500amol_R2.mzML" --size 733797294
& $python tools/scripts/benchmark/peptide_double_competition/download_ranges.py "$data/UPS1_50000amol_R1.mzML" "$work/raw/UPS1_50000amol_R1.mzML" --size 767043386
& $python tools/scripts/benchmark/peptide_double_competition/download_ranges.py "$data/UPS1_50000amol_R2.mzML" "$work/raw/UPS1_50000amol_R2.mzML" --size 767355283
Invoke-WebRequest "$data/uniprot_yeast_reviewed_isoforms_ups1_crap.fasta_td.fasta" -OutFile "$work/uniprot_yeast_reviewed_isoforms_ups1_crap.fasta_td.fasta"
& $python tools/scripts/benchmark/peptide_double_competition/prepare_target_fasta.py "$work/uniprot_yeast_reviewed_isoforms_ups1_crap.fasta_td.fasta" "$work/target.fasta"

# Generate Crux's explicit peptide pairing from the searched target database.
& $crux tide-index --overwrite T --output-dir "$work/index_metadata" --peptide-list T --missed-cleavages 1 --enzyme trypsin --seed 1 --mods-spec C+57.02146 --decoy-format shuffle --keep-terminal-aminos NC "$work/target.fasta" "$work/ups48"

# Search the same concatenated target/decoy database; top-match 1 performs the PSM competition.
foreach ($run in 'UPS1_2500amol_R1','UPS1_2500amol_R2','UPS1_50000amol_R1','UPS1_50000amol_R2') {
  & $crux tide-search --output-dir "$work/search/$run" --overwrite T --concat T --top-match 1 --num-threads 8 --precursor-window 25 --precursor-window-type ppm --mzid-output T "$work/raw/$run.mzML" "$work/ups48"
  if ($LASTEXITCODE -ne 0) { throw "Crux search failed for $run" }
  & $python tools/scripts/benchmark/peptide_double_competition/crux_to_idxml.py "$work/search/$run/tide-search.txt" "$work/index_metadata/tide-index.peptides.txt" "$work/converted/$run.idXML" "$work/converted/$run.pairs.tsv"
  if ($LASTEXITCODE -ne 0) { throw "Crux conversion failed for $run" }
}

$env:OPENMS_DATA_PATH = (Resolve-Path '.\share\OpenMS').Path
foreach ($run in 'UPS1_2500amol_R1','UPS1_2500amol_R2','UPS1_50000amol_R1','UPS1_50000amol_R2') {
  $input = "$work/converted/$run.idXML"
  # Legacy peptide-level FDR; defaults are unchanged.
  & $fdr -in $input -out "$work/converted/$run.legacy.idXML" -PSM true -peptide true -protein false -algorithm:add_decoy_peptides true
  # New method. Include decoy winners here so the summary can report the decoy counts at each q-value threshold.
  & $fdr -in $input -out "$work/converted/$run.double.idXML" -PSM true -peptide true -protein false -peptide_method PSM_and_peptide -peptide_pairs "$work/converted/$run.pairs.tsv" -algorithm:add_decoy_peptides true
}
& $python tools/scripts/benchmark/peptide_double_competition/summarize_idxml.py `
  "$work/converted/UPS1_2500amol_R1.legacy.idXML" "$work/converted/UPS1_2500amol_R1.double.idXML" `
  "$work/converted/UPS1_2500amol_R2.legacy.idXML" "$work/converted/UPS1_2500amol_R2.double.idXML" `
  "$work/converted/UPS1_50000amol_R1.legacy.idXML" "$work/converted/UPS1_50000amol_R1.double.idXML" `
  "$work/converted/UPS1_50000amol_R2.legacy.idXML" "$work/converted/UPS1_50000amol_R2.double.idXML"
```

For runtime and peak-memory comparison, build the `FalseDiscoveryRate` target in Release mode and repeat the `FalseDiscoveryRate` commands without `-algorithm:add_decoy_peptides true`; that is the default-output mode measured below. Run each method three times per input on the same machine and report the median `FalseDiscoveryRate took ... (wall)` and `Peak Memory Usage` values from its standard output. The decoy-output run is only for counting decoy sequences at the reported thresholds.

Crux generates 436,658 target/decoy peptide pairs for this FASTA; 77 target peptides have no shuffled decoy and are not pairable. The converter excludes PSMs for those unpaired targets, verifies each remaining decoy PSM against Crux's `original target sequence` and `tide-index.peptides.txt`, and writes the exact pairs to the TSV passed to OpenMS. This prevents a proxy pairing rule from changing the method being evaluated.

## Results

Counts are unique peptide sequences at the indicated peptide q-value thresholds. `UPS+` is the subset of target sequences assigned to one of the 48 human UPS proteins. `Decoy` is the decoy sequence count emitted by `-algorithm:add_decoy_peptides true`; for double competition these are peptide-competition winners. The decoy counts and q-value thresholds are the target-decoy calibration diagnostic, not an independent measurement of empirical false discoveries. Runtime and peak memory are medians across three Release runs per method/input in the normal default-output mode. Measurements used a Windows 10 Pro x64 machine (build 26300) with an Intel Core i7-13700HX (16 cores / 24 logical processors) and the Visual Studio 18 2026 generator.

| Run | Unpaired PSMs excluded | Paired PSMs | q ≤ 1%: target (UPS+) / decoy, legacy → double | q ≤ 5%: target (UPS+) / decoy, legacy → double | q ≤ 10%: target (UPS+) / decoy, legacy → double | Runtime, legacy → double | Peak memory, legacy → double |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2500 amol R1 | 2 | 36,913 | 31 (0) / 1 → 0 (0) / 0 | 147 (0) / 7 → 111 (0) / 4 | 284 (0) / 28 → 281 (0) / 27 | 1.06 s → 2.12 s | 88 MB → 216 MB |
| 2500 amol R2 | 0 | 36,947 | 79 (0) / 1 → 0 (0) / 0 | 102 (0) / 5 → 96 (0) / 3 | 251 (1) / 26 → 245 (1) / 23 | 0.97 s → 1.84 s | 88 MB → 217 MB |
| 50000 amol R1 | 3 | 37,902 | 45 (2) / 1 → 0 (0) / 0 | 161 (16) / 9 → 143 (14) / 6 | 397 (34) / 40 → 381 (30) / 37 | 0.89 s → 1.87 s | 90 MB → 218 MB |
| 50000 amol R2 | 2 | 37,564 | 17 (3) / 1 → 0 (0) / 0 | 67 (9) / 4 → 0 (0) / 0 | 321 (23) / 33 → 321 (23) / 31 | 0.86 s → 1.85 s | 89 MB → 217 MB |

On these four inputs, double competition did not increase target peptide yield at any reported threshold. It returned no target peptide at 1% in any run. At 5%, target yield was lower in every run, including none at 5% for the 50000 amol R2 run. At 10%, the results were closer but remained equal or lower for double competition. UPS-positive recovery followed the same pattern. Across the four runs, double competition took 1.90–2.15× as long and used 2.42–2.47× peak memory. These ratios describe this dataset and implementation; they are not a general performance guarantee.

This benchmark reports the estimator's q-value/decoy-winner behavior and UPS-positive recovery on the same spectra. It does not claim an independent empirical FDP: the sample also contains a yeast digest, so non-UPS target peptides cannot be labeled false merely because they are not in the UPS standard. The observed lack of a yield gain is reported as measured; it is not evidence that the paper's reported average gain is incorrect, since this is one instrument/data-processing setup and four runs from one spike-in series.

## Scripts

- `download_ranges.py`: parallel, size-checked byte-range download using `curl.exe`.
- `prepare_target_fasta.py`: removes the already-present `DECOY_` entries before Crux generates its own shuffled pairs.
- `crux_to_idxml.py`: validates Crux pairing and converts one-hit TSV results to idXML plus the OpenMS pair TSV.
- `summarize_idxml.py`: counts target, decoy, and UPS-positive peptide hits at 1%, 5%, and 10% q-value.
