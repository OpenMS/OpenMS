#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib.sh
source "$SCRIPT_DIR/lib.sh"

DATASET="${1:-/home/sing/Documents/datasets/proteobench-astral-entrapment}"
OUT="${2:-$PWD/opendia-benchmark/astral}"
shift $(( $# >= 2 ? 2 : $# ))

require_executable "$OPEN_DIA"
mapfile -t RAW_FILES < <(find "$DATASET/raw" -maxdepth 1 -type f -name '*.raw' | sort)
(( ${#RAW_FILES[@]} == 3 )) || { echo "ERROR: expected 3 Astral .raw files, found ${#RAW_FILES[@]}" >&2; exit 2; }

if [[ -n "${FASTA:-}" ]]; then
  FASTA="$FASTA"
else
  mapfile -t FASTA_FILES < <(find "$DATASET" -maxdepth 2 -type f -name '*.fasta' | sort)
  (( ${#FASTA_FILES[@]} == 1 )) || { echo "ERROR: set FASTA when the dataset contains zero or multiple FASTA files" >&2; exit 2; }
  FASTA="${FASTA_FILES[0]}"
fi
require_file "$FASTA"
prepare_output "$OUT"
mkdir -p "$OUT/results" "$OUT/intermediates" "$OUT/tmp"
printf '%s\n' "${RAW_FILES[@]}" > "$OUT/meta/raw_files.txt"
sha256sum "${RAW_FILES[@]}" > "$OUT/meta/raw_files.sha256"
sha256sum "$FASTA" > "$OUT/meta/fasta.sha256"

MAX_PEPTIDE_LENGTH="${MAX_PEPTIDE_LENGTH:-$(awk '
  /^>/ { if (length(sequence) > max) max=length(sequence); sequence=""; next }
  { sequence=sequence $0 }
  END { if (length(sequence) > max) max=length(sequence); print max }
' "$FASTA")}"
PREDICTION_THREADS="${PREDICTION_THREADS:-$THREADS}"
(( PREDICTION_THREADS > 8 )) && PREDICTION_THREADS=8
VARIABLE_MOD_ARGS=()
if [[ "${OXIDATION:-true}" == "true" ]]; then
  VARIABLE_MOD_ARGS=(-PeptideQueryParameters:LibraryPrediction:variable_modifications "Oxidation (M)")
fi

command=(
  env "LD_LIBRARY_PATH=$BUILD_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
  "$OPEN_DIA"
  -in "${RAW_FILES[@]}"
  -tr "$FASTA"
  -tr_type fasta
  -workflow:library_mode predicted
  -workflow:working_format sqlite
  -workflow:keep_intermediate_files true
  -workflow:intermediate_dir "$OUT/intermediates"
  -out_dir "$OUT/results"
  -threads "$THREADS"
  -PeptideQueryParameters:LibraryPrediction:enzyme "no cleavage"
  -PeptideQueryParameters:LibraryPrediction:missed_cleavages 0
  -PeptideQueryParameters:LibraryPrediction:min_peptide_length 7
  -PeptideQueryParameters:LibraryPrediction:max_peptide_length "$MAX_PEPTIDE_LENGTH"
  -PeptideQueryParameters:LibraryPrediction:precursor_charges 1 2 3 4 5
  -PeptideQueryParameters:LibraryPrediction:fixed_modifications "Carbamidomethyl (C)"
  "${VARIABLE_MOD_ARGS[@]}"
  -PeptideQueryParameters:LibraryPrediction:max_variable_modifications 1
  -PeptideQueryParameters:LibraryPrediction:clip_nterm_methionine false
  -PeptideQueryParameters:LibraryPrediction:prediction_batch_size 5000
  -PeptideQueryParameters:LibraryPrediction:inference_threads "$PREDICTION_THREADS"
  -PeptideQueryParameters:LibraryPrediction:nce 25
  -PeptideQueryParameters:LibraryPrediction:instrument_index 4
  -PeptideQueryParameters:LibraryPrediction:predict_ccs false
  -PeptideQueryParameters:AssayGenerator:min_transitions 6
  -PeptideQueryParameters:AssayGenerator:max_transitions 6
  -PeptideQueryParameters:AssayGenerator:precursor_lower_mz_limit 400
  -PeptideQueryParameters:AssayGenerator:precursor_upper_mz_limit 1000
  -PeptideQueryParameters:AssayGenerator:product_lower_mz_limit 100
  -PeptideQueryParameters:AssayGenerator:product_upper_mz_limit 1800
  -PeptideQueryParameters:DecoyGenerator:method pseudo-reverse
  -PeptideQueryParameters:DecoyGenerator:aim_decoy_fraction 1
  -TargetedDataExtraction:enable_ipf false
  -TargetedDataExtraction:ion_mobility_mode disabled
  "${OPEN_DIA_SPEED_ARGS[@]}"
  "${OPEN_DIA_ANALYSIS_ARGS[@]}"
  -TargetedDataExtraction:tempDirectory "$OUT/tmp"
  -TargetedDataExtraction:write_debug_files false
  -TargetedDataExtraction:Scoring:stop_report_after_feature 5
  -TargetedDataExtraction:Calibration:auto_irt:enabled true
  -TargetedDataExtraction:Calibration:auto_irt:irt_seed 5489
  -TargetedDataExtraction:Calibration:auto_irt:irt_bins 100
  -TargetedDataExtraction:Calibration:auto_irt:irt_peptides_per_bin 10
  -TargetedDataExtraction:Calibration:auto_irt:linear_top_fraction 0.10
  -TargetedDataExtraction:Calibration:auto_irt:irt_peptides_per_bin_nonlinear 0
  -TargetedDataExtraction:Calibration:linear:outlier_detection iter_residual
  -TargetedDataExtraction:Calibration:RTNormalization:estimateBestPeptides
  -TargetedDataExtraction:Calibration:qc:min_rsq 0.95
  -TargetedDataExtraction:Calibration:qc:min_coverage 0.60
  -TargetedDataExtraction:Calibration:MassIMCorrection:mz_correction_function quadratic_regression_delta_ppm
  -Rescoring:level ms1ms2
  -Rescoring:percolator:subset_max_train 100000
  -Export:results:run true
  -Export:results:format tsv
  "$@"
)
run_measured "$OUT" "${command[@]}"
python3 "$SCRIPT_DIR/metrics_astral.py" "$OUT/results/OpenDIA.results.tsv" \
  --output "$OUT/metrics.json"
