#!/bin/bash
#SBATCH --job-name=afull-evaluate
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -euo pipefail
export PATH=/ebi/slurm/codon/bin:$PATH
cd /nfs/research/juan/singj/opendia/astral-full-20261123
mkdir -p results
sacct -j "$(paste -sd, job_ids.txt)" --parsable2 --units=K --format=JobID,JobName,State,ExitCode,Submit,Start,End,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS,NodeList > results/slurm_accounting.tsv
"$HOME/.local/bin/uv" run -q --no-project --with pandas python eval_full.py library/search6.pqp runs results > results/evaluation.log
if test -s results/calibration_curve.tsv; then
 "$HOME/.local/bin/uv" run -q --no-project --with pandas --with matplotlib python plot_full.py results
fi
printf 'raw\tseed\tpath\n' > results/quant_manifest.tsv
complete=0
for raw in 01 02 03; do
 for seed in 0 1 2 3 4; do
  directory="runs/raw${raw}_s${seed}"
  file="$directory/LFQ_Astral_DIA_15min_50ng_Human_${raw}.precursors.tsv"
  if test -f "$directory/exit_code.txt" && test "$(cat "$directory/exit_code.txt")" = 0 && test -f "$file"; then
   printf '%s\t%s\t../%s\n' "$raw" "$seed" "$file" >> results/quant_manifest.tsv
   complete=$((complete+1))
  fi
 done
done
if test "$complete" = 15; then
 for q in 0.005 0.01 0.02; do
  for population in all human; do
   "$HOME/.local/bin/uv" run -q --no-project python eval_quant.py results/quant_manifest.tsv --q-cutoff "$q" --library library/search6.pqp --population "$population" > "results/quant_${population}_q${q}.json"
  done
 done
else
 printf 'Quantification incomplete: %s/15 successful runs. No substitute seeds.\n' "$complete" > results/QUANT_INCOMPLETE.txt
fi
find library runs -type f -printf '%s\t%p\n' > results/artifact_bytes.tsv
printf 'evaluation_complete\n' > results/evaluation.complete
