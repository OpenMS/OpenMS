#!/usr/bin/env bash
# Re-run the key diagnosis summaries over all 10 dev runs (raw01/02 seeds 0-4). New output files only (*_all10.txt).
set -euo pipefail
cd /home/singjc/.jcode/scratch/fdr-20261125
L=/home/singjc/.jcode/scratch/perf/lib/search.pqp
py() { ~/.local/bin/uv run -q --no-project --with pandas --with numpy --with scipy python "$@"; }
py summary_table.py "$L" raw0[12]_s[0-4]/*.features.tsv > summary_table_all10.txt
py summary_ov0.py "$L" raw0[12]_s[0-4]/*.features.tsv > summary_ov0_all10.txt
py sim_combo.py "$L" raw0[12]_s[0-4]/*.features.tsv > sim_combo_all10.txt
py diag_overlap2.py "$L" raw0[12]_s[0-4]/*.features.tsv > diag_overlap2_all10.txt
py diag_sweep.py "$L" raw0[12]_s[0-4]/*.features.tsv > diag_sweep_all10.txt
tail -n 30 summary_table_all10.txt summary_ov0_all10.txt sim_combo_all10.txt diag_overlap2_all10.txt diag_sweep_all10.txt
