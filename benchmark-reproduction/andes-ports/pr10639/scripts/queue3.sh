cd /home/user/bench
D=frozen/z/prose-andes-reproduction/broad_benchmark/data
for r in 2 3; do for sl in 10 40 150; do scripts/run_andes.sh andes_runs/phos2_E_sl${sl}_r$r $sl db/hye_E.fasta $D/hfx_A2.mzML andes_runs/mods_phos2.txt 2; done; done
python3 -I scripts/bench.py experiments.json phos2_E_chunk 3
echo QUEUE3_DONE
