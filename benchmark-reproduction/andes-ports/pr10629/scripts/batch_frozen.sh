# frozen-database arms: Velos natural entrapment, shuffled r = 1 entrapment (frozen make_entrapment.py) for HF-X/Astral/TMTpro
S=$1; B=$2; E=$S/work/broad_benchmark/entrapment
cd $S
for d in velos_125_R1 velos_5000_R2 velos_25000_R3; do python3 -I scripts/run.py $B $d $E/hy_natural_td.fasta runs/$d/natural; done
for d in hfx_A2 hfx_B1 hfx_B3 astral_A2 astral_B1 astral_B3; do python3 -I scripts/run.py $B $d $E/hye_entrap_td.fasta runs/$d/frozenE; done
for d in eclipse_tmtpro_10855 eclipse_tmtpro_10858 eclipse_tmtpro_10863; do python3 -I scripts/run.py $B $d $E/human_entrap_td.fasta runs/$d/frozenE; done
echo BATCH_DONE
