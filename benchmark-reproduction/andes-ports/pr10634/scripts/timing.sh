#!/bin/bash
# Interleaved timing: PR build, option off vs on, 3 repetitions each, 4 threads; search only (no Percolator),
# then one Percolator 3.09 run (seed 42) per arm and repetition on the PIN of that repetition.
R=/home/user/bench; D=$R/frozen/broad_benchmark/data
run() { # name params mzml fasta
  for rep in 1 2 3; do for arm in off on; do
    X=(); [ $arm = on ] && X=(--set Search:report:collapse_positional_isomers=true)
    o=$R/runs/timing/$1/$arm.$rep; rm -rf $o
    python3 -I $R/scripts/run_arm.py --bin /opt/openms-pr --params $2 --mzml $3 --fasta $4 --out $o --seeds "${X[@]}" > /dev/null
    s=$(date +%s.%N); $R/percolator -Y -U --seed 42 --results-psms $o/t.tsv --decoy-results-psms $o/d.tsv $o/input.pin > $o/perc.log 2>&1
    e=$(date +%s.%N)
    python3 -c "import json;m=json.load(open('$o/metrics.json'));print('$1','$arm',$rep,m['wall_s'],m['max_rss_kb'],round($e-$s,2),sep='\t')"
    rm -f $o/input.pin $o/t.tsv $o/d.tsv
  done; done
}
echo -e "dataset\tarm\trep\tsearch_wall_s\tsearch_max_rss_kb\tpercolator_wall_s"
run hfx_A2 $R/cfg/hfx_A2.json $D/hfx_A2.mzML $D/hye_td.fasta
run phospho $R/cfg/phospho.json $R/phospho/control2_phospho.mzML $R/phospho_db/mouse_entrap_td.fasta
