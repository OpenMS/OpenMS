#!/bin/bash
# #10693 step 2: extraction evidence (FFID quality, M0 width ratio, mass traces) as extra SVM predictors
# (PIPECHO_EXTRA_PREDICTORS), prototype v3 = develop 5edf8394 + predictors_v3.patch. Four jobs in parallel.
S=/home/user/bench/scripts; D=/home/user/data; R=/home/user/runs; P=/opt/openms-proto3
export OPENMS_DATA_PATH=/home/user/OpenMS/share/OpenMS PIPECHO_EXTRA_PREDICTORS=1
H="$D/db/HYE_td.fasta $D/ids_hye"; MX="-targeted_only false -pip_echo true"
mkdir -p $R/pred
job() {
  local name=$1; shift
  [ -s $R/pred/$name/out_msstats.csv ] && { echo "skip $name"; return; }
  case $name in
    ctrl_x_q05)   env PIPECHO_PROBE_DECOYS=120 PIPECHO_PROBE_FIXED=1 PIPECHO_SCORE_DUMP=$R/pred/$name.scores.tsv $S/plfq.sh $P $H $R/seed/hye/feat_mx $R/pred/$name $MX -PipEcho:max_training_points 15000 ;;
    ctrl_x_q01)   env PIPECHO_PROBE_DECOYS=120 PIPECHO_PROBE_FIXED=1 $S/plfq.sh $P $H $R/seed/hye/feat_mx $R/pred/$name $MX -PipEcho:max_training_points 15000 -PipEcho:fdr 0.01 ;;
    proto_x_q05)  env PLFQ_EXTRA_SEEDS_DIR=$R/proto/seeds_hye PIPECHO_PROBE_DECOYS=120 PIPECHO_PROBE_FIXED=1 PIPECHO_SCORE_DUMP=$R/pred/$name.scores.tsv $S/plfq.sh $P $H $R/proto/hye/feat_proto $R/pred/$name $MX -PipEcho:max_training_points 15000 ;;
    proto_x_q01)  env PLFQ_EXTRA_SEEDS_DIR=$R/proto/seeds_hye PIPECHO_PROBE_DECOYS=120 PIPECHO_PROBE_FIXED=1 $S/plfq.sh $P $H $R/proto/hye/feat_proto $R/pred/$name $MX -PipEcho:max_training_points 15000 -PipEcho:fdr 0.01 ;;
    rand_x)       env PIPECHO_PSEUDO=15 PIPECHO_PSEUDO_OUT=$R/pred/rand_x.tsv PIPECHO_SCORE_DUMP=$R/pred/$name.scores.tsv $S/plfq.sh $P $H $R/seed/hye/feat_mx $R/pred/$name $MX ;;
    rand_x_absent) env PIPECHO_PSEUDO=15 PIPECHO_PSEUDO_MODE=absent PIPECHO_PSEUDO_OUT=$R/pred/rand_x_absent.tsv PIPECHO_SCORE_DUMP=$R/pred/$name.scores.tsv $S/plfq.sh $P $H $R/seed/hye/feat_mx $R/pred/$name $MX ;;
    probe_x)      env PIPECHO_PSEUDO=15 PIPECHO_PROBE_DECOYS=120 PIPECHO_PSEUDO_OUT=$R/pred/probe_x.tsv $S/plfq.sh $P $H $R/seed/hye/feat_mx $R/pred/$name $MX ;;
    probe_x_absent) env PIPECHO_PSEUDO=15 PIPECHO_PSEUDO_MODE=absent PIPECHO_PROBE_DECOYS=120 PIPECHO_PSEUDO_OUT=$R/pred/probe_x_absent.tsv $S/plfq.sh $P $H $R/seed/hye/feat_mx $R/pred/$name $MX ;;
  esac && echo "done $name" || echo "FAILED $name"
}
export -f job; export S D R P H MX
printf '%s\n' proto_x_q05 ctrl_x_q05 rand_x rand_x_absent proto_x_q01 ctrl_x_q01 probe_x probe_x_absent | xargs -P4 -I{} bash -c 'job {}'
echo PRED_DONE
