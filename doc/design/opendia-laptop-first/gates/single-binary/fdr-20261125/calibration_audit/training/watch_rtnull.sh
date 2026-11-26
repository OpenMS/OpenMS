#!/usr/bin/env bash
# Follow-through only. No submission or parameter changes.
set -euo pipefail
BASE=/home/singjc/.jcode/scratch/rtnull-20261126-monitor
DOC=/home/singjc/Documents/github/OpenMS-wp00-design-docs/doc/design/opendia-laptop-first/gates/single-binary/fdr-20261125
JOBS=44052409,44052410,44052411,44052412,44052413,44052414,44052415,44052416,44052417,44052418
mkdir -p "$BASE"
remote() {
 ssh -n -o BatchMode=yes -o ControlMaster=no -o ControlPath=none codon "become pst_prd bash -c '$*'"
}
previous=
while true; do
 if state=$(remote "export PATH=/ebi/slurm/codon/bin:\$PATH; sacct -n -X -j $JOBS --format=JobIDRaw,State,ExitCode -P"); then
  printf '%s\n' "$state" > "$BASE/status.latest"
  if [[ $state != "$previous" ]]; then
   printf 'JCODE_PROGRESS {"message":"RTNULL Slurm state changed"}\n'
   date -u +%Y-%m-%dT%H:%M:%SZ
   printf '%s\n' "$state"
   previous=$state
  fi
  if printf '%s\n' "$state" | awk -F '|' 'NF>=3 {n++; if($2 !~ /^(COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|BOOT_FAIL|DEADLINE)/) bad=1} END {exit !(n==10 && !bad)}'; then break; fi
 fi
 sleep 60
done
FETCH=$BASE/fetch
mkdir "$FETCH"
remote 'cd /nfs/research/juan/singj/opendia/rtnull-20261126; find runs logs -type f \( -name "*.precursors.tsv" -o -name "*.log" -o -name "*.exit" -o -name "*.out" \) -print0 | tar --null -czf - --files-from=-' > "$BASE/outputs.tgz"
tar -xzf "$BASE/outputs.tgz" -C "$FETCH"
remote 'cd /nfs/research/juan/singj/opendia/rtnull-20261126; find runs logs -type f \( -name "*.precursors.tsv" -o -name "*.log" -o -name "*.exit" -o -name "*.out" \) -print0 | sort -z | xargs -0 sha256sum' > "$BASE/fetched.sha256"
(cd "$FETCH" && sha256sum -c "$BASE/fetched.sha256")
bash "$DOC/rtnull/rtnull_eval.sh" "$FETCH" > "$BASE/EVAL_RTNULL.txt" 2>&1
cat "$BASE/EVAL_RTNULL.txt"
