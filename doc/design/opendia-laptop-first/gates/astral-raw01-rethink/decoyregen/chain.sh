#!/bin/bash
# After extraction: decode, compute features and score each decoy family, per run.
cd "$(dirname "$0")"; PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
export NUMBA_NUM_THREADS=12
for RUN in 02 01; do
  while [ ! -f extract$RUN.EXIT ]; do sleep 30; done
  [ "$(cat extract$RUN.EXIT)" = "0" ] || { echo "extract $RUN failed" > chain$RUN.EXIT; continue; }
  OPENDIA_RUN=$RUN $PY build_xic2.py > build$RUN.log 2>&1 && \
  OPENDIA_RUN=$RUN $PY feats2.py > feats$RUN.log 2>&1 && \
  OPENDIA_RUN=R$RUN $PY score2.py 2>/dev/null > score$RUN.txt
  echo $? > chain$RUN.EXIT
done
