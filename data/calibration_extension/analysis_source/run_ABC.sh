#!/usr/bin/env bash
# Waves A, B, C in protocol order; only launched after analyze_calext.py anchor reports ANCHOR_PASS=true.
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1
PY=/c/Users/admin/miniforge3/python.exe
grep -q '"ANCHOR_PASS": true' results/analysis_anchor.json || { echo "anchor not passed; refusing"; exit 2; }
for w in A B C; do
  $PY calext_studies.py $w 24 > logs/$w.log 2> logs/$w.err || echo "wave $w exited nonzero" >> logs/$w.log
  # resume pass for any failed block
  $PY calext_studies.py $w 24 >> logs/$w.log 2>> logs/$w.err
done
echo ALL_WAVES_DONE >> logs/C.log
