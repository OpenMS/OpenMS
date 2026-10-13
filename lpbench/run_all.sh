#!/bin/bash
# Build OpenMS once per LP backend, run the full ctest suite and the LPWrapper benchmark.
# Usage: lpbench/run_all.sh [SOLVER...]   (default: HIGHS COIN GLPK)
set -u
SRC=/home/user/OpenMS
OUT=$SRC/lpbench/results
mkdir -p "$OUT"
SOLVERS=${*:-HIGHS COIN GLPK}
for S in $SOLVERS; do
  B=$SRC/OpenMS-build-${S,,}
  log=$OUT/$S
  echo "=== $S $(date -Is)" >> "$OUT/progress.log"
  if [ ! -f "$B/build.ninja" ]; then
    cmake -S "$SRC" -B "$B" -G Ninja -DCMAKE_BUILD_TYPE=Release -DLP_SOLVER=$S \
      -DCMAKE_PREFIX_PATH="/opt/arrow;/opt/highs" \
      -DENABLE_DOCS=OFF -DWITH_THERMO_RAW=OFF -DWITH_OPENTIMS=OFF -DENABLE_UPDATE_CHECK=OFF \
      -DGIT_TRACKING=OFF -DHAS_XSERVER=OFF -DMT_REQUIRE_OPENMP=ON > "$log.configure.log" 2>&1 \
      || { echo "$S configure FAILED" >> "$OUT/progress.log"; continue; }
  fi
  /usr/bin/time -f "%e s wall" -o "$log.buildtime" cmake --build "$B" -j4 -k 0 > "$log.build.log" 2>&1
  echo "$S build rc=$? $(cat $log.buildtime) $(date -Is)" >> "$OUT/progress.log"
  /usr/bin/time -f "%e s wall" -o "$log.testtime" ctest --test-dir "$B" -j4 --timeout 1500 \
    --output-junit "$log.junit.xml" > "$log.ctest.log" 2>&1
  echo "$S ctest rc=$? $(tail -n 30 $log.ctest.log | grep 'tests passed') $(date -Is)" >> "$OUT/progress.log"
  # benchmark binary against this build's libOpenMS
  g++ -O2 -std=c++20 -fopenmp "$SRC/lpbench/lpbench.cpp" -o "$B/lpbench" \
    -I"$SRC/src/openms/include" -I"$B/src/openms/include" -I"$SRC/src/openms/extern" \
    $(find "$B/src/openms" -maxdepth 3 -type d -name include -printf '-I%p ') \
    -L"$B/lib" -lOpenMS -Wl,-rpath,"$B/lib:/opt/highs/lib:/opt/arrow/lib" > "$log.bench_build.log" 2>&1 \
    && "$B/lpbench" > "$log.bench.csv" 2> "$log.bench.err"
  echo "$S bench rc=$? $(date -Is)" >> "$OUT/progress.log"
done
echo "ALL DONE $(date -Is)" >> "$OUT/progress.log"
