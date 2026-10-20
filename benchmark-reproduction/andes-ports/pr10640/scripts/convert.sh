#!/bin/bash
# convert finished raw downloads to MS2-only mzML (marked <name>.mzML.ready) and delete the raw files
export LD_LIBRARY_PATH=/opt/openms-pr/lib OPENMS_TOOL_REGISTRY_PATH=/opt/openms-pr/registry
B=/opt/openms-pr/bin
D=/home/user/data/pxd022787; R=/home/user/bench/results; L=/home/user/bench/logs
while true; do
  while [ $(ls $D/mzml/*.ready 2>/dev/null | wc -l) -ge 20 ]; do sleep 20; done
  for done in $D/raw/*.raw.done; do
    [ -e "$done" ] || continue
    raw=${done%.done}; [ -e $raw ] || continue
    name=$(basename ${raw%.raw})
    $B/FileConverter -in $raw -out $D/mzml/tmp_$name.mzML -process_lowmemory -no_progress > $L/conv_$name.log 2>&1 &&
    $B/FileFilter -in $D/mzml/tmp_$name.mzML -out $D/mzml/$name.mzML -peak_options:level 2 -peak_options:zlib_compression true -no_progress >> $L/conv_$name.log 2>&1 &&
    rm -f $raw $D/mzml/tmp_$name.mzML && touch $D/mzml/$name.mzML.ready || echo "conversion failed: $name"
  done
  grep -q FETCH_DONE /home/user/bench/fetch.log && [ -z "$(ls $D/raw/*.raw 2>/dev/null)" ] && { echo CONVERT_DONE; exit 0; }
  sleep 20
done
