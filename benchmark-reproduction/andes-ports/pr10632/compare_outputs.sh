#!/bin/bash
# Hash the outputs of two arms. MSstats CSVs are compared as they are; consensusXML and mzTab without the lines
# that carry the OpenMS version/revision, dates, run ids or file paths of the run.
a=$1; b=$2
for f in out_msstats.csv norm_msstats.csv; do
  echo "$f: $(sha256sum < $a/$f | cut -c1-16) $(sha256sum < $b/$f | cut -c1-16)"
done
filt='feature_id|parameter: |version|revision|date|Date|time|UID|id="|unique_id|file_path|location|software\[|Spectra_Filepath|ms_run\[|filename|document_id'
for f in out.consensusXML out.mzTab; do
  echo "$f (filtered): $(grep -Ev "$filt" $a/$f | sha256sum | cut -c1-16) $(grep -Ev "$filt" $b/$f | sha256sum | cut -c1-16)"
done
