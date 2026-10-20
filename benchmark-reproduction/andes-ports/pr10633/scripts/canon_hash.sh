#!/bin/bash
# sha256 of an output with run-specific noise removed, for the base/PR identity check:
#  consensusXML: unique ids, feature_id references, dates and the run's own output paths/parameter list
#  mzTab: metadata (MTD) lines, random unique feature ids (>= 10 digits), rows sorted
#  MSstats CSV: as is
f=$1
case $f in
  *.consensusXML) sed -E 's/ id="[a-zA-Z_]*[0-9]+"//g; s/ (date|document_id|completion_time)="[^"]*"//g; s/unique_id="[0-9]+"//g; s/name="feature_id" value="[0-9]+"/name="feature_id"/; /name="parameter: /d' "$f" | sha256sum | cut -d" " -f1 ;;
  *.mzTab) grep -v "^MTD" "$f" | sed -E 's/\b[0-9]{10,}\b/ID/g' | sort | sha256sum | cut -d" " -f1 ;;
  *) sha256sum "$f" | cut -d" " -f1 ;;
esac
