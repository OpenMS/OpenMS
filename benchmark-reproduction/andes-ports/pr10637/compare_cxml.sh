#!/bin/bash
# Compare two ProteomicsLFQ consensusXML(.gz) files up to what differs between any two processes: unique ids and
# feature_id values (process-random), map/ref/protein_refs ids, completion time and the 'parameter: ...' provenance
# lines (output paths, and the new PipEcho:shape_features entry). Prints the number of remaining differing lines.
n() { zcat -f "$1" | sed -E 's/ id="[^"]*"//g; s/(identification_run_ref|protein_refs|unique_id|completion_time)="[^"]*"//g; s/ (map|ref)="[^"]*"//g; s/name="feature_id" value="[^"]*"/name="feature_id"/' | grep -v 'name="parameter: '; }
diff <(n "$1") <(n "$2") | grep -c '^<'
