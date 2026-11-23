"""Summarize preserved library-only audit records without any feature/candidate claims.
Usage: summarize_asymmetry.py AUDIT.json OUTPUT.json
"""
import collections
import hashlib
import json
import sys

source, output = sys.argv[1:]
with open(source) as handle:
    audit = json.load(handle)
pairs = audit['own_decoy_pairs']
deltas = collections.Counter(row['delta_library_siblings'] for row in pairs)
sides = {}
for label in ['target', 'decoy']:
    counts = collections.Counter(row[label + '_library_siblings'] for row in pairs)
    sides[label] = dict(distribution=dict(sorted(counts.items())),
                        mean=sum(k*v for k, v in counts.items()) / max(len(pairs), 1))
for row in pairs:
    assert row['target_candidate_windows'] is None
    assert row['decoy_candidate_windows'] is None
with open(source, 'rb') as handle:
    digest = hashlib.file_digest(handle, 'sha256').hexdigest()
result = dict(source=source, source_sha256=digest, provenance=audit['provenance'],
              candidate_population=audit['candidate_population'], library_precursors=len(audit['precursors']),
              own_decoy_pairs=len(pairs), unpaired_precursors=len(audit['unpaired_precursor_ids']),
              equal=sum(v for k,v in deltas.items() if k == 0),
              target_more=sum(v for k,v in deltas.items() if k > 0),
              decoy_more=sum(v for k,v in deltas.items() if k < 0),
              target_minus_decoy_siblings=dict(sorted(deltas.items())), sides=sides,
              interpretation='Structural opportunities only, not proof of null asymmetry or FDR inflation')
with open(output, 'x') as handle:
    json.dump(result, handle, indent=2)
    handle.write('\n')
print(json.dumps(result, indent=2))
