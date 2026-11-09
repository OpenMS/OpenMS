"""Remap OpenDIACalibrate anchors (targets only) from the old library's precursor IDs to the new library's, by TRAML_ID.
Usage: remap_anchors.py OLD.pqp NEW.pqp anchors.tsv out.tsv"""
import sys, csv, sqlite3
old, new, src, dst = sys.argv[1:5]
q = 'select ID, TRAML_ID, DECOY from PRECURSOR'
o = {str(i): (t, d) for i, t, d in sqlite3.connect(f'file:{old}?mode=ro', uri=True).execute(q)}
n = {t: str(i) for i, t, d in sqlite3.connect(f'file:{new}?mode=ro', uri=True).execute(q) if d == 0}
kept = miss = 0
with open(src, newline='') as fi, open(dst, 'w', newline='') as fo:
    r = csv.DictReader(fi, delimiter='\t'); w = csv.DictWriter(fo, r.fieldnames, delimiter='\t', lineterminator='\n'); w.writeheader()
    for row in r:
        t, d = o[row['precursor_id']]
        if d != 0 or t not in n:
            miss += 1; continue
        row['precursor_id'] = n[t]; w.writerow(row); kept += 1
print(f'anchors remapped {kept}, dropped {miss}')
if miss: sys.exit(1)
