"""PLAN_DECOYGEN amendment 1: give each arm the PAIRLIB run calibration (as k11, astral-commatched-decoys-20261109).
cal.json is copied unchanged: it holds the RT map, m/z and peak-boundary models and counts only (no precursor IDs; the
target library RTs are identical across arms). anchors.tsv (target anchors only) is remapped from PAIRLIB precursor IDs to
the arm's IDs by TRAML_ID; every anchor must map. Usage: reuse_cal.py <arm>"""
import csv, os, shutil, sqlite3, sys
ARM = sys.argv[1]
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
P = B + '/pairlib-20261113'; H = f'{B}/decoygen-20261113/arm{ARM}/home'
q = 'select ID, TRAML_ID, DECOY from PRECURSOR'
old = {str(i): (t, d) for i, t, d in sqlite3.connect(f'file:{P}/lib/library.pqp?mode=ro', uri=True).execute(q)}
new = {t: str(i) for i, t, d in sqlite3.connect(f'file:{H}/lib/library.pqp?mode=ro', uri=True).execute(q) if d == 0}
for raw in ('raw01', 'raw02'):
    os.makedirs(f'{H}/{raw}/cal', exist_ok=True)
    shutil.copyfile(f'{P}/{raw}/cal/cal.json', f'{H}/{raw}/cal/cal.json')
    kept = 0
    with open(f'{P}/{raw}/cal/anchors.tsv', newline='') as fi, open(f'{H}/{raw}/cal/anchors.tsv', 'w', newline='') as fo:
        r = csv.DictReader(fi, delimiter='\t'); w = csv.DictWriter(fo, r.fieldnames, delimiter='\t', lineterminator='\n'); w.writeheader()
        for row in r:
            t, d = old[row['precursor_id']]
            assert d == 0 and t in new, (raw, row['precursor_id'])
            row['precursor_id'] = new[t]; w.writerow(row); kept += 1
    print(f'arm {ARM} {raw}: cal.json copied, anchors remapped {kept}')
