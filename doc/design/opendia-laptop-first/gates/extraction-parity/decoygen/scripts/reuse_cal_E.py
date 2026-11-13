"""PLAN_DECOYGEN addendum (arm E): the PAIRLIB run calibration for the C++-built library. cal.json is copied unchanged;
target anchors are remapped from PAIRLIB precursor IDs to arm E's by (modified sequence, charge), because the C++ PQP writer
renumbers TRAML IDs. Anchors whose target is not in arm E are dropped and counted. Usage: reuse_cal_E.py"""
import csv, os, shutil, sqlite3
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
P = B + '/pairlib-20261113'; H = B + '/decoygen-20261113/armE/home'
q = """select p.ID, pe.MODIFIED_SEQUENCE || '/' || p.CHARGE, p.DECOY from PRECURSOR p
  join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID"""
old = {str(i): k for i, k, d in sqlite3.connect(f'file:{P}/lib/library.pqp?mode=ro', uri=True).execute(q) if d == 0}
new = {k: str(i) for i, k, d in sqlite3.connect(f'file:{H}/lib/library.pqp?mode=ro', uri=True).execute(q) if d == 0}
for raw in ('raw01', 'raw02'):
    os.makedirs(f'{H}/{raw}/cal', exist_ok=True)
    shutil.copyfile(f'{P}/{raw}/cal/cal.json', f'{H}/{raw}/cal/cal.json')
    kept = dropped = 0
    with open(f'{P}/{raw}/cal/anchors.tsv', newline='') as fi, open(f'{H}/{raw}/cal/anchors.tsv', 'w', newline='') as fo:
        r = csv.DictReader(fi, delimiter='\t'); w = csv.DictWriter(fo, r.fieldnames, delimiter='\t', lineterminator='\n'); w.writeheader()
        for row in r:
            k = old.get(row['precursor_id'])
            if k is None or k not in new:
                dropped += 1; continue
            row['precursor_id'] = new[k]; w.writerow(row); kept += 1
    print(f'arm E {raw}: cal.json copied, anchors remapped {kept}, dropped {dropped}')
