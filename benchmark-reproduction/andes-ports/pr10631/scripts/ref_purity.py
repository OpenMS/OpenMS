#!/usr/bin/env python3
"""Reference precursor purity, a Python transcription of andes crates/quant/src/purity.rs (bigbio/andes@40f3ee1):
window [target - lower, target + upper] centred on the isolation target, isotope ladder anchored on the selected
ion m/z, only peaks inside the window count, 10 ppm isotope tolerance, survey scan = the MS2's spectrumRef.
Used to check OpenMS' purity values for MS2 precursors whose selected ion lies outside the isolation window.
usage: ref_purity.py <mzML> <scans.tsv> <out.tsv>  (computes it for every MS2)"""
import sys, csv, base64, zlib, struct
import numpy as np
import xml.etree.ElementTree as ET

ISO = 1.00335483507

def decode(b64, is64, comp):
    raw = base64.b64decode(b64)
    if comp: raw = zlib.decompress(raw)
    return np.frombuffer(raw, dtype='<f8' if is64 else '<f4').astype(float)

def read_ms1(path, wanted):
    out = {}
    for ev, el in ET.iterparse(path, events=('end',)):
        if el.tag.rsplit('}', 1)[-1] != 'spectrum': continue
        sid = el.get('id')
        if sid in wanted:
            arrs = {}
            for bda in el.iter():
                if bda.tag.rsplit('}', 1)[-1] != 'binaryDataArray': continue
                accs = {c.get('accession') for c in bda.iter() if c.tag.rsplit('}', 1)[-1] == 'cvParam'}
                b = [x for x in bda.iter() if x.tag.rsplit('}', 1)[-1] == 'binary'][0].text or ''
                v = decode(b, 'MS:1000523' in accs, 'MS:1000574' in accs)
                arrs['mz' if 'MS:1000514' in accs else 'int' if 'MS:1000515' in accs else 'x'] = v
            o = np.argsort(arrs['mz']); out[sid] = (arrs['mz'][o], arrs['int'][o])
        el.clear()
    return out

def best_within(mz, it, x, tol):
    lo, hi = np.searchsorted(mz, x - tol), np.searchsorted(mz, x + tol, 'right')
    if hi <= lo: return None
    j = lo + int(np.argmax(it[lo:hi])); return mz[j], it[j]

def purity(mz, it, prec, target, z, lo_off, hi_off, ppm=10.0):
    c = target if target else prec
    lo, hi = c - lo_off, c + hi_off
    if not hi > lo or len(mz) == 0: return None
    a, b = np.searchsorted(mz, lo), np.searchsorted(mz, hi, 'right')
    total = it[a:b].sum()
    if total <= 0: return None
    sp = ISO / max(z, 1); tol = lambda m: max(m * ppm * 1e-6, 0.001); inw = lambda m: lo <= m <= hi
    anc = best_within(mz, it, prec, tol(prec))
    if anc is None: return 0.0
    env = anc[1] if inw(anc[0]) else 0.0
    for d, edge in ((-1, lo), (1, hi)):
        exp = anc[0] + d * sp; inside = steps = 0
        while d * (edge - exp) >= 0 and inside < 10 and steps < 64:
            was = inw(exp); h = best_within(mz, it, exp, tol(exp))
            if h is not None:
                if inw(h[0]): env += h[1]
                exp = h[0] + d * sp
            else:
                exp += d * sp
            inside += int(was); steps += 1
    return min(env / total, 1.0)

def main(mzml, scans, out):
    rows = list(csv.DictReader(open(scans), delimiter='\t'))
    ms2 = [r for r in rows if r['ms_level'] == '2' and r['selected_mz'] and r['iso_lo']]
    ms1 = read_ms1(mzml, {r['spectrum_ref'] for r in ms2})
    w = csv.writer(open(out, 'w', newline=''), delimiter='\t'); w.writerow(['ms2_id', 'ref_purity', 'ref_purity_selected_centre'])
    for r in ms2:
        if r['spectrum_ref'] not in ms1: continue
        mz, it = ms1[r['spectrum_ref']]
        z = int(r['charge'] or 0); s = float(r['selected_mz']); t = float(r['iso_target'] or 0) or None
        p1 = purity(mz, it, s, t, z, float(r['iso_lo']), float(r['iso_hi']))
        p0 = purity(mz, it, s, None, z, float(r['iso_lo']), float(r['iso_hi']))
        w.writerow([r['native_id'], '' if p1 is None else round(p1, 6), '' if p0 is None else round(p0, 6)])

if __name__ == '__main__':
    main(*sys.argv[1:4])
