"""Own b/y ion-series features at each candidate apex (h11). Own implementation; nothing from OSW's scoring code.
For every candidate (precursor, apex RT) in xicR{RUN}, take the MS2 spectra of its isolation window whose RT is the
nearest to the apex and its two neighbours (apex +- 1 cycle). In each, match all theoretical b and y ions (fragment charge
1 and 2, >= 150.06 Th, <= 2000 Th, carbamidomethyl C, oxidised M) within +-17 ppm. Features, averaged over the 3 spectra:
  by_nb, by_ny: fraction of theoretical b / y ions matched (charge 1 or 2);
  by_run: longest run of consecutive cleavage sites covered by b or y (fraction of length - 1);
  by_int: matched-ion intensity / total spectrum intensity (TIC share).
Spectra are streamed from mzML once; candidates are bucketed by isolation window.
Output: $JCODE_SCRATCH_DIR/xicR{RUN}/by_feats.npy, aligned with oswlike_X.npy rows.
v2 also writes by_counts.npy: integer counts of matched b and y ions, fragment charge 1 only, nearest peak within 10 ppm
and intensity > 300, in the apex spectrum only (count variant; by_vs_osw02.txt motivates it).
Usage: OPENDIA_RUN=R02."""
import os, re, sys, base64, sqlite3, time, numpy as np, pandas as pd, numba as nb
RUN = os.environ['OPENDIA_RUN']; S = os.environ['JCODE_SCRATCH_DIR'] + '/xic' + RUN
HERE = os.path.dirname(os.path.abspath(__file__)); MZML = os.environ['JCODE_SCRATCH_DIR'] + f'/mzml/raw{RUN[-2:]}.mzML'
m = pd.read_parquet(S + '/oswlike_meta.parquet'); apex = m.apex.values.astype(np.float64); pid = m.precursor_id.values
c = sqlite3.connect('file:' + HERE + '/merged.pqp?mode=ro', uri=True)
pp = pd.read_sql('select p.ID pid, p.PRECURSOR_MZ mz, pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING mm on mm.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=mm.PEPTIDE_ID', c).drop_duplicates('pid').set_index('pid')
AA = dict(G=57.02146, A=71.03711, S=87.03203, P=97.05276, V=99.06841, T=101.04768, C=103.00919 + 57.02146, L=113.08406,
          I=113.08406, N=114.04293, D=115.02694, Q=128.05858, K=128.09496, E=129.04259, M=131.04049, H=137.05891,
          F=147.06841, R=156.10111, Y=163.06333, W=186.07931)
PROTON, H2O = 1.007276, 18.010565
def residues(ms):
    out = []
    for a, mod in re.findall(r'([A-Z])(\(UniMod:\d+\))?', ms):
        v = AA[a]
        if mod == '(UniMod:35)': v += 15.994915
        out.append(v)
    return np.array(out)
# per unique precursor: theoretical ions, with series id (0 b, 1 y) and cleavage index
up = np.unique(pid); ION_OFF = np.zeros(len(up) + 1, np.int64); ions = []; LEN = np.zeros(len(up), np.int64); pmz = np.zeros(len(up))
for k, p in enumerate(up):
    r = residues(pp.ms.loc[p]); n = len(r); LEN[k] = n; pmz[k] = pp.mz.loc[p]
    b = np.cumsum(r)[:-1]; y = np.cumsum(r[::-1])[:-1] + H2O; idx = np.arange(1, n)
    rows = []
    for z in (1, 2):
        rows.append(np.c_[(b + z * PROTON) / z, np.zeros(n - 1), idx, np.full(n - 1, z)]); rows.append(np.c_[(y + z * PROTON) / z, np.ones(n - 1), n - idx, np.full(n - 1, z)])
    a = np.vstack(rows); a = a[(a[:, 0] >= 150.06) & (a[:, 0] <= 2000)]
    ions.append(a); ION_OFF[k + 1] = ION_OFF[k] + len(a)
ION = np.vstack(ions); kof = np.searchsorted(up, pid)
print(f'run {RUN}: candidates {len(pid)}, precursors {len(up)}, theoretical ions {len(ION)}', flush=True)
# ---- pass 1: MS2 spectrum index (RT, isolation target) by streaming headers only
hdr = re.compile(rb'name="ms level" value="(\d)".*?name="scan start time" value="([0-9.eE+-]+)".*?(?:isolation window target m/z" value="([0-9.]+)")?', re.S)
t0 = time.time(); spec_rt, spec_tg, spec_pos = [], [], []
with open(MZML, 'rb') as fh:
    buf = b''; pos0 = 0
    while True:
        chunk = fh.read(1 << 26)
        if not chunk: break
        buf += chunk
        last = buf.rfind(b'</spectrum>')
        if last < 0: continue
        block = buf[:last + 11]
        for sm in re.finditer(rb'<spectrum [^>]*>(.*?)</spectrum>', block, re.S):
            s = sm.group(1)
            if b'name="ms level" value="2"' not in s[:2000]: continue
            rt = float(re.search(rb'scan start time" value="([0-9.eE+-]+)"', s).group(1))
            tg = float(re.search(rb'isolation window target m/z" value="([0-9.eE+-]+)"', s).group(1))
            spec_rt.append(rt); spec_tg.append(tg)
        buf = buf[last + 11:]
spec_rt = np.array(spec_rt); spec_tg = np.array(spec_tg)
print(f'MS2 spectra {len(spec_rt)} in {time.time() - t0:.0f} s', flush=True)
# windows: unique isolation targets (rounded), and per window the sorted spectrum ordinals
tg_r = np.round(spec_tg, 2); wins = np.unique(tg_r); wof = np.searchsorted(wins, tg_r)
# candidate -> window: the window whose target is nearest the precursor m/z (DIA windows +-1 Th)
cw = np.searchsorted(wins, pmz[kof]); cw = np.clip(cw, 1, len(wins) - 1); cw = np.where(np.abs(wins[cw - 1] - pmz[kof]) < np.abs(wins[cw] - pmz[kof]), cw - 1, cw)
# per window: spectrum ordinals sorted by RT; per candidate: ordinals of nearest and +-1
need = {}  # spectrum ordinal -> list of (candidate, slot)
for w in range(len(wins)):
    so = np.where(wof == w)[0]; so = so[np.argsort(spec_rt[so])]; cand = np.where(cw == w)[0]
    if len(cand) == 0 or len(so) == 0: continue
    j = np.clip(np.searchsorted(spec_rt[so], apex[cand]), 1, len(so) - 1)
    j = np.where(np.abs(spec_rt[so][j - 1] - apex[cand]) < np.abs(spec_rt[so][j] - apex[cand]), j - 1, j)
    for d in (-1, 0, 1):
        jj = np.clip(j + d, 0, len(so) - 1)
        for ci, s in zip(cand, so[jj]): need.setdefault(int(s), []).append((int(ci), d == 0))
print(f'spectra needed {len(need)}', flush=True)

@nb.njit(cache=True)
def score_spec(mz, inten, cands, kof, ION, ION_OFF, LEN, acc, centre, cnt):
    tic = inten.sum() + 1e-9
    for q in range(len(cands)):
        ci = cands[q]
        k = kof[ci]; a, b = ION_OFF[k], ION_OFF[k + 1]; n = LEN[k]
        covered = np.zeros(n + 1, np.bool_); nb_ = 0; ny_ = 0; tb = 0; ty = 0; mi = 0.0; cb = 0; cy = 0
        for t in range(a, b):
            f = ION[t, 0]; tol = 17e-6 * f
            j = np.searchsorted(mz, f - tol); hit = False; best = 0.0; bestppm = 1e9; bestint = 0.0
            while j < len(mz) and mz[j] <= f + tol:
                hit = True; best = max(best, inten[j])
                pp = abs(mz[j] - f) / f * 1e6
                if pp < bestppm: bestppm = pp; bestint = inten[j]
                j += 1
            if ION[t, 1] == 0: tb += 1
            else: ty += 1
            if hit:
                if ION[t, 1] == 0: nb_ += 1
                else: ny_ += 1
                covered[int(ION[t, 2])] = True; mi += best
                # count variant: fragment charge 1 only (ION[t, 3]), nearest peak within 10 ppm, intensity floor
                if centre[q] and ION[t, 3] == 1 and bestppm <= 10.0 and bestint > 300.0:
                    if ION[t, 1] == 0: cb += 1
                    else: cy += 1
        run = 0; br = 0
        for i in range(1, n):
            if covered[i]: run += 1; br = max(br, run)
            else: run = 0
        acc[ci, 0] += nb_ / max(tb, 1); acc[ci, 1] += ny_ / max(ty, 1); acc[ci, 2] += br / max(n - 1, 1); acc[ci, 3] += mi / tic; acc[ci, 4] += 1
        if centre[q]: cnt[ci, 0] = cb; cnt[ci, 1] = cy
acc = np.zeros((len(pid), 5)); cnt = np.zeros((len(pid), 2)); t0 = time.time(); ordn = -1
def decode(s, tag):
    blk = re.search(rb'name="' + tag + rb' array".*?name="(64|32)-bit float".*?<binary>([^<]*)</binary>', s, re.S)
    return np.frombuffer(base64.b64decode(blk.group(2)), dtype=np.float64 if blk.group(1) == b'64' else np.float32).astype(np.float64)
with open(MZML, 'rb') as fh:
    buf = b''
    while True:
        chunk = fh.read(1 << 26)
        if not chunk: break
        buf += chunk
        last = buf.rfind(b'</spectrum>')
        if last < 0: continue
        for sm in re.finditer(rb'<spectrum [^>]*>(.*?)</spectrum>', buf[:last + 11], re.S):
            s = sm.group(1)
            if b'name="ms level" value="2"' not in s[:2000]: continue
            ordn += 1
            cl = need.get(ordn)
            if cl is None: continue
            mz = decode(s, b'm/z'); it = decode(s, b'intensity'); o = np.argsort(mz)
            score_spec(mz[o], it[o], np.array([x[0] for x in cl], np.int64), kof, ION, ION_OFF, LEN, acc, np.array([x[1] for x in cl], np.bool_), cnt)
        buf = buf[last + 11:]
F = acc[:, :4] / np.maximum(acc[:, 4:5], 1)
np.save(S + '/by_feats.npy', F.astype(np.float32)); np.save(S + '/by_counts.npy', cnt.astype(np.float32))
dec = m.decoy.values.astype(bool); right = m.right.values
print(f'done in {time.time() - t0:.0f} s | mean by_ny: right {F[right, 1].mean():.3f}, other targets {F[(~dec) & (~right), 1].mean():.3f}, decoys {F[dec, 1].mean():.3f}', flush=True)
