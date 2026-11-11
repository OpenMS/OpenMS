"""Stage-D examples PDF (raw 01, old population, D11FWT). Descriptive only; DIA-NN is used as a reference only.
Per precursor:
  (a) MS2 XICs of the 6 library fragments from the OpenDIA .xipm export (pyopenms XIPMParquetFile), summed per scan within
      +-10 ppm of the calibrated target m/z, apex +-40 s;
  (b) MS1 XICs M+0..M+3 (+-10 ppm after the calibrated offset);
  (c) the MS2 scan nearest our apex whose isolation window holds the precursor m/z (raw01.mzML, OnDiscMSExperiment),
      mirrored against the library intensities; peaks within 10 ppm of a calibrated library fragment are annotated.
Lines: our winning apex (solid) and DIA-NN RT with RT.Start..RT.Stop (dashed, shaded) when DIA-NN reports the precursor.
Groups are chosen by fixed rules with a fixed seed (see pick())."""
import json, sqlite3, numpy as np, pandas as pd, pyarrow.parquet as pq, pyopenms as oms
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; V = B + '/astral-voxel-pilot-20261110/v2'; OUT = B + '/dd-plots-20261111'
RAW = 'raw01'; PPM = 10.0; HALF = 40.0; RR = 148983 / 151239
OFF = float(json.load(open(B + '/astral-design-raw01-20261107/window300k/cal.json'))['ms2_mz']['offset_ppm'])
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
prec = pd.read_sql("""select p.ID precursor_id, p.PRECURSOR_MZ pmz, p.CHARGE z, p.LIBRARY_RT lrt, pe.MODIFIED_SEQUENCE seq from PRECURSOR p
  join PRECURSOR_PEPTIDE_MAPPING pp on pp.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=pp.PEPTIDE_ID""", c).drop_duplicates('precursor_id').set_index('precursor_id')
tr = pd.read_sql("""select m.PRECURSOR_ID precursor_id, t.PRODUCT_MZ mz, t.ANNOTATION ann, t.LIBRARY_INTENSITY li from TRANSITION_PRECURSOR_MAPPING m
  join TRANSITION t on t.ID=m.TRANSITION_ID""", c)
DN = pq.read_table('/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet',
                   columns=['Run', 'Stripped.Sequence', 'Precursor.Charge', 'RT', 'RT.Start', 'RT.Stop', 'Precursor.Quantity', 'Q.Value', 'Evidence', 'Decoy']).to_pandas()
DN = DN[(DN.Decoy == 0) & DN.Run.str.endswith('_01')]; DN['key'] = DN['Stripped.Sequence'] + '/' + DN['Precursor.Charge'].astype(str)
DN = DN.drop_duplicates('key').set_index('key')

# Winners with the evaluator's per-fold decoy median/MAD normalisation; the FDP_c 1 % cut as in misses.py.
pop = pd.read_csv(f'{V}/{RAW}/population.tsv', sep='\t'); s = np.load(f'{V}/{RAW}/scores_D11FWT.npy').astype(np.float64)
w = pop.assign(s=s).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
for f_, g in w[w.is_decoy == 1].groupby('fold'):
    med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
t = w[w.is_decoy == 0].assign(e=lambda d: d.precursor_id.map(ent).values).sort_values('s', ascending=False)
fdp = np.cumsum(t.e.values) * (1 + 1 / RR) / np.arange(1, len(t) + 1); CUT = t.s.values[int(np.where(fdp <= 0.01)[0].max())]
w = w.set_index('precursor_id')
st = pd.read_csv(f'{V}/v4/dd_diann_{RAW}.tsv', sep='\t').set_index('precursor_id')  # E / D with DIA-NN columns and evidence measures
band = pd.read_csv(f'{V}/v4/dd_band_{RAW}.tsv', sep='\t').set_index('precursor_id')  # D / DEC / ENT winners above the decoy q90

def pick():
    rng = np.random.default_rng(11); e, d = st[st.grp == 'E'], st[st.grp == 'D']
    t1, t2 = st.lq.quantile([1 / 3, 2 / 3])
    take = lambda df, k: list(df.index[rng.choice(len(df), min(k, len(df)), replace=False)]) if len(df) else []
    dsc = w.loc[d.index, 's']; lo, hi = dsc.quantile([0.25, 0.75])
    near = lambda g: band[(band.grp == g)].index[(w.loc[band[band.grp == g].index, 's'] >= lo).values & (w.loc[band[band.grp == g].index, 's'] <= hi).values]
    P = [('E accepted, high abundance', take(e[e.lq >= t2], 2)), ('E accepted, low abundance', take(e[e.lq < t1], 2)),
         ('D missed, low abundance, DIA-NN q <= 0.001', take(d[(d.lq < t1) & (d['Q.Value'] <= 0.001)], 3)),
         ('D missed, low abundance, DIA-NN q > 0.003', take(d[(d.lq < t1) & (d['Q.Value'] > 0.003)], 2)),
         ('D missed, mid abundance', take(d[(d.lq >= t1) & (d.lq < t2)], 2)), ('D missed, high abundance', take(d[d.lq >= t2], 1)),
         ('Decoy winner, score in D interquartile range', take(pd.DataFrame(index=near('DEC')), 3)),
         ('Entrapment winner, score in D interquartile range', take(pd.DataFrame(index=near('ENT')), 2))]
    return P, (t1, t2), (lo, hi)

def xics(xf, pid):
    d = xf.getPeakMaps(precursor_id=int(pid), explode=False); ms2, ms1 = [], []
    for k in range(len(d['ms_level'])):
        rt, mz, it = np.asarray(d['rt'][k]), np.asarray(d['mz'][k]), np.asarray(d['intensity'][k]); tg = d['target_mz'][k]
        if len(rt) == 0: continue
        if d['ms_level'][k] == 1:
            iso = int(d['annotation'][k].rsplit('_i', 1)[1])
            if iso > 3: continue
            m = np.abs((mz - tg) / tg * 1e6 - OFF) <= PPM; lab = f'M+{iso}'; dst = ms1
        else:
            m = np.abs((mz - tg) / tg * 1e6) <= PPM; lab = d['annotation'][k]; dst = ms2
        df = pd.DataFrame({'rt': rt[m], 'int': it[m]}).groupby('rt', as_index=False).int.sum()
        allrt = pd.DataFrame({'rt': np.unique(rt)}).merge(df, how='left').fillna(0.0)  # zero where a scan had no point in tolerance
        dst.append(allrt.assign(ion=lab, target_mz=tg))
    return (pd.concat(ms2) if ms2 else pd.DataFrame(columns=['rt', 'int', 'ion'])), (pd.concat(ms1) if ms1 else pd.DataFrame(columns=['rt', 'int', 'ion'])), d

def apex_spectrum(od, ms2idx, pmz, rt):
    c_ = ms2idx[(ms2idx.lo <= pmz) & (ms2idx.hi >= pmz)]
    j = int(c_.iloc[(c_.rt - rt).abs().argsort().values[0]].idx); sp = od.getSpectrum(j); mz, it = sp.get_peaks()
    return pd.DataFrame({'mz': mz, 'int': it}), sp.getRT()

def page(pdf, title, pid, xf, od, ms2idx):
    p = prec.loc[pid]; seq = p.seq; key = ''.join(ch for ch in seq.split('(')[0]) if '(' not in seq else None
    import re; stripped = re.sub(r'\(.*?\)', '', seq); key = f'{stripped}/{int(p.z)}'; dn = DN.loc[key] if key in DN.index else None
    apex = float(w.loc[pid, 'rt_apex']); sc = float(w.loc[pid, 's'])
    ms2, ms1, _ = xics(xf, pid)
    fig = plt.figure(figsize=(11.7, 8.3)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], hspace=0.38, wspace=0.28)
    ax1, ax2, ax3 = fig.add_subplot(gs[0, :2]), fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[1, :])
    for ax, df, name in ((ax1, ms2, 'MS2 fragment XICs (6 library fragments)'), (ax2, ms1, 'MS1 isotopes')):
        sub = df[(df.rt > apex - HALF) & (df.rt < apex + HALF)].copy()
        if len(sub):
            sub.plot(kind='chromatogram', x='rt', y='int', by='ion', backend='ms_matplotlib', canvas=ax, show_plot=False,
                     legend_config={'show': True, 'loc': 'right', 'fontsize': 7}, grid=False, xlabel='RT (s)', ylabel='intensity', title=name)
        ax.axvline(apex, color='k', lw=1, label='our apex')
        if dn is not None:
            ax.axvline(dn.RT * 60, color='r', ls='--', lw=1); ax.axvspan(dn['RT.Start'] * 60, dn['RT.Stop'] * 60, color='r', alpha=0.08)
        ax.set_xlim(apex - HALF, apex + HALF); ax.set_title(name, fontsize=9)
    spec, srt = apex_spectrum(od, ms2idx, p.pmz, apex)
    lib = tr[tr.precursor_id == pid].copy(); lib['cmz'] = lib.mz * (1 + OFF * 1e-6)
    lo_, hi_ = lib.cmz.min() - 50, lib.cmz.max() + 50; spec = spec[(spec.mz > max(lo_, 100)) & (spec.mz < hi_)].copy()
    spec['ion'] = ''
    for _, r in lib.iterrows():
        m = np.abs(spec.mz - r.cmz) / r.cmz * 1e6 <= PPM
        if m.any(): j = spec[m].int.idxmax(); spec.loc[j, 'ion'] = r.ann
    spec['kind'] = np.where(spec.ion != '', 'library fragment', 'other')
    ref = pd.DataFrame({'mz': lib.cmz, 'int': lib.li, 'ion': lib.ann, 'kind': 'library fragment'})
    spec.plot(kind='spectrum', x='mz', y='int', by='kind', backend='ms_matplotlib', canvas=ax3, reference_spectrum=ref, mirror_spectrum=True,
              ion_annotation='ion', annotate_mz=False, show_plot=False, grid=False, xlabel='m/z', ylabel='relative intensity (top: scan, bottom: library)',
              title=f'MS2 scan at RT {srt:.1f} s (nearest to our apex), isolation window holds precursor m/z {p.pmz:.4f}', legend_config={'show': True, 'fontsize': 7})
    ax3.set_title(ax3.get_title(), fontsize=9)
    n_hit = int((spec.ion != '').sum()); dec = 'decoy' if w.loc[pid, 'is_decoy'] == 1 else ('entrapment target' if ent.get(pid, 0) == 1 else 'human target')
    ex = st.loc[pid] if pid in st.index else (band.loc[pid] if pid in band.index else None)
    meas = '' if ex is None else f"our lib_cos {ex.lib_cos:.2f}, corr_med {ex.corr_med:.2f}, log10 central area {ex.log_area:.2f}"
    dtxt = 'not reported by DIA-NN' if dn is None else (f"DIA-NN: RT {dn.RT*60:.1f} s, log10 Precursor.Quantity {np.log10(dn['Precursor.Quantity']):.2f}, "
                                                          f"q {dn['Q.Value']:.4f}, Evidence {dn.Evidence:.2f}")
    fig.suptitle(f"{title}\n{seq} {int(p.z)}+ (precursor {pid}, {dec}); our score {sc:.1f} vs cut {CUT:.1f}; {meas}\n{dtxt}; "
                 f"{n_hit} of {len(lib)} library fragments found in the apex scan", fontsize=9.5, x=0.01, ha='left')
    pdf.savefig(fig); plt.close(fig)

def summary(pdf, tert, iqr):
    fig, ax = plt.subplots(1, 3, figsize=(11.7, 8.3)); fig.subplots_adjust(top=0.70, wspace=0.3)
    e, d = st[st.grp == 'E'], st[st.grp == 'D']
    allq = np.log10(DN['Precursor.Quantity'])
    bins = np.linspace(2.5, 9, 40)
    ax[0].hist(e.lq, bins, alpha=0.6, label=f'E accepted (n {len(e)})'); ax[0].hist(d.lq, bins, alpha=0.6, label=f'D missed (n {len(d)})')
    for x in tert: ax[0].axvline(x, color='k', ls=':')
    ax[0].set_xlabel('DIA-NN log10 Precursor.Quantity'); ax[0].set_ylabel('precursors'); ax[0].legend(fontsize=8); ax[0].set_title('abundance (dotted: tertiles)', fontsize=9)
    ax[1].scatter(e.lq, e.log_area, s=4, alpha=0.4, label='E'); ax[1].scatter(d.lq, d.log_area, s=4, alpha=0.4, label='D')
    ax[1].set_xlabel('DIA-NN log10 Precursor.Quantity'); ax[1].set_ylabel('our log10 central fragment area'); ax[1].legend(fontsize=8); ax[1].set_title('same signal seen by both', fontsize=9)
    g = st.assign(qb=pd.qcut(st.lq, 10, labels=False)).groupby('qb').agg(lq=('lq', 'median'), D=('grp', lambda v: (v == 'D').mean()))
    ax[2].plot(g.lq, g.D, 'o-'); ax[2].set_xlabel('DIA-NN log10 Precursor.Quantity (decile median)'); ax[2].set_ylabel('share missed (D)'); ax[2].set_ylim(0, 1)
    ax[2].set_title('miss share by abundance decile', fontsize=9)
    txt = (f"Stage-D examples, raw 01, D11FWT (old population). Accepted at FDP_c 1 %: score cut {CUT:.1f} (normalised winner score).\n"
           f"'Low abundance' = lowest DIA-NN quantity tertile: log10 Precursor.Quantity < {tert[0]:.2f} (about {10**tert[0]:.0f}); high tertile > {tert[1]:.2f}.\n"
           f"All DIA-NN 1 % precursors in this run: median log10 quantity {allq.median():.2f}. D score interquartile range {iqr[0]:.1f}..{iqr[1]:.1f}.\n"
           f"Per page: (top left) 6 library-fragment XICs from the .xipm export, (top right) MS1 M+0..M+3, (bottom) apex MS2 scan mirrored with the library.\n"
           f"Black line: our winning apex. Red dashed line and band: DIA-NN RT and RT.Start..RT.Stop.")
    fig.text(0.01, 0.97, txt, va='top', fontsize=9.5); pdf.savefig(fig); plt.close(fig)

if __name__ == '__main__':
    P, tert, iqr = pick()
    xf = oms.XIPMParquetFile(f'{V}/export01/v2.xipm')
    od = oms.OnDiscMSExperiment(); od.openFile(f'{OUT}/raw01.mzML', False); meta = od.getMetaData()
    rows = []
    for i in range(meta.getNrSpectra()):
        sp = meta.getSpectrum(i)
        if sp.getMSLevel() != 2 or not sp.getPrecursors(): continue
        pr = sp.getPrecursors()[0]; rows.append((i, sp.getRT(), pr.getMZ() - pr.getIsolationWindowLowerOffset(), pr.getMZ() + pr.getIsolationWindowUpperOffset()))
    ms2idx = pd.DataFrame(rows, columns=['idx', 'rt', 'lo', 'hi'])
    sel = [(t_, pid) for t_, ids in P for pid in ids]
    pd.DataFrame(sel, columns=['group', 'precursor_id']).to_csv(f'{OUT}/examples.tsv', sep='\t', index=False)
    with PdfPages(f'{OUT}/stageD_examples_raw01.pdf') as pdf:
        summary(pdf, tert, iqr)
        for t_, pid in sel:
            print(t_, pid, flush=True); page(pdf, t_, pid, xf, od, ms2idx)
