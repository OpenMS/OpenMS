"""Voxel v2 diagnosis part 4: why does V0 rank entrapment winners above decoy winners?
Hypotheses tested (no training, entrapment labels used only to define the groups):
  H1 borrowed signal: the window sits on a real (base-accepted, human) peptide that co-elutes (|dRT| <= 5 s) with a
     close precursor m/z (<= 2 Th, a proxy for the same DIA isolation window), and the window's library fragments
     match that peptide's fragments (<= 20 ppm).
  H2 MS1: the window has a stronger MS1 isotope trace (from ms1.npy).
  H3 sequence relation: entrapment peptides are shuffled human peptides. Their same-composition human parent has the
     same precursor m/z, so the parent's MS1 and some shared fragments (y1 = C-terminal K/R) are present.
Groups (raw 01 and raw 02, V0 and FS): E_top = entrapment target winners in the arm's top 1500 targets,
D_top = decoy winners with a score >= the 1500th target, E_all / D_all = random 3000 of all entrapment / decoy winners."""
import sqlite3, sys, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = B + '/astral-voxel-pilot-20261110/v2'
R = 148983 / 151239
BASEF = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion/features.tsv',
         'raw02': B + '/astral-raw02-byion-design-20261109/byion/features.tsv'}
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pre = pd.read_sql("""select p.ID precursor_id, p.PRECURSOR_MZ mz, p.CHARGE z, p.DECOY d, pe.UNMODIFIED_SEQUENCE seq,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING pp on pp.PRECURSOR_ID=p.ID
  join PEPTIDE pe on pe.ID=pp.PEPTIDE_ID join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by p.ID""", c).set_index('precursor_id')
tr = pd.read_sql("""select m.PRECURSOR_ID precursor_id, t.PRODUCT_MZ fmz, t.ANNOTATION ann from TRANSITION_PRECURSOR_MAPPING m join TRANSITION t on t.ID=m.TRANSITION_ID""", c)
frag = tr.groupby('precursor_id').fmz.apply(np.sort).to_dict()
ann = tr.groupby('precursor_id').ann.apply(list).to_dict()
comp = pre.seq.map(lambda s: ''.join(sorted(s)))
human_comp = set(zip(comp[(pre.d == 0) & (pre.e == 0)], pre.z[(pre.d == 0) & (pre.e == 0)]))
rng = np.random.default_rng(3)

def winners(pop, s):
    w = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in w[w.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median()
        w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
    return w.sort_values('s', ascending=False)

for raw in sys.argv[1:] or ['raw01', 'raw02']:
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t')
    ms1 = np.load(f'{O}/{raw}/ms1.npy', mmap_mode='r'); xic = np.load(f'{O}/{raw}/xic.npy', mmap_mode='r')
    base = pd.read_csv(BASEF[raw], sep='\t', usecols=['sample_id', 'stage2_score']).set_index('sample_id').stage2_score.loc[pop.sample_id].values
    wb = winners(pop, base); tb = wb[wb.is_decoy == 0].copy(); tb['e'] = tb.precursor_id.map(pre.e).values
    f = np.cumsum(tb.e.values) * (1 + 1 / R) / np.arange(1, len(tb) + 1); n1 = int(np.where(f <= 0.01)[0].max() + 1)
    real = tb.iloc[:n1]; real = real[real.e == 0]  # base-accepted human hits = "real peptides present"
    real = real.assign(mz=real.precursor_id.map(pre.mz).values).sort_values('rt_apex')
    r_rt, r_mz, r_pid = real.rt_apex.values, real.mz.values, real.precursor_id.values
    print(f'\n===== {raw}: {len(real)} base-accepted human hits used as "real signal"')
    for arm in ('V0', 'FS'):
        s = np.load(f'{O}/{raw}/scores_{arm}.npy'); w = winners(pop, s)
        t = w[w.is_decoy == 0].reset_index(drop=True); t['e'] = t.precursor_id.map(pre.e).values; cut = t.s.iloc[1499]
        d = w[w.is_decoy == 1]; real_set = set(r_pid)
        grp = {'E_top': t.iloc[:1500][t.e.iloc[:1500] == 1], 'D_top': d[d.s >= cut],
               'H_top(not base hit)': t.iloc[:1500][(t.e.iloc[:1500] == 0) & ~t.precursor_id.iloc[:1500].isin(real_set)],
               'E_all': t[t.e == 1].sample(3000, random_state=1), 'D_all': d.sample(3000, random_state=1)}
        print(f'--- {arm}')
        print(f'    {"group":22s} {"n":>5s} {"H1 neighbour":>13s} {"shared frags":>13s} {"MS1 max":>9s} {"MS2 sum":>9s} {"H3 parent":>10s} {"same-pep":>9s}')
        for g, x in grp.items():
            nb, sh, m1, m2, par, same = [], [], [], [], [], []
            for pid, rt, row in zip(x.precursor_id.values, x.rt_apex.values, x.row.values):
                lo, hi = np.searchsorted(r_rt, rt - 5), np.searchsorted(r_rt, rt + 5)
                cand = [(r_pid[j]) for j in range(lo, hi) if abs(r_mz[j] - pre.mz[pid]) <= 2.0 and r_pid[j] != pid]
                best = 0
                for q in cand:
                    fq = frag[q]; fp = frag[pid]
                    i = np.searchsorted(fq, fp); i = np.clip(i, 1, len(fq) - 1)
                    dd = np.minimum(np.abs(fq[i] - fp), np.abs(fq[i - 1] - fp)) / fp * 1e6
                    best = max(best, int((dd <= 20).sum()))
                nb.append(len(cand) > 0); sh.append(best)
                same.append(any(abs(pre.mz[q] - pre.mz[pid]) / pre.mz[pid] * 1e6 <= 10 and pre.z[q] == pre.z[pid] for q in cand))
                m1.append(float(np.asarray(ms1[row, 0]).max())); m2.append(float(np.asarray(xic[row]).sum()))
                seqp = pre.seq[pid] if pre.d[pid] == 0 else pre.seq[pid]
                par.append((''.join(sorted(seqp)), pre.z[pid]) in human_comp)
            print(f'    {g:22s} {len(x):5d} {np.mean(nb):13.2f} {np.mean(sh):13.2f} {np.median(m1):9.2g} {np.median(m2):9.2g} {np.mean(par):10.2f} {np.mean(same):9.2f}')
