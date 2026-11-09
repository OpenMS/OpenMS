"""Raw 02 gap anatomy on the re-extracted XICs (h8, diagnostic). OSW+Percolator reaches 1287 at entrapment FDP_c 1% on
raw 02, our 5-seed MLP 1083. Set G = human (non-entrapment) targets in OSW's export at m_score<=0.01 (reference only)
that are NOT in our accepted set at entrapment FDP_c 1%. For each, using OSW's reported RT:
  - no candidate: none of our candidates lies within 10 s of OSW's RT (candidate generation or extraction);
  - wrong winner: a candidate is within 10 s but our winner is another peak;
  - right winner, low score: our winner is within 10 s but below our FDP_c 1% cut.
For the last two: the median rank of the OSW-matching candidate within the precursor, and its score gap to the cut in
units of the decoy-winner MAD. Also: the same split for the accepted human set (as a reference), and the share of G
that DIA-NN also identifies.
Usage: OPENDIA_RUN=R02 (R01 compares against OSW raw 01)."""
import os, sys, sqlite3, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../gap2')
import oswlike_score as S
RUN = os.environ['OPENDIA_RUN']; RUNNO = RUN[-2:]
_sel = S.select
def select(idx, score, fdr):
    pos, neg = _sel(idx, score, fdr)
    if len(pos) < 200:
        w = S.winners(idx, score); tw = w[~S.dec[w]]; pos = tw[np.argsort(-score[tw], kind='stable')[:200]]
    return pos, neg
S.select = select
meta = S.m; fam = meta.family.values; keep = (fam == 'T') | (fam == 'O')
for k_, v in {k: getattr(S, k).copy() for k in ('X0', 'prec', 'dec', 'fold', 'pair', 'entr', 'right')}.items(): setattr(S, k_, v[keep])
m = meta[keep].reset_index(drop=True); apex = m.apex.values; dnid = m.dn_rt.notna().values
sc = np.mean([S.crossfit(S.ALL, kind='mlp', seed=s) for s in (0, 10, 20, 30, 40)], axis=0)
wn = S.winners(np.arange(len(sc)), sc); tg = wn[~S.dec[wn]]; o = tg[np.argsort(-sc[tg], kind='stable')]
fd = np.cumsum(S.entr[o]) * (1 + 1 / S.R_ENT) / np.arange(1, len(o) + 1); nacc = int(np.where(fd <= 0.01)[0].max() + 1)
acc = o[:nacc]; cut = sc[o[nacc - 1]]
dwin = wn[S.dec[wn]]; mad = 1.4826 * np.median(np.abs(sc[dwin] - np.median(sc[dwin])))
# map merged-PQP target precursor -> modified sequence/charge, to match OSW's export
c = sqlite3.connect('file:' + os.path.dirname(os.path.abspath(__file__)) + '/merged.pqp?mode=ro', uri=True)
pp = pd.read_sql('select p.ID pid, p.CHARGE z, pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING mm on mm.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=mm.PEPTIDE_ID where p.DECOY=0', c).drop_duplicates('pid')
key2pid = dict(zip(pp.ms + '/' + pp.z.astype(str), pp.pid))
osw = pd.read_csv('/home/singjc/Documents/datasets/proteobench-astral-entrapment/100K_sampled_fasta/opendia_results/OpenDIA.results.tsv', sep='\t',
                  usecols=['run_name', 'FullPeptideName', 'Charge', 'RT', 'decoy', 'm_score', 'ProteinName', 'peak_group_rank'])
osw = osw[osw.run_name.str.endswith('_' + RUNNO) & (osw.decoy == 0) & (osw.m_score <= 0.01) & (osw.peak_group_rank == 1)]
osw['ent'] = osw.ProteinName.str.split(';').map(lambda a: all('_p_target' in x for x in a))
osw = osw[~osw.ent]; osw['pid'] = (osw.FullPeptideName + '/' + osw.Charge.astype(str)).map(key2pid)
print(f'run {RUN}: our N at FDP_c 1% {nacc}; OSW human m_score<=0.01 {len(osw)}, mapped {osw.pid.notna().mean():.3f}', flush=True)
osw = osw.dropna(subset=['pid']); osw['pid'] = osw.pid.astype(int)
byp = pd.Series(np.arange(len(sc))).groupby(S.prec).apply(np.array).to_dict()
win_of = dict(zip(S.prec[wn], wn)); acc_set = set(S.prec[acc])
rows = []
for pid, rt in zip(osw.pid, osw.RT):
    idx = byp.get(pid)
    if idx is None: rows.append(('not extracted', np.nan, np.nan, False)); continue
    near = idx[np.abs(apex[idx] - rt) <= 10]
    w = win_of[pid]; d = dnid[idx[0]]
    if pid in acc_set: cat = 'accepted'
    elif len(near) == 0: cat = 'no candidate near OSW RT'
    elif abs(apex[w] - rt) <= 10: cat = 'right winner, below cut'
    else: cat = 'wrong winner'
    if len(near):
        b = near[np.argmax(sc[near])]; rank = int((sc[idx] > sc[b]).sum()) + 1; gap = (cut - sc[b]) / mad
    else: rank, gap = np.nan, np.nan
    rows.append((cat, rank, gap, d))
r = pd.DataFrame(rows, columns=['cat', 'rank', 'gap', 'dn'])
g = r.groupby('cat').agg(n=('cat', 'size'), dn_share=('dn', 'mean'), rank_med=('rank', 'median'), gap_med=('gap', 'median'))
print(g.to_string(), flush=True)
miss = r[r.cat != 'accepted']
print(f'missed OSW human IDs: {len(miss)}; of these with the OSW peak at our rank 1: {(miss["rank"] == 1).mean():.3f}; within 1 MAD of the cut: {(miss.gap <= 1).mean():.3f}', flush=True)
