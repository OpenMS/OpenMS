"""V18 reporting (evaluation only): q<=0.01 accepted target winners of D11FWT vs D11FSWT per raw, human / entrapment
counts, paired decoys above the cut, gained / lost precursors and peptides, and historical D / E movement."""
import sys, sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; O = B + '/astral-voxel-pilot-20261110/v2'
sys.path.insert(0, O); 
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
pepmap = pd.read_sql("select PRECURSOR_ID precursor_id, PEPTIDE_ID peptide_id from PRECURSOR_PEPTIDE_MAPPING", c).set_index('precursor_id').peptide_id
def tdc(pop, s):
    w = pop.assign(s=s).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for f_, g in w[w.is_decoy == 1].groupby('fold'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median(); w.loc[w.fold == f_, 's'] = (w.loc[w.fold == f_, 's'] - med) / mad
    best = w.sort_values(['s', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('s', ascending=False)
    T = np.cumsum(1 - best.is_decoy.values); D = np.cumsum(best.is_decoy.values)
    best['q'] = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
    return best
CTL, ARM = sys.argv[1:3]
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t')
    st = pd.read_csv(f'{O}/v4/stageD_v16/score/reference_winners_{raw}.tsv.gz', sep='\t', usecols=['precursor_id', 'stage']).set_index('precursor_id').stage
    acc = {}
    for a in (CTL, ARM):
        b = tdc(pop, np.load(f'{O}/{raw}/scores_{a}.npy')); k = b[b.q <= 0.01]
        acc[a] = (set(k[k.is_decoy == 0].precursor_id), int((k.is_decoy == 1).sum()))
    (tc, dc), (ta, da) = acc[CTL], acc[ARM]
    E = lambda S: sum(int(ent.get(p, 0)) for p in S)
    pep = lambda S: {pepmap[p] for p in S if not ent.get(p, 0)}
    g, l = ta - tc, tc - ta
    print(f'{raw} q<=0.01 {CTL}: targets {len(tc)} human {len(tc) - E(tc)} entrap {E(tc)} paired decoys {dc} human peptides {len(pep(tc))}')
    print(f'{raw} q<=0.01 {ARM}: targets {len(ta)} human {len(ta) - E(ta)} entrap {E(ta)} paired decoys {da} human peptides {len(pep(ta))}')
    print(f'{raw} gained: {len(g)} (human {len(g) - E(g)}, entrap {E(g)}) | lost: {len(l)} (human {len(l) - E(l)}, entrap {E(l)}) | '
          f'new human peptides {len(pep(ta) - pep(tc))}, lost human peptides {len(pep(tc) - pep(ta))}')
    sg = st.reindex(list(g)).value_counts().to_dict(); sl = st.reindex(list(l)).value_counts().to_dict()
    print(f'{raw} gained by V16 stage: {sg} | lost by V16 stage: {sl}', flush=True)
