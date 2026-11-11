"""A+B fix check on the extended-calibration search (rtx/<raw>/search.features.tsv); DIA-NN and entrapment used only to
measure. For the old and new search tables: candidate windows, population (pair max transition_support >= 3) size by
target/decoy, and the DIA-NN reference funnel up to 'a population window within 10 s of the DIA-NN RT' (stages A, B
and 'has right window'). Also the entrapment share among population target precursors (a symmetric-change check)."""
import sys, sqlite3, numpy as np, pandas as pd
R = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; V = R + '/astral-voxel-pilot-20261110/v2'
OLD = {'raw01': R + '/astral-design-raw01-20261107/pipeline', 'raw02': R + '/astral-entrapment-raw02-20261108'}
c = sqlite3.connect('file:' + R + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pre = pd.read_sql('select ID precursor_id, TRAML_ID traml, DECOY decoy from PRECURSOR', c)
pre['pair_key'] = np.where(pre.decoy == 1, pre.traml.str.replace('^DECOY_', '', regex=True), pre.traml)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c).set_index('precursor_id').e
ref_all = pd.read_parquet(R + '/astral-diann-gap-20261108/ref_exact.parquet')
for raw in sys.argv[1:] or ['raw01', 'raw02']:
    ref = ref_all[ref_all.Run.str.endswith('_' + raw[-2:]) & ~ref_all.entrap].drop_duplicates('precursor_id').copy(); ref['rt_s'] = ref.RT * 60
    print(f'\n== {raw}: DIA-NN reference {len(ref)}')
    for name, path in (('old', OLD[raw]), ('extended', f'{V}/rtx/{raw}')):
        cand = pd.read_csv(path + '/search.features.tsv', sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'transition_support'])
        cand = cand.merge(pre[['precursor_id', 'pair_key']], on='precursor_id')
        pm = cand.groupby('pair_key').transition_support.max(); pop = cand[cand.pair_key.isin(pm[pm >= 3].index)]
        pp = pop.drop_duplicates('precursor_id'); tp = pp[pp.is_decoy == 0]
        r = ref.merge(pop[['precursor_id', 'rt_apex']], on='precursor_id', how='left')
        right = set(r[(r.rt_apex - r.rt_s).abs() <= 10].precursor_id)
        inpop = set(pop.precursor_id)
        A = int((~ref.precursor_id.isin(inpop)).sum()); Bn = int((ref.precursor_id.isin(inpop) & ~ref.precursor_id.isin(right)).sum())
        print(f'{name:9s} windows {len(cand):8d} | population windows {len(pop):8d}, precursors T {len(tp):6d} D {int((pp.is_decoy == 1).sum()):6d}, '
              f'entrapment share of T {tp.precursor_id.map(ent).mean():.4f} | A {A} B {Bn} has right window {len(right)}')
