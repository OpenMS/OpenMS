"""Direct test of the parent-co-elution mechanism (evaluation only).
A winner has an eluting parent if the ProteoBench DIA-NN reference (same raw, q <= 0.01) identified a peptide with the
same amino-acid composition and charge (so the same precursor m/z, within 10 ppm) whose RT is within 10 s of the
winner's apex. Entrapments are FDRBench shuffles of human peptides; their parents are in the full proteome and are
often identified by the reference, but are mostly absent from our 100K library. Decoys are mass shifted, so for a decoy
the same test asks whether its own target's composition elutes at the decoy's apex with the decoy's m/z (should be rare).
Also: rho among winners with and without such a parent, and the background rate in the whole population."""
import sys, numpy as np, pandas as pd
sys.argv = [sys.argv[0], sys.argv[1] if len(sys.argv) > 1 else 'D11FWT', '../pairlib-20261113', 'xic']  # PLAN_PAIRLIB
exec(open('diagnose.py').read().split("for raw in ('raw01', 'raw02'):")[0])
ref = pd.read_parquet('/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/input_file.parquet',
                      columns=['Run', 'Stripped.Sequence', 'Precursor.Charge', 'Precursor.Mz', 'RT', 'Q.Value', 'Decoy', 'Protein.Ids'])
ref = ref[(ref['Q.Value'] <= 0.01) & (ref.Decoy == 0)]
seq = pd.read_sql("""select p.ID precursor_id, pe.UNMODIFIED_SEQUENCE seq from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING pp on pp.PRECURSOR_ID=p.ID
  join PEPTIDE pe on pe.ID=pp.PEPTIDE_ID""", c).drop_duplicates('precursor_id').set_index('precursor_id').seq
comp = lambda s: ''.join(sorted(s))
names = {0: 'human', 1: 'entrap', 2: 'dec(hum)', 3: 'dec(ent)'}
for raw, run in (('raw01', '_01'), ('raw02', '_02')):
    rr = ref[ref.Run.str.endswith(run)].copy()
    rr['key'] = rr['Stripped.Sequence'].map(comp) + '/' + rr['Precursor.Charge'].astype(str)
    rr['rts'] = rr.RT * 60.0
    parents = rr.groupby('key').agg(pmz=('Precursor.Mz', 'first'), prt=('rts', list), pseq=('Stripped.Sequence', lambda x: set(x)))
    pop = pd.read_csv(f'{ROOT}/{raw}/population.tsv', sep='\t')
    s = np.load(f'{ROOT}/{raw}/scores_{ARM}.npy').astype(np.float64)
    tgt_of_pair = pop[pop.is_decoy == 0].drop_duplicates('pair_key').set_index('pair_key').precursor_id
    pop['cls'] = np.where(pop.is_decoy == 0, pop.precursor_id.map(ent), 2 + pop.pair_key.map(tgt_of_pair).map(ent).fillna(-9))
    w, best = tdc(pop, s)
    best = best.join(pre[['charge', 'mz']], on='precursor_id')
    best['seq'] = best.precursor_id.map(seq)
    best['key'] = best.seq.map(comp) + '/' + best.charge.astype(str)
    j = best.join(parents, on='key')

    def has_parent(r, drt=10.0):
        if not isinstance(r.prt, list):
            return 0
        if abs(r.pmz - r.mz) > r.mz * 10e-6:   # mass shifted (decoys): no co-isolated parent
            return 0
        if r.seq in r.pseq:                     # the winner itself is the reference ID, not a different parent
            others = True
        return int(any(abs(t - r.rt_apex) <= drt for t in r.prt))

    def self_is_ref(r):
        return int(isinstance(r.pseq, set) and r.seq in r.pseq)

    j['parent'] = [has_parent(r) for r in j.itertuples()]
    j['self_ref'] = [self_is_ref(r) for r in j.itertuples()]
    # A parent that is a DIFFERENT peptide (shuffle partner), eluting at the apex.
    j['other_parent'] = ((j.parent == 1) & (j.self_ref == 0)).astype(int)
    print(f'==== {raw}')
    for q in (0.01, 0.05, 0.10, 0.20):
        cut = j[j.q <= q].s.min(); top = j[j.s >= cut]
        line = f'  q <= {q:.2f}: winners with a DIFFERENT reference peptide of the same composition+charge eluting within 10 s:'
        print(line)
        for k_, g in top.groupby('cls'):
            if k_ in names:
                print(f'    {names[k_]:9s} n {len(g):5d}  {g.other_parent.sum():4d} ({g.other_parent.mean():.3f})')
        for flag, g in top.groupby('other_parent'):
            E = int((g.cls == 1).sum()); D = int((g.is_decoy == 1).sum()); N = int((g.is_decoy == 0).sum())
            print(f'    other_parent={flag}:' + fmt(E, D, N))
    # Background over the population: a parent of the same composition+charge+m/z identified anywhere in the run.
    allp = pop.drop_duplicates('precursor_id').join(pre[['charge', 'mz']], on='precursor_id')
    allp['seq'] = allp.precursor_id.map(seq); allp['key'] = allp.seq.map(comp) + '/' + allp.charge.astype(str)
    a = allp.join(parents, on='key')
    a['any_parent'] = [int(isinstance(r.prt, list) and abs(r.pmz - r.mz) <= r.mz * 10e-6 and not (r.seq in r.pseq)) for r in a.itertuples()]
    print('  background: population precursors with a different same-composition reference peptide anywhere in the run:',
          {names.get(k, k): round(v, 4) for k, v in a.groupby('cls').any_parent.mean().items()})
    ex = j[(j.q <= 0.01) & (j.cls == 1)][['precursor_id', 'seq', 'charge', 'rt_apex', 'q', 'other_parent', 'pseq', 'prt']]
    print('  q <= 0.01 entrapment winners:')
    for r in ex.itertuples():
        near = [round(t - r.rt_apex, 1) for t in (r.prt if isinstance(r.prt, list) else [])]
        print(f'    {r.seq}/{r.charge} apex {r.rt_apex:.0f} parent {sorted(r.pseq) if isinstance(r.pseq, set) else "-"} dRT {near[:3]}')
