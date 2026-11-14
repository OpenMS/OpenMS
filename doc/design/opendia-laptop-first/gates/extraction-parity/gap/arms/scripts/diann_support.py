"""Post-hoc independent check (evaluation only): are the IDs that block I newly accepts supported by DIA-NN?
For each raw: accepted sets at the N@FDP_c 1 % cut of the control D11FDWT and of D11FDIWT (frozen seed). For the
human targets gained by I (accepted by I, not by the control) and lost, report DIA-NN's evidence: share present in
DIA-NN's report at all, DIA-NN q <= 0.01 / 0.05, and |apex - DIA-NN RT| <= 6 s for the winning window. Entrapment gains
are counted separately (they are false by construction). A gain supported by DIA-NN at least as often as the control's
marginal accepts (the last 200 above its cut) indicates real peptides rather than ranking noise."""
import os, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114'
SRC = open(B + '/funnel.py').read().split("pp = pop[pop.precursor_id")[0]
REF = {'raw01': '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/diann-pairlib-20261113/report_raw01.parquet',
       'raw02': '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/diann-decoys-20261114/own/raw02/report.parquet'}
for raw in ('raw01', 'raw02'):
    os.environ['RAW'] = raw
    res = {}
    for arm in ('D11FDWT', 'D11FDIWT'):
        g = {}
        exec(SRC.replace("np.load(f'{D}/scores_D11FDWT.npy')", f"np.load('{B}/run/{raw}/scores_{arm}.npy')"), g)
        res[arm] = g
    g = res['D11FDWT']; tg = g['tg'].reset_index().set_index('precursor_id'); ent = g['ent']
    rep = pd.read_parquet(REF[raw], columns=['Modified.Sequence', 'Precursor.Charge', 'RT', 'Q.Value'])
    rep['key'] = rep['Modified.Sequence'] + '/' + rep['Precursor.Charge'].astype(str)
    rep = rep.sort_values('Q.Value').drop_duplicates('key').set_index('key')
    def support(ids, arm):
        ids = [p for p in ids if ent.get(p, 0) == 0]
        w = res[arm]['w'].set_index('precursor_id')
        keys = tg.loc[ids, 'key']
        r = rep.reindex(keys.values)
        present = r['Q.Value'].notna().values
        near = present & (np.abs(w.loc[ids, 'rt_apex'].values - r.RT.values * 60) <= 6)
        return (f'n {len(ids):4d} | in DIA-NN report {present.mean():.2f} | q<=0.01 {(r["Q.Value"] <= 0.01).mean():.2f} '
                f'| q<=0.05 {(r["Q.Value"] <= 0.05).mean():.2f} | same peak (<=6 s) {near.mean():.2f}')
    ac, ai = res['D11FDWT']['accepted'], res['D11FDIWT']['accepted']
    gained, lost = ai - ac, ac - ai
    ent_g = sum(ent.get(p, 0) for p in gained); ent_l = sum(ent.get(p, 0) for p in lost)
    marg = list(res['D11FDWT']['tw'].precursor_id.iloc[res['D11FDWT']['n1'] - 200:res['D11FDWT']['n1']])
    print(f'{raw}: control accepts {len(ac)}, I accepts {len(ai)}; gained {len(gained)} (entrapment {ent_g}), lost {len(lost)} (entrapment {ent_l})')
    print('  gained by I (human)    ', support(gained, 'D11FDIWT'))
    print('  lost by I (human)      ', support(lost, 'D11FDWT'))
    print('  control marginal 200   ', support(marg, 'D11FDWT'))
    print('  control all accepted   ', support(ac, 'D11FDWT'))
