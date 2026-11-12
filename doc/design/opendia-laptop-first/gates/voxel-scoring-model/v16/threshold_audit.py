"""V16 descriptive threshold semantics and historical-arm churn. No training."""
import json
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path('/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107')
V = ROOT / 'astral-voxel-pilot-20261110/v2'
OUT = V / 'v4/stageD_v16'
R = 148983 / 151239
ARMS = ['D11FWT', 'D11FPWT', 'D11FXWT', 'D11FNWT', 'D11FEWT', 'D11FUWT', 'D11FIWT']


def tdc(pop, scores):
    """Match evaluate_v2.py, including fold normalization and decoy tie priority."""
    winners = pop.assign(s=scores, row=np.arange(len(pop))).sort_values(
        ['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').copy()
    for fold, decoys in winners[winners.is_decoy == 1].groupby('fold'):
        med = decoys.s.median()
        mad = 1.4826 * (decoys.s - med).abs().median()
        if mad <= 0:
            raise ValueError('Degenerate decoy scale')
        mask = winners.fold == fold
        winners.loc[mask, 's'] = (winners.loc[mask, 's'] - med) / mad
    best = winners.sort_values(['s', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('s', ascending=False).copy()
    targets = np.cumsum(1 - best.is_decoy.to_numpy())
    decoys = np.cumsum(best.is_decoy.to_numpy())
    best['q'] = np.minimum.accumulate(np.minimum(1, (decoys + 1) / np.maximum(targets, 1))[::-1])[::-1]
    return winners, best


def n_at(scores, entrap, fdp):
    order = np.argsort(-scores, kind='stable')
    curve = np.cumsum(entrap[order]) * (1 + 1 / R) / np.arange(1, len(entrap) + 1)
    valid = np.flatnonzero(curve <= fdp)
    return int(valid[-1] + 1) if len(valid) else 0


def tail_percentile(values, null):
    """Midrank CDF with half credit for ties."""
    null = np.sort(np.asarray(null))
    return (np.searchsorted(null, values, 'left') + np.searchsorted(null, values, 'right')) / (2 * len(null))


def main():
    with sqlite3.connect(f'file:{ROOT}/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True) as con:
        ent = pd.read_sql('''select pp.PRECURSOR_ID precursor_id,
          min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
          from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID
          join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID''', con).set_index('precursor_id').e
    summary, churn, bands, coverage = [], [], [], []
    for raw in ['raw01', 'raw02']:
        pop = pd.read_csv(V / raw / 'population.tsv', sep='\t')
        stages = pd.read_csv(V / 'v3seq' / f'misses_D11FWT_{raw}.tsv', sep='\t')
        d_ids = set(stages.loc[stages.stage.str.startswith('D'), 'precursor_id'])
        e_ids = set(stages.loc[stages.stage.str.startswith('E'), 'precursor_id'])
        assert len(d_ids) == {'raw01': 866, 'raw02': 929}[raw]
        assert len(e_ids) == {'raw01': 1413, 'raw02': 1362}[raw]
        base, best = tdc(pop, np.load(V / raw / 'scores_D11FWT.npy').astype(float))
        base['e'] = base.precursor_id.map(ent)
        assert base.loc[base.is_decoy == 0, 'e'].notna().all()
        q_index = best.set_index('precursor_id').q
        base['q'] = base.precursor_id.map(q_index)
        base['D'] = base.precursor_id.isin(d_ids)
        base['E'] = base.precursor_id.isin(e_ids)
        base['null_pct'] = np.nan
        for fold, group in base.groupby('fold'):
            base.loc[group.index, 'null_pct'] = tail_percentile(group.s, group.loc[group.is_decoy == 1, 's'])
        targets = base[base.is_decoy == 0].sort_values('s', ascending=False)
        q_acc = base[(base.is_decoy == 0) & (base.q <= .01)]
        n1 = n_at(targets.s.to_numpy(), targets.e.to_numpy(), .01)
        assert n1 == {'raw01': 1496, 'raw02': 1426}[raw]
        assert len(q_acc) == {'raw01': 1614, 'raw02': 1536}[raw]
        empirical_ids = set(targets.head(n1).precursor_id)
        assert not (d_ids & empirical_ids)
        d = base[base.D].copy()
        diann = pd.read_csv(V / 'v4' / f'dd_diann_{raw}.tsv', sep='\t').drop_duplicates('precursor_id').set_index('precursor_id')
        d['diann_q'] = d.precursor_id.map(diann['Q.Value'])
        d['rt_diff'] = d.precursor_id.map(diann['rt_diff'])
        d['high_ref_confidence'] = d.diann_q <= .001
        d['q_band'] = pd.cut(d.q, [0, .01, .05, .15, .5, 1], include_lowest=True).astype(str)
        d.loc[d.q.isna(), 'q_band'] = 'lost_pair_competition'
        d.to_csv(OUT / f'threshold_D_{raw}.tsv', sep='\t', index=False)
        for name, group in d.groupby('q_band'):
            bands.append(dict(raw=raw, band=name, n=len(group), high_diann_confidence=int(group.high_ref_confidence.sum()),
                              median_null_percentile=float(group.null_pct.median()), median_score=float(group.s.median())))
        summary.append(dict(raw=raw, stage_D=len(d), empirical_n1=n1, empirical_cut=float(targets.iloc[n1-1].s),
                            q01_targets=len(q_acc), q01_human=int((q_acc.e == 0).sum()), q01_entrap=int(q_acc.e.sum()),
                            q01_pair_decoys=int(((best.is_decoy == 1) & (best.q <= .01)).sum()),
                            D_q01=int((d.q <= .01).sum()), D_lost_pair=int(d.q.isna().sum()),
                            D_null_below90=int((d.null_pct < .9).sum()), D_null_above999=int((d.null_pct >= .999).sum()),
                            D_null_percentiles=np.quantile(d.null_pct, [.1,.25,.5,.75,.9]).tolist()))
        # Diagnostic cost of merely lowering the unchanged baseline cut. This is not a proposed acceptance rule.
        for fraction in [.25, .5, .75, 1.]:
            cutoff = float(d.s.sort_values(ascending=False).iloc[int(np.ceil(len(d) * fraction)) - 1])
            accepted = targets[targets.s >= cutoff]
            coverage.append(dict(raw=raw, D_coverage_fraction=fraction, cutoff=cutoff, D_above=int((d.s >= cutoff).sum()),
                                 human=int((accepted.e == 0).sum()), entrap=int(accepted.e.sum()),
                                 decoy_precursor_winners=int(((base.is_decoy == 1) & (base.s >= cutoff)).sum()),
                                 empirical_fdp=float(accepted.e.sum() * (1+1/R) / len(accepted))))
        base_acc = set(q_acc.precursor_id)
        for arm in ARMS:
            w, b = tdc(pop, np.load(V / raw / f'scores_{arm}.npy').astype(float))
            b['e'] = b.precursor_id.map(ent)
            acc = b[(b.is_decoy == 0) & (b.q <= .01)]
            ids = set(acc.precursor_id)
            gained = acc[~acc.precursor_id.isin(base_acc)]
            lost = q_acc[~q_acc.precursor_id.isin(ids)]
            churn.append(dict(raw=raw, arm=arm, targets=len(acc), entrap=int(acc.e.sum()),
                              decoy_competitors=int(((b.is_decoy == 1) & (b.q <= .01)).sum()),
                              D_accepted=len(d_ids & ids), D_new_vs_base=len((d_ids & ids) - base_acc),
                              baseline_q_targets_lost=len(base_acc-ids), baseline_E_lost=len(e_ids-ids),
                              gained_human=int((gained.e == 0).sum()), gained_entrap=int(gained.e.sum()),
                              lost_human=int((lost.e == 0).sum()), lost_entrap=int(lost.e.sum()),
                              empirical_fdp=float(acc.e.sum()*(1+1/R)/max(len(acc),1))))
        print(raw, json.dumps(summary[-1]), flush=True)
    pd.DataFrame(summary).to_json(OUT/'threshold_summary.json', orient='records', indent=2)
    pd.DataFrame(bands).to_csv(OUT/'threshold_bands.tsv', sep='\t', index=False)
    pd.DataFrame(coverage).to_csv(OUT/'threshold_lowering_cost.tsv', sep='\t', index=False)
    pd.DataFrame(churn).to_csv(OUT/'historical_churn.tsv', sep='\t', index=False)
    print(pd.DataFrame(bands).to_string(index=False))
    print(pd.DataFrame(coverage).to_string(index=False))
    print(pd.DataFrame(churn).to_string(index=False))


if __name__ == '__main__':
    main()
