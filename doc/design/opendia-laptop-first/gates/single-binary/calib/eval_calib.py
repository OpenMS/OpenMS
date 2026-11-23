"""Descriptive calibration curves. No independence or FDR-control claim.
Usage: eval_calib.py LIBRARY RUN_LAYOUT SET OUTPUT_PREFIX
Layout SET_{never,always}_s{0..4}/*.precursors.tsv[.gz].
Pairing by composition is NOT an authenticated generator pairing map. Until one
is recovered, report lower/combined only. Repeated seeds reuse the same raw:
pooled counts are descriptive, not independent Poisson replicates.
"""
import glob
import os
import sqlite3
import sys
import pandas as pd

ALPHAS = [0.001, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.05]


def evaluate(library, layout, dataset):
    with sqlite3.connect(f'file:{library}?mode=ro', uri=True) as conn:
        labels = pd.read_sql_query("""SELECT p.ID precursor_id,
          min(CASE WHEN pr.PROTEIN_ACCESSION LIKE '%_p_target%' THEN 1 ELSE 0 END) e
          FROM PRECURSOR p JOIN PRECURSOR_PEPTIDE_MAPPING m ON m.PRECURSOR_ID=p.ID
          JOIN PEPTIDE_PROTEIN_MAPPING pp ON pp.PEPTIDE_ID=m.PEPTIDE_ID
          JOIN PROTEIN pr ON pr.ID=pp.PROTEIN_ID WHERE p.DECOY=0 GROUP BY p.ID""", conn)
    ent = labels.set_index('precursor_id').e
    ratio = int(labels.e.sum()) / int((labels.e == 0).sum())
    rows = []
    for arm in ['never', 'always']:
        for seed in range(5):
            directory = os.path.join(layout, f'{dataset}_{arm}_s{seed}')
            files = sorted(glob.glob(directory + '/*.precursors.tsv*'))
            files = [f for f in files if os.path.basename(f) != 'precursors.tsv']
            if len(files) != 1:
                raise ValueError(f'{directory}: expected exactly one precursor TSV, got {files}')
            frame = pd.read_csv(files[0], sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
            frame = frame[(frame.is_decoy == 0) & (frame.score > -1e29)]
            frame = frame.sort_values('score', ascending=False).drop_duplicates('precursor_id')
            frame['e'] = frame.precursor_id.map(ent)
            if frame.e.isna().any():
                raise ValueError('unmapped precursor IDs')
            if not frame.q_value.between(0, 1).all():
                raise ValueError('invalid q values')
            for alpha in ALPHAS:
                accepted = frame[frame.q_value <= alpha]
                n, e = len(accepted), int(accepted.e.sum())
                rows.append(dict(dataset=dataset, arm=arm, seed=seed, alpha=alpha,
                                 N=n, E=e, r=ratio, lower=e/n if n else float('nan'),
                                 combined=e*(1+1/ratio)/n if n else float('nan'),
                                 paired_status='BLOCKED_NO_AUTHENTICATED_PAIR_MAP'))
    return pd.DataFrame(rows)


if __name__ == '__main__':
    lib, layout, dataset, output = sys.argv[1:]
    per = evaluate(lib, layout, dataset)
    per.to_csv(output + '.per_run.tsv', sep='\t', index=False)
    summary = per.groupby(['dataset', 'arm', 'alpha']).agg(
        N_sum=('N', 'sum'), E_sum=('E', 'sum'), N_median=('N', 'median'),
        lower_min=('lower', 'min'), lower_median=('lower', 'median'), lower_max=('lower', 'max'),
        combined_min=('combined', 'min'), combined_median=('combined', 'median'),
        combined_max=('combined', 'max'), r=('r', 'first')).reset_index()
    summary['lower_pooled_descriptive'] = summary.E_sum / summary.N_sum.replace(0, float('nan'))
    summary['combined_pooled_descriptive'] = summary.lower_pooled_descriptive * (1 + 1 / summary.r)
    summary.to_csv(output + '.curve.tsv', sep='\t', index=False)
    print(summary.to_string(index=False))
