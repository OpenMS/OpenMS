"""Frozen descriptive full-library ID/FDP/resource evaluation.
Usage: eval_full.py SEARCH.pqp RUNS RESULTS_DIRECTORY
Missing/failed runs remain failures, never substituted or dropped from success counts.
"""
from pathlib import Path
import re
import sqlite3
import sys
import pandas as pd

ALPHAS = [.001,.0025,.005,.0075,.01,.015,.02,.03,.05]


def evaluate(library, runs):
    with sqlite3.connect(Path(library).resolve().as_uri() + '?mode=ro', uri=True) as conn:
        labels = pd.read_sql_query("""SELECT p.ID precursor_id,
        min(CASE WHEN pr.PROTEIN_ACCESSION LIKE '%_p_target%' THEN 1 ELSE 0 END) e
        FROM PRECURSOR p JOIN PRECURSOR_PEPTIDE_MAPPING m ON m.PRECURSOR_ID=p.ID
        JOIN PEPTIDE_PROTEIN_MAPPING pp ON pp.PEPTIDE_ID=m.PEPTIDE_ID
        JOIN PROTEIN pr ON pr.ID=pp.PROTEIN_ID WHERE p.DECOY=0 GROUP BY p.ID""", conn)
    ent = labels.set_index('precursor_id').e
    ratio = int(labels.e.sum()) / int((labels.e == 0).sum())
    counts, status = [], []
    for raw in ['01','02','03']:
        for seed in range(5):
            directory = Path(runs) / f'raw{raw}_s{seed}'
            exit_file = directory / 'exit_code.txt'
            state = dict(raw=raw, seed=seed, exit_code=exit_file.read_text().strip() if exit_file.exists() else 'MISSING')
            time_file = directory / 'time.txt'
            text = time_file.read_text() if time_file.exists() else ''
            rss = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', text)
            state['max_rss_kib'] = int(rss[1]) if rss else None
            wall = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(.+)', text)
            state['wall_time'] = wall[1].strip() if wall else None
            state['within_32GiB_observed_rss'] = state['max_rss_kib'] <= 32*1024**2 if rss else None
            file = directory / f'LFQ_Astral_DIA_15min_50ng_Human_{raw}.precursors.tsv'
            state['status'] = 'SUCCESS' if state['exit_code'] == '0' and file.exists() else 'FAILED_OR_MISSING'
            status.append(state)
            if state['status'] != 'SUCCESS':
                continue
            frame = pd.read_csv(file, sep='\t')
            frame = frame[(frame.is_decoy == 0) & (frame.score > -1e29)]
            if frame.precursor_id.duplicated().any():
                raise ValueError('duplicate target winners in single-run precursor output')
            frame['e'] = frame.precursor_id.map(ent)
            if frame.e.isna().any() or not frame.q_value.between(0,1).all():
                raise ValueError('invalid precursor IDs or q values')
            for alpha in ALPHAS:
                accepted = frame[frame.q_value <= alpha]
                n, e = len(accepted), int(accepted.e.sum())
                human = accepted[accepted.e == 0]
                counts.append(dict(raw=raw, seed=seed, alpha=alpha, N=n, E=e, human_precursors=n-e,
                                   human_modified_peptides=human.sequence.nunique(), r=ratio,
                                   lower=e/n if n else None, combined=e*(1+1/ratio)/n if n else None))
    return pd.DataFrame(counts), pd.DataFrame(status)


if __name__ == '__main__':
    library, runs, output = sys.argv[1:]
    root = Path(output)
    root.mkdir(parents=True, exist_ok=True)
    counts, status = evaluate(library, runs)
    status.to_csv(root/'run_status.tsv', sep='\t', index=False)
    counts.to_csv(root/'calibration_per_run.tsv', sep='\t', index=False)
    if not counts.empty:
        summary = counts.groupby(['raw','alpha']).agg(successful_seeds=('seed','nunique'),
          N_median=('N','median'), human_median=('human_precursors','median'), E_sum=('E','sum'), N_sum=('N','sum'),
          lower_median=('lower','median'), lower_min=('lower','min'), lower_max=('lower','max'),
          combined_median=('combined','median'), combined_min=('combined','min'), combined_max=('combined','max')).reset_index()
        summary.to_csv(root/'calibration_curve.tsv', sep='\t', index=False)
    print(status.to_string(index=False))
    print('Descriptive only. Seed spread is not an independent-replicate confidence interval.')
