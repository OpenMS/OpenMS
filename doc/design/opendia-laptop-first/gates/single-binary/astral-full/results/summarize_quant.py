"""Compact descriptive summary of the frozen eval_quant JSONs (drops per-feature rows, adds CV medians)."""
import glob
import json
import os
import statistics as st

out = {}
for path in sorted(glob.glob('results/quant_*.json')):
    with open(path) as f:
        d = json.load(f)
    name = os.path.basename(path)[:-5]
    s = dict(q_cutoff=d['q_cutoff'], population=d['population'],
             audit=[{k: v for k, v in a.items() if k != 'path'} for a in d['audit']])
    for block in ('across_raws', 'within_raw_seeds'):
        s[block] = {}
        for level, groups in d[block].items():
            s[block][level] = {}
            for key, g in groups.items():
                rows = g.pop('features')
                cv = [r['cv_raw'] for r in rows if r['cv_raw'] is not None]
                cvn = [r['cv_median_normalized'] for r in rows if r['cv_median_normalized'] is not None]
                g['cv_raw_median'] = st.median(cv) if cv else None
                g['cv_median_normalized_median'] = st.median(cvn) if cvn else None
                g['cv_n'] = len(cv)
                g['pairs'] = [{k: p[k] for k in ('a', 'b', 'jaccard', 'positive_pairs', 'log2_pearson')} for p in g['pairs']]
                s[block][level][key] = g
    out[name] = s
    print(name, flush=True)
with open('diagnostics/quant_summary.json', 'w') as f:
    json.dump(out, f, indent=1)
