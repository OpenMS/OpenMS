#!/usr/bin/env python3
"""Paired evaluation of PR #10631 (baseline vs PR) on isobaric quantification outputs.

For each run, tool (IA = IsobaricAnalyzer + PSMs joined by scan, IW = IsobaricWorkflow) and arm:
  features, PSM-quantified features (1% FDR targets), precursor purity present/median/quantiles,
  values changed vs the other arm, reporter scans whose scan window starts above the lowest reporter,
  and for the TKO design (PXD017803, TMTpro16) the knockout-channel interference:
  median over PSMs of the knockout proteins of mean(KO channels)/mean(other channels),
  overall, on scans whose scan window covers 126 Th, and at purity >= 0.9 (as in andes docs/benchmarks).
Arms: base (merge-base), pr (PR head), fix (PR head + in-window fix). "changed" counts values that differ from base.
Runs given as name[:tko][@group]; runs sharing a group are also pooled (KO medians over the pooled PSMs).
usage: evaluate.py <data_dir> <out_tsv> run1[:tko][@group] ...
"""
import sys, csv, statistics, os

ARMS = [a for a in os.environ.get('ARMS', 'base pr fix').split()]

KO = {'P05694': range(0, 4), 'P00815': range(4, 8), 'P07259': range(8, 12)}  # MET6, HIS4, URA2 (126-128N, 128C-130N, 130C-132N)

def read(p):
    return list(csv.DictReader(open(p), delimiter='\t'))

def q(v, f):
    v = sorted(v); return v[min(len(v) - 1, int(f * len(v)))] if v else float('nan')

def ko_protein(acc):
    ids = {a.split('|')[1] if '|' in a else a for a in acc.split(';') if a}
    hits = [k for k in KO if k in ids]
    return hits[0] if len(hits) == 1 and len(ids) == 1 else None

def ko_ratio(row, prot, nch):
    ints = [float(row[f'ch{i}']) for i in range(nch)]
    ko = [ints[i] for i in KO[prot]]
    other = [ints[i] for i in range(12) if i not in KO[prot]]  # the other knockout strains' channels
    other += ints[12:16]  # and the parental strain
    mo = statistics.mean(other)
    return statistics.mean(ko) / mo if mo > 0 else None

def load(data, run, tool, arm, scans, ids):
    rows = read(f'{data}/quant/{arm}/{run}.{tool}.tsv')
    nch = len([k for k in rows[0] if k.startswith('ch')]) if rows else 0
    feat = {}
    for r in rows:
        key = r['id_scan_id'] or r['scan_id']
        if tool == 'IA':
            idr = ids.get(key)
            if idr is None: continue
            r['sequence'], r['accessions'] = idr['sequence'], idr['accessions']
        elif r['target_decoy'] != 'target':
            continue
        r['_pur'] = float(r['precursor_purity']) if r['precursor_purity'] else 0.0  # not written when 0 / unavailable
        sc = scans.get(r['scan_id'], {})
        r['_covers'] = not (sc.get('sw_lo') and float(sc['sw_lo']) > 126.0)
        r['_nch'] = nch
        feat[key] = r
    return feat, len(rows)

def summarize(run, tool, arm, feats, nall, base, tko):
    pur = [r['_pur'] for r in feats.values()]
    shared = set(feats) & set(base)
    rec = dict(run=run, tool=tool, arm=arm, features=nall, psms_quantified=len(feats),
               psms_scan_window_misses_reporters=sum(1 for r in feats.values() if not r['_covers']),
               purity_present=sum(1 for p in pur if p > 0), purity_median=round(q(pur, 0.5), 3),
               purity_q25=round(q(pur, 0.25), 3), purity_ge_0_9=sum(1 for p in pur if p >= 0.9),
               purity_changed_vs_base=sum(1 for k in shared if abs(feats[k]['_pur'] - base[k]['_pur']) > 1e-6))
    if tko:
        def med(sel):
            v = [x for x in (ko_ratio(r, ko_protein(r['accessions']), r['_nch']) for r in sel) if x is not None]
            return (round(statistics.median(v), 3) if v else float('nan')), len(v)
        kos = [r for r in feats.values() if ko_protein(r['accessions'])]
        rec['ko_interference_all'], rec['ko_n_all'] = med(kos)
        cov = [r for r in kos if r['_covers']]
        rec['ko_interference_covering'], rec['ko_n_covering'] = med(cov)
        rec['ko_interference_pur_ge_0_9'], rec['ko_n_pur_ge_0_9'] = med([r for r in cov if r['_pur'] >= 0.9])
    return rec

def main():
    data, out = sys.argv[1], sys.argv[2]
    rows_out = []
    groups = {}
    for spec in sys.argv[3:]:
        spec, _, group = spec.partition('@')
        run, _, tko = spec.partition(':')
        sub = 'pxd017803' if os.path.exists(f'{data}/pxd017803/{run}.scans.tsv') else 'arch'
        scans = {r['native_id']: r for r in read(f'{data}/{sub}/{run}.scans.tsv')}
        ids = {r['spectrum_reference']: r for r in read(f'{data}/search/{run}.ids.tsv') if r['target_decoy'] == 'target'}
        for tool in ('IA', 'IW'):
            loaded = {arm: load(data, run, tool, arm, scans, ids) for arm in ARMS}
            for arm in ARMS:
                feats, nall = loaded[arm]
                rows_out.append(summarize(run, tool, arm, feats, nall, loaded['base'][0], tko))
                if group:
                    g = groups.setdefault((group, tool, arm), [{}, 0, {}, tko])
                    g[0].update({run + '|' + k: v for k, v in feats.items()}); g[1] += nall
                    g[2].update({run + '|' + k: v for k, v in loaded['base'][0].items()})
    for (group, tool, arm), (feats, nall, base, tko) in groups.items():
        rows_out.append(summarize(group, tool, arm, feats, nall, base, tko))
    keys = []
    for r in rows_out:
        for k in r:
            if k not in keys: keys.append(k)
    w = csv.DictWriter(open(out, 'w', newline=''), fieldnames=keys, delimiter='\t', restval='')
    w.writeheader(); w.writerows(rows_out)

if __name__ == '__main__':
    main()
