#!/usr/bin/env python3
"""Score a --precursor-mono-dump against the pGlyco2 reference for one run.

    score_dump.py <dump.tsv> <truth.tsv.gz> <run> [testset.tsv] [--min-fit F --min-gain G --min-snr S --backoff B]

For every reference scan, the TRUE offset is the integer k such that
    recorded_neutral - (peptide + 57.02146*C + 15.99491*ox + glycan) ~= k * 1.00335
minimised over ox in 0..#M (pGlyco2 searched Carbamidomethyl-C fixed, Oxidation-M
variable). Scans whose residual is not within 0.15 Da of an integer are dropped
("unresolved").

Reports the confusion of best_shift vs k_true, and, for a given decision rule
(re-applied here from the dump's fit columns so thresholds can be swept without
re-running andes), how many reference scans end up REACHABLE by the default 0..2
sweep from the corrected precursor: k_true - applied in {0,1,2}.
"""
import csv, gzip, sys, argparse
from collections import Counter, defaultdict

AA = {'G':57.021464,'A':71.037114,'S':87.032028,'P':97.052764,'V':99.068414,'T':101.047679,
      'C':103.009185,'L':113.084064,'I':113.084064,'N':114.042927,'D':115.026943,'Q':128.058578,
      'K':128.094963,'E':129.042593,'M':131.040485,'H':137.058912,'F':147.068414,'R':156.101111,
      'Y':163.06332,'W':186.079313}
H2O = 18.010565; PROTON = 1.007276; ISO = 1.003355; CAM = 57.02146; OX = 15.99491
MONO = {'HexNAc':203.079373,'Hex':162.052824,'Fuc':146.057909,'NeuAc':291.095417,'NeuGc':307.090331}

def glycan_mass(g):
    import re
    m = 0.0
    for name, n in re.findall(r'(HexNAc|Hex|Fuc|NeuAc|NeuGc)(\d+)', g):
        m += MONO[name]*int(n)
    return m

def pep_mass(p):
    return sum(AA[c] for c in p) + H2O + CAM*p.count('C')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dump'); ap.add_argument('truth'); ap.add_argument('run')
    ap.add_argument('testset', nargs='?')
    ap.add_argument('--min-fit', type=float, default=0.90)
    ap.add_argument('--min-gain', type=float, default=0.15)
    ap.add_argument('--min-snr', type=float, default=3.0)
    ap.add_argument('--backoff', type=int, default=1)
    ap.add_argument('--sweep', type=int, default=2, help='upper isotope-error of the search sweep')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args()

    dump = {}
    with open(a.dump) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        fitcols = [c for c in r.fieldnames if c.startswith('fit') and c[3:].isdigit()]
        for row in r:
            dump[int(row['scan'])] = row
    truth = {}
    with gzip.open(a.truth, 'rt') as fh:
        for row in csv.DictReader((l for l in fh if not l.startswith('#')), delimiter='\t'):
            if row['run'] == a.run:
                truth[int(row['scan'])] = row
    testset = {}
    if a.testset:
        with open(a.testset) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                testset[int(row['scan'])] = row

    def decide(row):
        """Re-apply the decision rule on the dump's per-k fits."""
        if row['shift'] == 'NA':
            return None, None
        fits = [float(row[c]) for c in fitcols]
        best = max(range(len(fits)), key=lambda k: (fits[k], -k))
        # snr in the dump is that of the SEARCHED hypothesis; for re-deciding we need
        # the best hypothesis' SNR, which the dump gives only when it was applied.
        # Use the dump's own decision when thresholds are the defaults; otherwise the
        # SNR gate is approximated by the dump's snr column when best was applied.
        snr_ok = True
        if best > 0 and int(row['best_shift']) == best and int(row['shift']) > 0:
            snr_ok = float(row['snr']) >= a.min_snr
        elif best > 0 and int(row['shift']) == 0 and int(row['best_shift']) == best:
            # not applied in the run: either fit/gain/snr failed. We cannot see the
            # best hypothesis' SNR, so trust the run's SNR gate implicitly only when
            # fit and gain would have passed (then SNR was the reason).
            pass
        apply = best > 0 and fits[best] >= a.min_fit and fits[best]-fits[0] >= a.min_gain and snr_ok
        applied = max(best - a.backoff, 0) if apply else 0
        return best, applied

    conf = Counter(); reach = Counter(); tot_by_k = Counter(); unresolved = 0; nolink = 0
    applied_hist = Counter(); false_shift = Counter()
    per_scan = {}
    for scan, t in truth.items():
        d = dump.get(scan)
        if d is None or d['shift'] == 'NA':
            nolink += 1; continue
        z = int(d['charge']); rec = float(d['recorded_mz'])
        neutral = (rec - PROTON) * z
        pep = t['peptide']; nM = pep.count('M')
        base = pep_mass(pep) + glycan_mass(t['glycan'])
        best_res = None
        for ox in range(nM+1):
            delta = neutral - (base + OX*ox)
            k = round(delta/ISO); res = abs(delta - k*ISO)
            if best_res is None or res < best_res[0]:
                best_res = (res, k, ox)
        res, ktrue, ox = best_res
        if res > 0.15 or ktrue < -1 or ktrue > 8:
            unresolved += 1; continue
        best, applied = decide(d)
        tot_by_k[ktrue] += 1
        conf[(ktrue, best)] += 1
        applied_hist[(ktrue, applied)] += 1
        reachable = 0 <= (ktrue - applied) <= a.sweep
        reach[(ktrue, reachable)] += 1
        if applied > 0 and not reachable:
            false_shift[ktrue] += 1
        per_scan[scan] = (ktrue, best, applied, reachable, d)

    print(f"run {a.run}: {len(truth)} reference scans; {nolink} without MS1 link/charge; {unresolved} unresolved mass")
    ks = sorted(tot_by_k)
    print("\nbest_shift (columns) vs true offset (rows):")
    cols = sorted({b for (_, b) in conf})
    print("k_true\tn\t" + "\t".join(f"b={c}" for c in cols))
    for k in ks:
        print(f"{k}\t{tot_by_k[k]}\t" + "\t".join(str(conf[(k, c)]) for c in cols))
    print(f"\napplied shift (rule: min_fit {a.min_fit}, min_gain {a.min_gain}, backoff {a.backoff}); reachable = k_true-applied in 0..{a.sweep}")
    acols = sorted({b for (_, b) in applied_hist})
    print("k_true\tn\treachable\tlost\t" + "\t".join(f"a={c}" for c in acols))
    tr = tl = 0
    for k in ks:
        r_ = reach[(k, True)]; l_ = reach[(k, False)]; tr += r_; tl += l_
        print(f"{k}\t{tot_by_k[k]}\t{r_}\t{l_}\t" + "\t".join(str(applied_hist[(k, c)]) for c in acols))
    print(f"TOTAL\t{tr+tl}\t{tr}\t{tl}")
    base_reach = sum(tot_by_k[k] for k in ks if 0 <= k <= a.sweep)
    print(f"reachable without correction (k_true in 0..{a.sweep}): {base_reach}; with: {tr}; net {tr-base_reach:+d}")
    print(f"wrongly shifted (applied>0 and truth no longer reachable): {sum(false_shift.values())} {dict(false_shift)}")

    if testset:
        print("\n84-scan test set (bigbio/andes#64):")
        print("scan\tz\tissue_offset\tk_true\tbest\tapplied\treachable\tfit0\tfit_best\tsnr")
        hit = 0; n = 0
        for scan in sorted(testset):
            ts = testset[scan]
            if scan in per_scan:
                ktrue, best, applied, reachable, d = per_scan[scan]
                n += 1; hit += reachable
                print(f"{scan}\t{d['charge']}\t{ts['offset_Da']}\t{ktrue}\t{best}\t{applied}\t{int(reachable)}\t{d['fit_recorded']}\t{d['fit_best']}\t{d['snr']}")
            else:
                d = dump.get(scan)
                print(f"{scan}\t{ts['charge']}\t{ts['offset_Da']}\t?\t{d['best_shift'] if d else '-'}\t{d['shift'] if d else '-'}\t?\t{d['fit_recorded'] if d else '-'}\t{d['fit_best'] if d else '-'}\t-")
        print(f"test set: {hit} of {n} resolvable targets reachable after correction")

    # Overall run-level histogram of applied shifts (all MS2, not only reference)
    allhist = Counter()
    for scan, d in dump.items():
        if d['shift'] == 'NA': continue
        best, applied = decide(d)
        allhist[applied] += 1
    print("\nall MS2: applied-shift histogram under this rule:", dict(sorted(allhist.items())))

if __name__ == '__main__':
    main()
