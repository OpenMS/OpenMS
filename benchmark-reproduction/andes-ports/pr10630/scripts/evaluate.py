"""Evaluate the PR #10630 benchmark (scoring:candidate_tie_break none vs intensity) on the frozen #10364 suite.

Reads runs/<arm>/<dataset>/summary.json, s<seed>/target.tsv.gz and search.pin.tsv.gz; writes the tables into
the output folder (default: results/). Arms: base (merge-base build), none / int (PR build, option off / on),
none_E / int_E (PR build, shuffled paired-mass entrapment databases, r = 1).

- Yield: target PSMs at Percolator q <= 0.01, mean over seeds 1/42/137 per file, then over files per group.
  z = (mean_int - mean_none) / sqrt(sd_int^2/3 + sd_none^2/3); a file is up/down when |z| >= 2 (#10364 rule).
- Entrapment: shuffled (r = 1) combined FDP = 2 * sum N_E / sum N (seed means summed over a group's files);
  Velos natural entrapment (UPS1 in yeast vs all human proteins, as in pr-prototypes-20261001/entrapment_velos.py)
  combined FDP = N_E (1 + 1/r) / N, summed over 3 files x 3 seeds.
- Top-1 changes: the rank-1 PSM (highest 'score' = ln hyperscore) per spectrum in the PIN, none vs int.
"""
import csv, gzip, json, math, re, statistics, sys
from collections import defaultdict
from pathlib import Path

BENCH = Path('/home/user/bench')
RUNS = BENCH / 'runs'
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else BENCH / 'results'
OUT.mkdir(parents=True, exist_ok=True)
DATA = BENCH / 'frozen/broad_benchmark/data'
SEEDS = (1, 42, 137)
GROUPS = [('velos', 'Velos CID'), ('hfx', 'HF-X HCD'), ('astral', 'Astral HCD'), ('lumos_lfq', 'Lumos HCD LFQ'),
          ('lumos_tmt', 'Lumos CID TMT6plex'), ('eclipse', 'Exploris 480 TMTpro'), ('tims', 'timsTOF HT')]
DATASETS = list(json.load(open(BENCH / 'frozen/inputs.json')))


def group(d):
    return next(name for p, name in GROUPS if d.startswith(p))


def summary(arm, d):
    p = RUNS / arm / d / 'summary.json'
    return json.load(open(p)) if p.exists() else None


def accepted(arm, d, seed):
    """scan -> (peptide, proteins) of the target PSMs at q <= 0.01"""
    out = {}
    with gzip.open(RUNS / arm / d / f's{seed}' / 'target.tsv.gz', 'rt') as f:
        next(f)
        for line in f:
            x = line.rstrip('\n').split('\t')
            if float(x[2]) <= 0.01: out[x[0]] = (x[4], [p for p in x[5:] if p])
    return out


def write(name, rows):
    if not rows: return
    with open(OUT / name, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)


def zscore(a, b):
    se = math.sqrt(statistics.variance(a) / len(a) + statistics.variance(b) / len(b))
    d = statistics.mean(b) - statistics.mean(a)
    return d / se if se > 0 else (0.0 if d == 0 else math.copysign(math.inf, d))


# ---------- yields ----------
def yields(ref, arm, tag):
    per_file, per_group = [], []
    gm = defaultdict(lambda: defaultdict(list))
    for d in DATASETS:
        a, b = summary(ref, d), summary(arm, d)
        if not a or not b: continue
        z = zscore(a['seed_psms'], b['seed_psms'])
        ma, mb = statistics.mean(a['seed_psms']), statistics.mean(b['seed_psms'])
        per_file.append(dict(dataset=d, group=group(d), ref=ref, arm=arm, ref_seeds='/'.join(map(str, a['seed_psms'])),
                             arm_seeds='/'.join(map(str, b['seed_psms'])), ref_mean=round(ma, 1), arm_mean=round(mb, 1),
                             delta_pct=round(100 * (mb - ma) / ma, 2), z=round(z, 2),
                             call='up' if z >= 2 else 'down' if z <= -2 else 'flat'))
        gm[group(d)]['a'].append(ma); gm[group(d)]['b'].append(mb)
    for _, g in GROUPS:
        if g not in gm: continue
        ma, mb = statistics.mean(gm[g]['a']), statistics.mean(gm[g]['b'])
        rows = [r for r in per_file if r['group'] == g]
        per_group.append(dict(group=g, n=len(gm[g]['a']), ref=ref, arm=arm, ref_mean=round(ma, 1), arm_mean=round(mb, 1),
                              delta_pct=round(100 * (mb - ma) / ma, 2),
                              up_flat_down='/'.join(str(sum(r['call'] == c for r in rows)) for c in ('up', 'flat', 'down'))))
    ta = sum(r['ref_mean'] for r in per_file); tb = sum(r['arm_mean'] for r in per_file)
    per_group.append(dict(group='All files (sum)', n=len(per_file), ref=ref, arm=arm, ref_mean=round(ta, 1), arm_mean=round(tb, 1),
                          delta_pct=round(100 * (tb - ta) / ta, 2) if ta else '',
                          up_flat_down='/'.join(str(sum(r['call'] == c for r in per_file)) for c in ('up', 'flat', 'down'))))
    write(f'per_file_{tag}.tsv', per_file); write(f'per_group_{tag}.tsv', per_group)
    return per_file, per_group


# ---------- entrapment ----------
UPS1 = """P00915 P00918 P01031 P69905 P68871 P41159 P02768 P0CG48 P0CG47 P62987 P62979 P04040 P00167 P01133 P02144 P15559
P62937 Q06830 P63165 P00709 P06730 P12081 P61626 Q15843 P02753 P16083 P63279 P01008 P61769 P55957 O76070 P08263 P01344
P01127 P10599 P99999 P06396 P09211 P01112 P01579 P02787 P02788 P05413 P10636 P10145 P02741 O00762 P01375 P51965 P08758""".split()


def read_fasta(p):
    entries, hdr, seq = [], None, []
    for line in open(p):
        if line.startswith('>'):
            if hdr: entries.append((hdr, ''.join(seq)))
            hdr, seq = line[1:].strip(), []
        else: seq.append(line.strip())
    if hdr: entries.append((hdr, ''.join(seq)))
    return [(h, s) for h, s in entries if not h.startswith('DECOY_')]


def digest(s):
    sites = [0] + [i + 1 for i in range(len(s) - 1) if s[i] in 'KR' and s[i + 1] != 'P'] + [len(s)]
    out = set()
    for i in range(len(sites) - 1):
        for j in range(i + 1, min(i + 4, len(sites))):
            p = s[sites[i]:sites[j]]
            if 7 <= len(p) <= 40: out.add(p)
    return out


def stripped(peptide):
    plain = re.sub(r'\[[^\]]*\]|\([^)]*\)', '', peptide)
    core = plain.split('.')[1] if plain.count('.') == 2 else plain
    return re.sub(r'[^A-Z]', '', core).replace('I', 'L')


def velos_entrapment(arms):
    entries = read_fasta(DATA / 'hy_td.fasta')
    sample, ups, human = [], [], []
    sample_peps, ent_peps = set(), set()
    for h, s in entries:
        acc = h.split('|')[1] if '|' in h else h.split()[0]
        s2 = s.replace('I', 'L'); first = h.split()[0]
        if 'Cont_' in h or '_YEAST' in first: sample.append(s2)
        elif acc in UPS1: ups.append(s2)
        elif '_HUMAN' in first: human.append(s2)
        else: sample.append(s2)
        if '_HUMAN' in first and 'Cont_' not in h and acc not in UPS1: ent_peps |= digest(s2)
        else: sample_peps |= digest(s2)
    r = len(ent_peps - sample_peps) / len(sample_peps)
    db = {k: '\n'.join(v) for k, v in dict(sample=sample, ups1=ups, human=human).items()}
    cache = {}

    def cat(p):
        if p not in cache:
            cache[p] = 'sample' if p in db['sample'] else 'ups1' if p in db['ups1'] else 'entrapment' if p in db['human'] else 'unknown'
        return cache[p]
    rows = []
    acc = {(arm, d, s): accepted(arm, d, s) for arm in arms for d in DATASETS if d.startswith('velos') and summary(arm, d) for s in SEEDS}
    for arm in arms:
        n = ne = gained = gained_e = 0
        for (a, d, s), psms in acc.items():
            if a != arm: continue
            for scan, (pep, _) in psms.items():
                c = cat(stripped(pep)); n += 1; ne += c == 'entrapment'
            if arm != arms[0]:
                ref = acc[(arms[0], d, s)]
                for scan, (pep, _) in psms.items():
                    if scan not in ref: gained += 1; gained_e += cat(stripped(pep)) == 'entrapment'
        if n:
            rows.append(dict(group='Velos CID (natural, r=%.3f)' % r, arm=arm, psms_sum_3files_3seeds=n, entrapment=ne,
                             fdp_combined_pct=round(100 * ne * (1 + 1 / r) / n, 3), gained_vs_ref=gained, gained_entrapment=gained_e))
    return rows


def shuffled_entrapment(arms):
    rows, per_file = [], []
    for _, g in GROUPS:
        if g == 'Velos CID': continue
        for arm in arms:
            N = NE = 0.0; nf = 0
            gained = gained_e = 0
            for d in DATASETS:
                if group(d) != g or not summary(arm, d): continue
                n, ne = [], []
                for s in SEEDS:
                    psms = accepted(arm, d, s)
                    n.append(len(psms)); ne.append(sum(all(p.startswith('ENTRAP_') for p in prot) for _, prot in psms.values()))
                    if arm != arms[0] and summary(arms[0], d):
                        ref = accepted(arms[0], d, s)
                        for scan, (_, prot) in psms.items():
                            if scan not in ref: gained += 1; gained_e += all(p.startswith('ENTRAP_') for p in prot)
                per_file.append(dict(dataset=d, arm=arm, psms=round(statistics.mean(n), 1), entrapment=round(statistics.mean(ne), 1),
                                     fdp_combined_pct=round(200 * sum(ne) / sum(n), 2), entrapment_seeds='/'.join(map(str, ne))))
                N += statistics.mean(n); NE += statistics.mean(ne); nf += 1
            if nf:
                rows.append(dict(group=g, arm=arm, files=nf, psms=round(N, 1), entrapment=round(NE, 1),
                                 fdp_combined_pct=round(200 * NE / N, 2), gained_vs_ref=gained, gained_entrapment=gained_e,
                                 gained_fdp_combined_pct=round(200 * gained_e / gained, 1) if gained else ''))
    return rows, per_file


# ---------- top-1 changes in the PIN ----------
def top1(arm, d):
    best = {}
    with gzip.open(RUNS / arm / d / 'search.pin.tsv.gz', 'rt') as f:
        hdr = next(f).rstrip('\n').split('\t')
        i_score, i_pep, i_lab = hdr.index('score'), hdr.index('Peptide'), hdr.index('Label')
        for line in f:
            x = line.rstrip('\n').split('\t')
            sc = float(x[i_score])
            cur = best.get(x[0])
            if cur is None or sc > cur[0]: best[x[0]] = (sc, x[i_pep], x[i_lab])
    return best


def top1_changes():
    rows = []
    for d in DATASETS:
        if not summary('none', d) or not summary('int', d): continue
        a, b = top1('none', d), top1('int', d)
        both = set(a) & set(b)
        changed = [k for k in both if a[k][1] != b[k][1]]
        rows.append(dict(dataset=d, group=group(d), spectra_none=len(a), spectra_int=len(b),
                         only_none=len(set(a) - set(b)), only_int=len(set(b) - set(a)), top1_changed=len(changed),
                         int_higher_score=sum(b[k][0] > a[k][0] for k in changed),
                         int_lower_score=sum(b[k][0] < a[k][0] for k in changed),
                         equal_score=sum(b[k][0] == a[k][0] for k in changed),
                         changed_target_to_decoy=sum(a[k][2] == '1' and b[k][2] == '-1' for k in changed),
                         changed_decoy_to_target=sum(a[k][2] == '-1' and b[k][2] == '1' for k in changed)))
    write('top1_changes.tsv', rows)
    return rows


def accepted_changes():
    rows = []
    for d in DATASETS:
        if not summary('none', d) or not summary('int', d): continue
        g = l = ch = 0
        for s in SEEDS:
            a, b = accepted('none', d, s), accepted('int', d, s)
            g += len(set(b) - set(a)); l += len(set(a) - set(b))
            ch += sum(stripped(a[k][0]) != stripped(b[k][0]) for k in set(a) & set(b))
        rows.append(dict(dataset=d, group=group(d), gained_scans=round(g / 3, 1), lost_scans=round(l / 3, 1),
                         same_scan_other_peptide=round(ch / 3, 1)))
    write('accepted_changes.tsv', rows)
    return rows


def identity():
    """base (merge-base build) vs none (PR build, option off): PIN of the main run, idXML of the timing replicates
    (hash without the date attributes and the search_engine_version string, which carries the branch name)."""
    rows = []
    for d in DATASETS:
        a, b = summary('base', d), summary('none', d)
        if not a or not b: continue
        reps = {arm: [json.load(open(p)) for p in sorted((RUNS / arm / d).glob('rep*/summary.json'))] for arm in ('base', 'none', 'int')}
        idx = {arm: {r['idxml_norm_sha256'] for r in v} for arm, v in reps.items()}
        pins = {arm: {r['pin_sha256'] for r in v} for arm, v in reps.items()}
        rows.append(dict(dataset=d, pin_identical=a['pin_sha256'] == b['pin_sha256'],
                         idxml_identical=(len(idx['base']) == 1 and idx['base'] == idx['none']) if reps['base'] and reps['none'] else '',
                         reps_deterministic=all(len(pins[x]) == 1 for x in pins if reps[x]),
                         seeds_base='/'.join(map(str, a['seed_psms'])), seeds_none='/'.join(map(str, b['seed_psms']))))
    write('identity_base_vs_none.tsv', rows)
    return rows


def resources():
    rows = []
    for d in DATASETS:
        for arm in ('base', 'none', 'int'):
            reps = [json.load(open(p)) for p in sorted((RUNS / arm / d).glob('rep*/summary.json'))]
            if not reps: continue
            rows.append(dict(dataset=d, arm=arm, reps=len(reps), wall_s_median=statistics.median(r['wall_s'] for r in reps),
                             cpu_s_median=statistics.median(r['cpu_s'] for r in reps),
                             max_rss_mib_median=round(statistics.median(r['max_rss_kb'] for r in reps) / 1024, 1)))
    if rows:
        for arm in ('base', 'none', 'int'):
            r = [x for x in rows if x['arm'] == arm]
            rows.append(dict(dataset='SUM / median', arm=arm, reps=min(x['reps'] for x in r), wall_s_median=round(sum(x['wall_s_median'] for x in r), 1),
                             cpu_s_median=round(sum(x['cpu_s_median'] for x in r), 1),
                             max_rss_mib_median=statistics.median(x['max_rss_mib_median'] for x in r)))
    write('resources.tsv', rows)
    return rows


def diag():
    """OPENMS_TIEBREAK_DIAG output (diagnostic build): one row per fragment-index query."""
    rows = []
    for d in DATASETS:
        p = RUNS / 'int_diag' / d / 'tiebreak_diag.tsv'
        if not p.exists(): continue
        n = over = tie = sw = swp = q = 0; tie_sizes = []
        for line in open(p):
            x = line.rstrip('\n').split('\t')
            q += 1
            nel, bt, v, at_v, slots, s, sp = int(x[1]), int(x[2]), int(x[3]), int(x[4]), int(x[5]), int(x[6]), int(x[7])
            over += nel > 50; tie += bt; sw += s > 0; swp += sp > 0
            if bt: tie_sizes.append(at_v)
        rows.append(dict(dataset=d, group=group(d), queries=q, over_cap=over, boundary_tie=tie,
                         boundary_tie_pct=round(100 * tie / q, 1) if q else '',
                         median_tie_size=statistics.median(tie_sizes) if tie_sizes else '',
                         candidate_set_changed=sw, peptide_set_changed=swp))
    write('vote_ties.tsv', rows)
    return rows


if __name__ == '__main__':
    yields('none', 'int', 'int_vs_none')
    yields('base', 'none', 'none_vs_base')
    identity()
    ve = velos_entrapment(['none', 'int', 'intiso'])
    # exploratory variant (local patch variant_intensity_after_isotope.patch, not in the PR): intensity key after the isotope keys
    yields('none', 'intiso', 'intiso_vs_none')
    yields('none_E', 'intiso_E', 'intiso_vs_none_entrapment_db')
    ie, ief = shuffled_entrapment(['none_E', 'intiso_E'])
    yields('none_E', 'int_E', 'int_vs_none_entrapment_db')
    se, sef = shuffled_entrapment(['none_E', 'int_E'])
    # like-for-like andes reference: andes 695acf7, --fragment-index-intensity-tiebreak off/on, same inputs
    yields('andes_off_E', 'andes_on_E', 'andes_on_vs_off_entrapment_db')
    ae, aef = shuffled_entrapment(['andes_off_E', 'andes_on_E'])
    write('entrapment_velos_natural.tsv', ve); write('entrapment_shuffled.tsv', se + [r for r in ie if r['arm'] != 'none_E'] + ae); write('entrapment_shuffled_per_file.tsv', sef + [r for r in ief if r['arm'] != 'none_E'] + aef)
    top1_changes(); accepted_changes(); resources(); diag()
    for f in sorted(OUT.glob('*.tsv')):
        print('==', f.name); print(f.read_text())
