"""Site accuracy on the MHCquant spike-in (PXD011628): 66 synthetic HLA peptides, each with exactly one heavy V/P/L.

  mhc_site.py <spiked_peptides_labels.tsv> <arm_dir> [<arm_dir> ...]

For each arm dir (run_arm.py output: input.pin, s<seed>/target.tsv) and seed: accepted PSMs (q <= 0.01) whose
residue sequence is a spiked peptide, split by the label they carry: exactly one heavy label at the true position
(correct), exactly one at another position (wrong site), none (light) or several. Also ProSE's own rank-1 row of the
same spectra (best 'score' in the PIN) and, where Percolator accepted another row than ProSE's rank 1, which of the two
had the correct site.
"""
import collections, re, sys
from pathlib import Path

HEAVY = {'V': 6.0138, 'P': 6.0138, 'L': 7.0172}


def parse(pep):
    """'-.AIV[+6.013809]DKVPSV.-' -> ('AIVDKVPSV', [(pos1, residue, delta), ...])"""
    core = pep.split('.', 1)[1].rsplit('.', 1)[0] if pep.count('.') >= 2 else pep
    res, mods = [], []
    for m in re.finditer(r'([A-Z])((?:\[[^\]]*\]|\([^)]*\))*)', core):
        res.append(m.group(1))
        for d in re.findall(r'\[([+-]?[0-9.]+)\]', m.group(2)):
            mods.append((len(res), m.group(1), float(d)))
    return ''.join(res), mods


def classify(pep, truth):
    seq, mods = parse(pep)
    if seq not in truth:
        return None
    pos, aa = truth[seq]
    heavy = [(p, a) for p, a, d in mods if a in HEAVY and abs(d - HEAVY[a]) < 0.01]
    if not heavy:
        return 'light'
    if len(heavy) > 1:
        return 'multi'
    return 'correct' if heavy[0] == (pos, aa) else 'wrong_site'


def main():
    truth = {}
    for l in open(sys.argv[1]).read().splitlines()[1:]:
        c = l.split('\t'); truth[c[0]] = (int(c[3]), c[2])
    out = ['arm\tseed\taccepted\tspiked_psms\tcorrect\twrong_site\tlight\tmulti\tsite_accuracy_pct\t'
           'spiked_peptides_correct\tprose_rank1_correct\tprose_rank1_wrong_site\tswitched\tswitched_percolator_correct\t'
           'switched_prose_correct']
    for d in map(Path, sys.argv[2:]):
        best = {}
        with open(d / 'input.pin') as f:
            h = f.readline().rstrip('\n').split('\t'); i_id, i_sc, i_p = h.index('SpecId'), h.index('score'), h.index('Peptide')
            for line in f:
                c = line.rstrip('\n').split('\t'); s = float(c[i_sc])
                if c[i_id] not in best or s > best[c[i_id]][0]:
                    best[c[i_id]] = (s, c[i_p])
        for t in sorted(d.glob('s*/target.tsv'), key=lambda p: int(p.parent.name[1:])):
            n = 0; cnt = collections.Counter(); peps = set(); r1 = collections.Counter(); sw = collections.Counter()
            with open(t) as f:
                f.readline()
                for line in f:
                    c = line.rstrip('\n').split('\t')
                    if float(c[2]) > 0.01: continue
                    n += 1
                    k = classify(c[4], truth)
                    if k is None: continue
                    cnt[k] += 1
                    if k == 'correct': peps.add(parse(c[4])[0])
                    kr = classify(best[c[0]][1], truth)
                    r1[kr] += 1
                    if best[c[0]][1] != c[4]:
                        sw['n'] += 1; sw['perc'] += k == 'correct'; sw['prose'] += kr == 'correct'
            tot = sum(cnt.values()); lab = cnt['correct'] + cnt['wrong_site']
            out.append('\t'.join(map(str, [d.name, t.parent.name[1:], n, tot, cnt['correct'], cnt['wrong_site'], cnt['light'],
                       cnt['multi'], f"{100 * cnt['correct'] / lab:.2f}" if lab else '', len(peps), r1['correct'],
                       r1['wrong_site'], sw['n'], sw['perc'], sw['prose']])))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
