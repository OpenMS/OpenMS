"""Tables for PR #10634 from runs/<suite>/<dataset>/<arm>/metrics.json.

  evaluate.py <runs_dir> <out_dir> [--require-complete]

Per file: mean over Percolator seeds 1/42/137 of target PSMs at q <= 0.01; per group: mean over files.
z = delta / sqrt(var_on/3 + var_off/3) over the seeds; a file is up/down when |z| >= 2 (#10364 rule).
Entrapment: N and N_E are seed means summed over a group's files; combined FDP = 2 N_E / N (r = 1).
"""
import json, statistics as st, sys
from pathlib import Path

GROUPS = [('Velos CID', ['velos_125_R1', 'velos_5000_R2', 'velos_25000_R3']),
          ('HF-X HCD', ['hfx_A2', 'hfx_B1', 'hfx_B3']),
          ('Astral HCD', ['astral_A2', 'astral_B1', 'astral_B3']),
          ('Lumos HCD LFQ', ['lumos_lfq_5192', 'lumos_lfq_5194', 'lumos_lfq_5199']),
          ('Lumos CID TMT6plex', ['lumos_tmt_5058', 'lumos_tmt_5059', 'lumos_tmt_5066']),
          ('Exploris 480 TMTpro', ['eclipse_tmtpro_10855', 'eclipse_tmtpro_10858', 'eclipse_tmtpro_10863']),
          ('timsTOF HT plasma', ['tims_plasma_30min', 'tims_plasma_50min'])]
SEEDS = ['1', '42', '137']


def load(runs, d, arm):
    p = runs / d / arm / 'metrics.json'
    return json.loads(p.read_text()) if p.exists() else None


def seeds(m, key):
    return [m['seeds'][s][key] for s in SEEDS]


def z(a, b):
    se = (st.variance(a) / 3 + st.variance(b) / 3) ** 0.5
    d = st.mean(a) - st.mean(b)
    return d / se if se > 0 else (0.0 if d == 0 else float('inf') * (1 if d > 0 else -1))


def main():
    runs, out = Path(sys.argv[1]), Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
    complete = '--require-complete' in sys.argv
    files = [d for _, fs in GROUPS for d in fs]
    pf = ['dataset\tbase_psms\toff_psms\ton_psms\tdelta\tz\tcall\toff_seeds\ton_seeds\toff_peptides\ton_peptides\t'
          'off_rows\ton_rows\tisomer_extra_rows_target\tisomer_extra_rows_decoy\tbase_eq_off_pin\tbase_eq_off_idxml\t'
          'off_wall_s\ton_wall_s']
    per = {}
    for d in files:
        b, o, n = (load(runs, d, a) for a in ('base', 'off', 'on'))
        if None in (o, n):
            if complete: sys.exit(f'missing {d}')
            continue
        so, sn = seeds(o, 'psms'), seeds(n, 'psms')
        zz = z(sn, so); call = 'up' if zz >= 2 else 'down' if zz <= -2 else 'flat'
        per[d] = dict(off=st.mean(so), on=st.mean(sn), base=st.mean(seeds(b, 'psms')) if b else None, call=call,
                      off_pep=st.mean(seeds(o, 'peptides')), on_pep=st.mean(seeds(n, 'peptides')))
        pf.append('\t'.join(map(str, [d, f"{per[d]['base']:.1f}" if b else '', f'{st.mean(so):.1f}', f'{st.mean(sn):.1f}',
                  f'{st.mean(sn) - st.mean(so):+.1f}', f'{zz:+.2f}', call, '/'.join(map(str, so)), '/'.join(map(str, sn)),
                  f"{per[d]['off_pep']:.1f}", f"{per[d]['on_pep']:.1f}",
                  o['pin']['rows_target'] + o['pin']['rows_decoy'], n['pin']['rows_target'] + n['pin']['rows_decoy'],
                  o['pin']['isomer_extra_rows_target'], o['pin']['isomer_extra_rows_decoy'],
                  (b['pin_sha256'] == o['pin_sha256']) if b else '', (b['idxml_sha256_normalized'] == o['idxml_sha256_normalized']) if b else '',
                  o['wall_s'], n['wall_s']])))
    (out / 'per_file.tsv').write_text('\n'.join(pf) + '\n')
    pg = ['group\tn\tbase\toff\ton\tdelta_pct\tup/flat/down\toff_peptides\ton_peptides\tpeptides_delta_pct']
    tot = dict(base=0, off=0, on=0, op=0, np=0)
    for g, fs in GROUPS:
        fs = [d for d in fs if d in per]
        if not fs: continue
        m = {k: st.mean(per[d][k] for d in fs) for k in ('off', 'on', 'off_pep', 'on_pep')}
        bb = [per[d]['base'] for d in fs if per[d]['base'] is not None]
        calls = [per[d]['call'] for d in fs]
        for d in fs:
            tot['off'] += per[d]['off']; tot['on'] += per[d]['on']; tot['base'] += per[d]['base'] or 0
            tot['op'] += per[d]['off_pep']; tot['np'] += per[d]['on_pep']
        pg.append(f"{g}\t{len(fs)}\t{st.mean(bb):.1f}\t{m['off']:.1f}\t{m['on']:.1f}\t{100 * (m['on'] / m['off'] - 1):+.2f}\t"
                  f"{calls.count('up')}/{calls.count('flat')}/{calls.count('down')}\t{m['off_pep']:.1f}\t{m['on_pep']:.1f}\t"
                  f"{100 * (m['on_pep'] / m['off_pep'] - 1):+.2f}")
    calls = [per[d]['call'] for d in per]
    pg.append(f"All files (sum)\t{len(per)}\t{tot['base']:.1f}\t{tot['off']:.1f}\t{tot['on']:.1f}\t"
              f"{100 * (tot['on'] / tot['off'] - 1):+.2f}\t{calls.count('up')}/{calls.count('flat')}/{calls.count('down')}\t"
              f"{tot['op']:.1f}\t{tot['np']:.1f}\t{100 * (tot['np'] / tot['op'] - 1):+.2f}")
    (out / 'per_group.tsv').write_text('\n'.join(pg) + '\n')
    # entrapment
    en = ['group\tarm\tN\tN_E\tcombined_fdp_pct\tper_seed_N_E']
    for g, fs in GROUPS:
        for arm in ('off_entrap', 'on_entrap', 'off_frozenE', 'on_frozenE'):
            ms = [load(runs, d, arm) for d in fs]
            if any(m is None for m in ms): continue
            N = sum(st.mean(seeds(m, 'psms')) for m in ms); E = sum(st.mean(seeds(m, 'entrapment_psms')) for m in ms)
            per_seed = [sum(m['seeds'][s]['entrapment_psms'] for m in ms) for s in SEEDS]
            en.append(f'{g}\t{arm}\t{N:.1f}\t{E:.1f}\t{200 * E / N:.2f}\t{"/".join(map(str, per_seed))}')
    (out / 'entrapment.tsv').write_text('\n'.join(en) + '\n')
    for f in ('per_group.tsv', 'entrapment.tsv'):
        print((out / f).read_text())


if __name__ == '__main__':
    main()
