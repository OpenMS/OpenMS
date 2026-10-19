"""One ProSE search + external Percolator (3 seeds) for one dataset and arm.

  run_arm.py --bin /opt/openms-pr --params frozen_params.json --mzml x.mzML --fasta db.fasta --out DIR
             [--set Search:key=value ...] [--seeds 1 42 137] [--threads 4]

The frozen #10364 parameter dictionary (jobs.json, native_dedup:<dataset>:default) is passed as
-Search:<key> <value>. ProSE writes the PIN (-out_pin, no in-process rescoring), Percolator 3.09 rescoring
runs per seed with -Y -U, and metrics.json holds per-seed target PSMs / peptidoforms / peptides at q <= 0.01,
entrapment PSMs (all proteins ENTRAP_) and the PIN row statistics.
"""
import argparse, collections, hashlib, json, os, re, subprocess, sys, time
from pathlib import Path

PERC = os.environ.get('PERCOLATOR', '/home/user/bench/percolator')


def sha(p):
    with open(p, 'rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def val(v):
    if isinstance(v, list):
        return [str(x) for x in v]
    return [str(v)]


def strip(pep):
    """Percolator/PIN peptide 'K.PEP(Phospho)TIDE.R' -> residues only."""
    core = pep.split('.', 1)[1].rsplit('.', 1)[0] if pep.count('.') >= 2 else pep
    core = re.sub(r'\([^()]*\)|\[[^\]]*\]', '', core)
    return re.sub(r'[^A-Z]', '', core)


def mods(pep):
    core = pep.split('.', 1)[1].rsplit('.', 1)[0] if pep.count('.') >= 2 else pep
    return sorted(re.findall(r'\(([^()]*)\)|\[([^\]]*)\]', core))


def pin_stats(pin):
    """Rows, targets/decoys, and positional-isomer duplicates (same scan, charge, residues, mod multiset)."""
    rows = collections.Counter(); groups = collections.defaultdict(list)
    with open(pin) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        i_id, i_lab = hdr.index('SpecId'), hdr.index('Label')
        i_pep = hdr.index('Peptide')
        i_scan = hdr.index('ScanNr')
        ch_cols = [(i, int(h[6:])) for i, h in enumerate(hdr) if re.fullmatch(r'[Cc]harge\d+', h)]
        for line in f:
            c = line.rstrip('\n').split('\t')
            if len(c) < len(hdr) - 1 or c[i_lab] not in ('1', '-1'):
                continue
            rows[c[i_lab]] += 1
            z = next((zz for i, zz in ch_cols if c[i] == '1'), 0)
            pep = c[i_pep]
            groups[(c[i_scan], z, c[i_lab], strip(pep), tuple(mods(pep)))].append(pep)
    dup = collections.Counter(); sets = collections.Counter()
    for k, v in groups.items():
        if len(set(v)) > 1:
            sets[k[2]] += 1
            dup[k[2]] += len(v) - 1
    return dict(rows_target=rows['1'], rows_decoy=rows['-1'],
                isomer_sets_target=sets['1'], isomer_sets_decoy=sets['-1'],
                isomer_extra_rows_target=dup['1'], isomer_extra_rows_decoy=dup['-1'],
                spectra=len({k[0] for k in groups}))


def perc_counts(target_tsv, q=0.01):
    psms = 0; entrap = 0; forms = set(); peps = set(); phospho_psms = 0
    with open(target_tsv) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        iq, ip, ipr = hdr.index('q-value'), hdr.index('peptide'), hdr.index('proteinIds')
        for line in f:
            c = line.rstrip('\n').split('\t')
            if float(c[iq]) > q:
                continue
            psms += 1
            prots = [x for x in c[ipr:] if x]
            if prots and all('ENTRAP_' in x for x in prots):
                entrap += 1
            pep = c[ip]
            forms.add(pep.split('.', 1)[1].rsplit('.', 1)[0] if pep.count('.') >= 2 else pep)
            peps.add(strip(pep))
            if 'Phospho' in pep or '79.96' in pep:
                phospho_psms += 1
    return dict(psms=psms, peptidoforms=len(forms), peptides=len(peps), entrapment_psms=entrap,
                phospho_psms=phospho_psms,
                fdp=(2.0 * entrap / psms) if psms else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bin', required=True); ap.add_argument('--params', required=True)
    ap.add_argument('--mzml', required=True); ap.add_argument('--fasta', required=True)
    ap.add_argument('--out', required=True); ap.add_argument('--set', action='append', default=[])
    ap.add_argument('--seeds', nargs='*', type=int, default=[1, 42, 137])
    ap.add_argument('--threads', type=int, default=4); ap.add_argument('--keep-idxml', action='store_true')
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    params = json.loads(Path(a.params).read_text())
    valid = set(json.loads(Path(a.bin, 'valid_params.json').read_text()))
    cmd = [f'{a.bin}/bin/ProSE', '-in', a.mzml, '-database', a.fasta, '-out_pin', str(out / 'input.pin'),
           '-out_idxml', str(out / 'search.idXML'), '-threads', str(a.threads), '-no_summary']
    ignored = []
    for k, v in params.items():
        if k not in valid:
            ignored.append(k); continue
        cmd += ['-Search:' + k] + val(v)
    for s in a.set:
        k, v = s.split('=', 1)
        cmd += [k if k.startswith('-') else '-' + k] + v.split(',')
    env = dict(os.environ, OPENMS_DATA_PATH=f'{a.bin}/share/OpenMS', LD_LIBRARY_PATH=f'{a.bin}/lib',
               OMP_NUM_THREADS=str(a.threads))
    t0 = time.monotonic()
    with open(out / 'search.log', 'w') as log:
        r = subprocess.run(['/usr/bin/time', '-v'] + cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
    wall = time.monotonic() - t0
    if r.returncode:
        sys.exit(f'ProSE failed, see {out}/search.log')
    rss = int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', (out / 'search.log').read_text()).group(1))
    res = dict(cmd=cmd, ignored_parameters=ignored, wall_s=round(wall, 2), max_rss_kb=rss,
               pin_sha256=sha(out / 'input.pin'), idxml_sha256=None, pin=pin_stats(out / 'input.pin'), seeds={})
    # idXML carries the date and file paths; hash it without those lines
    norm = hashlib.sha256()
    with open(out / 'search.idXML', 'rb') as f:
        for line in f:
            if b'date=' in line or b'file_origin' in line or b'SearchParameters' in line and b'db=' in line:
                continue
            norm.update(line)
    res['idxml_sha256_normalized'] = norm.hexdigest()
    if not a.keep_idxml:
        (out / 'search.idXML').unlink()
    for s in a.seeds:
        d = out / f's{s}'; d.mkdir(exist_ok=True)
        with open(d / 'percolator.log', 'w') as log:
            subprocess.run([PERC, '-Y', '-U', '--seed', str(s), '--results-psms', str(d / 'target.tsv'),
                            '--decoy-results-psms', str(d / 'decoy.tsv'), str(out / 'input.pin')],
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        res['seeds'][s] = perc_counts(d / 'target.tsv')
        (d / 'decoy.tsv').unlink()
    (out / 'metrics.json').write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ('wall_s', 'max_rss_kb', 'pin')}), {s: v['psms'] for s, v in res['seeds'].items()})


if __name__ == '__main__':
    main()
