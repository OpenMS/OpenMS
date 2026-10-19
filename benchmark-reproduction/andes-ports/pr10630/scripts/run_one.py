"""One ProSE search of the frozen #10364 suite plus external Percolator (seeds 1/42/137).

usage: run_one.py <build> <arm> <dataset> [--entrap] [--rep N] [--keep-idxml]
  build   : directory /opt/openms-<build> (bin/, lib/, share/) holding the ProSE binary
  arm     : free label; 'int' arms add -Search:scoring:candidate_tie_break intensity
  --entrap: search the shuffled paired-mass entrapment database (<db>_E_td.fasta, r = 1)
  --rep N : timing replicate N (> 0: search only, no Percolator, outputs removed)

Parameters are the frozen native_dedup:<dataset>:default parameters of the #10364 bundle
(jobs.json), mapped onto the current ProSE parameter names; the two removed ones
(precursor:isotopes, snes_enabled) are dropped. Everything else stays at the ProSE default.
Rescoring: -out_pin (no in-process -rescore), Percolator -Y -U --seed s, q <= 0.01.
"""
import hashlib, json, os, re, subprocess, sys, time
from pathlib import Path

BENCH = Path('/home/user/bench')
DATA = BENCH / 'frozen/broad_benchmark/data'
PERCOLATOR = BENCH / 'THIRDPARTY/Linux/x86_64/Percolator/percolator'
SEEDS = (1, 42, 137)
DROPPED = {'precursor:isotopes', 'snes_enabled'}


def sha(p):
    with open(p, 'rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


IDB = BENCH / 'idb/data'


def idb_params():
    """Second dataset (id_bruker, MHC class I, timsTOF Pro): the frozen tims parameters with the immunopeptidomics
    search space (unspecific cleavage, length 8-12, no fixed modification, <= 1 Ox(M))."""
    p, _ = frozen('tims_plasma_30min')
    p = dict(p)
    p.update({'enzyme': 'unspecific cleavage', 'peptide:enzyme_specificity': 'none', 'peptide:min_size': 8, 'peptide:max_size': 12,
              'peptide:missed_cleavages': 0, 'modifications:fixed': [], 'database:chunk_size': CHUNK})
    return p, 'human_td.fasta'


CHUNK = 5000


def frozen(dataset):
    for j in json.load(open(BENCH / 'frozen/jobs.json')):
        if j['stage'] == 'native_dedup' and j['arm'] == 'default' and j['dataset'] == dataset:
            fasta = next(x for x in j['original_command'] if x.endswith('.fasta')).split('/')[-1]
            return j['parameters'], fasta
    raise KeyError(dataset)


def value(v):
    if isinstance(v, list): return [str(x) for x in v]
    return [str(v)]


def main():
    a = sys.argv[1:]
    build, arm, dataset = a[:3]
    entrap = '--entrap' in a
    rep = int(a[a.index('--rep') + 1]) if '--rep' in a else 0
    keep_idxml = '--keep-idxml' in a
    params, fasta = idb_params() if dataset.startswith('idb_') else frozen(dataset)
    data = IDB if dataset.startswith('idb_') else DATA
    if entrap: fasta = fasta.replace('_td.fasta', '_E_td.fasta')
    out = BENCH / 'runs' / (arm + ('_E' if entrap else '')) / dataset / (f'rep{rep}' if rep else '')
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'summary.json').exists():
        print('done already', out); return
    prefix = Path('/opt/openms-' + build)
    cmd = [str(prefix / 'bin/ProSE'), '-in', str(data / (dataset + '.mzML')), '-database', str(data / fasta),
           '-out_pin', str(out / 'search.pin.tsv'), '-out_idxml', str(out / 'search.idXML'), '-threads', '4', '-no_summary']
    for k, v in params.items():
        if k in DROPPED: continue
        cmd += ['-Search:' + k] + (value(v) if v != [] else [])
    if arm.startswith('int'):
        cmd += ['-Search:scoring:candidate_tie_break', 'intensity']
    env = dict(os.environ, OPENMS_DATA_PATH=str(prefix / 'share/OpenMS'), LD_LIBRARY_PATH=str(prefix / 'lib'), OMP_NUM_THREADS='4')
    if build == 'diag': env['OPENMS_TIEBREAK_DIAG'] = str(out / 'tiebreak_diag.tsv')
    (out / 'command.json').write_text(json.dumps(cmd, indent=1))
    t0 = time.monotonic()
    with open(out / 'search.log', 'w') as log:
        subprocess.run(['/usr/bin/time', '-v'] + cmd, stdout=log, stderr=subprocess.STDOUT, env=env, check=True)
    wall = time.monotonic() - t0
    logtxt = (out / 'search.log').read_text()
    rss_kb = int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', logtxt).group(1))
    cpu = sum(float(x) for x in re.findall(r'(?:User|System) time \(seconds\): ([\d.]+)', logtxt))
    res = dict(build=build, arm=arm, dataset=dataset, entrap=entrap, rep=rep, fasta=fasta, fasta_sha256=sha(data / fasta),
               mzml_sha256=sha(data / (dataset + '.mzML')), wall_s=round(wall, 2), cpu_s=round(cpu, 2), max_rss_kb=rss_kb)
    res['pin_sha256'] = sha(out / 'search.pin.tsv')
    res['idxml_sha256'] = sha(out / 'search.idXML')
    # idXML without the date/time attribute and the version string, for identity checks
    txt = (out / 'search.idXML').read_text()
    txt = re.sub(r'date="[^"]*"', '', txt)
    txt = re.sub(r'search_engine_version="[^"]*"', '', txt)  # the build's branch name
    res['idxml_norm_sha256'] = hashlib.sha256(txt.encode()).hexdigest()
    if rep:
        for f in ('search.pin.tsv', 'search.idXML'): (out / f).unlink()
        (out / 'summary.json').write_text(json.dumps(res, indent=1)); print(json.dumps(res)); return
    pin = out / 'search.pin.tsv'
    res['pin_rows'] = sum(1 for _ in open(pin)) - 1
    res['seed_psms'] = []
    for s in SEEDS:
        d = out / f's{s}'; d.mkdir(exist_ok=True)
        with open(d / 'percolator.log', 'w') as log:
            subprocess.run([str(PERCOLATOR), '-Y', '-U', '--seed', str(s), '--results-psms', str(d / 'target.tsv'),
                            '--decoy-results-psms', str(d / 'decoy.tsv'), str(pin)], stdout=log, stderr=subprocess.STDOUT, check=True)
        n = 0
        with open(d / 'target.tsv') as f:
            next(f)
            for line in f:
                if float(line.split('\t')[2]) <= 0.01: n += 1
        res['seed_psms'].append(n)
        for f in ('target.tsv', 'decoy.tsv'):
            subprocess.run(['gzip', '-f', str(d / f)], check=True)
    if not keep_idxml: (out / 'search.idXML').unlink()
    subprocess.run(['gzip', '-f', str(pin)], check=True)
    (out / 'summary.json').write_text(json.dumps(res, indent=1)); print(json.dumps(res), flush=True)


if __name__ == '__main__':
    main()
