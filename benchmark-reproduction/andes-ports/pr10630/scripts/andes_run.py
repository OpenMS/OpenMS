"""andes 695acf7 (parent of bigbio/andes#116, which made the tie-break the default) with
--fragment-index-intensity-tiebreak off/on, frozen #10364 andes search space (jobs.json; --candidate-index
and --mmap-window-cache-candidates/--fragment-index-slice-da left at the andes defaults of this commit),
on the shuffled entrapment databases; Percolator -Y -U seeds 1/42/137 on andes' native PIN.
usage: andes_run.py <off|on> <dataset>"""
import gzip, hashlib, json, re, subprocess, sys, time
from pathlib import Path
BENCH = Path('/home/user/bench'); DATA = BENCH / 'frozen/broad_benchmark/data'
ANDES = BENCH / 'andes/target/release/andes'; PERC = BENCH / 'THIRDPARTY/Linux/x86_64/Percolator/percolator'
arm, d = sys.argv[1], sys.argv[2]
fasta = {'astral': 'hye', 'hfx': 'hye', 'lumos': 'hy'}[d.split('_')[0]] + '_E_td.fasta'
out = BENCH / 'runs' / f'andes_{arm}_E' / d; out.mkdir(parents=True, exist_ok=True)
if (out / 'summary.json').exists(): sys.exit()
cmd = [str(ANDES), '--spectrum', str(DATA / f'{d}.mzML'), '--database', str(DATA / fasta), '--output-pin', str(out / 'native.pin'),
       '--output-tsv', str(out / 'native.tsv'), '--threads', '4', '--decoy-strategy', 'none', '--decoy-prefix', 'DECOY',
       '--precursor-tol', '20ppm', '--charge', '2..5', '--top-n', '10', '--max-missed-cleavages', '2', '--min-length', '7',
       '--max-length', '40', '--max-mods', '1', '--mods', str(BENCH / 'andes_runs/mods_none.txt'), '--precursor-cal', 'off',
       '--isotope-error=-1..2', '--protocol', 'standard']
if arm == 'on': cmd.append('--fragment-index-intensity-tiebreak')
(out / 'command.json').write_text(json.dumps(cmd, indent=1))
t0 = time.monotonic()
with open(out / 'search.log', 'w') as log: subprocess.run(['/usr/bin/time', '-v'] + cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
res = dict(arm=arm, dataset=d, fasta=fasta, wall_s=round(time.monotonic() - t0, 1),
           max_rss_kb=int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', (out / 'search.log').read_text()).group(1)),
           pin_sha256=hashlib.sha256((out / 'native.pin').read_bytes()).hexdigest(), seed_psms=[])
for s in (1, 42, 137):
    sd = out / f's{s}'; sd.mkdir(exist_ok=True)
    with open(sd / 'percolator.log', 'w') as log:
        subprocess.run([str(PERC), '-Y', '-U', '--seed', str(s), '--results-psms', str(sd / 'target.tsv'),
                        '--decoy-results-psms', str(sd / 'decoy.tsv'), str(out / 'native.pin')], stdout=log, stderr=subprocess.STDOUT, check=True)
    with open(sd / 'target.tsv') as f:
        next(f); res['seed_psms'].append(sum(float(l.split('\t')[2]) <= 0.01 for l in f))
    for f in ('target.tsv', 'decoy.tsv'): subprocess.run(['gzip', '-f', str(sd / f)], check=True)
for f in ('native.pin', 'native.tsv'): subprocess.run(['gzip', '-f', str(out / f)], check=True)
(out / 'summary.json').write_text(json.dumps(res, indent=1)); print(json.dumps(res), flush=True)
