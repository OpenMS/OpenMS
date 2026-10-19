"""Interleaved paired runs: for rep in 1..R, for arm in arms: one ProSE search (run1.sh, /usr/bin/time -v).
usage: bench.py <experiments.json> <exp-name> <reps>
experiments.json: {exp: {"db": ..., "in": [...], "threads": N, "arms": {arm: {"prefix": ..., "ini": ..., "args": [...]}}}}
Appends one row per run to results/<exp>.tsv: exp arm rep threads wall_s maxrss_kb idxml_sha256 pin_sha256 status"""
import hashlib, json, re, subprocess, sys
from pathlib import Path
B = Path('/home/user/bench'); exps = json.load(open(sys.argv[1])); name = sys.argv[2]; reps = int(sys.argv[3])
e = exps[name]; res = B / 'results'; res.mkdir(exist_ok=True); tsv = res / f'{name}.tsv'
if not tsv.exists(): tsv.write_text('exp\tarm\trep\tthreads\twall_s\tmaxrss_kb\tidxml_sha256\tpin_sha256\tstatus\n')
def norm_idxml(p):  # drop the run date, engine version and paths, which differ between builds and output folders
    h = hashlib.sha256()
    for line in open(p, 'rb'):
        line = re.sub(rb' date="[^"]*"| search_engine_version="[^"]*"| db="[^"]*"| location="[^"]*"', b'', line)
        if b'name="ProSE:' in line or b'spectra_data' in line: continue
        h.update(line)
    return h.hexdigest()
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
done = {tuple(l.split('\t')[1:3]) for l in tsv.read_text().splitlines()[1:]}
for rep in range(1, reps + 1):
    for arm, a in e['arms'].items():
        if (arm, str(rep)) in done: continue
        out = B / 'runs' / name / f'{arm}_r{rep}'
        subprocess.run(['rm', '-rf', str(out)])
        r = subprocess.run([str(B / 'scripts/run1.sh'), a['prefix'], a['ini'], str(out), str(e['threads']), e['db'], *e['in'], '--', *a.get('args', [])],
                           cwd=B, capture_output=True, text=True)
        t = (out / 'time.txt').read_text() if (out / 'time.txt').exists() else ''
        m = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)', t)
        wall = sum(float(x) * 60 ** i for i, x in enumerate(reversed(m.group(1).split(':')))) if m else float('nan')
        rss = re.search(r'Maximum resident set size \(kbytes\): (\d+)', t); rss = rss.group(1) if rss else ''
        ok = r.returncode == 0
        ih = ','.join(norm_idxml(p) for p in sorted(out.glob('*.idXML'))) if ok else ''
        ph = ','.join(sha(p) for p in sorted(out.glob('*.pin'))) if ok else ''
        status = 'ok' if ok else ('killed' if 'signal 9' in t else 'failed')
        with open(tsv, 'a') as f: f.write(f'{name}\t{arm}\t{rep}\t{e["threads"]}\t{wall:.2f}\t{rss}\t{ih}\t{ph}\t{status}\n')
        print(name, arm, rep, f'{wall:.1f}s', rss, status, flush=True)
        if rep > 1 and ok:  # keep the outputs of rep 1 only
            for p in list(out.glob('*.idXML')) + list(out.glob('*.pin')): p.unlink()
