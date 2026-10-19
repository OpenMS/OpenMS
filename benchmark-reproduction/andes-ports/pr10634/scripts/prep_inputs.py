"""Recreate the frozen #10364 inputs (8,000 MS2 per file, three target/decoy FASTAs) with the bundle's own
prepare.py / prepare_databases.py; checks every hash against inputs.json / databases.json."""
import hashlib, json, subprocess, sys
from pathlib import Path
W = Path(sys.argv[1]); B = W/'broad_benchmark'; PY = sys.argv[2]
INPUTS = json.loads((W/'inputs.json').read_text())
def digest(p):
    with open(p, 'rb') as f: return hashlib.file_digest(f, 'sha256').hexdigest()
def fetch(url, path, sha, size):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and digest(path) == sha: return
    tmp = Path(str(path)+'.part')
    for _ in range(200):  # PRIDE FTP resets long transfers: resume until complete
        if tmp.exists() and tmp.stat().st_size >= size: break
        subprocess.run(['curl','-sSL','-C','-','--connect-timeout','30','--speed-limit','50000','--speed-time','60',url,'-o',tmp])
    assert tmp.stat().st_size == size, f'size {path}'
    assert digest(tmp) == sha, f'sha {path}'
    tmp.rename(path)
(B/'data').mkdir(exist_ok=True); (B/'raw').mkdir(exist_ok=True)
f = json.load(open(B/'provenance/fasta_integrity_correction.json'))
fetch(f['source'], B/'data/hye_targets.fasta', f['complete_sha256'], f['complete_bytes'])
if not (B/'data/hy_td.fasta').exists(): subprocess.run([PY, B/'prepare_databases.py'], check=True)
for name, info in json.load(open(W/'databases.json')).items():
    print(name, digest(B/'data'/(name+'_td.fasta')) == info['sha256'], flush=True)
ids = sys.argv[3:] or list(INPUTS)
for d in ids:
    out = B/'data'/(d+'.mzML')
    if not out.exists():
        s = INPUTS[d]['source']; raw = B/'raw'/(d+('.mzML.gz' if s['url'].endswith('.gz') else '.mzML'))
        fetch(s['url'], raw, s['sha256'], s['bytes'])
        subprocess.run([PY, B/'prepare.py', d], check=True)
        raw.unlink()
    print(d, 'OK' if digest(out) == INPUTS[d]['selected_sha256'] else 'MISMATCH', flush=True)
