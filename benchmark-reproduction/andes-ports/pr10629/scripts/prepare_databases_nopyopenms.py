"""FASTA part of the frozen broad_benchmark/prepare_databases.py without pyopenms (same output, checked by sha256)."""
import hashlib, sys
from pathlib import Path
data = Path(sys.argv[1])
entries, hdr, seq = [], None, []
for line in (data / 'hye_targets.fasta').read_text().splitlines():
    if line.startswith('>'):
        if hdr is not None: entries.append((hdr, ''.join(seq)))
        hdr, seq = line[1:], []
    else: seq.append(line.strip())
entries.append((hdr, ''.join(seq)))
def split(h):
    h = h.strip(); parts = h.split(None, 1)
    return parts[0], (parts[1] if len(parts) > 1 else '')
E = [(*split(h), s) for h, s in entries]
for db in ['hy', 'hye', 'human']:
    es = E
    if db == 'hy': es = [e for e in E if e[0].endswith(('_HUMAN', '_YEAST')) or 'Cont_' in e[0]]
    if db == 'human': es = [e for e in E if e[0].endswith('_HUMAN') or 'Cont_' in e[0]]
    out = data / (db + '_td.fasta')
    with out.open('w') as f:
        for i, d, s in es: f.write('>' + i + ' ' + d + '\n' + s + '\n')
        for i, d, s in es: f.write('>DECOY_' + i + '\n' + s[::-1] + '\n')
    print(db, len(es), hashlib.sha256(out.read_bytes()).hexdigest())
