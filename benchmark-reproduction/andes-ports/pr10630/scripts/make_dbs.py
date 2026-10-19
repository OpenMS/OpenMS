"""pyOpenMS-free equivalent of the frozen prepare_databases.py (FASTA part); output hashes are checked."""
import hashlib, json, sys
from pathlib import Path
B = Path(sys.argv[1])  # frozen/broad_benchmark
entries, hdr, seq = [], None, []
for line in (B / 'data/hye_targets.fasta').read_text().splitlines():
    if line.startswith('>'):
        if hdr is not None: entries.append((hdr, ''.join(seq)))
        hdr, seq = line[1:], []
    else: seq.append(line.strip())
if hdr is not None: entries.append((hdr, ''.join(seq)))
def split(h):
    h = h.strip(); i = h.find(' ')
    return (h, '') if i < 0 else (h[:i], h[i+1:].strip())
entries = [(*split(h), s) for h, s in entries]
expected = json.load(open(B.parent / 'databases.json'))
for db in ['hy', 'hye', 'human']:
    e = entries
    if db == 'hy': e = [x for x in e if x[0].endswith(('_HUMAN', '_YEAST')) or 'Cont_' in x[0]]
    if db == 'human': e = [x for x in e if x[0].endswith('_HUMAN') or 'Cont_' in x[0]]
    out = B / 'data' / (db + '_td.fasta')
    with out.open('w') as f:
        for i, d, s in e: f.write('>' + i + ' ' + d + '\n' + s + '\n')
        for i, d, s in e: f.write('>DECOY_' + i + '\n' + s[::-1] + '\n')
    h = hashlib.sha256(out.read_bytes()).hexdigest()
    print(db, len(e), h, 'OK' if h == expected[db]['sha256'] else 'MISMATCH')
