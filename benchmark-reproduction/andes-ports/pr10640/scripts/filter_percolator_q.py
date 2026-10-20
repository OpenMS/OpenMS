#!/usr/bin/env python3
"""Drop PeptideIdentifications whose best hit has Percolator q-value (MS:1001491) above a threshold.
Usage: filter_percolator_q.py in.idXML out.idXML max_q"""
import re, sys
max_q = float(sys.argv[3]); kept = dropped = 0
pat = re.compile(r'name="MS:1001491" value="([^"]+)"')
with open(sys.argv[1]) as f, open(sys.argv[2], 'w') as o:
    block = None
    for line in f:
        if block is None:
            if '<PeptideIdentification ' in line:
                block = [line]
            else:
                o.write(line)
            continue
        block.append(line)
        if '</PeptideIdentification>' in line:
            q = next((float(m.group(1)) for l in block for m in [pat.search(l)] if m), 0.0)
            if q <= max_q:
                o.writelines(block); kept += 1
            else:
                dropped += 1
            block = None
print(f'kept {kept}, dropped {dropped}', file=sys.stderr)
