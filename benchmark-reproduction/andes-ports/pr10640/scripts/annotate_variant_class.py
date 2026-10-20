#!/usr/bin/env python3
"""Add the PSM meta value 'variant_class' (variant / canonical) for FalseDiscoveryRate -PSM_class (PR #10636).
A hit is 'variant' if every protein it maps to is a variant entry (<acc>_VAR_<ref><pos><alt>, decoys included).
Usage: annotate_variant_class.py in.idXML out.idXML"""
import re, sys
VAR = re.compile(r'_VAR_[A-Z]\d+[A-Z]$')
acc_of = {}
n = {'variant': 0, 'canonical': 0}
with open(sys.argv[1]) as f, open(sys.argv[2], 'w') as o:
    for line in f:
        if '<ProteinHit ' in line:
            m_id, m_acc = re.search(r' id="([^"]+)"', line), re.search(r' accession="([^"]+)"', line)
            if m_id and m_acc: acc_of[m_id.group(1)] = m_acc.group(1)
        elif '<PeptideHit ' in line:
            refs = re.search(r' protein_refs="([^"]*)"', line)
            accs = [acc_of.get(r, r) for r in (refs.group(1).split() if refs else [])]
            cls = 'variant' if accs and all(VAR.search(a) for a in accs) else 'canonical'
            n[cls] += 1
            indent = re.match(r'\s*', line).group(0)
            param = f'{indent}\t<UserParam type="string" name="variant_class" value="{cls}"/>\n'
            if line.rstrip().endswith('/>'):
                o.write(line.rstrip()[:-2].rstrip() + ' >\n' + param + indent + '</PeptideHit>\n')
            else:
                o.write(line + param)
            continue
        o.write(line)
print(n, file=sys.stderr)
