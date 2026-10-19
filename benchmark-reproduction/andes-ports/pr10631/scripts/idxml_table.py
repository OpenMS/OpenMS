#!/usr/bin/env python3
"""idXML (FDR-filtered) -> TSV: spectrum_reference, rank-1 sequence, protein accessions, target_decoy, score."""
import sys, csv
import xml.etree.ElementTree as ET

def main(idxml, out):
    prot = {}
    w = csv.writer(open(out, 'w', newline=''), delimiter='\t')
    w.writerow(['spectrum_reference', 'sequence', 'accessions', 'target_decoy', 'score'])
    for ev, el in ET.iterparse(idxml, events=('end',)):
        if el.tag == 'ProteinHit':
            prot[el.get('id')] = el.get('accession')
        elif el.tag == 'PeptideIdentification':
            hits = el.findall('PeptideHit')
            if hits:
                h = hits[0]
                hup = {u.get('name'): u.get('value') for u in h.findall('UserParam')}
                acc = ';'.join(prot.get(r, r) for r in (h.get('protein_refs') or '').split())
                w.writerow([el.get('spectrum_reference'), h.get('sequence'), acc, hup.get('target_decoy', ''), h.get('score')])
            el.clear()

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
