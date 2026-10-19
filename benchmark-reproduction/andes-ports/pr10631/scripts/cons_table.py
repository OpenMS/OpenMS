#!/usr/bin/env python3
"""Flatten an isobaric consensusXML (IsobaricAnalyzer+IDMapper or IsobaricWorkflow) to TSV:
one row per consensus feature with scan ids, purities, channel intensities and the best peptide hit."""
import sys, csv
import xml.etree.ElementTree as ET

def main(cxml, out):
    prot = {}
    nmaps = 0
    w = None
    for ev, el in ET.iterparse(cxml, events=('end',)):
        tag = el.tag
        if tag == 'ProteinHit':
            prot[el.get('id')] = el.get('accession')
        elif tag == 'map':
            nmaps = max(nmaps, int(el.get('id')) + 1)
        elif tag == 'consensusElement':
            if w is None:
                w = csv.writer(open(out, 'w', newline=''), delimiter='\t')
                w.writerow(['scan_id', 'id_scan_id', 'precursor_purity', 'quant_precursor_purity', 'sequence', 'accessions', 'target_decoy', 'score'] + [f'ch{i}' for i in range(nmaps)])
            ints = [0.0] * nmaps
            for e in el.find('groupedElementList'):
                ints[int(e.get('map'))] = float(e.get('it'))
            up = {u.get('name'): u.get('value') for u in el.findall('UserParam')}
            seq = acc = td = score = ''
            pid = el.find('PeptideIdentification')
            if pid is not None:
                hits = pid.findall('PeptideHit')
                if hits:
                    h = hits[0]
                    seq = h.get('sequence')
                    acc = ';'.join(prot.get(r, r) for r in (h.get('protein_refs') or '').split())
                    hup = {u.get('name'): u.get('value') for u in h.findall('UserParam')}
                    td = hup.get('target_decoy', '')
                    score = h.get('score')
            w.writerow([up.get('scan_id', ''), up.get('id_scan_id', ''), up.get('precursor_purity', ''), up.get('quant_precursor_purity', ''), seq, acc, td, score] + ints)
            el.clear()
    if w is None:
        csv.writer(open(out, 'w', newline=''), delimiter='\t').writerow(['empty'])

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
