#!/usr/bin/env python3
"""Per-spectrum table from an mzML: native id, ms level, scan window, FAIMS CV,
isolation target / selected ion m/z, spectrumRef. Streams the file (iterparse)."""
import sys, gzip, csv
import xml.etree.ElementTree as ET

def opener(p):
    return gzip.open(p, 'rb') if p.endswith('.gz') else open(p, 'rb')

def main(mzml, out):
    w = csv.writer(open(out, 'w', newline=''), delimiter='\t')
    w.writerow(['native_id', 'ms_level', 'rt', 'sw_lo', 'sw_hi', 'faims_cv', 'iso_target', 'iso_lo', 'iso_hi', 'selected_mz', 'charge', 'spectrum_ref'])
    for ev, el in ET.iterparse(opener(mzml), events=('end',)):
        tag = el.tag.rsplit('}', 1)[-1]
        if tag != 'spectrum':
            if tag == 'chromatogramList': break
            continue
        rec = dict(native_id=el.get('id'), ms_level='', rt='', sw_lo='', sw_hi='', faims_cv='', iso_target='', iso_lo='', iso_hi='', selected_mz='', charge='', spectrum_ref='')
        first_prec = None
        for sub in el.iter():
            t = sub.tag.rsplit('}', 1)[-1]
            if t == 'precursor' and first_prec is None:
                first_prec = sub
                rec['spectrum_ref'] = sub.get('spectrumRef', '')
            if t != 'cvParam': continue
            acc, val = sub.get('accession'), sub.get('value')
            if acc == 'MS:1000511': rec['ms_level'] = val
            elif acc == 'MS:1000016': rec['rt'] = val
            elif acc == 'MS:1000501' and rec['sw_lo'] == '': rec['sw_lo'] = val
            elif acc == 'MS:1000500' and rec['sw_hi'] == '': rec['sw_hi'] = val
            elif acc == 'MS:1001581': rec['faims_cv'] = val
        if first_prec is not None:
            for sub in first_prec.iter():
                if sub.tag.rsplit('}', 1)[-1] != 'cvParam': continue
                acc, val = sub.get('accession'), sub.get('value')
                if acc == 'MS:1000827': rec['iso_target'] = val
                elif acc == 'MS:1000828': rec['iso_lo'] = val
                elif acc == 'MS:1000829': rec['iso_hi'] = val
                elif acc == 'MS:1000744' and rec['selected_mz'] == '': rec['selected_mz'] = val
                elif acc == 'MS:1000041' and rec['charge'] == '': rec['charge'] = val
        w.writerow([rec[k] for k in ['native_id', 'ms_level', 'rt', 'sw_lo', 'sw_hi', 'faims_cv', 'iso_target', 'iso_lo', 'iso_hi', 'selected_mz', 'charge', 'spectrum_ref']])
        el.clear()

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
