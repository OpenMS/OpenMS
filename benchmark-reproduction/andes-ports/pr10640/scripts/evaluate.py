#!/usr/bin/env python3
"""Count accepted PSMs in a filtered idXML (output of FalseDiscoveryRate) and label variant PSMs with the
PXD022787 donor genotypes.

A PSM (best hit per spectrum, target) is
  variant   - every protein is a variant entry (<acc>_VAR_<ref><pos><alt>),
  knockoff  - every protein is a knockoff entry (<acc>_KO<i>_..._VAR_...),
  other     - otherwise (canonical; also PSMs shared between variant/knockoff and canonical entries).
A variant/knockoff PSM is 'carried' if its peptide, extended to the substitution site if the site is next to it
(I = L), occurs in a protein sequence of the same gene that a
donor haplotype carries (personal FASTA tag with xHy or A), 'absent' if not but the peptide with the reference
residue does, else 'unlabelled' (also when the gene has no personal sequence). True FDP = absent / (absent + carried).
Usage: evaluate.py <variant db fasta> <personal fasta> <Swiss-Prot fasta> <label> <out.json> <in.idXML>..."""
import json, re, sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

VAR = re.compile(r'^(.+)_VAR_([A-Z])(\d+)([A-Z])$')
KO = re.compile(r'^(.+)_KO\d+_([A-Z])(\d+)([A-Z])_VAR_[A-Z]\d+[A-Z]$')

def read_fasta_full(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None: yield name, ''.join(seq)
                name, seq = line[1:].rstrip('\n'), []
            else:
                seq.append(line.strip())
    if name is not None: yield name, ''.join(seq)

def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None: yield name, ''.join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None: yield name, ''.join(seq)

class Truth:
    """donor-carried protein sequences per gene (personal FASTA gname:: = UniProt GN=)"""
    def __init__(self, personal, sprot):
        carried = defaultdict(list)
        for header, seq in read_fasta_full(personal):
            tag = header.split()[0]
            gene = re.search(r'gname::(\S+)', header)
            parts = tag.split('_')[1:]
            if gene and any(p == 'A' or re.fullmatch(r'\dH\d', p) for p in parts):
                carried[gene.group(1)].append(seq.replace('I', 'L'))
        self.text = {g: '|'.join(v) for g, v in carried.items()}
        self.gene = {}
        for header, seq in read_fasta_full(sprot):
            gn = re.search(r' GN=(\S+)', header)
            if gn: self.gene[header.split()[0].split('|')[1]] = gn.group(1)
    def gene_of(self, acc):
        """gene of a variant/knockoff accession 'sp|P04217|A1BG_HUMAN_VAR_H52R'"""
        parts = acc.split('|')
        return self.gene.get(parts[1]) if len(parts) > 2 else None
    def has(self, pep, gene):
        return gene in self.text and pep.replace('I', 'L') in self.text[gene]

def strip(seq):
    return re.sub(r'\([^)]*\)', '', seq).replace('.', '')

def label(pep, entries, seqs, truth):
    """carried/absent/unlabelled for a peptide on variant or knockoff entries (site from the accession). The
    peptide is extended to the substitution site when the site lies next to it (a new or lost cleavage site),
    so that the label is about the substitution, not about the peptide alone."""
    genes = {truth.gene_of(a) for a in entries} - {None}
    if not any(g in truth.text for g in genes): return 'unlabelled'
    labels = []
    for acc in entries:
        m = KO.match(acc) or VAR.match(acc)
        ref, pos = m.group(2), int(m.group(3)) - 1
        s, gene = seqs.get(acc), truth.gene_of(acc)
        if s is None or gene not in truth.text: continue
        start = s.find(pep)
        while start != -1:
            a, b = min(start, pos), max(start + len(pep), pos + 1)
            if b - a <= len(pep) + 1:  # the site is in the peptide or next to it
                alt_win = s[a:b]
                ref_win = s[a:pos] + ref + s[pos + 1:b]
                labels.append('carried' if truth.has(alt_win, gene) else
                              'absent' if truth.has(ref_win, gene) else 'unlabelled')
            start = s.find(pep, start + 1)
    for t in ('carried', 'absent'):
        if t in labels: return t
    return 'unlabelled'

def main():
    db, personal, sprot, arm, out = sys.argv[1:6]
    seqs = {a: s for a, s in read_fasta(db) if '_VAR_' in a}
    truth = Truth(personal, sprot)
    res = {'arm': arm, 'files': sys.argv[6:]}
    counts = Counter(); variant_ids = defaultdict(set); strata = defaultdict(Counter); per_file = defaultdict(Counter)
    for path in sys.argv[6:]:
        prot = {}
        for ev, el in ET.iterparse(path, events=('end',)):
            tag = el.tag
            if tag == 'ProteinHit':
                prot[el.get('id')] = el.get('accession')
            elif tag == 'PeptideIdentification':
                hits = el.findall('PeptideHit')
                if not hits: el.clear(); continue
                h = hits[0]
                accs = [prot.get(r, r) for r in (h.get('protein_refs') or '').split()]
                meta = {u.get('name'): u.get('value') for u in h.findall('UserParam')}
                if meta.get('target_decoy', 'target') == 'decoy' or all(a.startswith('DECOY_') for a in accs):
                    counts['decoy'] += 1; el.clear(); continue
                accs = [a for a in accs if not a.startswith('DECOY_')]
                counts['target'] += 1
                origin = next((u.get('value') for u in el.findall('UserParam') if u.get('name') == 'file_origin'), path)
                run = origin.rsplit('/', 1)[-1].replace('.idXML', '')
                per_file[run]['target'] += 1
                pep = strip(h.get('sequence'))
                if accs and all(KO.match(a) for a in accs):
                    cls = 'knockoff'
                elif accs and all(VAR.match(a) and not KO.match(a) for a in accs):
                    cls = 'variant'
                else:
                    cls = 'other'
                counts[cls] += 1
                per_file[run][cls] += 1
                if cls != 'other':
                    t = label(pep, accs, seqs, truth)
                    counts[f'{cls}_{t}'] += 1
                    per_file[run][f'{cls}_{t}'] += 1
                    if cls == 'variant':
                        variant_ids[t].add(tuple(sorted(accs)))
                        st = meta.get('knockoff_stratum', '')
                        strata[st][t] += 1
                el.clear()
    res['counts'] = dict(counts)
    for cls in ('variant', 'knockoff'):
        a, c = counts[f'{cls}_absent'], counts[f'{cls}_carried']
        res[f'{cls}_fdp'] = a / (a + c) if a + c else None
    va, vc = len(variant_ids['absent']), len(variant_ids['carried'])
    res['variants'] = {'carried': vc, 'absent': va, 'unlabelled': len(variant_ids['unlabelled']),
                       'fdp': va / (va + vc) if va + vc else None}
    res['strata'] = {k: dict(v) for k, v in strata.items()}
    res['per_file'] = {k: dict(v) for k, v in sorted(per_file.items())}
    json.dump(res, open(out, 'w'), indent=1)
    print(json.dumps({k: res[k] for k in ('arm', 'counts', 'variant_fdp', 'variants')}))

if __name__ == '__main__':
    main()
