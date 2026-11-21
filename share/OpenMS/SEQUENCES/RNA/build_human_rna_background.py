#!/usr/bin/env python3
"""Build the human background RNA database of NucleicAcidSearchEngine ('background:enable').

usage: build_human_rna_background.py <download dir> <output FASTA>

The download dir must contain (the script prints the commands to fetch them if they are missing):
  hg38-tRNAs.tar.gz  GtRNAdb, http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz
  NC_012920.gb       NCBI nuccore NC_012920.1 (mitochondrial genome, GenBank with parts)
  rrna.fasta         NCBI nuccore NR_046235.3, NR_023363.1 (FASTA)
  smallrna.fasta     NCBI nuccore SMALL_RNAS below (FASTA)

Output entries (accessions start with 'BG_'):
  - every mature tRNA of hg38-mature-tRNAs.fa, with 5' phosphate ('p') and each 3' end CCA, CC, C and none
  - every mitochondrial tRNA annotated in NC_012920.1, with 5' phosphate and 3' CCA
  - the mitochondrial genome on both strands (both are transcribed)
  - 45S pre-rRNA and 5S rRNA
  - the small RNAs in SMALL_RNAS
"""
import re
import sys
import tarfile

RRNAS = ['NR_046235.3', 'NR_023363.1']
SMALL_RNAS = ['NR_002715.1', 'NR_001445.2', 'NR_004430.2', 'NR_002716.3', 'NR_003925.1', 'NR_002756.2', 'NR_004394.1',
              'NR_004391.1', 'NR_004392.1', 'NR_004393.1', 'NR_001571.2', 'NR_026703.1', 'NR_002312.1', 'NR_003051.3',
              'NR_002999.2']
EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore'


def fasta(lines):
    entries, name = [], None
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            entries.append([line[1:], '']); name = line
        elif name is not None and line:
            entries[-1][1] += line.upper().replace('T', 'U')
    return entries


def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans('ACGU', 'UGCA'))


def main(download, out):
    try:
        with tarfile.open(download + '/hg38-tRNAs.tar.gz') as tar:
            member = next(m for m in tar.getmembers() if m.name.endswith('hg38-mature-tRNAs.fa'))
            trnas = fasta(tar.extractfile(member).read().decode().splitlines())
        genbank = open(download + '/NC_012920.gb').read()
        rrnas = fasta(open(download + '/rrna.fasta'))
        small = fasta(open(download + '/smallrna.fasta'))
    except (OSError, StopIteration) as error:
        sys.exit(f'{error}\nFetch the inputs with:\n'
                 f'  curl -L -o hg38-tRNAs.tar.gz http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz\n'
                 f'  curl -o NC_012920.gb "{EFETCH}&id=NC_012920.1&rettype=gbwithparts&retmode=text"\n'
                 f'  curl -o rrna.fasta "{EFETCH}&id={",".join(RRNAS)}&rettype=fasta&retmode=text"\n'
                 f'  curl -o smallrna.fasta "{EFETCH}&id={",".join(SMALL_RNAS)}&rettype=fasta&retmode=text"')

    entries = []
    for header, seq in trnas:
        name = header.split()[0]
        full = seq if seq.endswith('CCA') else seq + 'CCA'
        for cut, end in ((0, 'CCA'), (1, 'CC'), (2, 'C'), (3, 'none')):
            entries.append((f'BG_{name}_3end_{end} GtRNAdb hg38 mature tRNA, 3\' end {end}', 'p' + (full[:-cut] if cut else full)))

    origin = genbank.split('\nORIGIN')[1].split('//')[0]
    # rCRS keeps an 'N' at position 3107 as a placeholder for a deletion; dropping it gives the actual sequence (the
    # annotated coordinates refer to the sequence with the placeholder, so tRNAs are cut out before it is removed)
    chrm_n = re.sub(r'[\d\s]', '', origin).upper().replace('T', 'U')
    chrm = chrm_n.replace('N', '')
    for match in re.finditer(r'\n     tRNA\s+(complement\()?<?(\d+)\.\.>?(\d+)\)?\n(.*?)(?=\n     \S)', genbank, re.S):
        complement, start, end, body = match.groups()
        seq = chrm_n[int(start) - 1:int(end)]
        if complement: seq = reverse_complement(seq)
        product = re.search(r'/product="([^"]+)"', body).group(1).replace(' ', '_')
        entries.append((f'BG_MT_{product}_{start} NC_012920.1 mitochondrial tRNA {start}..{end}{" (minus)" if complement else ""}, 3\' CCA',
                        'p' + seq + 'CCA'))
    entries.append(('BG_NC_012920.1_plus NC_012920.1 mitochondrial genome, plus strand', chrm))
    entries.append(('BG_NC_012920.1_minus NC_012920.1 mitochondrial genome, minus strand', reverse_complement(chrm)))

    for header, seq in rrnas + small:
        accession, _, description = header.partition(' ')
        entries.append((f'BG_{accession} {description}', seq))

    with open(out, 'w') as handle:
        for header, seq in entries:
            if not re.fullmatch(r'p?[ACGU]+', seq):
                sys.exit(f'unexpected characters in {header.split()[0]}')
            handle.write(f'>{header}\n')
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + '\n')
    print(f'{out}: {len(entries)} sequences ({len(trnas)} tRNAs x 4 3\' ends)')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
