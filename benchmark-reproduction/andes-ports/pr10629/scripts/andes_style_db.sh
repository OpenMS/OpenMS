# andes 1:1 shuffled entrapment (docs/benchmarks/glyco/build_shuffled_entrap.py, verbatim, bigbio/andes cab36f6) on a target FASTA,
# accessions moved to prefix form (ENTRAP_<header>) for EntrapmentEvaluator, then #10364 full-reversal DECOY_ decoys.
# usage: andes_style_db.sh <andes script> <targets.fasta> <seed> <out.fasta> <swap_decoys.py>
set -e
tmp=$(mktemp -d)
python3 -I "$1" "$2" $tmp/e.fasta "$3" > /dev/null
python3 -I -c "
import sys
with open(sys.argv[2], 'w') as o:
    for l in open(sys.argv[1]):
        if l.startswith('>') and 'ENTRAP_' in l:
            p = l[1:].rstrip('\n').split('|')
            if len(p) >= 3: p[1] = p[1].replace('ENTRAP_', '', 1); p[2] = p[2].replace('ENTRAP_', '', 1)
            l = '>ENTRAP_' + '|'.join(p).replace('ENTRAP_', '', 1 if not l[1:].startswith('ENTRAP_') else 0) + '\n'
        o.write(l)
" $tmp/e.fasta $tmp/p.fasta
python3 -I "$5" $tmp/p.fasta "$4"
rm -r $tmp
