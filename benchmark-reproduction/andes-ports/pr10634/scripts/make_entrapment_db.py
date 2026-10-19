"""Target/decoy FASTA -> targets + 1:1 shuffled entrapment twins (andes build_shuffled_entrap.py, seed 1)
+ DECOY_ reversals of both. usage: make_entrapment_db.py <td.fasta> <out.fasta>"""
import random, sys
def read(p):
    h=None; s=[]
    for l in open(p):
        l=l.rstrip('\n')
        if l.startswith('>'):
            if h is not None: yield h, ''.join(s)
            h=l[1:]; s=[]
        else: s.append(l.strip())
    if h is not None: yield h, ''.join(s)
rng=random.Random(1)
t=[(h,s) for h,s in read(sys.argv[1]) if not h.startswith('DECOY_')]
e=[]
for h,s in t:
    c=list(s); rng.shuffle(c); e.append(('ENTRAP_'+h.split()[0], ''.join(c)))
with open(sys.argv[2],'w') as f:
    for h,s in t+e: f.write(f'>{h}\n{s}\n')
    for h,s in t+e: f.write(f'>DECOY_{h.split()[0]}\n{s[::-1]}\n')
print(len(t), 'targets', len(e), 'entrapment')
