"""Accepted entrapments (q <= 0.01) of D seed 0: is their parent also a target winner, and at the same apex?
(Diagnostic for the VERDICT reading; needs rt_apex, so reads the full per-run table on codon output if present.)"""
import sqlite3, glob, pandas as pd, numpy as np, sys
L='hela_entrap_dec.pqp'; c=sqlite3.connect(f'file:{L}?mode=ro',uri=True)
m=pd.read_sql("""select p.ID pid, pe.UNMODIFIED_SEQUENCE seq, p.CHARGE z, pr.PROTEIN_ACCESSION acc from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING pm on pm.PRECURSOR_ID=p.ID
 join PEPTIDE pe on pe.ID=pm.PEPTIDE_ID join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID where p.DECOY=0""",c)
m=m.drop_duplicates('pid')
# parent of an entrapment: the target peptide whose shuffle it is (from the FASTA pairs)
pairs={}; prev=None
for l in open('hela_entrap.fasta'):
    if l.startswith('>'): name=l[1:].split('|')[1]
    else:
        s=l.strip()
        if name.endswith('_p_target'): pairs[s]=prev
        else: prev=s
m['ent']=m.acc.str.contains('_p_target')
m['parent']=np.where(m.ent,m.seq.map(pairs),None)
key={(s,z):p for p,s,z,e in zip(m.pid,m.seq,m.z,m.ent) if not e}
for arm in sys.argv[1:]:
    f=glob.glob(f'full/D_{arm}_s0/*.precursors.tsv')[0]
    t=pd.read_csv(f,sep='\t',usecols=['precursor_id','is_decoy','score','q_value','rt_apex'])
    w=t[t.score>-1e29].sort_values('score',ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
    M=m.set_index('pid')
    acc=w[(w.is_decoy==0)&(w.q_value<=0.01)]
    e=acc[acc.index.map(M.ent).fillna(False).astype(bool)]
    par=[key.get((M.parent[i],M.z[i])) for i in e.index]
    has=[p is not None and p in w.index for p in par]
    acc_par=[p is not None and p in acc.index for p in par]
    dt=[abs(w.rt_apex[p]-e.rt_apex[i]) if (p is not None and p in w.index) else np.nan for p,i in zip(par,e.index)]
    dt=np.array(dt)
    print(f'{arm}: accepted entrapments {len(e)}; parent (same charge) also accepted {np.mean(acc_par):.2f}; '
          f'parent apex within 5 s {np.nanmean(dt<=5):.2f} (of those with a parent winner {np.mean(has):.2f})')
