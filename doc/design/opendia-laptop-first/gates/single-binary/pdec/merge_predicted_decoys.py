"""Assemble a PQP with *predicted* decoys: targets.pqp (PGENT2 targets + entrapments, decoys none) + decoys.pqp (DIALibGen
prediction of the pseudo-reversed sequences in decoy_seqs.fasta, same settings, no cleavage). Every decoy precursor gets
DECOY=1 and TRAML_ID 'DECOY_<target TRAML_ID>' of the target precursor with the same charge and the reversed sequence
(exact reversed modification positions), so OpenDIA pairs them. Decoy precursors without a target partner are dropped; so are targets
without a decoy (and their entrapment status is kept). Usage: merge_predicted_decoys.py targets.pqp decoys.pqp out.pqp"""
import sqlite3, shutil, sys, re
T, Dp, OUT = sys.argv[1:4]
shutil.copy(T, OUT)
o=sqlite3.connect(OUT); d=sqlite3.connect(Dp)
# Pair by the exact pseudo-reversed *modified* sequence (modifications move with their residue). Pairing by oxidation
# count alone (v1) gave one decoy to only one of several same-count positional isoforms, leaving 3.4 % unpaired.
tok=lambda s: re.findall(r'([A-Z])((?:\([^)]*\))?)', s)
def rev(ms):
    t=tok(ms); t=t[:-1][::-1]+t[-1:]; return ''.join(a+m for a,m in t)
tgt={}
for pid,tid,seq,ms,z in o.execute("select p.ID,p.TRAML_ID,pe.UNMODIFIED_SEQUENCE,pe.MODIFIED_SEQUENCE,p.CHARGE from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID where p.DECOY=0"):
    tgt[(rev(ms),z)]=tid
off={t:o.execute(f"select coalesce(max(ID),0)+1 from {t}").fetchone()[0] for t in ('PRECURSOR','PEPTIDE','TRANSITION','PROTEIN','COMPOUND') if o.execute(f"select count(*) from sqlite_master where name='{t}'").fetchone()[0]}
keep=set(); kept_t=set()
for pid,seq,ms,z in d.execute("select p.ID,pe.UNMODIFIED_SEQUENCE,pe.MODIFIED_SEQUENCE,p.CHARGE from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID"):
    k=(ms,z)
    if k in tgt and tgt[k] not in kept_t: keep.add((pid,tgt[k])); kept_t.add(tgt[k])
P=off['PRECURSOR']; PE=off['PEPTIDE']; TR=off['TRANSITION']; PO=off['PROTEIN']
idmap={}
for pid,tid in keep: idmap[pid]=(pid+P,'DECOY_'+tid)
cols=lambda t:[r[1] for r in d.execute(f'pragma table_info({t})')]
def copy(t,rows,cl):
    o.executemany(f"insert into {t}({','.join(cl)}) values({','.join('?'*len(cl))})",rows)
pc=cols('PRECURSOR')
rows=[]
for r in d.execute("select * from PRECURSOR"):
    r=dict(zip(pc,r))
    if r['ID'] not in idmap: continue
    r['ID'],r['TRAML_ID']=idmap[r['ID']]; r['DECOY']=1; rows.append([r[c] for c in pc])
copy('PRECURSOR',rows,pc)
pp=[(a+P,b+PE) for a,b in d.execute("select PRECURSOR_ID,PEPTIDE_ID from PRECURSOR_PEPTIDE_MAPPING") if a in idmap]
copy('PRECURSOR_PEPTIDE_MAPPING',pp,['PRECURSOR_ID','PEPTIDE_ID'])
peps={b for _,b in pp}
pec=cols('PEPTIDE'); rows=[]
for r in d.execute("select * from PEPTIDE"):
    r=dict(zip(pec,r))
    if r['ID']+PE not in peps: continue
    r['ID']+=PE; r['DECOY']=1; rows.append([r[c] for c in pec])
copy('PEPTIDE',rows,pec)
prc=cols('PROTEIN'); rows=[]
for r in d.execute("select * from PROTEIN"):
    r=dict(zip(prc,r)); r['ID']+=PO; r['DECOY']=1; rows.append([r[c] for c in prc])
copy('PROTEIN',rows,prc)
copy('PEPTIDE_PROTEIN_MAPPING',[(a+PE,b+PO) for a,b in d.execute("select PEPTIDE_ID,PROTEIN_ID from PEPTIDE_PROTEIN_MAPPING") if a+PE in peps],['PEPTIDE_ID','PROTEIN_ID'])
tpm=[(a+TR,b+P) for a,b in d.execute("select TRANSITION_ID,PRECURSOR_ID from TRANSITION_PRECURSOR_MAPPING") if b in idmap]
trs={a for a,_ in tpm}
tc=cols('TRANSITION'); rows=[]
for r in d.execute("select * from TRANSITION"):
    r=dict(zip(tc,r))
    if r['ID']+TR not in trs: continue
    r['ID']+=TR; r['TRAML_ID']='DECOY_'+str(r['TRAML_ID']); r['DECOY']=1; rows.append([r[c] for c in tc])
copy('TRANSITION',rows,tc)
copy('TRANSITION_PRECURSOR_MAPPING',tpm,['TRANSITION_ID','PRECURSOR_ID'])
if o.execute("select count(*) from sqlite_master where name='TRANSITION_PEPTIDE_MAPPING'").fetchone()[0]:
    copy('TRANSITION_PEPTIDE_MAPPING',[(a+TR,b+PE) for a,b in d.execute("select TRANSITION_ID,PEPTIDE_ID from TRANSITION_PEPTIDE_MAPPING") if a+TR in trs],['TRANSITION_ID','PEPTIDE_ID'])
o.commit()
# The PQP reader INNER JOINs PEPTIDE_GENE_MAPPING when GENE has rows: map every decoy peptide to the decoy peptides' gene
# (DIALibGen's own gene table rows, re-numbered), else to the first target gene.
if o.execute("select count(*) from sqlite_master where name='GENE'").fetchone()[0] and o.execute("select count(*) from GENE").fetchone()[0]:
    gid=o.execute("select min(ID) from GENE").fetchone()[0]
    have={r[0] for r in o.execute("select PEPTIDE_ID from PEPTIDE_GENE_MAPPING")}
    o.executemany("insert into PEPTIDE_GENE_MAPPING(PEPTIDE_ID,GENE_ID) values(?,?)",[(p,gid) for p in peps if p not in have])
    o.commit()
nt=o.execute("select count(*) from PRECURSOR where DECOY=0").fetchone()[0]
print(f'targets {nt}, predicted decoys paired {len(keep)}, targets with a decoy {len(kept_t)}')
