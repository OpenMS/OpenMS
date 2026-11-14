"""Rejected-vs-accepted signal diagnostics (evaluation only). Expected false among the rejected (sum of DIA-NN PEP);
fragments per precursor in DIA-NN's library; sparsity of OpenDIA's binned fragment arrays at the winner."""
import numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = B + '/gap-20261114/armE_raw01'
exec(open(B + '/gap-20261114/funnel.py').read().split("# Why stage 1 loses")[0].split('print(\'\\nFunnel')[0])
ra = pd.read_parquet(B + '/diann-pairlib-20261113/report_raw01.parquet', columns=['Modified.Sequence', 'Precursor.Charge', 'PEP', 'Q.Value', 'FWHM', 'RT.Start', 'RT.Stop'])
ra['key'] = ra['Modified.Sequence'] + '/' + ra['Precursor.Charge'].astype(str); ra = ra.drop_duplicates('key').set_index('key')
hp = human.set_index('precursor_id')
wr = w.set_index('precursor_id')
rej = sorted(near_is_winner - accepted); acc = sorted(near_is_winner & accepted)
for nm, ids in (('accepted', acc), ('rejected', rej)):
    r = ra.loc[hp.loc[ids, 'key']]
    print(f'{nm}: DIA-NN sum PEP {r.PEP.sum():.1f} (expected false), PEP>0.05 {int((r.PEP > 0.05).sum())}, FWHM p50 {r.FWHM.median()*60:.2f} s, width p50 {(r["RT.Stop"]-r["RT.Start"]).median()*60:.1f} s')
lt = pd.read_csv(B + '/diann-pairlib-20261113/lib_targets.tsv', sep='\t', nrows=200000)
print('DIA-NN library columns', list(lt.columns)[:12])
pc = [c for c in lt.columns if c.lower() in ('precursorid', 'modifiedpeptide', 'modified.sequence', 'modifiedsequence', 'transition_group_id')]
print('fragments per precursor in DIA-NN library (first 200k rows):', lt.groupby(pc[0]).size().describe()[['mean', 'min', 'max']].to_dict() if pc else 'n/a')
XIC = np.load(f'{D}/xic.npy', mmap_mode='r')
for nm, ids in (('accepted', acc), ('rejected', rej)):
    rows = np.sort(wr.loc[ids, 'row'].values)
    x = np.asarray(XIC[rows], np.float64)[:, :, 10:22]
    nz = (x > 0)
    print(f'{nm}: central 12 bins nonzero share per fragment p50 {np.median(nz.mean(2)):.2f}; fragments with >=3 nonzero bins p50 {np.median((nz.sum(2) >= 3).sum(1)):.0f}; '
          f'summed-profile nonzero bins p50 {np.median(nz.any(1).sum(1)):.0f}/12; apex/median-of-others ratio p50 '
          f'{np.median(x.sum(1).max(1) / (np.median(np.asarray(XIC[rows], np.float64).sum(1), 1) + 1)):.1f}')
