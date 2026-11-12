"""Pairwise rank inversions against the SAME individual decoy competitors."""
import numpy as np

def reversals(qfull,qlin,dfull,dlin):
    qfull,qlin,dfull,dlin=map(np.asarray,(qfull,qlin,dfull,dlin))
    grid=np.unique(dlin)
    did=np.searchsorted(grid,dlin)+1
    order=np.argsort(dfull,kind='stable')
    qo=np.argsort(qfull,kind='stable')
    def scan(inclusive):
        bit=np.zeros(len(grid)+1,dtype=np.int64)
        out=np.zeros(len(qfull),dtype=np.int64)
        pos=0
        def prefix(i):
            total=0
            while i>0:
                total+=bit[i]
                i-=i&-i
            return total
        for qi in qo:
            while pos<len(order) and (dfull[order[pos]]<=qfull[qi] if inclusive else dfull[order[pos]]<qfull[qi]):
                j=did[order[pos]]
                while j<len(bit):
                    bit[j]+=1
                    j+=j&-j
                pos+=1
            if inclusive:
                # linear wins that full no longer strictly loses to
                out[qi]=prefix(np.searchsorted(grid,qlin[qi],'left'))
            else:
                # full strictly wins but linear strictly loses
                out[qi]=pos-prefix(np.searchsorted(grid,qlin[qi],'right'))
        return out
    up=scan(False)
    down=np.searchsorted(np.sort(dlin),qlin,'left')-scan(True)
    return down/len(dlin),up/len(dlin)

if __name__=='__main__':
    from pathlib import Path
    import pandas as pd
    # Brute-force helper test includes ties and a query also present in the null set.
    qf=np.array([1,2,0,1]);ql=np.array([0,1,2,1]);df=np.array([0,1,1,2]);dl=np.array([2,1,0,1])
    d,u=reversals(qf,ql,df,dl)
    np.testing.assert_allclose(d,((qf[:,None]<df)&(ql[:,None]>dl)).mean(1))
    np.testing.assert_allclose(u,((qf[:,None]>df)&(ql[:,None]<dl)).mean(1))
    O=Path(__file__).resolve().parent
    results=[]
    for raw in ['raw01','raw02']:
        W=pd.read_csv(O/f'reference_winners_{raw}.tsv.gz',sep='\t')
        for f in range(3):
            ix=W.index[W.fold==f]
            dc=W[(W.fold==f)&W.group.eq('DEC')]
            d,u=reversals(W.loc[ix,'full'],W.loc[ix,'linear'],dc.full,dc.linear)
            W.loc[ix,'pairwise_down']=d
            W.loc[ix,'pairwise_up']=u
        masks={'D':W.stage.eq('D right window, below cut'),'E':W.stage.eq('E accepted'),'allH':W.group.eq('H'),'ENT':W.group.eq('ENT'),'DEC':W.group.eq('DEC')}
        for g,m in masks.items():
            G=W[m]
            results.append(dict(raw=raw,group=g,n=len(G),mean_down_fraction=G.pairwise_down.mean(),mean_up_fraction=G.pairwise_up.mean(),rows_any_down=int((G.pairwise_down>0).sum()),rows_any_up=int((G.pairwise_up>0).sum()),rows_more_down_than_up=int((G.pairwise_down>G.pairwise_up).sum())))
        W[['row','pairwise_down','pairwise_up']].to_csv(O/f'pairwise_reversals_{raw}.tsv.gz',sep='\t',index=False)
    pd.DataFrame(results).to_csv(O/'pairwise_reversal_summary.tsv',sep='\t',index=False)
    print('Brute-force reversal test passed, both raws complete')
