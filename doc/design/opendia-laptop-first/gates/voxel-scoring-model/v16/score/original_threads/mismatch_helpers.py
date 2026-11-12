"""Label-free selector-comparison helper; stage labels never choose rows."""
import numpy as np
import pandas as pd

def compare_rows(pop,reference_rows,alternative_rows):
    """Align exactly one alternative row per precursor to each reference row."""
    reference_rows=np.asarray(reference_rows,dtype=np.int64)
    alternative_rows=np.asarray(alternative_rows,dtype=np.int64)
    p=pop.precursor_id.to_numpy()
    assert len(np.unique(p[alternative_rows]))==len(alternative_rows)
    alt=pd.Series(alternative_rows,index=p[alternative_rows]).reindex(p[reference_rows])
    assert alt.notna().all()
    a=alt.to_numpy(dtype=np.int64)
    rt=pop.rt_apex.to_numpy()
    return a,a!=reference_rows,np.abs(rt[a]-rt[reference_rows])
