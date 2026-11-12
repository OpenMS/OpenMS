"""Frozen diagnostic helpers, no evaluation labels used in reconstruction."""
import numpy as np
import pandas as pd

def winners(pop, linear):
    return pop.assign(value=linear, row=np.arange(len(pop))).sort_values(['value','sample_id'], ascending=[False,True]).drop_duplicates('precursor_id').row.to_numpy()

def percentile(values, null):
    null=np.sort(np.asarray(null))
    return (np.searchsorted(null,values,'left')+np.searchsorted(null,values,'right'))/(2*len(null))

def auc(a,b):
    return float(np.mean(percentile(np.asarray(a),np.asarray(b)))) if len(a) and len(b) else np.nan

def membership_category(bits):
    bits=np.asarray(bits,dtype=bool)
    return np.where(~bits.any(axis=1),'never_admitted',np.where(bits[:,-1],'retained_r2','admitted_then_lost'))
