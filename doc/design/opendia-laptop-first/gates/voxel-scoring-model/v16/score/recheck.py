"""Reporting-only recheck of preserved failed reconstruction, without evaluation labels."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from helpers import winners
O=Path(__file__).resolve().parent
P=pd.read_csv(O.parents[2]/'raw01/population.tsv',sep='\t')
A=np.load(O/'raw01/reconstructed.npz')
ref=np.load(O.parents[2]/'raw01/lin_D11F.npy')
a=winners(P,A['linear']);b=winners(P,ref)
x=json.loads((O/'raw01/verification.json').read_text())
x['W_row_agreement']=float(np.isin(a,b).mean())
x['note']='Corrected row-overlap reporting only. Original gate failure is unchanged.'
(O/'raw01/verification_corrected.json').write_text(json.dumps(x,indent=2))
print(json.dumps(x,indent=2))
