"""Is NEW's raw 02 FDP_c excess beyond what one run's entrapment counts can resolve? Seeds of one raw share the data,
so the sampling unit is the run, not the seed. Under exact calibration at q <= 0.01, entrapments E ~ Poisson(0.01 * N / k),
k = 1 + 1/R ~ 2. Reports P(E >= observed | nominal 1 %) per run (seed 0 and seed median) and pooled over raws (seed 0)."""
from scipy.stats import poisson, chi2
k = 2.0018
runs = {  # (N at q<=0.01, entrapments) per seed, from Q4_FDP_POOLED / EVAL files
 'raw01 new': ([1879,1862,1908,1838,1848],[4,4,7,5,5]), 'raw01 cur': ([1772,1709,1702,1686,1744],[7,7,4,8,4]),
 'raw02 new': ([1920,1949,1965,1943,2010],[12,10,7,7,13]), 'raw02 cur': ([1814,1811,1811,1751,1671],[8,4,3,4,1]),
 'raw03 new': ([1980,1890,1930,1902,1947],[5,7,5,6,7]), 'raw03 cur': ([1876,1732,1757,1902,1849],[9,5,6,8,7])}
import numpy as np
for name,(N,E) in runs.items():
    n=np.median(N); e=np.median(E); mu=0.01*n/k
    lo,hi=chi2.ppf(.05,2*e)/2*k/n, chi2.ppf(.95,2*e+2)/2*k/n
    print(f'{name}: median N {n:.0f}, median E {e:.0f} -> FDP_c {e*k/n:.4f} (one-run 90 % {lo:.4f}-{hi:.4f}); '
          f'expected E at exactly 1 % = {mu:.1f}, P(E >= {e:.0f}) = {poisson.sf(e-1,mu):.2f}')
for arm in ('new','cur'):
    N=sum(runs[f'raw0{r} {arm}'][0][0] for r in (1,2,3)); E=sum(runs[f'raw0{r} {arm}'][1][0] for r in (1,2,3))
    print(f'seed 0 pooled over raws 01-03 ({arm}): E {E}, N {N}, FDP_c {E*k/N:.4f} (90 % {chi2.ppf(.05,2*E)/2*k/N:.4f}-{chi2.ppf(.95,2*E+2)/2*k/N:.4f})')
# Power: with one run, which true FDR could the plan's 'median FDP_c <= 0.010' rule tell apart?
n=1950
for f in (0.005,0.007,0.010,0.013):
    mu=f*n/k; print(f'true FDR {f:.3f}: P(observed FDP_c > 0.010) = {poisson.sf(int(0.010*n/k),mu):.2f}')
