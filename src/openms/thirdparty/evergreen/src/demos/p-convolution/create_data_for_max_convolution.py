import numpy as np
from scipy.signal import fftconvolve
import pylab as P
norm = lambda a : a / float(max(a))

N=4096

padded_N = N+9
i=np.arange(padded_N)
x= fftconvolve(np.random.uniform(0.5,1.0,N+9) * ( np.exp( -((i-400)/2500.0)**2 ) + 0.3*np.exp( -((i-padded_N/2.0)/100.0)**2 ) + 0.7*np.exp( -((i-padded_N)/400.0)**2 ) ), [1,2,3,4,4,3,2,1])[8:-8]
y= fftconvolve(np.random.uniform(0.5,1.0,padded_N) * ( 3.0*np.exp( -((i**0.9)/1000.0)**2 ) + np.exp( -((i-padded_N)/600.0)**2 ) ), [1,2,3,4,4,3,2,1])[8:-8]
x = norm(x)
y = norm(y)

outfile=open('x_and_y.txt', 'w')
outfile.write(str(N) + '\n')

for v in x:
  outfile.write(str(v) + ' ')
outfile.write('\n')
for v in y:
  outfile.write(str(v) + ' ')
outfile.write('\n')

