import numpy as np
from scipy.signal import fftconvolve
from time import time
import sys

def main(argv):
    if len(argv) != 1:
        print 'Usage: numpy_benchmark.py <LOG_N>'
        return;

    log_n = int(argv[0])
    n = 2**log_n

    x=np.arange(n)*(1+1j)
    y=np.arange(n)*(-1-1j)
    t1=time()
    z = fftconvolve(x,y)
    t2=time()

if __name__ == '__main__':
    main(sys.argv[1:])

