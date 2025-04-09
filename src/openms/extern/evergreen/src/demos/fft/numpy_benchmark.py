import numpy as np
from time import time
import sys

def main(argv):
    if len(argv) != 1:
        print 'Usage: numpy_benchmark.py <LOG_N>'
        return;

    logN = int(argv[0])
    N = 2**logN

    x=np.arange(N)*(1+1j)
    t1=time()
    y=np.fft.fftn(x)
    t2=time()

    print N, t2-t1

if __name__ == '__main__':
    main(sys.argv[1:])

