Performance Measurements
========================

FloatDataArray
--------------

Raw conversions of floating point data (arrays):

.. code-block:: bash

    python -m timeit -s 'import pyopenms as p, numpy as np;d=np.random.rand(100000).astype(np.float32);f=p.FloatDataArray();f.set_data(d)' 'np.sum(f.get_data())'
    10000 loops, best of 5: 27.9 usec per loop

    python -m timeit -s 'import pyopenms as p, numpy as np;d=np.random.rand(100000).astype(np.float32);f=p.FloatDataArray();f.set_data(d)' 'np.sum(d)'
    10000 loops, best of 5: 26.4 usec per loop

    python -m timeit -s 'import pyopenms, numpy; d = numpy.random.rand( 100000).astype(numpy.float32); f = pyopenms.FloatDataArray();' 'f.set_data(d)'    
    100000 loops, best of 3: 13.6 usec per loop

    python -m timeit -s 'import pyopenms, numpy; d = numpy.random.rand( 100000).astype(numpy.float32); f = pyopenms.FloatDataArray();' 'q = d.copy()'
    20000 loops, best of 5: 14.5 usec per loop

    python -m timeit -s 'import pyopenms, numpy; d = numpy.random.rand( 100000).astype(numpy.float32); f = pyopenms.FloatDataArray();' 'for val in d: f.push_back(val)'
    100 loops, best of 3: 8.5 msec per loop

OpenSwath Spectrum
------------------

Spectral conversion to OpenSwath Spectrum:

.. code-block:: bash

  python -m timeit -s 'import pyopenms as p, numpy as np;d=np.random.rand(100000);s=p.MSSpectrum();s.set_peaks((d,d));sa=p.OpenSwathDataAccessHelper()' 'sa.convertToSpectrumPtr(s)'

  1000 loops, best of 3: 413 usec per loop -- version 2.3.0.4
  1000 loops, best of 3: 321 usec per loop

Spectral conversion from OpenSwath Spectrum:

.. code-block:: bash

  python -m timeit -s 'import pyopenms as p, numpy as np;d=np.random.rand(100000);s=p.MSSpectrum();s.set_peaks((d,d));sa=p.OpenSwathDataAccessHelper();sptr=sa.convertToSpectrumPtr(s)' 's.clear(False); sa.convertToOpenMSSpectrum(sptr, s)'

  1000 loops, best of 3: 1.45 msec per loop -- version 2.3.0.4
  1000 loops, best of 3: 1.21 msec per loop

Extract data from spectrum

.. code-block:: bash

  python -m timeit -r 3 -s 'import pyopenms, numpy; d = numpy.random.rand( 100000).astype(numpy.float64); s=pyopenms.OSSpectrum();s.setMZArray(d.tolist()); rr= numpy.sum(d)'  \
  'q = s.getMZArray(); assert numpy.sum(q) == rr'
  50 loops, best of 3: 6.19 msec per loop -- version 2.7.0
  2000 loops, best of 3: 111 usec per loop

  'q = s.getMZArray_mv(); assert numpy.sum(q) == rr'
  10000 loops, best of 3: 32.6 usec per loop

  'q = s.getMZArray_mv(); res = q[0] + q[-1]'
  200000 loops, best of 3: 994 nsec per loop

  'arr = s.getDataArrays(); assert numpy.sum(arr[0].getData_mv()) == rr'
  10000 loops, best of 3: 32.7 usec per loop

  'arr = s.getDataArrays(); tmp = arr[0].getData_mv(); res = tmp[0] + tmp[-1]'
  200000 loops, best of 3: 1.16 usec per loop


MSSpectrum
----------

Spectral data setting:

.. code-block:: bash

  python -m timeit -s 'import pyopenms, numpy; d = numpy.random.rand(100000); s=pyopenms.MSSpectrum(); p=pyopenms.Peak1D()' 'for val in d: p.setMZ(val); p.setIntensity(val); s.push_back(p)'
  100 loops, best of 3: 16.8 msec per loop

  python -m timeit -s 'import pyopenms, numpy; d = numpy.random.rand( 100000).astype(numpy.float32); s=pyopenms.MSSpectrum();d=list(d)' 's.set_peaks([d,d])'
  100 loops, best of 3: 2.7 msec per loop -- version 2.3.0.4
  100 loops, best of 3: 2.48 msec per loop

  python -m timeit -s 'import pyopenms, numpy; d = numpy.random.rand( 100000).astype(numpy.float64); df = d.astype(numpy.float32); s=pyopenms.MSSpectrum();' 's._set_peaks_fast_dd( d,d )' 
  1000 loops, best of 5: 367 usec per loop

  python -m timeit -s 'import pyopenms, numpy; d = numpy.random.rand( 100000).astype(numpy.float64); df = d.astype(numpy.float32); s=pyopenms.MSSpectrum();' 's._set_peaks_fast_df( d,df )' 
  1000 loops, best of 5: 369 usec per loop


Strangely, using numpy arrays to set data is considerably slower when using the list interface:


.. code-block:: bash

  python -m timeit -s 'import pyopenms, numpy; d = numpy.random.rand( 100000).astype(numpy.float64); df = d.astype(numpy.float32); s=pyopenms.MSSpectrum();' 's._set_peaks_orig(d,df)'

  100 loops, best of 3: 6.46 msec per loop -- version 2.3.0.4 [use 's.set_peaks([d,df])']
  100 loops, best of 3: 6.23 msec per loop


Spectral data access:


.. code-block:: bash

  python -m timeit -r 3 -s 'import pyopenms, numpy; d = numpy.random.rand( 100000).astype(numpy.float32); s=pyopenms.MSSpectrum(); s.set_peaks([d,d])' 's.get_peaks()'

  1000 loops, best of 3: 207 usec per loop -- version 2.3.0.4
  2000 loops, best of 3: 109 usec per loop



