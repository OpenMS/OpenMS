Performance Measurements
========================

FloatDataArray
--------------

Raw conversions of floating point data (arrays):

.. code-block:: bash

    python -m timeit -s 'import pyopenms as p, numpy as np;d=np.random.rand(100000).astype(np.float32);f=p.FloatDataArray();f.set_data(d)' 'f.get_data()'
    1000000 loops, best of 3: 1.42 usec per loop

    python -m timeit -s 'import pyopenms as p, numpy as np;d=np.random.rand(100000).astype(np.float32);f=p.FloatDataArray();f.set_data(d)' 'q = [val for val in d]'
    100 loops, best of 3: 3.93 msec per loop

    python -m timeit -s 'import pyopenms, numpy; d = numpy.random.rand( 100000).astype(numpy.float32); f = pyopenms.FloatDataArray();' 'f.set_data(d)'    
    100000 loops, best of 3: 13.6 usec per loop

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
  1000 loops, best of 3: 792 usec per loop

  python -m timeit -s 'import pyopenms, numpy; d = numpy.random.rand( 100000).astype(numpy.float64); df = d.astype(numpy.float32); s=pyopenms.MSSpectrum();' 's._set_peaks_fast_df( d,df )' 
  1000 loops, best of 3: 662 usec per loop

Strangely, using numpy arrays to set data is considerably slower when using the list interface:


.. code-block:: bash

  python -m timeit -s 'import pyopenms, numpy; d = numpy.random.rand( 100000).astype(numpy.float64); df = d.astype(numpy.float32); s=pyopenms.MSSpectrum();' 's._set_peaks_orig(d,df)'

  100 loops, best of 3: 6.46 msec per loop -- version 2.3.0.4 [use 's.set_peaks([d,df])']
  100 loops, best of 3: 6.23 msec per loop


Spectral data access:


.. code-block:: bash

  python -m timeit -r 3 -s 'import pyopenms, numpy; d = numpy.random.rand( 100000).astype(numpy.float32); s=pyopenms.MSSpectrum(); s.set_peaks([d,d])' 's.get_peaks()'

  1000 loops, best of 3: 207 usec per loop -- version 2.3.0.4
  10000 loops, best of 3: 168 usec per loop


