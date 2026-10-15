Threads and Parallel Processing
===============================

pyOpenMS can use several CPU cores in two ways: many OpenMS algorithms spread
their work over several threads by themselves, and you can run your own code in
several Python threads or processes. This page explains how the two interact.

Parallel Work Inside pyOpenMS
*****************************

Many OpenMS classes use OpenMP to spread a single call over several threads;
for example, :py:meth:`~.MzMLFile.load` decodes the spectra of a file in
parallel, and parts of :term:`feature` detection and peptide search run in
parallel, too. :py:class:`~.OpenMSBuildInfo` tells you whether pyOpenMS was
built with OpenMP and sets the number of threads:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    print(oms.OpenMSBuildInfo.isOpenMPEnabled())
    oms.OpenMSBuildInfo.setOpenMPNumThreads(2)
    print(oms.OpenMSBuildInfo.getOpenMPMaxNumThreads())

.. code-block:: output

    True
    2

:py:meth:`~.OpenMSBuildInfo.getOpenMPMaxNumThreads` returns how many threads a
parallel section started from the calling thread may use (functions with their
own ``threads`` argument can differ). Initially, this is the value of the
environment variable ``OMP_NUM_THREADS`` if it is set before pyOpenMS is
imported, and otherwise usually the number of logical CPU cores.
:py:meth:`~.OpenMSBuildInfo.setOpenMPNumThreads` changes it only for the
calling thread; a thread started later begins with the initial value.

Python Threads and the GIL
**************************

Python's global interpreter lock (GIL) lets only one thread at a time run
Python code, and a C++ function called from Python holds the GIL unless it
releases it. pyOpenMS releases the GIL during some long-running calls: file
input and output (for example :py:meth:`~.MzMLFile.load` and
:py:meth:`~.MzMLFile.store`), the Arrow export and import, sorting an
:py:class:`~.MSExperiment`, and some algorithms for :term:`feature` detection,
deconvolution, peptide indexing and search. Other Python threads run meanwhile;
several threads that load files work in parallel. All other pyOpenMS calls hold
the GIL.

Releasing the GIL does **not** make a call safe to run at the same time as
other calls on the same data. Unless the documentation of a class states
otherwise:

* Do not use a pyOpenMS object from two threads at the same time, whether it
  holds data (:py:class:`~.MSExperiment`) or is an algorithm or file object
  (:py:class:`~.PeakPickerHiRes`, :py:class:`~.MzMLFile`).
* Do not read or change an object, or anything inside it, while another thread
  runs a call that can change it; for example, do not iterate over an
  :py:class:`~.MSExperiment` while another thread loads a file into it.
* Give each thread its own algorithm and data objects. To hand data to another
  thread, pass a copy (``copy.deepcopy(obj)``) or guard it with a
  ``threading.Lock``.
* Views (results of methods ending in ``_view``, ``_views`` or ``_struct``, and
  Arrow tables or buffers that point into OpenMS storage) share memory with
  their object: treat them like that object for as long as they exist.

Avoiding Oversubscription
*************************

If several workers run OpenMP-parallel code at the same time, the thread
counts multiply: on 8 cores, 4 workers with 8 OpenMP threads each run 32
threads. Give each worker about as many OpenMP threads as there are cores per
worker, and set this in each worker (for example in the ``initializer`` of the
pool), since the setting applies only to the calling thread:

.. code-block:: python
    :linenos:

    import os
    from concurrent.futures import ThreadPoolExecutor
    from urllib.request import urlretrieve

    N_WORKERS = 2
    OMP_THREADS = max(1, (os.cpu_count() or 1) // N_WORKERS)


    def init_worker():
        oms.OpenMSBuildInfo.setOpenMPNumThreads(OMP_THREADS)


    def count_spectra(path):
        # each call creates its own objects: nothing is shared between threads
        exp = oms.MSExperiment()
        oms.MzMLFile().load(path, exp)
        return path, exp.getNrSpectra()


    gh = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/doc/pyopenms"
    files = ["BSA1.mzML", "FeatureFinderCentroided_1_input.mzML", "tiny.mzML"]
    for name in files:
        urlretrieve(gh + "/src/data/" + name, name)

    with ThreadPoolExecutor(max_workers=N_WORKERS, initializer=init_worker) as pool:
        for path, n_spectra in pool.map(count_spectra, files):
            print(path, n_spectra)

.. code-block:: output

    BSA1.mzML 1684
    FeatureFinderCentroided_1_input.mzML 112
    tiny.mzML 4

Threads suit this example because :py:meth:`~.MzMLFile.load` releases the GIL.
For work that holds the GIL, use a ``concurrent.futures.ProcessPoolExecutor``
with the same ``initializer``. On Linux, pass
``mp_context=multiprocessing.get_context("spawn")``: before Python 3.14 the
default there is ``fork``, and a forked worker can hang in its first
OpenMP-parallel call if the parent process made such calls before. A ``spawn``
worker imports your script, so start the pool from a script and put what the
workers must not repeat (downloads, the pool itself) under
``if __name__ == "__main__":``.

A program that uses a single thread needs none of this. The contract behind
these rules, with an inventory of the calls that release the GIL, is
`THREAD_SAFETY.md <https://github.com/OpenMS/OpenMS/blob/develop/src/pyOpenMS/THREAD_SAFETY.md>`_.
