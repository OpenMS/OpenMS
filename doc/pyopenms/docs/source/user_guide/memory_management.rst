Memory Management
==================

On order to save memory, we can avoid loading the whole file into memory and
use the :py:class:`~.OnDiscMSExperiment` for reading data.

.. code-block:: python
  :linenos:

  import pyopenms as oms

  od_exp = oms.OnDiscMSExperiment()
  od_exp.openFile("test.mzML")

  e = oms.MSExperiment()
  for k in range(od_exp.getNrSpectra()):
      s = od_exp.getSpectrum(k)
      if s.getNativeID().startswith("scan="):
          e.addSpectrum(s)

  oms.MzMLFile().store("test_filtered.mzML", e)

Note that using the approach the output data ``e`` is still completely in
memory and may end up using a substantial amount of memory. We can avoid that
by using

.. code-block:: python
  :linenos:

  od_exp = oms.OnDiscMSExperiment()
  od_exp.openFile("test.mzML")

  consumer = oms.PlainMSDataWritingConsumer("test_filtered.mzML")

  e = oms.MSExperiment()
  for k in range(od_exp.getNrSpectra()):
      s = od_exp.getSpectrum(k)
      if s.getNativeID().startswith("scan="):
          consumer.consumeSpectrum(s)

  del consumer

Make sure you do not forget ``del consumer`` since otherwise the final part of
the :term:`mzML` may not get written to disk (and the consumer is still waiting for new
data).

