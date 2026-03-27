Reading Raw MS Data
===========================

:term:`mzML` Files in Memory
****************************

As discussed in the last section, the most straight forward way to load :term:`mass
spectrometry` data is using the :py:class:`~.MzMLFile` class:

.. code-block:: python

    import pyopenms as oms
    from urllib.request import urlretrieve

    gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-docs/master"
    urlretrieve(gh + "/src/data/tiny.mzML", "test.mzML")
    exp = oms.MSExperiment()
    oms.MzMLFile().load("test.mzML", exp)

which will load the content of the "test.mzML" file into the ``exp``
variable of type :py:class:`~.MSExperiment`. We can access the raw data and spectra through:

.. code-block:: python

    spectrum_data = exp.getSpectrum(0).get_peaks()
    chromatogram_data = exp.getChromatogram(0).get_peaks()

Which will allow us to compute on spectra and chromatogram data. We can
manipulate the spectra in the file for example as follows:

.. code-block:: python

    spec = []
    for s in exp.getSpectra():
        if s.getMSLevel() != 1:
            spec.append(s)

    exp.setSpectra(spec)

Which will only keep :term:`MS2` spectra in the :py:class:`~.MSExperiment`. We can then store the modified data structure on disk:

.. code-block:: python

    oms.MzMLFile().store("filtered.mzML", exp)

Putting this together, a small filtering program would look like this:

.. code-block:: python

    """
    Script to read mzML data and filter out all MS1 spectra
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load("test.mzML", exp)

    spec = []
    for s in exp.getSpectra():
        if s.getMSLevel() != 1:
            spec.append(s)

    exp.setSpectra(spec)

    oms.MzMLFile().store("filtered.mzML", exp)

Indexed :term:`mzML` Files
**************************

Since pyOpenMS 2.4, you can open, read and inspect files that use the
indexedMzML standard. This allows users to read MS data without loading all
data into memory:

.. code-block:: python

    od_exp = oms.OnDiscMSExperiment()
    od_exp.openFile("test.mzML")
    meta_data = od_exp.getMetaData()
    meta_data.getNrChromatograms()
    od_exp.getNrChromatograms()

    # data is not present in meta_data experiment
    sum(meta_data.getChromatogram(0).get_peaks()[1])  # no data!
    sum(od_exp.getChromatogram(0).get_peaks()[1])  # data is here!

    # meta data is present and identical in both data structures:
    meta_data.getChromatogram(0).getNativeID()  # fast
    od_exp.getChromatogram(0).getNativeID()  # slow

Note that the :py:class:`~.OnDiscMSExperiment` allows users to access meta data through
the :py:meth:`~.OnDiscMSExperiment.getMetaData` function, which allows easy selection and filtering on meta
data attributes (such as MS level, precursor m/z, retention time etc.) in
order to select spectra and chromatograms for analysis. Only once selection on
the meta data has been performed, will actual data be loaded into memory using
the :py:meth:`~.OnDiscMSExperiment.getChromatogram` and :py:meth:`~.OnDiscMSExperiment.getSpectrum` functions.

This approach is memory efficient in cases where computation should only occur
on part of the data or the whole data may not fit into memory.

:term:`mzML` Files as Streams
*****************************

In some instances it is impossible or inconvenient to load all data from an
mzML file directly into memory. OpenMS offers streaming-based access to mass
spectrometric data which uses a callback object that receives spectra and
chromatograms as they are read from the disk. A simple implementation could look like

.. code-block:: python

    class MSCallback:
        def setExperimentalSettings(self, s):
            pass

        def setExpectedSize(self, a, b):
            pass

        def consumeChromatogram(self, c):
            print("Read a chromatogram")

        def consumeSpectrum(self, s):
            print("Read a spectrum")


which can the be used as follows:

.. code-block:: python

    filename = b"test.mzML"
    consumer = MSCallback()
    oms.MzMLFile().transform(filename, consumer)
	
.. code-block:: output
	
    Read a spectrum
    Read a spectrum
    Read a spectrum
    Read a spectrum
    Read a chromatogram
    Read a chromatogram

which provides an intuition on how the callback object works: whenever a
spectrum or chromatogram is read from disk, the function ``consumeSpectrum`` or
``consumeChromatogram`` is called and a specific action is performed. We can
use this to implement a simple filtering function for mass spectra:

.. code-block:: python

	import os
	import pyopenms as oms
	from urllib.request import urlretrieve

	gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-docs/master"
	urlretrieve(gh + "/src/data/tiny.mzML", "test.mzML")

	print("Current Working Directory where all files are stored:", os.getcwd())

	class FilteringConsumer:
		"""
		Consumer that forwards all calls the internal consumer (after
		filtering out spectra with less than 'min_spec_size' peaks)
		"""

		def __init__(self, consumer, min_spec_size = 0):
			self._internal_consumer = consumer
			self._min_spec_size = min_spec_size

		def setExperimentalSettings(self, s):
			self._internal_consumer.setExperimentalSettings(s)

		def setExpectedSize(self, a, b):
			self._internal_consumer.setExpectedSize(a, b)

		def consumeChromatogram(self, c):
			# just forward; do nothing to chromatograms
			self._internal_consumer.consumeChromatogram(c)

		def consumeSpectrum(self, s):
			print("Spec has size: ", s.size())
			if s.size() >= self._min_spec_size:
				print(" --> keep it")
				self._internal_consumer.consumeSpectrum(s)
			else:
				print(" --> discard it")


	min_spec_size = 11  ## we will keep spectra with 11 or more peaks and discard the others
	inputfile = "test.mzML"
	outputfile = "out.mzML"

	consumer = oms.PlainMSDataWritingConsumer(outputfile)
	consumer = FilteringConsumer(consumer, min_spec_size)

	oms.MzMLFile().transform(inputfile, consumer)

, where the spectra are filtered by their size. It is
similarly trivial to implement filtering by other attributes. Note how the data
are written to disk using the :py:class:`~.PlainMSDataWritingConsumer` which is one of
multiple available consumer classes -- this specific class will simply take the
:py:class:`~.MSSpectrum` ``s`` or :py:class:`~.MSChromatogram` ``c`` and write it to disk (the location of the
output file is given by the ``outputfile`` variable).

Note that this approach is memory efficient in cases where computation should
only occur on part of the data or the whole data may not fit into memory.


Cached :term:`mzML` Files
*************************

In addition, since pyOpenMS 2.4 the user can efficiently cache :term:`mzML` files to disk which
provides very fast access with minimal overhead in memory. Basically the data
directly mapped into memory when requested. You can use this feature as follows:

.. code-block:: python

    # First load data and cache to disk
    exp = oms.MSExperiment()
    oms.MzMLFile().load("test.mzML", exp)
    oms.CachedmzML().store("myCache.mzML", exp)

    # Now load data
    cfile = oms.CachedmzML()
    oms.CachedmzML().load("myCache.mzML", cfile)

    meta_data = cfile.getMetaData()
    cfile.getNrChromatograms()
    cfile.getNrSpectra()

    # data is not present in meta_data experiment
    sum(meta_data.getChromatogram(0).get_peaks()[1])  # no data!
    sum(cfile.getChromatogram(0).get_peaks()[1])  # data is here!

    # meta data is present and identical in both data structures:
    meta_data.getChromatogram(0).getNativeID()  # fast
    cfile.getChromatogram(0).getNativeID()  # slow

Note that the :py:class:`~.CachedmzML` allows users to access meta data through
the :py:meth:`~.CachedmzML.getMetaData` function, which allows easy selection and filtering on meta
data attributes (such as MS level, precursor m/z, retention time etc.) in
order to select spectra and chromatograms for analysis. Only once selection on
the meta data has been performed, will actual data be loaded into memory using
the :py:meth:`~.CachedmzML.getChromatogram` and :py:meth:`~.CachedmzML.getSpectrum` functions.

Note that in the example above all data is loaded into memory first and then
cached to disk. This is not very efficient and we can use the
:py:class:`~.MSDataCachedConsumer` to directly cache to disk (without loading any data
into memory):

.. code-block:: python

    # First cache to disk
    # Note: writing meta data to myCache2.mzML is required
    cacher = oms.MSDataCachedConsumer("myCache2.mzML.cached")
    exp = oms.MSExperiment()
    oms.MzMLFile().transform(b"test.mzML", cacher, exp)
    oms.CachedMzMLHandler().writeMetadata(exp, "myCache2.mzML")
    del cacher

    # Now load data
    cfile = oms.CachedmzML()
    oms.CachedmzML().load("myCache2.mzML", cfile)

    meta_data = cfile.getMetaData()
    # data is not present in meta_data experiment
    sum(meta_data.getChromatogram(0).get_peaks()[1])  # no data!
    sum(cfile.getChromatogram(0).get_peaks()[1])  # data is here!

This approach is now memory efficient in cases where computation should only occur
on part of the data or the whole data may not fit into memory.
