Vendor Formats
==============

pyOpenMS reads two vendor formats directly, without a conversion to :term:`mzML`: Thermo Fisher
``.raw`` files and Bruker timsTOF ``.d`` directories. The pyOpenMS wheels include both readers.
A pyOpenMS built from source has them only if OpenMS was built with ``WITH_THERMO_RAW`` and
``WITH_OPENTIMS`` (the default), so code that may run on such a build should check for them:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    print(hasattr(oms, "ThermoRawFile"), hasattr(oms, "BrukerTimsFile"))

The examples below read ``my_run.raw`` and ``my_run.d``; use one of your own files instead.

Thermo Fisher .raw
******************

The reader uses Thermo's RawFileReader library, which runs on .NET. It needs the **.NET 8
runtime** (or newer), which you can install from the `.NET download page
<https://dotnet.microsoft.com/download>`_ or with your package manager. If the runtime is
installed in a non-standard location, set the environment variable ``DOTNET_ROOT`` to the
directory that contains the ``dotnet`` executable and the ``shared`` sub-directory. Everything
else the reader needs comes with the wheel.

:py:class:`~.FileHandler` reads ``.raw`` files like any other format. It applies Thermo's peak
picking, as the TOPP tool FileConverter does by default, so it returns centroided spectra:

.. code-block:: python
    :linenos:

    exp = oms.MSExperiment()
    oms.FileHandler().loadExperiment("my_run.raw", exp)
    print(exp[0].getMSLevel(), exp[0].getType())

.. code-block:: output

    1 SpectrumType.CENTROID

:py:class:`~.ThermoRawFile` returns the spectra as they were acquired, which for most
:term:`Orbitrap` methods means profile :term:`MS1` spectra. Its options set, among other
things, whether it applies Thermo's peak picking:

.. code-block:: python
    :linenos:

    reader = oms.ThermoRawFile()
    acquired = oms.MSExperiment()
    reader.load("my_run.raw", acquired)
    print(acquired[0].getMSLevel(), acquired[0].getType())

    options = reader.getOptions()
    options.centroid = True
    reader.setOptions(options)
    centroided = oms.MSExperiment()
    reader.load("my_run.raw", centroided)
    print(centroided[0].getMSLevel(), centroided[0].getType())

.. code-block:: output

    1 SpectrumType.PROFILE
    1 SpectrumType.CENTROID

The other fields of :py:class:`~.ThermoRawFileOptions` add the charges the instrument assigned
to centroided peaks (``charge_data``), the noise data (``noise_data``) and the traces of other
detectors such as UV (``all_detectors``). ``preserve_trailers``, ``instrument_methods`` and
``checksum``, which are on by default, keep the scan trailers, the instrument method and a SHA-1
checksum of the file.

Bruker timsTOF .d
*****************

:py:class:`~.BrukerTimsFile` reads the ``.d`` directories of timsTOF instruments. It also reads
a ``.d`` directory packed into a ZIP archive (``.d.zip``), which it unpacks into a temporary
directory. ``load`` returns a new :py:class:`~.MSExperiment`:

.. code-block:: python
    :linenos:

    reader = oms.BrukerTimsFile()
    exp = reader.load("my_run.d")

For DDA-PASEF data, the experiment holds one :term:`MS1` spectrum per frame, which stores the ion
mobility (1/K0) of every peak in a float data array, and one :term:`MS2` spectrum per precursor with a
single ion mobility value. For DIA-PASEF data, it holds one :term:`MS2` spectrum per frame and isolation
window, again with the ion mobility of every peak.

.. code-block:: python
    :linenos:

    ms1 = next(s for s in exp if s.getMSLevel() == 1)
    print(ms1.containsIMData(), ms1.getFloatDataArrays()[0].getName())
    ms2 = next(s for s in exp if s.getMSLevel() == 2)
    print(ms2.getDriftTimeUnit())

.. code-block:: output

    True raw inverse reduced ion mobility array
    DriftTimeUnit.VSSC

A :py:class:`~.BrukerTimsFile.Config` changes how the data is read. For a database search, for
example, the :term:`MS1` spectra are not needed:

.. code-block:: python
    :linenos:

    config = oms.BrukerTimsFile.Config()
    config.load_ms1 = False
    ms2_only = reader.load("my_run.d", config)

``export_mode`` selects the layout: ``AUTO`` (the default) detects DDA or DIA, ``SPECTRUM``
forces the per-precursor layout and ``FRAME`` returns the raw frames. Other fields control m/z
recalibration, frame aggregation and centroiding along the ion mobility axis.
``ms2_centroid_algo`` (:py:class:`~.BrukerTimsFile.Config.CentroidAlgo`: ``OFF``, ``GREEDY2D`` or
``HILL_BASED``) selects the MS2 centroiding of DIA-PASEF and DDA-PASEF. ``GREEDY2D`` needs RT-neighbour
aggregation (``dia_ms2_n_neighbors > 0``); ``HILL_BASED`` centroids across the ion mobility scans of each
frame and also works with ``dia_ms2_n_neighbors = 0``.
:py:class:`~.FileHandler` reads ``.d`` directories with the default configuration.

Bruker's SDK is not needed: m/z and ion mobility are computed from the calibration stored in the
file. If Bruker's SDK library (``timsdata.dll`` or ``libtimsdata.so``) is installed, set
``config.bruker_sdk_path``, or the environment variable ``OPENMS_BRUKER_SDK_PATH``, to its path
to convert ion mobility with Bruker's own code.
