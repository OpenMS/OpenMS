Mass Spectrometry Imaging
=========================

pyOpenMS reads and writes imzML, the standard format for mass spectrometry imaging data. An imzML data set
consists of two files with the same name: the ``.imzML`` file holds the metadata and the position of every pixel,
the ``.ibd`` file the peak data. Keep both in the same directory.

Loading an imzML file
*********************

:py:class:`~.ImzMLFile` loads a data set into an :py:class:`~.MSImagingExperiment`, which holds one spectrum per
pixel and the geometry of the image. ``FileHandler().loadImagingExperiment()`` loads imaging data too and picks
the reader by the type of the file; ``FileHandler().loadExperiment()`` does not read imzML.

The geometry is two-dimensional and covers the first plane (``z = 1``) of a data set. Spectra of other planes are
loaded too, and ``getMSExperiment()`` returns them, but pixel access, regions and ion images leave them out.

.. code-block:: python
    :linenos:

    import pyopenms as oms
    from urllib.request import urlretrieve

    url = "https://raw.githubusercontent.com/OpenMS/OpenMS/develop/src/tests/class_tests/openms/data/"
    for ext in (".imzML", ".ibd"):
        urlretrieve(url + "ImzMLFile_1_Example_Continuous" + ext, "example" + ext)

    img = oms.MSImagingExperiment()
    oms.ImzMLFile().load("example.imzML", img)
    geometry = img.getGeometry()
    print(geometry.getWidth(), geometry.getHeight(), img.getNumberOfSpectra())

.. code-block:: output

    3 3 9

Pixel coordinates start at 0. The spectrum of the top-left pixel is:

.. code-block:: python
    :linenos:

    spectrum = img.getSpectrum(0, 0)
    print(spectrum.size())

.. code-block:: output

    8399

The experiment is also a container of its pixels: ``len(img)`` is the number of pixels, ``img[i]`` the spectrum of
pixel ``i`` (negative indices count from the end), and a loop visits the spectrum of every pixel:

.. code-block:: python
    :linenos:

    print(len(img))
    sizes = [spectrum.size() for spectrum in img]

.. code-block:: output

    9

The index ``i`` follows the pixels of the geometry, not the spectra of the file. The two differ where the geometry
leaves out or reorders spectra, e.g. those of other planes. Use ``img.experiment_view()[k]`` for the spectrum
with index ``k``, e.g. one of those that ``getRegionSpectrumIndices()`` returns.

Copies and views
****************

``img[i]``, ``getSpectrum(x, y)``, ``getGeometry()`` and ``getMSExperiment()`` return copies. A change to a copy
has to be written back with ``img[i] = spectrum``, ``setSpectrum(x, y, spectrum)``, ``setGeometry()`` or
``setMSExperiment()``. The methods whose name ends in ``_view`` return the object inside the experiment instead,
so a change takes effect at once:

.. code-block:: python
    :linenos:

    img.spectrum_view(0, 0).setMSLevel(1)  # the pixel at x = 0, y = 0
    img.spectrum_view(-1).setMSLevel(1)  # the last pixel
    for spectrum in img.iter_spectrum_views():
        spectrum.setMSLevel(1)

``spectrum_views()`` returns the views of all pixels as a list, ``experiment_view()`` the
:py:class:`~.MSExperiment` and ``geometry_view()`` the geometry. A view of a spectrum is valid until spectra are
added to or removed from the experiment, or the experiment is replaced.

Ion images
**********

``extractIonImage`` sums, for every pixel, the intensities within a tolerance (in ppm) around an m/z value. The
result is an :py:class:`~.IonImage`, whose ``get_data()`` returns a copy of the image as a numpy array of shape
(height, width), for example to show it with matplotlib's ``imshow``:

.. code-block:: python
    :linenos:

    ion = img.extractIonImage(450.0, 10.0)
    data = ion.get_data()
    print(data.shape)

.. code-block:: output

    (3, 3)

``data_view()`` returns an array that shares its memory with the image, so writing to it changes the image, e.g.
``ion.data_view()[:] /= ion.data_view().max()``. It does not change which pixels count as present (see
``get_mask()`` below); ``setIntensity(x, y, value)`` does.

.. note::

    In pyOpenMS 3.6.0, ``get_data()`` returned such a shared array. Code that wrote to the result of
    ``get_data()`` has to use ``data_view()`` now.

A region limits the image to part of the pixels. Regions belong to the geometry, so add them to a copy of it and
set the copy:

.. code-block:: python
    :linenos:

    left = oms.MSImagingRegion.rectangle(1, "left", 0, 0, 1, 2)  # id, name, min x, min y, max x, max y
    geometry.addRegion(left)
    img.setGeometry(geometry)
    ion_left = img.extractIonImage(450.0, 10.0, 1)
    print(ion_left.get_mask().sum())

.. code-block:: output

    6

``get_mask()`` tells which pixels of the image belong to the region. ``img.geometry_view().addRegion(left)`` adds
the region without the copy of the geometry.

Large data sets
***************

An :py:class:`~.OnDiscImzMLExperiment` keeps the peak data in the ``.ibd`` file and reads a spectrum only when it
is needed, so it can open data sets larger than the memory. It extracts ion images the same way:

.. code-block:: python
    :linenos:

    on_disc = oms.OnDiscImzMLExperiment()
    on_disc.open("example.imzML")
    ion = on_disc.extractIonImage(450.0, 10.0)
    spectrum = on_disc.getSpectrumAtCoord(1, 1)

``getSpectrumAtCoord`` takes the coordinates as imzML stores them, starting at 1, so ``(1, 1)`` is the pixel
that :py:class:`~.MSImagingExperiment` calls ``(0, 0)``.

Like :py:class:`~.MSImagingExperiment`, it maps only the plane ``z = 1`` to pixels: ``getSpectrumAtCoord`` raises
an error for other values of ``z``, and ``extractIonImage`` leaves those spectra out. ``getSpectrum(i)`` reads any
spectrum by its index, and ``getIndex(i)`` returns its coordinates (``x``, ``y`` and ``z``) without reading the
peaks.

``len(on_disc)``, ``on_disc[i]`` and a loop follow the pixels, as they do for an
:py:class:`~.MSImagingExperiment`, and read each spectrum from the file. ``getNrSpectra()`` and ``getSpectrum(i)``
follow the spectra in the order of the file instead, so ``on_disc[i]`` and ``on_disc.getSpectrum(i)`` are
different spectra where the geometry leaves out or reorders spectra.

Bruker timsTOF MALDI imaging
****************************

:py:class:`~.BrukerTimsImagingFile` reads the ``.d`` directory of a timsTOF MALDI imaging run, or a ``.d.zip``
archive of it, into an :py:class:`~.MSImagingExperiment`. The positions of the pixels come from the
``MaldiFrameInfo`` table of the data set. The reader is part of the pyOpenMS wheels, but not of every build from
source; ``hasattr(oms, "BrukerTimsImagingFile")`` tells whether it is available. The examples read
``my_image.d``; use one of your own data sets instead.

.. code-block:: python
    :linenos:

    print(oms.BrukerTimsImagingFile.isImagingDataset("my_image.d"))
    img = oms.BrukerTimsImagingFile().load("my_image.d")

``load`` raises an error for a ``.d`` directory that is not a MALDI imaging run and for a data set with more than
one plane. A :py:class:`~.BrukerTimsImagingFile.Config` changes how the frames are read: its ``inner_config`` is a
:py:class:`~.BrukerTimsFile.Config` (see `Bruker timsTOF .d <vendor_formats.html>`_), of which ``export_mode``,
``load_ms1`` and ``load_ms2`` are fixed.

.. code-block:: python
    :linenos:

    config = oms.BrukerTimsImagingFile.Config()
    config.inner_config.calibrate = True  # recalibrate m/z
    img = oms.BrukerTimsImagingFile().load("my_image.d", config)

If ``inner_config`` limits the frames, e.g. to a range of frame ids, the pixels of the other frames are missing
from the image.

``FileHandler().loadImagingExperiment("my_image.d", img)`` reads a MALDI imaging run with the default
configuration into an existing :py:class:`~.MSImagingExperiment`.

Writing imzML
*************

``ImzMLFile().store()`` writes an :py:class:`~.MSImagingExperiment` as a pair of ``.imzML`` and ``.ibd`` files:

.. code-block:: python
    :linenos:

    oms.ImzMLFile().store("copy.imzML", img)

The :term:`TOPP tools` and :term:`TOPPView` do not read imzML files yet.
