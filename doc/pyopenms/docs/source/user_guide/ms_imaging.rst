Mass Spectrometry Imaging
=========================

pyOpenMS reads and writes imzML, the standard format for mass spectrometry imaging data. An imzML data set
consists of two files with the same name: the ``.imzML`` file holds the metadata and the position of every pixel,
the ``.ibd`` file the peak data. Keep both in the same directory. pyOpenMS also reads the MALDI imaging runs of
Bruker timsTOF instruments (see `Bruker timsTOF MALDI imaging`_).

Loading an imzML file
*********************

:py:class:`~.ImzMLFile` loads a data set into an :py:class:`~.MSImagingExperiment`, which holds one spectrum per
pixel and the geometry of the image. ``FileHandler().loadImagingExperiment()`` does the same and detects the
format from the file; ``FileHandler().loadExperiment()`` does not read imaging data.

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

Pixels
******

An :py:class:`~.MSImagingExperiment` is a container of its pixels: ``len()``, indexing and iteration follow the
pixels of the geometry and return a copy of each pixel's spectrum.

.. code-block:: python
    :linenos:

    print(len(img), img[0].size())
    tic = [sum(spectrum.get_peaks()[1]) for spectrum in img]

.. code-block:: output

    9 8399

A change to such a copy does not reach the experiment. To change a spectrum in place, take a live view of it
with ``spectrum_view``, which accepts a pixel index or the coordinates of a pixel; ``spectrum_views()`` returns
the views of all pixels:

.. code-block:: python
    :linenos:

    img.spectrum_view(0, 0).setMetaValue("note", "top left")
    print(img.getSpectrum(0, 0).getMetaValue("note"))

.. code-block:: output

    top left

In the same way, ``geometry_view()`` and ``experiment_view()`` return the geometry and the
:py:class:`~.MSExperiment` of the imaging experiment itself, whereas ``getGeometry()`` and ``getMSExperiment()``
return copies. A view is only valid as long as the spectra it belongs to are not replaced or resized.

The index of a pixel is not always the index of its spectrum: a geometry can leave spectra out (for example those
of other planes) or order them differently. ``experiment_view()[k]`` is the spectrum with index ``k``, and an
experiment without a geometry has no pixels (``len()`` is 0) although it may hold spectra.

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

``data_view()`` returns the same array without copying it. Values written to it change the
:py:class:`~.IonImage`, for example to normalize an image in place (``ion.data_view()[:] /= norm``).

A region limits the image to part of the pixels. Regions belong to the geometry, so add them to the geometry of
the experiment with ``geometry_view()``:

.. code-block:: python
    :linenos:

    left = oms.MSImagingRegion.rectangle(1, "left", 0, 0, 1, 2)  # id, name, min x, min y, max x, max y
    img.geometry_view().addRegion(left)
    ion_left = img.extractIonImage(450.0, 10.0, 1)
    print(ion_left.get_mask().sum())

.. code-block:: output

    6

``get_mask()`` tells which pixels of the image belong to the region.

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

``len()``, indexing and iteration follow the pixels here too, so ``on_disc[i]`` reads the spectrum of the pixel
with index ``i``, while ``getSpectrum(i)`` takes the index of a spectrum in the file and ``getNrSpectra()`` counts
the spectra. The two differ if the geometry leaves spectra out.

Bruker timsTOF MALDI imaging
****************************

:py:class:`~.BrukerTimsImagingFile` loads a MALDI imaging run of a timsTOF instrument (a ``.d`` directory, also
zipped as ``.d.zip``) into an :py:class:`~.MSImagingExperiment`. Like :py:class:`~.BrukerTimsFile` (see
:doc:`vendor_formats`), it is part of the pyOpenMS wheels and needs no Bruker software; a pyOpenMS built from
source has it only if OpenMS was built with ``WITH_OPENTIMS`` (the default). The example reads ``my_image.d``;
use one of your own runs instead.

.. code-block:: python
    :linenos:

    if oms.BrukerTimsImagingFile.isImagingDataset("my_image.d"):
        img = oms.BrukerTimsImagingFile().load("my_image.d")
        print(img.getGeometry().getWidth(), img.getGeometry().getHeight(), len(img))

Every pixel holds the spectrum of one frame. The pixel coordinates start at 0, as for imzML, whatever the
positions stored in the run. Only two-dimensional images are read: a run with several sections (values of ``z``)
raises an error.

A :py:class:`~.BrukerTimsImagingFile.Config` changes how the run is read. Its ``inner_config`` is a
:py:class:`~.BrukerTimsFile.Config`, for example to recalibrate m/z or to centroid along the ion mobility axis.
The frames are always read one by one (``export_mode`` is set to ``FRAME``). If the configuration restricts the
frames, for example with ``frame_id_min`` and ``frame_id_max``, the pixels of the other frames are missing from
the image.

.. code-block:: python
    :linenos:

    config = oms.BrukerTimsImagingFile.Config()
    config.inner_config.frame_id_max = 1000
    part = oms.BrukerTimsImagingFile().load("my_image.d", config)

``readMaldiFrameInfo`` returns the frame and the position of every pixel as stored in the run, without loading
the spectra.

Writing imzML
*************

``ImzMLFile().store()`` writes an :py:class:`~.MSImagingExperiment` as a pair of ``.imzML`` and ``.ibd`` files:

.. code-block:: python
    :linenos:

    oms.ImzMLFile().store("copy.imzML", img)

The :term:`TOPP tools` and :term:`TOPPView` do not read imzML files yet.
