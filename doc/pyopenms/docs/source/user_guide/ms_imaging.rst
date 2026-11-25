Mass Spectrometry Imaging
=========================

pyOpenMS reads and writes imzML, the standard format for mass spectrometry imaging data. An imzML data set
consists of two files with the same name: the ``.imzML`` file holds the metadata and the position of every pixel,
the ``.ibd`` file the peak data. Keep both in the same directory.

Loading an imzML file
*********************

:py:class:`~.ImzMLFile` loads a data set into an :py:class:`~.MSImagingExperiment`, which holds one spectrum per
pixel and the geometry of the image. ``FileHandler().loadExperiment()`` does not read imaging data;
``FileHandler().loadImagingExperiment(filename, img)`` does, and picks the reader from the file type.

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

Pixels as a container
*********************

An :py:class:`~.MSImagingExperiment` behaves like a list of its pixels: ``len(img)`` is the number of pixels,
``img[i]`` the spectrum of the i-th pixel (a negative ``i`` counts from the end), and a loop visits the spectra
of all pixels. The order is that of the pixels in the geometry.

.. code-block:: python
    :linenos:

    print(len(img), img[0].size())
    total = sum(spectrum.size() for spectrum in img)

.. code-block:: output

    9 8399

A pixel index is not a spectrum index: the two differ as soon as the geometry leaves out or reorders spectra,
for example those of other planes. ``getNumberOfSpectra()`` counts the spectra, and an experiment without a
geometry has no pixels (``len(img)`` is 0) even if it holds spectra. To reach a spectrum by its spectrum index,
such as one that ``getRegionSpectrumIndices()`` returns, use ``img.experiment_view()[k]``.

``img[i]`` and ``getSpectrum(x, y)`` return copies, as element access does everywhere in pyOpenMS, so a change
to such a spectrum does not reach the experiment unless you write it back (``img[i] = spectrum``). To change the
data in place, use the views, whose names end in ``_view``:

.. code-block:: python
    :linenos:

    img.spectrum_view(0).setMetaValue("checked", "yes")  # by pixel index
    img.spectrum_view(2, 1).setMetaValue("checked", "yes")  # by pixel coordinates
    for spectrum in img.iter_spectrum_views():
        spectrum.setMSLevel(1)

``experiment_view()`` and ``geometry_view()`` are the views of the whole :py:class:`~.MSExperiment` and of the
geometry; ``getMSExperiment()`` and ``getGeometry()`` return copies of them. A spectrum view stays valid only
until the list of spectra is resized or reordered, for example by ``setMSExperiment()`` or by sorting through
``experiment_view()``.

Ion images
**********

``extractIonImage`` sums, for every pixel, the intensities within a tolerance (in ppm) around an m/z value. The
result is an :py:class:`~.IonImage`, whose ``get_data()`` returns the image as a numpy array of shape (height,
width), for example to show it with matplotlib's ``imshow``. The array is a copy; ``data_view()`` returns the
intensities of the image itself without copying them, so that writing to that array changes the image:

.. code-block:: python
    :linenos:

    ion = img.extractIonImage(450.0, 10.0)
    data = ion.get_data()
    print(data.shape)
    ion.data_view()[:] /= max(data.max(), 1.0)  # scales the image in place

.. code-block:: output

    (3, 3)

.. note::

    In pyOpenMS 3.6, ``get_data()`` returned the array that ``data_view()`` returns now. Code that changed an
    image by writing to the result of ``get_data()`` has to use ``data_view()`` instead.

A region limits the image to part of the pixels. Regions belong to the geometry, so add them to a copy of it and
set the copy (or add them in place with ``img.geometry_view().addRegion(left)``):

.. code-block:: python
    :linenos:

    left = oms.MSImagingRegion.rectangle(1, "left", 0, 0, 1, 2)  # id, name, min x, min y, max x, max y
    geometry.addRegion(left)
    img.setGeometry(geometry)
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
spectrum by its index in the file (there are ``getNrSpectra()`` of them), and ``getIndex(i)`` returns its
coordinates (``x``, ``y`` and ``z``) without reading the peaks.

It is a container of pixels as well: ``len(on_disc)``, ``on_disc[i]`` and a loop follow the pixels of the
geometry and read each spectrum from the file when it is asked for. ``on_disc[i]`` therefore takes a pixel index
and ``getSpectrum(i)`` a spectrum index; the two differ if the file holds spectra that are not pixels of the
image.

.. code-block:: python
    :linenos:

    print(len(on_disc), on_disc.getNrSpectra())
    tic = [spectrum.calculateTIC() for spectrum in on_disc]

.. code-block:: output

    9 9

Bruker timsTOF MALDI imaging
****************************

:py:class:`~.BrukerTimsImagingFile` reads the ``.d`` directories (or ``.d.zip`` archives) of timsTOF MALDI
imaging runs. Like the other `vendor format <vendor_formats.html>`_ readers, it is part of the pyOpenMS wheels,
but of a pyOpenMS built from source only if OpenMS was built with ``WITH_OPENTIMS`` (the default);
``hasattr(oms, "BrukerTimsImagingFile")`` tells. The example reads ``my_image.d``; use one of your own runs
instead.

.. code-block:: python
    :linenos:

    if oms.BrukerTimsImagingFile.isImagingDataset("my_image.d"):
        img = oms.BrukerTimsImagingFile().load("my_image.d")
        print(img.getGeometry().getWidth(), img.getGeometry().getHeight(), len(img))

``load`` returns an :py:class:`~.MSImagingExperiment` with one spectrum per frame. The pixel positions come from
the run's metadata and are shifted so that they start at 0. A run with more than one section (several values of
``z``) is rejected.

A :py:class:`~.BrukerTimsImagingFile.Config` changes how the run is read. Its ``inner_config`` is the
:py:class:`~.BrukerTimsFile.Config` of the frame reader (m/z recalibration, centroiding, and so on) and can be
edited in place. If it limits the frames that are read, the pixels of the other frames are left out of the
image. ``strict_imaging_only`` (on by default) rejects runs that were not acquired as imaging runs.

.. code-block:: python
    :linenos:

    config = oms.BrukerTimsImagingFile.Config()
    config.inner_config.calibrate = True
    img = oms.BrukerTimsImagingFile().load("my_image.d", config)

Writing imzML
*************

``ImzMLFile().store()`` writes an :py:class:`~.MSImagingExperiment` as a pair of ``.imzML`` and ``.ibd`` files:

.. code-block:: python
    :linenos:

    oms.ImzMLFile().store("copy.imzML", img)

The :term:`TOPP tools` and :term:`TOPPView` do not read imzML files yet.
