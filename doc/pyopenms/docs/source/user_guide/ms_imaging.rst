Mass Spectrometry Imaging
=========================

pyOpenMS reads and writes imzML, the standard format for mass spectrometry imaging data. An imzML data set
consists of two files with the same name: the ``.imzML`` file holds the metadata and the position of every pixel,
the ``.ibd`` file the peak data. Keep both in the same directory.

Loading an imzML file
*********************

:py:class:`~.ImzMLFile` loads a data set into an :py:class:`~.MSImagingExperiment`, which holds one spectrum per
pixel and the geometry of the image. :py:meth:`.FileHandler.loadImagingExperiment` does the same for any imaging
format that OpenMS can detect (imzML, and Bruker timsTOF MALDI imaging ``.d`` folders if OpenMS was built with
OpenTIMS support), while :py:meth:`.FileHandler.loadExperiment` rejects imaging files.

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
    oms.ImzMLFile().load("example.imzML", img)  # or: oms.FileHandler().loadImagingExperiment("example.imzML", img)
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

Ion images
**********

``extractIonImage`` sums, for every pixel, the intensities within a tolerance (in ppm) around an m/z value. The
result is an :py:class:`~.IonImage`, whose ``get_data()`` returns the image as a numpy array of shape (height,
width), for example to show it with matplotlib's ``imshow``:

.. code-block:: python
    :linenos:

    ion = img.extractIonImage(450.0, 10.0)
    data = ion.get_data()
    print(data.shape)

.. code-block:: output

    (3, 3)

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

Writing imzML
*************

``ImzMLFile().store()`` writes an :py:class:`~.MSImagingExperiment` as a pair of ``.imzML`` and ``.ibd`` files:

.. code-block:: python
    :linenos:

    oms.ImzMLFile().store("copy.imzML", img)

The :term:`TOPP tools` and :term:`TOPPView` do not read imzML files yet.
