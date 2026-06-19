Mass Spectrometry Imaging (MSI)
================================

pyOpenMS provides classes for working with Mass Spectrometry Imaging (MSI)
datasets, where each pixel of a 2-D sample surface is associated with one
mass spectrum.

Core classes
------------

* :py:class:`~.MSImagingExperiment` — an :py:class:`~.MSExperiment` extended
  with a 2-D pixel grid.
* :py:class:`~.MSImagingGeometry` — the pixel coordinate system: dimensions,
  pixel size, and the mapping from (x, y) to spectrum index.
* :py:class:`~.IonImage` — a 2-D ion-intensity image returned by
  :py:meth:`~.MSImagingExperiment.extractIonImage`.
* :py:class:`~.MSImagingRegion` — a spatial region (rectangle or arbitrary
  mask) that can be annotated on the geometry.

Building an imaging dataset
---------------------------

An MSI dataset is constructed from an :py:class:`~.MSExperiment` (spectra in
acquisition order) paired with a geometry that maps each pixel to its
spectrum:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    # 1. Load spectra from disk (any format accepted by MzMLFile, etc.)
    exp = oms.MSExperiment()
    oms.MzMLFile().load("imaging_run.mzML", exp)

    # 2. Define the pixel grid (width x height) and map pixels to spectra
    geom = oms.MSImagingGeometry()
    geom.setDimensions(10, 10)             # 10 x 10 pixels
    geom.setPixelSize(50.0, 50.0, "micrometer")

    spectrum_index = 0
    for y in range(10):
        for x in range(10):
            geom.addPixel(x, y, spectrum_index)
            spectrum_index += 1

    # 3. Combine into an MSImagingExperiment
    mie = oms.MSImagingExperiment(exp)
    mie.setGeometry(geom)
    print("pixels:", mie.getNumberOfPixels())
    print("spectra:", mie.getNumberOfSpectra())

Extracting ion images
---------------------

:py:meth:`~.MSImagingExperiment.extractIonImage` sums all peaks within a
mass window (centre m/z + ppm tolerance) into a 2-D image:

.. code-block:: python
    :linenos:

    # Extract ion image for m/z 760.6 ± 10 ppm
    img = mie.extractIonImage(760.6, 10.0)

    # Get intensity data as a (height, width) numpy float64 array
    data = img.get_data()   # zero-copy view
    print(data.shape)       # (10, 10)

    # Boolean mask: 1 where a pixel was acquired, 0 where it is absent
    mask = img.get_mask()   # copy
    print(mask.dtype)       # uint8

    # Individual pixel access
    if img.hasPixel(3, 2):
        print(img.getIntensity(3, 2))

    # m/z window actually applied
    mz_range = img.getMzRange()
    print(mz_range.getMinMZ(), mz_range.getMaxMZ())

Because ``get_data()`` returns a zero-copy NumPy view, you can write into it
and the change is reflected in the C++ object immediately.

Defining spatial regions
------------------------

:py:class:`~.MSImagingRegion` lets you annotate named sub-areas of the
sample.  Two shapes are supported: an axis-aligned rectangle and an arbitrary
per-pixel bitmask.

**Rectangle region**

.. code-block:: python
    :linenos:

    import pyopenms as oms

    # rectangle(id, name, min_x, min_y, max_x, max_y)  — all inclusive, zero-based
    region = oms.MSImagingRegion.rectangle(1, "tumour", 2, 3, 6, 8)

    print(region.getId())           # 1
    print(region.getName())         # "tumour"
    print(region.getBBoxWidth())    # 5  (columns 2..6)
    print(region.getBBoxHeight())   # 6  (rows 3..8)
    print(region.area())            # 30

    print(region.contains(4, 5))   # True  — inside
    print(region.contains(0, 0))   # False — outside

**Mask region (arbitrary shape)**

.. code-block:: python
    :linenos:

    import numpy as np
    import pyopenms as oms

    # 2-D bool array — shape is (height, width); True marks pixels inside the region
    footprint = np.array([
        [True,  False],
        [False, True ],
    ], dtype=bool)

    # fromMask(id, name, origin_x, origin_y, mask_2d)
    region = oms.MSImagingRegion.fromMask(2, "stroma", 5, 10, footprint)
    print(region.area())                # 2 (two True cells)
    print(region.contains(5, 10))      # True  — (origin_x, origin_y) is set
    print(region.contains(6, 10))      # False — bit [0,1] is False

Adding regions to the geometry
------------------------------

Regions are attached to the :py:class:`~.MSImagingGeometry`.  The geometry
rejects duplicate IDs and overlapping footprints:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    geom = mie.getGeometry()                       # from the experiment above

    r1 = oms.MSImagingRegion.rectangle(1, "region_A", 0, 0, 4, 4)
    r2 = oms.MSImagingRegion.rectangle(2, "region_B", 5, 0, 9, 4)
    geom.addRegion(r1)
    geom.addRegion(r2)

    print(geom.getNumberOfRegions())               # 2
    print(geom.regionOf(2, 2))                     # 1  — belongs to region_A
    print(geom.regionOf(7, 3))                     # 2  — belongs to region_B
    print(geom.regionOf(0, 9))                     # MSImagingGeometry.NO_REGION

    # Retrieve spectrum indices that fall inside a region
    indices = geom.getRegionSpectrumIndices(1)
    print(list(indices))                           # spectrum indices for region_A

    # Remove a region (leaves acquired pixels untouched)
    geom.removeRegion(1)

Extracting ion images for a region
-----------------------------------

Pass the region ID as a third argument to
:py:meth:`~.MSImagingExperiment.extractIonImage` to restrict the image to
that region.  Pixels outside the region are marked absent in the mask but the
image dimensions still match the full dataset so coordinates stay consistent:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    # Attach the region if not done already
    r = oms.MSImagingRegion.rectangle(1, "roi", 0, 0, 4, 4)
    mie.getGeometry().addRegion(r)

    # Two-argument form: whole dataset
    img_full   = mie.extractIonImage(760.6, 10.0)

    # Three-argument form: region 1 only
    img_region = mie.extractIonImage(760.6, 10.0, 1)

    mask = img_region.get_mask()
    # Pixels outside region 1 are 0 even if they were acquired
    data = img_region.get_data()
