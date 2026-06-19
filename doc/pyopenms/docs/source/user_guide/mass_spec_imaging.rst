Mass Spectrometry Imaging
=========================

Mass spectrometry imaging (MSI) maps the spatial distribution of molecules
across a tissue section or other surface. pyOpenMS provides a dedicated object
model that keeps the spatial geometry decoupled from the raw spectra, making it
easy to extract ion images, define sub-regions, and integrate with NumPy.

Core Classes
------------

* :py:class:`~.MSImagingExperiment` — the top-level container wrapping an
  :py:class:`~.MSExperiment` with its associated :py:class:`~.MSImagingGeometry`.
* :py:class:`~.MSImagingGeometry` — stores the pixel grid dimensions, physical
  pixel size, the (x, y) → spectrum-index mapping, and an optional set of
  named :py:class:`~.MSImagingRegion` overlays.
* :py:class:`~.IonImage` — a 2D intensity buffer returned by
  :py:meth:`~.MSImagingExperiment.extractIonImage`.
* :py:class:`~.MSImagingRegion` — a spatial sub-region footprint (rectangle or
  arbitrary pixel mask) that can be registered on the geometry.

Loading imzML Data
------------------

The standard interchange format for MSI data is imzML. Load it directly into
an :py:class:`~.MSImagingExperiment`:

.. code-block:: python

    import pyopenms as oms

    # ImzMLFile sets the geometry automatically
    exp = oms.MSExperiment()
    oms.ImzMLFile().load("tissue.imzML", exp)

    mie = oms.MSImagingExperiment(exp)
    geom = mie.getGeometry()
    print("Grid:", geom.getWidth(), "×", geom.getHeight())
    print("Acquired pixels:", mie.getNumberOfPixels())

Building a Geometry Manually
-----------------------------

You can also build the geometry from scratch, for example when creating
synthetic datasets or wrapping data from another format:

.. code-block:: python

    import pyopenms as oms

    # Three spectra at known pixel positions
    exp = oms.MSExperiment()
    for mz_vals, intensities in [
        ([499.9, 500.1], [10.0, 20.0]),
        ([490.0],        [5.0]),
        ([500.0],        [100.0]),
    ]:
        s = oms.MSSpectrum()
        s.set_peaks((mz_vals, intensities))
        exp.addSpectrum(s)

    geom = oms.MSImagingGeometry()
    geom.setDimensions(2, 2)            # 2×2 grid
    geom.setPixelSize(25.0, 25.0, "micrometer")
    geom.addPixel(0, 0, 0)             # (col, row, spectrum_index)
    geom.addPixel(1, 0, 1)
    geom.addPixel(0, 1, 2)
    # pixel (1, 1) is intentionally absent

    mie = oms.MSImagingExperiment(exp)
    mie.setGeometry(geom)

Extracting Ion Images
---------------------

:py:meth:`~.MSImagingExperiment.extractIonImage` sums all peaks within a
*±ppm* window around a target m/z across all acquired pixels and returns an
:py:class:`~.IonImage`:

.. code-block:: python

    # 200 ppm window around m/z 500 (covers 499.9 – 500.1)
    img = mie.extractIonImage(500.0, 200.0)

    print("Image size:", img.getWidth(), "×", img.getHeight())
    print("Pixel (0,0) acquired?", img.hasPixel(0, 0))   # True
    print("Pixel (1,1) acquired?", img.hasPixel(1, 1))   # False (gap)

Accessing Ion Image Data with NumPy
-------------------------------------

:py:class:`~.IonImage` exposes its intensity buffer and validity mask as NumPy
arrays:

.. code-block:: python

    import numpy as np

    data = img.get_data()   # shape (height, width), dtype float64, zero-copy
    mask = img.get_mask()   # shape (height, width), dtype uint8, 1 = acquired

    # Apply the mask before downstream analysis
    valid = data * mask.astype(float)

The intensity array is a *zero-copy* view into the C++ buffer — writes are
reflected in the :py:class:`~.IonImage` object immediately.

Pixel Metadata
--------------

The geometry exposes all acquired pixels as a structured NumPy array for
efficient bulk access:

.. code-block:: python

    pix = mie.getGeometry().get_pixels_struct()
    # dtype has fields: x (uint32), y (uint32), spectrum_index (uint64)
    print(pix["x"], pix["y"], pix["spectrum_index"])

Named Sub-Regions
-----------------

:py:class:`~.MSImagingRegion` defines a spatial footprint — either an
axis-aligned rectangle or a per-pixel bitmask — that can be registered on
the geometry. Regions must not overlap; adding an overlapping region raises
an exception. Coordinates are zero-based.

**Rectangle region**

.. code-block:: python

    from pyopenms import MSImagingRegion

    # region id=1, spanning columns 0–1, rows 0–2
    r = MSImagingRegion.rectangle(1, "sample_area", min_x=0, min_y=0, max_x=1, max_y=2)
    print("Shape:", r.getShape())       # Shape.Rectangle
    print("Area (pixels):", r.area())  # 6

**Mask region (from Python list or NumPy array)**

.. code-block:: python

    import numpy as np

    # 2D bool array — height × width, infers dimensions from shape
    footprint = np.array([[True, False],
                           [False, True]], dtype=bool)
    r_mask = MSImagingRegion.fromMask(2, "sparse_area", origin_x=5, origin_y=3,
                                      mask=footprint)
    print("Area (pixels):", r_mask.area())   # 2

**Registering regions and region-based extraction**

.. code-block:: python

    import pyopenms as oms
    from pyopenms import MSImagingRegion

    # (re-use mie from above)
    region = MSImagingRegion.rectangle(1, "left_column", 0, 0, 0, 1)
    mie.getGeometry().addRegion(region)

    # extractIonImage with a region id — only pixels inside the region are
    # included; the image canvas still has the full grid dimensions so spatial
    # context is preserved
    img_region = mie.extractIonImage(500.0, 200.0, 1)

    data = img_region.get_data()
    mask = img_region.get_mask()
    # pixels outside the region are masked out (mask == 0)

**Querying region membership**

.. code-block:: python

    geom = mie.getGeometry()

    NO_REGION = oms.MSImagingGeometry.NO_REGION

    print(geom.regionOf(0, 0))   # 1  (inside region 1)
    print(geom.regionOf(1, 0))   # NO_REGION  (outside all regions)

    # Spectrum indices that fall inside a region
    indices = geom.getRegionSpectrumIndices(1)
    print("Spectra in region 1:", list(indices))

Checking Overlaps
-----------------

Before adding a region you can test for overlaps explicitly:

.. code-block:: python

    r1 = MSImagingRegion.rectangle(1, "A", 0, 0, 3, 3)
    r2 = MSImagingRegion.rectangle(2, "B", 2, 2, 5, 5)  # overlaps r1
    r3 = MSImagingRegion.rectangle(3, "C", 10, 10, 12, 12)  # disjoint

    print(r1.intersects(r2))   # True
    print(r1.intersects(r3))   # False
