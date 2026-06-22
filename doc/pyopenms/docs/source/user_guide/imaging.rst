MS Imaging
==========

Mass spectrometry imaging (MSI) records a full mass spectrum at each spatial
position on a tissue or sample surface, producing a three-dimensional dataset
(x, y, m/z). pyOpenMS represents this as an
:py:class:`~.MSImagingExperiment`, which combines an
:py:class:`~.MSExperiment` (holding one spectrum per pixel) with an
:py:class:`~.MSImagingGeometry` (describing the spatial layout and any defined
regions of interest).

Building an MSImagingExperiment Programmatically
*************************************************

Each acquired pixel is represented by one :py:class:`~.MSSpectrum` in an
:py:class:`~.MSExperiment`. The spatial mapping from pixel coordinates
``(x, y)`` to spectrum index is held by :py:class:`~.MSImagingGeometry` via
:py:meth:`~.MSImagingGeometry.addPixel`. The two are combined into a single
:py:class:`~.MSImagingExperiment`:

.. code-block:: python
    :linenos:

    import pyopenms as oms
    import numpy as np

    width, height = 3, 2   # 3 columns (x), 2 rows (y) -> 6 pixels

    # --- spectra (one per pixel, in any order) ---
    exp = oms.MSExperiment()
    mz_axis = np.array([100.0, 200.0, 300.0], dtype=np.float64)
    for y in range(height):
        for x in range(width):
            spec = oms.MSSpectrum()
            intensities = np.array([float(x + 1), float(y + 1), float(x + y + 1)],
                                   dtype=np.float64)
            spec.set_peaks((mz_axis, intensities))
            exp.addSpectrum(spec)

    # --- geometry (pixel grid + pixel-to-spectrum index mapping) ---
    geom = oms.MSImagingGeometry()
    geom.setDimensions(width, height)
    geom.setPixelSize(25.0, 25.0, "micrometer")   # optional physical size

    spectrum_index = 0
    for y in range(height):
        for x in range(width):
            geom.addPixel(x, y, spectrum_index)
            spectrum_index += 1

    # --- assemble ---
    msi = oms.MSImagingExperiment(exp)
    msi.setGeometry(geom)

    print("Spectra:", msi.getNumberOfSpectra())   # 6
    print("Pixels: ", msi.getNumberOfPixels())    # 6
    print("Grid:   ", geom.getWidth(), "x", geom.getHeight())


Loading imzML Files
-------------------

Real MSI data is typically stored as `imzML
<https://imzml.org>`_ and loaded directly from disk.
:py:class:`~.ImzMLFile` loads directly into an
:py:class:`~.MSImagingExperiment` (geometry is built automatically):

.. code-block:: python
    :linenos:

    import pyopenms as oms

    msi = oms.MSImagingExperiment()
    oms.ImzMLFile().load("sample.imzML", msi)
    print("Loaded", msi.getNumberOfPixels(), "pixels")


Extracting Ion Images
*********************

Full-Grid Ion Image
-------------------

:py:meth:`~.MSImagingExperiment.extractIonImage` collects, for every acquired
pixel, the summed intensity within ``tolerance_ppm`` of the target m/z and
returns an :py:class:`~.IonImage` object. Call ``get_data()`` for a 2-D
``float64`` numpy array with shape ``(height, width)``, or ``get_mask()`` for a
``uint8`` array where ``1`` marks acquired pixels:

.. code-block:: python
    :linenos:

    import pyopenms as oms
    import numpy as np

    # Build 3x2 experiment from the code above (msi variable)

    ion_img = msi.extractIonImage(200.0, 10.0)   # m/z = 200, 10 ppm window

    data = ion_img.get_data()    # shape (2, 3), float64
    mask = ion_img.get_mask()    # shape (2, 3), uint8 (1 = acquired pixel)
    print("Image shape:", data.shape)
    print("Ion image:\n", data)

    # Access a single pixel
    print("Pixel (1,0) intensity:", ion_img.getIntensity(1, 0))

Region-Filtered Ion Image
--------------------------

When a region of interest has been defined (see below), you can restrict
extraction to pixels inside that region by supplying its integer id as a third
argument.  Pixels outside the region are zero in the data array and zero in the
mask:

.. code-block:: python
    :linenos:

    # Region 1 must be added to the geometry before extraction (see below)
    ion_img_roi = msi.extractIonImage(200.0, 10.0, 1)
    data_roi = ion_img_roi.get_data()
    mask_roi = ion_img_roi.get_mask()   # 0 outside the region, 1 inside
    print("ROI image:\n", data_roi)


Accessing Pixel Spectra
***********************

:py:meth:`~.MSImagingExperiment.getSpectrum` retrieves the spectrum at pixel
coordinates ``(x, y)``:

.. code-block:: python
    :linenos:

    spec = msi.getSpectrum(0, 0)   # pixel at column 0, row 0
    mz, intensity = spec.get_peaks()
    print("Spectrum at (0,0):", mz, intensity)

    # Check whether a pixel was acquired
    print("Has pixel (2, 1):", msi.hasPixel(2, 1))


Defining Regions of Interest
*****************************

An :py:class:`~.MSImagingRegion` describes a spatial sub-area of the sample.
Two shapes are supported: axis-aligned rectangles and arbitrary bitmasks.
Regions are managed through :py:class:`~.MSImagingGeometry` and may not
overlap each other.

Rectangular Regions
-------------------

Use the :py:meth:`~.MSImagingRegion.rectangle` factory to define a rectangular
region from its bounding-box coordinates (all coordinates are inclusive,
zero-based):

.. code-block:: python
    :linenos:

    import pyopenms as oms

    # rectangle(id, name, min_x, min_y, max_x, max_y)
    rect_region = oms.MSImagingRegion.rectangle(1, "tumor_core", 1, 0, 2, 1)

    print("Shape:", rect_region.getShape())   # MSImagingRegion.Shape.Rectangle
    print("BBox:", rect_region.getMinX(), rect_region.getMinY(),
          "->", rect_region.getMaxX(), rect_region.getMaxY())
    print("Width x Height:", rect_region.getBBoxWidth(), "x", rect_region.getBBoxHeight())
    print("Area:", rect_region.area())                 # 4 pixels (2 x 2)
    print("Contains (2,1):", rect_region.contains(2, 1))   # True
    print("Contains (0,0):", rect_region.contains(0, 0))   # False

Bitmask Regions
---------------

Irregular shapes are described via a boolean mask.  The mask covers the
bounding box; a ``True`` / ``1`` entry marks pixels belonging to the region.
``origin_x`` / ``origin_y`` give the top-left corner of the bounding box in
global coordinates.

:py:meth:`~.MSImagingRegion.fromMask` accepts a Python list, a 1-D numpy
``bool`` array (with explicit ``width`` and ``height``), or a 2-D numpy
``bool`` array (width and height are inferred from the shape):

.. code-block:: python
    :linenos:

    import pyopenms as oms
    import numpy as np

    # L-shaped region as a 3x3 bitmask (row-major: rows first, then columns)
    mask_2d = np.array([
        [True, True,  True],
        [True, False, False],
        [True, False, False],
    ], dtype=bool)

    # fromMask(id, name, origin_x, origin_y, mask_2d)
    # width and height are inferred from mask_2d.shape
    mask_region = oms.MSImagingRegion.fromMask(2, "margin", 0, 0, mask_2d)

    print("Shape:", mask_region.getShape())   # MSImagingRegion.Shape.Mask
    print("BBox:", mask_region.getBBoxWidth(), "x", mask_region.getBBoxHeight())
    print("Area:", mask_region.area())        # 5 pixels (L-shape)

    # Retrieve the mask as a 2-D uint8 numpy array (height x width)
    retrieved = mask_region.get_mask()
    print("Mask array (1=inside):\n", retrieved)


Registering Regions with the Geometry
--------------------------------------

Regions are added to the :py:class:`~.MSImagingGeometry` attached to the
experiment.  Two regions may not share the same integer id, and their
footprints must not overlap; both violations raise an exception:

.. code-block:: python
    :linenos:

    import pyopenms as oms
    import numpy as np

    # Re-use the msi / geom built earlier
    region_a = oms.MSImagingRegion.rectangle(1, "left_half",  0, 0, 1, 1)
    region_b = oms.MSImagingRegion.rectangle(2, "right_half", 2, 0, 2, 1)

    geom.addRegion(region_a)
    geom.addRegion(region_b)

    print("Regions:", geom.getNumberOfRegions())              # 2
    print("Names:", [r.getName() for r in geom.getRegions()])

    # Look up a region by id
    r = geom.getRegion(1)
    print("Region 1 name:", r.getName())

    # Check which region a pixel belongs to
    # Returns the region id, or MSImagingGeometry.NO_REGION when outside all regions
    rid = geom.regionOf(0, 0)
    if rid != oms.MSImagingGeometry.NO_REGION:
        print("Pixel (0,0) is in region:", rid)        # 1

    rid2 = geom.regionOf(1, 0)
    print("Pixel (1,0) is in region:", rid2)           # 1

    rid3 = geom.regionOf(2, 0)
    print("Pixel (2,0) is in region:", rid3)           # 2


Querying Pixels and Spectra of a Region
*****************************************

Once regions are registered you can retrieve the pixel coordinates or spectrum
indices that fall inside a region.

:py:meth:`~.MSImagingExperiment.getRegionSpectrumIndices` returns spectrum
indices directly, suitable for accessing spectra via
:py:meth:`~.MSExperiment.getSpectrum`:

.. code-block:: python
    :linenos:

    # Spectrum indices belonging to region 1
    indices = msi.getRegionSpectrumIndices(1)
    print("Spectrum indices for region 1:", list(indices))

    for idx in indices:
        spec = msi.getMSExperiment().getSpectrum(idx)
        mz, intensity = spec.get_peaks()
        print(f"  spectrum {idx}: total intensity = {intensity.sum():.1f}")

:py:meth:`~.MSImagingGeometry.getRegionPixels` returns indices into the
geometry's pixel list (from :py:meth:`~.MSImagingGeometry.get_pixels_struct`),
which lets you recover the ``(x, y)`` coordinates:

.. code-block:: python
    :linenos:

    # Pixel positions belonging to region 1
    pix_all = geom.get_pixels_struct()     # structured numpy array: x, y, spectrum_index
    pix_indices = geom.getRegionPixels(1)  # indices into pix_all
    for i in pix_indices:
        x, y = int(pix_all['x'][i]), int(pix_all['y'][i])
        spec = msi.getSpectrum(x, y)
        mz, intensity = spec.get_peaks()
        print(f"  pixel ({x},{y}) total intensity = {intensity.sum():.1f}")

    # Remove a region (spectral data is unaffected)
    geom.removeRegion(2)
    print("Regions after removal:", geom.getNumberOfRegions())   # 1

    # Remove all regions at once
    geom.clearRegions()
    print("Regions after clearRegions:", geom.getNumberOfRegions())  # 0


Checking Region Relationships
*******************************

:py:class:`~.MSImagingRegion` provides geometry helpers for working with
multiple regions:

.. code-block:: python
    :linenos:

    import pyopenms as oms

    r1 = oms.MSImagingRegion.rectangle(1, "A", 0, 0, 3, 3)
    r2 = oms.MSImagingRegion.rectangle(2, "B", 2, 2, 5, 5)
    r3 = oms.MSImagingRegion.rectangle(3, "C", 4, 4, 7, 7)

    print("r1 intersects r2:", r1.intersects(r2))   # True  (overlap at corner)
    print("r1 intersects r3:", r1.intersects(r3))   # False (no overlap)
    print("r1 area:", r1.area())                     # 16 pixels (4 x 4)
