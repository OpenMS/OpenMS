"""
Tests for IMAGING bindings: IonImage, MSImagingGeometry, MSImagingExperiment.
"""

import pytest
import numpy as np


def _make_fixture():
    """2x2 grid; pixel (1, 1) intentionally missing."""
    from pyopenms import MSExperiment, MSSpectrum, MSImagingExperiment, MSImagingGeometry

    exp = MSExperiment()
    s0 = MSSpectrum(); s0.set_peaks(([499.95, 500.05], [10.0, 20.0])); exp.addSpectrum(s0)
    s1 = MSSpectrum(); s1.set_peaks(([490.0],          [5.0]));        exp.addSpectrum(s1)
    s2 = MSSpectrum(); s2.set_peaks(([500.0],          [100.0]));      exp.addSpectrum(s2)

    geom = MSImagingGeometry()
    geom.setDimensions(2, 2)
    geom.addPixel(0, 0, 0)
    geom.addPixel(1, 0, 1)
    geom.addPixel(0, 1, 2)

    mie = MSImagingExperiment(exp)
    mie.setGeometry(geom)
    return mie


def test_ion_image_basic():
    from pyopenms import IonImage

    img = IonImage(3, 2)
    assert img.getWidth() == 3
    assert img.getHeight() == 2
    assert img.hasPixel(0, 0) is False

    img.setIntensity(1, 0, 7.5)
    assert img.hasPixel(1, 0) is True
    assert img.getIntensity(1, 0) == pytest.approx(7.5)


def test_ion_image_get_data_zero_copy():
    from pyopenms import IonImage

    img = IonImage(3, 2)
    img.setIntensity(1, 0, 7.5)

    data = img.get_data()
    assert data.shape == (2, 3)
    assert data.dtype == np.float64
    assert data[0, 1] == pytest.approx(7.5)

    # Zero-copy: writes through the view are reflected in the C++ object
    data[1, 2] = 42.0
    assert img.getIntensity(2, 1) == pytest.approx(42.0)


def test_ion_image_get_mask_copy():
    from pyopenms import IonImage

    img = IonImage(2, 2)
    img.setIntensity(0, 0, 1.0)
    img.setIntensity(1, 1, 2.0)

    mask = img.get_mask()
    assert mask.shape == (2, 2)
    assert mask.dtype == np.uint8
    assert mask[0, 0] == 1
    assert mask[1, 1] == 1
    assert mask[0, 1] == 0
    assert mask[1, 0] == 0


def test_ion_image_oob_raises():
    from pyopenms import IonImage

    img = IonImage(3, 2)
    with pytest.raises(Exception):
        img.getIntensity(5, 0)


def test_ms_imaging_geometry_basic():
    from pyopenms import MSImagingGeometry

    g = MSImagingGeometry()
    g.setDimensions(4, 3)
    g.setPixelSize(25.0, 25.0, "micrometer")
    g.addPixel(0, 0, 100)
    g.addPixel(1, 0, 101)
    assert g.getWidth() == 4
    assert g.getHeight() == 3
    assert g.getNumberOfPixels() == 2
    assert g.hasPixel(0, 0) is True
    assert g.hasPixel(2, 0) is False
    assert g.getSpectrumIndex(1, 0) == 101


def test_ms_imaging_geometry_pixels_struct():
    g = _make_fixture().getGeometry()
    pix = g.pixels_struct()
    assert pix.dtype.names == ("x", "y", "spectrum_index")
    assert pix.dtype.itemsize == 16
    assert pix.shape == (3,)
    assert pix["x"].tolist() == [0, 1, 0]
    assert pix["y"].tolist() == [0, 0, 1]
    assert pix["spectrum_index"].tolist() == [0, 1, 2]


def test_ms_imaging_geometry_rejects_oob():
    from pyopenms import MSImagingGeometry

    g = MSImagingGeometry()
    g.setDimensions(2, 2)
    g.addPixel(1, 1, 0)
    with pytest.raises(Exception):
        g.addPixel(2, 0, 1)
    with pytest.raises(Exception):
        g.addPixel(1, 1, 99)  # duplicate


def test_ms_imaging_experiment_extract_ion_image():
    mie = _make_fixture()
    assert mie.getNumberOfPixels() == 3
    assert mie.getNrSpectra() == 3
    assert mie.hasPixel(1, 1) is False

    # 200 ppm @ 500 Da -> window [499.9, 500.1]
    img = mie.extractIonImage(500.0, 200.0)
    assert img.getWidth() == 2
    assert img.getHeight() == 2

    data = img.get_data()
    assert data[0, 0] == pytest.approx(30.0)   # 10 + 20
    assert data[0, 1] == pytest.approx(0.0)    # spectrum 1 has only 490
    assert data[1, 0] == pytest.approx(100.0)  # spectrum 2 has 500
    assert data[1, 1] == pytest.approx(0.0)    # pixel missing, but data buffer reads 0

    valid = img.get_mask()
    assert valid[1, 1] == 0  # pixel (1, 1) absent -> masked out

    assert img.getMzRange().getMinMZ() == pytest.approx(499.9)
    assert img.getMzRange().getMaxMZ() == pytest.approx(500.1)


def test_ms_imaging_experiment_extract_rejects_bad_inputs():
    mie = _make_fixture()
    with pytest.raises(Exception):
        mie.extractIonImage(-1.0, 100.0)
    with pytest.raises(Exception):
        mie.extractIonImage(500.0, -1.0)


def test_ms_imaging_experiment_get_spectrum():
    mie = _make_fixture()
    sp = mie.getSpectrum(0, 1)
    mz, intens = sp.get_peaks()
    assert list(mz) == pytest.approx([500.0])
    assert list(intens) == pytest.approx([100.0])

    with pytest.raises(Exception):
        mie.getSpectrum(1, 1)


def test_ms_imaging_experiment_assignment_clears_geometry():
    from pyopenms import MSExperiment, MSSpectrum

    mie = _make_fixture()
    fresh = MSExperiment()
    fresh.addSpectrum(MSSpectrum())
    mie.setMSExperiment(fresh)
    # setMSExperiment keeps the geometry
    assert mie.getNumberOfPixels() == 3
    # operator= is exposed in C++ but we can re-create via constructor + reassign:
    mie2 = type(mie)(fresh)
    assert mie2.getNrSpectra() == 1
    assert mie2.getNumberOfPixels() == 0

def test_ms_imaging_region_rectangle():
    from pyopenms import MSImagingRegion

    r = MSImagingRegion.rectangle(7, "region", 1,2,4,7)
    assert r.getId() == 7
    assert r.getName() == "region"
    assert r.getShape() == MSImagingRegion.Shape.Rectangle
    assert r.getBBoxWidth() == 4
    assert r.getBBoxHeight() == 6
    assert r.area() == 24
    assert r.contains(1,2) is True
    assert r.contains(0,2) is False

    mask = r.get_mask()
    assert mask.shape == (6,4)
    assert mask.dtype == np.uint8
    assert mask.min() == 1  #check we have fully true values

    # intersects: r covers x in [1,4], y in [2,7]
    assert r.intersects(MSImagingRegion.rectangle(8, "overlap", 3, 3, 9, 9)) is True
    assert r.intersects(MSImagingRegion.rectangle(9, "disjoint", 10, 10, 12, 12)) is False


def test_ms_imaging_region_from_mask():
    from pyopenms import MSImagingRegion

    r = MSImagingRegion.fromMask(3, "m", 10, 20, 2, 2, [True, False, False, True])
    assert r.getShape() == MSImagingRegion.Shape.Mask
    assert r.area() == 2
    assert r.contains(10, 20) is True
    assert r.contains(11, 20) is False   # dropout bit

    mask = r.get_mask()
    assert mask.shape == (2, 2)
    assert mask.dtype == np.uint8
    assert mask[0, 0] == 1
    assert mask[0, 1] == 0
    assert mask[1, 1] == 1


def test_ms_imaging_region_from_mask_numpy():
    from pyopenms import MSImagingRegion

    # 1D numpy bool array (explicit width/height)
    m1d = np.array([True, False, False, True], dtype=bool)
    r = MSImagingRegion.fromMask(1, "m1d", 10, 20, 2, 2, m1d)
    assert r.getShape() == MSImagingRegion.Shape.Mask
    assert r.area() == 2
    assert r.contains(10, 20) is True    # bit 0
    assert r.contains(11, 20) is False   # bit 1
    assert r.contains(11, 21) is True    # bit 3

    # 2D numpy bool array (height, width); width/height inferred from shape
    m2d = np.array([[True, False],
                    [False, True]], dtype=bool)
    r2 = MSImagingRegion.fromMask(2, "m2d", 10, 20, m2d)   # no width/height args
    assert r2.getShape() == MSImagingRegion.Shape.Mask
    assert r2.getBBoxWidth() == 2
    assert r2.getBBoxHeight() == 2
    assert r2.area() == 2
    assert r2.contains(10, 20) is True    # [row 0, col 0] -> (x=10, y=20)
    assert r2.contains(11, 20) is False   # [0, 1]
    assert r2.contains(11, 21) is True    # [1, 1]

    # non-contiguous (strided) array must be normalized correctly by c_contig
    full = np.array([[True, False, False],
                     [False, True, False]], dtype=bool)   # shape (2, 3)
    strided = full[:, ::2]   # shape (2, 2), columns 0 and 2 -> non-contiguous
    assert not strided.flags["C_CONTIGUOUS"]
    r3 = MSImagingRegion.fromMask(3, "strided", 0, 0, strided)
    # logical strided cells: [[True, False], [False, False]] -> only (0,0) set
    assert r3.area() == 1
    assert r3.contains(0, 0) is True
    assert r3.contains(1, 1) is False


def test_ms_imaging_geometry_regions():
    from pyopenms import MSImagingGeometry, MSImagingRegion

    g = _make_fixture().getGeometry()
    g.addRegion(MSImagingRegion.rectangle(1, "col_0", 0,0,0,1)) # cover column 0 of our data

    assert g.getNumberOfRegions() == 1
    assert g.regionOf(0,0) == 1
    assert g.regionOf(1,0) == MSImagingGeometry.NO_REGION
    assert g.regionOf(10,10) == MSImagingGeometry.NO_REGION

    assert g.getRegion(1).getId() == 1
    spec = g.getRegionSpectrumIndices(1)
    assert list(spec) == [0,2]

    # duplicate id is rejected
    with pytest.raises(Exception):
        g.addRegion(MSImagingRegion.rectangle(1, "duplicate", 5,5,6,6))

    # a distinct id whose footprint overlaps region 1's column 0 is rejected as overlapping
    with pytest.raises(Exception):
        g.addRegion(MSImagingRegion.rectangle(2, "overlapping", 0,0,1,1))

    g.removeRegion(1)
    assert g.getNumberOfRegions() == 0
    assert g.getNumberOfPixels() == 3

def test_ms_imaging_experiment_extract_region():
    from pyopenms import MSImagingRegion

    mie = _make_fixture()
    geom = mie.getGeometry()
    geom.addRegion(MSImagingRegion.rectangle(1, "col0", 0, 0, 0, 1))
    mie.setGeometry(geom)
    img = mie.extractIonImage(500.0, 200.0, 1)   # region overload (3 args)

    assert img.getWidth() == 2 and img.getHeight() == 2   # global dims
    data = img.get_data()
    assert data[0, 0] == pytest.approx(30.0)
    assert data[1, 0] == pytest.approx(100.0)   # (0,1) in region

    mask = img.get_mask()
    assert mask[0, 0] == 1
    assert mask[0, 1] == 0   # (1,0) acquired but OUTSIDE region -> masked

    with pytest.raises(Exception):
        mie.extractIonImage(500.0, 200.0, 99)   # unknown region id

# ---------------------------------------------------------------------------
# Coordinates live on the spectra; the geometry is a rebuildable index.
# ---------------------------------------------------------------------------

def test_ms_imaging_experiment_set_geometry_records_coordinates():
    from pyopenms import MSImagingExperiment

    mie = _make_fixture()
    assert MSImagingExperiment.META_PIXEL_X == "imzml:x"
    # 1-based on the spectrum, 0-based in the geometry
    assert mie[1].getMetaValue("imzml:x") == 2 and mie[1].getMetaValue("imzml:y") == 1
    assert mie[2].getMetaValue("imzml:x") == 1 and mie[2].getMetaValue("imzml:y") == 2
    assert MSImagingExperiment.getPixelCoordinate(mie[2]) == (0, 1)
    mie.validate()


def test_ms_imaging_experiment_index_access():
    mie = _make_fixture()
    assert mie.getNrSpectra() == 3
    assert len(mie) == 3
    assert mie.size() == 3 and mie.empty() is False
    assert mie[2].get_peaks()[0].tolist() == pytest.approx([500.0])
    assert mie.getSpectrum(0).size() == 2
    assert [s.size() for s in mie] == [2, 1, 1]
    with pytest.raises(IndexError):
        mie[3]
    with pytest.raises(IndexError):
        mie.getSpectrum(3)

    # getters copy; the write-back route is __setitem__ / setSpectrum(i, ...)
    s = mie[1]
    s.setRT(7.0)
    assert mie[1].getRT() == -1.0                      # default RT: the copy was edited, not the experiment
    mie[1] = s
    assert mie[1].getRT() == 7.0
    mie.setSpectrum(1, mie[0])
    assert mie[1].size() == 2


def test_ms_imaging_experiment_bind_pixel():
    from pyopenms import MSExperiment, MSSpectrum, MSImagingExperiment

    exp = MSExperiment()
    for _ in range(3):
        exp.addSpectrum(MSSpectrum())
    mie = MSImagingExperiment(exp)
    mie.geometry_view().setDimensions(2, 2)
    mie.bindPixel(1, 1, 2)
    assert mie.hasPixel(1, 1) is True
    assert mie[2].getMetaValue("imzml:x") == 2 and mie[2].getMetaValue("imzml:y") == 2
    assert MSImagingExperiment.getPixelCoordinate(mie[2]) == (1, 1)
    assert MSImagingExperiment.getPixelCoordinate(mie[0]) is None
    with pytest.raises(Exception):
        mie.bindPixel(0, 0, 3)     # no such spectrum
    with pytest.raises(Exception):
        mie.bindPixel(1, 1, 0)     # duplicate pixel
    mie.validate()


def test_ms_imaging_experiment_rebuild_geometry_after_reorder():
    mie = _make_fixture()
    exp = mie.msexperiment_view()
    # reorder the spectra behind the geometry's back: the index is stale ...
    s0, s2 = exp[0], exp[2]
    exp[0], exp[2] = s2, s0
    with pytest.raises(Exception):
        mie.validate()
    # ... and is repaired from the coordinates the spectra carry
    mie.rebuildGeometry()
    mie.validate()
    assert mie.getNumberOfPixels() == 3
    assert mie.geometry_view().getSpectrumIndex(0, 0) == 2
    assert mie.geometry_view().getSpectrumIndex(0, 1) == 0
    assert mie.getSpectrum(0, 0).get_peaks()[0].tolist() == pytest.approx([499.95, 500.05])

    # a dropped spectrum drops its pixel and re-indexes the rest
    exp.setSpectra([exp[0], exp[2]])
    mie.rebuildGeometry()
    mie.validate()
    assert mie.getNumberOfPixels() == 2
    assert mie.hasPixel(1, 0) is False


def test_ms_imaging_experiment_validate_detects_stale_index():
    from pyopenms import MSExperiment

    mie = _make_fixture()
    exp = mie.getMSExperiment()           # a copy
    exp.setSpectra([exp[1], exp[2]])      # spectrum 0 gone: every index shifts by one
    mie.setMSExperiment(exp)              # geometry is kept, now stale
    with pytest.raises(Exception):
        mie.validate()
    mie.rebuildGeometry()
    mie.validate()
    assert mie.getNumberOfPixels() == 2


def test_ms_imaging_experiment_static_coordinate_helpers():
    from pyopenms import MSSpectrum, MSImagingExperiment, MSImagingGeometry, MSExperiment

    s = MSSpectrum()
    assert MSImagingExperiment.getPixelCoordinate(s) is None
    MSImagingExperiment.setPixelCoordinate(s, 3, 0)
    assert s.getMetaValue("imzml:x") == 4 and s.getMetaValue("imzml:y") == 1
    assert s.getMetaValue("imzml:z") == 1
    assert MSImagingExperiment.getPixelCoordinate(s) == (3, 0)

    exp = MSExperiment()
    exp.addSpectrum(s)
    geom = MSImagingGeometry()
    MSImagingExperiment.bindPixelsFromSpectra(exp, geom)
    assert geom.getNumberOfPixels() == 1
    assert geom.getSpectrumIndex(3, 0) == 0
    assert (geom.getWidth(), geom.getHeight()) == (4, 1)


# ---------------------------------------------------------------------------
# Views: msexperiment_view / geometry_view / spectrum_view alias the storage.
# ---------------------------------------------------------------------------

def test_ms_imaging_experiment_msexperiment_view_aliases():
    mie = _make_fixture()
    exp = mie.msexperiment_view()
    assert exp.getNrSpectra() == 3
    exp.spectrum_view(1).setRT(12.5)         # lands without any write-back
    assert mie.getSpectrum(1, 0).getRT() == 12.5
    # the view follows setMSExperiment(): the member is replaced in place
    fresh = mie.getMSExperiment()
    fresh.addSpectrum(fresh[0])
    mie.setMSExperiment(fresh)
    assert exp.getNrSpectra() == 4
    # zero-copy iteration of every spectrum through the view
    for spec in mie.iter_spectrum_views():
        spec.setMSLevel(2)
    assert all(s.getMSLevel() == 2 for s in mie)
    assert len(mie.spectrum_views()) == 4


def test_ms_imaging_experiment_geometry_view_aliases():
    from pyopenms import MSImagingRegion

    mie = _make_fixture()
    g = mie.geometry_view()
    g.addRegion(MSImagingRegion.rectangle(1, "col0", 0, 0, 0, 1))
    assert mie.getRegionSpectrumIndices(1) == [0, 2]   # landed, no setGeometry() needed
    pix = g.pixels_struct()                             # zero-copy end to end
    assert pix["spectrum_index"].tolist() == [0, 1, 2]
    assert mie.getGeometry().getNumberOfRegions() == 1


def test_ms_imaging_experiment_spectrum_view_by_pixel_and_index():
    mie = _make_fixture()
    mie.spectrum_view(0, 1).setRT(3.0)
    assert mie.getSpectrum(0, 1).getRT() == 3.0
    assert mie.spectrum_view(2).getRT() == 3.0          # pixel (0, 1) is spectrum 2
    mie.spectrum_view(2).setRT(4.0)
    assert mie[2].getRT() == 4.0
    with pytest.raises(Exception):
        mie.spectrum_view(1, 1)                          # pixel not in the geometry
    with pytest.raises(IndexError):
        mie.spectrum_view(3)


def test_ms_imaging_experiment_views_keep_parent_alive():
    import gc

    exp = _make_fixture().msexperiment_view()
    spec = _make_fixture().spectrum_view(0, 0)
    gc.collect()
    assert exp.getNrSpectra() == 3                      # parent kept alive by the view
    assert spec.size() == 2


def test_ms_imaging_experiment_setters_keep_pixel_coordinate():
    from pyopenms import MSSpectrum, MSImagingExperiment

    mie = _make_fixture()
    # a spectrum built from scratch inherits the slot's coordinate ...
    mie[1] = MSSpectrum()
    assert MSImagingExperiment.getPixelCoordinate(mie[1]) == (1, 0)
    mie.setSpectrum(2, MSSpectrum())
    assert MSImagingExperiment.getPixelCoordinate(mie[2]) == (0, 1)
    # ... one that carries its own keeps it (the caller may be moving spectra deliberately)
    moved = MSSpectrum()
    MSImagingExperiment.setPixelCoordinate(moved, 1, 1)
    mie[1] = moved
    assert MSImagingExperiment.getPixelCoordinate(mie[1]) == (1, 1)
    with pytest.raises(Exception):
        mie.validate()                                      # index is now stale, as it should report
    # the pixel-addressed setter records the pixel it was given
    mie.setSpectrum(1, 0, MSSpectrum())
    assert MSImagingExperiment.getPixelCoordinate(mie[1]) == (1, 0)
    mie.validate()


def test_ms_imaging_experiment_setters_keep_other_plane_coordinates():
    """A z != 1 spectrum (kept by index, never mapped by the 2D geometry) keeps its own
    acquisition location when swapped into a mapped slot."""
    from pyopenms import MSSpectrum, MSImagingExperiment

    mie = _make_fixture()
    exp = mie.msexperiment_view()
    deep = MSSpectrum()
    deep.setMetaValue("imzml:x", 3)
    deep.setMetaValue("imzml:y", 4)
    deep.setMetaValue("imzml:z", 2)
    exp.addSpectrum(deep)                      # index 3, plane 2: not a geometry pixel
    assert MSImagingExperiment.getPixelCoordinate(mie[3]) is None
    assert mie.getNumberOfPixels() == 3

    # swap the plane-2 spectrum with the plane-1 spectrum of pixel (1, 0)
    s1, s3 = mie[1], mie[3]
    mie[1], mie[3] = s3, s1
    assert (mie[1].getMetaValue("imzml:x"), mie[1].getMetaValue("imzml:y"), mie[1].getMetaValue("imzml:z")) == (3, 4, 2)
    assert (mie[3].getMetaValue("imzml:x"), mie[3].getMetaValue("imzml:y"), mie[3].getMetaValue("imzml:z")) == (2, 1, 1)
    with pytest.raises(Exception):
        mie.validate()                         # pixel (1, 0) now points at a plane-2 spectrum
    mie.rebuildGeometry()                      # ... and the index follows the spectra
    assert mie.geometry_view().getSpectrumIndex(1, 0) == 3
    mie.validate()


def test_ms_imaging_experiment_bind_pixel_rejects_second_binding():
    from pyopenms import MSExperiment, MSSpectrum, MSImagingExperiment

    exp = MSExperiment()
    for _ in range(2):
        exp.addSpectrum(MSSpectrum())
    mie = MSImagingExperiment(exp)
    mie.geometry_view().setDimensions(2, 2)
    mie.bindPixel(0, 0, 0)
    with pytest.raises(Exception):
        mie.bindPixel(1, 0, 0)                 # spectrum 0 already has a pixel
    assert mie.getNumberOfPixels() == 1
    assert MSImagingExperiment.getPixelCoordinate(mie[0]) == (0, 0)
    mie.validate()
    # moving a binding is explicit: unbind, then bind
    mie.unbindPixel(0, 0)
    assert mie.hasPixel(0, 0) is False
    assert MSImagingExperiment.getPixelCoordinate(mie[0]) is None
    mie.bindPixel(1, 0, 0)
    assert mie.geometry_view().getSpectrumIndex(1, 0) == 0
    mie.validate()
    # a geometry binding one spectrum twice is rejected as a whole
    from pyopenms import MSImagingGeometry
    twice = MSImagingGeometry()
    twice.addPixel(0, 0, 1)
    twice.addPixel(0, 1, 1)
    with pytest.raises(Exception):
        mie.setGeometry(twice)
    assert mie.getNumberOfPixels() == 1
