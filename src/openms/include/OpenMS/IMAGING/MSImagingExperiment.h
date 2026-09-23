// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Patrick Boschmann $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/IMAGING/IonImage.h>
#include <OpenMS/IMAGING/MSImagingGeometry.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
/**
  @brief In-memory model for a 2D imaging mass spectrometry dataset.

  Owns an MSExperiment (spectra) and an MSImagingGeometry (pixel grid +
  pixel -> spectrum index mapping). Provides pixel-based spectrum access
  and a simple sum-based ion image extraction.

  @section msimagingexperiment_model Where the pixel coordinates live

  Like every other acquisition attribute in OpenMS (RT, drift time, native
  ID), a spectrum's pixel coordinate is stored on the spectrum itself, as
  the meta values @ref META_PIXEL_X / @ref META_PIXEL_Y (1-based, imzML
  convention; @ref META_PIXEL_Z = 1 for the single plane modeled here). The
  geometry is a spatial @em index derived from them: it answers
  (x, y) -> spectrum_index in O(1) and additionally carries the dataset-level
  grid dimensions, pixel size and regions.

  Consequently:
  - bindPixel() and setGeometry() record the coordinates on the bound spectra,
    so the loaders (ImzMLFile, BrukerTimsImagingFile) and hand-assembled
    experiments all leave the spectra self-describing.
  - Reordering, erasing or subsetting spectra in the owned MSExperiment (e.g.
    MSExperiment::sortSpectra(), dropping empty spectra, any tool that
    rewrites the spectrum list) leaves the geometry's indices stale but
    loses no information: call rebuildGeometry() afterwards.
  - validate() detects a stale index by comparing every pixel against the
    coordinate its spectrum carries, not just by checking the index range.

  A 3D / serial-section experiment is modeled as a collection of
  MSImagingExperiment objects, one per section.
*/
class OPENMS_DLLAPI MSImagingExperiment final
{
public:
  /// @name Per-spectrum pixel coordinate meta value keys (1-based, imzML convention)
  //@{
  /// Column (x) coordinate of the pixel a spectrum was acquired at.
  static constexpr const char* META_PIXEL_X = "imzml:x";
  /// Row (y) coordinate of the pixel a spectrum was acquired at.
  static constexpr const char* META_PIXEL_Y = "imzml:y";
  /// Plane (z) coordinate; only z == 1 is represented by the 2D geometry.
  static constexpr const char* META_PIXEL_Z = "imzml:z";
  //@}

  typedef MSExperiment::Iterator Iterator;
  typedef MSExperiment::ConstIterator ConstIterator;

  MSImagingExperiment() = default;

  /**
    @brief Constructs an MSImagingExperiment wrapping @p exp with an empty geometry.
    @param[in] exp MSExperiment to take ownership of (moved in).
  */
  explicit MSImagingExperiment(MSExperiment exp);

  /**
    @brief Replaces the owned MSExperiment and resets the geometry.

    The previous geometry's spectrum indices are tied to the previous
    experiment's spectrum order and become meaningless after assignment,
    so the geometry is cleared. Use setMSExperiment() if you want to
    keep the existing geometry, or rebuildGeometry() to derive a new one
    from the coordinates the spectra carry.

    @param[in] exp MSExperiment to take ownership of (moved in).
    @return Reference to *this.
  */
  MSImagingExperiment& operator=(MSExperiment exp);

  /// @name Access to the owned MSExperiment and geometry
  //@{
  /// @brief Mutable access to the owned MSExperiment.
  ///
  /// Operations that reorder or resize its spectrum list (sortSpectra(),
  /// erase, resize, clear, ...) invalidate the geometry's spectrum indices;
  /// call rebuildGeometry() afterwards (validate() reports the mismatch).
  /// @return Reference to the owned experiment.
  MSExperiment& getMSExperiment();

  /// @brief Read access to the owned MSExperiment.
  /// @return Const reference to the owned experiment.
  const MSExperiment& getMSExperiment() const;

  /**
    @brief Replaces the owned MSExperiment, keeping the geometry.

    The geometry is not touched and the spectra are not re-annotated: if the
    new spectra carry coordinates that disagree with the geometry, validate()
    reports it and rebuildGeometry() re-derives the index from the spectra.

    @param[in] exp New experiment (moved in).
  */
  void setMSExperiment(MSExperiment exp);

  /// @brief Mutable access to the owned geometry (e.g. to add regions).
  ///
  /// Pixels added directly through the geometry are not recorded on their
  /// spectra; prefer bindPixel() for that.
  /// @return Reference to the owned geometry.
  MSImagingGeometry& getGeometry();

  /// @brief Read access to the owned geometry.
  /// @return Const reference to the owned geometry.
  const MSImagingGeometry& getGeometry() const;

  /**
    @brief Replaces the owned geometry and records its pixel coordinates on the bound spectra.

    Every pixel whose spectrum_index lies inside the owned experiment gets
    META_PIXEL_X / META_PIXEL_Y (and META_PIXEL_Z = 1, unless already set)
    written to that spectrum, overwriting whatever coordinate it carried.
    Pixels referencing a missing spectrum are kept in the geometry (validate()
    reports them).

    @param[in] geom New geometry (moved in).
  */
  void setGeometry(MSImagingGeometry geom);
  //@}

  /// @name Spectrum access by index (mirrors OnDiscImzMLExperiment)
  //@{
  /// @brief Number of spectra in the underlying experiment.
  /// @return Spectrum count.
  Size getNrSpectra() const;

  /// @brief Number of spectra in the underlying experiment (alias of getNrSpectra()).
  Size size() const;

  /// @brief true if the underlying experiment holds no spectra.
  bool empty() const;

  /**
    @brief Mutable access to the spectrum at linear @p index.
    @param[in] index Index into the underlying experiment.
    @return Reference to the spectrum.
    @throws Exception::IndexOverflow if @p index >= getNrSpectra().
  */
  MSSpectrum& getSpectrum(Size index);

  /**
    @brief Read access to the spectrum at linear @p index.
    @param[in] index Index into the underlying experiment.
    @return Const reference to the spectrum.
    @throws Exception::IndexOverflow if @p index >= getNrSpectra().
  */
  const MSSpectrum& getSpectrum(Size index) const;

  /// @brief Bounds-checked shorthand for getSpectrum(index).
  MSSpectrum& operator[](Size index);

  /// @brief Bounds-checked shorthand for getSpectrum(index).
  const MSSpectrum& operator[](Size index) const;

  /// @brief Iterator over the spectra of the underlying experiment (acquisition order).
  Iterator begin();
  /// @brief End iterator over the spectra of the underlying experiment.
  Iterator end();
  /// @brief Const iterator over the spectra of the underlying experiment.
  ConstIterator begin() const;
  /// @brief Const end iterator over the spectra of the underlying experiment.
  ConstIterator end() const;
  //@}

  /// @name Pixel-based access (zero-based geometry coordinates)
  //@{
  /// @brief Number of pixels in the geometry.
  /// @return Geometry pixel count.
  Size getNumberOfPixels() const;

  /**
    @brief Tests pixel presence at (@p x, @p y).
    @param[in] x Column index.
    @param[in] y Row index.
    @return true if a pixel exists at that coordinate.
  */
  bool hasPixel(UInt x, UInt y) const;

  /**
    @brief Binds the pixel at (@p x, @p y) to the spectrum at @p spectrum_index.

    Adds the pixel to the geometry and records the coordinate on the spectrum
    (META_PIXEL_X / META_PIXEL_Y = x + 1 / y + 1, META_PIXEL_Z = 1 unless set),
    so the spectrum stays self-describing and rebuildGeometry() can restore
    the binding later.

    @param[in] x              Column index (zero-based).
    @param[in] y              Row index (zero-based).
    @param[in] spectrum_index Index into the underlying experiment.
    @throws Exception::IndexOverflow if @p spectrum_index >= getNrSpectra().
    @throws Exception::InvalidValue on duplicate coordinates, or if the
            geometry has dimensions and (@p x, @p y) lies outside them.
  */
  void bindPixel(UInt x, UInt y, Size spectrum_index);

  /**
    @brief Mutable access to the spectrum bound to the pixel at (@p x, @p y).
    @param[in] x Column index.
    @param[in] y Row index.
    @return Reference to the bound spectrum.
    @throws Exception::ElementNotFound if no pixel exists at that coordinate.
    @throws Exception::InvalidValue if the pixel references a spectrum_index
            outside the underlying experiment.
  */
  MSSpectrum& getSpectrum(UInt x, UInt y);

  /**
    @brief Read access to the spectrum bound to the pixel at (@p x, @p y).
    @param[in] x Column index.
    @param[in] y Row index.
    @return Const reference to the bound spectrum.
    @throws Exception::ElementNotFound if no pixel exists at that coordinate.
    @throws Exception::InvalidValue if the pixel references a spectrum_index
            outside the underlying experiment.
  */
  const MSSpectrum& getSpectrum(UInt x, UInt y) const;
  //@}

  /// @name Consistency between spectra and geometry
  //@{
  /**
    @brief Re-derives the pixel -> spectrum_index mapping from the coordinates the spectra carry.

    Drops all pixels from the geometry (dimensions, pixel size and regions are
    kept) and re-adds one per spectrum that carries META_PIXEL_X / META_PIXEL_Y,
    bound to that spectrum's current index. Call this after reordering,
    erasing or subsetting spectra in the owned MSExperiment.

    Follows the loaders' tolerance rules: spectra without coordinates, with a
    coordinate < 1, or on a plane other than z == 1 are skipped silently or
    with a warning; a spectrum whose pixel lies outside the geometry's declared
    dimensions is excluded with a warning; when several spectra share a pixel,
    only the first one is mapped (a summary warning lists the others). Every
    spectrum stays reachable by index. If the geometry has no dimensions yet,
    they are set to the bounding box of the mapped pixels.
  */
  void rebuildGeometry();

  /**
    @brief Re-derives a geometry's pixels from the coordinates carried by @p exp's spectra.

    Static building block of rebuildGeometry() and ImzMLFile::buildImagingGeometry():
    appends one pixel per coordinate-carrying spectrum to @p geom (which is expected
    to hold no pixels yet; dimensions and pixel size are honoured if set, and the
    dimensions are grown to the mapped bounding box if unset). Tolerance rules as
    described for rebuildGeometry().

    @param[in]     exp  Experiment whose spectra carry META_PIXEL_X / META_PIXEL_Y.
    @param[in,out] geom Geometry receiving the pixels.
  */
  static void bindPixelsFromSpectra(const MSExperiment& exp, MSImagingGeometry& geom);

  /**
    @brief Validates that the geometry is a consistent index of the spectra.

    Checks, for every pixel, that its spectrum_index lies inside the experiment
    and — when that spectrum carries META_PIXEL_X / META_PIXEL_Y — that the
    coordinate it carries is the pixel's own (and META_PIXEL_Z, if present, is 1).
    A mismatch means the spectrum list was reordered or replaced after the
    geometry was built; rebuildGeometry() repairs it. Spectra without
    coordinates are accepted as-is.

    @throws Exception::InvalidValue on the first dangling or mismatching pixel.
  */
  void validate() const;

  /**
    @brief Reads the pixel coordinate a spectrum carries.
    @param[in]  spec Spectrum to inspect.
    @param[out] x    Zero-based column index (META_PIXEL_X - 1).
    @param[out] y    Zero-based row index (META_PIXEL_Y - 1).
    @return true if the spectrum carries META_PIXEL_X and META_PIXEL_Y, both >= 1,
            and lies on plane z == 1 (META_PIXEL_Z absent or 1); false otherwise.
  */
  static bool getPixelCoordinate(const MSSpectrum& spec, UInt& x, UInt& y);

  /**
    @brief Records a pixel coordinate on a spectrum (META_PIXEL_X/Y = x + 1 / y + 1; META_PIXEL_Z = 1 unless set).
    @param[in,out] spec Spectrum to annotate.
    @param[in]     x    Zero-based column index.
    @param[in]     y    Zero-based row index.
  */
  static void setPixelCoordinate(MSSpectrum& spec, UInt x, UInt y);
  //@}

  /**
    @brief Extracts an ion image by summing peak intensities inside
           [mz - dm, mz + dm], with dm = mz * tolerance_ppm * 1e-6.

    Preconditions: each referenced spectrum must be sorted by m/z.
    The ImzMLFile loader guarantees this; callers assembling the
    experiment by hand must ensure it themselves.

    Pixels not present in the geometry remain invalid in the returned
    image. Pixels with a spectrum but no peaks in the window are marked
    valid with intensity 0. The returned image's m/z range is set to
    [mz - dm, mz + dm].

    @param[in] mz             m/z center of the extraction window (>= 0).
    @param[in] tolerance_ppm  Half-window width in ppm (>= 0).
    @return Image of the same dimensions as the geometry.
    @throws Exception::InvalidValue if @p mz or @p tolerance_ppm is
            negative or non-finite, or if any pixel references a missing
            spectrum_index.
    @throws Exception::IndexOverflow if a pixel coordinate falls outside
            the geometry's declared dimensions.
  */
  IonImage extractIonImage(double mz, double tolerance_ppm) const;

  /**
    @brief Extracts an ion image by summing peak intensities inside
           [mz - dm, mz + dm], with dm = mz * tolerance_ppm * 1e-6.
           Overloaded with region id to extract a specific region only.

    Preconditions: each referenced spectrum must be sorted by m/z.
    The ImzMLFile loader guarantees this; callers assembling the
    experiment by hand must ensure it themselves.

    Pixels not present in the geometry remain invalid in the returned
    image. Pixels with a spectrum but no peaks in the window are marked
    valid with intensity 0. The returned image's m/z range is set to
    [mz - dm, mz + dm].

    @param[in] mz             m/z center of the extraction window (>= 0).
    @param[in] tolerance_ppm  Half-window width in ppm (>= 0).
    @param[in] region_id region of pixels to be extracted
    @return Image of the same dimensions as the geometry.
    @throws Exception::InvalidValue if @p mz or @p tolerance_ppm is
            negative or non-finite, or if any pixel references a missing
            spectrum_index.
    @throws Exception::IndexOverflow if a pixel coordinate falls outside
            the geometry's declared dimensions.
    @throws Exception::ElementNotFound if an unknown region id is supplied.
  */
  IonImage extractIonImage(double mz, double tolerance_ppm, Size region_id) const;

  /**
    @brief Spectrum indices of the acquired pixels belonging to a region.

    Delegates to the geometry; the returned values are spectrum_index values
    (indices into the bound MSExperiment), not pixel positions.
    @param[in] region_id Region id.
    @return The spectrum_index values of the region's member pixels.
    @throws Exception::ElementNotFound if @p region_id is unknown.
  */
  std::vector<Size> getRegionSpectrumIndices(Size region_id) const;


private:
  /// Writes every in-range pixel's coordinate onto its spectrum (see setGeometry()).
  void stampPixelCoordinates_();

  MSExperiment experiment_;
  MSImagingGeometry geometry_;
};

} // namespace OpenMS
