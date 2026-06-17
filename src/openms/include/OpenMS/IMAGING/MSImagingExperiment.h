// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
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

    A 3D / serial-section experiment is modeled as a collection of
    MSImagingExperiment objects, one per section.
  */
  class OPENMS_DLLAPI MSImagingExperiment final
  {
  public:
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
      keep the existing geometry.

      @param[in] exp MSExperiment to take ownership of (moved in).
      @return Reference to *this.
    */
    MSImagingExperiment& operator=(MSExperiment exp);

    /// @brief Mutable access to the owned MSExperiment.
    /// @return Reference to the owned experiment.
    MSExperiment& getMSExperiment();

    /// @brief Read access to the owned MSExperiment.
    /// @return Const reference to the owned experiment.
    const MSExperiment& getMSExperiment() const;

    /**
      @brief Replaces the owned MSExperiment.
      @param[in] exp New experiment (moved in).
    */
    void setMSExperiment(MSExperiment exp);

    /// @brief Mutable access to the owned geometry.
    /// @return Reference to the owned geometry.
    MSImagingGeometry& getGeometry();

    /// @brief Read access to the owned geometry.
    /// @return Const reference to the owned geometry.
    const MSImagingGeometry& getGeometry() const;

    /**
      @brief Replaces the owned geometry.
      @param[in] geom New geometry (moved in).
    */
    void setGeometry(MSImagingGeometry geom);

    /// @brief Number of pixels in the geometry.
    /// @return Geometry pixel count.
    Size getNumberOfPixels() const;

    /// @brief Number of spectra in the underlying experiment.
    /// @return Spectrum count.
    Size getNumberOfSpectra() const;

    /**
      @brief Tests pixel presence at (@p x, @p y).
      @param[in] x Column index.
      @param[in] y Row index.
      @return true if a pixel exists at that coordinate.
    */
    bool hasPixel(UInt x, UInt y) const;

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

    /**
      @brief Extracts an ion image by summing peak intensities inside
             [mz - dm, mz + dm], with dm = mz * tolerance_ppm * 1e-6.

      Preconditions: each referenced spectrum must be sorted by m/z.
      Phase 2's ImzMLFile loader will guarantee this; Phase 1 callers
      must ensure it manually.

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
      @brief Validates that every pixel references an in-range spectrum_index.
      @throws Exception::InvalidValue on the first dangling reference.
    */
    void validate() const;

  private:
    MSExperiment experiment_;
    MSImagingGeometry geometry_;
  };

} // namespace OpenMS
