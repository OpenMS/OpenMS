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
    @brief In-memory model for an imaging mass spectrometry dataset.

    Owns an MSExperiment (spectra) and an MSImagingGeometry (pixel grid +
    pixel -> spectrum index mapping). Provides pixel-based spectrum access
    and a simple sum-based ion image extraction.
  */
  class OPENMS_DLLAPI MSImagingExperiment final
  {
  public:
    MSImagingExperiment() = default;

    /// @name MSExperiment accessors (owned by value)
    ///@{
    MSExperiment& getMSExperiment();
    const MSExperiment& getMSExperiment() const;
    void setMSExperiment(MSExperiment exp);
    ///@}

    /// @name MSImagingGeometry accessors (owned by value)
    ///@{
    MSImagingGeometry& getGeometry();
    const MSImagingGeometry& getGeometry() const;
    void setGeometry(MSImagingGeometry geom);
    ///@}

    /// @brief Number of pixels in the geometry.
    Size getNumberOfPixels() const;

    /// @brief Number of spectra in the underlying experiment.
    Size getNumberOfSpectra() const;

    /// @brief True if a pixel exists at (x, y, z).
    bool hasPixel(UInt x, UInt y, UInt z = 0) const;

    /// @brief Spectrum bound to the pixel at (x, y, z).
    /// @throws Exception::ElementNotFound if no pixel exists at that coordinate.
    MSSpectrum& getSpectrum(UInt x, UInt y, UInt z = 0);

    /// @brief Spectrum bound to the pixel at (x, y, z).
    /// @throws Exception::ElementNotFound if no pixel exists at that coordinate.
    const MSSpectrum& getSpectrum(UInt x, UInt y, UInt z = 0) const;

    /**
      @brief Extract an ion image by summing peak intensities inside
             [mz - dm, mz + dm], with dm = mz * tolerance_ppm * 1e-6.

      Preconditions: each referenced spectrum must be sorted by m/z.
      Phase 2's ImzMLFile loader will guarantee this; Phase 1 callers
      must ensure it manually.

      Pixels not present in the geometry remain invalid in the returned
      image. Pixels with a spectrum but no peaks in the window are marked
      valid with intensity 0.

      The returned image's m/z range is set to [mz - dm, mz + dm].

      @throws Exception::InvalidValue if any pixel references a missing
              spectrum_index.
      @throws Exception::IndexOverflow if a pixel coordinate falls outside
              the geometry's declared dimensions.
    */
    IonImage extractIonImage(double mz, double tolerance_ppm) const;

    /// @brief Throws if any pixel references a spectrum_index that is out
    ///        of bounds for the underlying experiment.
    /// @throws Exception::InvalidValue
    void validate() const;

  private:
    MSExperiment experiment_;
    MSImagingGeometry geometry_;
  };

} // namespace OpenMS
