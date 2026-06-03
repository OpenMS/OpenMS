// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/IMAGING/MSImagingExperiment.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <cmath>
#include <utility>

namespace OpenMS
{

  namespace
  {
    /// Sums peak intensities of @p spec within the closed m/z window [mz_lo, mz_hi].
    /// Precondition: @p spec is sorted by m/z.
    double sumIntensityInWindow(const MSSpectrum& spec, double mz_lo, double mz_hi)
    {
      double sum = 0.0;
      auto it_hi = spec.MZEnd(mz_hi);
      for (auto it = spec.MZBegin(mz_lo); it != it_hi; ++it)
      {
        sum += static_cast<double>(it->getIntensity());
      }
      return sum;
    }
  } // namespace

  MSImagingExperiment::MSImagingExperiment(MSExperiment exp) : experiment_(std::move(exp))
  {
  }

  MSImagingExperiment& MSImagingExperiment::operator=(MSExperiment exp)
  {
    experiment_ = std::move(exp);
    geometry_.clear();
    return *this;
  }

  MSExperiment& MSImagingExperiment::getMSExperiment()
  {
    return experiment_;
  }

  const MSExperiment& MSImagingExperiment::getMSExperiment() const
  {
    return experiment_;
  }

  void MSImagingExperiment::setMSExperiment(MSExperiment exp)
  {
    experiment_ = std::move(exp);
  }

  MSImagingGeometry& MSImagingExperiment::getGeometry()
  {
    return geometry_;
  }

  const MSImagingGeometry& MSImagingExperiment::getGeometry() const
  {
    return geometry_;
  }

  void MSImagingExperiment::setGeometry(MSImagingGeometry geom)
  {
    geometry_ = std::move(geom);
  }

  Size MSImagingExperiment::getNumberOfPixels() const
  {
    return geometry_.getNumberOfPixels();
  }

  Size MSImagingExperiment::getNumberOfSpectra() const
  {
    return experiment_.getNrSpectra();
  }

  bool MSImagingExperiment::hasPixel(UInt x, UInt y) const
  {
    return geometry_.hasPixel(x, y);
  }

  MSSpectrum& MSImagingExperiment::getSpectrum(UInt x, UInt y)
  {
    const Size idx = geometry_.getSpectrumIndex(x, y);
    if (idx >= experiment_.getNrSpectra())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Pixel references missing spectrum",
                                    String(idx));
    }
    return experiment_[idx];
  }

  const MSSpectrum& MSImagingExperiment::getSpectrum(UInt x, UInt y) const
  {
    const Size idx = geometry_.getSpectrumIndex(x, y);
    if (idx >= experiment_.getNrSpectra())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Pixel references missing spectrum",
                                    String(idx));
    }
    return experiment_[idx];
  }

  IonImage MSImagingExperiment::extractIonImage(double mz, double tolerance_ppm) const
  {
    if (!std::isfinite(mz) || !std::isfinite(tolerance_ppm)
        || mz < 0.0 || tolerance_ppm < 0.0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "mz and tolerance_ppm must be finite and non-negative",
                                    "mz=" + String(mz) + ", tolerance_ppm=" + String(tolerance_ppm));
    }
    const double dm = mz * tolerance_ppm * 1e-6;
    const double mz_lo = mz - dm;
    const double mz_hi = mz + dm;

    IonImage image(geometry_.getWidth(), geometry_.getHeight());
    image.setMzRange(RangeMZ(mz_lo, mz_hi));

    const Size n_spectra = experiment_.getNrSpectra();
    for (const auto& p : geometry_.getPixels())
    {
      if (p.spectrum_index >= n_spectra)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Pixel references missing spectrum",
                                      String(p.spectrum_index));
      }
      const MSSpectrum& spec = experiment_[p.spectrum_index];
      image.setIntensity(p.x, p.y, sumIntensityInWindow(spec, mz_lo, mz_hi));
    }
    return image;
  }

  IonImage MSImagingExperiment::extractIonImage(double mz, double tolerance_ppm, const MSImagingRegion& region) const
  {
    if (!std::isfinite(mz) || !std::isfinite(tolerance_ppm)
        || mz < 0.0 || tolerance_ppm < 0.0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "mz and tolerance_ppm must be finite and non-negative",
                                    "mz=" + String(mz) + ", tolerance_ppm=" + String(tolerance_ppm));
    }
    const double dm = mz * tolerance_ppm * 1e-6;
    const double mz_lo = mz - dm;
    const double mz_hi = mz + dm;

    // Crop to the region's bounding box; local (0,0) == region (min_x, min_y).
    IonImage image(region.getBBoxWidth(), region.getBBoxHeight());
    image.setMzRange(RangeMZ(mz_lo, mz_hi));

    const Size n_spectra = experiment_.getNrSpectra();
    for (const auto& p : geometry_.getPixelsInRegion(region))
    {
      if (p.spectrum_index >= n_spectra)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Pixel references missing spectrum",
                                      String(p.spectrum_index));
      }
      const MSSpectrum& spec = experiment_[p.spectrum_index];
      image.setIntensity(region.toLocalX(p.x), region.toLocalY(p.y),
                         sumIntensityInWindow(spec, mz_lo, mz_hi));
    }
    return image;
  }

  void MSImagingExperiment::validate() const
  {
    const Size n_spectra = experiment_.getNrSpectra();
    for (const auto& p : geometry_.getPixels())
    {
      if (p.spectrum_index >= n_spectra)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Pixel references missing spectrum",
                                      String(p.spectrum_index));
      }
    }
  }

} // namespace OpenMS
