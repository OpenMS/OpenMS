// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/IMAGING/MSImagingExperiment.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <utility>

namespace OpenMS
{

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

  bool MSImagingExperiment::hasPixel(UInt x, UInt y, UInt z) const
  {
    return geometry_.hasPixel(x, y, z);
  }

  MSSpectrum& MSImagingExperiment::getSpectrum(UInt x, UInt y, UInt z)
  {
    return experiment_[geometry_.getSpectrumIndex(x, y, z)];
  }

  const MSSpectrum& MSImagingExperiment::getSpectrum(UInt x, UInt y, UInt z) const
  {
    return experiment_[geometry_.getSpectrumIndex(x, y, z)];
  }

  IonImage MSImagingExperiment::extractIonImage(double mz, double tolerance_ppm) const
  {
    const double dm = mz * tolerance_ppm * 1e-6;
    const double mz_lo = mz - dm;
    const double mz_hi = mz + dm;

    IonImage image(geometry_.getWidth(), geometry_.getHeight(), geometry_.getDepth());
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
      double sum = 0.0;
      auto it_lo = spec.MZBegin(mz_lo);
      auto it_hi = spec.MZEnd(mz_hi);
      for (auto it = it_lo; it != it_hi; ++it)
      {
        sum += static_cast<double>(it->getIntensity());
      }
      image.setIntensity(p.x, p.y, p.z, sum);
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
