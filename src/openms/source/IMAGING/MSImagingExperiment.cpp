// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Patrick Boschmann $
// --------------------------------------------------------------------------

#include <OpenMS/IMAGING/MSImagingExperiment.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/IMAGING/IonImageExtraction.h>
#include <cmath>
#include <numeric>
#include <utility>

namespace OpenMS
{

MSImagingExperiment::MSImagingExperiment(MSExperiment exp): experiment_(std::move(exp))
{
}

MSImagingExperiment& MSImagingExperiment::operator=(MSExperiment exp)
{
  experiment_ = std::move(exp);
  geometry_.clear();
  return *this;
}

MSExperiment& MSImagingExperiment::getMSExperiment()
{ return experiment_; }

const MSExperiment& MSImagingExperiment::getMSExperiment() const
{ return experiment_; }

void MSImagingExperiment::setMSExperiment(MSExperiment exp)
{ experiment_ = std::move(exp); }

MSImagingGeometry& MSImagingExperiment::getGeometry()
{ return geometry_; }

const MSImagingGeometry& MSImagingExperiment::getGeometry() const
{ return geometry_; }

void MSImagingExperiment::setGeometry(MSImagingGeometry geom)
{ geometry_ = std::move(geom); }

Size MSImagingExperiment::getNumberOfPixels() const
{ return geometry_.getNumberOfPixels(); }

Size MSImagingExperiment::getNumberOfSpectra() const
{ return experiment_.getNrSpectra(); }

bool MSImagingExperiment::hasPixel(UInt x, UInt y) const
{ return geometry_.hasPixel(x, y); }

  MSSpectrum& MSImagingExperiment::getSpectrum(UInt x, UInt y)
  {
    const Size idx = geometry_.getSpectrumIndex(x, y);
    if (idx >= experiment_.getNrSpectra())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Pixel references missing spectrum",
                                    StringUtils::toStr(idx));
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
                                    StringUtils::toStr(idx));
    }
    return experiment_[idx];
  }

  IonImage MSImagingExperiment::extractIonImage(double mz, double tolerance_ppm) const
  {
    std::vector<Size> all(geometry_.getNumberOfPixels());
    std::iota(all.begin(), all.end(), Size(0));
    return Internal::extractIonImage(geometry_, mz, tolerance_ppm, all, experiment_.getNrSpectra(),
                                     [this](Size i) -> const MSSpectrum& { return experiment_[i]; });
  }

  IonImage MSImagingExperiment::extractIonImage(double mz, double tolerance_ppm, Size region_id) const
  {
    const auto& indices = geometry_.getRegionPixels(region_id);
    return Internal::extractIonImage(geometry_, mz, tolerance_ppm, indices, experiment_.getNrSpectra(),
                                     [this](Size i) -> const MSSpectrum& { return experiment_[i]; });
  }

  std::vector<Size> MSImagingExperiment::getRegionSpectrumIndices(Size region_id) const
  {
    return geometry_.getRegionSpectrumIndices(region_id);
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
                                      StringUtils::toStr(p.spectrum_index));
      }
    }
  }

} // namespace OpenMS
