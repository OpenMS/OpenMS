// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/COMPARISON/SpectraSTSimilarityScore.h>

#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

namespace OpenMS
{
  SpectraSTSimilarityScore::SpectraSTSimilarityScore() :
    PeakSpectrumCompareFunctor()
  {
    setName("SpectraSTSimilarityScore");
  }

  SpectraSTSimilarityScore::SpectraSTSimilarityScore(const SpectraSTSimilarityScore & source) = default;

  SpectraSTSimilarityScore::~SpectraSTSimilarityScore() = default;

  SpectraSTSimilarityScore & SpectraSTSimilarityScore::operator=(const SpectraSTSimilarityScore & source)
  {
    if (this != &source)
    {
      PeakSpectrumCompareFunctor::operator=(source);
    }
    return *this;
  }

  double SpectraSTSimilarityScore::operator()(const PeakSpectrum & spec) const
  {
    return operator()(spec, spec);
  }

  double SpectraSTSimilarityScore::operator()(const PeakSpectrum & s1, const PeakSpectrum & s2) const
  {
    // TODO: check if this operator makes sense (as it doesn't allow to fine tune resolution)
    BinnedSpectrum bin1(s1, 1, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES);
    BinnedSpectrum bin2(s2, 1, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES);

    // normalized dot product
    *bin1.getBins() /= bin1.getBins()->norm();
    *bin2.getBins() /= bin2.getBins()->norm();
    return bin1.getBins()->dot(*bin2.getBins());
  }

  double SpectraSTSimilarityScore::operator()(const BinnedSpectrum & bin1, const BinnedSpectrum & bin2) const
  {
    return bin1.getBins()->dot(*bin2.getBins());
  }

  bool SpectraSTSimilarityScore::preprocess(PeakSpectrum & spec,
                                            float remove_peak_intensity_threshold,
                                            UInt cut_peaks_below,
                                            Size min_peak_number,
                                            Size max_peak_number)
  {
    double min_high_intensity = 0.;
    if (!spec.empty())
    {
      double max_el = std::max_element(spec.begin(),spec.end(),Peak1D::IntensityLess())->getIntensity();
      min_high_intensity = (1.0 / cut_peaks_below) * max_el;
    }

    spec.sortByPosition();

    PeakSpectrum tmp;
    Size s = 0;
    for (PeakSpectrum::iterator k = spec.begin(); k < spec.end() && s < max_peak_number; ++k, ++s)
    {
      Peak1D peak;
      if (k->getIntensity() >  remove_peak_intensity_threshold && k->getIntensity() > min_high_intensity)
      {
        peak.setIntensity(sqrt(k->getIntensity()));
        peak.setMZ(k->getMZ());
        peak.setPosition(k->getPosition());
        tmp.push_back(peak);
      }

    }
    spec = tmp;
    //if not enough peaks in the spectrum pass that one out
    return spec.size() >= min_peak_number;
  }

  BinnedSpectrum SpectraSTSimilarityScore::transform(const PeakSpectrum & spec)
  {
    // TODO: resolution seems rather low. Check with current original implementations.
    BinnedSpectrum bin(spec, 1, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES);
    *bin.getBins() /= bin.getBins()->norm();
    return bin;
  }

  double SpectraSTSimilarityScore::dot_bias(const BinnedSpectrum & bin1, const BinnedSpectrum & bin2, double dot_product) const
  {
    double numerator = (bin1.getBins()->cwiseProduct(*bin2.getBins())).norm();
    
    if (dot_product != 0)
    {
      return (double)numerator / dot_product;
    }
    else
    {
      return (double)numerator / (*this)(bin1, bin2);
    }
  }

  double SpectraSTSimilarityScore::delta_D(double top_hit, double runner_up)
  {
    if (top_hit == 0)
    {
      throw Exception::DivisionByZero(__FILE__, __LINE__, __FUNCTION__);
    }
    else
    {
      return (double)(top_hit - runner_up) / top_hit;
    }
  }

  double SpectraSTSimilarityScore::compute_F(double dot_product, double delta_D, double dot_bias)
  {
    double b(0);
    if (dot_bias < 0.1 || (0.35 < dot_bias && dot_bias <= 0.4))
    {
      b = 0.12;
    }
    else if (0.4 < dot_bias && dot_bias <= 0.45)
    {
      b = 0.18;
    }
    else if (dot_bias > 0.45)
    {
      b  = 0.24;
    }
    return 0.6 * dot_product + 0.4 * delta_D - b;
  }

}
