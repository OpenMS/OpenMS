// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/BinnedSpectralContrastAngle.h>

#include <Eigen/Sparse>

using namespace std;

namespace OpenMS
{
  BinnedSpectralContrastAngle::BinnedSpectralContrastAngle() :
    BinnedSpectrumCompareFunctor()
  {
    setName("BinnedSpectralContrastAngle");
    defaultsToParam_();
  }

  BinnedSpectralContrastAngle::BinnedSpectralContrastAngle(const BinnedSpectralContrastAngle& source) :
    BinnedSpectrumCompareFunctor(source)
  {
  }

  BinnedSpectralContrastAngle::~BinnedSpectralContrastAngle() = default;

  BinnedSpectralContrastAngle& BinnedSpectralContrastAngle::operator=(const BinnedSpectralContrastAngle& source)
  {
    if (this != &source)
    {
      BinnedSpectrumCompareFunctor::operator=(source);
    }
    return *this;
  }

  double BinnedSpectralContrastAngle::operator()(const BinnedSpectrum& spec) const
  {
    return operator()(spec, spec);
  }

  void BinnedSpectralContrastAngle::updateMembers_()
  {
  }

  double BinnedSpectralContrastAngle::operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const
  {
    OPENMS_PRECONDITION(BinnedSpectrum::isCompatible(spec1, spec2), "Binned spectra have different bin size or spread");

    // resulting score standardized to interval [0,1]
    const double sum1 = spec1.getBins()->dot(*spec1.getBins());
    const double sum2 = spec2.getBins()->dot(*spec2.getBins());
    const double numerator = spec1.getBins()->dot(*spec2.getBins());
    const double score = numerator / (sqrt(sum1 * sum2));

    return score;
  }
}

