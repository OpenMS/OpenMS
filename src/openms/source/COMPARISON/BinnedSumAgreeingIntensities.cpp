// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/BinnedSumAgreeingIntensities.h>

#include <Eigen/Sparse>

using namespace std;

namespace OpenMS
{
  BinnedSumAgreeingIntensities::BinnedSumAgreeingIntensities() :
    BinnedSpectrumCompareFunctor()
  {
    setName("BinnedSumAgreeingIntensities");
    defaultsToParam_();
  }

  BinnedSumAgreeingIntensities::BinnedSumAgreeingIntensities(const BinnedSumAgreeingIntensities& source) :
    BinnedSpectrumCompareFunctor(source)
  {
  }

  BinnedSumAgreeingIntensities::~BinnedSumAgreeingIntensities() = default;

  BinnedSumAgreeingIntensities& BinnedSumAgreeingIntensities::operator=(const BinnedSumAgreeingIntensities& source)
  {
    if (this != &source)
    {
      BinnedSpectrumCompareFunctor::operator=(source);
    }
    return *this;
  }

  double BinnedSumAgreeingIntensities::operator()(const BinnedSpectrum& spec) const
  {
    return operator()(spec, spec);
  }

  void BinnedSumAgreeingIntensities::updateMembers_()
  {
  }

  double BinnedSumAgreeingIntensities::operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const
  {
    OPENMS_PRECONDITION(BinnedSpectrum::isCompatible(spec1, spec2), "Binned spectra have different bin size or spread");

    const double sum1 = spec1.getBins()->sum();
    const double sum2 = spec2.getBins()->sum();

    // For maximum speed, keep in single expression
    // 1. calculate mean minus difference: x = mean(a,b) - abs(a-b)
    // 2. truncate negative values:        y = max(0, x)
    // 3. calculate sum of entries:   sum_nn = y.sum()
    BinnedSpectrum::SparseVectorType s = ((*spec1.getBins() + *spec2.getBins()) * 0.5) - ((*spec1.getBins() - *spec2.getBins()).cwiseAbs());
    double sum_nn = s.coeffs().cwiseMax(0).sum();

    // resulting score normalized to interval [0,1]
    return min(sum_nn / ((sum1 + sum2) / 2.0), 1.0);
  }
}

