// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/BinnedSharedPeakCount.h>

using namespace std;

#include <Eigen/Sparse>

namespace OpenMS
{
  BinnedSharedPeakCount::BinnedSharedPeakCount() :
    BinnedSpectrumCompareFunctor()
  {
    setName("BinnedSharedPeakCount");
    defaultsToParam_();
  }

  BinnedSharedPeakCount::BinnedSharedPeakCount(const BinnedSharedPeakCount& source) :
    BinnedSpectrumCompareFunctor(source)
  {
  }

  BinnedSharedPeakCount::~BinnedSharedPeakCount() = default;

  BinnedSharedPeakCount& BinnedSharedPeakCount::operator=(const BinnedSharedPeakCount& source)
  {
    if (this != &source)
    {
      BinnedSpectrumCompareFunctor::operator=(source);
    }
    return *this;
  }

  double BinnedSharedPeakCount::operator()(const BinnedSpectrum& spec) const
  {
    return operator()(spec, spec);
  }

  void BinnedSharedPeakCount::updateMembers_()
  {
  }

  double BinnedSharedPeakCount::operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const
  {
    OPENMS_PRECONDITION(BinnedSpectrum::isCompatible(spec1, spec2), "Binned spectra have different bin size or spread");

    size_t denominator(max(spec1.getBins()->nonZeros(), spec2.getBins()->nonZeros()));

    // Note: keep in single expression for faster computation via expression templates
    // Calculate coefficient-wise product and count non-zero entries
    BinnedSpectrum::SparseVectorType s = spec1.getBins()->cwiseProduct(*spec2.getBins());
    
    // resulting score normalized to interval [0,1]
    return static_cast<double>(s.nonZeros()) / denominator;
  }

}
