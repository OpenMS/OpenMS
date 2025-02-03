// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/SpectrumPrecursorComparator.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace std;

namespace OpenMS
{
  SpectrumPrecursorComparator::SpectrumPrecursorComparator() :
    PeakSpectrumCompareFunctor()
  {
    setName("SpectrumPrecursorComparator");
    defaults_.setValue("window", 2, "Allowed deviation between precursor peaks.");
    defaultsToParam_();
  }

  SpectrumPrecursorComparator::SpectrumPrecursorComparator(const SpectrumPrecursorComparator & source) = default;

  SpectrumPrecursorComparator::~SpectrumPrecursorComparator() = default;

  SpectrumPrecursorComparator & SpectrumPrecursorComparator::operator=(const SpectrumPrecursorComparator & source)
  {
    if (this != &source)
    {
      PeakSpectrumCompareFunctor::operator=(source);
    }
    return *this;
  }

  double SpectrumPrecursorComparator::operator()(const PeakSpectrum & spec) const
  {
    return operator()(spec, spec);
  }

  double SpectrumPrecursorComparator::operator()(const PeakSpectrum & x, const PeakSpectrum & y) const
  {
    double window = (double)param_.getValue("window");

    double mz1 = 0.0;
    if (!x.getPrecursors().empty())
    {
      mz1 = x.getPrecursors()[0].getMZ();
    }
    double mz2 = 0.0;
    if (!y.getPrecursors().empty())
    {
      mz2 = y.getPrecursors()[0].getMZ();
    }

    if (fabs(mz1 - mz2) > window)
    {
      return 0;
    }

    return window - fabs(mz1 - mz2);
  }

}
