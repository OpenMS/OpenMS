// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>

#include <OpenMS/COMPARISON/SPECTRA/SpectrumCheapDPCorr.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumPrecursorComparator.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SteinScottImproveScore.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakAlignment.h>
#include <OpenMS/CONCEPT/Factory.h>

using namespace std;

namespace OpenMS
{
  PeakSpectrumCompareFunctor::PeakSpectrumCompareFunctor() :
    DefaultParamHandler("PeakSpectrumCompareFunctor")
  {
  }

  PeakSpectrumCompareFunctor::PeakSpectrumCompareFunctor(const PeakSpectrumCompareFunctor & source) = default;

  PeakSpectrumCompareFunctor::~PeakSpectrumCompareFunctor() = default;

  PeakSpectrumCompareFunctor & PeakSpectrumCompareFunctor::operator=(const PeakSpectrumCompareFunctor & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

}
