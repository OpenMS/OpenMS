// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/PeakSpectrumCompareFunctor.h>

#include <OpenMS/COMPARISON/SpectrumCheapDPCorr.h>
#include <OpenMS/COMPARISON/SpectrumPrecursorComparator.h>
#include <OpenMS/COMPARISON/ZhangSimilarityScore.h>
#include <OpenMS/COMPARISON/SpectrumAlignmentScore.h>
#include <OpenMS/COMPARISON/SteinScottImproveScore.h>
#include <OpenMS/COMPARISON/PeakAlignment.h>

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
