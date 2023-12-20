// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
    DefaultParamHandler(PeakSpectrumCompareFunctor::getProductName())
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

  void PeakSpectrumCompareFunctor::registerChildren()
  {
    Factory<PeakSpectrumCompareFunctor>::registerProduct(SpectrumCheapDPCorr::getProductName(), &SpectrumCheapDPCorr::create);
    Factory<PeakSpectrumCompareFunctor>::registerProduct(SpectrumPrecursorComparator::getProductName(), &SpectrumPrecursorComparator::create);
    Factory<PeakSpectrumCompareFunctor>::registerProduct(ZhangSimilarityScore::getProductName(), &ZhangSimilarityScore::create);
    Factory<PeakSpectrumCompareFunctor>::registerProduct(SpectrumAlignmentScore::getProductName(), &SpectrumAlignmentScore::create);
    Factory<PeakSpectrumCompareFunctor>::registerProduct(SteinScottImproveScore::getProductName(), &SteinScottImproveScore::create);
    Factory<PeakSpectrumCompareFunctor>::registerProduct(PeakAlignment::getProductName(), &PeakAlignment::create);
  }

}
