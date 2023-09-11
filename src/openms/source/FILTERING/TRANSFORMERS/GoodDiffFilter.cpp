// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>

using namespace std;

namespace OpenMS
{
  GoodDiffFilter::GoodDiffFilter() :
    FilterFunctor()
  {
    setName(GoodDiffFilter::getProductName());
    //values from kinter sherman


    // TODO from CHEMISTRY!
    aamass_.insert(make_pair(57.02, 'G'));
    aamass_.insert(make_pair(71.04, 'A'));
    aamass_.insert(make_pair(87.03, 'S'));
    aamass_.insert(make_pair(97.05, 'P'));
    aamass_.insert(make_pair(99.07, 'V'));
    aamass_.insert(make_pair(101.05, 'T'));
    aamass_.insert(make_pair(103.01, 'C'));
    aamass_.insert(make_pair(113.08, 'L')); //and also I, but the chars are fillers, anyway
    aamass_.insert(make_pair(114.04, 'N'));
    aamass_.insert(make_pair(115.03, 'D'));
    aamass_.insert(make_pair(128.06, 'Q'));
    aamass_.insert(make_pair(128.09, 'K'));
    aamass_.insert(make_pair(129.04, 'E'));
    aamass_.insert(make_pair(131.04, 'M'));
    aamass_.insert(make_pair(137.06, 'H'));
    aamass_.insert(make_pair(147.07, 'F'));
    aamass_.insert(make_pair(156.10, 'R'));
    aamass_.insert(make_pair(163.06, 'Y'));
    aamass_.insert(make_pair(186.06, 'W'));

    //value from Bioinformatics, Bern 2004
    defaults_.setValue("tolerance", 0.37, "Tolerance value as defined by Bern et al.");
    defaultsToParam_();
  }

  GoodDiffFilter::GoodDiffFilter(const GoodDiffFilter & source)  = default;

  GoodDiffFilter & GoodDiffFilter::operator=(const GoodDiffFilter & source)
  {
    if (this != &source)
    {
      FilterFunctor::operator=(source);
      aamass_ = source.aamass_;
    }
    return *this;
  }

  GoodDiffFilter::~GoodDiffFilter() = default;

}
