// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/PROCESSING/SCALING/Normalizer.h>

using namespace std;
namespace OpenMS
{
  Normalizer::Normalizer() :
    DefaultParamHandler("Normalizer")
  {
    defaults_.setValue("method", "to_one", "Normalize via dividing by TIC ('to_TIC') per spectrum (i.e. all peaks sum to 1) or normalize to max. intensity to one ('to_one') per spectrum.");
    defaults_.setValidStrings("method", {"to_one","to_TIC"});
    defaultsToParam_();
  }

  Normalizer::~Normalizer() = default;

  Normalizer::Normalizer(const Normalizer & source) :
    DefaultParamHandler(source)
  {
  }

  Normalizer & Normalizer::operator=(const Normalizer & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  void Normalizer::filterPeakSpectrum(PeakSpectrum& spectrum) const
  {
    filterSpectrum(spectrum);
  }

  void Normalizer::filterPeakMap(PeakMap& exp) const
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }

  void Normalizer::updateMembers_()
  {
    method_ = param_.getValue("method").toString();
  }

}
