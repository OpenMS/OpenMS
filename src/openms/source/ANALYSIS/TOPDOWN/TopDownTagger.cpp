// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/TopDownTagger.h>

namespace OpenMS
{
  TopDownTagger::TopDownTagger() : DefaultParamHandler("TopDownTagger")
  {
    setDefaultParams_();
  }

  TopDownTagger::TopDownTagger(const TopDownTagger& other) : DefaultParamHandler(other)
  {
  }

  TopDownTagger& TopDownTagger::operator=(const TopDownTagger& rhs)
  {
    if (this == &rhs)
      return *this;

    DefaultParamHandler::operator=(rhs);
    return *this;
  }

  void TopDownTagger::setDefaultParams_()
  {
    defaults_.setValue("min_length", 3, "Minimum length of the tags.");
    defaults_.setMinInt("min_length", 3);

    defaults_.setValue("max_length", 1000, "Maximum length of the tags.");
    defaults_.setMinInt("max_length", 1000);

    defaults_.setValue("tol", DoubleList {10.0, 10.0}, "ppm tolerance for tag generation.");
    defaultsToParam_();
  }

  void TopDownTagger::updateMembers_()
  {
    min_tag_length_ = param_.getValue("min_length");
    max_tag_length_ = param_.getValue("max_length");
    max_tag_length_ = max_tag_length_ < min_tag_length_ ? min_tag_length_ : max_tag_length_;
    ppm_ = param_.getValue("tol");
  }

  void TopDownTagger::run(DeconvolvedSpectrum& dspec, std::vector<std::string>& tags)
  {
    double tol = ppm_[dspec.getOriginalSpectrum().getMSLevel() - 1];
    auto spec = dspec.toSpectrum(1, tol);
    auto tagger = Tagger(min_tag_length_, tol, max_tag_length_ , 1, 1);
    tagger.setUseAbsoluteMzForTol();
    tagger.getTag(spec, tags);
  }
}