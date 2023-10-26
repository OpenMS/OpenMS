// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeon $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/Tagger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <iomanip>
#include <iostream>

namespace OpenMS
{
  /**
  @brief
  @ingroup Topdown
  */

  class OPENMS_DLLAPI TopDownTagger : public DefaultParamHandler
  {
  public:

    /// constructor
    TopDownTagger();

    /// destructor
    ~TopDownTagger() override = default;

    /// copy constructor
    TopDownTagger(const TopDownTagger&);

    /// move constructor
    TopDownTagger(TopDownTagger&& other) = default;

    /// assignment operator
    TopDownTagger& operator=(const TopDownTagger& other);

    void run(DeconvolvedSpectrum& dspec, std::vector<std::string>& tags);

  protected:
    void updateMembers_() override;
    /// implemented for DefaultParamHandler
    void setDefaultParams_();

  private:
    uint min_tag_length_;
    uint max_tag_length_;
    DoubleList ppm_;
  };
} // namespace OpenMS