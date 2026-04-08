// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------
 
#pragma once
 
#include <OpenMS/CONCEPT/Types.h>
#include <functional>
#include <sstream>
 
namespace OpenMS
{

  enum HydrophobicityScaleNumber
    {
      KYTEDOOLITTLE = 0,
      EISENBERG = 1,
      HOPPWOODS = 2,
      BULLBREESE = 3,
      BLACKMOULD = 4,
      GUY = 5
    }; 


  /* class HydrophobicityScale
  {
public:

    enum HydrophobicityScaleNumber
    {
      KYTEDOOLITTLE = 0,
      EISENBERG = 1,
      HOPPWOODS = 2,
      BULLBREESE = 3,
      BLACKMOULD = 4,
      GUY = 5
    }; 

private:
    std::vector<double> HydrophobicityTable_;
  } */
  
} // namespace OpenMS