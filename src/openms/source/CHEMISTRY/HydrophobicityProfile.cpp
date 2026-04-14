// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------
 
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/CHEMISTRY/HydrophobicityProfile.h>

namespace OpenMS
{
   double HydrophobicityProfile::computeGRAVY(const AASequence& seq)
   {
      if(seq.empty())
      {
         return 0;
      }      
      double sum = 0;
      for (const auto& residue : seq) //std accumulate
      {     
         sum += residue.getHydrophobicity(HydrophobicityScaleMethod::KYTE_DOOLITTLE);
      }
      return sum/seq.size(); 
   } 
      
   int HydrophobicityProfile::testfunktion()
   {
      Residue res;
      res.setOneLetterCode("A");
      return 1;
   } 
} // namespace OpenMS