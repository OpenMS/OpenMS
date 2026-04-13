// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------
 
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/CHEMISTRY/HydrophobicityProfile.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

namespace OpenMS
{
   double HydrophobicityProfile::computeGRAVY(const AASequence& seq)
   {
      if(seq.size() == 0)
      {
         return 0;
      }      
      double sum = 0;
      for (auto residue : seq) 
      {     
         sum += residue.getHydrophobicity(KYTEDOOLITTLE);
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