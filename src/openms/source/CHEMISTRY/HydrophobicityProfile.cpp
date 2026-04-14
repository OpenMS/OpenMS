// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------
 
/// @todo make the calculation of average hydrophobicity reusable
/// hydrophobic moment explainnnnn

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

  std::vector<double> HydrophobicityProfile::computeProfile
  (
    const AASequence& seq, 
    const HydrophobicityScaleMethod scale
  )
  {
    if (seq.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Can not compute profile of an empty sequence", "");
    }
    std::vector<double> profile;
    for (const auto& residue : seq) 
    {  
      profile.push_back(residue.getHydrophobicity(scale));
    }
    return profile;
  }

  std::vector<double> HydrophobicityProfile::computeWindowedProfile
  (
    const AASequence& seq,
    Size window_size,
    const HydrophobicityScaleMethod scale
  )
  {
    if (window_size>seq.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,window_size, "Window size can not be larger than the sequence");
    }
    if (seq.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Can not compute profile of an empty sequence", "");
    }
    std::vector<double> profile;
    for (auto i = 0; i <= seq.size()-window_size; i++)
    {
      AASequence subsequence = seq.getSubsequence(i,window_size);
      double sum = 0;
      for (const auto& residue : subsequence)
      {     
        sum += residue.getHydrophobicity(scale);
      }
      sum = sum / window_size;
      profile.push_back(sum);
    }
    return profile;
  }

  std::vector<double> HydrophobicityProfile::computeHydrophobicMoment
  (
    const AASequence& seq,
    Size window_size,
    double angle  
  )
  {
    if (window_size>seq.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,window_size, "Window size can not be larger than the sequence");
    }
    if (seq.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Can not compute profile of an empty sequence", "");
    }
    std::vector<double> profile;
    for (auto i = 0; i <= seq.size()-window_size; i++)
    {
      AASequence subsequence = seq.getSubsequence(i,window_size);
      double sum_1 = 0;
      double sum_2 = 0;
      for (const auto& residue : subsequence)
      {
        sum_1 += residue.getHydrophobicity(HydrophobicityScaleMethod::EISENBERG)*std::sin((100*window_size*3.14159)/180);
        sum_1 += residue.getHydrophobicity(HydrophobicityScaleMethod::EISENBERG)*std::cos((100*window_size*3.14159)/180);
      }
      sum_1 = std::pow(sum_1,2);
      profile.push_back(std::pow(sum_1+sum_2,0.5));
    }
  }


  
      
  int HydrophobicityProfile::testfunktion()
  {
    Residue res;
    res.setOneLetterCode("A");
    return 1;
  } 
} // namespace OpenMS