// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------
 
#include <OpenMS/CHEMISTRY/HydrophobicityProfile.h>
#include <OpenMS/KERNEL/DPeak.h>


namespace OpenMS
{
  double HydrophobicityProfile::computeGRAVY(const AASequence& seq)
  {
    if (seq.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Can not compute GRAVY score of an empty sequence", "");
    }  
    double sum = 0;
    for (const auto& residue : seq)
    {     
      sum += residue.getHydrophobicity(HydrophobicityScaleMethod::KYTE_DOOLITTLE);
    }
    return sum / seq.size(); 
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
    if (seq.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Can not compute profile of an empty sequence", "");
    }
    if (window_size<=0)
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,window_size,"window size can not be 0");
    }
    if (window_size > seq.size())
    {
      OPENMS_LOG_WARN << "Warning: window size is larger than sequence size. Using seequence size as window size: " << seq.size() << "\n";
    }
    std::vector<double> profile;
    AASequence subsequence = seq.getSubsequence(0, std::min(window_size, seq.size())); // first window
    double sum = 0;
    for (const auto& residue : subsequence) // calculate score for first window
    {
      sum += residue.getHydrophobicity(scale);
    }
    profile.push_back(sum / std::min(window_size, seq.size()));
    double prev_value = seq[0].getHydrophobicity(scale); // value for amino acid that gets removed from next window
    for (auto i = std::min(window_size, seq.size()); i < seq.size(); i++) // window slides by one amino acid in each iteration
    {
      sum -= prev_value;
      sum += seq[i].getHydrophobicity(scale); // value for new amino acid gets added
      profile.push_back(sum / std::min(window_size, seq.size()));
      prev_value = seq[i - std::min(window_size, seq.size()) + 1].getHydrophobicity(scale);
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
    if (seq.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Can not compute profile of an empty sequence", "");
    }
    if (window_size<=0)
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,window_size,"window size can not be 0");
    }
    if (window_size > seq.size())
    {
      OPENMS_LOG_WARN << "Warning: window size is larger than sequence size. Setting window size to sequence size: " << seq.size() << "\n";
    }
    std::vector<double> profile;
    for (auto i = 0; i <= seq.size() - std::min(window_size, seq.size()); i++) // window slides by one amino acid in each iteration
    {
      AASequence subsequence = seq.getSubsequence(i, std::min(window_size, seq.size()));
      double sum_sin = 0;
      double sum_cos = 0;
      double curr_position = 0; // position in window
      for (const auto& residue : subsequence) // sum the values in the window
      {
        sum_sin += residue.getHydrophobicity(HydrophobicityScaleMethod::EISENBERG)*std::sin((angle*curr_position*3.14159265359)/180);
        sum_cos += residue.getHydrophobicity(HydrophobicityScaleMethod::EISENBERG)*std::cos((angle*curr_position*3.14159265359)/180);
        curr_position++;
      }
      sum_sin = std::pow(sum_sin,2);
      sum_cos = std::pow(sum_cos,2);
      profile.push_back(std::pow(sum_sin+sum_cos,0.5) / std::min(window_size, seq.size()));
    }
    return profile;
  }
} // namespace OpenMS