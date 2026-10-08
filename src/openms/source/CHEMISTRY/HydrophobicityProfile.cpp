// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------
 
#include <OpenMS/CHEMISTRY/HydrophobicityProfile.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/DPeak.h>


namespace OpenMS
{
  double HydrophobicityProfile::computeGRAVY(const AASequence& seq) const
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
  const
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
    const Size window_size,
    const HydrophobicityScaleMethod scale
  )
  const
  {
    if (seq.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Can not compute profile of an empty sequence", "");
    }
    if (window_size == 0)
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,window_size,"window size can not be 0");
    }
    Size effective_window = window_size;
    if (window_size > seq.size())
    {
      OPENMS_LOG_WARN << "Warning: window size (" << window_size << ") is larger than sequence length. Window size clamped to sequence length: " << seq.size() << "\n";
      effective_window = seq.size();
    }
    std::vector<double> profile;
    AASequence subsequence = seq.getSubsequence(0, effective_window); // first window
    double sum = 0;
    for (const auto& residue : subsequence) // calculate score for first window
    {
      sum += residue.getHydrophobicity(scale);
    }
    profile.push_back(sum / effective_window);
    double prev_value = seq[0].getHydrophobicity(scale); // value for amino acid that gets removed from next window
    for (auto i = effective_window; i < seq.size(); i++) // window slides by one amino acid in each iteration
    {
      sum -= prev_value;
      sum += seq[i].getHydrophobicity(scale); // value for new amino acid gets added
      profile.push_back(sum / effective_window);
      prev_value = seq[i - effective_window + 1].getHydrophobicity(scale);
    }
    return profile;
  }

  std::vector<double> HydrophobicityProfile::computeHydrophobicMoment
  (
    const AASequence& seq,
    const Size window_size,
    const double angle  
  )
  const
  {
    if (seq.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Can not compute profile of an empty sequence", "");
    }
    if (window_size == 0)
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,window_size,"window size can not be 0");
    }
    Size effective_window = window_size;
    if (window_size > seq.size())
    {
      OPENMS_LOG_WARN << "Warning: window size (" << window_size << ") is larger than sequence length. Window size clamped to sequence length: " << seq.size() << "\n";
      effective_window = seq.size();
    }
    std::vector<double> profile;
    for (Size i = 0; i + effective_window <= seq.size(); i++) // window slides by one amino acid in each iteration
    {
      AASequence subsequence = seq.getSubsequence(i, effective_window);
      double sum_sin = 0;
      double sum_cos = 0;
      double curr_position = 0; // position in window
      double angle_rad = angle * Constants::PI / 180;
      for (const auto& residue : subsequence) // sum the values in the window
      {
        double hydrophobicity = residue.getHydrophobicity(HydrophobicityScaleMethod::EISENBERG);
        sum_sin += hydrophobicity * std::sin(angle_rad * curr_position);
        sum_cos += hydrophobicity * std::cos(angle_rad * curr_position);
        curr_position++;
      }
      sum_sin *= sum_sin; //square
      sum_cos *= sum_cos; //square
      profile.push_back(std::sqrt(sum_sin + sum_cos) / effective_window);
    }
    return profile;
  }
} // namespace OpenMS