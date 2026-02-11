// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Development Team $
// $Authors: Srikanth K N $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/SequenceCoverage.h>

#include <vector>

namespace OpenMS
{

  double SequenceCoverage::getCoverage(
    const AASequence& protein,
    const std::vector<AASequence>& peptides)
  {
    const String protein_str = protein.toUnmodifiedString();
    const Size protein_length = protein.size();

    if (protein_length == 0 || peptides.empty())
    {
      return 0.0;
    }

    std::vector<bool> covered(protein_length, false);

    for (const auto& pep : peptides)
    {
      const String peptide_str = pep.toUnmodifiedString();

      if (peptide_str.empty())
      {
        continue;
      }

      Size pos = protein_str.find(peptide_str);

      while (pos != String::npos)
      {
        for (Size i = pos; i < pos + peptide_str.size(); ++i)
        {
          covered[i] = true;
        }

        pos = protein_str.find(peptide_str, pos + 1);
      }
    }

    Size covered_count = 0;

    for (bool c : covered)
    {
      if (c)
      {
        ++covered_count;
      }
    }

    return static_cast<double>(covered_count) * 100.0 / protein_length;
  }

} // namespace OpenMS
