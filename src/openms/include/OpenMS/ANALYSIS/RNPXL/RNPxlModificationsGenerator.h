// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <vector>
#include <map>
#include <set>

namespace OpenMS
{  
  class AASequence;

  struct OPENMS_DLLAPI RNPxlModificationMassesResult
  {
    std::map<String, double> mod_masses; ///< empirical formula -> mass
    std::map<String, std::set<String> > mod_combinations; ///< empirical formula -> nucleotide formula(s) (formulas if modifications lead to ambiguities)
  };

  class OPENMS_DLLAPI RNPxlModificationsGenerator
  {
    public:
      static RNPxlModificationMassesResult initModificationMassesRNA(const StringList& target_nucleotides,
                                                                     const StringList& nt_groups,
                                                                     const std::set<char>& can_xl,
                                                                     const StringList& mappings,
                                                                     const StringList& modifications,
                                                                     String sequence_restriction,
                                                                     bool cysteine_adduct,
                                                                     Int max_length = 4);
    private:
      static bool notInSeq(const String& res_seq, const String& query);
      static void generateTargetSequences(const String& res_seq, Size param_pos, const std::map<char, std::vector<char> >& map_source2target, StringList& target_sequences);
    };
}


