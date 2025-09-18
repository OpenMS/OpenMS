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
#include <iostream>

namespace OpenMS
{  
  class AASequence;

  /*
      formula2mass holds the map from empirical formula to mass

      mod_combinations holds the map from empirical formula to (potentially ambigious) nucleotide formulae
      e.g.,: 
       C10H14N5O7P   -> {A}
       C10H14N5O8P   -> {G}
       C18H22N4O16P2 -> { CU-H3N1, UU-H2O1 }
   */
  struct OPENMS_DLLAPI NuXLModificationMassesResult
  {
    struct MyStringLengthCompare
    {
      bool operator () (const std::string & p_lhs, const std::string & p_rhs) const
      {
        const size_t lhsLength = p_lhs.length() ;
        const size_t rhsLength = p_rhs.length() ;
        if(lhsLength == rhsLength)
        {
            return (p_lhs < p_rhs) ; // when two strings have the same
                                     // length, defaults to the normal
                                     // string comparison
        }
        return (lhsLength < rhsLength) ; // compares with the length
      }
    };
    std::map<String, double> formula2mass; ///< empirical formula -> mass

    using NucleotideFormulas = std::set<String, MyStringLengthCompare>;
    using MapSumFormulaToNucleotideFormulas = std::map<String, NucleotideFormulas>;
    MapSumFormulaToNucleotideFormulas mod_combinations; ///< empirical formula -> nucleotide formula(s) (formulas if modifications lead to ambiguities)
  };

  class OPENMS_DLLAPI NuXLModificationsGenerator
  {
    public:
      /* @brief generate all combinations of precursor adducts
         @param target_nucleotides the list of nucleotides: e.g., "U", "C", "G", "A" or "U", "T", "G", "A"
         @param can_xl the set of cross-linkable nucleotides
         @param mappings
         @param modifications additional losses associated with the precursor adduct: e.g., "-H2O"
         @param sequence_restriction only precursor adducts that are substrings of this NA sequence are generated
         @param cysteine_adduct special DTT adduct
         @param max_length maximum oligo length
      */
      static NuXLModificationMassesResult initModificationMassesNA(const StringList& target_nucleotides,
                                                                     const StringList& nt_groups,
                                                                     const std::set<char>& can_xl,
                                                                     const StringList& mappings,
                                                                     const StringList& modifications,
                                                                     String sequence_restriction = "",
                                                                     bool cysteine_adduct = false,
                                                                     Int max_length = 4);
    private:
      /// return true if qery is not in sequence
      static bool notInSeq(const String& res_seq, const String& query);

      static void generateTargetSequences(const String& res_seq, Size param_pos, const std::map<char, std::vector<char> >& map_source2target, StringList& target_sequences);
    };
}

