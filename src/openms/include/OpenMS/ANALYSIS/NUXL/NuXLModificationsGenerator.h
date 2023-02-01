// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

