// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/Ribonucleotide.h>

#include <vector>
#include <map>
#include <set>
#include <functional>

namespace OpenMS
{
  /*
   * @brief This class applies fixed and variable modifications to (unmodified)
   * nucleic acid sequences, combinatorially generating modified sequences.
   *
   */
  class OPENMS_DLLAPI ModifiedNASequenceGenerator
  {
  public:
    using ConstRibonucleotidePtr = const Ribonucleotide*;

    /// Applies fixed modifications to a single NASequence
    static void applyFixedModifications(
      const std::set<ConstRibonucleotidePtr>& fixed_mods,
      NASequence& sequence);

    /**
       @brief Applies variable modifications to a single NASequence.

       @param var_mods Set of variable modifications to consider
       @param seq Input sequence to modify
       @param max_var_mods Max. number of var. mods that can be added
       @param all_modified_seqs Output list of modified sequences
       @param keep_unmodified Include the original (unmodified) sequence in the results?
       @param max_missed_cleavages Max. number of missed cleavages that can be introduced (by adding inosine)

       @p max_missed_cleavages is ignored if the value is negative.
    */
    static void applyVariableModifications(
      const std::set<ConstRibonucleotidePtr>& var_mods,
      const NASequence& seq, Size max_var_mods,
      std::vector<NASequence>& all_modified_seqs,
      bool keep_unmodified = true, Int max_missed_cleavages = -1);

  protected:
    /// Sequence with modification information
    struct ModSeqInfo
    {
      NASequence seq; //< sequence (unmodified or partially modified)
      Size var_mods_left; //< number of var. mods that can still be applied
      Int missed_cleavages_left; //< for inosine/RNase T1 special case

      /// C'tor for convenience
      ModSeqInfo(NASequence seq, Size var_mods_left,
                 Int missed_cleavages_left = -1):
        seq(seq), var_mods_left(var_mods_left),
        missed_cleavages_left(missed_cleavages_left)
      {
      }
    };

    /**
       @brief Add a modification to a set of sequences using a functor

       @param temp_seqs List of sequences that can be modified
       @param n_temp_seqs Number of sequences from the list that should be modified
       @param finished_seqs Output list of sequences that cannot be further modified
       @param applyMod Functor that modifies a sequence (and optionally checks/updates missed cleavages), returns success status

       Only the first @p n_temp_seqs from @p temp_seqs will be modified.
       After the modification is applied, each new sequence is appended either to @p temp_seqs or @p finished_seqs, depending on whether the max. number of variable mods has been reached.
    */
    static void addModToSequences_(
      std::vector<ModSeqInfo>& temp_seqs, Size n_temp_seqs,
      std::vector<NASequence>& finished_seqs,
      const std::function<bool(NASequence&, Int&)>& applyMod);
  };
}
