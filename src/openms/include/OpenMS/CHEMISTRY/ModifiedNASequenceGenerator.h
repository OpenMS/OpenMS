// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

namespace OpenMS
{
  /*
   * @brief This class applies fixed and variable modifications to (unmodified)
   * nucleic acid sequences, combinatorically generating modified sequences.
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

    /// Applies variable modifications to a single NASequence. If keep_original is set the original (e.g. unmodified version) is also returned
    static void applyVariableModifications(
      const std::set<ConstRibonucleotidePtr>& var_mods,
      const NASequence& seq, Size max_variable_mods_per_NASequence,
      std::vector<NASequence>& all_modified_NASequences,
      bool keep_original = true);

  protected:
    /// Recursively generate all combinatorial placements at compatible sites
    static void recurseAndGenerateVariableModifiedSequences_(
      const std::vector<int>& subset_indices,
      const std::map<int, std::vector<ConstRibonucleotidePtr>>& map_compatibility,
      int depth,
      const NASequence& current_NASequence,
      std::vector<NASequence>& modified_NASequences);

    /// Fast implementation of modification placement. No combinatorial placement is needed in this case
    /// - just every site is modified once by each compatible modification.
    /// Already modified residues are skipped
    static void applyAtMostOneVariableModification_(
      const std::set<ConstRibonucleotidePtr>& var_mods,
      const NASequence& seq,
      std::vector<NASequence>& all_modified_NASequences,
      bool keep_original = true);
  };
}

