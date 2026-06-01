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
  /**
    @brief Generator of modified nucleic-acid sequences for fixed and variable modification placement.

    Counterpart of @ref ModifiedPeptideGenerator for @ref NASequence. The class is a namespace
    in disguise: all entry points are static and operate either on a single sequence in place
    (fixed modifications) or by enumerating every legal combination of placements (variable
    modifications) into a caller-provided output vector.

    Compatibility of a modification with a sequence position is determined by:
      - The residue's one-letter code matching the modification's @c origin character.
      - The modification's term-specificity (@c FIVE_PRIME / @c THREE_PRIME / @c ANYWHERE)
        being satisfied: @c ANYWHERE applies to any internal nucleotide,
        @c FIVE_PRIME / @c THREE_PRIME only at the corresponding sequence terminus.
      - The residue / terminal slot not already carrying a modification — pre-existing
        modifications are preserved and never overwritten.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI ModifiedNASequenceGenerator
  {
  public:
    using ConstRibonucleotidePtr = const Ribonucleotide*;

    /**
      @brief Apply all compatible fixed modifications from @p fixed_mods to @p sequence in place.

      Two passes:
        - Terminal modifications (@c FIVE_PRIME / @c THREE_PRIME) are installed on the
          corresponding terminal slot only when that slot is still empty. Pre-existing
          terminal modifications are preserved.
        - Each nucleotide is then visited: residues already carrying a modification are
          skipped, and residues whose one-letter code matches the modification's @c origin
          receive that modification provided its term-specificity is @c ANYWHERE.

      When multiple fixed modifications match the same residue, the iteration order of
      @p fixed_mods (a pointer-keyed @c std::set) decides which placement persists.

      @param[in]     fixed_mods Fixed modifications to apply.
      @param[in,out] sequence   Sequence modified in place.
    */
    static void applyFixedModifications(
      const std::set<ConstRibonucleotidePtr>& fixed_mods,
      NASequence& sequence);

    /**
      @brief Enumerate all sequence variants obtained by combinatorially placing up to @p max_variable_mods_per_NASequence variable modifications.

      Produces every legal combination of @em compatible variable modifications drawn from
      @p var_mods, with the number of placements per variant ranging from @c 1 to
      @c min(max_variable_mods_per_NASequence, @c N) where @c N is the number of compatible sites. 5' / 3' terminal
      modifications participate in the enumeration as virtual sites and are considered only
      when the corresponding terminal slot of @p seq is empty; nucleotides that are already
      modified are skipped during compatibility analysis.

      The output @p all_modified_NASequences is @em appended to (existing entries are
      preserved). The unmodified @p seq is included in the appended block iff @p keep_original
      is @c true. If @p var_mods is empty, @p max_variable_mods_per_NASequence is @c 0, or no
      compatible sites are found, only the unmodified sequence is appended (and only when
      @p keep_original is @c true). A fast specialisation is used when
      @p max_variable_mods_per_NASequence equals @c 1.

      @param[in]  var_mods                          Variable modifications to consider.
      @param[in]  seq                               Source sequence; not modified.
      @param[in]  max_variable_mods_per_NASequence  Maximum number of variable modifications per variant; @c 0 disables enumeration.
      @param[out] all_modified_NASequences          Generated variants are appended here.
      @param[in]  keep_original                     If @c true, also emit @p seq unchanged.
    */
    static void applyVariableModifications(
      const std::set<ConstRibonucleotidePtr>& var_mods,
      const NASequence& seq, Size max_variable_mods_per_NASequence,
      std::vector<NASequence>& all_modified_NASequences,
      bool keep_original = true);

  protected:
    /**
      @brief Recursive helper for @ref applyVariableModifications — enumerate every assignment of compatible modifications to a fixed subset of sequence sites.

      For the position at @c subset_indices[depth], iterate over every modification listed
      in @p map_compatibility for that position, apply it to a copy of @p current_NASequence,
      and recurse with @c depth+1. When @p depth reaches @c subset_indices.size() the current
      sequence is appended to @p modified_NASequences. Sentinel indices @c -1 / @c -2 select
      the 5' / 3' terminal slot instead of an internal residue.

      @param[in]     subset_indices       Selected positions (and 5'/3' sentinels) at which to place modifications.
      @param[in]     map_compatibility    Per-position list of compatible modifications.
      @param[in]     depth                Current recursion depth (caller passes @c 0).
      @param[in]     current_NASequence   Sequence accumulated so far.
      @param[in,out] modified_NASequences Output vector to which fully enumerated variants are appended.
    */
    static void recurseAndGenerateVariableModifiedSequences_(
      const std::vector<int>& subset_indices,
      const std::map<int, std::vector<ConstRibonucleotidePtr>>& map_compatibility,
      int depth,
      const NASequence& current_NASequence,
      std::vector<NASequence>& modified_NASequences);

    /**
      @brief Fast specialisation of variable-modification placement for at most one modification per sequence.

      Emits one variant per (site, compatible modification) pair: no combinations are
      enumerated. Already-modified residues are skipped. The unmodified sequence is emitted
      first when @p keep_original is @c true.

      @param[in]  var_mods                 Variable modifications to consider.
      @param[in]  seq                      Source sequence; not modified.
      @param[out] all_modified_NASequences Generated variants are appended here.
      @param[in]  keep_original            If @c true, also emit @p seq unchanged before any modified variants.
    */
    static void applyAtMostOneVariableModification_(
      const std::set<ConstRibonucleotidePtr>& var_mods,
      const NASequence& seq,
      std::vector<NASequence>& all_modified_NASequences,
      bool keep_original = true);
  };
}

