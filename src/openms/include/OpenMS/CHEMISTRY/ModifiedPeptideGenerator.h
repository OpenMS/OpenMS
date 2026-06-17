// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>

#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>

namespace OpenMS
{
  /**
    @brief Generate fixed- and variable-modification variants of an @ref AASequence, lock-free.

    The class implements the modification-placement stage of database-search engines and assay
    generators:

      -# getModifications() — resolve a list of UniMod names (e.g. @c "Oxidation (M)") into a
         cached @ref MapToResidueType keyed by @ref ResidueModification. This step performs the
         only @ref ResidueDB lookups; the apply-helpers below are lock-free.
      -# applyFixedModifications() — apply mandatory modifications in place to one peptide.
      -# applyVariableModifications() — enumerate all peptide variants obtained by placing up
         to @c max_variable_mods_per_peptide compatible variable modifications.

    Why pre-caching matters: @ref ResidueDB is process-wide and locked on every getResidue /
    getModifiedResidue call. In a peptide-search workload that means contention proportional
    to the candidate count. getModifications() does the lookups once up front; the apply
    methods then reuse the cached @ref Residue pointers and never touch the DB.

    All methods are static; the class is a namespace in disguise.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI ModifiedPeptideGenerator
  {
  public:
    /**
      @brief Cached mapping ResidueModification* -> already-instantiated modified Residue*.

      Built once by getModifications() (or createResidueModificationToResidueMap_()) and
      consumed by the apply* methods. For modifications that have no associated residue
      (e.g. strict @c "Protein N-term" without an amino-acid origin), the residue pointer
      is @c nullptr.

      The wrapping struct (rather than a bare @c unordered_map) exists so that pyOpenMS can
      bind this aggregate as a single type without having to template-instantiate the map.
    */
    struct MapToResidueType { std::unordered_map<const ResidueModification*, const Residue*> val; };

    /**
      @brief Resolve modification names from UniMod terms and pre-cache the corresponding modified @ref Residue pointers.

      Performs all @ref ResidueDB / @ref ModificationsDB lookups in one place so that the
      downstream apply methods can run lock-free. For modifications without an amino-acid
      origin (strict terminal modifications) the cached residue pointer is @c nullptr.

      @param[in] modNames List of UniMod modification names (e.g. @c "Oxidation (M)").
      @return Mapping ResidueModification* -> modified Residue* suitable for the apply* methods.
    */
    static MapToResidueType getModifications(const StringList& modNames);

    /**
      @brief Apply all compatible fixed modifications in @p fixed_mods to @p peptide in place.

      Two passes:
        - Strict terminal modifications (@c N_TERM / @c C_TERM with no specific residue
          origin) are installed on the peptide's terminal slots if those are still empty.
        - Each residue is then visited: residues already carrying a modification are skipped
          (fixed mods do not overwrite), residues whose one-letter code matches a fixed mod
          receive that mod. Term-specificity @c ANYWHERE applies anywhere; @c C_TERM /
          @c N_TERM only when the residue is the last / first residue of the peptide.

      @param[in]     fixed_mods Cached fixed-modification table (typically from getModifications()).
      @param[in,out] peptide    Peptide to modify in place.
    */
    static void applyFixedModifications(
      const MapToResidueType& fixed_mods,
      AASequence& peptide);

    /**
      @brief Enumerate all peptide variants obtained by combinatorially placing up to @p max_variable_mods_per_peptide variable modifications.

      Treats each modification site independently and produces every legal combination of
      compatible variable mods up to the cap. The original (unmodified) peptide is included
      in the output iff @p keep_original is @c true. The output @p all_modified_peptides is
      appended to (existing contents are preserved).

      Site compatibility is determined by:
        - Residue one-letter code matching the modification's @c origin.
        - Term-specificity (@c ANYWHERE / @c N_TERM / @c C_TERM) being satisfied by the residue's
          position in the peptide.
        - The residue not already carrying a fixed modification (those slots are skipped).

      The implementation uses a fast specialisation (@ref applyAtMostOneVariableModification_)
      when @p max_variable_mods_per_peptide == 1.

      @param[in]  var_mods                      Cached variable-modification table (from getModifications()).
      @param[in]  peptide                       Source peptide; fixed modifications already applied if desired.
      @param[in]  max_variable_mods_per_peptide Maximum number of variable mods to place; 0 disables enumeration.
      @param[out] all_modified_peptides         Generated variants are appended here (existing entries preserved).
      @param[in]  keep_original                 If @c true, also emit @p peptide unchanged as the first entry.
    */
    static void applyVariableModifications(
     const MapToResidueType& var_mods,
     const AASequence& peptide,
     Size max_variable_mods_per_peptide,
     std::vector<AASequence>& all_modified_peptides,
     bool keep_original=true);

  protected:
    /// Sentinel index used internally to mark a strict N_TERM-only modification placed at the N-terminal residue (distinguishes it from an @c ANYWHERE mod that happens to land there)
    static const int N_TERM_MODIFICATION_INDEX;
    /// Sentinel index used internally to mark a strict C_TERM-only modification placed at the C-terminal residue (distinguishes it from an @c ANYWHERE mod that happens to land there)
    static const int C_TERM_MODIFICATION_INDEX;

    /**
      @brief Build the ResidueModification -> Residue cache used by the apply* methods (the lock-free trick — see class brief).

      For each modification: looks up the corresponding modified @ref Residue in @ref ResidueDB
      and stores the pointer. Strict terminal modifications without an amino-acid origin
      (e.g. @c "Protein N-term" with @c origin == 'X') map to @c nullptr because no residue
      can carry them — the apply* code handles those via the peptide's terminal slots.

      @param[in] mods Modifications to cache (typically the resolved list inside getModifications()).
      @return Cache ready for use by applyFixedModifications() / applyVariableModifications().
    */
    static MapToResidueType createResidueModificationToResidueMap_(const std::vector<const ResidueModification*>& mods);

    /**
      @brief Fast path for the common case where exactly one variable modification is placed per peptide.

      Avoids the combinatoric enumeration of @ref applyVariableModifications(): every
      compatible (residue, modification) pair yields exactly one output peptide; already-modified
      residues are skipped. The original peptide is emitted iff @p keep_original is @c true.

      @param[in]  var_mods              Cached variable-modification table.
      @param[in]  peptide               Source peptide.
      @param[out] all_modified_peptides Generated variants appended here.
      @param[in]  keep_original         If @c true, also emit @p peptide unchanged.
    */
    static void applyAtMostOneVariableModification_(
      const MapToResidueType& var_mods,
      const AASequence& peptide,
      std::vector<AASequence>& all_modified_peptides,
      bool keep_original=true);

  private:
    /// For each modification in @p mods, extend @p original_sequences with a copy carrying that mod at @p idx_to_modify; the first mod in @p mods is applied in-place to the existing entries (avoids one copy)
    static void applyAllModsAtIdxAndExtend_(std::vector<AASequence>& original_sequences, int idx_to_modify, const std::vector<const ResidueModification*>& mods, const MapToResidueType& var_mods);
    /// Install modification @p m on @p current_peptide at residue index @p current_index, looking up the pre-cached modified @ref Residue via @p var_mods; overwrites an existing modification if any
    static void applyModToPep_(AASequence& current_peptide, int current_index, const ResidueModification* m, const MapToResidueType& var_mods);
  };
}
