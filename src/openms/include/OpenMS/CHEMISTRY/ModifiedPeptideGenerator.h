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
  class OPENMS_DLLAPI ModifiedPeptideGenerator
  {
   /*
    * @brief Modifications can be generated and applied to AASequences. 
    */

  public:
    // struct needed to wrap the template for pyOpenMS
    struct MapToResidueType { std::unordered_map<const ResidueModification*, const Residue*> val; };

      /**
      * @brief Retrieve modifications from strings
      * 
      * @param modNames The list of modification names
      * @return A map of modifications and associated residue
      * ResidueModifications are referenced by Residues in AASequence objects. Every time an AASequence object
      * with modifications is constructed, it needs to query if the (modified) Residue is already
      * registered in ResidueDB. This implies a lock of the whole db. To make modified peptide generation lock-free, we
      * query and cache all modified residues once so we can directly apply them without further queries.
      */
    static MapToResidueType getModifications(const StringList& modNames);

    // Applies fixed modifications to a single peptide
    static void applyFixedModifications(
      const MapToResidueType& fixed_mods, 
      AASequence& peptide);

    // Applies variable modifications to a single peptide. If keep_original is set the original (e.g. unmodified version) is also returned
    static void applyVariableModifications(
     const MapToResidueType& var_mods, 
     const AASequence& peptide, 
     Size max_variable_mods_per_peptide, 
     std::vector<AASequence>& all_modified_peptides, 
     bool keep_original=true);

  protected:
    static const int N_TERM_MODIFICATION_INDEX; // magic constant to distinguish N_TERM only modifications from ANYWHERE modifications placed at N-term residue
    static const int C_TERM_MODIFICATION_INDEX; // magic constant to distinguish C_TERM only modifications from ANYWHERE modifications placed at C-term residue

    // Lookup datastructure to allow lock-free generation of modified peptides. Modifications without origin (e.g., "Protein N-term") set the residue to nullptr.
    static MapToResidueType createResidueModificationToResidueMap_(const std::vector<const ResidueModification*>& mods);

    // Fast implementation of modification placement. No combinatoric placement is needed in this case - just every site is modified once by each compatible modification. Already modified residues are skipped
    static void applyAtMostOneVariableModification_(
      const MapToResidueType& var_mods, 
      const AASequence& peptide, 
      std::vector<AASequence>& all_modified_peptides, 
      bool keep_original=true);

  private:
    /// take a vector of AASequences @p original_sequences, and for each mod in @p mods, add a version with mod at index @p idx_to_modify. In-place, with the original sequences recieving the first mod in @p mods.
    static void applyAllModsAtIdxAndExtend_(std::vector<AASequence>& original_sequences, int idx_to_modify, const std::vector<const ResidueModification*>& mods, const MapToResidueType& var_mods);
    /// applies a modification @p m to the @p current_peptide at @p current_index. Overwrites mod if it exists. Looks up in var_mods for existing modified Residue pointers.
    static void applyModToPep_(AASequence& current_peptide, int current_index, const ResidueModification* m, const MapToResidueType& var_mods);
  };
}
