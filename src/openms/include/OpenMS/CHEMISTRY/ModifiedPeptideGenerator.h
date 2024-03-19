// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

#include <OpenMS/CHEMISTRY/ProteaseDigestion.h> // for digestandmodify
#include <OpenMS/DATASTRUCTURES/StringView.h> // for digestandmodify

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

    /**
    * @brief Apply variable modifications to a peptide.
    *
    * This function applies variable modifications to a given peptide based on the provided mapping
    * of modifications to residue types. It generates all possible modified peptide variants
    * considering the specified constraints such as the maximum number of variable modifications
    * allowed per peptide and whether unmodified peptides should be retained.
    *
    * @param var_mods A mapping of modifications to residue types.
    * @param peptide The original peptide to which modifications will be applied.
    * @param max_variable_mods_per_peptide The maximum number of variable modifications allowed per peptide.
    * @param all_modified_peptides Output container. Generated modified peptides are appended.
    * @param keep_unmodified Flag indicating whether unmodified peptides should be retained.
    */    
    static void applyVariableModifications(
     const MapToResidueType& var_mods, 
     const AASequence& peptide, 
     Size max_variable_mods_per_peptide, 
     std::vector<AASequence>& all_modified_peptides, 
     bool keep_original=true);

    struct SequenceMassPair
    {
      SequenceMassPair(const AASequence& s, double m) noexcept :
        mass(m), 
        sequence(s)
        {}
    
      SequenceMassPair(AASequence&& s, double m) noexcept :
        mass(m),
        sequence(std::move(s))
        {}

      double mass{};
      AASequence sequence;
    };

    /**
    * @brief Apply variable modifications to a peptide to obtain the modified sequence and mass.
    *
    * This function applies variable modifications to a given peptide based on the provided mapping
    * of modifications to residue types. It generates all possible modified peptide variants
    * considering the specified constraints such as the maximum number of variable modifications
    * allowed per peptide and whether unmodified peptides should be retained. Additionally, it
    * calculates and stores the monoisotopic mass of each modified peptide.
    *
    * @param var_mods A mapping of modifications to residue types.
    * @param peptide The original peptide to which modifications will be applied.
    * @param max_variable_mods_per_peptide The maximum number of variable modifications allowed per peptide.
    * @param peptide_is_n_terminal Whether the peptide is N-terminal in the protein sequence. Needed to correctly apply terminal modifications.
    * @param peptide_is_c_terminal Whether the peptide is C-terminal. Needed to correctly apply terminal modifications.          
    * @param all_modified_peptides Output container. Generated modified peptides with their masses are appended.
    * @param keep_unmodified Flag indicating whether unmodified peptides should be retained.
    */
    static void generateVariableModifiedPeptidesWithMasses(
     const MapToResidueType& var_mods, 
     const AASequence& peptide,
     bool peptide_is_n_terminal,
     bool peptide_is_c_terminal,          
     Size max_variable_mods_per_peptide,
     std::vector<SequenceMassPair>& all_modified_peptides,
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
    static void applyAllModsAtIdxAndExtend_(std::vector<SequenceMassPair>& original_sequences, int idx_to_modify, const std::vector<const ResidueModification*>& mods, const MapToResidueType& var_mods);
    /// applies a modification @p m to the @p current_peptide at @p current_index. Overwrites mod if it exists. Looks up in var_mods for existing modified Residue pointers.
    static void applyModToPep_(AASequence& current_peptide, int current_index, const ResidueModification* m, const MapToResidueType& var_mods);
  };

  class OPENMS_DLLAPI DigestAndModify
  {
    private:
      ProteaseDigestion digestor_;
      Size max_variable_mods_per_peptide_;
      Size min_length_;
      Size max_length_;
      ModifiedPeptideGenerator::MapToResidueType variable_mods_;
      ModifiedPeptideGenerator::MapToResidueType fixed_mods_;

    public:
    DigestAndModify(const String& enzyme, const Size& missed_cleavages, const StringList& fixed_modifications, const StringList& variable_modifications, Size max_variable_mods_per_peptide = 2, Size min_length = 1, Size max_length = 0) :
       max_variable_mods_per_peptide_(max_variable_mods_per_peptide),
       min_length_(min_length),
       max_length_(max_length)
    { 
      // set up digestion
      digestor_.setEnzyme(enzyme);
      digestor_.setMissedCleavages(missed_cleavages);

      // set up modifications
      fixed_mods_ = ModifiedPeptideGenerator::getModifications(fixed_modifications);
      variable_mods_ = ModifiedPeptideGenerator::getModifications(variable_modifications);
    }

    void getPeptides(const StringView& protein_sequence, std::vector<ModifiedPeptideGenerator::SequenceMassPair>& modified_peptides)
    {
      std::vector<std::pair<Size, Size>> output;
      digestor_.digestUnmodified(protein_sequence, output, min_length_, max_length_);

      for (const auto&[start, end] : output)
      {
        AASequence peptide = AASequence::fromString(protein_sequence.substr(start, end - start).getString());
        ModifiedPeptideGenerator::applyFixedModifications(fixed_mods_, peptide); // this could be done on protein level. not sure if worth the rewrite

        // determine if peptide is C-terminal, N-terminal, both (or none of them) and pass it to ensure modifications are correctly applied
        bool peptide_is_n_terminal = (start == 0);
        bool peptide_is_c_terminal = (end == peptide.size() - 1);
        ModifiedPeptideGenerator::generateVariableModifiedPeptidesWithMasses(variable_mods_, 
          peptide,
          peptide_is_c_terminal,
          peptide_is_n_terminal,
          max_variable_mods_per_peptide_, 
          modified_peptides, 
          true);
      }
      return;
    }

    // Getter for digestor_
    const ProteaseDigestion& getDigestor() const 
    {
      return digestor_;
    }

    // Getter for min_length_
    Size getMinLength() const 
    {
      return min_length_;
    }

    // Getter for max_length_
    Size getMaxLength() const 
    {
      return max_length_;
    }

    // Getter for max_variable_mods_per_peptide_
    Size getMaxVariableModsPerPeptide() const 
    {
      return max_variable_mods_per_peptide_;
    }

    // Getter for variable_mods_
    const ModifiedPeptideGenerator::MapToResidueType& getVariableMods() const 
    {
      return variable_mods_;
    }

    // Getter for fixed_mods_
    const ModifiedPeptideGenerator::MapToResidueType& getFixedMods() const 
    {
      return fixed_mods_;
    }    

  };
}
