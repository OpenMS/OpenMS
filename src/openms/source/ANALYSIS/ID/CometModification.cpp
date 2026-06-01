// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Copilot $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/CometModification.h>

namespace OpenMS
{

  CometModification::CometModification(const ResidueModification* mod, int binary_grp, int max_mods)
    : mass(mod->getDiffMonoMass()), binary_group(binary_grp), max_mods_per_peptide(max_mods)
  {
    residues = String(mod->getOrigin());

    auto term_spec = mod->getTermSpecificity();
    if (term_spec == ResidueModification::C_TERM)
    {
      if (mod->getOrigin() == 'X')
      {
        residues = "c";
      }
      term_distance = 0;
      nc_term = 3;
    }
    else if (term_spec == ResidueModification::N_TERM)
    {
      if (mod->getOrigin() == 'X')
      {
        residues = "n";
      }
      term_distance = 0;
      nc_term = 2;
    }
    else if (term_spec == ResidueModification::PROTEIN_N_TERM)
    {
      if (mod->getOrigin() == 'X')
      {
        residues = "n";
      }
      term_distance = 0;
      nc_term = 0;
    }
    else if (term_spec == ResidueModification::PROTEIN_C_TERM)
    {
      if (mod->getOrigin() == 'X')
      {
        residues = "c";
      }
      term_distance = 0;
      nc_term = 1;
    }
  }

  bool CometModification::isMergeableWith(const CometModification& other) const
  {
    // Must have same mass (within floating point tolerance)
    if (std::abs(mass - other.mass) > MASS_TOLERANCE)
    {
      return false;
    }

    // Must have same binary group
    if (binary_group != other.binary_group)
    {
      return false;
    }

    // Protein terminal modifications (nc_term 0 or 1 with term_distance 0)
    // must remain separate — they cannot be merged with anything else.
    // Two protein terminal mods of the same type can merge (same residues combined).
    bool this_is_protein_term = (term_distance == 0 && (nc_term == 0 || nc_term == 1));
    bool other_is_protein_term = (other.term_distance == 0 && (other.nc_term == 0 || other.nc_term == 1));

    if (this_is_protein_term || other_is_protein_term)
    {
      return (nc_term == other.nc_term && term_distance == other.term_distance);
    }

    // For non-protein-terminal modifications:
    // - Regular amino acid mods (term_distance=-1) can merge with each other
    // - Peptide N-term (nc_term=2) can merge with amino acids
    // - Peptide C-term (nc_term=3) can merge with amino acids
    // But peptide N-term cannot merge with peptide C-term
    bool this_is_peptide_nterm = (term_distance == 0 && nc_term == 2);
    bool this_is_peptide_cterm = (term_distance == 0 && nc_term == 3);
    bool other_is_peptide_nterm = (other.term_distance == 0 && other.nc_term == 2);
    bool other_is_peptide_cterm = (other.term_distance == 0 && other.nc_term == 3);

    if ((this_is_peptide_nterm && other_is_peptide_cterm) ||
        (this_is_peptide_cterm && other_is_peptide_nterm))
    {
      return false;
    }

    return true;
  }

  void CometModification::merge(const CometModification& other)
  {
    // Add residues from other (avoiding duplicates)
    for (char c : other.residues)
    {
      if (residues.find(c) == String::npos)
      {
        residues += c;
      }
    }

    // Take the maximum max_mods_per_peptide
    max_mods_per_peptide = std::max(max_mods_per_peptide, other.max_mods_per_peptide);

    // When merging terminal with non-terminal modifications, the terminal
    // specificity is already encoded by 'n'/'c' in the residue string.
    // Per Comet convention (e.g., "42.010565 nK 0 3 -1 0 0 0.0"), set
    // term_distance=-1 and nc_term=0 since nc_term is unused when term_distance=-1.
    bool this_is_terminal = (term_distance == 0);
    bool other_is_terminal = (other.term_distance == 0);

    if (this_is_terminal != other_is_terminal)
    {
      // Mixed terminal + non-terminal: terminal info is in residue chars
      term_distance = -1;
      nc_term = 0;
    }
    // If both are terminal (same type, validated by isMergeableWith): keep as-is
    // If both are non-terminal: keep as-is

    required = required || other.required;
  }

  String CometModification::toCometString(Size index) const
  {
    // Format: variable_modXX = <mass> <residues> <binary_group> <max_mods> <term_distance> <nc_term> <required> <neutral_loss>
    std::ostringstream os;
    os // use default precision (6 significant digits) to match previous CometAdapter behavior
       << "variable_mod"
       << std::setw(2) << std::setfill('0') << index
       << " = "
       << mass << " " << residues << " "
       << binary_group << " "
       << max_mods_per_peptide << " "
       << term_distance << " "
       << nc_term << " "
       << (required ? 1 : 0) << " "
       << "0.0";  // neutral loss — not currently supported
    return String(os.str());
  }

  std::vector<CometModification> CometModification::mergeModifications(const std::vector<CometModification>& mods)
  {
    std::vector<CometModification> merged;
    merged.reserve(mods.size());

    for (const auto& mod : mods)
    {
      bool was_merged = false;
      for (auto& existing : merged)
      {
        if (existing.isMergeableWith(mod))
        {
          existing.merge(mod);
          was_merged = true;
          break;
        }
      }
      if (!was_merged)
      {
        merged.push_back(mod);
      }
    }
    return merged;
  }

} // namespace OpenMS
