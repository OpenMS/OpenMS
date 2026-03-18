// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Copilot $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace OpenMS
{

/**
 * @brief Helper struct to represent a Comet variable modification entry.
 *
 * This struct facilitates merging of modifications with the same mass and compatible
 * term specificity. Comet supports specifying multiple residues in a single modification
 * entry (e.g., "KR" for lysine and arginine), which significantly improves search performance.
 *
 * Merging rules based on Comet's internal logic and documentation:
 * - Modifications with the same mass and binary group can be merged
 * - Amino acid mods can merge with each other (e.g., "K" + "R" -> "KR")
 * - Peptide N-term and amino acids can be merged (e.g., "n" + "RST" -> "nRST")
 * - Peptide C-term and amino acids can be merged (e.g., "c" + "KR" -> "cKR")
 * - Protein terminal mods (nc_term 0 or 1) must remain separate
 * - Peptide N-term and C-term cannot merge with each other
 *
 * When terminal and amino acid mods merge, term_distance becomes -1 and nc_term
 * becomes 0, because the terminal specificity is encoded by 'n'/'c' in the residue
 * string (as per Comet documentation: e.g., "42.010565 nK 0 3 -1 0 0 0.0").
 */
struct OPENMS_DLLAPI CometModification
{
  /// Mass tolerance for comparing modification masses
  static constexpr double MASS_TOLERANCE = 1e-6;

  double mass{0.0};              ///< Modification mass difference
  String residues;               ///< Residue(s) this modification applies to (e.g., "K", "KR", "n", "nKR")
  int binary_group{0};           ///< Binary modification group (for SILAC etc.)
  int max_mods_per_peptide{5};   ///< Maximum number of this modification per peptide
  int term_distance{-1};         ///< Terminal distance (-1 = no constraint, 0 = terminal only)
  int nc_term{0};                ///< Terminal specificity: 0=protein N-term, 1=protein C-term, 2=peptide N-term, 3=peptide C-term
  bool required{false};          ///< Whether this modification is required

  /// Default constructor
  CometModification() = default;

  /// Constructor from ResidueModification
  CometModification(const ResidueModification* mod, int binary_grp, int max_mods);

  /// Check if this modification can be merged with another
  bool isMergeableWith(const CometModification& other) const;

  /// Merge another modification into this one
  void merge(const CometModification& other);

  /// Generate the Comet param file line format for the given 1-based index
  String toCometString(Size index) const;

  /// Merge a vector of CometModifications, combining compatible entries
  static std::vector<CometModification> mergeModifications(const std::vector<CometModification>& mods);
};

} // namespace OpenMS
