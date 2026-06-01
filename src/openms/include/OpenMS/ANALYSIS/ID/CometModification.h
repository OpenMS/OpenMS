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
  @brief Helper struct that represents one Comet @c variable_modXX entry and supports merging compatible entries.

  Comet supports specifying multiple residues in a single modification
  entry (e.g. @c "KR" for lysine and arginine), which significantly
  improves search performance. This struct mirrors the @c variable_modXX
  fields and bundles the rules used to collapse compatible entries.

  Merging rules (verified against the cpp):
    - Two entries must have the same mass (within @ref MASS_TOLERANCE) and the same @ref binary_group.
    - Plain amino-acid mods (@c term_distance @c == @c -1) can merge with each other (e.g. @c "K" + @c "R" -> @c "KR").
    - Peptide N-term (@c nc_term @c == @c 2, @c term_distance @c == @c 0) and amino acids can merge (e.g. @c "n" + @c "RST" -> @c "nRST").
    - Peptide C-term (@c nc_term @c == @c 3, @c term_distance @c == @c 0) and amino acids can merge (e.g. @c "c" + @c "KR" -> @c "cKR").
    - Protein terminal mods (@c nc_term @c == @c 0 or @c 1 with @c term_distance @c == @c 0) only merge with another protein terminal mod of the same kind.
    - Peptide N-term and peptide C-term cannot merge with each other.

  When a terminal and an amino-acid mod merge, @c term_distance is set
  to @c -1 and @c nc_term to @c 0 because the terminal specificity is
  encoded by the @c 'n' / @c 'c' character in the residue string (per
  Comet convention, e.g. @c "42.010565 nK 0 3 -1 0 0 0.0").

  @ingroup Analysis_ID
*/
struct OPENMS_DLLAPI CometModification
{
  /// Absolute tolerance used by @ref isMergeableWith when comparing the @ref mass field.
  static constexpr double MASS_TOLERANCE = 1e-6;

  double mass{0.0};              ///< Modification mass difference (in Da).
  String residues;               ///< Residue(s) this modification applies to (e.g. @c "K", @c "KR", @c "n", @c "nKR").
  int binary_group{0};           ///< Comet binary modification group (used e.g. for SILAC).
  int max_mods_per_peptide{5};   ///< Maximum number of occurrences of this modification per peptide.
  int term_distance{-1};         ///< Terminal distance constraint: @c -1 = no constraint, @c 0 = terminal only.
  int nc_term{0};                ///< Terminal specificity: @c 0 = protein N-term, @c 1 = protein C-term, @c 2 = peptide N-term, @c 3 = peptide C-term. Only meaningful when @c term_distance @c == @c 0.
  bool required{false};          ///< Whether this modification is required to appear in every peptide.

  /// Default constructor; produces a zeroed entry with @c mass @c == @c 0.
  CometModification() = default;

  /**
    @brief Build an entry from an OpenMS @ref ResidueModification.

    Copies the modification's diff-monoisotopic mass into @ref mass and
    the residue origin character into @ref residues, then sets
    @ref term_distance and @ref nc_term according to @p mod 's term
    specificity:
      - @c PEPTIDE_C_TERM -> @c term_distance @c = @c 0, @c nc_term @c = @c 3.
      - @c PEPTIDE_N_TERM -> @c term_distance @c = @c 0, @c nc_term @c = @c 2.
      - @c PROTEIN_N_TERM -> @c term_distance @c = @c 0, @c nc_term @c = @c 0.
      - @c PROTEIN_C_TERM -> @c term_distance @c = @c 0, @c nc_term @c = @c 1.
      - Anywhere -> @c term_distance and @c nc_term are left at their defaults.
    Terminal mods whose origin is @c 'X' are translated into a residue
    string of @c "n" or @c "c" so Comet recognises them; non-@c 'X'
    terminal mods keep the original residue character.

    @param[in] mod         Source OpenMS modification.
    @param[in] binary_grp  Value for @ref binary_group.
    @param[in] max_mods    Value for @ref max_mods_per_peptide.
  */
  CometModification(const ResidueModification* mod, int binary_grp, int max_mods);

  /**
    @brief Whether two entries can be combined into a single Comet @c variable_modXX line.

    Implements the merging rules listed in the struct brief. Returns
    @c false when the masses differ by more than @ref MASS_TOLERANCE,
    the @ref binary_group values differ, when only one of the two is a
    protein-terminal mod (or they are different protein-terminal kinds),
    or when one is a peptide N-term and the other a peptide C-term.

    @param[in] other Candidate merge partner.
    @return @c true if the entries can be merged.
  */
  bool isMergeableWith(const CometModification& other) const;

  /**
    @brief Merge another compatible entry into this one.

    Unions @p other 's residue characters into @ref residues (skipping
    duplicates), keeps the larger @ref max_mods_per_peptide, and ORs
    the @ref required flag. When a terminal mod is merged with a
    non-terminal mod the terminal information moves into the residue
    string and @ref term_distance / @ref nc_term are reset to
    @c -1 / @c 0; same-kind terminal merges leave both fields
    unchanged.

    Behaviour is only defined for pairs that @ref isMergeableWith
    accepts; the caller is responsible for that check.

    @param[in] other Compatible entry to absorb.
  */
  void merge(const CometModification& other);

  /**
    @brief Render this entry as a Comet @c variable_modXX parameter line.

    The output is one whitespace-separated line of the form
    @verbatim
    variable_mod<II> = <mass> <residues> <binary_group> <max_mods_per_peptide> <term_distance> <nc_term> <required> 0.0
    @endverbatim
    with the @c II suffix being @p index zero-padded to two digits. The last
    column is hard-coded to @c "0.0" (neutral loss; OpenMS does not
    currently surface it). Numeric formatting uses @c std::ostream
    default precision (six significant digits).

    @param[in] index 1-based modification index used to form the @c variable_modXX key.
    @return The fully-formed @c variable_modXX = ... line.
  */
  String toCometString(Size index) const;

  /**
    @brief Greedy first-fit merge of a vector of entries.

    Walks @p mods in input order; each entry is merged into the first
    compatible existing result entry (via @ref isMergeableWith and
    @ref merge), or appended verbatim if no compatible partner is found.
    The relative order of the kept entries follows the input order; no
    entry is dropped.

    @param[in] mods Input modifications.
    @return Merged modifications.
  */
  static std::vector<CometModification> mergeModifications(const std::vector<CometModification>& mods);
};

} // namespace OpenMS
