// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace OpenMS
{
/**
  @brief Localizes a cross-link by the chemistry of the cross-linker, within the residues the fragment ions support.

  Used for cross-linkers that bind the protein through an NHS ester (e.g. SDA): these react with the
  protein N-terminus and lysines, and much less with serine, threonine and tyrosine (mainly next to a
  histidine, see Mädler et al., J. Mass Spectrom. 2009, 44, 694-706).

  Localization works in two steps:
  1. @ref supportedRange determines the residues the fragment ions can't tell apart (the ambiguity
     stretch, like the lower-case range in the MSFragger localization).
  2. @ref localizeNHS picks the site in that stretch by the NHS chemistry.

  The fragment evidence is passed per residue, as summed intensities of the a/b- and y-ions:
  index @c i of @p n_shifts / @p n_noshifts is the prefix ion that ends at residue @c i (b<sub>i+1</sub>),
  index @c i of @p c_shifts / @p c_noshifts is the suffix ion that starts at residue @c i
  (y<sub>length-i</sub>). "shifts" are ions carrying a cross-linker fragment adduct, "noshifts" ions
  without. All residue positions are 0-based.

  @ingroup Analysis_ID
*/
class NuXLSiteLocalization
{
public:
  /// Result of localizing one cross-link
  struct Site
  {
    int position = -1;                      ///< localized residue, -1 if unlocalized
    int range_begin = -1;                   ///< first residue supported by the fragment ions, -1 without fragment evidence
    int range_end = -1;                     ///< last residue supported by the fragment ions (inclusive), -1 without fragment evidence
    std::string method = "unlocalized";     ///< how the site was chosen: "N-term", "K", "STY+H" or "unlocalized"
  };

  /**
    @brief Residues at which a cross-link explains the fragment ions best.

    For every residue, sums the intensity of the ions consistent with a cross-link there (shifted ions
    that contain the residue, unshifted ions that don't) minus the intensity of the inconsistent ones.
    The returned range spans the residues with the highest sum: residues no ion separates get the
    same sum (e.g., the first two residues if neither b<sub>1</sub> nor y<sub>length-1</sub> is observed).

    @param[in] n_shifts intensity of the shifted prefix (a/b) ions ending at each residue
    @param[in] n_noshifts intensity of the unshifted prefix (a/b) ions ending at each residue
    @param[in] c_shifts intensity of the shifted suffix (y) ions starting at each residue
    @param[in] c_noshifts intensity of the unshifted suffix (y) ions starting at each residue
    @return 0-based first and last residue (inclusive); (-1, -1) if no shifted ion was observed
  */
  static std::pair<int, int> supportedRange(const std::vector<double>& n_shifts,
                                            const std::vector<double>& n_noshifts,
                                            const std::vector<double>& c_shifts,
                                            const std::vector<double>& c_noshifts);

  /**
    @brief Localizes a cross-link attached through an NHS ester.

    Within the @ref supportedRange, in this order:
    1. the peptide N-terminus, if @p reactive_protein_n_term (the peptide starts the protein and its
       N-terminus is free) and the range starts at the first residue;
    2. the lysine closest to the middle of the range;
    3. the serine, threonine or tyrosine closest to the middle of the range that has a histidine within
       four residues in the peptide;
    4. otherwise the cross-link stays unlocalized.

    Ties go to the residue closer to the N-terminus.

    @param[in] sequence unmodified peptide sequence
    @param[in] reactive_protein_n_term whether the peptide N-terminus is the free N-terminus of a protein
    @param[in] n_shifts intensity of the shifted prefix (a/b) ions ending at each residue
    @param[in] n_noshifts intensity of the unshifted prefix (a/b) ions ending at each residue
    @param[in] c_shifts intensity of the shifted suffix (y) ions starting at each residue
    @param[in] c_noshifts intensity of the unshifted suffix (y) ions starting at each residue
  */
  static Site localizeNHS(const std::string& sequence,
                          bool reactive_protein_n_term,
                          const std::vector<double>& n_shifts,
                          const std::vector<double>& n_noshifts,
                          const std::vector<double>& c_shifts,
                          const std::vector<double>& c_noshifts);

  /**
    @brief Sorted protein starts to look up protein N-terminal peptides with @ref isProteinNTerminal.

    For every protein, its sequence and, if it starts with methionine, its sequence without the
    initiator methionine.
  */
  static std::vector<std::string> proteinNTermini(const std::vector<std::string_view>& protein_sequences);

  /// Whether @p peptide starts a protein, i.e. is a prefix of one of the sorted @p protein_n_termini
  static bool isProteinNTerminal(const std::vector<std::string>& protein_n_termini, std::string_view peptide);
};
} // namespace OpenMS
