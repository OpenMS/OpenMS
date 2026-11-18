// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

#include <string>
#include <string_view>
#include <vector>

namespace OpenMS
{
/**
  @brief Localizes a cross-link by the chemistry of the cross-linker, within the residues the fragment ions support.

  Used for cross-linkers that bind the protein through an NHS ester (e.g. SDA): these react with the
  protein N-terminus and lysines, and much less with serine, threonine and tyrosine (mainly next to a
  histidine, see Mädler et al., J. Mass Spectrom. 2009, 44, 694-706).

  Localization works in two steps:
  1. @ref supportedResidues determines the residues the fragment ions can't tell apart (the ambiguity
     stretch, like the lower-case range in the MSFragger localization).
  2. @ref localizeNHS picks the site in that stretch by the NHS chemistry.

  For cross-linkers without such a preference (e.g. UV cross-linking), @ref consistencyScores and
  @ref localizeByConsistency localize by the fragment ions alone, allowing for ions that lost the
  nucleic acid.

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
  /// How cross-links are localized (option NuXL:localization)
  enum class Method
  {
    FRAGMENTS,  ///< at the residue the shifted fragment ions point to (NuXLAnnotateAndLocate)
    NHS,        ///< by the NHS-ester chemistry within the residues the fragment ions support (@ref localizeNHS)
    CONSISTENCY ///< at the residue most consistent with all fragment ions, allowing for loss of the nucleic acid (@ref localizeByConsistency)
  };

  /// Default fraction of the fragment intensity by which a residue's support may fall short of the best and still count as supported
  static constexpr double DEFAULT_TOLERANCE = 0.1;

  /// Default weight of an unshifted ion on a rung without a shifted ion, relative to a shifted ion (see @ref consistencyScores)
  static constexpr double DEFAULT_UNSHIFTED_WEIGHT = 0.5;

  /// Result of localizing one cross-link
  struct Site
  {
    int position = -1;                      ///< localized residue, -1 if unlocalized
    int range_begin = -1;                   ///< first residue supported by the fragment ions, -1 without fragment evidence
    int range_end = -1;                     ///< last residue supported by the fragment ions, -1 without fragment evidence
    std::string method = "unlocalized";     ///< how the site was chosen: "N-term", "K", "STY+H" (@ref localizeNHS), "consistency" (@ref localizeByConsistency) or "unlocalized"
    double delta = 0.0;                     ///< @ref localizeByConsistency only: score of the site minus the best score of another candidate residue; 0 if they tie or there is no other candidate
  };

  /**
    @brief Residues at which a cross-link explains the fragment ions best.

    For every residue, sums the intensity of the ions consistent with a cross-link there (shifted ions
    that contain the residue, unshifted ions that don't) minus the intensity of the inconsistent ones.
    Returns the residues whose sum is at most @p tolerance times the total fragment intensity below the
    highest one, so that a single weak or misassigned ion doesn't decide. They are usually contiguous;
    contradicting ions can support residues on both sides of unsupported ones. Residues no
    ion separates get the same sum (e.g., the first two residues without a<sub>1</sub>/b<sub>1</sub> and
    y<sub>length-1</sub> ions). Unshifted ions alone also limit the range. Unshifted a<sub>1</sub>/b<sub>1</sub>
    ions are ignored: b<sub>1</sub> of lysine coincides with common lysine ions of any peptide.

    @param[in] n_shifts intensity of the shifted prefix (a/b) ions ending at each residue
    @param[in] n_noshifts intensity of the unshifted prefix (a/b) ions ending at each residue
    @param[in] c_shifts intensity of the shifted suffix (y) ions starting at each residue
    @param[in] c_noshifts intensity of the unshifted suffix (y) ions starting at each residue
    @param[in] tolerance fraction of the total fragment intensity by which a residue's sum may fall short of the highest
    @return supported residues (0-based, ascending); empty without any fragment evidence
  */
  static std::vector<int> supportedResidues(const std::vector<double>& n_shifts,
                                            const std::vector<double>& n_noshifts,
                                            const std::vector<double>& c_shifts,
                                            const std::vector<double>& c_noshifts,
                                            double tolerance = DEFAULT_TOLERANCE);

  /**
    @brief Localizes a cross-link attached through an NHS ester.

    Among the @ref supportedResidues, in this order:
    1. the peptide N-terminus, if @p reactive_protein_n_term (the peptide starts the protein and its
       N-terminus is free) and the first residue is supported;
    2. the lysine closest to the middle of the supported residues;
    3. the serine, threonine or tyrosine closest to the middle of the supported residues that has a
       histidine within four residues in the peptide;
    4. otherwise the cross-link stays unlocalized.

    Ties go to the residue closer to the N-terminus. Residues flagged in @p blocked are skipped.

    @param[in] sequence unmodified peptide sequence
    @param[in] reactive_protein_n_term whether the peptide N-terminus is the free N-terminus of a protein
    @param[in] n_shifts intensity of the shifted prefix (a/b) ions ending at each residue
    @param[in] n_noshifts intensity of the unshifted prefix (a/b) ions ending at each residue
    @param[in] c_shifts intensity of the shifted suffix (y) ions starting at each residue
    @param[in] c_noshifts intensity of the unshifted suffix (y) ions starting at each residue
    @param[in] tolerance see @ref supportedResidues
    @param[in] blocked per residue, whether a modification blocks its reactive group (e.g. an acetylated lysine); empty: none
  */
  static Site localizeNHS(const std::string& sequence,
                          bool reactive_protein_n_term,
                          const std::vector<double>& n_shifts,
                          const std::vector<double>& n_noshifts,
                          const std::vector<double>& c_shifts,
                          const std::vector<double>& c_noshifts,
                          double tolerance = DEFAULT_TOLERANCE,
                          const std::vector<bool>& blocked = {});

  /**
    @brief For every residue, how consistent the fragment ions are with a cross-link there, allowing for ions that lost the nucleic acid.

    Each a/b- or y-ion counts with its intensity:
    - a shifted ion for the residues it contains and against the others;
    - an unshifted ion on a rung (ion index) that also has a shifted ion not at all: often only some of the ions of
      a rung lose the nucleic acid, so the unshifted one doesn't contradict a cross-link on the residues it contains;
    - an unshifted ion on a rung without a shifted ion against the residues it contains, weighted by
      @p unshifted_weight, as it, too, may have lost the nucleic acid.

    Unshifted a<sub>1</sub>/b<sub>1</sub> ions are ignored (see @ref supportedResidues). Unlike @ref supportedResidues,
    an unshifted ion doesn't count for the residues it doesn't contain: these residues are not separated from
    each other by it.

    @param[in] n_shifts intensity of the shifted prefix (a/b) ions ending at each residue
    @param[in] n_noshifts intensity of the unshifted prefix (a/b) ions ending at each residue
    @param[in] c_shifts intensity of the shifted suffix (y) ions starting at each residue
    @param[in] c_noshifts intensity of the unshifted suffix (y) ions starting at each residue
    @param[in] unshifted_weight weight of an unshifted ion on a rung without a shifted ion, relative to a shifted ion
    @return per residue, the summed intensity of the consistent minus the inconsistent ions, as a fraction of the
            (weighted) intensity of all ions counted: 1 if all of them are consistent, -1 if all contradict it;
            empty without any fragment evidence
  */
  static std::vector<double> consistencyScores(const std::vector<double>& n_shifts,
                                               const std::vector<double>& n_noshifts,
                                               const std::vector<double>& c_shifts,
                                               const std::vector<double>& c_noshifts,
                                               double unshifted_weight = DEFAULT_UNSHIFTED_WEIGHT);

  /**
    @brief Localizes a cross-link at the residue with the highest @ref consistencyScores.

    Ties go to the residue closer to the N-terminus: residues that no ion separates get exactly the same score.
    The range spans the residues whose score is at most @p tolerance below that of the site, so that the residues
    weak or misassigned ions alone separate count as ambiguous. Residues flagged in @p blocked are no candidates
    (e.g. a C-terminal lysine or arginine the enzyme would have cleaved after).

    @param[in] scores per residue, from @ref consistencyScores
    @param[in] tolerance by how much a residue's score may fall short of the best to be in the range
    @param[in] blocked per residue, whether it can't carry the cross-link; empty: none
    @return the site, with method "consistency"; unlocalized without @p scores or candidate
  */
  static Site localizeByConsistency(const std::vector<double>& scores,
                                    double tolerance = DEFAULT_TOLERANCE,
                                    const std::vector<bool>& blocked = {});

  /**
    @brief The peptides of a protein as searched: its digest and, optionally, its N-terminal peptides after removal of the initial methionine.

    @param[in] digestion the enzyme and missed cleavages of the search
    @param[in] protein protein sequence
    @param[in] min_length minimum peptide length
    @param[in] max_length maximum peptide length (0: no limit)
    @param[in] clip_nterm_methionine also the N-terminal products of the protein without its initial methionine, if it starts with M
    @param[out] peptides views into @p protein
  */
  static void digestProtein(const EnzymaticDigestion& digestion,
                            std::string_view protein,
                            size_t min_length,
                            size_t max_length,
                            bool clip_nterm_methionine,
                            std::vector<std::string_view>& peptides);

  /// Protein N- and C-terminal peptides of a digest, sorted, to look them up with @ref isProteinNTerminal and @ref isProteinCTerminal
  struct ProteinTermini
  {
    std::vector<std::string> n_terminal; ///< peptides that start a protein (also right behind a clipped initial methionine)
    std::vector<std::string> c_terminal; ///< peptides that end a protein
  };

  /**
    @brief The protein-terminal peptides among the peptides @ref digestProtein yields for @p proteins.

    Only digestion products count: a peptide that merely is a prefix or suffix of a protein, without being produced from its terminus, is not protein-terminal.
  */
  static ProteinTermini proteinTermini(const std::vector<std::string_view>& proteins,
                                       const EnzymaticDigestion& digestion,
                                       size_t min_length,
                                       size_t max_length,
                                       bool clip_nterm_methionine);

  /// Whether @p peptide is a protein N-terminal digestion product (see @ref proteinTermini)
  static bool isProteinNTerminal(const ProteinTermini& termini, std::string_view peptide);

  /// Whether @p peptide is a protein C-terminal digestion product (see @ref proteinTermini)
  static bool isProteinCTerminal(const ProteinTermini& termini, std::string_view peptide);
};
} // namespace OpenMS
