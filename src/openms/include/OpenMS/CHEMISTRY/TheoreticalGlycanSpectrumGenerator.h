// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/GlycanStructure.h>
#include <map>

namespace OpenMS
{
class AASequence;
class MSSpectrum;

/**
  @brief Generate positive-mode glycan and singly glycosylated peptide fragments.

  Composition fragments are explicitly labelled as such: their compositions are
  possibilities, not claims about topology. Structural fragments are connected
  subtrees. A B fragment carries a root cleavage and, for internal fragments,
  additional branch cleavages. Y fragments retain the reducing end. Node indices
  identify cleaved bonds by their non-reducing-side node (zero is the attachment).

  Residue formulas exclude water: B = residue sum, C = B + H2O,
  Y = reducing end + retained residues, Z = Y - H2O per branch cleavage.
  Free glycans use H2O as the reducing end; glycopeptides use the peptide.
  Positive charges add protons. Intensities are uniformly 1, not predicted.

  mzPAF does not define glycan series. getAnnotation() therefore uses its named
  compound syntax, keeping topology/composition and attachment in the name.
  Resource limits fail explicitly rather than returning a truncated spectrum.

  @ingroup Chemistry
  @see glycan_fragmentation
*/
class OPENMS_DLLAPI TheoreticalGlycanSpectrumGenerator
{
public:
  using Composition = ProForma::GlycanComposition;

  /// Glycan series and peptide backbone fragments are kept distinct.
  enum class IonType
  {
    DIAGNOSTIC,
    B,
    C,
    Y,
    Z,
    PEPTIDE
  };
  /// EThcD combines b/y and c/z ions with intact glycan retention.
  enum class FragmentationMethod
  {
    HCD,
    ETD,
    ETHCD
  };

  /// Independent choices of glycan material retained by a backbone ion.
  struct OPENMS_DLLAPI PeptideRetention
  {
    bool intact = false;            ///< Include the complete glycan
    bool stripped = true;           ///< Include the peptide with no glycan
    std::vector<Composition> stubs; ///< Explicit retained compositions, each a subset of the input
  };

  /// Generation settings; all size and charge ranges are inclusive.
  struct OPENMS_DLLAPI Options
  {
    bool add_diagnostic_ions = true;
    bool add_b_ions = true;
    bool add_y_ions = true;
    bool add_c_ions = false;
    bool add_z_ions = false;
    bool add_internal_fragments = true;
    bool allow_structural = true;  ///< False selects the composition fallback even for a tree
    Size min_composition_size = 1; ///< Retained residues for both composition B and Y
    Size max_composition_size = 3; ///< Y0 is included separately when Y is enabled
    Size max_cleavages = 2;        ///< Includes the root cleavage of B/internal fragments
    Size max_fragments = 10000;    ///< Includes charge states, losses and peptide fragments
    Size max_states = 100000;      ///< Also bounds enumeration work, including rejected candidates
    Int min_charge = 1;
    Int max_charge = 2;
    Int min_oxonium_charge = 1;
    Int max_oxonium_charge = 1;
    std::vector<EmpiricalFormula> neutral_losses;                                 ///< Single losses, not a combinatorial loss ladder
    std::map<std::string, std::vector<EmpiricalFormula>> specific_neutral_losses; ///< Applied only if that residue is retained
    std::map<char, PeptideRetention> peptide_retention;                           ///< Optional overrides for b/y/c/z
  };

  /// One interpretation of a fragment, including its charge and attachment.
  struct OPENMS_DLLAPI Fragment
  {
    IonType ion_type = IonType::DIAGNOSTIC;
    Composition composition;
    double neutral_mass = 0.0; ///< Before protonation; includes any neutral loss
    Int charge = 1;
    std::optional<Size> attachment_position; ///< Zero-based peptide index
    std::string attachment_residue;
    std::optional<Size> root_cleavage;
    std::vector<Size> branch_cleavages;
    std::string name; ///< Nonempty uncharged name without whitespace, square brackets or NUL characters
    double getMZ() const;
    /// Serialize as an mzPAF named compound; throws InvalidParameter for an invalid name.
    std::string getAnnotation() const;
  };

  /// Construct with default options.
  TheoreticalGlycanSpectrumGenerator();
  /// Construct with validated options.
  explicit TheoreticalGlycanSpectrumGenerator(const Options& options);
  /// Set options; throws InvalidParameter for invalid ranges, loss formulas or rules.
  void setOptions(const Options& options);
  /// Return the current options.
  const Options& getOptions() const;

  /**
    @brief Generate diagnostic and bounded composition B/Y (optionally C/Z) ions.
    @param[in] composition ProForma composition, with nonnegative counts
    @return Fragments sorted by m/z; empty input produces no ions
    @throws Exception::ElementNotFound if a named monosaccharide is unknown to MonosaccharideDB
    @throws Exception::InvalidParameter for invalid chemistry or exhausted resource limits
  */
  std::vector<Fragment> getFragments(const Composition& composition) const;

  /**
    @brief Generate structural B/Y/C/Z and internal B fragments from a tree.
    @param[in] structure Rooted glycan tree
    @return Fragments sorted by m/z, including topology-preserving isobaric interpretations
  */
  std::vector<Fragment> getFragments(const GlycanStructure& structure) const;

  /**
    @brief Generate glycan and peptide fragments for one localized glycan.

    The peptide must exclude the glycan modification; other peptide modifications
    are preserved. HCD defaults to b/y ions stripped or with one HexNAc (if present).
    ETD defaults to c/z ions with the full glycan. EThcD generates b/y/c/z
    ions with the full glycan; per-series overrides can add stripped/stub ions.
    Backbone ions named peptide:zN use the radical z+1 form (Residue::Zp1Ion)
    for both ETD and EThcD; this is distinct from the glycan Z series.
    Only backbone fragments containing the attachment receive retained glycan mass.
    Explicit stubs allow, for example, HexNAc+Fuc retention without assuming a core
    topology from a composition. ETD omits glycan cleavage ions by default.

    @param[in] peptide Peptide without its glycan modification
    @param[in] composition Glycan composition
    @param[in] attachment_position Zero-based residue index in the peptide
    @param[in] method Fragmentation method
    @return Fragments sorted by m/z
    @throws Exception::ElementNotFound if a named monosaccharide is unknown to MonosaccharideDB
    @throws Exception::InvalidParameter for an empty glycan, invalid site, or impossible stub
  */
  std::vector<Fragment>
  getGlycopeptideFragments(const AASequence& peptide, const Composition& composition, Size attachment_position, FragmentationMethod method) const;
  /// As above, using topology for glycan fragments.
  std::vector<Fragment>
  getGlycopeptideFragments(const AASequence& peptide, const GlycanStructure& structure, Size attachment_position, FragmentationMethod method) const;

  /**
    @brief Convert fragment interpretations to a fresh annotated MS2 spectrum.
    @param[in] fragments Fragments to convert
    @return Sorted spectrum with aligned IonNames (mzPAF) and Charges arrays
  */
  static MSSpectrum toSpectrum(const std::vector<Fragment>& fragments);

private:
  Options options_;
  std::vector<Fragment>
  generate_(const Composition& composition, const GlycanStructure* structure, const AASequence* peptide, Size site, FragmentationMethod method) const;
};
} // namespace OpenMS
