// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/ID/FragmentIonLikelihoodModel.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <array>
#include <cstdint>
#include <vector>

namespace OpenMS
{
  /**
    @brief Self-trained fragment ion likelihood for protein-nucleic acid cross-links (after ProSE's ion priors).

    Every b and y ion of a hit is matched against its spectrum (closest peak within the fragment tolerance) and gets
    an outcome: absent, or intensity-rank bin x mass-error bin (see FragmentIonLikelihoodModel, RICH contexts).
    The ions fall into four groups, each with its own FragmentIonLikelihoodModel:
    - LINEAR: ions of a hit without nucleic acid adduct
    - XL_UNSHIFTED: ions of a localized cross-link that do not contain the cross-linked residue
    - XL_SHIFTED: ions that contain it; they carry one of the fragment adducts of the cross-linked nucleotide, and
      the best outcome over these adducts counts
    - XL_UNLOCALIZED: ions of a cross-link without localization; the best outcome over no shift and every fragment
      adduct counts

    The noise hypothesis is the hit's sequence reversed except for its C-terminal residue (modifications move with
    their residues, so does the cross-linked residue), matched with the same shift alternatives against the same
    spectrum. Models are trained on the confident hits of one half of the spectra and score the other half.
  */
  class NuXLIonLikelihood
  {
  public:
    using Model = FragmentIonLikelihoodModel;

    enum Group : std::uint8_t { LINEAR = 0, XL_UNSHIFTED = 1, XL_SHIFTED = 2, XL_UNLOCALIZED = 3, GROUPS = 4 };

    /// Matched ions of one hypothesis (hit or its reversed sequence)
    struct Outcomes
    {
      std::vector<Model::Ion> ions;
      std::vector<std::uint8_t> groups; ///< Group of each ion
    };

    /// Matched ions of a hit and of its noise hypothesis
    struct Evidence
    {
      Outcomes signal;
      Outcomes noise; ///< empty if the reversed sequence equals the sequence
    };

    /// Features of one hit
    struct Features
    {
      double llr = 0.0;                    ///< summed log-likelihood ratio of all ions
      double llr_shifted = 0.0;            ///< the part of llr from XL_SHIFTED and XL_UNLOCALIZED ions
      double explained_presence = 0.0;     ///< sum of P(present) of the matched ions / sum over all ions
      double top_predicted_observed = 0.0; ///< fraction of the 6 ions with the highest P(present) that matched
    };

    /// @param tolerance fragment tolerance, @param ppm whether it is in ppm, @param max_fragment_charge highest fragment charge
    NuXLIonLikelihood(double tolerance, bool ppm, int max_fragment_charge = 2);

    /// Compact copies of the spectra the hits refer to (by scan index)
    void setSpectra(const PeakMap& spectra);

    /**
      @brief Match the ions of a hit and of its reversed sequence.

      @param scan Scan index of the hit's spectrum.
      @param peptide The peptide, with fixed and variable modifications but without nucleic acid adduct.
      @param precursor_charge Precursor charge of the hit.
      @param shifts Masses of the fragment adducts of the cross-linked nucleotide; empty for hits without cross-link.
      @param site Position of the cross-linked residue, -1 if not localized (ignored without @p shifts).
    */
    Evidence match(Size scan, const AASequence& peptide, int precursor_charge, const std::vector<double>& shifts, int site) const;

    /// Clear the models of @p fold
    void resetFold(Size fold);

    /// Add a confident hit to the models of @p fold
    void addTraining(Size fold, const Evidence& evidence);

    /// Finalize the models of both folds
    void finalize();

    /// Features of a hit scored by the models of @p fold
    Features score(Size fold, const Outcomes& outcomes) const;

    /// Reversed sequence of the noise hypothesis (all but the C-terminal residue)
    static AASequence reversedSequence(const AASequence& peptide);

  private:
    Outcomes matchHypothesis_(Size scan, const AASequence& peptide, int precursor_charge, const std::vector<double>& shifts, int site) const;

    double tolerance_;
    bool ppm_;
    int max_fragment_charge_;
    Model::PeakLists peaks_;
    TheoreticalSpectrumGenerator tsg_;
    std::array<std::array<Model, GROUPS>, 2> models_;
  };
} // namespace OpenMS
