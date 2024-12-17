// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ITransition.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <OpenMS/KERNEL/RangeManager.h>

namespace OpenMS
{
  class TheoreticalSpectrumGenerator;

  /**
    @brief Scoring of an spectrum at the peak apex of an chromatographic elution peak.

    In DIA (data independent acquisition) / SWATH analysis, at each
    chromatographic point a full MS2 spectrum is recorded. This class allows to
    compute a number of scores based on the full MS2 spectrum available. The scores are the following:

    - isotope scores:
      - isotope_corr: computes the correlation of each fragment ion with the
         theoretical isotope distribution. This is the pearson correlation to
         the theoretical isotope pattern weighted by the relative intensity of
         the transition (more is better).
      - isotope_overlap: checks whether a signal at position (mz - 1) / charge
         exists and how strong it is. This would be an indication that the current
         peak is an isotopic signal of another peak. This simply counts how
         often a peak was observed that is higher than the current peak, thus
         number is then weighted by the relative intensity of the transition
         (thus less is better here).

    - massdiff score: computes the difference in ppm of the experimental signal to the expected signal (thus less is better). 
      - Equation: sum(ppm_difference) / # transitions
      - Notes: 
        - Divide by the total number of transitions and is thus quite punishing if a transition is missing
        - Also outputs a list of all the ppm differences, if signal is not found output -1.0

    - b/y ion score: checks for the presence of b/y ions of the peptide in question

    - theoretical spectrum: a dotproduct and a manhattan score with a theoretical spectrum

    This class expects spectra objects that implement the OpenSWATH Spectrum
    interface. Transitions are expected to be in the light transition format
    (defined in OPENSWATHALGO/DATAACCESS/TransitionExperiment.h).

  @htmlinclude OpenMS_DIAScoring.parameters

  */
  class OPENMS_DLLAPI DIAScoring :
    public DefaultParamHandler
  {
    ///Type definitions
    //@{
    /// Spectrum type, see Spectrum interface
    typedef OpenSwath::SpectrumPtr SpectrumPtrType;
    /// Transition interface (Transition, Peptide, Protein)
    typedef OpenSwath::LightTransition TransitionType;
    //@}

public:

    ///@name Constructors and Destructor
    //@{
    /// Default constructor
    DIAScoring();

    /// Destructor
    ~DIAScoring() override;
    //@}

    ///////////////////////////////////////////////////////////////////////////
    // DIA / SWATH scoring

    ///@name DIA Scores
    //@{
    /// Isotope scores, see class description
    void dia_isotope_scores(const std::vector<TransitionType>& transitions,
                            SpectrumSequence& spectrum,
                            OpenSwath::IMRMFeature* mrmfeature,
                            const RangeMobility& im_range,
                            double& isotope_corr,
                            double& isotope_overlap) const;

    /// Massdiff scores, see class description
    void dia_massdiff_score(const std::vector<TransitionType>& transitions,
                            const SpectrumSequence& spectrum,
                            const std::vector<double>& normalized_library_intensity,
                            const RangeMobility& im_range,
                            double& ppm_score,
                            double& ppm_score_weighted,
                            std::vector<double>& diff_ppm) const;

    /**
      Precursor massdifference score

      @param precursor_mz Exact m/z of the precursor to be evaluated
      @param spectrum MS1 spectrum to be evaluated
      @param im_range Ion mobility range to keep (filter data); can be empty
      @param ppm_score Resulting score
      @return False if no signal was found (and no sensible score calculated), true otherwise
    */
    bool dia_ms1_massdiff_score(double precursor_mz, const SpectrumSequence& spectrum, const RangeMobility& im_range,
                                double& ppm_score) const;

    /// Precursor isotope scores for precursors (peptides and metabolites)
    void dia_ms1_isotope_scores_averagine(double precursor_mz, const SpectrumSequence& spectrum, int charge_state, RangeMobility& im_range,
                                          double& isotope_corr, double& isotope_overlap) const;
    void dia_ms1_isotope_scores(double precursor_mz, const std::vector<SpectrumPtrType>& spectrum, RangeMobility& im_range,
                                double& isotope_corr, double& isotope_overlap, const EmpiricalFormula& sum_formula) const;

    /// b/y ion scores
    void dia_by_ion_score(const SpectrumSequence& spectrum, AASequence& sequence,
                          int charge, const RangeMobility& im_range, double& bseries_score, double& yseries_score) const;

    /// Dotproduct / Manhattan score with theoretical spectrum
    void score_with_isotopes(SpectrumSequence& spectrum,
                             const std::vector<TransitionType>& transitions,
                             const RangeMobility& im_range,
                             double& dotprod,
                             double& manhattan) const;
    //@}

private:

    /// Copy constructor (algorithm class)
    DIAScoring(const DIAScoring& rhs);

    /// Assignment operator (algorithm class)
    DIAScoring& operator=(const DIAScoring& rhs);

    /// Synchronize members with param class
    void updateMembers_() override;

    /// Subfunction of dia_isotope_scores
    void diaIsotopeScoresSub_(const std::vector<TransitionType>& transitions,
                              const SpectrumSequence& spectrum,
                              std::map<std::string, double>& intensities,
                              const RangeMobility& im_range,
                              double& isotope_corr,
                              double& isotope_overlap) const;

    /// retrieves intensities from MRMFeature
    /// computes a vector of relative intensities for each feature (output to intensities)
    void getFirstIsotopeRelativeIntensities_(const std::vector<TransitionType>& transitions,
                                            OpenSwath::IMRMFeature* mrmfeature,
                                            std::map<std::string, double>& intensities //experimental intensities of transitions
                                            ) const;

private:

    /**
      @brief Determine whether the current m/z value is a monoisotopic peak

      This function will try to determine whether the current peak is a
      monoisotopic peak or not. It will do so by searching for an intense peak
      at a lower m/z that could explain the current peak as part of a isotope
      pattern.

      @param spectrum The spectrum (MS1 or MS2)
      @param mono_mz The m/z value where a monoisotopic is expected
      @param mono_int The intensity of the monoisotopic peak (peak at mono_mz)
      @param nr_occurrences Will contain the count of how often a peak is found at lower m/z than mono_mz with an intensity higher than mono_int. Multiple charge states are tested, see class parameter dia_nr_charges_
      @param max_ratio Will contain the maximum ratio of a peaks intensity compared to the monoisotopic peak intensity how often a peak is found at lower m/z than mono_mz with an intensity higher than mono_int. Multiple charge states are tested, see class parameter dia_nr_charges_
      @param im_range Ion mobility subrange to consider (used as filter); can be empty (i.e. no IM filtering)
    */
    void largePeaksBeforeFirstIsotope_(const SpectrumSequence& spectrum, double mono_mz, double mono_int, int& nr_occurrences, double& max_ratio, const RangeMobility& im_range) const;

    /**
      @brief Compare an experimental isotope pattern to a theoretical one

      This function will take an array of isotope intensities @p isotopes_int and compare them
      (by order only; no m/z matching) to the theoretically expected ones for the given @p product_mz using an averagine
      model. The returned value is a Pearson correlation between the
      experimental and theoretical pattern.
    */
    double scoreIsotopePattern_(const std::vector<double>& isotopes_int,
                                double product_mz,
                                int putative_fragment_charge) const;

    /**
    @brief Compare an experimental isotope pattern to a theoretical one

    This function will take an array of isotope intensities and compare them
    (by order only; no m/z matching) to the theoretically expected ones for the given @p sum_formula.
    The returned value is a Pearson correlation between the experimental and theoretical pattern.
    */
    double scoreIsotopePattern_(const std::vector<double>& isotopes_int,
                                const EmpiricalFormula& sum_formula) const;

    /**
    @brief Compare an experimental isotope pattern to a theoretical one

    This function will take an array of isotope intensities and compare them
    (by order only; no m/z matching) to the theoretically expected ones given by @p isotope_dist.
    The returned value is a Pearson correlation between the experimental and theoretical pattern.
    */
    double scoreIsotopePattern_(const std::vector<double>& isotopes_int,
                                const IsotopeDistribution& isotope_dist) const;

    /// Get the intensities of isotopes around @p precursor_mz in experimental @p spectrum
    /// and fill @p isotopes_int.
    void getIsotopeIntysFromExpSpec_(double precursor_mz, const SpectrumSequence& spectrum, int charge_state, const RangeMobility& im_range,
                                     std::vector<double>& isotopes_int) const;

    // Parameters
    double dia_extract_window_;
    double dia_byseries_intensity_min_;
    double dia_byseries_ppm_diff_;
    double dia_nr_isotopes_;
    double dia_nr_charges_;
    double peak_before_mono_max_ppm_diff_;
    bool dia_extraction_ppm_;
    bool dia_centroided_;

    TheoreticalSpectrumGenerator * generator;
  };
}
