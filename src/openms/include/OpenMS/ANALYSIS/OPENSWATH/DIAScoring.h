// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

  @ingroup TargetedQuantitation
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
    /**
      @brief Compute the isotope-pattern correlation and isotope-overlap scores for a transition group (see class description).

      Delegates to the internal @ref diaIsotopeScoresSub_ after zero-initialising the two
      output arguments. @p mrmfeature is consulted for per-transition apex intensities used
      to weight the correlation.

      @param[in]  transitions      Transitions whose fragment-ion isotope patterns are scored.
      @param[in]  spectrum         DIA MS2 spectrum sequence (one or more spectra around the apex) the patterns are read from.
      @param[in]  mrmfeature       MRM feature providing per-transition apex intensities used to weight the correlation.
      @param[in]  im_range         Ion-mobility range filter; pass an empty range to disable IM filtering.
      @param[out] isotope_corr     Intensity-weighted Pearson correlation against the averagine isotope distribution (higher = better).
      @param[out] isotope_overlap  Intensity-weighted count of peaks at lower m/z than each transition's monoisotope (lower = better).
    */
    void dia_isotope_scores(const std::vector<TransitionType>& transitions,
                            SpectrumSequence& spectrum,
                            OpenSwath::IMRMFeature* mrmfeature,
                            const RangeMobility& im_range,
                            double& isotope_corr,
                            double& isotope_overlap) const;

    /**
      @brief Compute fragment-ion mass-error scores in ppm against the theoretical product m/z of each transition.

      For every transition, integrates the spectrum in a window around the theoretical
      product m/z (window width and ppm-vs-Th controlled by parameter @c "dia_extraction_window"
      and @c "dia_extraction_unit") and records the signed ppm error between the measured
      and theoretical m/z. Transitions where no signal was found contribute @c -1.0 to
      @p diff_ppm and are excluded from the averages.

      @param[in]  transitions                 Transitions whose theoretical product m/z values are the targets.
      @param[in]  spectrum                    DIA MS2 spectrum sequence.
      @param[in]  normalized_library_intensity Per-transition library intensities (parallel to @p transitions) used to weight @p ppm_score_weighted.
      @param[in]  im_range                    Ion-mobility range filter; pass an empty range to disable IM filtering.
      @param[out] ppm_score                   Mean of the absolute ppm errors over observed transitions (0 if none observed).
      @param[out] ppm_score_weighted          Sum of @c |ppm| * @p normalized_library_intensity over observed transitions.
      @param[out] diff_ppm                    Per-transition signed ppm error; @c -1.0 entries mark transitions with no signal in the extraction window.
    */
    void dia_massdiff_score(const std::vector<TransitionType>& transitions,
                            const SpectrumSequence& spectrum,
                            const std::vector<double>& normalized_library_intensity,
                            const RangeMobility& im_range,
                            double& ppm_score,
                            double& ppm_score_weighted,
                            std::vector<double>& diff_ppm) const;

    /**
      @brief Compute the precursor-mass-error score on an MS1 spectrum.

      Integrates the spectrum in a window around @p precursor_mz and reports the absolute
      ppm difference between the measured and theoretical m/z.

      @param[in]  precursor_mz Exact m/z of the precursor to be evaluated.
      @param[in]  spectrum     MS1 spectrum sequence.
      @param[in]  im_range     Ion-mobility range filter; pass an empty range to disable IM filtering.
      @param[out] ppm_score    On success, the absolute ppm error of the measured peak. On failure (no signal found), set to the absolute ppm width of the extraction window as a worst-case bound — not @c -1.
      @return @c true if a signal was integrated successfully, @c false if no peak was found inside the extraction window (in which case @p ppm_score still carries the worst-case bound described above).
    */
    bool dia_ms1_massdiff_score(double precursor_mz, const SpectrumSequence& spectrum, const RangeMobility& im_range,
                                double& ppm_score) const;

    /**
      @brief Precursor isotope-pattern scores using an averagine model.

      Drives the same scoring as @ref dia_ms1_isotope_scores but builds the theoretical
      isotope pattern via an averagine model parameterised by @p precursor_mz and
      @p charge_state (suitable for peptides where no exact sum formula is available).

      @param[in]  precursor_mz    Exact m/z of the precursor to be evaluated.
      @param[in]  spectrum        MS1 spectrum sequence.
      @param[in]  charge_state    Precursor charge state used to space the averagine isotope envelope.
      @param[in,out] im_range     Ion-mobility range filter; may be narrowed during the call.
      @param[out] isotope_corr    Pearson correlation against the averagine isotope distribution.
      @param[out] isotope_overlap Count of peaks at lower m/z that could explain @p precursor_mz as part of a larger isotope envelope.
    */
    void dia_ms1_isotope_scores_averagine(double precursor_mz, const SpectrumSequence& spectrum, int charge_state, RangeMobility& im_range,
                                          double& isotope_corr, double& isotope_overlap) const;
    /**
      @brief Precursor isotope-pattern scores using an explicit sum formula.

      Builds the theoretical pattern from @p sum_formula (and uses its charge), then places
      it at @p precursor_mz (the actual precursor m/z, which may differ slightly from what
      @p sum_formula would predict — e.g. when adducts are present).

      @param[in]  precursor_mz    Exact m/z of the precursor to be evaluated.
      @param[in]  spectrum        MS1 spectrum sequence.
      @param[in,out] im_range     Ion-mobility range filter; may be narrowed during the call.
      @param[out] isotope_corr    Pearson correlation against the formula-derived isotope distribution.
      @param[out] isotope_overlap Count of peaks at lower m/z that could explain @p precursor_mz as part of a larger isotope envelope.
      @param[in]  sum_formula     Neutral sum formula of the precursor compound; its @c getCharge() is used for the isotope spacing.
    */
    void dia_ms1_isotope_scores(double precursor_mz, const std::vector<SpectrumPtrType>& spectrum, RangeMobility& im_range,
                                double& isotope_corr, double& isotope_overlap, const EmpiricalFormula& sum_formula) const;

    /**
      @brief Count the b- and y-fragment ions of a peptide that have evidence in a DIA spectrum.

      For every b/y ion of @p sequence at the given @p charge, a peak is
      considered a match if its ppm error is below @c "dia_byseries_ppm_diff"
      and its intensity is above @c "dia_byseries_intensity_min". The
      returned scores are the number of matching b- and y-ions.

      @param[in]  spectrum       DIA MS2 spectrum sequence.
      @param[in]  sequence       Peptide sequence.
      @param[in]  charge         Fragment-ion charge state. Must be > 0; passing a non-positive value is a programming error (checked in debug builds).
      @param[in]  im_range       Ion-mobility range filter; pass an empty range to disable IM filtering.
      @param[out] bseries_score  Number of matching b-ions.
      @param[out] yseries_score  Number of matching y-ions.
    */
    void dia_by_ion_score(const SpectrumSequence& spectrum, AASequence& sequence,
                          int charge, const RangeMobility& im_range, double& bseries_score, double& yseries_score) const;

    /**
      @brief Dot-product and Manhattan-distance scores between the spectrum
             and the transition group's theoretical isotope pattern.

      @param[in]  spectrum     DIA MS2 spectrum sequence.
      @param[in]  transitions  Transition group whose theoretical pattern is scored against @p spectrum.
      @param[in]  im_range     Ion-mobility range filter; pass an empty range to disable IM filtering.
      @param[out] dotprod      Dot product of the two intensity vectors (higher = better).
      @param[out] manhattan    Manhattan distance of the two intensity vectors (lower = better).
    */
    void score_with_isotopes(SpectrumSequence& spectrum,
                             const std::vector<TransitionType>& transitions,
                             const RangeMobility& im_range,
                             double& dotprod,
                             double& manhattan) const;

    /**
      @brief Return the DIA extraction window.

      The value is expressed in Thomson or ppm, depending on isDIAExtractionPPM().
    */
    double getDIAExtractionWindow() const
    {
      return dia_extract_window_;
    }

    /**
      @brief Return whether the DIA extraction window is interpreted in ppm.

      Returns false when the extraction window is interpreted in Thomson.
    */
    bool isDIAExtractionPPM() const
    {
      return dia_extraction_ppm_;
    }
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
                              OpenSwath::IMRMFeature* mrmfeature,
                              const RangeMobility& im_range,
                              double& isotope_corr,
                              double& isotope_overlap) const;

private:

    /**
      @brief Determine whether the current m/z value is a monoisotopic peak

      This function will try to determine whether the current peak is a
      monoisotopic peak or not. It will do so by searching for an intense peak
      at a lower m/z that could explain the current peak as part of a isotope
      pattern.

      @param[in] spectrum The spectrum (MS1 or MS2)
      @param[in] mono_mz The m/z value where a monoisotopic is expected
      @param[in] mono_int The intensity of the monoisotopic peak (peak at mono_mz)
      @param[out] nr_occurrences Will contain the count of how often a peak is found at lower m/z than mono_mz with an intensity higher than mono_int. Multiple charge states are tested, see class parameter dia_nr_charges_
      @param[out] max_ratio Will contain the maximum ratio of a peaks intensity compared to the monoisotopic peak intensity how often a peak is found at lower m/z than mono_mz with an intensity higher than mono_int. Multiple charge states are tested, see class parameter dia_nr_charges_
      @param[in] im_range Ion mobility subrange to consider (used as filter); can be empty (i.e. no IM filtering)
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
