// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_DIASCORING_H
#define OPENMS_ANALYSIS_OPENSWATH_DIASCORING_H

#include <boost/math/special_functions/fpclassify.hpp> // for isnan
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ITransition.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

namespace OpenMS
{
  class TheoreticalSpectrumGenerator;

  /**
    @brief Scoring of an spectrum at the peak apex of an chromatographic elution peak.

    In DIA (data independent acquisition) / SWATH analysis, at each
    chromatographic point a full MS2 spectrum is recorded. This class allows to
    compute a number of scores based on the full MS2 spectrum available. The scores are the following:

    - isotope scores:
      -- isotope_corr: computes the correlation of each fragment ion with the
         theoretical isotope distribution. This is the pearson correlation to
         the theoretical isotope pattern weighted by the relative intensity of
         the transition (more is better).
      -- isotope_overlap: checks whether a signal at position (mz - 1) / charge
         exists and how strong it is. This would be an indication that the current
         peak is an isotopic signal of another peak. This simply counts how
         often a peak was observed that is higher than the current peak, thus
         number is then weighted by the relative intensity of the transition
         (thus less is better here).

    - massdiff score: computes the difference in ppm of the experimental signal
         to the expected signal (thus less is better)

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
                            SpectrumPtrType spectrum, OpenSwath::IMRMFeature* mrmfeature, double& isotope_corr,
                            double& isotope_overlap);

    /// Massdiff scores, see class description
    void dia_massdiff_score(const std::vector<TransitionType>& transitions,
                            SpectrumPtrType spectrum, const std::vector<double>& normalized_library_intensity,
                            double& ppm_score, double& ppm_score_weighted);

    /**
      Precursor massdifference score

      @param precursor_mz Exact m/z of the precursor to be evaluated
      @param spectrum MS1 spectrum to be evaluated
      @param ppm_score Resulting score
      @return False if no signal was found (and no sensible score calculated), true otherwise
    */
    bool dia_ms1_massdiff_score(double precursor_mz, SpectrumPtrType spectrum,
                                double& ppm_score);

    /// Precursor isotope scores for precursors (peptides and metabolites)
    void dia_ms1_isotope_scores(double precursor_mz, SpectrumPtrType spectrum, size_t charge_state, 
                                double& isotope_corr, double& isotope_overlap, std::string sum_formula = "");

    /// b/y ion scores
    void dia_by_ion_score(SpectrumPtrType spectrum, AASequence& sequence,
                          int charge, double& bseries_score, double& yseries_score);

    /// Dotproduct / Manhatten score with theoretical spectrum
    void score_with_isotopes(SpectrumPtrType spectrum, const std::vector<TransitionType>& transitions,
                             double& dotprod, double& manhattan);
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
                                SpectrumPtrType spectrum, std::map<std::string, double>& intensities,
                                double& isotope_corr, double& isotope_overlap);

    /// retrieves intensities from MRMFeature
    /// computes a vector of relative intensities for each feature (output to intensities)
    void getFirstIsotopeRelativeIntensities_(const std::vector<TransitionType>& transitions,
                                            OpenSwath::IMRMFeature* mrmfeature,
                                            std::map<std::string, double>& intensities //experimental intensities of transitions
                                            );

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
      @param nr_occurrences Will contain the maximum ratio of a peaks intensity compared to the monoisotopic peak intensity how often a peak is found at lower m/z than mono_mz with an intensity higher than mono_int. Multiple charge states are tested, see class parameter dia_nr_charges_

    */
    void largePeaksBeforeFirstIsotope_(SpectrumPtrType spectrum, double mono_mz, double mono_int, int& nr_occurrences, double& max_ratio);

    /**
      @brief Compare an experimental isotope pattern to a theoretical one

      This function will take an array of isotope intensities and compare them
      to the theoretically expected ones for the given m/z using an averagine
      model. The returned value is a Pearson correlation between the
      experimental and theoretical pattern.
    */
    double scoreIsotopePattern_(double product_mz, const std::vector<double>& isotopes_int, 
                                int putative_fragment_charge, std::string sum_formula = "");

    // Parameters
    double dia_extract_window_;
    double dia_centroided_;
    double dia_byseries_intensity_min_;
    double dia_byseries_ppm_diff_;
    double dia_nr_isotopes_;
    double dia_nr_charges_;
    double peak_before_mono_max_ppm_diff_;
    bool dia_extraction_ppm_;

    TheoreticalSpectrumGenerator * generator;
  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_DIASCORING_H

