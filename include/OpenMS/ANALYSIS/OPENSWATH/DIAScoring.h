// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_DIASCORING_H_
#define OPENMS_ANALYSIS_OPENSWATH_DIASCORING_H_

#include <boost/math/special_functions/fpclassify.hpp> // for isnan
#include <OpenMS/CHEMISTRY/AASequence.h>

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ITransition.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h"

namespace OpenSwath
{
  using namespace OpenMS;
  /**
    @brief Scoring of an spectrum at the peak apex of an chromatographic elution peak.

    In DIA (data independent acquisition) / SWATH analysis, at each
    chromatographic point a full MS2 spectrum is recorded. This class allows to
    compute a number of scores based on the full MS2 spectrum available. The scores are the following:

    - isotope scores:
      -- isotope_corr: computes the correlation of each fragment ion with the
         theoretical isotope distribution.
      -- isotope_overlap: checks whether a signal at position (mz - 1) / charge
         exists and how strong it is. This would be an indication that the current
         peak is an isotopic signal of another peak.

    - massdiff score: computes the difference in ppm of the experimental signal to the exepcted signal

    - b/y ion score: checks for the presence of b/y ions of the peptide in question

    - theoretical spectrum: a dotproduct and a manhattan score with a theoretical spectrum

    This class expects spectra objects that implement the OpenSWATH Spectrum
    interface. Transitions are expected to be in the light transition format
    (defined in OPENSWATHALGO/DATAACCESS/TransitionExperiment.h).

  */
  class OPENMS_DLLAPI DIAScoring
  {
    ///Type definitions
    //@{
    /// Spectrum type, see Spectrum interface
    typedef OpenSwath::SpectrumPtr SpectrumType;
    /// Transition interface (Transition, Peptide, Protein)
    typedef OpenSwath::LightTransition TransitionType;
    typedef OpenSwath::LightPeptide PeptideType;
    typedef OpenSwath::LightProtein ProteinType;
    //@}

public:

    ///@name Constructors and Destructor
    //@{
    /// Default constructor
    DIAScoring() {}

    /// Destructor
    virtual ~DIAScoring() {}
    //@}

    ///@name Accessors
    //@{
    /// set parameters for the algorithm
    void set_dia_parameters(double dia_extract_window, double dia_centroided,
      double dia_byseries_intensity_min, double dia_byseries_ppm_diff, double dia_nr_isotopes, double dia_nr_charges)
    {
      dia_extract_window_ = dia_extract_window;
      dia_centroided_ = dia_centroided;
      dia_byseries_intensity_min_ = dia_byseries_intensity_min;
      dia_byseries_ppm_diff_ = dia_byseries_ppm_diff;
      dia_nr_isotopes_ = dia_nr_isotopes;
      dia_nr_charges_ = dia_nr_charges;
    }
    //@}

    ///////////////////////////////////////////////////////////////////////////
    // DIA / SWATH scoring

    ///@name DIA Scores
    //@{
    /// Isotope scores, see class description
    void dia_isotope_scores(const std::vector<TransitionType> & transitions,
      SpectrumType spectrum, OpenSwath::IMRMFeature * mrmfeature, double & isotope_corr,
      double & isotope_overlap);

    /// Massdiff scores, see class description
    void dia_massdiff_score(const std::vector<TransitionType> & transitions,
      SpectrumType spectrum, const std::vector<double> & normalized_library_intensity,
      double & ppm_score, double & ppm_score_weighted);

    /// b/y ion scores
    void dia_by_ion_score(SpectrumType spectrum, AASequence & sequence,
      int charge, double & bseries_score, double & yseries_score);

    /// Dotproduct / Manhatten score with theoretical spectrum
    void score_with_isotopes(SpectrumType spectrum, const std::vector<TransitionType> & transitions,
      double & dotprod, double & manhattan);
    //@}

private:

    /// Copy constructor (algorithm class)
    DIAScoring(const DIAScoring & rhs);

    /// Assignment operator (algorithm class)
    DIAScoring & operator=(const DIAScoring & rhs);

    /// Subfunction of dia_isotop_scores
    void dia_isotope_scores_sub(const std::vector<TransitionType> & transitions,
      SpectrumType spectrum, std::map<std::string, double> & intensities,
      double & isotope_corr, double & isotope_overlap);

    /// retrieves intensities from MRMFeature
    /// computes a vector of relative intensities for each feature (output to intensities) 
    void getFirstIsotopeRelativeIntensities(const std::vector<TransitionType> & transitions,
      OpenSwath::IMRMFeature * mrmfeature,
      std::map<std::string, double> & intensities     //experimental intensities of transitions
      );

#if 0
    /// TODO (wolski): what is this doing here? is this code dead?
    void dia_isotope_scores(const std::vector<TransitionType> & transitions,
      SpectrumType spectrum, int putative_fragment_charge,
      double & isotope_corr, double & isotope_overlap);

    /* computes apex mz and area for a given spectrum fragment */
    /// TODO (wolski): what is this doing here? is this code dead?
    void getSpectrumIntensities(const std::vector<TransitionType> & transitions,
      SpectrumType spectrum, double extractWindow,
      std::vector<double> & mzv, std::vector<double> & intensityv);


    // TODO (wolski) what does this method do, is it called somewhere? 
    void getFirstIsotopeRelativeIntensities(
      const std::vector<TransitionType> & transitions,
      SpectrumType spectrum, std::map<std::string, double> & intensities     //experimental intensities of transitions
      );

    // TODO (wolski) this method is dead? where is the implementation?
    void getFirstIsotopeRelativeIntensities(
      const std::vector<TransitionType> & transitions,
      SpectrumType spectrum, std::vector<double> & intensities     //experimental intensities of transitions
      );
#endif

private:

    /**
      @brief Integrate intensity in a spectrum from start to end 

      This function will integrate the intensity in a spectrum between mz_start
      and mz_end, returning the total intensity and an intensity-weighted m/z
      value.

      @note If there is no signal, mz will be set to -1 and intensity to 0
    */
    void getIntensePeakInWindow(const SpectrumType spectrum, double mz_start,
      double mz_end, double & mz, double & intensity, bool centroided);

    /**
      @brief Search for a large peak _before_ (lower m/z) the current peak 

      This function will try to determine whether the current peak is part of
      an isotopic pattern that does NOT have the current peak as monoisotopic
      peak.
    */
    DoubleReal largePeaksBeforeFirstIsotope(double product_mz,
      SpectrumType & spectrum, double max_ppm_diff, double main_peak);

    /**
      @brief Compare an experimental isotope pattern to a theoretical one 

      This function will take an array of isotope intensities and compare them
      to the theoritcally expected ones using pearson correlation.
    */
    DoubleReal scoreIsotopePattern(double product_mz,
      const std::vector<double> & isotopes_int, int putative_fragment_charge);

    // Parameters
    double dia_extract_window_;
    double dia_centroided_;
    double dia_byseries_intensity_min_;
    double dia_byseries_ppm_diff_;
    double dia_nr_isotopes_;
    double dia_nr_charges_;

  };
}

#endif
