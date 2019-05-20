// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Hannes Roest$
// --------------------------------------------------------------------------

#pragma once

#include <string>
#include <boost/math/special_functions/fpclassify.hpp> // for isnan
#include <boost/numeric/conversion/cast.hpp>

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>

#include "OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h"
#include "OpenMS/OPENSWATHALGO/DATAACCESS/ITransition.h"
#include "OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h"
#include "OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h"
#include "OpenMS/OPENSWATHALGO/ALGO/Scoring.h"
//#include "OpenMS/OPENSWATHALGO/ALGO/DIAHelpers.h"

namespace OpenSwath
{

  /**
      @brief This class implements different scores for peaks found in SRM/MRM.

      It uses scores based on different parameters of the peaks from the
      individual transitions and stores them individually. The idea and the
      scores are based on the following paper:
        Reiter L, Rinner O, Picotti P, Huettenhain R, Beck M, Brusniak MY,
        Hengartner MO, Aebersold R.  mProphet: automated data processing and
        statistical validation for large-scale SRM experiments.  Nat Methods.
        2011 May;8(5):430-5. Epub 2011 Mar 20.

      The currently implemented scores include:
      - xcorr_coelution: Cross-correlation of the different transitions
      - xcorr_shape: Cross-correlation shape score (whether the maximal
                     Cross-correlation coincides with the maximal intensity)
      - library_rmsd: normalized RMSD of the measured intensities to the expected intensities
      - library_correlation: correlation of the measured intensities to the expected intensities
      - rt_score: deviation from the expected retention time
      - elution_fit_score: how well the elution profile fits a theoretical elution profile

  */
  class OPENSWATHALGO_DLLAPI MRMScoring
  {

public:

    ///Type definitions
    //@{
    /// Cross Correlation array
    typedef OpenSwath::Scoring::XCorrArrayType XCorrArrayType;
    /// Cross Correlation matrix
    typedef std::vector<std::vector<XCorrArrayType> > XCorrMatrixType;

    typedef std::string String;

    typedef OpenSwath::SpectrumPtr SpectrumType;
    typedef OpenSwath::LightTransition TransitionType;
    typedef OpenSwath::LightCompound PeptideType;
    typedef OpenSwath::LightProtein ProteinType;

    typedef boost::shared_ptr<OpenSwath::IFeature> FeatureType;
    //@}

    /** @name Accessors */
    //@{
    /// non-mutable access to the cross-correlation matrix
    const XCorrMatrixType& getXCorrMatrix() const;
    //@}

    /// non-mutable access to the cross-correlation contrast matrix
    const XCorrMatrixType& getXCorrContrastMatrix() const;
    //@}

    /// non-mutable access to the cross-correlation precursor contrast matrix
    const XCorrMatrixType& getXCorrPrecursorContrastMatrix() const;
    //@}

    /// non-mutable access to the cross-correlation precursor contrast matrix
    const XCorrMatrixType& getXCorrPrecursorCombinedMatrix() const;
    //@}

    /** @name Scores */
    //@{
    /// Initialize the scoring object and building the cross-correlation matrix
    void initializeXCorrMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& native_ids);

    /// Initialize the scoring object and building the cross-correlation matrix of chromatograms of set1 (e.g. identification transitions) vs set2 (e.g. detection transitions)
    void initializeXCorrContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& native_ids_set1, const std::vector<String>& native_ids_set2);

    /// Initialize the scoring object and building the cross-correlation matrix
    void initializeXCorrPrecursorMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& precursor_ids);

    /// Initialize the scoring object and building the cross-correlation matrix of chromatograms of precursor isotopes vs transitions
    void initializeXCorrPrecursorContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& precursor_ids, const std::vector<String>& native_ids);

    /// Initialize the scoring object and building the cross-correlation matrix of chromatograms of precursor isotopes and transitions
    void initializeXCorrPrecursorCombinedMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& precursor_ids, const std::vector<String>& native_ids);

    /// calculate the cross-correlation score
    double calcXcorrCoelutionScore();

    /// calculate the weighted cross-correlation score
    double calcXcorrCoelutionWeightedScore(const std::vector<double>& normalized_library_intensity);

    /// calculate the cross-correlation contrast score
    double calcXcorrContrastCoelutionScore();

    /// calculate the separate cross-correlation contrast score
    std::vector<double> calcSeparateXcorrContrastCoelutionScore();

    /// calculate the precursor cross-correlation contrast score
    double calcXcorrPrecursorCoelutionScore();

    /// calculate the precursor cross-correlation contrast score against the transitions
    double calcXcorrPrecursorContrastCoelutionScore();

    /// calculate the precursor cross-correlation coelution score including the transitions
    double calcXcorrPrecursorCombinedCoelutionScore();

    /// calculate the cross-correlation shape score
    double calcXcorrShapeScore();

    /// calculate the weighted cross-correlation shape score
    double calcXcorrShapeWeightedScore(const std::vector<double>& normalized_library_intensity);

    /// calculate the cross-correlation contrast shape score
    double calcXcorrContrastShapeScore();

    /// calculate the separate cross-correlation contrast shape score
    std::vector<double> calcSeparateXcorrContrastShapeScore();

    /// calculate the precursor cross-correlation shape score
    double calcXcorrPrecursorShapeScore();

    /// calculate the precursor cross-correlation shape score against the transitions
    double calcXcorrPrecursorContrastShapeScore();

    /// calculate the precursor cross-correlation shape score including the transitions
    double calcXcorrPrecursorCombinedShapeScore();

    /// calculate the library correlation score
    static void calcLibraryScore(OpenSwath::IMRMFeature* mrmfeature,
      const std::vector<TransitionType>& transitions, double& correlation, 
      double& norm_manhattan, double& manhattan, double& dotprod, 
      double& spectral_angle, double& rmsd);

    /// calculate the retention time correlation score
    static double calcRTScore(const PeptideType& peptide, double normalized_experimental_rt);

    /// calculate the Signal to Noise ratio
    //  using a vector of SignalToNoiseEstimatorMedian that were calculated for
    //  each chromatogram of the transition_group.
    static double calcSNScore(OpenSwath::IMRMFeature* mrmfeature, 
        std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators);

    static std::vector<double> calcSeparateSNScore(OpenSwath::IMRMFeature* mrmfeature, 
        std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators);

    /// non-mutable access to the MI matrix
    const std::vector< std::vector<double> > & getMIMatrix() const;
    //@}

    /// non-mutable access to the MI contrast matrix
    const std::vector< std::vector<double> > & getMIContrastMatrix() const;
    //@}

    /// non-mutable access to the MI precursor contrast matrix
    const std::vector< std::vector<double> > & getMIPrecursorContrastMatrix() const;
    //@}

    /// non-mutable access to the MI precursor combined matrix
    const std::vector< std::vector<double> > & getMIPrecursorCombinedMatrix() const;
    //@}

    /// Initialize the scoring object and building the MI matrix
    void initializeMIMatrix(OpenSwath::IMRMFeature* mrmfeature, std::vector<String> native_ids);

    /// Initialize the scoring object and building the MI matrix of chromatograms of set1 (e.g. identification transitions) vs set2 (e.g. detection transitions)
    void initializeMIContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, std::vector<String> native_ids_set1, std::vector<String> native_ids_set2);

    /// Initialize the scoring object and building the MI matrix
    void initializeMIPrecursorMatrix(OpenSwath::IMRMFeature* mrmfeature, std::vector<String> precursor_ids);

    /// Initialize the mutual information vector against the MS1 trace
    void initializeMIPrecursorContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& precursor_ids, const std::vector<String>& native_ids);

    /// Initialize the mutual information vector with the MS1 trace
    void initializeMIPrecursorCombinedMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<String>& precursor_ids, const std::vector<String>& native_ids);

    double calcMIScore();
    double calcMIWeightedScore(const std::vector<double>& normalized_library_intensity);
    double calcMIPrecursorScore();
    double calcMIPrecursorContrastScore();
    double calcMIPrecursorCombinedScore();
    std::vector<double> calcSeparateMIContrastScore();

    //@}

private:

    /** @name Members */
    //@{
    /// the precomputed cross correlation matrix
    XCorrMatrixType xcorr_matrix_;

    /// the precomputed contrast cross correlation
    XCorrMatrixType xcorr_contrast_matrix_;
    //@}

    /// the precomputed cross correlation matrix of the MS1 trace
    XCorrMatrixType xcorr_precursor_matrix_;

    /// the precomputed cross correlation against the MS1 trace
    XCorrMatrixType xcorr_precursor_contrast_matrix_;
    //@}

    /// the precomputed cross correlation with the MS1 trace
    XCorrMatrixType xcorr_precursor_combined_matrix_;
    //@}

    /// the precomputed mutual information matrix
    std::vector< std::vector<double> > mi_matrix_;

    /// the precomputed contrast mutual information matrix
    std::vector< std::vector<double> > mi_contrast_matrix_;

    /// the precomputed mutual information matrix of the MS1 trace
    std::vector< std::vector<double> > mi_precursor_matrix_;

    /// the precomputed contrast mutual information matrix against the MS1 trace
    std::vector< std::vector<double> > mi_precursor_contrast_matrix_;
    //@}

    /// the precomputed contrast mutual information matrix with the MS1 trace
    std::vector< std::vector<double> > mi_precursor_combined_matrix_;
    //@}

  };
}

