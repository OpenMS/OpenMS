// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $a
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <string>

#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ITransition.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h>
#include <OpenMS/OPENSWATHALGO/ALGO/Scoring.h>

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
    class OPENMS_DLLAPI MRMScoring
    {

    public:
        ///Type definitions
        //@{
        /// Cross Correlation array
        typedef OpenSwath::Scoring::XCorrArrayType XCorrArrayType;
        /// Cross Correlation matrix
        typedef OpenMS::Matrix<XCorrArrayType> XCorrMatrixType;

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
        void initializeXCorrMatrix(const std::vector< std::vector< double > >& data);

        /// Initialize the scoring object and building the cross-correlation matrix
        void initializeXCorrMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& native_ids);

        /// Initialize the scoring object and building the cross-correlation matrix of chromatograms of set1 (e.g. identification transitions) vs set2 (e.g. detection transitions)
        void initializeXCorrContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& native_ids_set1, const std::vector<std::string>& native_ids_set2);

        /// Initialize the scoring object and building the cross-correlation matrix
        void initializeXCorrPrecursorMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids);

        /// Initialize the scoring object and building the cross-correlation matrix of chromatograms of precursor isotopes vs transitions
        void initializeXCorrPrecursorContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids, const std::vector<std::string>& native_ids);

        /// Initialize the scoring object and building the cross-correlation matrix of chromatograms of precursor isotopes vs transitions
        void initializeXCorrPrecursorContrastMatrix(const std::vector< std::vector< double > >& data_precursor, const std::vector< std::vector< double > >& data_fragments);

        /// Initialize the scoring object and building the cross-correlation matrix of chromatograms of precursor isotopes and transitions
        void initializeXCorrPrecursorCombinedMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids, const std::vector<std::string>& native_ids);

        /**
           @brief Calculate the cross-correlation coelution score

           The score is a distance where zero indicates perfect coelution.
        */
        double calcXcorrCoelutionScore();

        /**
           @brief Calculate the weighted cross-correlation coelution score

           The score is a distance where zero indicates perfect coelution. The
           score is weighted by the transition intensities, non-perfect coelution
           in low-intensity transitions should thus become less important.
        */
        double calcXcorrCoelutionWeightedScore(const std::vector<double>& normalized_library_intensity);

        /// calculate the separate cross-correlation contrast score
        std::vector<double> calcSeparateXcorrContrastCoelutionScore();

        /// calculate the precursor cross-correlation contrast score
        double calcXcorrPrecursorCoelutionScore();

        /**
           @brief Calculate the precursor cross-correlation contrast score against the transitions

           The score is a distance where zero indicates perfect coelution.
        */
        double calcXcorrPrecursorContrastCoelutionScore();

        /**
           @brief Calculate the precursor cross-correlation contrast score against the sum of transitions
	   implemented the same as calcXcorrPrecursorCoelutionScore(), however assertion is different.

           The score is a distance where zero indicates perfect coelution.
        */
        double calcXcorrPrecursorContrastSumFragCoelutionScore();

        /// calculate the precursor cross-correlation coelution score including the transitions
        double calcXcorrPrecursorCombinedCoelutionScore();

        /**
           @brief Calculate the cross-correlation shape score

           The score is a correlation measure where 1 indicates perfect correlation
           and 0 means no correlation.
        */
        double calcXcorrShapeScore();

        /**
           @brief Calculate the weighted cross-correlation shape score

           The score is a correlation measure where 1 indicates perfect correlation
           and 0 means no correlation. The score is weighted by the transition
           intensities, non-perfect coelution in low-intensity transitions should
           thus become less important.
        */
        double calcXcorrShapeWeightedScore(const std::vector<double>& normalized_library_intensity);

        /// calculate the cross-correlation contrast shape score
        double calcXcorrContrastShapeScore();

        /// calculate the separate cross-correlation contrast shape score
        std::vector<double> calcSeparateXcorrContrastShapeScore();

        /// calculate the precursor cross-correlation shape score
        double calcXcorrPrecursorShapeScore();

        /// calculate the precursor cross-correlation shape score against the transitions
        double calcXcorrPrecursorContrastShapeScore();

        /**
           @brief Calculate the precursor cross-correlation contrast score against the sum of transitions
	   implemented the same as calcXcorrPrecursorContrastShapeScore(), however assertion is different.

           The score is a distance where zero indicates perfect coelution.
        */
        double calcXcorrPrecursorContrastSumFragShapeScore();

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
        const OpenMS::Matrix<double> & getMIMatrix() const;
        //@}

        /// non-mutable access to the MI contrast matrix
        const OpenMS::Matrix<double> & getMIContrastMatrix() const;
        //@}

        /// non-mutable access to the MI precursor contrast matrix
        const OpenMS::Matrix<double> & getMIPrecursorContrastMatrix() const;
        //@}

        /// non-mutable access to the MI precursor combined matrix
        const OpenMS::Matrix<double> & getMIPrecursorCombinedMatrix() const;
        //@}

        /// Initialize the scoring object and building the MI matrix
        void initializeMIMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& native_ids);

        /// Initialize the scoring object and building the MI matrix of chromatograms of set1 (e.g. identification transitions) vs set2 (e.g. detection transitions)
        void initializeMIContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& native_ids_set1, const std::vector<std::string>& native_ids_set2);

        /// Initialize the scoring object and building the MI matrix
        void initializeMIPrecursorMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids);

        /// Initialize the mutual information vector against the MS1 trace
        void initializeMIPrecursorContrastMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids, const std::vector<std::string>& native_ids);

        /// Initialize the mutual information vector with the MS1 trace
        void initializeMIPrecursorCombinedMatrix(OpenSwath::IMRMFeature* mrmfeature, const std::vector<std::string>& precursor_ids, const std::vector<std::string>& native_ids);

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

        /// contains max Peaks from xcorr_matrix_
        OpenMS::Matrix<int> xcorr_matrix_max_peak_;
        OpenMS::Matrix<double> xcorr_matrix_max_peak_sec_;

        /// the precomputed contrast cross correlation
        XCorrMatrixType xcorr_contrast_matrix_;
        //@}

        /// contains max Peaks from xcorr_contrast_matrix_
        OpenMS::Matrix<double> xcorr_contrast_matrix_max_peak_sec_;

        /// the precomputed cross correlation matrix of the MS1 trace
        XCorrMatrixType xcorr_precursor_matrix_;

        /// the precomputed cross correlation against the MS1 trace
        XCorrMatrixType xcorr_precursor_contrast_matrix_;
        //@}

        /// the precomputed cross correlation with the MS1 trace
        XCorrMatrixType xcorr_precursor_combined_matrix_;
        //@}

        /// the precomputed mutual information matrix

        OpenMS::Matrix<double> mi_matrix_;
        /// the precomputed contrast mutual information matrix
        OpenMS::Matrix<double> mi_contrast_matrix_;

        /// the precomputed mutual information matrix of the MS1 trace
        OpenMS::Matrix<double> mi_precursor_matrix_;

        /// the precomputed contrast mutual information matrix against the MS1 trace
        OpenMS::Matrix<double> mi_precursor_contrast_matrix_;
        //@}

        /// the precomputed contrast mutual information matrix with the MS1 trace
        OpenMS::Matrix<double> mi_precursor_combined_matrix_;
        //@}
    };
}
