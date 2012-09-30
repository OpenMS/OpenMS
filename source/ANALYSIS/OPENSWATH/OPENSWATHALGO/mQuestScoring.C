// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MRMFEATURESCORING_C
#define OPENMS_KERNEL_MRMFEATURESCORING_C

//#define MRMSCORING_TESTING
#include <algorithm>
#include <algorithm>
#include <iterator>
#include <iostream>

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/mQuestScoring.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/meanAndSd.h"

//#include <gsl/gsl_statistics_int.h>
//#include <gsl/gsl_statistics_double.h>

namespace OpenMS
{

	const MRMFeatureScoring::XCorrMatrixType & MRMFeatureScoring::getXCorrMatrix() const
	{
		return xcorr_matrix_;
	}

	// see /IMSB/users/reiterl/bin/code/biognosys/trunk/libs/mrm_libs/MRM_pgroup.pm
	// _calc_xcorr_coelution_score
	//
	//   for each i,j get xcorr_matrix array => find max of the crosscorrelation
	//   store the delta to the retention time
	// return $deltascore_mean + $deltascore_stdev
	double MRMFeatureScoring::calcXcorrCoelutionScore()
	{
		//TODO OPENMS_PRECONDITION(xcorr_matrix_.size() > 1, "Expect cross-correlation matrix of at least 2x2");

		std::vector<int> deltas;
		for (std::size_t i = 0; i < xcorr_matrix_.size(); i++) {
			for (std::size_t  j = i; j < xcorr_matrix_.size(); j++) {
				// first is the X value (RT), should be an int
				deltas.push_back( std::fabs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first));
#ifdef MRMSCORING_TESTING
				std::cout << "&&_xcoel append " << std::fabs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first) << std::endl;
#endif
			}
		}

		OpenSwath::mean_and_stddev msc;
		msc = std::for_each(deltas.begin(),deltas.end(),msc);
		double deltas_mean = msc.mean();
		double deltas_stdv = msc.sample_stddev();
		//double deltas_mean = gsl_stats_int_mean(&deltas[0], 1, deltas.size());
		//double deltas_stdv = gsl_stats_int_sd(&deltas[0], 1, deltas.size());

		double xcorr_coelution_score = deltas_mean + deltas_stdv;
		return xcorr_coelution_score;
	}

	double MRMFeatureScoring::calcXcorrCoelutionScore_weighted(
			const std::vector<double> & normalized_library_intensity)
	{
		//TODO OPENMS_PRECONDITION(xcorr_matrix_.size() > 1, "Expect cross-correlation matrix of at least 2x2");

#ifdef MRMSCORING_TESTING
		double weights = 0;
#endif
		std::vector<double> deltas;
		for (std::size_t i = 0; i < xcorr_matrix_.size(); i++) {
			deltas.push_back(
					std::fabs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][i])->first)
							* normalized_library_intensity[i]
							* normalized_library_intensity[i]);
#ifdef MRMSCORING_TESTING
			std::cout << "_xcoel_weighted " << i << " " << i << " " << Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][i])->first << " weight " <<
			normalized_library_intensity[i] * normalized_library_intensity[i] << std::endl;
			weights += normalized_library_intensity[i] * normalized_library_intensity[i];
#endif
			for (std::size_t j = i + 1; j < xcorr_matrix_.size(); j++) {
				// first is the X value (RT), should be an int
				deltas.push_back(
						std::fabs(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first)
								* normalized_library_intensity[i]
								* normalized_library_intensity[j] * 2);
#ifdef MRMSCORING_TESTING
				std::cout << "_xcoel_weighted " << i << " " << j << " " << Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->first << " weight " <<
				normalized_library_intensity[i] * normalized_library_intensity[j] * 2 << std::endl;
				weights += normalized_library_intensity[i] * normalized_library_intensity[j];
#endif

			}
		}

#ifdef MRMSCORING_TESTING
		std::cout << " all weights sum " << weights << std::endl;
#endif

		/*
		 double deltas_mean = gsl_stats_int_mean(&deltas[0],1,deltas.size() );
		 double deltas_stdv = gsl_stats_int_sd(  &deltas[0],1,deltas.size() );

		 double xcorr_coelution_score = deltas_mean + deltas_stdv;
		 return xcorr_coelution_score;
		 */
		return std::accumulate(deltas.begin(), deltas.end(), 0.0);
	}

	// see /IMSB/users/reiterl/bin/code/biognosys/trunk/libs/mrm_libs/MRM_pgroup.pm
	// _calc_xcorr_shape_score
	//
	//   for each i,j get xcorr_matrix array => find max of the crosscorrelation
	//   calculate whether the maximal crosscorrelation coincides with the maximal intensity
	double MRMFeatureScoring::calcXcorrShape_score()
	{
		//TODO OPENMS_PRECONDITION(xcorr_matrix_.size() > 1, "Expect cross-correlation matrix of at least 2x2");

		std::vector<double> intensities;
		for (std::size_t i = 0; i < xcorr_matrix_.size(); i++) {
			for (std::size_t j = i; j < xcorr_matrix_.size(); j++) {
				// second is the Y value (intensity)
				intensities.push_back(Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->second);
			}
		}
		OpenSwath::mean_and_stddev msc;
		msc = std::for_each(intensities.begin(),intensities.end(),msc);
		return msc.mean();
	}

	double MRMFeatureScoring::calcXcorrShape_score_weighted(
			const std::vector<double> & normalized_library_intensity)
	{
		// TODO OPENMS_PRECONDITION(xcorr_matrix_.size() > 1, "Expect cross-correlation matrix of at least 2x2");

		// TODO : see _calc_weighted_xcorr_shape_score in MRM_pgroup.pm
		//         -- they only multiply up the intensity once
		std::vector<double> intensities;
		for (std::size_t i = 0; i < xcorr_matrix_.size(); i++) {
			intensities.push_back(
					Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][i])->second
							* normalized_library_intensity[i]
							* normalized_library_intensity[i]);
#ifdef MRMSCORING_TESTING
			std::cout << "_xcorr_weighted " << i << " " << i << " " << Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][i])->second << " weight " <<
			normalized_library_intensity[i] * normalized_library_intensity[i] << std::endl;
#endif
			for (std::size_t j = i + 1; j < xcorr_matrix_.size(); j++) {
				intensities.push_back(
						Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->second
								* normalized_library_intensity[i]
								* normalized_library_intensity[j] * 2);
#ifdef MRMSCORING_TESTING
				std::cout << "_xcorr_weighted " << i << " " << j << " " << Scoring::xcorrArrayGetMaxPeak(xcorr_matrix_[i][j])->second << " weight " <<
				normalized_library_intensity[i] * normalized_library_intensity[j] * 2 << std::endl;
#endif
			}
		}
		return std::accumulate(intensities.begin(), intensities.end(), 0.0);
	}


}

#endif
