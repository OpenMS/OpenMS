// -*- C++: make; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MARRWAVELETSEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MARRWAVELETSEEDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSweepSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>


#include <OpenMS/CONCEPT/Exception.h>

#include <gsl/gsl_cdf.h>

namespace OpenMS
{
	/**
		@brief Seeding module based on the Marr wavelet transform to detect (poorly resolved) isotopic pattern.

 		Uses the continuous wavelet transform (and the Marr mother wavelet) to detect isotopic pattern
		in each scan. Patterns that occur in several consecutive scans are joined to seeding regions
		for the extension phase.

		The algorithm considers local maxima in the wavelet transform signal and checks for maxima
		with a distance corresponding to isotopic pattern (e.g. 1 Th, 1/2 Th etc).

		Regions with local maxima a scored based on a F-statistic (compares variance of intervals in cwt).

		@improvement The same thing is called avg_signal_factor in MarrWaveletSeeder and signal_avg_factor in IsotopeWaveletSeeder.  Unify this.
		 
		@ref MarrWaveletSeeder_Parameters are explained on a separate page.

		@ingroup FeatureFinder
	*/
  class MarrWaveletSeeder
    : public BaseSweepSeeder
  {
  	public:

			/// intensity of a peak
			typedef FeaFiTraits::IntensityType IntensityType;
			/// coordinate ( in rt or m/z )
		  typedef FeaFiTraits::CoordinateType CoordinateType;
			/// score
		  typedef DoubleReal ProbabilityType;

			/// a single MS spectrum
			typedef BaseSweepSeeder::SpectrumType SpectrumType;

			/// charge state estimate with associated score
			typedef BaseSweepSeeder::ScoredChargeType ScoredChargeType;
			/// m/z position in spectrum with charge estimate and score
			typedef BaseSweepSeeder::ScoredMZType ScoredMZType;
			/// container of scored m/z positions
			typedef BaseSweepSeeder::ScoredMZVector ScoredMZVector;

		  /// default constructor
	    MarrWaveletSeeder();

	    /// destructor
	    virtual ~MarrWaveletSeeder();

	    /// copy constructor
	    MarrWaveletSeeder(const MarrWaveletSeeder& rhs);

	    /// assignment operator
	    MarrWaveletSeeder& operator= (const MarrWaveletSeeder& rhs);

			/// Creates an instance of this class
	    static BaseSeeder* create()
	    {
	      return new MarrWaveletSeeder();
	    }

			/// Well....
	    static const String getProductName()
	    {
	      return "MarrWaveletSeeder";
	    }

			protected:

			/// keeps member and param entries in synchrony
	  	virtual void updateMembers_();

			/// detects an isotopic pattern in a scan
			ScoredMZVector detectIsotopicPattern_(SpectrumType& scan );

			/// Finds local maxima in cwt
			void getMaxPositions_( const SpectrumType::const_iterator& first,
														 const SpectrumType::const_iterator& last,
														 std::vector<Int>& localmax
														#ifdef DEBUG_FEATUREFINDER
														 ,CoordinateType curr_peak
														#endif
														 );

			/// Compute local variance and test for significance
			ProbabilityType testLocalVariance_(const std::vector<Int>& local_maxima, const UInt max_index);

	    /// estimate charge state
	    UInt distanceToCharge_(CoordinateType dist);

		  /// indicates whether this module has been initialized
		  bool is_initialized_;

			/// upper bound for distance between charge 1 peaks
			CoordinateType charge1_ub_;
			/// lower bound for distance between charge 1 peaks
			CoordinateType charge1_lb_;

			/// upper bound for distance between charge 2 peaks
			CoordinateType charge2_ub_;
			/// lower bound for distance between charge 2 peaks
			CoordinateType charge2_lb_;

			/// upper bound for distance between charge 3 peaks
			CoordinateType charge3_ub_;
			/// lower bound for distance between charge 3 peaks
			CoordinateType charge3_lb_;

			/// upper bound for distance between charge 4 peaks
			CoordinateType charge4_ub_;
			/// lower bound for distance between charge 4 peaks
			CoordinateType charge4_lb_;

			/// upper bound for distance between charge 5 peaks
			CoordinateType charge5_ub_;
			/// lower bound for distance between charge 4 peaks
			CoordinateType charge5_lb_;

			/// computes the wavelet transform for a given scan
			ContinuousWaveletTransformNumIntegration cwt_;

			/// intensity threshold for spectrum
			IntensityType avg_signal_factor_;

	   	/// intensity threshold for cwt
	   	IntensityType avg_cwt_factor_;

			/// Marr wavelet scale
			CoordinateType cwt_scale_;

			/// Minimum number of local maxima in cwt for an isotopic pattern
			UInt min_peaks_;

  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_MARRWAVELETSEEDER_H
