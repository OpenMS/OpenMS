// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_PICKEDPEAKSEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_PICKEDPEAKSEEDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSweepSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>

#include <gsl/gsl_cdf.h>

namespace OpenMS
{
	/** 
		@brief Seeding module which tries to find seeds by looking at 

		This extender module sweeps through the scans and classifies cluster 
 		of peaks as candidate peptides if the distance between successive peaks 
 		is 1 Da (charge 1) , 0.5 Da (charge 2) or 0.3 Da (charge 3) etc.
		 
		@ref PickedPeakSeeder_Parameters are explained on a separate page.
		  	
		@ingroup FeatureFinder
	*/
	class PickedPeakSeeder
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
		
	    /// Default constructor
	    PickedPeakSeeder();
	
	    /// destructor
	    virtual ~PickedPeakSeeder();

	    /// Copy constructor
	    PickedPeakSeeder(const PickedPeakSeeder& rhs);
	    
	    /// Assignment operator
	    PickedPeakSeeder& operator= (const PickedPeakSeeder& rhs);
	
			/// Creates an instance of this class
	    static BaseSeeder* create()
	    {
	        return new PickedPeakSeeder();
	    }
	
			/// Well....
	    static const String getProductName()
	    {
	        return "PickedPeakSeeder";
	    }
		
		protected:
		
			/// keeps member and param entries in synchrony
			virtual void updateMembers_();
	
			/// detects isotopic pattern in a scan
			ScoredMZVector detectIsotopicPattern_(SpectrumType& scan );	

	    /// Find out which charge state belongs to this distance
	    UInt distanceToCharge_(CoordinateType dist);
			
			/// Scores a group of peaks
			ProbabilityType scorePattern_(std::vector<IntensityType>& data, std::vector<IntensityType>& model);
		  		
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
	    /// lower bound for distance between charge 5 peaks
	    CoordinateType charge5_lb_;
			
			/// Averagine model used to determine the significance of a isotopic pattern
			IsotopeModel isomodel_;
			
			/// Minimum number of peaks per pattern
			UInt min_peaks_;
	};
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_PICKEDPEAKSEEDER_H
