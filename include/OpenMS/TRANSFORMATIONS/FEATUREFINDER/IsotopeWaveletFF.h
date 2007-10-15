// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Rene Hussong$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELET_FF_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_ISOTOPEWAVELET_FF_H

#ifndef NULL
#define NULL 0
#endif

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletTransform.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CHEMISTRY/Averagine.h>


namespace OpenMS
{

template <typename MapType>

/** @brief Implements the isotope wavelet feature finder.
 *
 * 	The IsotopeWaveletFF class has been designed for finding features in 1D or 2D MS data sets using the isotope wavelet.
 * 	In the case of two dimensional data, the class provides additionally the sweep line algorithm. Please not that in
 * 	its current implementation the istope wavelet feature finder is only applicable to raw data (not to picked data). 
 *
 * 	Before you start the algorithm by calling runFF, you have to set up the class by a call to initializeFF. 
 * 	Please note that this class features a singleton implementation, i.e. yo cannot directly intantiate this
 * 	class by a call to its default constructor. */
class IsotopeWaveletFF
{

	public:

		/** @brief Destructor. */		
		virtual ~IsotopeWaveletFF() throw ();

		/** @brief Initializes the isotope wavelet feature finder. 
 			* @param experiment Your MS data.		
 			* @param max_charge Determines the maximal charge state considered by the algorithm (e.g. for MALDI data this will be 1). 
 			* @param threshold The parameter of the isotope wavelet transform. Although this parameter essentially behaves like 
 			* a normal cut off parameter, its meaning is a bit more complicated. Usually you have to play around with this parameter
 			* to determine a reasonable cut off between signal and noise. For peptide mass fingerprints of single scans, it is useful 
 			* to start with @p threshold=0 and to increase the parameter if too many false positives should be found. 
 			* Usually, the obtained result for @p threshold=0 is very significant in terms of MASCOT scores and coverage values. 
 			* For every other type of data, a good starting point might be: @p threshold=10..20. 
 			* @param RT_votes_cutoff The minimum number of scans a pattern should cover in order to be considered as significant. If a pattern
 			* occurs in less than @p RT_votes_cutoff scans, it will be dropped during the sweep line phase. 
 			* @param RT_interleave The number of scans a pattern is allowed to "skip". E.g. if a pattern spreads over 10 patterns, 
 			* but its signal is disrupted, s.t. the peptide occurs in 4 subsequent and later on again in 6 subsequent scans, 
 			* it will be dropped if the interleft section is larger than @p RT_interleave scans. Usually you do not have to adapt this
 			* parameter. 
 			* @param hash_precision. An internally used parameter. You do not have to set this parameter. In the case of 
 			* really highly resoluted spectra, you might increase this parameter to 4 or 5 to get better predictions w.r.t. the 
 			* m/z axis. */
		static IsotopeWaveletFF* initializeFF (MapType& experiment, const unsigned int max_charge, const double threshold, 
			const unsigned int RT_votes_cutoff, const unsigned int RT_interleave=2, const unsigned int hash_precision=3)
		{
		  if (me_ == NULL)
			{
				me_ = new IsotopeWaveletFF (experiment, max_charge, threshold, RT_votes_cutoff, RT_interleave, hash_precision);
				IsotopeWavelet::preComputeGammaFunction();
			}

		  return (me_);
		};

		/** @brief Runs the algorithm on the provided data set beginning at scan @p start_scan and ending (exclusively) at @p end_scan. */
		virtual FeatureMap<Feature> runFF (const int start_scan=0, const int end_scan=-1) throw ();
	

	protected:

		/** @brief Internally used data struture for the sweep line algorithm. */
		struct BoxElement
			{			
				double mz;
				unsigned int c; //Note, this is not the charge (it is charge-1!!!)
				double score;
				double intens;
				double RT; //The elution time (not the scan index)
			};				

		typedef std::map<unsigned int, BoxElement> Box; //Key: RT (index), value: BoxElement
		typedef DRawDataPoint<2> RawDataPoint2D; 

	
		/** @brief Default Constructor */
		IsotopeWaveletFF () throw ();
		
		/** @brief Constructor */
		IsotopeWaveletFF (MapType& experiment, const unsigned int max_charge, const double threshold, 
			const unsigned int RT_votes_cutoff, const unsigned int RT_interleave, const unsigned int hash_precision) throw ();

		/** Singleton pointer */
		static IsotopeWaveletFF* me_;		

		/** The experimentally recorded MS map */
		MapType& experiment_;

		unsigned int max_charge_;
		double threshold_;
		unsigned int RT_votes_cutoff_;
		unsigned int RT_interleave_;
		double hash_precision_;
};

#include "../../../../source/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletFF.C"

} //namespace

#endif 
