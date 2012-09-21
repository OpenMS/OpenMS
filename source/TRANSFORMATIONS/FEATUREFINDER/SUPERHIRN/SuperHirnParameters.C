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
// $Maintainer: Peter Kunszt $
// $Authors: Lukas Mueller, Markus Mueller Peter Kunszt $
// --------------------------------------------------------------------------
//

#include <map>
#include <stdio.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>

namespace OpenMS
{

	SuperHirnParameters* SuperHirnParameters::instance_ = NULL;
	bool SuperHirnParameters::haveInstance_ = false;

	SuperHirnParameters::SuperHirnParameters()
	{
		// Background Intensity static vars
		backgroundIntensityBinsTR_ = 2.0;
		backgroundIntensityBinsMZ_ = 50.0;
		backgroundIntensityBinsIntens_ = 50;
		backgroundIntensityBinsMinBinCount_ = 1;

		// these will be overwritten at config time
		minTR_ = 0;
		maxTR_ = 0;
		minFeatureMZ_ = 0;
		maxFeatureMZ_ = 0;
		minFeatureChrg_ = 0;
		maxFeatureChrg_ = 0;


		// minimal intensity level:  NEVER USED
		intensityThreshold_ = 0;

		// m/z tolerance value: NEVER CONFIGURED
		toleranceMZ_ = 10;

		// max_distance from next elution peak member in min.:
		maxInterScanRetentionTimeDistance_ = 0;

		// define minimal number of members in LC elution peaks cluster
		minNbClusterMembers_ = 0;

		/*
		 // to track detected monoistopic mass for debugging:
		 monoIsoDebugging_ = false;
		 debugMonoIsoMassMin_ = 318.00;
		 debugMonoIsoMassMax_ = 319.00;
		 */

		// if data are in centroid form or not:
		centroidDataModus_ = false;

		massTolPpm_ = 10.0;
		massTolDa_ = 0.01;
		minIntensity_ = 0.0;
		intensityFloor_ = 1.0;

		peptideProbabilityThreshold_ = 0.9; // this is hardcoded. it is never configured.
		storeAllLowProbabilityMS2Scans_ = false; // this is hardcoded.


		createFeatureElutionProfiles_ = false;
		/*
		elutionPeakDebugging_ = false;
		elutionPeakMassMin_ = -1;
		elutionPeakMassMax_ = -2;

		// if this option is on, then construct fake features for available MS2 features
		featureFakeInsertionBasedOnMS2Feature_ = true;
		*/

		lowIntensityMSSignalThreshold_ = 1.0; // never configured, but used
		initIsotopeDist_ = false;
	}

}
