// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/OfflinePrecursorIonSelection.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/ANALYSIS/ID/ILPWrapper.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ILPWrapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ILPWrapper* ptr = 0;
START_SECTION(ILPWrapper())
{
	ptr = new ILPWrapper();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~ILPWrapper())
{
	delete ptr;
}
END_SECTION


START_SECTION((template < typename InputPeakType > void encodeModelForKnownLCMSMapFeatureBased(FeatureMap<> &features, MSExperiment< InputPeakType > &experiment, std::vector< IndexTriple > &variable_indices, std::vector< std::vector< std::pair< Size, Size > > > &mass_ranges, std::set< Int > &charges_set, UInt ms2_spectra_per_rt_bin)))
{
	std::set<Int> charges_set;
	charges_set.insert(1);
	
	FeatureMap<> features;
	MSExperiment<> exp;
	std::vector<ILPWrapper::IndexTriple > variable_indices;
	std::vector<std::vector<std::pair<Size,Size> > > mass_ranges;
	ILPWrapper wrapper;
	FeatureMap<> map;

  std::vector<int> solution_indices;

  // test empty input
	ILPWrapper wrapper2;
  wrapper2.createAndSolveILPForKnownLCMSMapFeatureBased(features,exp,variable_indices,mass_ranges,charges_set,1,solution_indices);
	TEST_EQUAL(variable_indices.size(),0)
	TEST_EQUAL(solution_indices.size(),0)
	solution_indices.clear();

	// now with the same input as with the offline precursor ion selection (can't test them separately)
	FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("OfflinePrecursorIonSelection_features.featureXML"),map);
	MSExperiment<> raw_data;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OfflinePrecursorIonSelection_raw_data.mzML"),raw_data);
	mass_ranges.clear();
	OfflinePrecursorIonSelection ops;
	ops.getMassRanges(map,raw_data,mass_ranges);
	wrapper.createAndSolveILPForKnownLCMSMapFeatureBased(map,raw_data,variable_indices,mass_ranges,charges_set,1,solution_indices);
	TEST_EQUAL(variable_indices.size(),6)
	variable_indices.clear();
	TEST_EQUAL(solution_indices.size(),3)

		
// 	// test if min_peak_distance prohibits nearby peaks	
 // 	ILPWrapper wrapper3;	
// 	wrapper3.createAndSolveILPForKnownLCMSMapFeatureBased(map,raw_data,variable_indices,mass_ranges,charges_set,1,20.,solution_indices);
// 	TEST_EQUAL(variable_indices.size(),2)
//   TEST_EQUAL(solution_indices.size(),2)
	
}
END_SECTION
		
	      
START_SECTION((void solve(std::vector<int>& solution_indices)))
{
}
END_SECTION	      

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
