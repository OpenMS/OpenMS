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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TARGETED/OfflinePrecursorIonSelection.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(OfflinePrecursorIonSelection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OfflinePrecursorIonSelection* ptr = 0;
OfflinePrecursorIonSelection* nullPointer = 0;
START_SECTION(OfflinePrecursorIonSelection())
	ptr = new OfflinePrecursorIonSelection();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~OfflinePrecursorIonSelection())
	delete ptr;
END_SECTION

ptr = new OfflinePrecursorIonSelection();
std::vector<PeptideIdentification> pep_ids;
std::vector<ProteinIdentification> prot_ids;
//IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("OfflinePrecursorIonSelection_ids.IdXML"),prot_ids,pep_ids);

FeatureMap<> map;
FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("OfflinePrecursorIonSelection_features.featureXML"),map);
MSExperiment<> raw_data;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OfflinePrecursorIonSelection_raw_data.mzML"),raw_data);


START_SECTION((template < typename InputPeakType > void makePrecursorSelectionForKnownLCMSMap(const FeatureMap<> &features, const MSExperiment< InputPeakType > &experiment, MSExperiment< InputPeakType > &ms2, std::set< Int > &charges_set, bool feature_based)))
{
	MSExperiment<Peak1D> ms2;
	std::set<Int> charges_set;
	charges_set.insert(1);
	bool feature_based = true;
	Param param;
	param.setValue("ms2_spectra_per_rt_bin",1);
	ptr->setParameters(param);
	ptr->makePrecursorSelectionForKnownLCMSMap(map,raw_data,ms2,charges_set,feature_based);
	TEST_EQUAL(ms2.size(),3)
	TEST_REAL_SIMILAR(ms2[0].getRT(),45)
	TEST_REAL_SIMILAR(ms2[0].getPrecursors()[0].getMZ(),336.14)
	TEST_REAL_SIMILAR(ms2[1].getRT(),55)
	TEST_REAL_SIMILAR(ms2[1].getPrecursors()[0].getMZ(),319.19)
	TEST_REAL_SIMILAR(ms2[2].getRT(),65)
	TEST_REAL_SIMILAR(ms2[2].getPrecursors()[0].getMZ(),478.29)
		
	ms2.clear(true);
	feature_based = false;
	ptr->makePrecursorSelectionForKnownLCMSMap(map,raw_data,ms2,charges_set,feature_based);
	TEST_EQUAL(ms2.size(),3)
	TEST_REAL_SIMILAR(ms2[0].getRT(),45)
	TEST_REAL_SIMILAR(ms2[0].getPrecursors()[0].getMZ(),336.14)
	TEST_REAL_SIMILAR(ms2[1].getRT(),55)
	TEST_REAL_SIMILAR(ms2[1].getPrecursors()[0].getMZ(),336.14)
	TEST_REAL_SIMILAR(ms2[2].getRT(),65)
	TEST_REAL_SIMILAR(ms2[2].getPrecursors()[0].getMZ(),336.14)

	ms2.clear(true);
	feature_based = true;
	param.setValue("exclude_overlapping_peaks","true");
	param.setValue("min_peak_distance",40.);
	ptr->setParameters(param);
	ptr->makePrecursorSelectionForKnownLCMSMap(map,raw_data,ms2,charges_set,feature_based);
	TEST_EQUAL(ms2.size(),2)
	TEST_REAL_SIMILAR(ms2[0].getRT(),45)
	TEST_REAL_SIMILAR(ms2[0].getPrecursors()[0].getMZ(),336.14)
	TEST_REAL_SIMILAR(ms2[1].getRT(),65)
	TEST_REAL_SIMILAR(ms2[1].getPrecursors()[0].getMZ(),478.29)
		
}
END_SECTION	     

START_SECTION((template < typename InputPeakType > void getMassRanges(const FeatureMap<> &features, const MSExperiment< InputPeakType > &experiment, std::vector< std::vector< std::pair< Size, Size > > > &indices)))
{
	Param param;
	param.setValue("exclude_overlapping_peaks","false");
	ptr->setParameters(param);
	std::vector<std::vector<std::pair<Size,Size> > >  indices;
	ptr->getMassRanges(map,raw_data,indices);
	TEST_EQUAL(indices.size(),3)
	TEST_EQUAL(indices[0][0].first,0)
	TEST_EQUAL(indices[0][0].second,0)
	TEST_EQUAL(indices[0][1].second,0)
	TEST_EQUAL(indices[1][0].first,1)
	TEST_EQUAL(indices[1][0].second,0)
	TEST_EQUAL(indices[1][1].second,0)	
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



