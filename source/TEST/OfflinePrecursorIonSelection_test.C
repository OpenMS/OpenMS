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
START_SECTION(OfflinePrecursorIonSelection())
{
	ptr = new OfflinePrecursorIonSelection();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~OfflinePrecursorIonSelection())
{
	delete ptr;
}
END_SECTION

ptr = new OfflinePrecursorIonSelection();

std::vector<PeptideIdentification> pep_ids;
std::vector<ProteinIdentification> prot_ids;
//IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("OfflinePrecursorIonSelection_ids.IdXML"),prot_ids,pep_ids);

FeatureMap<> map;
FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("OfflinePrecursorIonSelection_features.featureXML"),map);
MSExperiment<> raw_data;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OfflinePrecursorIonSelection_raw_data.mzML"),raw_data);

START_SECTION((void computeOptimalSolution(std::vector<ProteinIdentification>& prot_ids,std::vector<PeptideIdentification>& pep_ids,FeatureMap<>& features,FeatureMap<>& optimal_set,bool filter)))
{
	FeatureMap<> optimal_map;
}
END_SECTION

START_SECTION((void makePrecursorSelectionForKnownLCMSMap(FeatureMap<>& features,MSExperiment< Peak1D > & experiment,MSExperiment< Peak1D > & ms2,std::set<Int>& charges_set,bool feature_based)))
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
	TEST_REAL_SIMILAR(ms2[0].getRT(),40)
	TEST_REAL_SIMILAR(ms2[0].getPrecursors()[0].getMZ(),336.14)
	TEST_REAL_SIMILAR(ms2[1].getRT(),50)
	TEST_REAL_SIMILAR(ms2[1].getPrecursors()[0].getMZ(),319.19)
	TEST_REAL_SIMILAR(ms2[2].getRT(),60)
	TEST_REAL_SIMILAR(ms2[2].getPrecursors()[0].getMZ(),478.29)
		
	ms2.clear();
	feature_based = false;
	ptr->makePrecursorSelectionForKnownLCMSMap(map,raw_data,ms2,charges_set,feature_based);
	TEST_EQUAL(ms2.size(),3)
	TEST_REAL_SIMILAR(ms2[0].getRT(),40)
	TEST_REAL_SIMILAR(ms2[0].getPrecursors()[0].getMZ(),336.14)
	TEST_REAL_SIMILAR(ms2[1].getRT(),50)
	TEST_REAL_SIMILAR(ms2[1].getPrecursors()[0].getMZ(),336.14)
	TEST_REAL_SIMILAR(ms2[2].getRT(),60)
	TEST_REAL_SIMILAR(ms2[2].getPrecursors()[0].getMZ(),336.14)
	
}
END_SECTION	     

START_SECTION((void getMassRanges(FeatureMap<>& features, MSExperiment<>& experiment,std::vector<std::vector<std::pair<Size,Size> > > & indices)))
{
	
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

START_SECTION((void computeOptimalSolution(std::vector<ProteinIdentification>& prot_ids,std::vector<PeptideIdentification>& pep_ids,MSExperiment<>& experiment,FeatureMap<>& features,FeatureMap<>& optimal_set,bool filter)))
{
	
}
END_SECTION	      

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



