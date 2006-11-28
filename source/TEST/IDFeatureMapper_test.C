// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <iostream>

#include <OpenMS/ANALYSIS/ID/IDFeatureMapper.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>

///////////////////////////

START_TEST(IDFeatureMapper, "$Id: IDSpectrumMapper_test.C 500 2006-09-06 16:39:00Z marc_sturm $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


IDFeatureMapper* ptr = 0;
CHECK(IDFeatureMapper())
	ptr = new IDFeatureMapper();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~IDFeatureMapper())
	delete ptr;
RESULT

CHECK((void annotate(DFeatureMap<2> fm, const std::vector<Identification>& identifications, const std::vector<float>& precursor_retention_times, const std::vector<float>& precursor_mz_values)))

//load id data
vector<Identification> identifications; 
vector<ProteinIdentification> protein_identifications; 
vector<float> precursor_retention_times;
vector<float> precursor_mz_values;
ContactPerson contact_person;
AnalysisXMLFile().load("data/IDFeatureMapper_test.analysisXML",
												protein_identifications, 
									   		identifications, 
												precursor_retention_times, 
												precursor_mz_values, 
												contact_person);

//load feature data
DFeatureMap<2> fm;
DFeatureMapFile().load("data/IDFeatureMapper_test.feat", fm);

//map
IDFeatureMapper().annotate(fm,identifications,precursor_retention_times,precursor_mz_values);

//TODO
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
