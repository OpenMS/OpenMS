// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

///////////////////////////

START_TEST(IDFeatureMapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


IDFeatureMapper* ptr = 0;
CHECK((IDFeatureMapper()))
	ptr = new IDFeatureMapper();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(([EXTRA]~IDFeatureMapper()))
	delete ptr;
RESULT

CHECK((void annotate(FeatureMap<> &fm, const std::vector< PeptideIdentification > &ids, const std::vector< ProteinIdentification > &protein_ids) ))
	//load id data
	vector<PeptideIdentification> identifications; 
	vector<ProteinIdentification> protein_identifications; 
	IdXMLFile().load("data/IDFeatureMapper_test.idXML", protein_identifications, identifications);
	
	//load feature data
	FeatureMap<> fm;
	FeatureXMLFile().load("data/IDFeatureMapper_test.featureXML", fm);
	
	//map
	IDFeatureMapper().annotate(fm,identifications,protein_identifications);
	
	//test protein ids
	TEST_EQUAL(fm.getProteinIdentifications().size(),1)
	TEST_EQUAL(fm.getProteinIdentifications()[0].getHits().size(),2)
	TEST_EQUAL(fm.getProteinIdentifications()[0].getHits()[0].getAccession(),"ABCDE")
	TEST_EQUAL(fm.getProteinIdentifications()[0].getHits()[1].getAccession(),"FGHIJ")
	
	//test peptide ids
	TEST_EQUAL(fm[0].getPeptideIdentifications().size(),5)
	TEST_EQUAL(fm[1].getPeptideIdentifications().size(),0)
	TEST_EQUAL(fm[2].getPeptideIdentifications().size(),0)
	TEST_EQUAL(fm[0].getPeptideIdentifications()[0].getHits().size(),1)
	TEST_EQUAL(fm[0].getPeptideIdentifications()[1].getHits().size(),1)
	TEST_EQUAL(fm[0].getPeptideIdentifications()[2].getHits().size(),1)
	TEST_EQUAL(fm[0].getPeptideIdentifications()[3].getHits().size(),1)
	TEST_EQUAL(fm[0].getPeptideIdentifications()[4].getHits().size(),1)
	TEST_EQUAL(fm[0].getPeptideIdentifications()[0].getHits()[0].getSequence(),"A")
	TEST_EQUAL(fm[0].getPeptideIdentifications()[1].getHits()[0].getSequence(),"B")
	TEST_EQUAL(fm[0].getPeptideIdentifications()[2].getHits()[0].getSequence(),"C")
	TEST_EQUAL(fm[0].getPeptideIdentifications()[3].getHits()[0].getSequence(),"D")
	TEST_EQUAL(fm[0].getPeptideIdentifications()[4].getHits()[0].getSequence(),"E")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
