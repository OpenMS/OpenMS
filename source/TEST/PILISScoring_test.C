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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/ID/PILISScoring.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

///////////////////////////

START_TEST(PILISScoring_test.C, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PILISScoring* ptr = 0;
String filename(OPENMS_GET_TEST_DATA_PATH("IDFilter_test2.idXML"));
START_SECTION(PILISScoring())
	ptr = new PILISScoring();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~PILISScoring())
	delete ptr;
END_SECTION

ptr = new PILISScoring();

START_SECTION(PILISScoring(const PILISScoring& source))
	PILISScoring copy(*ptr);
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(PILISScoring& operator = (const PILISScoring& source))
	PILISScoring copy;
	copy = *ptr;
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(void getScores(std::vector<PeptideIdentification>& ids))
	vector<PeptideIdentification> ids;
	vector<ProteinIdentification> prot_ids;
	String document_id;
	IdXMLFile().load(filename, prot_ids, ids, document_id);
	ptr->getScores(ids);
	for (vector<PeptideIdentification>::const_iterator it = ids.begin(); it != ids.end(); ++it)
	{
		TEST_EQUAL(it->getScoreType(), "PILIS-E-value")
	}
END_SECTION

START_SECTION(void getScore(PeptideIdentification& id))
	vector<PeptideIdentification> ids;
	vector<ProteinIdentification> prot_ids;
	String document_id;
	IdXMLFile().load(filename, prot_ids, ids, document_id);
	ptr->getScore(ids[0]);
	TEST_REAL_SIMILAR(ids[0].getHits().begin()->getScore(), 33.85)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
