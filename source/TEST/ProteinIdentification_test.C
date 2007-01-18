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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <string>

#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

///////////////////////////

START_TEST(ProteinIdentification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ProteinIdentification* ptr1 = 0;
ProteinIdentification* ptr2 = 0;
float protein_significance_threshold = 63.2;
std::vector<ProteinHit> protein_hits;
ProteinHit protein_hit;
ProteinIdentification protein_identification;
DateTime date;
MascotXMLFile xml_file;

date.now();

protein_hits.push_back(protein_hit);

CHECK((ProteinIdentification()))
	ptr1 = new ProteinIdentification();
	TEST_NOT_EQUAL(ptr1, 0)
	TEST_NOT_EQUAL(&(ptr1->getDateTime()), 0)
	TEST_EQUAL(ptr1->getProteinSignificanceThreshold(), 0)
	TEST_NOT_EQUAL(&(ptr1->getProteinHits()), 0)
RESULT

CHECK((ProteinIdentification(const ProteinIdentification& source)))
	ptr1 = new ProteinIdentification();
	ptr1->setDateTime(date);
	ptr1->setProteinSignificanceThreshold(protein_significance_threshold);
	ptr1->insertProteinHit(protein_hit);
	ptr2 = new ProteinIdentification(*ptr1);
	TEST_EQUAL(ptr1->getDateTime() == ptr2->getDateTime(), true)
	TEST_EQUAL(ptr1->getProteinSignificanceThreshold(), ptr2->getProteinSignificanceThreshold())
	TEST_EQUAL(ptr1->getProteinHits().size() == 1, true)
	TEST_EQUAL(ptr1->getProteinHits()[0].getSequence(), String(""))	
	TEST_EQUAL(ptr1->getProteinHits()[0] == protein_hit, true)	
RESULT

CHECK((~ProteinIdentification()))
	ptr1 = new ProteinIdentification();
	delete ptr1;
RESULT

CHECK((ProteinIdentification& operator=(const ProteinIdentification& source)))
	ptr1 = new ProteinIdentification();
	ptr1->setDateTime(date);
	ptr1->setProteinSignificanceThreshold(protein_significance_threshold);
	ptr1->insertProteinHit(protein_hit);
	*ptr2 = *ptr1;
	TEST_EQUAL(ptr1->getDateTime() == ptr2->getDateTime(), true)
	TEST_EQUAL(ptr1->getProteinSignificanceThreshold(), ptr2->getProteinSignificanceThreshold())
	TEST_EQUAL(ptr1->getProteinHits().size() == 1, true)
	TEST_EQUAL(*(ptr1->getProteinHits().begin()) == protein_hit, true)	
RESULT

CHECK((bool operator != (const ProteinIdentification& rhs) const))
	ProteinIdentification search1;
	ProteinIdentification search2;
	
	TEST_EQUAL(search1 != search2, false)
	search1.setDateTime(date);
	TEST_EQUAL(search1 != search2, true)	
	search2.setDateTime(date);
	TEST_EQUAL(search1 != search2, false)	
	search1.setProteinSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(search1 != search2, true)	
	search2.setProteinSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(search1 != search2, false)	
RESULT

CHECK((bool operator == (const ProteinIdentification& rhs) const))
	ProteinIdentification search1;
	ProteinIdentification search2;
	
	TEST_EQUAL(search1 == search2, true)
	search1.setDateTime(date);
	TEST_EQUAL(search1 == search2, false)	
	search2.setDateTime(date);
	TEST_EQUAL(search1 == search2, true)	
	search1.setProteinSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(search1 == search2, false)	
	search2.setProteinSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(search1 == search2, true)	
RESULT

CHECK((template<typename Arg> bool operator()(const Arg& a, const Arg& b)))
	ProteinHit hit_1;
	ProteinHit hit_2;
	ProteinHit hit_3;
	vector<ProteinHit> hits;
	ProteinIdentification id;
	
	hit_1.setScore(23);
	hit_2.setScore(11);
	hit_3.setScore(45);

	hit_1.setScoreType("Mascot");
	hit_2.setScoreType("Mascot");
	hit_3.setScoreType("Mascot");
	hit_1.setAccession("SECONDPROTEIN");
	hit_2.setAccession("THIRDPROTEIN");
	hit_3.setAccession("FIRSTPROTEIN");
	hits.push_back(hit_1);
	hits.push_back(hit_2);
	hits.push_back(hit_3);
	id.setProteinHits(hits);
	id.assignRanks();
	TEST_EQUAL(id.getProteinHits()[0].getAccession(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getProteinHits()[1].getAccession(), "SECONDPROTEIN")
	TEST_EQUAL(id.getProteinHits()[2].getAccession(), "THIRDPROTEIN")
	TEST_EQUAL(id.getProteinHits()[0].getRank(), 1)	
	TEST_EQUAL(id.getProteinHits()[1].getRank(), 2)
	TEST_EQUAL(id.getProteinHits()[2].getRank(), 3)

RESULT

CHECK((DateTime& getDateTime()))
	ptr1 = new ProteinIdentification();
	ptr1->setDateTime(date);
	TEST_EQUAL(ptr1->getDateTime() == date, true)  
RESULT

CHECK((const DateTime& getDateTime() const))
	ptr1 = new ProteinIdentification();
	ptr1->setDateTime(date);
	const DateTime& date_time = ptr1->getDateTime();
	TEST_EQUAL(date_time == date, true)  
RESULT

CHECK((float getProteinSignificanceThreshold() const))
	ptr1 = new ProteinIdentification();
	ptr1->setProteinSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(ptr1->getProteinSignificanceThreshold(), protein_significance_threshold)	
RESULT

CHECK((const std::vector<ProteinHit>& getProteinHits() const))
	ptr1 = new ProteinIdentification();
	ptr1->insertProteinHit(protein_hit);
	TEST_EQUAL(ptr1->getProteinHits().size() == 1, true)
	TEST_EQUAL(*(ptr1->getProteinHits().begin()) == protein_hit, true)	
RESULT

CHECK((void clear()))
	DateTime test_date;
	ptr1 = new ProteinIdentification();
	ptr1->setDateTime(date);
	ptr1->setProteinSignificanceThreshold(protein_significance_threshold);
	ptr1->clear();
	TEST_EQUAL(ptr1->getDateTime() == test_date, true)
	TEST_EQUAL(ptr1->getProteinSignificanceThreshold(), 0)
	TEST_EQUAL(ptr1->getProteinHits().size() == 0, true)
RESULT

CHECK((void insertProteinHit(const ProteinHit& input)))
	ptr1 = new ProteinIdentification();
	ptr1->insertProteinHit(protein_hit);
	TEST_EQUAL(ptr1->getProteinHits().size() == 1, true)
	TEST_EQUAL(*(ptr1->getProteinHits().begin()) == protein_hit, true)
RESULT

CHECK((void setDateTime(const DateTime& date)))
	ptr1 = new ProteinIdentification();
	ptr1->setDateTime(date);
	TEST_EQUAL(ptr1->getDateTime() == date, true)
RESULT

CHECK((void setProteinSignificanceThreshold(float value)))
	ptr1 = new ProteinIdentification();
	ptr1->setProteinSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(ptr1->getProteinSignificanceThreshold(), protein_significance_threshold)
RESULT

CHECK((bool empty() const))
	ptr1 = new ProteinIdentification();
	TEST_EQUAL(ptr1->empty(), true)
	ptr1->setProteinSignificanceThreshold(1);
	TEST_EQUAL(ptr1->empty(), false)
	ptr1->setProteinSignificanceThreshold(0);
	TEST_EQUAL(ptr1->empty(), true)
	ptr1->clear();
	ptr1->insertProteinHit(protein_hit);
	TEST_EQUAL(ptr1->empty(), false)
	ptr1->clear();
	ptr1->setDateTime(date);	
	TEST_EQUAL(ptr1->empty(), false)
	ptr1->clear();
RESULT

CHECK((void sort()))

	vector<ProteinIdentification> protein_identifications; 
	vector<IdentificationData> identifications; 
	ProteinHit hit;
	
	hit.setAccession("TESTPROTEIN");
	hit.setScore(80.4);
	hit.setScoreType("Mascot");

	AnalysisXMLFile().load("data/AnalysisXMLFile_test.analysisXML",
							protein_identifications, 
				   		identifications);
							
	protein_identifications[0].insertProteinHit(hit);
	protein_identifications[0].sort();

	TEST_EQUAL(protein_identifications.size(), 1)
	TEST_EQUAL(protein_identifications[0].getProteinHits()[1].getAccession(), "TESTPROTEIN")
	
RESULT

CHECK((void assignRanks()))
	ProteinHit hit_1;
	ProteinHit hit_2;
	ProteinHit hit_3;
	vector<ProteinHit> hits;
	ProteinIdentification id;
	
	hit_1.setScore(23);
	hit_2.setScore(11);
	hit_3.setScore(45);
	hit_1.setScoreType("Mascot");
	hit_2.setScoreType("Mascot");
	hit_3.setScoreType("Mascot");
	hit_1.setAccession("SECONDPROTEIN");
	hit_2.setAccession("THIRDPROTEIN");
	hit_3.setAccession("FIRSTPROTEIN");
	hits.push_back(hit_1);
	hits.push_back(hit_2);
	hits.push_back(hit_3);
	id.setProteinHits(hits);
	id.assignRanks();
	TEST_EQUAL(id.getProteinHits()[0].getAccession(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getProteinHits()[1].getAccession(), "SECONDPROTEIN")
	TEST_EQUAL(id.getProteinHits()[2].getAccession(), "THIRDPROTEIN")
	TEST_EQUAL(id.getProteinHits()[0].getRank(), 1)	
	TEST_EQUAL(id.getProteinHits()[1].getRank(), 2)
	TEST_EQUAL(id.getProteinHits()[2].getRank(), 3)
	
RESULT

CHECK((void setProteinHits(const std::vector<ProteinHit>& protein_hits)))
	ProteinHit hit_1;
	ProteinHit hit_2;
	ProteinHit hit_3;
	vector<ProteinHit> hits;
	ProteinIdentification id;
	
	hit_1.setScore(23);
	hit_2.setScore(11);
	hit_3.setScore(45);
	hit_1.setScoreType("Mascot");
	hit_2.setScoreType("Mascot");
	hit_3.setScoreType("Mascot");
	hit_1.setAccession("SECONDPROTEIN");
	hit_2.setAccession("THIRDPROTEIN");
	hit_3.setAccession("FIRSTPROTEIN");
	hits.push_back(hit_1);
	hits.push_back(hit_2);
	hits.push_back(hit_3);
	id.setProteinHits(hits);
	TEST_EQUAL(id.getProteinHits()[2].getAccession(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getProteinHits()[0].getAccession(), "SECONDPROTEIN")
	TEST_EQUAL(id.getProteinHits()[1].getAccession(), "THIRDPROTEIN")
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
