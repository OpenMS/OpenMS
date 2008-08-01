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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <string>

#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

///////////////////////////

START_TEST(Identification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

float protein_significance_threshold = 63.2;
std::vector<ProteinHit> protein_hits;
ProteinHit protein_hit;
ProteinIdentification protein_identification;
DateTime date;
MascotXMLFile xml_file;

date.now();

protein_hits.push_back(protein_hit);

ProteinIdentification* ptr = 0;
CHECK((ProteinIdentification()))
	ptr = new ProteinIdentification();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~ProteinIdentification()))
	ProteinIdentification hits;
	delete ptr;
RESULT

CHECK((ProteinIdentification(const ProteinIdentification &source)))
	ProteinIdentification hits;
	hits.setDateTime(date);
	hits.setSignificanceThreshold(protein_significance_threshold);
	hits.insertHit(protein_hit);
	hits.setMetaValue("label",17);
	hits.setIdentifier("id");
	hits.setScoreType("score_type");
	hits.setHigherScoreBetter(false);
	hits.setSearchEngine("Mascot");
	hits.setSearchEngineVersion("2.1");
	ProteinIdentification::SearchParameters param;
	param.db = "RefSeq";
	hits.setSearchParameters(param);
	
	ProteinIdentification hits2(hits);

	TEST_EQUAL(hits.getDateTime() == hits2.getDateTime(), true)
	TEST_EQUAL(hits.getSignificanceThreshold(), hits2.getSignificanceThreshold())
	TEST_EQUAL(hits.getHits().size() == 1, true)
	TEST_EQUAL(hits.getHits()[0].getSequence(), String(""))	
	TEST_EQUAL(hits.getHits()[0] == protein_hit, true)	
	TEST_EQUAL((UInt)hits.getMetaValue("label"),17)
	TEST_EQUAL(hits.getIdentifier(),"id")
	TEST_EQUAL(hits.getScoreType(),"score_type")
	TEST_EQUAL(hits.isHigherScoreBetter(),false)
	TEST_EQUAL(hits.getSearchEngine(), "Mascot")
	TEST_EQUAL(hits.getSearchEngineVersion(), "2.1")
	TEST_EQUAL(hits.getSearchParameters()==param, true)
RESULT


CHECK((ProteinIdentification& operator=(const ProteinIdentification& source)))
	ProteinIdentification hits;
	hits.setDateTime(date);
	hits.setSignificanceThreshold(protein_significance_threshold);
	hits.insertHit(protein_hit);
	hits.setIdentifier("id");
	hits.setScoreType("score_type");
	hits.setHigherScoreBetter(false);
	hits.setSearchEngine("Mascot");
	hits.setSearchEngineVersion("2.1");
	ProteinIdentification::SearchParameters param;
	param.db = "RefSeq";
	hits.setSearchParameters(param);
		
	ProteinIdentification hits2;
	hits2 = hits;
	
	TEST_EQUAL(hits.getDateTime() == hits2.getDateTime(), true)
	TEST_EQUAL(hits.getSignificanceThreshold(), hits2.getSignificanceThreshold())
	TEST_EQUAL(hits.getHits().size() == 1, true)
	TEST_EQUAL(*(hits.getHits().begin()) == protein_hit, true)
	TEST_EQUAL(hits.getIdentifier(),"id")
	TEST_EQUAL(hits.getScoreType(),"score_type")
	TEST_EQUAL(hits.isHigherScoreBetter(),false)
	TEST_EQUAL(hits.getSearchEngine(), "Mascot")
	TEST_EQUAL(hits.getSearchEngineVersion(), "2.1")
	TEST_EQUAL(hits.getSearchParameters()==param, true)
RESULT

CHECK((bool operator == (const ProteinIdentification& rhs) const))
	ProteinIdentification search1;
	ProteinIdentification search2;
	TEST_EQUAL(search1 == search2, true)
	
	search1.setDateTime(date);
	TEST_EQUAL(search1 == search2, false)	
	search1 = search2;
		
	search1.setSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(search1 == search2, false)	
	search1 = search2;

	search2.setIdentifier("id");
	TEST_EQUAL(search1 == search2, false)	
	search1 = search2;

	search2.setScoreType("score_type");
	TEST_EQUAL(search1 == search2, false)	
	search1 = search2;

	search2.setHigherScoreBetter(false);
	TEST_EQUAL(search1 == search2, false)	
	search1 = search2;

	search2.setSearchEngine("Mascot");
	TEST_EQUAL(search1 == search2, false)	
	search1 = search2;
	
	search2.setSearchEngineVersion("2.1");
	TEST_EQUAL(search1 == search2, false)	
	search1 = search2;
	
	ProteinIdentification::SearchParameters param;
	param.db = "RefSeq";
	search2.setSearchParameters(param);
	TEST_EQUAL(search1 == search2, false)	
	search1 = search2;

RESULT

CHECK((bool operator != (const ProteinIdentification& rhs) const))
	ProteinIdentification search1;
	ProteinIdentification search2;
	TEST_EQUAL(search1 != search2, false)
	
	search1.setDateTime(date);
	TEST_EQUAL(search1 != search2, true)	

	//rest does not need to be tested, as it is tested in the operator== test implicitly!
RESULT

CHECK((const DateTime& getDateTime() const))
	ProteinIdentification hits;
	hits.setDateTime(date);
	const DateTime& date_time = hits.getDateTime();
	TEST_EQUAL(date_time == date, true)  
RESULT

CHECK((Real getSignificanceThreshold() const))
	ProteinIdentification hits;
	hits.setSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(hits.getSignificanceThreshold(), protein_significance_threshold)	
RESULT

CHECK((const std::vector<ProteinHit>& getHits() const))
	ProteinIdentification hits;
	hits.insertHit(protein_hit);
	TEST_EQUAL(hits.getHits().size() == 1, true)
	TEST_EQUAL(*(hits.getHits().begin()) == protein_hit, true)	
RESULT

CHECK((void insertHit(const ProteinHit& input)))
	ProteinIdentification hits;
	hits.insertHit(protein_hit);
	TEST_EQUAL(hits.getHits().size() == 1, true)
	TEST_EQUAL(*(hits.getHits().begin()) == protein_hit, true)
RESULT

CHECK((void setDateTime(const DateTime& date)))
	ProteinIdentification hits;
	hits.setDateTime(date);
	TEST_EQUAL(hits.getDateTime() == date, true)
RESULT

CHECK((void setSignificanceThreshold(Real value)))
	ProteinIdentification hits;
	hits.setSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(hits.getSignificanceThreshold(), protein_significance_threshold)
RESULT

CHECK((void setHits(const std::vector< ProteinHit > &hits)))
	ProteinHit hit_1;
	ProteinHit hit_2;
	ProteinHit hit_3;
	vector<ProteinHit> hits;
	ProteinIdentification id;
	
	hit_1.setScore(23);
	hit_2.setScore(11);
	hit_3.setScore(45);
	hit_1.setAccession("SECONDPROTEIN");
	hit_2.setAccession("THIRDPROTEIN");
	hit_3.setAccession("FIRSTPROTEIN");
	hits.push_back(hit_1);
	hits.push_back(hit_2);
	hits.push_back(hit_3);
	id.setHits(hits);
	TEST_EQUAL(id.getHits()[2].getAccession(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getHits()[0].getAccession(), "SECONDPROTEIN")
	TEST_EQUAL(id.getHits()[1].getAccession(), "THIRDPROTEIN")
RESULT

CHECK((const String& getScoreType() const))
	ProteinIdentification hits;
	TEST_EQUAL(hits.getScoreType(),"")
RESULT

CHECK((void setScoreType(const String& type)))
	ProteinIdentification hits;
	hits.setScoreType("bla");
	TEST_EQUAL(hits.getScoreType(),"bla")
RESULT

CHECK((bool isHigherScoreBetter() const))
	ProteinIdentification hits;
	TEST_EQUAL(hits.isHigherScoreBetter(),true)
RESULT

CHECK((void setHigherScoreBetter(bool higher_is_better)))
	ProteinIdentification hits;
	hits.setHigherScoreBetter(false);
	TEST_EQUAL(hits.isHigherScoreBetter(),false)
RESULT

CHECK((const String& getIdentifier() const))
	ProteinIdentification hits;
	TEST_EQUAL(hits.getIdentifier(),"")
RESULT

CHECK((void setIdentifier(const String& id)))
	ProteinIdentification hits;
	hits.setIdentifier("bla");
	TEST_EQUAL(hits.getIdentifier(),"bla")
RESULT

CHECK((const String& getSearchEngine() const))
	ProteinIdentification hits;
	TEST_EQUAL(hits.getSearchEngine(),"")
RESULT

CHECK((void setSearchEngine(const String &search_engine)))
	ProteinIdentification hits;
	hits.setIdentifier("bla");
	TEST_EQUAL(hits.getIdentifier(),"bla")
RESULT

CHECK((const String& getSearchEngineVersion() const))
	ProteinIdentification hits;
	TEST_EQUAL(hits.getSearchEngineVersion(),"")
RESULT

CHECK((void setSearchEngineVersion(const String &search_engine_version)))
	ProteinIdentification hits;
	hits.setSearchEngineVersion("bla");
	TEST_EQUAL(hits.getSearchEngineVersion(),"bla")
RESULT

CHECK((const SearchParameters& getSearchParameters() const))
	ProteinIdentification hits;
	TEST_EQUAL(hits.getSearchParameters()==ProteinIdentification::SearchParameters(),true)
RESULT

CHECK((void setSearchParameters(const SearchParameters &search_parameters)))
	ProteinIdentification hits;
	ProteinIdentification::SearchParameters param;
	param.db="Mascot";
	hits.setSearchParameters(param);
	TEST_EQUAL(hits.getSearchParameters()==ProteinIdentification::SearchParameters(),false)
RESULT

CHECK((void sort()))
	ProteinIdentification id;
	ProteinHit hit;
	hit.setScore(23);
	hit.setAccession("SECONDPROTEIN");
	id.insertHit(hit);
	hit.setScore(45);
	hit.setAccession("FIRSTPROTEIN");
	id.insertHit(hit);
	hit.setScore(7);
	hit.setAccession("THIRDPROTEIN");
	id.insertHit(hit);
	
	//higher score is better
	id.sort();

	TEST_EQUAL(id.getHits()[0].getAccession(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getHits()[1].getAccession(), "SECONDPROTEIN")
	TEST_EQUAL(id.getHits()[2].getAccession(), "THIRDPROTEIN")
	TEST_EQUAL(id.getHits()[0].getScore(), 45)	
	TEST_EQUAL(id.getHits()[1].getScore(), 23)
	TEST_EQUAL(id.getHits()[2].getScore(), 7)

	//lower score is better
	id.setHigherScoreBetter(false);
	id.sort();

	TEST_EQUAL(id.getHits()[0].getAccession(), "THIRDPROTEIN")
	TEST_EQUAL(id.getHits()[1].getAccession(), "SECONDPROTEIN")
	TEST_EQUAL(id.getHits()[2].getAccession(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getHits()[0].getScore(), 7)	
	TEST_EQUAL(id.getHits()[1].getScore(), 23)
	TEST_EQUAL(id.getHits()[2].getScore(), 45)
RESULT

CHECK((void assignRanks()))
	ProteinIdentification id;
	ProteinHit hit;
	hit.setScore(23);
	hit.setAccession("SECONDPROTEIN");
	id.insertHit(hit);
	hit.setScore(45);
	hit.setAccession("FIRSTPROTEIN");
	id.insertHit(hit);
	hit.setScore(7);
	hit.setAccession("THIRDPROTEIN");
	id.insertHit(hit);

	id.assignRanks();

	TEST_EQUAL(id.getHits()[0].getAccession(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getHits()[1].getAccession(), "SECONDPROTEIN")
	TEST_EQUAL(id.getHits()[2].getAccession(), "THIRDPROTEIN")
	TEST_EQUAL(id.getHits()[0].getRank(), 1)	
	TEST_EQUAL(id.getHits()[1].getRank(), 2)
	TEST_EQUAL(id.getHits()[2].getRank(), 3)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
