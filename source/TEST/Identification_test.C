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
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

///////////////////////////

START_TEST(Identification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

Identification* ptr1 = 0;
Identification* ptr2 = 0;
float peptide_significance_threshold = 42.3;
float protein_significance_threshold = 63.2;
std::vector<PeptideHit> peptide_hits;
std::vector<ProteinHit> protein_hits;
ProteinHit protein_hit;
PeptideHit peptide_hit;
ProteinIdentification protein_identification;
vector<IdentificationData> identifications; 
DateTime date;
MascotXMLFile xml_file;
vector<PeptideHit>* hits;

date.now();

peptide_hits.push_back(peptide_hit);
protein_hits.push_back(protein_hit);

CHECK((Identification()))
	ptr1 = new Identification();
	TEST_NOT_EQUAL(ptr1, 0)
	TEST_NOT_EQUAL(&(ptr1->getDateTime()), 0)
	TEST_EQUAL(ptr1->getProteinSignificanceThreshold(), 0)
	TEST_EQUAL(ptr1->getPeptideSignificanceThreshold(), 0)
	TEST_NOT_EQUAL(&(ptr1->getProteinHits()), 0)
	TEST_NOT_EQUAL(&(ptr1->getPeptideHits()), 0)
RESULT

CHECK((Identification(const Identification& source)))
	ptr1 = new Identification();
	ptr1->setDateTime(date);
	ptr1->setPeptideSignificanceThreshold(peptide_significance_threshold);
	ptr1->setProteinSignificanceThreshold(protein_significance_threshold);
	ptr1->setPeptideAndProteinHits(peptide_hits, protein_hits);
	ptr2 = new Identification(*ptr1);
	TEST_EQUAL(ptr1->getDateTime() == ptr2->getDateTime(), true)
	TEST_EQUAL(ptr1->getProteinSignificanceThreshold(), ptr2->getProteinSignificanceThreshold())
	TEST_EQUAL(ptr1->getPeptideSignificanceThreshold(), ptr2->getPeptideSignificanceThreshold())
	TEST_EQUAL(ptr1->getPeptideHits().size() == 1, true)
	TEST_EQUAL(*(ptr1->getPeptideHits().begin()) == peptide_hit, true)	
	TEST_EQUAL(ptr1->getProteinHits().size() == 1, true)
	TEST_EQUAL(*(ptr1->getProteinHits().begin()) == protein_hit, true)	
RESULT

CHECK((~Identification()))
	ptr1 = new Identification();
	delete ptr1;
RESULT

CHECK((Identification& operator=(const Identification& source)))
	ptr1 = new Identification();
	ptr1->setDateTime(date);
	ptr1->setPeptideSignificanceThreshold(peptide_significance_threshold);
	ptr1->setProteinSignificanceThreshold(protein_significance_threshold);
	ptr1->setPeptideAndProteinHits(peptide_hits, protein_hits);
	*ptr2 = *ptr1;
	TEST_EQUAL(ptr1->getDateTime() == ptr2->getDateTime(), true)
	TEST_EQUAL(ptr1->getProteinSignificanceThreshold(), ptr2->getProteinSignificanceThreshold())
	TEST_EQUAL(ptr1->getPeptideSignificanceThreshold(), ptr2->getPeptideSignificanceThreshold())
	TEST_EQUAL(ptr1->getPeptideHits().size() == 1, true)
	TEST_EQUAL(*(ptr1->getPeptideHits().begin()) == peptide_hit, true)	
	TEST_EQUAL(ptr1->getProteinHits().size() == 1, true)
	TEST_EQUAL(*(ptr1->getProteinHits().begin()) == protein_hit, true)	
RESULT

CHECK((bool operator != (const Identification& rhs) const))
	Identification search1;
	Identification search2;
	
	TEST_EQUAL(search1 != search2, false)
	TEST_EQUAL(search1 != search2, false)	
	search1.setDateTime(date);
	TEST_EQUAL(search1 != search2, true)	
	search2.setDateTime(date);
	TEST_EQUAL(search1 != search2, false)	
	search1.setPeptideSignificanceThreshold(peptide_significance_threshold);
	TEST_EQUAL(search1 != search2, true)	
	search2.setPeptideSignificanceThreshold(peptide_significance_threshold);
	TEST_EQUAL(search1 != search2, false)	
	search1.setProteinSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(search1 != search2, true)	
	search2.setProteinSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(search1 != search2, false)	
	search1.setPeptideAndProteinHits(peptide_hits, protein_hits);
	TEST_EQUAL(search1 != search2, true)	
	search2.setPeptideAndProteinHits(peptide_hits, protein_hits);
	TEST_EQUAL(search1 != search2, false)	
RESULT

CHECK((bool operator == (const Identification& rhs) const))
	Identification search1;
	Identification search2;
	
	TEST_EQUAL(search1 == search2, true)
	search1.setDateTime(date);
	TEST_EQUAL(search1 == search2, false)	
	search2.setDateTime(date);
	TEST_EQUAL(search1 == search2, true)	
	search1.setPeptideSignificanceThreshold(peptide_significance_threshold);
	TEST_EQUAL(search1 == search2, false)	
	search2.setPeptideSignificanceThreshold(peptide_significance_threshold);
	TEST_EQUAL(search1 == search2, true)	
	search1.setProteinSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(search1 == search2, false)	
	search2.setProteinSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(search1 == search2, true)	
	search1.setPeptideAndProteinHits(peptide_hits, protein_hits);
	TEST_EQUAL(search1 == search2, false)	
	search2.setPeptideAndProteinHits(peptide_hits, protein_hits);
	TEST_EQUAL(search1 == search2, true)	
RESULT

CHECK((float getPeptideSignificanceThreshold() const))
	ptr1 = new Identification();
	ptr1->setPeptideSignificanceThreshold(peptide_significance_threshold);
	TEST_EQUAL(ptr1->getPeptideSignificanceThreshold(), peptide_significance_threshold)
RESULT

CHECK((std::vector<PeptideHit>& getPeptideHits()))
	ptr1 = new Identification();
	ptr1->insertPeptideHit(peptide_hit);
	TEST_EQUAL(ptr1->getPeptideHits().size() == 1, true)
	TEST_EQUAL(*(ptr1->getPeptideHits().begin()) == peptide_hit, true)	
RESULT

CHECK((const std::vector<PeptideHit>& getPeptideHits() const))
	ptr1 = new Identification();
	ptr1->insertPeptideHit(peptide_hit);
	TEST_EQUAL(ptr1->getPeptideHits().size() == 1, true)
	PeptideHit& temp_hit = *(ptr1->getPeptideHits().begin());
	TEST_EQUAL(temp_hit == peptide_hit, true)	
RESULT


CHECK((void clear()))
	DateTime test_date;
	ptr1 = new Identification();
	ptr1->setDateTime(date);
	ptr1->setPeptideSignificanceThreshold(peptide_significance_threshold);
	ptr1->setProteinSignificanceThreshold(protein_significance_threshold);
	ptr1->setPeptideAndProteinHits(peptide_hits, protein_hits);
	ptr1->clear();
	TEST_EQUAL(ptr1->getDateTime() == test_date, true)
	TEST_EQUAL(ptr1->getPeptideSignificanceThreshold(), 0)
	TEST_EQUAL(ptr1->getProteinSignificanceThreshold(), 0)
	TEST_EQUAL(ptr1->getProteinHits().size() == 0, true)
	TEST_EQUAL(ptr1->getPeptideHits().size() == 0, true)
RESULT

CHECK((void insertPeptideHit(const PeptideHit& input)))
	ptr1 = new Identification();
	ptr1->insertPeptideHit(peptide_hit);
	TEST_EQUAL(ptr1->getPeptideHits().size() == 1, true)
	TEST_EQUAL(*(ptr1->getPeptideHits().begin()) == peptide_hit, true)	
RESULT

CHECK((void insertProteinHit(const ProteinHit& input)))
	ptr1 = new Identification();
	ptr1->insertProteinHit(protein_hit);
	TEST_EQUAL(ptr1->getProteinHits().size() == 1, true)
	TEST_EQUAL(*(ptr1->getProteinHits().begin()) == protein_hit, true)
RESULT

CHECK((void setPeptideAndProteinHits(const std::vector<PeptideHit>& peptide_hits, const std::vector<ProteinHit>& protein_hits)))
	ptr1 = new Identification();
	ptr1->setPeptideAndProteinHits(peptide_hits, protein_hits);
RESULT

CHECK((void setPeptideSignificanceThreshold(float value)))
	ptr1 = new Identification();
	ptr1->setPeptideSignificanceThreshold(peptide_significance_threshold);
	TEST_EQUAL(ptr1->getPeptideSignificanceThreshold(), peptide_significance_threshold)
RESULT

CHECK((bool empty() const))
	ptr1 = new Identification();
	TEST_EQUAL(ptr1->empty(), true)
	ptr1->setPeptideSignificanceThreshold(1);
	TEST_EQUAL(ptr1->empty(), false)
	ptr1->setPeptideSignificanceThreshold(0);
	ptr1->setProteinSignificanceThreshold(1);
	TEST_EQUAL(ptr1->empty(), false)
	ptr1->setProteinSignificanceThreshold(0);
	ptr1->insertPeptideHit(peptide_hit);
	TEST_EQUAL(ptr1->empty(), false)
	ptr1->clear();
	ptr1->insertProteinHit(protein_hit);
	TEST_EQUAL(ptr1->empty(), false)
	ptr1->clear();
	ptr1->setDateTime(date);	
	TEST_EQUAL(ptr1->empty(), false)
	ptr1->clear();
RESULT

CHECK((std::vector<PeptideHit>* getReferencingHits(String date_time, String accession) const))
	
	xml_file.load("data/MascotXMLFile_test_1.mascotXML", protein_identification, identifications);
	hits = identifications[0].id.getReferencingHits(String("2006-03-09 11:31:52"), String("AAN17824"));
	TEST_EQUAL(hits->size(), 2)
	TEST_EQUAL((*hits)[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL((*hits)[1].getSequence(), "MRSLGYVAVISAVATDTDK")
	TEST_REAL_EQUAL((*hits)[0].getScore(), 33.85)
	TEST_REAL_EQUAL((*hits)[1].getScore(), 33.12)
	hits = identifications[0].id.getReferencingHits(String("2006-03-09 11:31:52"), String("BAN17824"));
	TEST_EQUAL(hits->size(), 0)	
RESULT

CHECK((template<class iteratorT> std::vector<PeptideHit>* getNonReferencingHits(iteratorT protein_hits_begin, iteratorT protein_hits_end, const String& date_time) const))
							
	xml_file.load("data/MascotXMLFile_test_1.mascotXML", protein_identification, identifications);
							
	hits = identifications[2].id.getNonReferencingHits(protein_identification.getProteinHits().begin(),
																									protein_identification.getProteinHits().end(),
																									String("2006-03-09 11:31:52"));
	TEST_EQUAL(hits->size(), 2)
	TEST_EQUAL((*hits)[0].getSequence(), "RASNSPQDPQSATAHSFR")
	TEST_EQUAL((*hits)[1].getSequence(), "MYSTVGPA")
	TEST_REAL_EQUAL((*hits)[0].getScore(), 5.41)
	TEST_REAL_EQUAL((*hits)[1].getScore(), 7.87)
	delete hits;
	hits = identifications[0].id.getNonReferencingHits(protein_identification.getProteinHits().begin(),
																									protein_identification.getProteinHits().end(),
																									String("2006-03-09 11:31:52"));
	TEST_EQUAL(hits->size(), 0)	
	delete hits;
	
RESULT

CHECK((std::vector<PeptideHit>* getNonReferencingHits(const std::multimap< String, ProteinHit >& protein_hits) const))
	multimap< String, ProteinHit > map;
	String date_time;
																																
	xml_file.load("data/MascotXMLFile_test_1.mascotXML", protein_identification, identifications);
	vector<ProteinHit> protein_hits = protein_identification.getProteinHits();

	protein_identification.getDateTime().get(date_time);														
	for(vector<ProteinHit>::iterator it = protein_hits.begin();
			it != protein_hits.end();
			it++)
	{
		map.insert(make_pair(date_time, *it));
	}																													
	hits = identifications[2].id.getNonReferencingHits(map);
	TEST_EQUAL(hits->size(), 2)
	TEST_EQUAL((*hits)[0].getSequence(), "RASNSPQDPQSATAHSFR")
	TEST_EQUAL((*hits)[1].getSequence(), "MYSTVGPA")
	TEST_REAL_EQUAL((*hits)[0].getScore(), 5.41)
	TEST_REAL_EQUAL((*hits)[1].getScore(), 7.87)
	delete hits;
	hits = identifications[0].id.getNonReferencingHits(map);
	TEST_EQUAL(hits->size(), 0)	
	delete hits;

RESULT

CHECK((void sort()))
	vector<ProteinIdentification> protein_identifications; 
	vector<IdentificationData> identifications; 
	PeptideHit hit;
	
	hit.setSequence("TESTPEPTIDE");
	hit.setScore(33.9);
	hit.setScoreType("Mascot");

	AnalysisXMLFile().load("data/AnalysisXMLFile_test.analysisXML", protein_identifications, identifications);
	TEST_EQUAL(identifications.size(), 3)

	identifications[0].id.insertPeptideHit(hit);
	identifications[0].id.sort();
	TEST_EQUAL(identifications[0].id.getPeptideHits().size(), 3)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getSequence(), "TESTPEPTIDE")

RESULT

CHECK((void assignRanks()))
	vector<ProteinIdentification> protein_identifications; 
	vector<IdentificationData> identifications;
	PeptideHit hit;
	
	hit.setSequence("TESTPEPTIDE");
	hit.setScore(33.9);
	hit.setScoreType("Mascot");

	AnalysisXMLFile().load("data/AnalysisXMLFile_test.analysisXML", protein_identifications, identifications);
	TEST_EQUAL(identifications.size(), 3)

	identifications[0].id.insertPeptideHit(hit);
	identifications[0].id.assignRanks();
	TEST_EQUAL(identifications[0].id.getPeptideHits().size(), 3)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getSequence(), "TESTPEPTIDE")
	TEST_EQUAL(identifications[0].id.getPeptideHits()[0].getRank(), 1)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[1].getRank(), 2)
	TEST_EQUAL(identifications[0].id.getPeptideHits()[2].getRank(), 3)

RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
