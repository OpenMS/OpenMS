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
// $Id: IDFilter_test.C,v 1.4 2006/06/09 23:47:35 nicopfeifer Exp $
// $Author: nicopfeifer $
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <string>

#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/FORMAT/MascotOutfile.h>

#include <vector>

///////////////////////////

START_TEST(IDFilter, "$Id: IDFilter_test.C,v 1.4 2006/06/09 23:47:35 nicopfeifer Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

IDFilter* ptr1;
IDFilter* ptr2;
Identification db_search;
MascotOutfile file("data/MascotOutfile.txt");
file >> db_search;

CHECK(IDFilter& operator = (const IDFilter& source))
	ptr1 = new IDFilter();
	ptr2 = new IDFilter();
  vector< pair<String, String> > proteins;
  
  proteins.push_back(pair<String, String>("Q824A5", "LHASGITVTEIPVTATNFK"));
  proteins.push_back(pair<String, String>("Q872T5", "THPYGHAIVAGIERYPSK"));
  ptr1->setProteinThresholdFraction(0.4);
  ptr1->setPeptideThresholdFraction(0.6);
  ptr1->setProteins(proteins);
  *ptr2 = *ptr1;
  TEST_EQUAL(ptr1->getProteinThresholdFraction(), ptr2->getProteinThresholdFraction());
  TEST_EQUAL(ptr1->getPeptideThresholdFraction(), ptr2->getPeptideThresholdFraction());
  TEST_EQUAL(ptr1->getProteins() == ptr2->getProteins(), true)
RESULT

CHECK(IDFilter())
	ptr1 = new IDFilter();
	TEST_NOT_EQUAL(ptr1, 0);
RESULT

CHECK(IDFilter(const IDFilter& source))
	ptr1 = new IDFilter();
  vector< pair<String, String> > proteins;
  proteins.push_back(pair<String, String>("Q824A5", "LHASGITVTEIPVTATNFK"));
  proteins.push_back(pair<String, String>("Q872T5", "THPYGHAIVAGIERYPSK"));
  ptr1->setProteinThresholdFraction(0.4);
  ptr1->setPeptideThresholdFraction(0.6);
  ptr1->setProteins(proteins);
  ptr2 = new IDFilter(*ptr1);
  TEST_EQUAL(ptr1->getProteinThresholdFraction(), ptr2->getProteinThresholdFraction());
  TEST_EQUAL(ptr1->getPeptideThresholdFraction(), ptr2->getPeptideThresholdFraction());
  TEST_EQUAL(ptr1->getProteins() == ptr2->getProteins(), true)
RESULT

CHECK(const Identification& filterIdentificationsByProteins(const Identification& db_search))
	Identification db_search2;
  vector< pair<String, String> > proteins;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

  proteins.push_back(pair<String, String>("Q824A5", "LHASGITVTEIPVTATNFK"));
  proteins.push_back(pair<String, String>("Q872T5", "THPYGHAIVAGIERYPSK"));
	ptr1 = new IDFilter();

	ptr1->setProteins(proteins);
	db_search2 = ptr1->filterIdentificationsByProteins(db_search);
	peptide_hits = db_search2.getPeptideHits();
	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_EQUAL(peptide_hits[0].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[0].getScoreType() , "Mascot")
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL(peptide_hits[1].getSequence() , "THPYGHAIVAGIERYPSK")
	TEST_REAL_EQUAL(peptide_hits[1].getScore() , 11.26)
	TEST_EQUAL(peptide_hits[1].getScoreType() , "Mascot")
	TEST_EQUAL(peptide_hits[1].getRank() , 2)	
	protein_hits = db_search2.getProteinHits();
	TEST_EQUAL(protein_hits.size(), 2)
	TEST_EQUAL(protein_hits[0].getAccession(), "Q824A5")
	TEST_EQUAL(protein_hits[1].getAccession(), "Q872T5")
RESULT

CHECK(const Identification& filterIdentificationsByProteins(const Identification& db_search, vector< pair<String, String>> proteins))
	Identification db_search2;
  vector< pair<String, String> > proteins;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

  proteins.push_back(pair<String, String>("Q824A5", "LHASGITVTEIPVTATNFK"));
  proteins.push_back(pair<String, String>("Q872T5", "THPYGHAIVAGIERYPSK"));
	ptr1 = new IDFilter();

	db_search2 = ptr1->filterIdentificationsByProteins(db_search, proteins);
	peptide_hits = db_search2.getPeptideHits();
	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_EQUAL(peptide_hits[0].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[0].getScoreType() , "Mascot")
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL(peptide_hits[1].getSequence() , "THPYGHAIVAGIERYPSK")
	TEST_REAL_EQUAL(peptide_hits[1].getScore() , 11.26)
	TEST_EQUAL(peptide_hits[1].getScoreType() , "Mascot")
	TEST_EQUAL(peptide_hits[1].getRank() , 2)	
	protein_hits = db_search2.getProteinHits();
	TEST_EQUAL(protein_hits.size(), 2)
	TEST_EQUAL(protein_hits[0].getAccession(), "Q824A5")
	TEST_EQUAL(protein_hits[1].getAccession(), "Q872T5")
    
RESULT

CHECK(const Identification& filterIdentificationsByThresholds(const Identification& db_search))
	Identification db_search2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

	ptr1 = new IDFilter();
	db_search2 = ptr1->filterIdentificationsByThresholds(db_search);
	peptide_hits = db_search2.getPeptideHits();
	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_EQUAL(peptide_hits[0].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[0].getScoreType() , "Mascot")
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL(peptide_hits[1].getSequence() , "MRSLGYVAVISAVATDTDK")
	TEST_REAL_EQUAL(peptide_hits[1].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[1].getRank() , 2)
	TEST_EQUAL(peptide_hits[1].getScoreType() , "Mascot")
	protein_hits = db_search2.getProteinHits();
	TEST_EQUAL(protein_hits.size(), 1)
	TEST_EQUAL(protein_hits[0].getAccession(), "Q824A5")
	
RESULT

CHECK((const Identification& filterIdentificationsByThresholds(const Identification& db_search, const double& peptide_threshold_fraction, const double& protein_threshold_fraction)))
	Identification db_search2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

	ptr1 = new IDFilter();
	db_search2 = ptr1->filterIdentificationsByThresholds(db_search, 1.3, 1);
	peptide_hits = db_search2.getPeptideHits();
	TEST_EQUAL(peptide_hits.size(), 0)
	db_search2 = ptr1->filterIdentificationsByThresholds(db_search, 1, 1);
	peptide_hits = db_search2.getPeptideHits();
	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_EQUAL(peptide_hits[0].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[0].getScoreType() , "Mascot")
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL(peptide_hits[1].getSequence() , "MRSLGYVAVISAVATDTDK")
	TEST_REAL_EQUAL(peptide_hits[1].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[1].getRank() , 2)
	TEST_EQUAL(peptide_hits[1].getScoreType() , "Mascot")
	protein_hits = db_search2.getProteinHits();
	TEST_EQUAL(protein_hits.size(), 1)
	TEST_EQUAL(protein_hits[0].getAccession(), "Q824A5")

RESULT

CHECK(const double& getPeptideThresholdFraction() const)
	ptr1 = new IDFilter();
  
  ptr1->setPeptideThresholdFraction(0.6);
  TEST_EQUAL(ptr1->getPeptideThresholdFraction(), 0.6);
RESULT

CHECK(const double& getProteinThresholdFraction() const)
	ptr1 = new IDFilter();
  
  ptr1->setProteinThresholdFraction(0.6);
  TEST_EQUAL(ptr1->getProteinThresholdFraction(), 0.6);
RESULT

CHECK((const vector< pair<String, String>& getProteins() const))
  vector< pair<String, String> > proteins;

  proteins.push_back(pair<String, String>("Q824A5", "LHASGITVTEIPVTATNFK"));
  proteins.push_back(pair<String, String>("Q872T5", "THPYGHAIVAGIERYPSK"));
	ptr1 = new IDFilter();
	ptr1->setProteins(proteins);
	TEST_EQUAL(ptr1->getProteins() == proteins, true)
RESULT

CHECK(void setPeptideThresholdFraction(const double& peptide_threshold_fraction))
	ptr1 = new IDFilter();
  
  ptr1->setPeptideThresholdFraction(0.6);
  TEST_EQUAL(ptr1->getPeptideThresholdFraction(), 0.6);
RESULT

CHECK(void setProteinThresholdFraction(const double& protein_threshold_fraction))
	ptr1 = new IDFilter();
  
  ptr1->setProteinThresholdFraction(0.6);
  TEST_EQUAL(ptr1->getProteinThresholdFraction(), 0.6);
RESULT

CHECK(void setProteins(const vector<String>& proteins))
  vector< pair<String, String> > proteins;

  proteins.push_back(pair<String, String>("Q824A5", "LHASGITVTEIPVTATNFK"));
  proteins.push_back(pair<String, String>("Q872T5", "THPYGHAIVAGIERYPSK"));
	ptr1 = new IDFilter();
	ptr1->setProteins(proteins);
	TEST_EQUAL(ptr1->getProteins() == proteins, true)
RESULT

CHECK(~IDFilter())
	ptr1 = new IDFilter();
	delete ptr1;	
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
