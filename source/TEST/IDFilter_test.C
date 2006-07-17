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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <string>

#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/FORMAT/MascotOutfile.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>

#include <vector>

///////////////////////////

START_TEST(IDFilter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

IDFilter* ptr1;
IDFilter* ptr2;
Identification identification;

AnalysisXMLFile xml_file;
std::vector< ProteinIdentification > protein_identifications;
std::vector< Identification > identifications;	
std::vector< Real > precursor_retention_times;
std::vector< Real > precursor_mz_values;
ContactPerson person;	
	
xml_file.load("data/IDFilter_test.analysisXML", 
							&protein_identifications, 
							&identifications, 
							&precursor_retention_times, 
							&precursor_mz_values,
							&person);

identification = identifications[0];							
	
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

CHECK(const Identification& filterIdentificationsByProteins(const Identification& identification))
	Identification identification2;
  vector< pair<String, String> > proteins;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

  proteins.push_back(pair<String, String>("Q824A5", "LHASGITVTEIPVTATNFK"));
  proteins.push_back(pair<String, String>("Q872T5", "THPYGHAIVAGIERYPSK"));
	ptr1 = new IDFilter();

	ptr1->setProteins(proteins);
	identification2 = ptr1->filterIdentificationsByProteins(identification);
	peptide_hits = identification2.getPeptideHits();
	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_EQUAL(peptide_hits[0].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[0].getScoreType() , "Mascot")
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL(peptide_hits[1].getSequence() , "THPYGHAIVAGIERYPSK")
	TEST_REAL_EQUAL(peptide_hits[1].getScore() , 11.26)
	TEST_EQUAL(peptide_hits[1].getScoreType() , "Mascot")
	TEST_EQUAL(peptide_hits[1].getRank() , 2)	
	protein_hits = identification2.getProteinHits();
	TEST_EQUAL(protein_hits.size(), 2)
	TEST_EQUAL(protein_hits[0].getAccession(), "Q824A5")
	TEST_EQUAL(protein_hits[1].getAccession(), "Q872T5")
RESULT

CHECK((const Identification& filterIdentificationsByProteins(const Identification& identification, std::vector< std::pair<String, String> > proteins)))
	Identification identification2;
  vector< pair<String, String> > proteins;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

  proteins.push_back(pair<String, String>("Q824A5", "LHASGITVTEIPVTATNFK"));
  proteins.push_back(pair<String, String>("Q872T5", "THPYGHAIVAGIERYPSK"));
	ptr1 = new IDFilter();

	identification2 = ptr1->filterIdentificationsByProteins(identification, proteins);
	peptide_hits = identification2.getPeptideHits();
	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_EQUAL(peptide_hits[0].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[0].getScoreType() , "Mascot")
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL(peptide_hits[1].getSequence() , "THPYGHAIVAGIERYPSK")
	TEST_REAL_EQUAL(peptide_hits[1].getScore() , 11.26)
	TEST_EQUAL(peptide_hits[1].getScoreType() , "Mascot")
	TEST_EQUAL(peptide_hits[1].getRank() , 2)	
	protein_hits = identification2.getProteinHits();
	TEST_EQUAL(protein_hits.size(), 2)
	TEST_EQUAL(protein_hits[0].getAccession(), "Q824A5")
	TEST_EQUAL(protein_hits[1].getAccession(), "Q872T5")
    
RESULT

CHECK((const Identification& filterIdentificationsByThresholds(const Identification& identification, bool strict = false)))
	Identification identification2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

	ptr1 = new IDFilter();
	identification2 = ptr1->filterIdentificationsByThresholds(identification);
	peptide_hits = identification2.getPeptideHits();
	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_EQUAL(peptide_hits[0].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[0].getScoreType() , "Mascot")
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL(peptide_hits[1].getSequence() , "MRSLGYVAVISAVATDTDK")
	TEST_REAL_EQUAL(peptide_hits[1].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[1].getRank() , 2)
	TEST_EQUAL(peptide_hits[1].getScoreType() , "Mascot")
	protein_hits = identification2.getProteinHits();
	TEST_EQUAL(protein_hits.size(), 1)
	TEST_EQUAL(protein_hits[0].getAccession(), "Q824A5")
	
RESULT

CHECK((const Identification& filterIdentificationsByThresholds(const Identification& identification, const double& peptide_threshold_fraction, const double& protein_threshold_fraction, bool strict = false)))
	Identification identification2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

	ptr1 = new IDFilter();
	identification2 = ptr1->filterIdentificationsByThresholds(identification, 1.3, 1);
	peptide_hits = identification2.getPeptideHits();
	TEST_EQUAL(peptide_hits.size(), 0)
	identification2 = ptr1->filterIdentificationsByThresholds(identification, 1, 1);
	peptide_hits = identification2.getPeptideHits();
	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_EQUAL(peptide_hits[0].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[0].getScoreType() , "Mascot")
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL(peptide_hits[1].getSequence() , "MRSLGYVAVISAVATDTDK")
	TEST_REAL_EQUAL(peptide_hits[1].getScore() , 33.85)
	TEST_EQUAL(peptide_hits[1].getRank() , 2)
	TEST_EQUAL(peptide_hits[1].getScoreType() , "Mascot")
	protein_hits = identification2.getProteinHits();
	TEST_EQUAL(protein_hits.size(), 1)
	TEST_EQUAL(protein_hits[0].getAccession(), "Q824A5")

RESULT

CHECK((const Identification& filterIdentificationsByExclusionPeptides(const Identification& identification, std::vector<String> peptides)))
	Identification identification2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;
	vector<String> peptides;
	
	peptides.push_back("LHASGITVTEIPVTATNFK");
	peptides.push_back("MRSLGYVAVISAVATDTDK");
	peptides.push_back("EGASTDFAALRTFLAEDGK");
	peptides.push_back("DLEPGTDYEVTVSTLFGR");
	peptides.push_back("FINFGVNVEVLSRFQTK");
	peptides.push_back("MSLLSNMISIVKVGYNAR");
	peptides.push_back("THPYGHAIVAGIERYPSK");
	peptides.push_back("AITSDFANQAKTVLQNFK");

// TGCDTWGQGTLVTVSSASTK
// TLCHHDATFDNLVWTPK

	ptr1 = new IDFilter();
	identification2 = ptr1->filterIdentificationsByExclusionPeptides(identification, peptides);
	peptide_hits = identification2.getPeptideHits();
	TEST_EQUAL(peptide_hits.size(), 2)
	TEST_EQUAL(peptide_hits[0].getSequence() , "TGCDTWGQGTLVTVSSASTK")
	TEST_REAL_EQUAL(peptide_hits[0].getScore() , 10.93)
	TEST_EQUAL(peptide_hits[0].getScoreType() , "Mascot")
	TEST_EQUAL(peptide_hits[0].getRank() , 1)
	TEST_EQUAL(peptide_hits[1].getSequence() , "TLCHHDATFDNLVWTPK")
	TEST_REAL_EQUAL(peptide_hits[1].getScore() , 10.37)
	TEST_EQUAL(peptide_hits[1].getRank() , 2)
	TEST_EQUAL(peptide_hits[1].getScoreType() , "Mascot")
	protein_hits = identification2.getProteinHits();
	TEST_EQUAL(protein_hits.size(), 50)

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

CHECK((const std::vector< std::pair<String, String>& getProteins() const))
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

CHECK((void setProteins(const std::vector< std::pair<String, String> >& proteins)))
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
