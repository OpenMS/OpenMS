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
ProteinIdentification protein_identification;

std::vector< ProteinIdentification > protein_identifications;
std::vector< IdentificationData > identifications;
std::vector< IdentificationData > identifications2;
IdentificationData data;
	
AnalysisXMLFile().load("data/IDFilter_test.analysisXML", 
							protein_identifications, 
							identifications);

identifications2.push_back(data);
identifications2[0].rt = identifications[0].rt;
identifications2[0].mz = identifications[0].mz;

identification = identifications[0].id;							
protein_identification = protein_identifications[0];							
	
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

CHECK(void filterIdentificationsByProteins(const Identification& identification, Identification& filtered_identification))
	Identification identification2;
	ProteinIdentification protein_identification2;
  vector< pair<String, String> > proteins;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

  proteins.push_back(pair<String, String>("Q824A5", "LHASGITVTEIPVTATNFK"));
  proteins.push_back(pair<String, String>("Q872T5", "THPYGHAIVAGIERYPSK"));
	ptr1 = new IDFilter();

	ptr1->setProteins(proteins);
	ptr1->filterIdentificationsByProteins(identification, identification2);
	ptr1->filterIdentificationsByProteins(protein_identification, protein_identification2);
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
	protein_hits = protein_identification2.getProteinHits();
	TEST_EQUAL(protein_hits.size(), 2)
	TEST_EQUAL(protein_hits[0].getAccession(), "Q824A5")
	TEST_EQUAL(protein_hits[1].getAccession(), "Q872T5")

RESULT

CHECK((void filterIdentificationsByProteins(const Identification& identification, std::vector< std::pair<String, String> > proteins, Identification& filtered_identification)))
	Identification identification2;
	ProteinIdentification protein_identification2;
  vector< pair<String, String> > proteins;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

  proteins.push_back(pair<String, String>("Q824A5", "LHASGITVTEIPVTATNFK"));
  proteins.push_back(pair<String, String>("Q872T5", "THPYGHAIVAGIERYPSK"));
	ptr1 = new IDFilter();

	ptr1->filterIdentificationsByProteins(identification, proteins, identification2);
	ptr1->filterIdentificationsByProteins(protein_identification, protein_identification2);
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
	protein_hits = protein_identification2.getProteinHits();
	TEST_EQUAL(protein_hits.size(), 2)
	TEST_EQUAL(protein_hits[0].getAccession(), "Q824A5")
	TEST_EQUAL(protein_hits[1].getAccession(), "Q872T5")
    
RESULT

CHECK((void filterIdentificationsByThresholds(const Identification& identification, Identification& filtered_identification)))
	Identification identification2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

	ptr1 = new IDFilter();
	ptr1->filterIdentificationsByThresholds(identification, identification2);
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
	
RESULT

CHECK((void filterIdentificationsByThresholds(const Identification& identification, const double& peptide_threshold_fraction, const double& protein_threshold_fraction, Identification& filtered_identification)))
	Identification identification2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;

	ptr1 = new IDFilter();
	TEST_EQUAL(identification.getPeptideHits().size(), 10)	
	ptr1->filterIdentificationsByThresholds(identification, 1.3, 1, identification2);
	peptide_hits = identification2.getPeptideHits();
	TEST_EQUAL(peptide_hits.size(), 0)
	ptr1->filterIdentificationsByThresholds(identification, 1, 1, identification2);
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
	protein_hits = protein_identification.getProteinHits();
	TEST_EQUAL(protein_hits.size(), 50)
	TEST_EQUAL(protein_hits[0].getAccession(), "Q824A5")

RESULT

CHECK((void filterIdentificationsByExclusionPeptides(const Identification& identification, std::vector<String> peptides, Identification& filtered_identification)))
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
	ptr1->filterIdentificationsByExclusionPeptides(identification, peptides, identification2);
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
	protein_hits = protein_identification.getProteinHits();
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

CHECK((template<class PeakT> void filterIdentificationsByProteins(MSExperiment< PeakT >& experiment, std::vector< std::pair<String, String> >proteins)))
	MSExperiment< DPeak<1> > experiment;
  vector< pair<String, String> > proteins;
  vector< Identification > ids;
	Identification identification2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;
  
  ids.push_back(identification);

  proteins.push_back(pair<String, String>("Q824A5", "LHASGITVTEIPVTATNFK"));
  proteins.push_back(pair<String, String>("Q872T5", "THPYGHAIVAGIERYPSK"));
	
	for(UnsignedInt i = 0; i < 5; ++i)
	{
		experiment.push_back(MSSpectrum< DPeak<1> >());
	}
	experiment[3].setMSLevel(2);
	experiment[3].setIdentifications(ids);
	
	ptr1->filterIdentificationsByProteins(experiment, proteins);

	identification2 = experiment[3].getIdentifications()[0];
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
	
RESULT

CHECK((template<class PeakT> void filterIdentificationsByThresholds(MSExperiment< PeakT >& experiment, double peptide_threshold_fraction, double protein_threshold_fraction)))
	MSExperiment< DPeak<1> > experiment;
  vector< Identification > ids;
	Identification identification2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;
  
  ids.push_back(identification);
	
	for(UnsignedInt i = 0; i < 5; ++i)
	{
		experiment.push_back(MSSpectrum< DPeak<1> >());
	}
	experiment[3].setMSLevel(2);
	experiment[3].setIdentifications(ids);

	ptr1->filterIdentificationsByThresholds(experiment, 1, 1);
	identification2 = experiment[3].getIdentifications()[0];
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
	protein_hits = protein_identification.getProteinHits();
	TEST_EQUAL(protein_hits.size(), 50)
	TEST_EQUAL(protein_hits[0].getAccession(), "Q824A5")
	ptr1->filterIdentificationsByThresholds(experiment, 1.3, 1);
	identification2 = experiment[3].getIdentifications()[0];
	peptide_hits = identification2.getPeptideHits();
	TEST_EQUAL(peptide_hits.size(), 0)

RESULT

CHECK((void filterIdentificationsByBestHits(const Identification& identification, Identification& filtered_identification, bool strict = false)))
	Identification identification2;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;
	bool strict = true;
	
	ptr1 = new IDFilter();

	ptr1->filterIdentificationsByBestHits(identification, identification2, strict);
	peptide_hits = identification2.getPeptideHits();
	TEST_EQUAL(peptide_hits.size(), 0)
	ptr1->filterIdentificationsByBestHits(identification, identification2);
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

	delete ptr1;		
RESULT

CHECK((void filterIdentificationsByRetentionTimes(const Identification& identification, const std::map<String, double>& predicted_retention_times, double measured_retention_time, double predicted_sigma, double allowed_deviation, double total_gradient_time, Identification& filtered_identification)))

RESULT

CHECK((void filterIdentificationsByProteins(const ProteinIdentification& identification, ProteinIdentification& filtered_identification)))
  vector< pair<String, String> > proteins;
	vector<PeptideHit> peptide_hits;
	vector<ProteinHit> protein_hits;
	ProteinIdentification id;
	ProteinIdentification filtered_id;
	IDFilter filter;
	
  proteins.push_back(pair<String, String>("Q824A5", "LHASGITVTEIPVTATNFK"));
  proteins.push_back(pair<String, String>("Q872T5", "THPYGHAIVAGIERYPSK"));

	protein_hits.push_back(ProteinHit(22, "Mascot", 1, "Q824A5", "MSDB", "LHASGITVTEIPVTATNFK"));
	protein_hits.push_back(ProteinHit(6, "Mascot", 5, "Q872T5", "MSDB", "THPYGHAIVAGIERYPSK"));
	protein_hits.push_back(ProteinHit(16, "Mascot", 3, "Q164A5", "MSDB", "LHASGIRYKKAATNFK"));
	protein_hits.push_back(ProteinHit(10, "Mascot", 4, "Q133A5", "MSDB", "LHASGGTRAYKPVTATNFK"));
	protein_hits.push_back(ProteinHit(19, "Mascot", 2, "Q255A5", "MSDB", "LHYRTKLLIVTATNFK"));
	protein_hits.push_back(ProteinHit(2, "Mascot", 6, "Q783A5", "MSDB", "LHAAELIIVTATNFK"));
	
	id.setProteinHits(protein_hits);
	filter.setProteins(proteins);
	filter.filterIdentificationsByProteins(id, filtered_id);
	protein_hits.clear();
	protein_hits = filtered_id.getProteinHits();
	TEST_EQUAL(protein_hits.size(), 2) 
	TEST_EQUAL(protein_hits[0].getAccession(), "Q824A5")
	TEST_EQUAL(protein_hits[1].getAccession(), "Q872T5")
	TEST_EQUAL(protein_hits[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(protein_hits[1].getSequence(), "THPYGHAIVAGIERYPSK")
RESULT

CHECK(~IDFilter())
	ptr1 = new IDFilter();
	delete ptr1;	
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
