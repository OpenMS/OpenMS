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
// $Id: MascotOutfile_test.C,v 1.16 2006/06/10 06:40:18 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

///////////////////////////

#include <OpenMS/FORMAT/MascotOutfile.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(String, "$Id: MascotOutfile_test.C,v 1.16 2006/06/10 06:40:18 marc_sturm Exp $")

/////////////////////////////////////////////////////////////

DateTime date;
date.set("27.01.2005 17:47:41");
MascotOutfile* ptr = 0;
CHECK(MascotOutfile(filename))
	ptr = new MascotOutfile("data/MascotOutfile.txt");
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~MascotOutfile())
	delete ptr;
RESULT

MascotOutfile file("data/MascotOutfile.txt");

CHECK( ok )
	TEST_EQUAL(file.ok(),true)
RESULT

CHECK( ok )
	TEST_EXCEPTION(Exception::ParseError,MascotOutfile file("data/TextFile_test_empty_infile.txt"))
RESULT

CHECK( operator >> PeptideHit )
	PeptideHit hit;

	file >> hit;
	TEST_EQUAL(hit.getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_EQUAL(hit.getScore() , 33.85)
	TEST_EQUAL(hit.getScoreType() , "Mascot")
	TEST_EQUAL(hit.getRank() , 1)

	file >> hit;
	TEST_EQUAL(hit.getSequence() , "MRSLGYVAVISAVATDTDK")
	TEST_REAL_EQUAL(hit.getScore() , 33.85)
	TEST_EQUAL(hit.getRank() , 2)

	file >> hit;
	TEST_EQUAL(hit.getSequence() , "EGASTDFAALRTFLAEDGK")
	TEST_REAL_EQUAL(hit.getScore() , 12.73)
	TEST_EQUAL(hit.getRank() , 3)

	file >> hit;
	TEST_EQUAL(hit.getSequence() , "DLEPGTDYEVTVSTLFGR")
	TEST_REAL_EQUAL(hit.getScore() , 11.79)
	TEST_EQUAL(hit.getRank() , 4)

	file >> hit;
	TEST_EQUAL(hit.getSequence() , "FINFGVNVEVLSRFQTK")
	TEST_REAL_EQUAL(hit.getScore() , 11.38)
	TEST_EQUAL(hit.getRank() , 5)
RESULT

CHECK( operator >> Identification )
	Identification search;
	file >> search;
	
	TEST_EQUAL(search.getDateTime() == date, true)
	TEST_EQUAL(search.getPeptideHits().size(),10)
	TEST_EQUAL(search.getProteinHits().size(),50)
	TEST_EQUAL(search.getCharge(),2)
RESULT

CHECK(MascotOutfile& operator>>(ProteinHit& protein_hit))
	ProteinHit hit;

	file >> hit;
	TEST_EQUAL(hit.getAccession() , "Q824A5")
	TEST_EQUAL(hit.getAccessionType() , "SwissProt")
//	TEST_EQUAL(hit.getSequence() , "LHASGITVTEIPVTATNFK")
	TEST_REAL_EQUAL(hit.getScore() , 32.3)
	TEST_EQUAL(hit.getScoreType() , "Mascot")
	TEST_EQUAL(hit.getRank() , 1)
	file >> hit;
	TEST_EQUAL(hit.getAccession() , "AAD30739")
	TEST_EQUAL(hit.getAccessionType() , "SwissProt")
//	TEST_EQUAL(hit.getSequence() , "TGCDTWGQGTLVTVSSASTK")
	TEST_REAL_EQUAL(hit.getScore() , 28.1)
	TEST_EQUAL(hit.getScoreType() , "Mascot")
	TEST_EQUAL(hit.getRank() , 2)

RESULT

CHECK(MascotOutfile& operator=(const MascotOutfile& source))
	MascotOutfile file2 = file;
	Identification db_search;
	std::vector<PeptideHit> peptide_hits;
	std::vector<ProteinHit> protein_hits;
	std::vector<uint> temp_peptide_indices;
	std::vector<uint> temp_protein_indices;
	PeptideHit temp_peptide_hit;
	ProteinHit temp_protein_hit;
		
	file2 >> db_search;

	protein_hits = db_search.getProteinHits();
	peptide_hits  = db_search.getPeptideHits();
	
	TEST_REAL_EQUAL(db_search.getPeptideSignificanceThreshold(), 31.862070)
	TEST_EQUAL(db_search.getDateTime() == date, true)
	TEST_REAL_EQUAL(db_search.getProteinSignificanceThreshold(), 0)	
	TEST_REAL_EQUAL(file2.getPrecursorRetentionTimes()[0], 25.379)	
RESULT

CHECK(MascotOutfile(const MascotOutfile& mascotoutfile))
	MascotOutfile file2 = MascotOutfile(file);
	Identification db_search;
	std::vector<PeptideHit> peptide_hits;
	std::vector<ProteinHit> protein_hits;
	std::vector<uint> temp_peptide_indices;
	std::vector<uint> temp_protein_indices;
	PeptideHit temp_peptide_hit;
	ProteinHit temp_protein_hit;
		
	file2 >> db_search;

	protein_hits = db_search.getProteinHits();
	peptide_hits  = db_search.getPeptideHits();
		
	TEST_REAL_EQUAL(db_search.getPeptideSignificanceThreshold(), 31.862070)
	TEST_EQUAL(db_search.getDateTime() == date, true)
	TEST_REAL_EQUAL(db_search.getProteinSignificanceThreshold(), 0)	
	TEST_REAL_EQUAL(file2.getPrecursorRetentionTimes()[0], 25.379)	
RESULT

CHECK(MascotOutfile(const std::string& filename) throw(Exception::ParseError))
	
	Identification db_search;
	std::vector<PeptideHit> peptide_hits;
	std::vector<ProteinHit> protein_hits;
	std::vector<uint> temp_peptide_indices;
	std::vector<uint> temp_protein_indices;
	PeptideHit temp_peptide_hit;
	ProteinHit temp_protein_hit;
		
	file >> db_search;

	protein_hits = db_search.getProteinHits();
	peptide_hits  = db_search.getPeptideHits();

	TEST_REAL_EQUAL(db_search.getPeptideSignificanceThreshold(), 31.862070)
	TEST_EQUAL(db_search.getDateTime() == date, true)
	TEST_REAL_EQUAL(db_search.getProteinSignificanceThreshold(), 0)
	TEST_REAL_EQUAL(file.getPrecursorRetentionTimes()[0], 25.379)	
RESULT

CHECK(bool ok() const)
  TEST_EQUAL(file.ok(), true)
RESULT

CHECK(const vector<float>& getPrecursorRetentionTimes() const)
	TEST_REAL_EQUAL(file.getPrecursorRetentionTimes()[0], 25.379)	
	ptr = new MascotOutfile("data/MascotOutfile2.txt");
  TEST_REAL_EQUAL(ptr->getPrecursorRetentionTimes().size(), 4)	
  TEST_REAL_EQUAL(ptr->getPrecursorRetentionTimes()[0], 88.3466f)
  TEST_REAL_EQUAL(ptr->getPrecursorRetentionTimes()[1], 96.9993f)
  TEST_REAL_EQUAL(ptr->getPrecursorRetentionTimes()[2], 105.615f)
  TEST_REAL_EQUAL(ptr->getPrecursorRetentionTimes()[3], 105.615f)
RESULT

CHECK(void setPrecursorRetentionTimes(const float& retention_time))
	vector<float> tmp;
	tmp.push_back(4.3f);
	file.setPrecursorRetentionTimes(tmp);
	TEST_REAL_EQUAL(file.getPrecursorRetentionTimes()[0], 4.3)	
RESULT

CHECK(const double& getPrecursorMZValues() const)
  TEST_REAL_EQUAL(file.getPrecursorMZValues()[0], 999.503912)
	ptr = new MascotOutfile("data/MascotOutfile2.txt");
	TEST_REAL_EQUAL(ptr->getPrecursorMZValues()[0], 508.119f)	
	TEST_REAL_EQUAL(ptr->getPrecursorMZValues()[1], 508.458f)	
	TEST_REAL_EQUAL(ptr->getPrecursorMZValues()[2], 517.267f)	
	TEST_REAL_EQUAL(ptr->getPrecursorMZValues()[3], 517.324f)	
RESULT

CHECK(void setPrecursorMZValues(const float& mz))
	vector<float> tmp;
	tmp.push_back(354.2f);
	file.setPrecursorMZValues(tmp);
  TEST_REAL_EQUAL(file.getPrecursorMZValues()[0], 354.2f)
RESULT

CHECK(void setIdentificationes(const std::vector<Identification>& db_searches))
	ptr = new MascotOutfile("data/MascotOutfile2.txt");
	vector<Identification> db_searches = ptr->getIdentifications();	
RESULT

CHECK(const std::vector<Identification>& getIdentificationes() const)
	ptr = new MascotOutfile("data/MascotOutfile2.txt");
	vector<Identification> db_searches = ptr->getIdentifications();
	TEST_EQUAL(db_searches.size(), 4)	
	TEST_EQUAL(db_searches[0].getCharge(), -2)	
	TEST_EQUAL(db_searches[1].getCharge(), 3)	
	TEST_EQUAL(db_searches[2].getCharge(), 3)	
	TEST_EQUAL(db_searches[3].getCharge(), 3)	
	TEST_EQUAL(db_searches[0].getPeptideHits().size(), 1)	
	TEST_EQUAL(db_searches[1].getPeptideHits().size(), 1)	
	TEST_EQUAL(db_searches[2].getPeptideHits().size(), 10)	
	TEST_EQUAL(db_searches[3].getPeptideHits().size(), 10)	
	TEST_REAL_EQUAL(db_searches[0].getPeptideHits()[0].getScore(), 19.1f)	
	TEST_EQUAL(db_searches[0].getPeptideHits()[0].getSequence(), "NSSEA")	
	TEST_EQUAL(db_searches[0].getPeptideHits()[0].getScoreType(), "Mascot")	
	TEST_EQUAL(db_searches[0].getPeptideHits()[0].getRank(), 1)	
	TEST_REAL_EQUAL(db_searches[1].getPeptideHits()[0].getScore(), 0.93f)	
	TEST_EQUAL(db_searches[1].getPeptideHits()[0].getSequence(), "FGASK")	
	TEST_EQUAL(db_searches[1].getPeptideHits()[0].getScoreType(), "Mascot")	
	TEST_EQUAL(db_searches[1].getPeptideHits()[0].getRank(), 1)	
	TEST_REAL_EQUAL(db_searches[2].getPeptideHits()[0].getScore(), 9.72f)	
	TEST_EQUAL(db_searches[2].getPeptideHits()[0].getSequence(), "AGGNAK")	
	TEST_EQUAL(db_searches[2].getPeptideHits()[0].getScoreType(), "Mascot")	
	TEST_EQUAL(db_searches[2].getPeptideHits()[0].getRank(), 1)	
	TEST_REAL_EQUAL(db_searches[2].getPeptideHits()[1].getScore(), 8.77f)	
	TEST_EQUAL(db_searches[2].getPeptideHits()[1].getSequence(), "KGANK")	
	TEST_EQUAL(db_searches[2].getPeptideHits()[1].getScoreType(), "Mascot")	
	TEST_EQUAL(db_searches[2].getPeptideHits()[1].getRank(), 2)	
	TEST_REAL_EQUAL(db_searches[2].getPeptideHits()[2].getScore(), 8.77f)	
	TEST_EQUAL(db_searches[2].getPeptideHits()[2].getSequence(), "KXANK")	
	TEST_EQUAL(db_searches[2].getPeptideHits()[2].getScoreType(), "Mascot")	
	TEST_EQUAL(db_searches[2].getPeptideHits()[2].getRank(), 3)	

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
