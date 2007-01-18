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

#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/DATASTRUCTURES/String.h>

///////////////////////////

START_TEST(ProteinHit, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ProteinHit* ptr1 = 0;
Real score = 4.4;
std::string score_type = "XCorr";
UnsignedInt rank = 3;
String sequence = "ARRAY";
String accession = "PROOE34";
std::string accession_type = "SWISSPROT";
	
CHECK(ProteinHit())
	ptr1 = new ProteinHit();
	TEST_NOT_EQUAL(ptr1, 0)
	TEST_EQUAL(ptr1->getScore(), 0)
	TEST_EQUAL(ptr1->getScoreType(), "")
	TEST_EQUAL(ptr1->getRank(), 0)
	TEST_EQUAL(ptr1->getAccession(), "")
	TEST_EQUAL(ptr1->getAccessionType(), "")
	TEST_EQUAL(ptr1->getSequence(), "")
RESULT

CHECK(~ProteinHit())
	ptr1 = new ProteinHit();	
  delete ptr1;
RESULT

CHECK(ProteinHit(const ProteinHit& source))
	ProteinHit source;
	source.setScore(score);
	source.setScoreType(score_type);
	source.setRank(rank);
	source.setAccession(accession);
	source.setAccessionType(accession_type);
	source.setSequence(sequence);

  ptr1 = new ProteinHit(source);
	TEST_EQUAL(ptr1->getScore(), source.getScore())
	TEST_EQUAL(ptr1->getScoreType(), source.getScoreType())
	TEST_EQUAL(ptr1->getRank(), source.getRank())
	TEST_EQUAL(ptr1->getAccession(), source.getAccession())
	TEST_EQUAL(ptr1->getAccessionType(), source.getAccessionType())
	TEST_EQUAL(ptr1->getSequence(), source.getSequence())		  
RESULT

CHECK((ProteinHit(DoubleReal score, std::string score_type, UnsignedInt rank, String accession, std::string accession_type, String sequence)))
	ptr1 = new ProteinHit(score, score_type, rank, accession, accession_type, sequence);
	TEST_EQUAL(ptr1->getScore(), score)
	TEST_EQUAL(ptr1->getScoreType(), score_type)
	TEST_EQUAL(ptr1->getRank(), rank)
	TEST_EQUAL(ptr1->getAccession(), accession)
	TEST_EQUAL(ptr1->getAccessionType(), accession_type)	
	TEST_EQUAL(ptr1->getSequence(), sequence)
RESULT

CHECK(ProteinHit& operator=(const ProteinHit& source))
	ProteinHit hit;
	ptr1 = new ProteinHit(score, score_type, rank, accession, accession_type, sequence);
	hit = *ptr1;
	TEST_EQUAL(hit.getScore(), ptr1->getScore())
	TEST_EQUAL(hit.getScoreType(), ptr1->getScoreType())
	TEST_EQUAL(hit.getRank(), ptr1->getRank())
	TEST_EQUAL(hit.getAccession(), ptr1->getAccession())
	TEST_EQUAL(hit.getAccessionType(), ptr1->getAccessionType())	
	TEST_EQUAL(hit.getSequence(), ptr1->getSequence())
RESULT

CHECK(bool operator != (const ProteinHit& rhs) const)
  ProteinHit hit, hit2;
  TEST_EQUAL(hit!=hit2,false);
  hit.setScore(score);
  TEST_EQUAL(hit!=hit2,true);
  hit2.setScore(score);
  TEST_EQUAL(hit!=hit2,false);
  hit.setScoreType(score_type);
  TEST_EQUAL(hit!=hit2,true);
  hit2.setScoreType(score_type);
  TEST_EQUAL(hit!=hit2,false);
	hit.setRank(rank);
  TEST_EQUAL(hit!=hit2,true);
	hit2.setRank(rank);
  TEST_EQUAL(hit!=hit2,false);
	hit.setAccession(accession);
  TEST_EQUAL(hit!=hit2,true);
	hit2.setAccession(accession);
	TEST_EQUAL(hit!=hit2,false);
	hit.setAccessionType(accession_type);
  TEST_EQUAL(hit!=hit2,true);
	hit2.setAccessionType(accession_type);
	TEST_EQUAL(hit!=hit2,false);
	hit.setSequence(sequence);
  TEST_EQUAL(hit!=hit2,true);
	hit2.setSequence(sequence);
	TEST_EQUAL(hit!=hit2,false);
RESULT

CHECK(bool operator == (const ProteinHit& rhs) const)
  ProteinHit hit, hit2;
  TEST_EQUAL(hit==hit2,true);
  hit.setScore(score);
  TEST_EQUAL(hit==hit2,false);
  hit2.setScore(score);
  TEST_EQUAL(hit==hit2,true);
  hit.setScoreType(score_type);
  TEST_EQUAL(hit==hit2,false);
  hit2.setScoreType(score_type);
  TEST_EQUAL(hit==hit2,true);
	hit.setRank(rank);
  TEST_EQUAL(hit==hit2,false);
	hit2.setRank(rank);
  TEST_EQUAL(hit==hit2,true);
	hit.setAccession(accession);
  TEST_EQUAL(hit==hit2,false);
	hit2.setAccession(accession);
	TEST_EQUAL(hit==hit2,true);
	hit.setAccessionType(accession_type);
  TEST_EQUAL(hit==hit2,false);
	hit2.setAccessionType(accession_type);
	TEST_EQUAL(hit==hit2,true);
	hit.setSequence(sequence);
  TEST_EQUAL(hit==hit2,false);
	hit2.setSequence(sequence);
	TEST_EQUAL(hit==hit2,true);
RESULT

CHECK(const String& getAccession() const)
	ptr1 = new ProteinHit(score, score_type, rank, accession, accession_type, sequence);
	TEST_EQUAL(ptr1->getAccession(), accession)
RESULT

CHECK(const String& getSequence() const)
	ptr1 = new ProteinHit(score, score_type, rank, accession, accession_type, sequence);
	TEST_EQUAL(ptr1->getSequence(), sequence)
RESULT

CHECK(Real getScore() const)
	ptr1 = new ProteinHit(score, score_type, rank, accession, accession_type, sequence);
	TEST_EQUAL(ptr1->getScore(), score)
RESULT

CHECK(const std::string& getScoreType() const)
	ptr1 = new ProteinHit(score, score_type, rank, accession, accession_type, sequence);
	TEST_EQUAL(ptr1->getScoreType(), score_type)
RESULT

CHECK(UnsignedInt getRank() const)
	ptr1 = new ProteinHit(score, score_type, rank, accession, accession_type, sequence);
	TEST_EQUAL(ptr1->getRank(), rank)
RESULT

CHECK(const std::string& getAccessionType() const)
	ptr1 = new ProteinHit(score, score_type, rank, accession, accession_type, sequence);
	TEST_EQUAL(ptr1->getAccessionType(), accession_type)
RESULT

CHECK(void clear())
	ProteinHit hit;	
	ptr1 = new ProteinHit(score, score_type, rank, accession, accession_type, sequence);
	ptr1->clear();
  TEST_EQUAL(hit==*ptr1,true);
RESULT

CHECK(void setRank(UnsignedInt newrank))
	ptr1 = new ProteinHit();
	ptr1->setRank(rank);
	TEST_EQUAL(ptr1->getRank(), rank)	
RESULT

CHECK(void setScore(const DoubleReal& score))
	ptr1->setScore(score);
	TEST_EQUAL(ptr1->getScore(), score);
RESULT

CHECK(void setScoreType(const std::string& score_type))
	ptr1->setScoreType(score_type);
	TEST_EQUAL(ptr1->getScoreType(), score_type)
RESULT

CHECK(void setSequence(const String& sequence))
	ptr1->setSequence(sequence);
	TEST_EQUAL(ptr1->getSequence(), sequence)
RESULT

CHECK(void setAccession(const String& accession))
	ptr1->setAccession(accession);
	TEST_EQUAL(ptr1->getAccession(), accession)
RESULT

CHECK(void setAccessionType(const std::string& accession_type))
	ptr1->setAccessionType(accession_type);
	TEST_EQUAL(ptr1->getAccessionType(), accession_type)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
