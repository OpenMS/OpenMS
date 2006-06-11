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
// $Id: PeptideHit_test.C,v 1.5 2006/06/09 23:47:35 nicopfeifer Exp $
// $Author: nicopfeifer $
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <string>

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/DATASTRUCTURES/String.h>

///////////////////////////

START_TEST(PeptideHit, "$Id: PeptideHit_test.C,v 1.5 2006/06/09 23:47:35 nicopfeifer Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

float score = 4.4;
std::string score_type = "XCorr";
uint rank = 3;
String sequence = "ARRAY";
std::string sequence2 = "  ARRAY  ";

PeptideHit* ptr1 = 0;

CHECK(PeptideHit())
	ptr1 = new PeptideHit();
	TEST_NOT_EQUAL(ptr1, 0)
	TEST_EQUAL(ptr1->getScore(), 0)
	TEST_EQUAL(ptr1->getScoreType(), "")
	TEST_EQUAL(ptr1->getRank(), 0)
	TEST_EQUAL(ptr1->getSequence(), "")
RESULT

CHECK(~PeptideHit())
	delete ptr1;
RESULT

CHECK(PeptideHit(score, score_type, rank, sequence))
	ptr1 = new PeptideHit(score, score_type, rank, sequence);
	TEST_EQUAL(ptr1->getScore(), score)
	TEST_EQUAL(ptr1->getScoreType(), score_type)
	TEST_EQUAL(ptr1->getRank(), rank)
	TEST_EQUAL(ptr1->getSequence(), sequence)
RESULT

CHECK(PeptideHit& operator=(const PeptideHit& source))
	PeptideHit hit;
	ptr1 = new PeptideHit(score, score_type, rank, sequence);
	hit = *ptr1;
	TEST_EQUAL(hit.getScore(), ptr1->getScore())
	TEST_EQUAL(hit.getScoreType(), ptr1->getScoreType())
	TEST_EQUAL(hit.getRank(), ptr1->getRank())
	TEST_EQUAL(hit.getSequence(), ptr1->getSequence())		
RESULT

CHECK(PeptideHit(const PeptideHit& source))
	PeptideHit source;
	source.setScore(score);
	source.setScoreType(score_type);
	source.setRank(rank);
	source.setSequence(sequence);

  ptr1 = new PeptideHit(source);
	TEST_EQUAL(ptr1->getScore(), source.getScore())
	TEST_EQUAL(ptr1->getScoreType(), source.getScoreType())
	TEST_EQUAL(ptr1->getRank(), source.getRank())
	TEST_EQUAL(ptr1->getSequence(), source.getSequence())		  
RESULT

CHECK(bool operator== (const PeptideHit& rhs))
  PeptideHit hit, hit2;
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
	hit.setSequence(sequence);
  TEST_EQUAL(hit==hit2,false);
	hit2.setSequence(sequence);
	TEST_EQUAL(hit==hit2,true);
RESULT

CHECK(bool operator!= (const Date& rhs))
  PeptideHit hit, hit2;
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
	hit.setSequence(sequence);
  TEST_EQUAL(hit!=hit2,true);
	hit2.setSequence(sequence);
	TEST_EQUAL(hit!=hit2,false);
RESULT

CHECK(const double& getScore() const)
	ptr1 = new PeptideHit(score, score_type, rank, sequence);
	TEST_EQUAL(ptr1->getScore(), score)
RESULT

CHECK(const std::string& getScoreType() const)
	ptr1 = new PeptideHit(score, score_type, rank, sequence);
	TEST_EQUAL(ptr1->getScoreType(), score_type)
RESULT

CHECK(const uint& getRank() const)
	ptr1 = new PeptideHit(score, score_type, rank, sequence);
	TEST_EQUAL(ptr1->getRank(), rank)
RESULT

CHECK(std::string getSequence() const)
	ptr1 = new PeptideHit(score, score_type, rank, sequence);
	TEST_EQUAL(ptr1->getSequence(), sequence)
RESULT

CHECK(void clear())
	ptr1 = new PeptideHit(score, score_type, rank, sequence);
	ptr1->clear();
	TEST_EQUAL(ptr1->getScore(), 0)
	TEST_EQUAL(ptr1->getScoreType(), "")
	TEST_EQUAL(ptr1->getRank(), 0)
	TEST_EQUAL(ptr1->getSequence(), "")
RESULT

CHECK(void serialize(PersistenceManager& pm))
  // ???
RESULT

CHECK(void setRank(uint newrank))
	ptr1 = new PeptideHit();
	ptr1->setRank(rank);
	TEST_EQUAL(ptr1->getRank(), rank)
RESULT

CHECK(void setScore(const double& score))
	ptr1 = new PeptideHit();
	ptr1->setScore(score);
	TEST_EQUAL(ptr1->getScore(), score)
RESULT

CHECK(void setScoreType(const std::string& score_type))
	ptr1 = new PeptideHit();
	ptr1->setScoreType(score_type);
	TEST_EQUAL(ptr1->getScoreType(), score_type)
RESULT

CHECK(void setSequence(const std::string& sequence))
	ptr1 = new PeptideHit();
	ptr1->setSequence(sequence);
	TEST_EQUAL(ptr1->getSequence(), sequence)
	ptr1->setSequence(sequence2);
	TEST_EQUAL(ptr1->getSequence(), sequence)	
RESULT

CHECK(~PeptideHit())
	ptr1 = new PeptideHit();
  delete ptr1;
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
