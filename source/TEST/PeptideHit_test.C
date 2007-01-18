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

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/DATASTRUCTURES/String.h>

///////////////////////////

START_TEST(PeptideHit, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

float score = 4.4;
std::string score_type = "XCorr";
uint rank = 3;
String sequence = "ARRAY";
std::string sequence2 = "  ARRAY  ";
SignedInt charge;

PeptideHit* ptr1 = 0;

CHECK((PeptideHit()))
	ptr1 = new PeptideHit();
	TEST_NOT_EQUAL(ptr1, 0)
	TEST_EQUAL(ptr1->getScore(), 0)
	TEST_EQUAL(ptr1->getScoreType(), "")
	TEST_EQUAL(ptr1->getRank(), 0)
	TEST_EQUAL(ptr1->getCharge(), 0)
	TEST_EQUAL(ptr1->getSequence(), "")
RESULT

CHECK((~PeptideHit()))
	delete ptr1;
RESULT

CHECK((PeptideHit(double score, std::string score_type, uint rank, SignedInt charge, String sequence)))
	ptr1 = new PeptideHit(score, score_type, rank, charge, sequence);
	TEST_EQUAL(ptr1->getScore(), score)
	TEST_EQUAL(ptr1->getScoreType(), score_type)
	TEST_EQUAL(ptr1->getRank(), rank)
	TEST_EQUAL(ptr1->getCharge(), charge)
	TEST_EQUAL(ptr1->getSequence(), sequence)
RESULT

CHECK((PeptideHit& operator=(const PeptideHit& source)))
	PeptideHit hit;
	ptr1 = new PeptideHit(score, score_type, rank, charge, sequence);
	hit = *ptr1;
	TEST_EQUAL(hit.getScore(), ptr1->getScore())
	TEST_EQUAL(hit.getScoreType(), ptr1->getScoreType())
	TEST_EQUAL(hit.getRank(), ptr1->getRank())
	TEST_EQUAL(hit.getSequence(), ptr1->getSequence())		
RESULT

CHECK((PeptideHit(const PeptideHit& source)))
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

CHECK((bool operator == (const PeptideHit& rhs) const))
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

CHECK((bool operator != (const PeptideHit& rhs) const))
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

CHECK((float getScore() const))
	ptr1 = new PeptideHit(score, score_type, rank, charge, sequence);
	TEST_EQUAL(ptr1->getScore(), score)
RESULT

CHECK((const std::string& getScoreType() const))
	ptr1 = new PeptideHit(score, score_type, rank, charge, sequence);
	TEST_EQUAL(ptr1->getScoreType(), score_type)
RESULT

CHECK((UnsignedInt getRank() const))
	ptr1 = new PeptideHit(score, score_type, rank, charge, sequence);
	TEST_EQUAL(ptr1->getRank(), rank)
RESULT

CHECK((String getSequence() const))
	ptr1 = new PeptideHit(score, score_type, rank, charge, sequence);
	TEST_EQUAL(ptr1->getSequence(), sequence)
RESULT

CHECK((void clear()))
	ptr1 = new PeptideHit(score, score_type, rank, charge, sequence);
	ptr1->clear();
	TEST_EQUAL(ptr1->getScore(), 0)
	TEST_EQUAL(ptr1->getScoreType(), "")
	TEST_EQUAL(ptr1->getRank(), 0)
	TEST_EQUAL(ptr1->getSequence(), "")
RESULT

CHECK((void setRank(UnsignedInt newrank)))
	ptr1 = new PeptideHit();
	ptr1->setRank(rank);
	TEST_EQUAL(ptr1->getRank(), rank)
RESULT

CHECK((void setScore(const double& score)))
	ptr1 = new PeptideHit();
	ptr1->setScore(score);
	TEST_EQUAL(ptr1->getScore(), score)
RESULT

CHECK((void setScoreType(const std::string& score_type)))
	ptr1 = new PeptideHit();
	ptr1->setScoreType(score_type);
	TEST_EQUAL(ptr1->getScoreType(), score_type)
RESULT

CHECK((void setSequence(const String& sequence)))
	ptr1 = new PeptideHit();
	ptr1->setSequence(sequence);
	TEST_EQUAL(ptr1->getSequence(), sequence)
	ptr1->setSequence(sequence2);
	TEST_EQUAL(ptr1->getSequence(), sequence)	
RESULT

CHECK((void addProteinIndex(const std::pair<String, String>& index)))
	String date;
	vector< pair<String, String> > indices;

	date = "2006-12-12 11:59:59";
	ptr1 = new PeptideHit();

	ptr1->addProteinIndex(make_pair(date, String("ACC392")));
	ptr1->addProteinIndex(make_pair(date, String("ACC392")));
	ptr1->addProteinIndex(make_pair(date, String("ACD392")));
	indices = ptr1->getProteinIndices();
	TEST_EQUAL(indices.size(), 2)
	TEST_EQUAL(indices[0].first == String("2006-12-12 11:59:59"), true)
	TEST_EQUAL(indices[0].second == String("ACC392"), true)
	TEST_EQUAL(indices[1].first == String("2006-12-12 11:59:59"), true)
	TEST_EQUAL(indices[1].second == String("ACD392"), true)

RESULT

CHECK((void addProteinIndex(const DateTime& date, const String& accession)))
	DateTime date;
	vector< pair<String, String> > indices;

	date.set("2006-12-12 11:59:59");
	ptr1 = new PeptideHit();

	ptr1->addProteinIndex(date, "ACC392");
	ptr1->addProteinIndex(date, "ACC392");
	ptr1->addProteinIndex(date, "ACD392");
	indices = ptr1->getProteinIndices();
	TEST_EQUAL(indices.size(), 2)
	TEST_EQUAL(indices[0].first == String("2006-12-12 11:59:59"), true)
	TEST_EQUAL(indices[0].second == String("ACC392"), true)
	TEST_EQUAL(indices[1].first == String("2006-12-12 11:59:59"), true)
	TEST_EQUAL(indices[1].second == String("ACD392"), true)

RESULT

CHECK((void setProteinIndices(const std::vector< std::pair<String, String> >& indices)))
	vector< pair<String, String> > indices1;
	vector< pair<String, String> > indices2;

	ptr1 = new PeptideHit();

	indices1.push_back(make_pair("2006-12-12 11:59:59", "ACC392"));
	indices1.push_back(make_pair("2006-12-12 11:59:59", "ACC392"));
	indices1.push_back(make_pair("2006-12-12 11:59:59", "ACD392"));
	ptr1->setProteinIndices(indices1);
	indices2 = ptr1->getProteinIndices();
	TEST_EQUAL(indices1 == indices2, true)

RESULT

CHECK((std::vector< std::pair<String, String>& getProteinIndices()))
	DateTime date;
	vector< pair<String, String> > indices;

	date.set("2006-12-12 11:59:59");
	ptr1 = new PeptideHit();

	ptr1->addProteinIndex(date, "ACC392");
	ptr1->addProteinIndex(date, "ACC392");
	ptr1->addProteinIndex(date, "ACD392");
	indices = ptr1->getProteinIndices();
	TEST_EQUAL(indices.size(), 2)
	TEST_EQUAL(indices[0].first == String("2006-12-12 11:59:59"), true)
	TEST_EQUAL(indices[0].second == String("ACC392"), true)
	TEST_EQUAL(indices[1].first == String("2006-12-12 11:59:59"), true)
	TEST_EQUAL(indices[1].second == String("ACD392"), true)

RESULT

CHECK((const std::vector< std::pair<String, String>& getProteinIndices() const))
	DateTime date;

	date.set("2006-12-12 11:59:59");
	ptr1 = new PeptideHit();

	ptr1->addProteinIndex(date, "ACC392");
	ptr1->addProteinIndex(date, "ACC392");
	ptr1->addProteinIndex(date, "ACD392");
	const vector< pair<String, String> >& indices = ptr1->getProteinIndices();
	TEST_EQUAL(indices.size(), 2)
	TEST_EQUAL(indices[0].first == String("2006-12-12 11:59:59"), true)
	TEST_EQUAL(indices[0].second == String("ACC392"), true)
	TEST_EQUAL(indices[1].first == String("2006-12-12 11:59:59"), true)
	TEST_EQUAL(indices[1].second == String("ACD392"), true)

RESULT

CHECK(SignedInt getCharge() const)
	PeptideHit hit;
	
	hit.setCharge(-43);
	TEST_EQUAL(-43, hit.getCharge())
RESULT

CHECK(void setCharge(SignedInt charge))
	PeptideHit hit;
	
	hit.setCharge(-43);
	TEST_EQUAL(-43, hit.getCharge())
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
