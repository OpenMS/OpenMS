// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

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


float score = 4.4f;
UInt rank = 3;
String sequence = "ARRAY";
String accession = "PROOE34";
String description = "class II antigen";

ProteinHit* ptr = nullptr;	
ProteinHit* nullPointer = nullptr;
START_SECTION(ProteinHit())
	ptr = new ProteinHit();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~ProteinHit())
  delete ptr;
END_SECTION

START_SECTION((ProteinHit(double score, UInt rank, String accession, String sequence)))
	ProteinHit hit(score, rank, accession, sequence);
	TEST_EQUAL(hit.getScore(), score)
	TEST_EQUAL(hit.getRank(), rank)
	TEST_EQUAL(hit.getAccession(), accession)
	TEST_EQUAL(hit.getSequence(), sequence)
	TEST_EQUAL(hit.getCoverage(), -1)
END_SECTION

START_SECTION(ProteinHit(const ProteinHit& source))
	ProteinHit source;
	source.setScore(score);
	source.setRank(rank);
	source.setAccession(accession);
	source.setDescription(description);
	source.setSequence(sequence);
	source.setMetaValue("label",17);
  source.setCoverage(123.123);
  
  ProteinHit hit(source);
	
	TEST_EQUAL(hit.getScore(), score)
	TEST_EQUAL(hit.getRank(), rank)
	TEST_EQUAL(hit.getAccession(), accession)
	TEST_EQUAL(hit.getDescription(), description)
	TEST_EQUAL(hit.getSequence(), sequence)
	TEST_EQUAL((UInt)hit.getMetaValue("label"),17)
	TEST_EQUAL(hit.getCoverage(), 123.123)
END_SECTION

START_SECTION(ProteinHit& operator=(const ProteinHit& source))
	ProteinHit hit;
	ProteinHit hit2(score, rank, accession, sequence);
	hit2.setMetaValue("label",17);
	hit2.setCoverage(123.123);
	hit2.setDescription(description);
	
	hit = hit2;
	
	TEST_EQUAL(hit.getScore(), score)
	TEST_EQUAL(hit.getRank(), rank)
	TEST_EQUAL(hit.getAccession(), accession)
	TEST_EQUAL(hit.getDescription(), description)
	TEST_EQUAL(hit.getSequence(), sequence)
	TEST_EQUAL((UInt)hit.getMetaValue("label"),17)
	TEST_EQUAL(hit.getCoverage(), 123.123)
END_SECTION

START_SECTION(ProteinHit& operator= (const MetaInfoInterface& source))
	ProteinHit hit(score, rank, accession, sequence);
	hit.setCoverage(123.123);
  MetaInfoInterface meta;
	meta.setMetaValue("label",17);
	
	hit = meta;
	
	TEST_EQUAL(hit.getScore(), score)
	TEST_EQUAL(hit.getRank(), rank)
	TEST_EQUAL(hit.getAccession(), accession)
	TEST_EQUAL(hit.getSequence(), sequence)
	TEST_EQUAL(hit.getCoverage(), 123.123)
	TEST_EQUAL((UInt)hit.getMetaValue("label"),17)
END_SECTION


START_SECTION(bool operator == (const ProteinHit& rhs) const)
  ProteinHit hit, hit2;
  TEST_EQUAL(hit==hit2,true);
  
  hit.setScore(score);
  TEST_EQUAL(hit==hit2,false);
	hit = hit2;
  
	hit.setRank(rank);
  TEST_EQUAL(hit==hit2,false);
	hit = hit2;
  
	hit.setAccession(accession);
  TEST_EQUAL(hit==hit2,false);
	hit = hit2;
  
	hit.setSequence(sequence);
  TEST_EQUAL(hit==hit2,false);
	hit = hit2;
	
	hit.setMetaValue("label",17);
  TEST_EQUAL(hit==hit2,false);
	hit = hit2;
	
	hit.setCoverage(123.123);
  TEST_EQUAL(hit==hit2,false);
	hit = hit2;	
END_SECTION

START_SECTION(bool operator != (const ProteinHit& rhs) const)
  ProteinHit hit, hit2;
  TEST_EQUAL(hit!=hit2,false);
  
  hit.setScore(score);
  TEST_EQUAL(hit!=hit2,true);
  hit = hit2;
  
  hit.setRank(rank);
  TEST_EQUAL(hit!=hit2,true);
	hit = hit2;
  
  hit.setAccession(accession);
  TEST_EQUAL(hit!=hit2,true);
	hit = hit2;
  
	hit.setSequence(sequence);
  TEST_EQUAL(hit!=hit2,true);
	hit = hit2;

	hit.setMetaValue("label",17);
  TEST_EQUAL(hit!=hit2,true);
	hit = hit2;

	hit.setCoverage(123.123);
  TEST_EQUAL(hit!=hit2,true);
	hit = hit2;		
END_SECTION

START_SECTION(const String& getAccession() const)
	ProteinHit hit(score, rank, accession, sequence);
	TEST_EQUAL(hit.getAccession(), accession)
END_SECTION

START_SECTION(const String& getDescription() const)
	ProteinHit hit(score, rank, accession, sequence);
  hit.setDescription(description);
	TEST_EQUAL(hit.getDescription(), description)
END_SECTION

START_SECTION(const String& getSequence() const)
	ProteinHit hit(score, rank, accession, sequence);
	TEST_EQUAL(hit.getSequence(), sequence)
END_SECTION

START_SECTION(float getScore() const)
	ProteinHit hit(score, rank, accession, sequence);
	TEST_EQUAL(hit.getScore(), score)
END_SECTION

START_SECTION(UInt getRank() const)
	ProteinHit hit(score, rank, accession, sequence);
	TEST_EQUAL(hit.getRank(), rank)
END_SECTION


START_SECTION(double getCoverage() const)
	ProteinHit hit(score, rank, accession, sequence);
	TEST_EQUAL(hit.getCoverage(), -1)
	hit.setCoverage(123.123);
	TEST_EQUAL(hit.getCoverage(), 123.123)
END_SECTION


START_SECTION(void setRank(UInt newrank))
	ProteinHit hit;
	hit.setRank(rank);
	TEST_EQUAL(hit.getRank(), rank)	
END_SECTION

START_SECTION(void setScore(const double score))
	ProteinHit hit;
	hit.setScore(score);
	TEST_EQUAL(hit.getScore(), score);
END_SECTION

START_SECTION(void setSequence(const String& sequence))
	ProteinHit hit;
	hit.setSequence(sequence);
	TEST_EQUAL(hit.getSequence(), sequence)
END_SECTION

START_SECTION(void setAccession(const String& accession))
	ProteinHit hit;
	hit.setAccession(accession);
	TEST_EQUAL(hit.getAccession(), accession)
END_SECTION

START_SECTION(void setDescription(const String& description))
	ProteinHit hit;
	hit.setDescription(description);
	TEST_EQUAL(hit.getDescription(), description)
END_SECTION

START_SECTION(void setCoverage(const double coverage))
	ProteinHit hit;
	hit.setCoverage(123.123);
	TEST_EQUAL(hit.getCoverage(), 123.123)
END_SECTION

START_SECTION(([ProteinHit::ScoreLess] template < typename Arg > bool operator()(const Arg &a, const Arg &b)))
{
  ProteinHit a,b;
  a.setScore(10);
  b.setScore(20);

  TEST_EQUAL(ProteinHit::ScoreLess().operator()(a,b), true)
  TEST_EQUAL(ProteinHit::ScoreLess().operator()(b,a), false)
  TEST_EQUAL(ProteinHit::ScoreLess().operator()(a,a), false)
}
END_SECTION

START_SECTION(([ProteinHit::ScoreMore] template < typename Arg > bool operator()(const Arg &a, const Arg &b)))
{
  ProteinHit a,b;
  a.setScore(20);
  b.setScore(10);

  TEST_EQUAL(ProteinHit::ScoreMore().operator()(a,b), true)
  TEST_EQUAL(ProteinHit::ScoreMore().operator()(b,a), false)
  TEST_EQUAL(ProteinHit::ScoreMore().operator()(a,a), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
