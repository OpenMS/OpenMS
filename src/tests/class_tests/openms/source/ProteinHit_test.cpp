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
#include <unordered_set>
#include <unordered_map>

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

START_SECTION((bool isDecoy() const))
{
  ProteinHit hit;
  
  // Test default behavior (no target_decoy meta value set)
  TEST_EQUAL(hit.isDecoy(), false);
  
  // Test with explicit "target" value
  hit.setMetaValue("target_decoy", "target");
  TEST_EQUAL(hit.isDecoy(), false);
  
  // Test with "decoy" value
  hit.setMetaValue("target_decoy", "decoy");
  TEST_EQUAL(hit.isDecoy(), true);
  
  // Test with "DECOY" (case insensitive)
  hit.setMetaValue("target_decoy", "DECOY");
  TEST_EQUAL(hit.isDecoy(), true);
  
  // Test after removing meta value
  hit.removeMetaValue("target_decoy");
  TEST_EQUAL(hit.isDecoy(), false);
}
END_SECTION

START_SECTION((void setTargetDecoyType(TargetDecoyType type)))
{
  ProteinHit hit;
  
  // Test setting TARGET
  hit.setTargetDecoyType(ProteinHit::TargetDecoyType::TARGET);
  TEST_EQUAL(hit.getMetaValue("target_decoy"), "target");
  
  // Test setting DECOY
  hit.setTargetDecoyType(ProteinHit::TargetDecoyType::DECOY);
  TEST_EQUAL(hit.getMetaValue("target_decoy"), "decoy");
  
  // Test setting UNKNOWN (should remove meta value)
  hit.setTargetDecoyType(ProteinHit::TargetDecoyType::UNKNOWN);
  TEST_EQUAL(hit.metaValueExists("target_decoy"), false);
  TEST_EQUAL(hit.getTargetDecoyType(), ProteinHit::TargetDecoyType::UNKNOWN);
}
END_SECTION

START_SECTION((TargetDecoyType getTargetDecoyType() const))
{
  ProteinHit hit;
  
  // Test default behavior (should return UNKNOWN when meta value doesn't exist)
  TEST_EQUAL(hit.getTargetDecoyType(), ProteinHit::TargetDecoyType::UNKNOWN);
  
  // Test with explicit "target" value
  hit.setMetaValue("target_decoy", "target");
  TEST_EQUAL(hit.getTargetDecoyType(), ProteinHit::TargetDecoyType::TARGET);
  
  // Test with "decoy" value
  hit.setMetaValue("target_decoy", "decoy");
  TEST_EQUAL(hit.getTargetDecoyType(), ProteinHit::TargetDecoyType::DECOY);
  
  // Test with "DECOY" (case insensitive)
  hit.setMetaValue("target_decoy", "DECOY");
  TEST_EQUAL(hit.getTargetDecoyType(), ProteinHit::TargetDecoyType::DECOY);
  
  // Test after removing meta value (should return UNKNOWN)
  hit.removeMetaValue("target_decoy");
  TEST_EQUAL(hit.getTargetDecoyType(), ProteinHit::TargetDecoyType::UNKNOWN);
}
END_SECTION

START_SECTION(([EXTRA] std::hash<ProteinHit>))
{
  // Test that equal ProteinHits have equal hashes
  ProteinHit hit1(4.5, 1, "ACC123", "PEPTIDE");
  hit1.setCoverage(45.5);

  ProteinHit hit2(4.5, 1, "ACC123", "PEPTIDE");
  hit2.setCoverage(45.5);

  TEST_EQUAL(hit1 == hit2, true)
  TEST_EQUAL(std::hash<ProteinHit>{}(hit1), std::hash<ProteinHit>{}(hit2))

  // Test that different ProteinHits have different hashes (not guaranteed but expected)
  ProteinHit hit3(4.5, 1, "ACC456", "PEPTIDE");
  hit3.setCoverage(45.5);
  TEST_NOT_EQUAL(std::hash<ProteinHit>{}(hit1), std::hash<ProteinHit>{}(hit3))

  // Test different scores
  ProteinHit hit4(5.5, 1, "ACC123", "PEPTIDE");
  hit4.setCoverage(45.5);
  TEST_NOT_EQUAL(std::hash<ProteinHit>{}(hit1), std::hash<ProteinHit>{}(hit4))

  // Test different rank
  ProteinHit hit5(4.5, 2, "ACC123", "PEPTIDE");
  hit5.setCoverage(45.5);
  TEST_NOT_EQUAL(std::hash<ProteinHit>{}(hit1), std::hash<ProteinHit>{}(hit5))

  // Test different sequence
  ProteinHit hit6(4.5, 1, "ACC123", "PROTEIN");
  hit6.setCoverage(45.5);
  TEST_NOT_EQUAL(std::hash<ProteinHit>{}(hit1), std::hash<ProteinHit>{}(hit6))

  // Test different coverage
  ProteinHit hit7(4.5, 1, "ACC123", "PEPTIDE");
  hit7.setCoverage(50.0);
  TEST_NOT_EQUAL(std::hash<ProteinHit>{}(hit1), std::hash<ProteinHit>{}(hit7))

  // Test use in unordered_set
  std::unordered_set<ProteinHit> hit_set;
  hit_set.insert(hit1);
  hit_set.insert(hit2); // Should not add duplicate
  hit_set.insert(hit3); // Should add (different accession)
  TEST_EQUAL(hit_set.size(), 2)
  TEST_EQUAL(hit_set.count(hit1), 1)
  TEST_EQUAL(hit_set.count(hit3), 1)

  // Test use in unordered_map
  std::unordered_map<ProteinHit, std::string> hit_map;
  hit_map[hit1] = "first";
  hit_map[hit3] = "second";
  TEST_EQUAL(hit_map.size(), 2)
  TEST_EQUAL(hit_map[hit1], "first")
  TEST_EQUAL(hit_map[hit2], "first") // hit2 == hit1, so should get same value
  TEST_EQUAL(hit_map[hit3], "second")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
