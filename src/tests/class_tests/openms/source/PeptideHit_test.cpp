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

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/DATASTRUCTURES/String.h>

///////////////////////////

START_TEST(PeptideHit, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

double score = 4.4;
UInt rank = 3;
AASequence sequence = AASequence::fromString("ARRAY");
std::string sequence2 = "  ARRAY  ";
Int charge = 2;

PeptideHit* ptr = nullptr;
PeptideHit* nullPointer = nullptr;
START_SECTION((PeptideHit()))
	ptr = new PeptideHit();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeptideHit()))
	delete ptr;
END_SECTION

START_SECTION((PeptideHit(double score, UInt rank, Int charge, const AASequence &sequence)))
	PeptideHit hit(score, rank, charge, sequence);
	TEST_EQUAL(hit.getScore(), score)
	TEST_EQUAL(hit.getRank(), rank)
	TEST_EQUAL(hit.getCharge(), charge)
	TEST_EQUAL(hit.getSequence(), sequence)
END_SECTION

START_SECTION((PeptideHit& operator=(const PeptideHit& source)))
	PeptideHit hit;
	PeptideHit hit2(score, rank, charge, sequence);
	hit2.setMetaValue("label",17);
	
	hit = hit2;
	
	TEST_EQUAL(hit.getScore(), score)
	TEST_EQUAL(hit.getRank(), rank)
	TEST_EQUAL(hit.getCharge(), charge)
	TEST_EQUAL(hit.getSequence(), sequence)
	TEST_EQUAL((UInt)hit.getMetaValue("label"),17)
END_SECTION

START_SECTION((PeptideHit(const PeptideHit& source)))
	PeptideHit source;
	source.setScore(score);
	source.setRank(rank);
	source.setSequence(sequence);
	source.setMetaValue("label",17);
	
  PeptideHit hit(source);
	
	TEST_EQUAL(hit.getScore(), source.getScore())
	TEST_EQUAL(hit.getRank(), source.getRank())
	TEST_EQUAL(hit.getSequence(), source.getSequence())
	TEST_EQUAL((UInt)hit.getMetaValue("label"),17) 
END_SECTION

START_SECTION((bool operator == (const PeptideHit& rhs) const))
  PeptideHit hit, hit2;
  TEST_EQUAL(hit==hit2,true);

  hit.setScore(score);
  TEST_EQUAL(hit==hit2,false);
	hit=hit2;
	
  hit.setRank(rank);
  TEST_EQUAL(hit==hit2,false);
	hit=hit2;
	
	hit.setSequence(sequence);
  TEST_EQUAL(hit==hit2,false);
	hit=hit2;
	
	hit.setMetaValue("label",17);
  TEST_EQUAL(hit==hit2,false);
	hit=hit2;
END_SECTION

START_SECTION((bool operator != (const PeptideHit& rhs) const))
  PeptideHit hit, hit2;
  TEST_EQUAL(hit!=hit2,false);

  hit.setScore(score);
  TEST_EQUAL(hit!=hit2,true);
	hit=hit2;
	
  hit.setRank(rank);
  TEST_EQUAL(hit!=hit2,true);
	hit=hit2;
	
	hit.setSequence(sequence);
  TEST_EQUAL(hit!=hit2,true);
	hit=hit2;
	
	hit.setMetaValue("label",17);
  TEST_EQUAL(hit!=hit2,true);
	hit=hit2;
END_SECTION

START_SECTION((double getScore() const ))
	PeptideHit hit(score, rank, charge, sequence);
	TEST_EQUAL(hit.getScore(), score)
END_SECTION

START_SECTION((UInt getRank() const))
	PeptideHit hit(score, rank, charge, sequence);
	TEST_EQUAL(hit.getRank(), rank)
END_SECTION

START_SECTION((const AASequence& getSequence() const))
	PeptideHit hit(score, rank, charge, sequence);
	TEST_EQUAL(hit.getSequence(), sequence)
END_SECTION

START_SECTION((void setRank(UInt newrank)))
	PeptideHit hit;
	hit.setRank(rank);
	TEST_EQUAL(hit.getRank(), rank)
END_SECTION

START_SECTION((void setScore(double score)))
	PeptideHit hit;
	hit.setScore(score);
	TEST_EQUAL(hit.getScore(), score)
END_SECTION

START_SECTION((void setSequence(const AASequence& sequence)))
	PeptideHit hit;
	hit.setSequence(sequence);
	TEST_EQUAL(hit.getSequence(), sequence)
	//hit.setSequence(sequence2);
	// @todo std::string interface?
	TEST_EQUAL(hit.getSequence(), sequence)	
END_SECTION

        ;
START_SECTION((void setPeptideEvidences(const vector<PeptideEvidence> & peptide_evidences)))
     PeptideHit hit;
     vector<PeptideEvidence> pes(2, PeptideEvidence());
     pes[0].setProteinAccession("ACC392");
     pes[1].setProteinAccession("ACD392");
     hit.setPeptideEvidences(pes);
    TEST_EQUAL(hit.getPeptideEvidences().size(), 2)
    TEST_EQUAL(hit.getPeptideEvidences()[0].getProteinAccession() == String("ACC392"), true)
    TEST_EQUAL(hit.getPeptideEvidences()[1].getProteinAccession() == String("ACD392"), true)
END_SECTION


START_SECTION((const std::set<String>& extractProteinAccessionsSet() const))
     PeptideHit hit;
     vector<PeptideEvidence> pes(2, PeptideEvidence());
     pes[0].setProteinAccession("ACC392");
     pes[1].setProteinAccession("ACD392");
     hit.setPeptideEvidences(pes);
     TEST_EQUAL(hit.extractProteinAccessionsSet().size(), 2)
     TEST_EQUAL(*hit.extractProteinAccessionsSet().begin(), "ACC392")
     TEST_EQUAL(*hit.extractProteinAccessionsSet().rbegin(), "ACD392")
END_SECTION

START_SECTION((Int getCharge() const))
	PeptideHit hit;
	
	hit.setCharge(-43);
	TEST_EQUAL(-43, hit.getCharge())
END_SECTION

START_SECTION((void setCharge(Int charge)))
	PeptideHit hit;
	
	hit.setCharge(-43);
	TEST_EQUAL(-43, hit.getCharge())
END_SECTION
/*
START_SECTION((void setAABefore(char acid)))
	PeptideHit hit;
	
	hit.setAABefore('R');
	TEST_EQUAL(hit.getAABefore(), 'R')
END_SECTION
START_SECTION((char getAABefore() const))
	PeptideHit hit;
	
	hit.setAABefore('R');
	TEST_EQUAL(hit.getAABefore(), 'R')
END_SECTION
START_SECTION((void setAAAfter(char acid)))
	PeptideHit hit;
	
	hit.setAAAfter('R');
	TEST_EQUAL(hit.getAAAfter(), 'R')
END_SECTION
START_SECTION((char getAAAfter() const))
	PeptideHit hit;
	
	hit.setAAAfter('R');
	TEST_EQUAL(hit.getAAAfter(), 'R')
END_SECTION
*/
START_SECTION(([PeptideHit::ScoreLess] template < typename Arg > bool operator()(const Arg &a, const Arg &b)))
{
  PeptideHit a,b;
  a.setScore(10);
  b.setScore(20);

  TEST_EQUAL(PeptideHit::ScoreLess().operator()(a,b), true)
  TEST_EQUAL(PeptideHit::ScoreLess().operator()(b,a), false)
  TEST_EQUAL(PeptideHit::ScoreLess().operator()(a,a), false)
}
END_SECTION

START_SECTION(([PeptideHit::ScoreMore] template < typename Arg > bool operator()(const Arg &a, const Arg &b)))
{
  PeptideHit a,b;
  a.setScore(20);
  b.setScore(10);

  TEST_EQUAL(PeptideHit::ScoreMore().operator()(a,b), true)
  TEST_EQUAL(PeptideHit::ScoreMore().operator()(b,a), false)
  TEST_EQUAL(PeptideHit::ScoreMore().operator()(a,a), false)
}
END_SECTION

START_SECTION((void setPeakAnnotations(const vector<PeptideHit::PeakAnnotation> & fragment_annotations)))
  PeptideHit hit;
  vector<PeptideHit::PeakAnnotation> frag_annos(2, PeptideHit::PeakAnnotation());
  frag_annos[0].annotation = "test string";
  frag_annos[0].charge = 2;
  frag_annos[0].mz = 1234.567;
  frag_annos[0].intensity = 1.0;
  frag_annos[1].annotation = "second test string";
  frag_annos[1].charge = 1;
  frag_annos[1].mz = 89.10;
  frag_annos[1].intensity = 0.5;
  hit.setPeakAnnotations(frag_annos);
  TEST_EQUAL(hit.getPeakAnnotations().size(), 2)
  TEST_EQUAL(hit.getPeakAnnotations()[0].annotation == "test string", true)
  TEST_EQUAL(hit.getPeakAnnotations()[0].charge == 2, true)
  TEST_EQUAL(hit.getPeakAnnotations()[0].mz == 1234.567, true)
  TEST_EQUAL(hit.getPeakAnnotations()[0].intensity == 1.0, true)
  TEST_EQUAL(hit.getPeakAnnotations()[1].annotation == "second test string", true)
  TEST_EQUAL(hit.getPeakAnnotations()[1].mz == 89.1, true)
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
