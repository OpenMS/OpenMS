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

START_SECTION((void setAnalysisResults(std::vector<PepXMLAnalysisResult> aresult)))
{
  PeptideHit hit;
  
  // Create some analysis results
  std::vector<PeptideHit::PepXMLAnalysisResult> results;
  
  PeptideHit::PepXMLAnalysisResult ar1;
  ar1.score_type = "peptideprophet";
  ar1.higher_is_better = true;
  ar1.main_score = 0.95;
  ar1.sub_scores["fval"] = 0.7114;
  ar1.sub_scores["ntt"] = 2.0;
  
  PeptideHit::PepXMLAnalysisResult ar2;
  ar2.score_type = "interprophet";
  ar2.higher_is_better = true;
  ar2.main_score = 0.98;
  ar2.sub_scores["nss"] = 0.0;
  ar2.sub_scores["nrs"] = 10.2137;
  
  results.push_back(ar1);
  results.push_back(ar2);
  
  // Set the analysis results
  hit.setAnalysisResults(results);
  
  // Check that the analysis results were set correctly
  TEST_EQUAL(hit.getAnalysisResults().size(), 2);
  
  TEST_EQUAL(hit.getAnalysisResults()[0].score_type, "peptideprophet");
  TEST_EQUAL(hit.getAnalysisResults()[0].higher_is_better, true);
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[0].main_score, 0.95);
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[0].sub_scores.at("fval"), 0.7114);
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[0].sub_scores.at("ntt"), 2.0);
  
  TEST_EQUAL(hit.getAnalysisResults()[1].score_type, "interprophet");
  TEST_EQUAL(hit.getAnalysisResults()[1].higher_is_better, true);
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[1].main_score, 0.98);
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[1].sub_scores.at("nss"), 0.0);
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[1].sub_scores.at("nrs"), 10.2137);
  
  // Check that the analysis results are stored as meta values
  TEST_EQUAL(hit.metaValueExists("_ar_0_score_type"), true);
  TEST_EQUAL(hit.metaValueExists("_ar_0_score"), true);
  TEST_EQUAL(hit.metaValueExists("_ar_0_higher_is_better"), true);
  TEST_EQUAL(hit.metaValueExists("_ar_0_subscore_fval"), true);
  TEST_EQUAL(hit.metaValueExists("_ar_0_subscore_ntt"), true);
  
  TEST_EQUAL(hit.metaValueExists("_ar_1_score_type"), true);
  TEST_EQUAL(hit.metaValueExists("_ar_1_score"), true);
  TEST_EQUAL(hit.metaValueExists("_ar_1_higher_is_better"), true);
  TEST_EQUAL(hit.metaValueExists("_ar_1_subscore_nss"), true);
  TEST_EQUAL(hit.metaValueExists("_ar_1_subscore_nrs"), true);
  
  TEST_EQUAL(hit.getMetaValue("_ar_0_score_type").toString(), "peptideprophet");
  TEST_REAL_SIMILAR(hit.getMetaValue("_ar_0_score"), 0.95);
  TEST_EQUAL(hit.getMetaValue("_ar_0_higher_is_better").toBool(), true);
  TEST_REAL_SIMILAR(hit.getMetaValue("_ar_0_subscore_fval"), 0.7114);
  TEST_REAL_SIMILAR(hit.getMetaValue("_ar_0_subscore_ntt"), 2.0);
  
  TEST_EQUAL(hit.getMetaValue("_ar_1_score_type").toString(), "interprophet");
  TEST_REAL_SIMILAR(hit.getMetaValue("_ar_1_score"), 0.98);
  TEST_EQUAL(hit.getMetaValue("_ar_1_higher_is_better").toBool(), true);
  TEST_REAL_SIMILAR(hit.getMetaValue("_ar_1_subscore_nss"), 0.0);
  TEST_REAL_SIMILAR(hit.getMetaValue("_ar_1_subscore_nrs"), 10.2137);
  
  // Test overwriting existing analysis results
  PeptideHit::PepXMLAnalysisResult ar3;
  ar3.score_type = "mascot";
  ar3.higher_is_better = true;
  ar3.main_score = 100.0;
  
  std::vector<PeptideHit::PepXMLAnalysisResult> new_results;
  new_results.push_back(ar3);
  
  hit.setAnalysisResults(new_results);
  
  TEST_EQUAL(hit.getAnalysisResults().size(), 1);
  TEST_EQUAL(hit.getAnalysisResults()[0].score_type, "mascot");
  TEST_EQUAL(hit.getAnalysisResults()[0].higher_is_better, true);
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[0].main_score, 100.0);
  
  // Check that the old meta values are gone
  TEST_EQUAL(hit.metaValueExists("_ar_0_score_type"), true);
  TEST_EQUAL(hit.metaValueExists("_ar_1_score_type"), false);
  
  TEST_EQUAL(hit.getMetaValue("_ar_0_score_type").toString(), "mascot");
  TEST_REAL_SIMILAR(hit.getMetaValue("_ar_0_score"), 100.0);
  TEST_EQUAL(hit.getMetaValue("_ar_0_higher_is_better").toBool(), true);
}
END_SECTION

START_SECTION((void addAnalysisResults(const PepXMLAnalysisResult& aresult)))
{
  PeptideHit hit;
  
  // Add first analysis result
  PeptideHit::PepXMLAnalysisResult ar1;
  ar1.score_type = "peptideprophet";
  ar1.higher_is_better = true;
  ar1.main_score = 0.95;
  ar1.sub_scores["fval"] = 0.7114;
  
  hit.addAnalysisResults(ar1);
  
  TEST_EQUAL(hit.getAnalysisResults().size(), 1);
  TEST_EQUAL(hit.getAnalysisResults()[0].score_type, "peptideprophet");
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[0].main_score, 0.95);
  
  // Add second analysis result
  PeptideHit::PepXMLAnalysisResult ar2;
  ar2.score_type = "interprophet";
  ar2.higher_is_better = true;
  ar2.main_score = 0.98;
  ar2.sub_scores["nrs"] = 10.2137;
  
  hit.addAnalysisResults(ar2);
  
  TEST_EQUAL(hit.getAnalysisResults().size(), 2);
  TEST_EQUAL(hit.getAnalysisResults()[1].score_type, "interprophet");
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[1].main_score, 0.98);
  
  // Check meta values
  TEST_EQUAL(hit.metaValueExists("_ar_0_score_type"), true);
  TEST_EQUAL(hit.metaValueExists("_ar_1_score_type"), true);
  
  TEST_EQUAL(hit.getMetaValue("_ar_0_score_type").toString(), "peptideprophet");
  TEST_EQUAL(hit.getMetaValue("_ar_1_score_type").toString(), "interprophet");
}
END_SECTION

START_SECTION((const std::vector<PepXMLAnalysisResult>& getAnalysisResults() const))
{
  PeptideHit hit;
  
  // Test empty analysis results
  TEST_EQUAL(hit.getAnalysisResults().size(), 0);
  
  // Add analysis results
  PeptideHit::PepXMLAnalysisResult ar1;
  ar1.score_type = "peptideprophet";
  ar1.higher_is_better = true;
  ar1.main_score = 0.95;
  
  hit.addAnalysisResults(ar1);
  
  TEST_EQUAL(hit.getAnalysisResults().size(), 1);
  TEST_EQUAL(hit.getAnalysisResults()[0].score_type, "peptideprophet");
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[0].main_score, 0.95);
}
END_SECTION

START_SECTION((PeptideHit(const PeptideHit& source) - with analysis results))
{
  PeptideHit source;
  
  // Add analysis results to source
  PeptideHit::PepXMLAnalysisResult ar1;
  ar1.score_type = "peptideprophet";
  ar1.higher_is_better = true;
  ar1.main_score = 0.95;
  ar1.sub_scores["fval"] = 0.7114;
  
  source.addAnalysisResults(ar1);
  
  // Copy construct
  PeptideHit hit(source);
  
  // Check that analysis results were copied
  TEST_EQUAL(hit.getAnalysisResults().size(), 1);
  TEST_EQUAL(hit.getAnalysisResults()[0].score_type, "peptideprophet");
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[0].main_score, 0.95);
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[0].sub_scores.at("fval"), 0.7114);
  
  // Check meta values
  TEST_EQUAL(hit.metaValueExists("_ar_0_score_type"), true);
  TEST_EQUAL(hit.getMetaValue("_ar_0_score_type").toString(), "peptideprophet");
}
END_SECTION

START_SECTION((PeptideHit& operator=(const PeptideHit& source) - with analysis results))
{
  PeptideHit source;
  
  // Add analysis results to source
  PeptideHit::PepXMLAnalysisResult ar1;
  ar1.score_type = "peptideprophet";
  ar1.higher_is_better = true;
  ar1.main_score = 0.95;
  ar1.sub_scores["fval"] = 0.7114;
  
  source.addAnalysisResults(ar1);
  
  // Assignment
  PeptideHit hit;
  hit = source;
  
  // Check that analysis results were copied
  TEST_EQUAL(hit.getAnalysisResults().size(), 1);
  TEST_EQUAL(hit.getAnalysisResults()[0].score_type, "peptideprophet");
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[0].main_score, 0.95);
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[0].sub_scores.at("fval"), 0.7114);
  
  // Check meta values
  TEST_EQUAL(hit.metaValueExists("_ar_0_score_type"), true);
  TEST_EQUAL(hit.getMetaValue("_ar_0_score_type").toString(), "peptideprophet");
}
END_SECTION

START_SECTION((bool operator==(const PeptideHit& rhs) const - with analysis results))
{
  PeptideHit hit1, hit2;
  
  // Empty hits should be equal
  TEST_EQUAL(hit1 == hit2, true);
  
  // Add analysis results to hit1
  PeptideHit::PepXMLAnalysisResult ar1;
  ar1.score_type = "peptideprophet";
  ar1.higher_is_better = true;
  ar1.main_score = 0.95;
  
  hit1.addAnalysisResults(ar1);
  
  // Now they should be different
  TEST_EQUAL(hit1 == hit2, false);
  
  // Add same analysis results to hit2
  hit2.addAnalysisResults(ar1);
  
  // Now they should be equal again
  TEST_EQUAL(hit1 == hit2, true);
  
  // Add different analysis results to hit2
  PeptideHit::PepXMLAnalysisResult ar2;
  ar2.score_type = "interprophet";
  ar2.higher_is_better = true;
  ar2.main_score = 0.98;
  
  hit2.addAnalysisResults(ar2);
  
  // Now they should be different again
  TEST_EQUAL(hit1 == hit2, false);
}
END_SECTION

START_SECTION((PeptideHit(PeptideHit&& source) noexcept - with analysis results))
{
  PeptideHit source;
  
  // Add analysis results to source
  PeptideHit::PepXMLAnalysisResult ar1;
  ar1.score_type = "peptideprophet";
  ar1.higher_is_better = true;
  ar1.main_score = 0.95;
  ar1.sub_scores["fval"] = 0.7114;
  
  source.addAnalysisResults(ar1);
  
  // Move construct
  PeptideHit hit(std::move(source));
  
  // Check that analysis results were moved
  TEST_EQUAL(hit.getAnalysisResults().size(), 1);
  TEST_EQUAL(hit.getAnalysisResults()[0].score_type, "peptideprophet");
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[0].main_score, 0.95);
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[0].sub_scores.at("fval"), 0.7114);
  
  // Check meta values
  TEST_EQUAL(hit.metaValueExists("_ar_0_score_type"), true);
  TEST_EQUAL(hit.getMetaValue("_ar_0_score_type").toString(), "peptideprophet");
  
  // Source should have no meta values after move
  TEST_EQUAL(source.metaValueExists("_ar_0_score_type"), false);
  TEST_EQUAL(source.getAnalysisResults().size(), 0);
}
END_SECTION

START_SECTION((PeptideHit& operator=(PeptideHit&& source) noexcept - with analysis results))
{
  PeptideHit source;
  
  // Add analysis results to source
  PeptideHit::PepXMLAnalysisResult ar1;
  ar1.score_type = "peptideprophet";
  ar1.higher_is_better = true;
  ar1.main_score = 0.95;
  ar1.sub_scores["fval"] = 0.7114;
  
  source.addAnalysisResults(ar1);
  
  // Move assignment
  PeptideHit hit;
  hit = std::move(source);
  
  // Check that analysis results were moved
  TEST_EQUAL(hit.getAnalysisResults().size(), 1);
  TEST_EQUAL(hit.getAnalysisResults()[0].score_type, "peptideprophet");
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[0].main_score, 0.95);
  TEST_REAL_SIMILAR(hit.getAnalysisResults()[0].sub_scores.at("fval"), 0.7114);
  
  // Check meta values
  TEST_EQUAL(hit.metaValueExists("_ar_0_score_type"), true);
  TEST_EQUAL(hit.getMetaValue("_ar_0_score_type").toString(), "peptideprophet");
  
  // Source should have no meta values after move
  TEST_EQUAL(source.metaValueExists("_ar_0_score_type"), false);
  TEST_EQUAL(source.getAnalysisResults().size(), 0);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
