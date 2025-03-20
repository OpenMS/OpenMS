// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <string>

#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

///////////////////////////

START_TEST(ProteinIdentification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

float protein_significance_threshold = 63.2f;
std::vector<ProteinHit> protein_hits;
ProteinHit protein_hit;
ProteinIdentification protein_identification;
MascotXMLFile xml_file;

protein_hits.push_back(protein_hit);

ProteinIdentification* ptr = nullptr;
ProteinIdentification* nullPointer = nullptr;

START_SECTION((ProteinIdentification()))
	ptr = new ProteinIdentification();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION


START_SECTION((virtual ~ProteinIdentification()))
	ProteinIdentification hits;
	delete ptr;
END_SECTION

START_SECTION((ProteinIdentification(const ProteinIdentification &source)))
	ProteinIdentification hits;
	hits.setDateTime(DateTime::now());
	hits.setSignificanceThreshold(protein_significance_threshold);
	hits.insertHit(protein_hit);
	hits.setMetaValue("label",17);
	hits.setIdentifier("id");
	hits.setScoreType("score_type");
	hits.setHigherScoreBetter(false);
	hits.setSearchEngine("Mascot");
	hits.setSearchEngineVersion("2.1");
	ProteinIdentification::SearchParameters param;
	param.db = "RefSeq";
	ProteinIdentification::ProteinGroup g;
	g.probability = 0.99;
	g.accessions.push_back("protein0");
	g.accessions.push_back("protein3");
	hits.insertProteinGroup(g);
	hits.insertProteinGroup(g);
	hits.setSearchParameters(param);

	ProteinIdentification hits2(hits);

	TEST_EQUAL(hits.getDateTime() == hits2.getDateTime(), true)
	TEST_EQUAL(hits.getSignificanceThreshold(), hits2.getSignificanceThreshold())
	TEST_EQUAL(hits.getHits().size() == 1, true)
	TEST_EQUAL(hits.getHits()[0].getSequence(), String(""))
	TEST_EQUAL(hits.getHits()[0] == protein_hit, true)
	TEST_EQUAL(hits.getProteinGroups().size() == 2, true)
	TEST_EQUAL(hits.getProteinGroups()[0] == g, true)
	TEST_EQUAL(hits.getProteinGroups()[1] == g, true)
	TEST_EQUAL((UInt)hits.getMetaValue("label"), 17)
	TEST_EQUAL(hits.getIdentifier(), "id")
	TEST_EQUAL(hits.getScoreType(), "score_type")
	TEST_EQUAL(hits.isHigherScoreBetter(), false)
	TEST_EQUAL(hits.getSearchEngine(), "Mascot")
	TEST_EQUAL(hits.getSearchEngineVersion(), "2.1")
	TEST_EQUAL(hits.getSearchParameters() == param, true)
END_SECTION


START_SECTION((ProteinIdentification& operator=(const ProteinIdentification& source)))
	ProteinIdentification hits;
	hits.setDateTime(DateTime::now());
	hits.setSignificanceThreshold(protein_significance_threshold);
	hits.insertHit(protein_hit);
	hits.setIdentifier("id");
	hits.setScoreType("score_type");
	hits.setHigherScoreBetter(false);
	hits.setSearchEngine("Mascot");
	hits.setSearchEngineVersion("2.1");
	ProteinIdentification::ProteinGroup g;
	g.probability = 0.99;
	g.accessions.push_back("protein0");
	g.accessions.push_back("protein3");
	hits.insertProteinGroup(g);
	ProteinIdentification::SearchParameters param;
	param.db = "RefSeq";
	hits.setSearchParameters(param);

	ProteinIdentification hits2;
	hits2 = hits;

	TEST_EQUAL(hits.getDateTime() == hits2.getDateTime(), true)
	TEST_EQUAL(hits.getSignificanceThreshold(), hits2.getSignificanceThreshold())
	TEST_EQUAL(hits2.getHits().size() == 1, true)
	TEST_EQUAL(*(hits2.getHits().begin()) == protein_hit, true)
	TEST_EQUAL(hits2.getProteinGroups().size() == 1, true)
	TEST_EQUAL(hits2.getProteinGroups()[0] == g, true)
	TEST_EQUAL(hits2.getIdentifier(), "id")
	TEST_EQUAL(hits2.getScoreType(), "score_type")
	TEST_EQUAL(hits2.isHigherScoreBetter(), false)
	TEST_EQUAL(hits2.getSearchEngine(), "Mascot")
	TEST_EQUAL(hits2.getSearchEngineVersion(), "2.1")
	TEST_EQUAL(hits2.getSearchParameters() == param, true)
END_SECTION


START_SECTION((bool operator==(const ProteinIdentification& rhs) const))
	ProteinIdentification search1;
	ProteinIdentification search2;
	TEST_TRUE(search1 == search2)

	search1.setDateTime(DateTime::now());
	TEST_EQUAL(search1 == search2, false)
	search1 = search2;

	search1.setSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(search1 == search2, false)
	search1 = search2;

	search2.setIdentifier("id");
	TEST_EQUAL(search1 == search2, false)
	search1 = search2;

	search2.setScoreType("score_type");
	TEST_EQUAL(search1 == search2, false)
	search1 = search2;

	search2.setHigherScoreBetter(false);
	TEST_EQUAL(search1 == search2, false)
	search1 = search2;

	search2.setSearchEngine("Mascot");
	TEST_EQUAL(search1 == search2, false)
	search1 = search2;

	search2.setSearchEngineVersion("2.1");
	TEST_EQUAL(search1 == search2, false)
	search1 = search2;

	ProteinIdentification::SearchParameters param;
	param.db = "RefSeq";
	search2.setSearchParameters(param);
	TEST_EQUAL(search1 == search2, false)
	search1 = search2;

	ProteinIdentification::ProteinGroup g;
	g.probability = 0.99;
	g.accessions.push_back("protein0");
	g.accessions.push_back("protein3");
	search2.insertProteinGroup(g);
	TEST_EQUAL(search1 == search2, false)
	search1 = search2;

END_SECTION


START_SECTION((bool operator!=(const ProteinIdentification& rhs) const))
	ProteinIdentification search1;
	ProteinIdentification search2;
	TEST_EQUAL(search1 != search2, false)

	search1.setDateTime(DateTime::now());
	TEST_FALSE(search1 == search2)

	//rest does not need to be tested, as it is tested in the operator== test implicitly!
END_SECTION


START_SECTION((const DateTime& getDateTime() const))
	ProteinIdentification hits;
	DateTime date = DateTime::now();
	hits.setDateTime(date);
	const DateTime& date_time = hits.getDateTime();
	TEST_TRUE(date_time == date)
END_SECTION


START_SECTION((double getSignificanceThreshold() const))
	ProteinIdentification hits;
	hits.setSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(hits.getSignificanceThreshold(), protein_significance_threshold)
END_SECTION


START_SECTION((const std::vector<ProteinHit>& getHits() const))
	ProteinIdentification hits;
	hits.insertHit(protein_hit);
	TEST_EQUAL(hits.getHits().size() == 1, true)
	TEST_EQUAL(*(hits.getHits().begin()) == protein_hit, true)
END_SECTION


START_SECTION((std::vector<ProteinHit>& getHits()))
	ProteinIdentification hits;
	hits.insertHit(protein_hit);
	TEST_EQUAL(hits.getHits().size() == 1, true)
	TEST_EQUAL(*(hits.getHits().begin()) == protein_hit, true)
END_SECTION


START_SECTION((void insertHit(const ProteinHit& input)))
	ProteinIdentification hits;
	hits.insertHit(protein_hit);
	TEST_EQUAL(hits.getHits().size() == 1, true)
	TEST_EQUAL(*(hits.getHits().begin()) == protein_hit, true)
END_SECTION


START_SECTION((void setDateTime(const DateTime& date)))
	ProteinIdentification hits;
	DateTime date = DateTime::now();
	hits.setDateTime(date);
	TEST_EQUAL(hits.getDateTime() == date, true)
END_SECTION


START_SECTION((void setSignificanceThreshold(double value)))
	ProteinIdentification hits;
	hits.setSignificanceThreshold(protein_significance_threshold);
	TEST_EQUAL(hits.getSignificanceThreshold(), protein_significance_threshold)
END_SECTION


START_SECTION((void setHits(const std::vector<ProteinHit> &hits)))
	ProteinHit hit_1;
	ProteinHit hit_2;
	ProteinHit hit_3;
	vector<ProteinHit> hits;
	ProteinIdentification id;

	hit_1.setScore(23);
	hit_2.setScore(11);
	hit_3.setScore(45);
	hit_1.setAccession("SECONDPROTEIN");
	hit_2.setAccession("THIRDPROTEIN");
	hit_3.setAccession("FIRSTPROTEIN");
	hits.push_back(hit_1);
	hits.push_back(hit_2);
	hits.push_back(hit_3);
	id.setHits(hits);
	TEST_EQUAL(id.getHits()[2].getAccession(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getHits()[0].getAccession(), "SECONDPROTEIN")
	TEST_EQUAL(id.getHits()[1].getAccession(), "THIRDPROTEIN")
END_SECTION


START_SECTION((const String& getScoreType() const))
	ProteinIdentification hits;
	TEST_EQUAL(hits.getScoreType(), "")
END_SECTION


START_SECTION((void setScoreType(const String& type)))
	ProteinIdentification hits;
	hits.setScoreType("bla");
	TEST_EQUAL(hits.getScoreType(), "bla")
END_SECTION


START_SECTION((bool isHigherScoreBetter() const))
	ProteinIdentification hits;
	TEST_EQUAL(hits.isHigherScoreBetter(), true)
END_SECTION


START_SECTION((void setHigherScoreBetter(bool higher_is_better)))
	ProteinIdentification hits;
	hits.setHigherScoreBetter(false);
	TEST_EQUAL(hits.isHigherScoreBetter(), false)
END_SECTION


START_SECTION((const String& getIdentifier() const))
	ProteinIdentification hits;
	TEST_EQUAL(hits.getIdentifier(), "")
END_SECTION


START_SECTION((void setIdentifier(const String& id)))
	ProteinIdentification hits;
	hits.setIdentifier("bla");
	TEST_EQUAL(hits.getIdentifier(), "bla")
END_SECTION


START_SECTION((const String& getSearchEngine() const))
	ProteinIdentification hits;
	TEST_EQUAL(hits.getSearchEngine(), "")
END_SECTION


START_SECTION((void setSearchEngine(const String& search_engine)))
	ProteinIdentification hits;
	hits.setIdentifier("bla");
	TEST_EQUAL(hits.getIdentifier(), "bla")
END_SECTION


START_SECTION((const String& getSearchEngineVersion() const))
	ProteinIdentification hits;
	TEST_EQUAL(hits.getSearchEngineVersion(), "")
END_SECTION


START_SECTION((void setSearchEngineVersion(const String &search_engine_version)))
	ProteinIdentification hits;
	hits.setSearchEngineVersion("bla");
	TEST_EQUAL(hits.getSearchEngineVersion(), "bla")
END_SECTION


START_SECTION((const SearchParameters& getSearchParameters() const))
	ProteinIdentification hits;
	TEST_EQUAL(hits.getSearchParameters() == ProteinIdentification::SearchParameters(), true)
END_SECTION


START_SECTION((void setSearchParameters(const SearchParameters& search_parameters)))
	ProteinIdentification hits;
	ProteinIdentification::SearchParameters param;
	param.db = "Mascot";
	hits.setSearchParameters(param);
	TEST_EQUAL(hits.getSearchParameters() == ProteinIdentification::SearchParameters(), false)
END_SECTION


START_SECTION((void sort()))
{
	ProteinIdentification id;
	ProteinHit hit;
	hit.setScore(23);
	hit.setAccession("SECONDPROTEIN");
	id.insertHit(hit);
	hit.setScore(45);
	hit.setAccession("FIRSTPROTEIN");
	id.insertHit(hit);
	hit.setScore(7);
	hit.setAccession("THIRDPROTEIN");
	id.insertHit(hit);

	//higher score is better
	id.sort();

	TEST_EQUAL(id.getHits()[0].getAccession(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getHits()[1].getAccession(), "SECONDPROTEIN")
	TEST_EQUAL(id.getHits()[2].getAccession(), "THIRDPROTEIN")
	TEST_EQUAL(id.getHits()[0].getScore(), 45)
	TEST_EQUAL(id.getHits()[1].getScore(), 23)
	TEST_EQUAL(id.getHits()[2].getScore(), 7)

	//lower score is better
	id.setHigherScoreBetter(false);
	id.sort();

	TEST_EQUAL(id.getHits()[0].getAccession(), "THIRDPROTEIN")
	TEST_EQUAL(id.getHits()[1].getAccession(), "SECONDPROTEIN")
	TEST_EQUAL(id.getHits()[2].getAccession(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getHits()[0].getScore(), 7)
	TEST_EQUAL(id.getHits()[1].getScore(), 23)
	TEST_EQUAL(id.getHits()[2].getScore(), 45)
}
{
	ProteinIdentification id;
	ProteinHit hit;
	hit.setScore(45);
	hit.setAccession("SECONDPROTEIN");
	id.insertHit(hit);
	hit.setScore(7);
	hit.setAccession("FOURTHPROTEIN");
	id.insertHit(hit);
	hit.setScore(23);
	hit.setAccession("THIRDPROTEIN");
	id.insertHit(hit);
	hit.setScore(99);
	hit.setAccession("FIRSTPROTEIN");
	id.insertHit(hit);

	ProteinIdentification::ProteinGroup g1, g2;
	g1.probability = 0.99;
	g1.accessions.push_back("FIRSTPROTEIN");
	g1.accessions.push_back("FOURTHPROTEIN");
	id.insertProteinGroup(g1);
	g2.probability = 0.96;
	g2.accessions.push_back("FIRSTPROTEIN");
	g2.accessions.push_back("SECONDPROTEIN");
	g2.accessions.push_back("THIRDPROTEIN");
	id.insertProteinGroup(g2);

	//higher score is better
	id.sort();

	TEST_EQUAL(id.getHits()[0].getAccession(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getHits()[1].getAccession(), "SECONDPROTEIN")
	TEST_EQUAL(id.getHits()[2].getAccession(), "THIRDPROTEIN")
	TEST_EQUAL(id.getHits()[0].getScore(), 99)
	TEST_EQUAL(id.getHits()[1].getScore(), 45)
	TEST_EQUAL(id.getHits()[2].getScore(), 23)

	TEST_EQUAL(id.getProteinGroups().size(), 2);
	TEST_EQUAL(id.getProteinGroups()[0] == g1, true);
	TEST_EQUAL(id.getProteinGroups()[1] == g2, true);
}

END_SECTION


START_SECTION((void assignRanks()))
	ProteinIdentification id;
	ProteinHit hit;
	hit.setScore(23);
	hit.setAccession("SECONDPROTEIN");
	id.insertHit(hit);
	hit.setScore(45);
	hit.setAccession("FIRSTPROTEIN");
	id.insertHit(hit);
	hit.setScore(7);
	hit.setAccession("THIRDPROTEIN");
	id.insertHit(hit);

	id.assignRanks();

	TEST_EQUAL(id.getHits()[0].getAccession(), "FIRSTPROTEIN")
	TEST_EQUAL(id.getHits()[1].getAccession(), "SECONDPROTEIN")
	TEST_EQUAL(id.getHits()[2].getAccession(), "THIRDPROTEIN")
	TEST_EQUAL(id.getHits()[0].getRank(), 1)
	TEST_EQUAL(id.getHits()[1].getRank(), 2)
	TEST_EQUAL(id.getHits()[2].getRank(), 3)
END_SECTION


START_SECTION((Size computeCoverage(const std::vector<PeptideIdentification>& pep_ids)))
	ProteinIdentification id;

  // prep hit
  ProteinHit hit;
  hit.setAccession("P1");
  hit.setSequence("MKQSTIALALLPLLFTPVTKARTPEMPVLENRAAQGDITAPGGARRLTGDQTAALRDSLS"
                  "DKPAKNIILLIGDGMGDSEITAARNYAEGAGGFFKGIDALPLTGQYTHYALNKKTGKPDY"
                  "VTDSAASATAWSTGVKTYNGALGVDIHEKDHPTILEMAKAAGLATGNVSTAELQDATPAA");
  id.insertHit(hit);
	hit.setAccession("P2");
  hit.setSequence("PEMPVLENRAAQGDITAPPGGARRLTGDQTAALRDSLS");
  id.insertHit(hit);

  // prep peptides
  std::vector<PeptideIdentification> pep_ids;
  PeptideIdentification pid;
  PeptideHit phit(0, 0, 1, AASequence::fromString(""));
  PeptideEvidence pe;
  pe.setProteinAccession("P1");
  pe.setStart(0);
  pe.setEnd(59);
  phit.addPeptideEvidence(pe);
  phit.setSequence(AASequence::fromString("MKQSTIALALLPLLFTPVTKARTPEMPVLENRAAQGDITAPGGARRLTGDQTAALRDSLS"));
  pid.insertHit(phit);
  pe.setStart(60);
  pe.setEnd(119);
  phit.addPeptideEvidence(pe);
  phit.setSequence(AASequence::fromString("DKPAKNIILLIGDGMGDSEITAARNYAEGAGGFFKGIDALPLTGQYTHYALNKKTGKPDY"));
  pid.insertHit(phit);
  pe.setStart(0);
  pe.setEnd(59);
  phit.addPeptideEvidence(pe);
  phit.setSequence(AASequence::fromString("MKQSTIALALLPLLFTPVTKARTPEMPVLENRAAQGDITAPGGARRLTGDQTAALRDSLS")); // should not count
  pid.insertHit(phit);
  pep_ids.push_back(pid);

  PeptideIdentification pid2;
  PeptideHit phit2(0, 0, 1, AASequence::fromString(""));
  pe.setStart(0);
  pe.setEnd(59);
  phit2.addPeptideEvidence(pe);
  phit2.setSequence(AASequence::fromString("MKQSTIALALLPLLFTPVTKARTPEMPVLENRAAQGDITAPGGARRLTGDQTAALRDSLS"));
  pid2.insertHit(phit2); // should not count
  pep_ids.push_back(pid2);

  id.computeCoverage(pep_ids);

  TEST_REAL_SIMILAR(id.getHits()[0].getCoverage(), 200.0 / 3.0);
  TEST_REAL_SIMILAR(id.getHits()[1].getCoverage(), 0.0);

  pe.setStart(120);
  pe.setEnd(179);
  phit2.addPeptideEvidence(pe);
  phit2.setSequence(AASequence::fromString("VTDSAASATAWSTGVKTYNGALGVDIHEKDHPTILEMAKAAGLATGNVSTAELQDATPAA"));
  pid2.insertHit(phit2);
  pep_ids.push_back(pid2);

  id.computeCoverage(pep_ids);

  TEST_REAL_SIMILAR(id.getHits()[0].getCoverage(), 100.0);
  TEST_REAL_SIMILAR(id.getHits()[1].getCoverage(), 0.0);

  pep_ids.clear();
  PeptideIdentification pid3;
  PeptideHit phit3(0, 0, 1, AASequence::fromString(""));
  PeptideEvidence pe2;
  pe2.setProteinAccession("P2");
  pe2.setStart(0);
  pe2.setEnd(18);
  phit3.addPeptideEvidence(pe2);
  phit3.setSequence(AASequence::fromString("PEMPVLENRAAQGDITAPP")); // 1st half
  pid3.insertHit(phit3);
  pe2.setStart(19);
  pe2.setEnd(37);
  phit3.addPeptideEvidence(pe2);
  phit3.setSequence(AASequence::fromString("GGARRLTGDQTAALRDSLS")); // 2nd half
  pid3.insertHit(phit3);
  pe2.setStart(8);
  pe2.setEnd(26);
  phit3.addPeptideEvidence(pe2);
  phit3.setSequence(AASequence::fromString("RAAQGDITAPPGGARRLTG")); // middle half
  pid3.insertHit(phit3);

  pep_ids.push_back(pid3);

  id.computeCoverage(pep_ids);

  TEST_REAL_SIMILAR(id.getHits()[0].getCoverage(), 0.0);
  TEST_REAL_SIMILAR(id.getHits()[1].getCoverage(), 100.0);

END_SECTION


START_SECTION(([ProteinIdentification::ProteinGroup] ProteinGroup()))
  ProteinIdentification::ProteinGroup p;
  TEST_EQUAL(p.probability, 0)
  TEST_EQUAL(p.accessions.size(), 0)
END_SECTION


START_SECTION(([ProteinIdentification::ProteinGroup] bool operator==(const ProteinGroup& rhs) const))
  ProteinIdentification::ProteinGroup p, p_c;

  p.probability = 0.5;
  TEST_NOT_EQUAL(p == p_c, true)

  p.probability = 0.0;
  p.accessions.push_back("bla");
  TEST_NOT_EQUAL(p == p_c, true)

  p_c = p;
  TEST_TRUE(p == p_c)
END_SECTION


START_SECTION(([ProteinIdentification::ProteinGroup] bool operator<(const ProteinGroup& rhs) const))
{
  ProteinIdentification::ProteinGroup p1, p2;

  // both are equal:
  TEST_EQUAL(p1 < p2, false);
  TEST_EQUAL(p2 < p1, false);

  // different probabilities:
  p1.probability = 0.1;
  p2.probability = 0.2;
  TEST_EQUAL(p2 < p1, true); // yes! (see documentation)
  TEST_EQUAL(p1 < p2, false);

  // equal again:
  p2.probability = 0.1;
  p1.accessions.push_back("bla");
  p2.accessions.push_back("bla");
  TEST_EQUAL(p1 < p2, false);
  TEST_EQUAL(p2 < p1, false);

  // different numbers of accessions:
  p2.accessions.push_back("blubb");
  TEST_EQUAL(p1 < p2, true);
  TEST_EQUAL(p2 < p1, false);

  // different accessions:
  p1.accessions.push_back("laber");
  TEST_EQUAL(p1 < p2, false);
  TEST_EQUAL(p2 < p1, true);
}
END_SECTION


START_SECTION(([ProteinIdentification::SearchParameters] SearchParameters()))
  ProteinIdentification::SearchParameters sp;

  TEST_EQUAL(sp.db.size(), 0)
  TEST_EQUAL(sp.db_version.size(), 0)
  TEST_EQUAL(sp.taxonomy.size(), 0)
  TEST_EQUAL(sp.charges.size(), 0)
  TEST_EQUAL(sp.mass_type, 0)
  TEST_EQUAL(sp.fixed_modifications.size(), 0)
  TEST_EQUAL(sp.variable_modifications.size(), 0)
  TEST_EQUAL(sp.digestion_enzyme.getName(), "unknown_enzyme")
  TEST_EQUAL(sp.missed_cleavages, 0)
  TEST_EQUAL(sp.fragment_mass_tolerance, 0.0)
  TEST_EQUAL(sp.fragment_mass_tolerance_ppm, false)
  TEST_EQUAL(sp.precursor_mass_tolerance, 0.0)
  TEST_EQUAL(sp.precursor_mass_tolerance_ppm, false)

END_SECTION


START_SECTION(([ProteinIdentification::SearchParameters] bool operator==(const SearchParameters& rhs) const))
  ProteinIdentification::SearchParameters sp, sp_n;
  sp_n.charges = "1,2,3";
  TEST_EQUAL(sp == sp_n, false)
END_SECTION


START_SECTION(([ProteinIdentification::SearchParameters] bool operator!=(const SearchParameters& rhs) const))
  ProteinIdentification::SearchParameters sp, sp_n;
  sp_n.charges = "1,2,3";
  TEST_FALSE(sp == sp_n)
END_SECTION


START_SECTION(([ProteinIdentification::SearchParameters] pair<int,int> getChargeRange() const))
{
  ProteinIdentification::SearchParameters sp;
  sp.charges = "1,2,3";
  auto range = sp.getChargeRange();
  TEST_EQUAL(range.first, 1);
  TEST_EQUAL(range.second, 3);
  sp.charges = "+2-+5";
  range = sp.getChargeRange();
  TEST_EQUAL(range.first, 2);
  TEST_EQUAL(range.second, 5);
  sp.charges = "-1,-2,-3";
  range = sp.getChargeRange();
  TEST_EQUAL(range.first, -3);
  TEST_EQUAL(range.second, -1);
  sp.charges = "2";
  range = sp.getChargeRange();
  TEST_EQUAL(range.first, 2);
  TEST_EQUAL(range.second, 2);
}
END_SECTION


START_SECTION((const vector<ProteinGroup>& getProteinGroups() const))
	ProteinIdentification id;
	ProteinIdentification::ProteinGroup g;
	g.probability = 0.99;
	g.accessions.push_back("protein0");
	g.accessions.push_back("protein3");
	id.insertProteinGroup(g);

	TEST_EQUAL(id.getProteinGroups().size(), 1);
	TEST_EQUAL(id.getProteinGroups()[0] == g, true);
END_SECTION


START_SECTION((vector<ProteinGroup>& getProteinGroups()))
	ProteinIdentification id;
	ProteinIdentification::ProteinGroup g;
	g.probability = 0.99;
	g.accessions.push_back("protein0");
	g.accessions.push_back("protein3");
	id.insertProteinGroup(g);

	TEST_EQUAL(id.getProteinGroups().size(), 1);
	TEST_EQUAL(id.getProteinGroups()[0] == g, true);
  TEST_EQUAL(id.getProteinGroups()[0].probability, 0.99);
END_SECTION


START_SECTION((void insertProteinGroup(const ProteinGroup& group)))
  NOT_TESTABLE
	//tested above
END_SECTION


START_SECTION((const vector<ProteinGroup>& getIndistinguishableProteins() const))
	ProteinIdentification id;
	ProteinIdentification::ProteinGroup g;
	g.accessions.push_back("protein0");
	g.accessions.push_back("protein3");
	id.insertIndistinguishableProteins(g);

	TEST_EQUAL(id.getIndistinguishableProteins().size(), 1);
	TEST_EQUAL(id.getIndistinguishableProteins()[0] == g, true);
END_SECTION


START_SECTION((vector<ProteinGroup>& getIndistinguishableProteins()))
	ProteinIdentification id;
	ProteinIdentification::ProteinGroup g;
	g.probability = 0.99;
	g.accessions.push_back("protein0");
	g.accessions.push_back("protein3");
	id.insertIndistinguishableProteins(g);

	TEST_EQUAL(id.getIndistinguishableProteins().size(), 1);
	TEST_EQUAL(id.getIndistinguishableProteins()[0] == g, true);
  id.getIndistinguishableProteins()[0].accessions.push_back("protein1");
	TEST_EQUAL(id.getIndistinguishableProteins()[0].accessions.size(), 3);
END_SECTION


START_SECTION((void insertIndistinguishableProteins(const ProteinGroup& group)))
	NOT_TESTABLE
	//tested above
END_SECTION


START_SECTION((vector<ProteinHit>::iterator findHit(const String& accession)))
{
	ProteinIdentification protein;
	ProteinHit hit;
	hit.setAccession("test1");
	protein.insertHit(hit);
	hit.setAccession("test2");
	protein.insertHit(hit);
	TEST_EQUAL(protein.findHit("test1")->getAccession(), "test1");
	TEST_EQUAL(protein.findHit("test2")->getAccession(), "test2");
	TEST_EQUAL(protein.findHit("test3") == protein.getHits().end(), true);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
