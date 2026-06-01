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

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/IdentifierMSRunMapper.h>
#include <OpenMS/METADATA/USI.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/String.h>

///////////////////////////

START_TEST(PeptideIdentification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


double peptide_significance_threshold = 42.3;
std::vector<PeptideHit> peptide_hits;
PeptideHit peptide_hit;
ProteinIdentification protein_identification;
vector<PeptideIdentification> identifications;
MascotXMLFile xml_file;

peptide_hits.push_back(peptide_hit);


PeptideIdentification* ptr = nullptr;
PeptideIdentification* nullPointer = nullptr;
START_SECTION((PeptideIdentification()))
  ptr = new PeptideIdentification();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeptideIdentification()))
  PeptideIdentification hits;
  delete ptr;
END_SECTION

// Copy Constructor
START_SECTION((PeptideIdentification(const PeptideIdentification& source)))
{
  PeptideIdentification hits;
  hits.setSignificanceThreshold(peptide_significance_threshold);
  hits.setHits(peptide_hits);
  hits.setMetaValue("label",17);
  hits.setIdentifier("id");
  hits.setScoreType("score_type");
  hits.setHigherScoreBetter(false);

  PeptideIdentification hits2(hits);

  TEST_EQUAL(hits.getSignificanceThreshold(), hits2.getSignificanceThreshold())
  TEST_TRUE(hits.getHits().size() == 1)
  TEST_TRUE(*(hits.getHits().begin()) == peptide_hit)
  TEST_EQUAL((UInt)hits.getMetaValue("label"),17)
  TEST_EQUAL(hits.getIdentifier(),"id")
  TEST_EQUAL(hits.getScoreType(),"score_type")
  TEST_FALSE(hits.isHigherScoreBetter())
}
END_SECTION

// Move Constructor
START_SECTION((PeptideIdentification(PeptideIdentification&& source) noexcept))
{
  // Ensure that PeptideIdentification has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_TRUE(noexcept(PeptideIdentification(std::declval<PeptideIdentification&&>())))

  PeptideIdentification hits;
  hits.setSignificanceThreshold(peptide_significance_threshold);
  hits.setHits(peptide_hits);
  hits.setMetaValue("label",17);
  hits.setIdentifier("id");
  hits.setScoreType("score_type");
  hits.setHigherScoreBetter(false);

  PeptideIdentification example(hits);

  PeptideIdentification hits2(std::move(example));

  TEST_EQUAL(hits.getSignificanceThreshold(), hits2.getSignificanceThreshold())
  TEST_TRUE(hits.getHits().size() == 1)
  TEST_TRUE(*(hits.getHits().begin()) == peptide_hit)
  TEST_EQUAL((UInt)hits.getMetaValue("label"),17)
  TEST_EQUAL(hits.getIdentifier(),"id")
  TEST_EQUAL(hits.getScoreType(),"score_type")
  TEST_FALSE(hits.isHigherScoreBetter())

  // the move source should be empty
  TEST_TRUE(example.getHits().empty())
  TEST_TRUE(example.isMetaEmpty())
  TEST_TRUE(example.getIdentifier().empty())
  TEST_TRUE(example.getScoreType().empty())
}
END_SECTION

START_SECTION((PeptideIdentification& operator=(const PeptideIdentification& source)))
  PeptideIdentification hits;
  hits.setSignificanceThreshold(peptide_significance_threshold);
  hits.setHits(peptide_hits);
  hits.setMetaValue("label",17);
  hits.setIdentifier("id");
  hits.setScoreType("score_type");
  hits.setHigherScoreBetter(false);

  PeptideIdentification hits2;
  hits2 = hits;

  TEST_EQUAL(hits.getSignificanceThreshold(), hits2.getSignificanceThreshold())
  TEST_TRUE(hits.getHits().size() == 1)
  TEST_TRUE(*(hits.getHits().begin()) == peptide_hit)
  TEST_EQUAL((UInt)hits.getMetaValue("label"),17)
  TEST_EQUAL(hits.getIdentifier(),"id")
  TEST_EQUAL(hits.getScoreType(),"score_type")
  TEST_FALSE(hits.isHigherScoreBetter())
END_SECTION

START_SECTION((bool operator == (const PeptideIdentification& rhs) const))
  PeptideIdentification search1, search2;
  TEST_TRUE(search1 == search2)

  search1.setSignificanceThreshold(peptide_significance_threshold);
  TEST_FALSE(search1 == search2)
  search1 = search2;

  search2.setMetaValue("label",17);
  TEST_FALSE(search1 == search2)
  search1 = search2;

  search2.setIdentifier("id");
  TEST_FALSE(search1 == search2)
  search1 = search2;

  search2.setScoreType("score_type");
  TEST_FALSE(search1 == search2)
  search1 = search2;

  search2.setHigherScoreBetter(false);
  TEST_FALSE(search1 == search2)
  search1 = search2;
END_SECTION


START_SECTION((bool operator != (const PeptideIdentification& rhs) const))
  PeptideIdentification search1, search2;
  TEST_FALSE(search1 != search2)

  search1.setSignificanceThreshold(peptide_significance_threshold);
  TEST_FALSE(search1 == search2)
  search1 = search2;

  //rest does not need to be tested, as it is tested in the operator== test implicitly!
END_SECTION


START_SECTION((double getRT() const))
	PeptideIdentification pi;
	TEST_FALSE(pi.hasRT());
	pi.setRT(1024.0);
	TEST_EQUAL(pi.getRT(), 1024.0);
END_SECTION

START_SECTION(void setRT(double mz))
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(bool hasRT())
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION((double getMZ() const))
	PeptideIdentification pi;
	TEST_FALSE(pi.hasMZ());
	pi.setMZ(1024.0);
	TEST_EQUAL(pi.getMZ(), 1024.0);
END_SECTION


START_SECTION(bool hasMZ())
  NOT_TESTABLE // tested above
END_SECTION


START_SECTION((double getSignificanceThreshold() const))
	PeptideIdentification hits;
	hits.setSignificanceThreshold(peptide_significance_threshold);
	TEST_EQUAL(hits.getSignificanceThreshold(), peptide_significance_threshold)
END_SECTION

START_SECTION((const std::vector<PeptideHit>& getHits() const))
  PeptideIdentification hits;
  hits.insertHit(peptide_hit);
  TEST_TRUE(hits.getHits().size() == 1)
  TEST_TRUE(hits.getHits()[0] == peptide_hit)
END_SECTION

START_SECTION((void insertHit(const PeptideHit &hit)))
  PeptideIdentification hits;
  hits.insertHit(peptide_hit);
  TEST_TRUE(hits.getHits().size() == 1)
  TEST_TRUE(*(hits.getHits().begin()) == peptide_hit)
END_SECTION

START_SECTION((void setHits(const std::vector< PeptideHit > &hits)))
  PeptideIdentification hits;
  hits.setHits(peptide_hits);
  TEST_TRUE(hits.getHits() == peptide_hits)
END_SECTION

START_SECTION((void setSignificanceThreshold(double value)))
  PeptideIdentification hits;
  hits.setSignificanceThreshold(peptide_significance_threshold);
  TEST_EQUAL(hits.getSignificanceThreshold(), peptide_significance_threshold)
END_SECTION

START_SECTION((void testSignificanceThresholdMetaValue()))
{
  PeptideIdentification id;
  
  // Test default value
  TEST_EQUAL(id.getSignificanceThreshold(), 0.0);
  TEST_EQUAL(id.metaValueExists(Constants::UserParam::SIGNIFICANCE_THRESHOLD), false);
  
  // Test setting a new value
  id.setSignificanceThreshold(42.0);
  TEST_EQUAL(id.getSignificanceThreshold(), 42.0);
  TEST_EQUAL(id.getMetaValue(Constants::UserParam::SIGNIFICANCE_THRESHOLD), 42.0);
  
  // Test that the meta value is copied
  PeptideIdentification id2(id);
  TEST_EQUAL(id2.getSignificanceThreshold(), 42.0);
  TEST_EQUAL(id2.getMetaValue(Constants::UserParam::SIGNIFICANCE_THRESHOLD), 42.0);
  
  // Test assignment operator
  PeptideIdentification id3;
  id3 = id;
  TEST_EQUAL(id3.getSignificanceThreshold(), 42.0);
  TEST_EQUAL(id3.getMetaValue(Constants::UserParam::SIGNIFICANCE_THRESHOLD), 42.0);
}
END_SECTION

START_SECTION((String& getScoreType() const))
	PeptideIdentification hits;
	TEST_EQUAL(hits.getScoreType(), "")
END_SECTION

START_SECTION((void setScoreType(const String& type)))
	PeptideIdentification hits;
	hits.setScoreType("bla");
	TEST_EQUAL(hits.getScoreType(), "bla")
END_SECTION

START_SECTION((bool isHigherScoreBetter() const))
	PeptideIdentification hits;
	TEST_TRUE(hits.isHigherScoreBetter())
END_SECTION

START_SECTION((void setHigherScoreBetter(bool value)))
  PeptideIdentification hits;
  hits.setHigherScoreBetter(false);
  TEST_FALSE(hits.isHigherScoreBetter())
END_SECTION

START_SECTION((const String& getIdentifier() const))
  PeptideIdentification hits;
  TEST_EQUAL(hits.getIdentifier(),"")
END_SECTION

START_SECTION((void setIdentifier(const String& id)))
  PeptideIdentification hits;
  hits.setIdentifier("bla");
  TEST_EQUAL(hits.getIdentifier(),"bla")
END_SECTION

START_SECTION((bool empty() const))
  PeptideIdentification hits;
  TEST_TRUE(hits.empty())
  hits.setSignificanceThreshold(1);
  TEST_FALSE(hits.empty())

  hits.setSignificanceThreshold(0);
  TEST_TRUE(hits.empty())

  hits.setBaseName("basename");
  TEST_FALSE(hits.empty())

  hits.setBaseName("");
  TEST_TRUE(hits.empty())

  hits.insertHit(peptide_hit);
  TEST_FALSE(hits.empty())
END_SECTION

START_SECTION((void sort()))
  PeptideIdentification id;
  PeptideHit hit;
  hit.setScore(23);
  hit.setSequence(AASequence::fromString("SECONDPROTEIN"));
  id.insertHit(hit);
  hit.setScore(45);
  hit.setSequence(AASequence::fromString("FIRSTPROTEIN"));
  id.insertHit(hit);
  hit.setScore(7);
  hit.setSequence(AASequence::fromString("THIRDPROTEIN"));
  id.insertHit(hit);

  //higher score is better
  id.sort();

  TEST_EQUAL(id.getHits()[0].getSequence(), AASequence::fromString("FIRSTPROTEIN"))
  TEST_EQUAL(id.getHits()[1].getSequence(), AASequence::fromString("SECONDPROTEIN"))
  TEST_EQUAL(id.getHits()[2].getSequence(), AASequence::fromString("THIRDPROTEIN"))
  TEST_EQUAL(id.getHits()[0].getScore(), 45)
  TEST_EQUAL(id.getHits()[1].getScore(), 23)
  TEST_EQUAL(id.getHits()[2].getScore(), 7)

  //lower score is better
  id.setHigherScoreBetter(false);
  id.sort();

  TEST_EQUAL(id.getHits()[0].getSequence(), AASequence::fromString("THIRDPROTEIN"))
  TEST_EQUAL(id.getHits()[1].getSequence(), AASequence::fromString("SECONDPROTEIN"))
  TEST_EQUAL(id.getHits()[2].getSequence(), AASequence::fromString("FIRSTPROTEIN"))
  TEST_EQUAL(id.getHits()[0].getScore(), 7)
  TEST_EQUAL(id.getHits()[1].getScore(), 23)
  TEST_EQUAL(id.getHits()[2].getScore(), 45)

END_SECTION


START_SECTION(static std::vector<PeptideHit> getReferencingHits(const std::vector<PeptideHit> & , const std::set<String> & accession))
{
  PeptideIdentification id;
  PeptideHit hit;
  vector< PeptideHit > peptide_hits;

  hit.setScore(23);
  hit.setSequence(AASequence::fromString("FIRSTPROTEIN"));
  PeptideEvidence pe;
  pe.setProteinAccession("TEST_PROTEIN1");
  hit.addPeptideEvidence(pe);
  id.insertHit(hit);

  hit = PeptideHit();
  hit.setScore(10);
  hit.setSequence(AASequence::fromString("SECONDPROTEIN"));
  pe.setProteinAccession("TEST_PROTEIN2");
  hit.addPeptideEvidence(pe);
  id.insertHit(hit);

  hit = PeptideHit();
  hit.setScore(11);
  hit.setSequence(AASequence::fromString("THIRDPROTEIN"));
  pe.setProteinAccession("TEST_PROTEIN2");
  hit.addPeptideEvidence(pe);
  id.insertHit(hit);

  set<String> query_accession;
  query_accession.insert("TEST_PROTEIN2");
  peptide_hits = PeptideIdentification::getReferencingHits(id.getHits(), query_accession);
  TEST_EQUAL(peptide_hits.size(), 2)
  TEST_EQUAL(peptide_hits[0].getSequence(), AASequence::fromString("SECONDPROTEIN"))
  TEST_EQUAL(peptide_hits[1].getSequence(), AASequence::fromString("THIRDPROTEIN"))

  query_accession.insert("TEST_PROTEIN3");
  peptide_hits = PeptideIdentification::getReferencingHits(id.getHits(), query_accession);
  TEST_EQUAL(peptide_hits.size(), 2)
  TEST_EQUAL(peptide_hits[0].getSequence(), AASequence::fromString("SECONDPROTEIN"))
  TEST_EQUAL(peptide_hits[1].getSequence(), AASequence::fromString("THIRDPROTEIN"))
}
{
  PeptideIdentification id;
  PeptideHit hit;
  vector< PeptideHit > peptide_hits;

  hit.setScore(23);
  hit.setSequence(AASequence::fromString("FIRSTPROTEIN"));
  PeptideEvidence pe;
  pe.setProteinAccession("TEST_PROTEIN1");
  hit.addPeptideEvidence(pe);
  id.insertHit(hit);

  hit = PeptideHit();
  hit.setScore(10);
  hit.setSequence(AASequence::fromString("SECONDPROTEIN"));
  pe.setProteinAccession("TEST_PROTEIN2");
  hit.addPeptideEvidence(pe);
  id.insertHit(hit);

  hit = PeptideHit();
  hit.setScore(11);
  hit.setSequence(AASequence::fromString("THIRDPROTEIN"));
  pe.setProteinAccession("TEST_PROTEIN3");
  hit.addPeptideEvidence(pe);
  id.insertHit(hit);

  set<String> query_accession;
  query_accession.insert("TEST_PROTEIN2");
  query_accession.insert("TEST_PROTEIN3");
  peptide_hits = PeptideIdentification::getReferencingHits(id.getHits(), query_accession);
  TEST_EQUAL(peptide_hits.size(), 2)
  TEST_EQUAL(peptide_hits[0].getSequence(), AASequence::fromString("SECONDPROTEIN"))
  TEST_EQUAL(peptide_hits[1].getSequence(), AASequence::fromString("THIRDPROTEIN"))
}
END_SECTION

START_SECTION((static std::multimap<String, std::pair<Size, Size>> buildUIDsFromAllPepIDs(const ConsensusMap &cmap)))
{
  ConsensusMap cmap;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_2.consensusXML"), cmap);
  std::multimap<String, std::pair<Size, Size>> map_of_UIDs = PeptideIdentification::buildUIDsFromAllPepIDs(cmap);

  auto b = map_of_UIDs.begin();
  TEST_EQUAL(b->first,
             "file:///C:/Users/bielow/AppData/Local/Temp/20190911_110348_8204_1/Untitled_workflow/002_FileFilter/out/ES-0014b_2.mzML|spectrum=112")
  TEST_EQUAL(b->second.first, Size(-1))
  TEST_EQUAL(b->second.second, 4)
  ++b;
  TEST_EQUAL(b->first,
             "file:///C:/Users/bielow/AppData/Local/Temp/20190911_110348_8204_1/Untitled_workflow/002_FileFilter/out/ES-0014b_2.mzML|spectrum=113")
  TEST_EQUAL(b->second.first, 11)
  TEST_EQUAL(b->second.second, 0)
  ++b;
  TEST_EQUAL(b->first,
             "file:///C:/Users/bielow/AppData/Local/Temp/20190911_110348_8204_1/Untitled_workflow/002_FileFilter/out/ES-0014b_2.mzML|spectrum=118")
  TEST_EQUAL(b->second.first, 4)
  TEST_EQUAL(b->second.second, 0)

  FeatureMap fmap;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_1.featureXML"), fmap);
  IdentifierMSRunMapper mapping(cmap.getProteinIdentifications());

  String uid_zero = PeptideIdentification::buildUIDFromPepID(cmap[0].getPeptideIdentifications()[0], mapping);
  String uid_one = PeptideIdentification::buildUIDFromPepID(cmap[0].getPeptideIdentifications()[1], mapping);
  String uid_two = PeptideIdentification::buildUIDFromPepID(cmap[0].getPeptideIdentifications()[2], mapping);

  TEST_EQUAL(uid_zero,"file:///C:/Users/bielow/AppData/Local/Temp/20190911_110348_8204_1/Untitled_workflow/002_FileFilter/out/ES-0014b_2.mzML|spectrum=219")
  TEST_EQUAL(uid_one, "file:///C:/Users/bielow/AppData/Local/Temp/20190911_110348_8204_1/Untitled_workflow/002_FileFilter/out/ES-0016_1.mzML|spectrum=33")
  TEST_EQUAL(uid_two, "file:///C:/Users/bielow/AppData/Local/Temp/20190911_110348_8204_1/Untitled_workflow/002_FileFilter/out/ES-0016_2_.mzML|spectrum=133")
}
END_SECTION

START_SECTION((static String buildUIDFromPepID(const PeptideIdentification& pep_id, const IdentifierMSRunMapper& mapping)))
{
      NOT_TESTABLE //Tested above
}
END_SECTION

START_SECTION((const String& getBaseName() const))
  PeptideIdentification id;
  TEST_EQUAL(id.getBaseName(), "")
  
  id.setBaseName("test_base_name");
  TEST_EQUAL(id.getBaseName(), "test_base_name")
  
  // Test that it's stored as a MetaValue
  TEST_TRUE(id.metaValueExists(Constants::UserParam::BASE_NAME))
  TEST_EQUAL(id.getMetaValue(Constants::UserParam::BASE_NAME), "test_base_name")
END_SECTION

START_SECTION((void setBaseName(const String& base_name)))
  PeptideIdentification id;
  
  // Test setting a non-empty base name
  id.setBaseName("test_base_name");
  TEST_EQUAL(id.getBaseName(), "test_base_name")
  TEST_TRUE(id.metaValueExists(Constants::UserParam::BASE_NAME))
  
  // Test setting an empty base name (should remove the meta value)
  id.setBaseName("");
  TEST_EQUAL(id.getBaseName(), "")
  TEST_FALSE(id.metaValueExists(Constants::UserParam::BASE_NAME))
  
  // Test that empty() method works correctly with base name as meta value
  TEST_TRUE(id.empty())
  id.setBaseName("test_base_name");
  TEST_FALSE(id.empty())
  id.setBaseName("");
  TEST_TRUE(id.empty())
END_SECTION

START_SECTION((std::hash<PeptideIdentification>))
{
  // Test 1: Equal objects have equal hashes
  PeptideIdentification id1, id2;
  id1.setRT(123.456);
  id1.setMZ(789.012);
  id1.setIdentifier("test_id");
  id1.setScoreType("Mascot");
  id1.setHigherScoreBetter(true);
  id1.setSignificanceThreshold(0.05);

  id2.setRT(123.456);
  id2.setMZ(789.012);
  id2.setIdentifier("test_id");
  id2.setScoreType("Mascot");
  id2.setHigherScoreBetter(true);
  id2.setSignificanceThreshold(0.05);

  TEST_EQUAL(id1 == id2, true)
  TEST_EQUAL(std::hash<PeptideIdentification>{}(id1), std::hash<PeptideIdentification>{}(id2))

  // Test 2: Different objects have different hashes (most likely)
  PeptideIdentification id3;
  id3.setRT(999.999);
  id3.setMZ(111.111);
  id3.setIdentifier("different_id");

  TEST_EQUAL(id1 == id3, false)
  TEST_NOT_EQUAL(std::hash<PeptideIdentification>{}(id1), std::hash<PeptideIdentification>{}(id3))

  // Test 3: Can use in unordered_set
  std::unordered_set<PeptideIdentification> id_set;
  id_set.insert(id1);
  id_set.insert(id2); // Should not increase size (equal to id1)
  id_set.insert(id3);
  TEST_EQUAL(id_set.size(), 2)

  // Test 4: Can use in unordered_map
  std::unordered_map<PeptideIdentification, std::string> id_map;
  id_map[id1] = "first";
  id_map[id3] = "third";
  TEST_EQUAL(id_map[id1], "first")
  TEST_EQUAL(id_map[id2], "first") // id2 == id1
  TEST_EQUAL(id_map[id3], "third")

  // Test 5: PeptideIdentification with peptide hits
  PeptideIdentification id4, id5;
  PeptideHit hit1;
  hit1.setScore(50.0);
  hit1.setSequence(AASequence::fromString("PEPTIDE"));
  hit1.setCharge(2);

  PeptideHit hit2;
  hit2.setScore(50.0);
  hit2.setSequence(AASequence::fromString("PEPTIDE"));
  hit2.setCharge(2);

  id4.insertHit(hit1);
  id5.insertHit(hit2);

  TEST_EQUAL(id4 == id5, true)
  TEST_EQUAL(std::hash<PeptideIdentification>{}(id4), std::hash<PeptideIdentification>{}(id5))

  // Test 6: Different peptide hits produce different hashes
  PeptideIdentification id6;
  PeptideHit hit3;
  hit3.setScore(99.0);
  hit3.setSequence(AASequence::fromString("DIFFERENT"));
  hit3.setCharge(3);
  id6.insertHit(hit3);

  TEST_EQUAL(id4 == id6, false)
  TEST_NOT_EQUAL(std::hash<PeptideIdentification>{}(id4), std::hash<PeptideIdentification>{}(id6))

  // Test 7: NaN handling - both NaN should be equal and have same hash
  PeptideIdentification id7, id8;
  // RT and MZ are NaN by default
  TEST_EQUAL(id7.hasRT(), false)
  TEST_EQUAL(id7.hasMZ(), false)
  TEST_EQUAL(id7 == id8, true)
  TEST_EQUAL(std::hash<PeptideIdentification>{}(id7), std::hash<PeptideIdentification>{}(id8))
}
END_SECTION

START_SECTION((USI buildUSI(const String& ms_run_name, const String& dataset_id, bool include_interpretation) const))
{
  // Test USI generation with scan number in native ID
  PeptideIdentification id1;
  id1.setSpectrumReference("controllerType=0 controllerNumber=1 scan=12345");
  id1.setRT(100.0);
  id1.setMZ(500.0);

  USI usi1 = id1.buildUSI("sample.mzML", "PXD000561", false);
  TEST_EQUAL(usi1.isValid(), true);
  TEST_STRING_EQUAL(usi1.getCollection(), "PXD000561");
  TEST_STRING_EQUAL(usi1.getMSRun(), "sample.mzML");
  TEST_STRING_EQUAL(usi1.getIndex(), "12345");
  TEST_EQUAL(usi1.getIndexType(), USI::IndexType::SCAN);
  TEST_STRING_EQUAL(usi1.getInterpretation(), "");

  // Test default dataset_id ("local")
  USI usi_local = id1.buildUSI("sample.mzML");
  TEST_EQUAL(usi_local.isValid(), true);
  TEST_STRING_EQUAL(usi_local.getCollection(), "local");
  TEST_STRING_EQUAL(usi_local.getMSRun(), "sample.mzML");

  // Test USI generation with peptide interpretation
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(2);
  hit.setScore(100.0);
  id1.insertHit(hit);
  id1.sort();

  USI usi2 = id1.buildUSI("sample.mzML", "PXD000561", true);
  TEST_EQUAL(usi2.isValid(), true);
  TEST_STRING_EQUAL(usi2.getInterpretation(), "PEPTIDEK/2");

  // Test USI generation with modified peptide interpretation (ProForma proteoform ion)
  PeptideIdentification id_mod;
  id_mod.setSpectrumReference("scan=12345");

  PeptideHit mod_hit;
  mod_hit.setSequence(AASequence::fromString("EM(Oxidation)K"));
  mod_hit.setCharge(2);
  mod_hit.setScore(100.0);
  id_mod.insertHit(mod_hit);
  id_mod.sort();

  USI usi_mod = id_mod.buildUSI("sample.mzML", "PXD000561", true);
  TEST_EQUAL(usi_mod.isValid(), true);
  TEST_STRING_EQUAL(usi_mod.getInterpretation(), "EM[UNIMOD:35]K/2");
  TEST_STRING_EQUAL(usi_mod.toString(), "mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2");

  // Test USI generation with simple scan= format
  PeptideIdentification id2;
  id2.setSpectrumReference("scan=54321");

  USI usi3 = id2.buildUSI("another.mzML", "PXD000562", false);
  TEST_EQUAL(usi3.isValid(), true);
  TEST_STRING_EQUAL(usi3.getIndex(), "54321");

  // Test USI with nativeId fallback (when scan can't be extracted)
  PeptideIdentification id3;
  id3.setSpectrumReference("custom_format_123");

  USI usi4 = id3.buildUSI("data.mzML", "PXD000563", false);
  TEST_EQUAL(usi4.isValid(), true);
  TEST_EQUAL(usi4.getIndexType(), USI::IndexType::NATIVEID);
  TEST_STRING_EQUAL(usi4.getIndex(), "custom_format_123");

  // Test empty spectrum reference returns invalid USI
  PeptideIdentification id4;
  USI usi5 = id4.buildUSI("test.mzML", "PXD000564", false);
  TEST_EQUAL(usi5.isValid(), false);
}
END_SECTION

START_SECTION(([EXTRA] buildUSI().toString() convenience pattern))
{
  PeptideIdentification id;
  id.setSpectrumReference("scan=12345");
  id.setRT(100.0);
  id.setMZ(500.0);

  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(2);
  id.insertHit(hit);

  // Test without interpretation
  String usi_str1 = id.buildUSI("sample.mzML", "PXD000561", false).toString();
  TEST_STRING_EQUAL(usi_str1, "mzspec:PXD000561:sample.mzML:scan:12345");

  // Test with interpretation
  String usi_str2 = id.buildUSI("sample.mzML", "PXD000561", true).toString();
  TEST_STRING_EQUAL(usi_str2, "mzspec:PXD000561:sample.mzML:scan:12345:PEPTIDEK/2");

  // Test default dataset_id ("local")
  String usi_str_local = id.buildUSI("sample.mzML").toString();
  TEST_STRING_EQUAL(usi_str_local, "mzspec:local:sample.mzML:scan:12345");

  // Test empty spectrum reference returns empty string
  PeptideIdentification id2;
  String usi_str3 = id2.buildUSI("sample.mzML", "PXD000561", false).toString();
  TEST_STRING_EQUAL(usi_str3, "");
}
END_SECTION

// Additional tests for file path handling in buildUSI
START_SECTION(([EXTRA] buildUSI with file paths and basenames))
{
  // Test buildUSI with full file path as MS run name
  PeptideIdentification id1;
  id1.setSpectrumReference("scan=12345");

  // Using full path as MS run (not recommended but should work)
  String full_path = "/path/to/data/sample.mzML";
  USI usi1 = id1.buildUSI(full_path, "PXD000561", false);
  TEST_EQUAL(usi1.isValid(), true);
  TEST_STRING_EQUAL(usi1.getMSRun(), full_path);

  // Test with basename extracted using USI::extractBasename
  String extracted_basename = USI::extractBasename(full_path);
  USI usi2 = id1.buildUSI(extracted_basename, "PXD000561", false);
  TEST_STRING_EQUAL(usi2.getMSRun(), "sample.mzML");
  TEST_STRING_EQUAL(usi2.toString(), "mzspec:PXD000561:sample.mzML:scan:12345");

  // Test with file:// URI (common in consensusXML)
  String file_uri = "file:///C:/Users/bielow/data/ES-0014b_2.mzML";
  USI usi3 = id1.buildUSI(USI::extractBasename(file_uri), "PXD000561", false);
  TEST_STRING_EQUAL(usi3.getMSRun(), "ES-0014b_2.mzML");

  // Test buildUSI().toString() with extracted basename
  String usi_str = id1.buildUSI(USI::extractBasename(file_uri), "PXD000561", false).toString();
  TEST_STRING_EQUAL(usi_str, "mzspec:PXD000561:ES-0014b_2.mzML:scan:12345");

  // Test with various native ID formats and basenames
  PeptideIdentification id2;
  id2.setSpectrumReference("spectrum=219");  // Another common format
  USI usi4 = id2.buildUSI("sample.mzML", "PXD000561", false);
  TEST_STRING_EQUAL(usi4.getIndex(), "219");
  TEST_EQUAL(usi4.getIndexType(), USI::IndexType::SCAN);

  // Test with Thermo native ID format
  PeptideIdentification id3;
  id3.setSpectrumReference("controllerType=0 controllerNumber=1 scan=7890");
  USI usi5 = id3.buildUSI("thermo_sample.raw", "PXD000561", false);
  TEST_STRING_EQUAL(usi5.getIndex(), "7890");
  TEST_STRING_EQUAL(usi5.getMSRun(), "thermo_sample.raw");
}
END_SECTION

// Tests for mapping-based buildUSI overloads (merged file support)
START_SECTION(([EXTRA] buildUSI with IdentifierMSRunMapper for merged files))
{
  // Create ProteinIdentifications to build a Mapping
  ProteinIdentification prot_id1;
  prot_id1.setIdentifier("run1");
  prot_id1.setPrimaryMSRunPath(StringList{"/data/file1.mzML"});

  ProteinIdentification prot_id2;
  prot_id2.setIdentifier("run2");
  prot_id2.setPrimaryMSRunPath(StringList{"/data/fileA.mzML", "/data/fileB.mzML"});  // merged

  std::vector<ProteinIdentification> prot_ids = {prot_id1, prot_id2};
  IdentifierMSRunMapper mapping(prot_ids);

  // Test single-file case (no merge index needed)
  PeptideIdentification id1;
  id1.setIdentifier("run1");
  id1.setSpectrumReference("scan=100");
  USI usi1 = id1.buildUSI(mapping, "PXD000561", false);
  TEST_EQUAL(usi1.isValid(), true);
  TEST_STRING_EQUAL(usi1.getMSRun(), "file1.mzML");
  TEST_STRING_EQUAL(usi1.getCollection(), "PXD000561");
  TEST_STRING_EQUAL(usi1.getIndex(), "100");

  // Test merged file case - default merge_index (0)
  PeptideIdentification id2;
  id2.setIdentifier("run2");
  id2.setSpectrumReference("scan=200");
  USI usi2 = id2.buildUSI(mapping, "PXD000561", false);
  TEST_EQUAL(usi2.isValid(), true);
  TEST_STRING_EQUAL(usi2.getMSRun(), "fileA.mzML");  // first file

  // Test merged file case - explicit merge_index = 1
  PeptideIdentification id3;
  id3.setIdentifier("run2");
  id3.setSpectrumReference("scan=300");
  id3.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, 1);
  USI usi3 = id3.buildUSI(mapping, "PXD000561", false);
  TEST_EQUAL(usi3.isValid(), true);
  TEST_STRING_EQUAL(usi3.getMSRun(), "fileB.mzML");  // second file

  // Test with default dataset_id ("local")
  USI usi4 = id1.buildUSI(mapping);
  TEST_EQUAL(usi4.isValid(), true);
  TEST_STRING_EQUAL(usi4.getCollection(), "local");

  // Test buildUSI().toString() with mapping
  String usi_str = id3.buildUSI(mapping, "PXD000561", false).toString();
  TEST_STRING_EQUAL(usi_str, "mzspec:PXD000561:fileB.mzML:scan:300");

  // Test invalid identifier returns invalid USI
  PeptideIdentification id_invalid;
  id_invalid.setIdentifier("unknown_run");
  id_invalid.setSpectrumReference("scan=999");
  USI usi_invalid = id_invalid.buildUSI(mapping, "PXD000561", false);
  TEST_EQUAL(usi_invalid.isValid(), false);

  // Test merge_index out of range returns invalid USI
  PeptideIdentification id_bad_index;
  id_bad_index.setIdentifier("run2");
  id_bad_index.setSpectrumReference("scan=999");
  id_bad_index.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, 99);  // out of range
  USI usi_bad = id_bad_index.buildUSI(mapping, "PXD000561", false);
  TEST_EQUAL(usi_bad.isValid(), false);

  // Test with interpretation
  PeptideIdentification id_with_hit;
  id_with_hit.setIdentifier("run1");
  id_with_hit.setSpectrumReference("scan=500");
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(2);
  id_with_hit.insertHit(hit);
  USI usi_interp = id_with_hit.buildUSI(mapping, "PXD000561", true);
  TEST_EQUAL(usi_interp.isValid(), true);
  TEST_STRING_EQUAL(usi_interp.getInterpretation(), "PEPTIDEK/2");
  TEST_STRING_EQUAL(usi_interp.toString(), "mzspec:PXD000561:file1.mzML:scan:500:PEPTIDEK/2");
}
END_SECTION

END_TEST

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
