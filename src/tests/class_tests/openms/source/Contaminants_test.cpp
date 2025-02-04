// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Dominik Schmitz, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/QC/Contaminants.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Contaminants, "$Id$")
FeatureMap fmap;
FeatureMap emptyFmap;
Feature f;

fmap.getProteinIdentifications().resize(1);
DigestionEnzymeProtein noenzyme("unknown_enzyme", "", set<String>(), "");
// set no digestion enzyme
fmap.getProteinIdentifications()[0].getSearchParameters().digestion_enzyme = noenzyme;
// set empty contaminants database
vector<FASTAFile::FASTAEntry> contaminantsFile;

// fill the featureMap of features with set sequence and intensity
{
  PeptideIdentification id;
  PeptideIdentification scnd_id;

  PeptideHit hit;
  PeptideHit scnd_hit;

  hit.setSequence(AASequence::fromString("AAAAAAAAAAK"));
  id.setHits({hit});
  f.setPeptideIdentifications({id});
  f.setIntensity(12.0);
  fmap.push_back(f);

  hit.setSequence(AASequence::fromString("R"));
  id.setHits({hit});
  f.setPeptideIdentifications({id});
  f.setIntensity(8.0);
  fmap.push_back(f);

  hit.setSequence(AASequence::fromString("R"));
  scnd_hit.setSequence(AASequence::fromString("QQQQQQQQQQ"));
  id.setHits({hit});
  scnd_id.setHits({scnd_hit});
  f.setPeptideIdentifications({id, scnd_id});
  f.setIntensity(10.0);
  fmap.push_back(f);

  hit.setSequence(AASequence::fromString("AAAAAAAAAAKR"));
  id.setHits({hit});
  f.setPeptideIdentifications({id});
  f.setIntensity(20.0);
  fmap.push_back(f);

  hit.setSequence(AASequence::fromString("AAAAAAAAAAKRAAAAAAAAAAKRCCCCCCCCCCKRCCCCCCCCCC"));
  id.setHits({hit});
  f.setPeptideIdentifications({id});
  f.setIntensity(10.0);
  fmap.push_back(f);

  f.setPeptideIdentifications({});
  fmap.push_back(f);
}

// set the unassigned peptideidentifications
std::vector<PeptideIdentification> ids2(3);
PeptideHit hit2;
hit2.setSequence(AASequence::fromString("AAAAAAAAAAK"));
ids2[0].setHits({hit2});
hit2.setSequence(AASequence::fromString("RCCCCCCCCCCK"));
ids2[1].setHits({hit2});
hit2.setSequence(AASequence::fromString("DDDDDDDDDD"));
ids2[2].setHits({hit2});

// check the constructor
Contaminants* ptr = nullptr;
Contaminants* nullPointer = nullptr;
START_SECTION(Contaminants())
ptr = new Contaminants();
TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~Contaminants())
delete ptr;
END_SECTION

START_SECTION((void compute(FeatureMap& features, const std::vector<FASTAFile::FASTAEntry>& contaminants)))
{
  Contaminants conts1, conts2, conts3, conts4, conts5, conts6, conts7;

  // test exception when the contaminants database is empty
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, conts1.compute(fmap, contaminantsFile), "No contaminants provided.");

  // set contaminant database "contaminantsFile"
  FASTAFile::FASTAEntry contaminantsProtein("test_protein", "protein consists of only Alanine or Cytosine", "AAAAAAAAAAKRAAAAAAAAAAKRCCCCCCCCCCKRCCCCCCCCCC");
  contaminantsFile.push_back(contaminantsProtein);

  // test exception when the proteinidentification in the FeatureMap is empty
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, conts2.compute(emptyFmap, contaminantsFile), "No proteinidentifications in FeatureMap.");

  // test exception when no digestion enzyme is given
  TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, conts6.compute(fmap, contaminantsFile), "No digestion enzyme in FeatureMap detected. No computation possible.");

  // tests without given missed cleavages and without given enzyme
  fmap.getProteinIdentifications()[0].getSearchParameters().digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("no cleavage");
  conts3.compute(fmap, contaminantsFile);
  std::vector<Contaminants::ContaminantsSummary> result3 = conts3.getResults();
  ABORT_IF(result3.size() != 1);
  TEST_REAL_SIMILAR(result3[0].assigned_contaminants_ratio, 1 / 6.0);
  TEST_REAL_SIMILAR(result3[0].assigned_contaminants_intensity_ratio, 1 / 7.0);
  TEST_REAL_SIMILAR(result3[0].all_contaminants_ratio, 1 / 6.0);
  TEST_EQUAL(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap[1].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap[2].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap[2].getPeptideIdentifications()[1].getHits()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap[3].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap[4].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 1);

  // set digestion enzyme to trypsin
  fmap.getProteinIdentifications()[0].getSearchParameters().digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme("trypsin");

  conts7.compute(fmap, contaminantsFile);
  std::vector<Contaminants::ContaminantsSummary> result7 = conts7.getResults();
  TEST_REAL_SIMILAR(result7[0].assigned_contaminants_ratio, 3 / 6.0);
  TEST_REAL_SIMILAR(result7[0].assigned_contaminants_intensity_ratio, 3 / 7.0);
  TEST_REAL_SIMILAR(result7[0].all_contaminants_ratio, 3 / 6.0);
  TEST_EQUAL(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[1].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[2].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[2].getPeptideIdentifications()[1].getHits()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap[3].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap[4].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 0);

  // fill the unassigned peptideidentifications
  fmap.setUnassignedPeptideIdentifications(ids2);

  // tests without given missed cleavages but with set enzyme
  conts4.compute(fmap, contaminantsFile);
  std::vector<Contaminants::ContaminantsSummary> result4 = conts4.getResults();
  ABORT_IF(result4.size() != 1);
  TEST_REAL_SIMILAR(result4[0].assigned_contaminants_ratio, 3 / 6.0);
  TEST_REAL_SIMILAR(result4[0].assigned_contaminants_intensity_ratio, 3 / 7.0);
  TEST_REAL_SIMILAR(result4[0].unassigned_contaminants_ratio, 1 / 3.0);
  TEST_REAL_SIMILAR(result4[0].all_contaminants_ratio, 4 / 9.0);
  TEST_EQUAL(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[1].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[2].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[2].getPeptideIdentifications()[1].getHits()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap[3].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap[4].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[1].getHits()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[2].getHits()[0].getMetaValue("is_contaminant"), 0);

  // set missed cleavages to 1
  fmap.getProteinIdentifications()[0].getSearchParameters().missed_cleavages = 1;

  // tests with set missed cleavages and set enzyme
  // also checks if the empty feature count is as expected
  conts5.compute(fmap, contaminantsFile);
  std::vector<Contaminants::ContaminantsSummary> result5 = conts5.getResults();
  ABORT_IF(result5.size() != 1);
  TEST_REAL_SIMILAR(result5[0].assigned_contaminants_ratio, 4 / 6.0);
  TEST_REAL_SIMILAR(result5[0].assigned_contaminants_intensity_ratio, 5 / 7.0);
  TEST_REAL_SIMILAR(result5[0].unassigned_contaminants_ratio, 2 / 3.0);
  TEST_REAL_SIMILAR(result5[0].all_contaminants_ratio, 6 / 9.0);
  TEST_EQUAL(fmap[0].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[1].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[2].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[2].getPeptideIdentifications()[1].getHits()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap[3].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap[4].getPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[0].getHits()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[1].getHits()[0].getMetaValue("is_contaminant"), 1);
  TEST_EQUAL(fmap.getUnassignedPeptideIdentifications()[2].getHits()[0].getMetaValue("is_contaminant"), 0);
  TEST_EQUAL(result5[0].empty_features.first, 1);
  TEST_EQUAL(result5[0].empty_features.second, 6);
}
END_SECTION


Contaminants temp;

START_SECTION(const String& getName() const override) {TEST_EQUAL(temp.getName(), "Contaminants")} END_SECTION


  START_SECTION(Status requirements() const override)
{
  TEST_EQUAL(temp.requirements() == (QCBase::Status(QCBase::Requires::POSTFDRFEAT) | QCBase::Requires::CONTAMINANTS), true);
}
END_SECTION
END_TEST
