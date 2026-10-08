// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/METADATA/ProteinRunReferences.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

///////////////////////////

START_TEST(ProteinRunReferences, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

auto run = [](const std::string& identifier)
{
  ProteinIdentification r;
  r.setIdentifier(identifier);
  return r;
};

auto peptide = [](const std::string& identifier)
{
  PeptideIdentification p;
  p.setIdentifier(identifier);
  return p;
};

const std::string missing_message =
  "Peptide identification has no matching protein run: 'other'. Every peptide identification needs the protein "
  "identification run (search run) with its identifier, which may have no protein hits.";

START_SECTION(static std::string missingRunMessage(const std::string& identifier))
{
  TEST_STRING_EQUAL(ProteinRunReferences::missingRunMessage("other"), missing_message)
}
END_SECTION

START_SECTION(static void check(const std::vector<ProteinIdentification>& runs, const std::string& identifier))
{
  const std::vector<ProteinIdentification> runs = {run("search"), run("")};
  ProteinRunReferences::check(runs, "search");
  ProteinRunReferences::check(runs, ""); // an empty identifier is an identifier like any other
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ProteinRunReferences::check(runs, "other"), missing_message)
  TEST_EXCEPTION(Exception::InvalidParameter, ProteinRunReferences::check(std::vector<ProteinIdentification>(), ""))
}
END_SECTION

START_SECTION(static void check(const std::vector<ProteinIdentification>& runs, const PeptideIdentificationList& peptides))
{
  // a run without protein hits (e.g. peptidomics or de novo) is a run
  const std::vector<ProteinIdentification> runs = {run("search")};
  PeptideIdentificationList peptides;
  ProteinRunReferences::check(std::vector<ProteinIdentification>(), peptides); // nothing to check
  peptides.push_back(peptide("search"));
  ProteinRunReferences::check(runs, peptides);
  peptides.push_back(peptide("other"));
  peptides.push_back(peptide("third"));
  // the first peptide identification without its run is reported
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ProteinRunReferences::check(runs, peptides), missing_message)
  TEST_EXCEPTION(Exception::InvalidParameter, ProteinRunReferences::check(std::vector<ProteinIdentification>(), peptides))
}
END_SECTION

START_SECTION(static void check(const std::vector<ProteinIdentification>& runs, const std::vector<const PeptideIdentification*>& peptides))
{
  const std::vector<ProteinIdentification> runs = {run("search")};
  const PeptideIdentification good = peptide("search"), bad = peptide("other");
  ProteinRunReferences::check(runs, std::vector<const PeptideIdentification*>{&good});
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ProteinRunReferences::check(runs, std::vector<const PeptideIdentification*>{&good, &bad}), missing_message)
}
END_SECTION

START_SECTION(static void check(const FeatureMap& map))
{
  FeatureMap map;
  map.setProteinIdentifications({run("search")});
  Feature feature;
  feature.getPeptideIdentifications().push_back(peptide("search"));
  map.push_back(feature);
  map.getUnassignedPeptideIdentifications().push_back(peptide("search"));
  ProteinRunReferences::check(map);

  // unassigned
  FeatureMap unassigned = map;
  unassigned.getUnassignedPeptideIdentifications().push_back(peptide("other"));
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ProteinRunReferences::check(unassigned), missing_message)

  // subordinate of a feature
  FeatureMap subordinate = map;
  Feature sub;
  sub.getPeptideIdentifications().push_back(peptide("other"));
  subordinate[0].getSubordinates().push_back(sub);
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ProteinRunReferences::check(subordinate), missing_message)

  // no runs at all
  FeatureMap no_runs = map;
  no_runs.getProteinIdentifications().clear();
  TEST_EXCEPTION(Exception::InvalidParameter, ProteinRunReferences::check(no_runs))
}
END_SECTION

START_SECTION(static void check(const ConsensusMap& map))
{
  ConsensusMap map;
  map.setProteinIdentifications({run("search")});
  ConsensusFeature feature;
  feature.getPeptideIdentifications().push_back(peptide("search"));
  map.push_back(feature);
  ProteinRunReferences::check(map);

  ConsensusMap assigned = map;
  assigned[0].getPeptideIdentifications().push_back(peptide("other"));
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ProteinRunReferences::check(assigned), missing_message)

  ConsensusMap unassigned = map;
  unassigned.getUnassignedPeptideIdentifications().push_back(peptide("other"));
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ProteinRunReferences::check(unassigned), missing_message)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
