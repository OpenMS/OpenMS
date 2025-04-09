// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/AbsoluteQuantitationStandards.h>
#include <OpenMS/FORMAT/AbsoluteQuantitationStandardsFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(AbsoluteQuantitationStandards, "$Id$")

/////////////////////////////////////////////////////////////

AbsoluteQuantitationStandards* ptr = 0;
AbsoluteQuantitationStandards* null_ptr = 0;

vector<AbsoluteQuantitationStandards::runConcentration> runs;
AbsoluteQuantitationStandards::runConcentration run;
for (Size i = 0; i < 10; ++i)
{
  run.sample_name = i < 5 ? String("sample1") : String("sample2");
  run.component_name = String("component") + i;
  run.IS_component_name = String("IS_component") + i;
  run.actual_concentration = i;
  run.IS_actual_concentration = i * 1.1;
  run.concentration_units = "uM";
  run.dilution_factor = 1;
  runs.push_back(run);
}
run.sample_name = "";
runs.push_back(run); // without sample_name
run.sample_name = "sample2";
run.component_name = "";
runs.push_back(run); // without component_name
run.component_name = "component10";
run.IS_component_name = "";
runs.push_back(run); // without IS_component_name
run.component_name = "component11";
runs.push_back(run); // without IS_component_name and no match for component_name
run.component_name = "component0";
runs.push_back(run); // with a component_name equal to one of those in sample1

vector<FeatureMap> fmaps;
FeatureMap fm;
Feature feature;
vector<Feature> subordinates;

fm.setPrimaryMSRunPath({"sample1.mzML"});
for (Size i = 0; i < 5; ++i)
{
  Feature f;
  f.setMetaValue("native_id", String("component") + i);
  subordinates.push_back(f);
  f.setMetaValue("native_id", String("IS_component") + i);
  subordinates.push_back(f);
}
feature.setSubordinates(subordinates);
fm.push_back(feature);
fmaps.push_back(fm);
// The first FeatureMap has sample_name: "sample1".  It contains 1 feature. This feature has 10 subordinates:
// 5 subordinates have native_id: "component0" to "component4", and the other 5 subordinates have native_id: "IS_component0" to "IS_component4".

fm.setPrimaryMSRunPath({"sample2.txt"});
Feature f;
f.setMetaValue("native_id", String("component10"));
subordinates.push_back(f);
f.setMetaValue("native_id", String("component0"));
subordinates.push_back(f);
feature.setSubordinates(subordinates);
fm.push_back(feature);
fmaps.push_back(fm);
// The second FeatureMap has sample_name: "sample2". It contains 1 feature. This feature has 2 subordinates. Their native_id are "component10" and "component0".

START_SECTION(AbsoluteQuantitationStandards())
{
  ptr = new AbsoluteQuantitationStandards();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~AbsoluteQuantitationStandards())
{
  delete ptr;
}
END_SECTION

START_SECTION(void mapComponentsToConcentrations(
  const std::vector<runConcentration>& run_concentrations,
  const std::vector<FeatureMap>& feature_maps,
  std::map<String, std::vector<featureConcentration>>& components_to_concentrations
) const)
{
  AbsoluteQuantitationStandards aqs;
  std::map<String, std::vector<AbsoluteQuantitationStandards::featureConcentration>> m;
  aqs.mapComponentsToConcentrations(runs, fmaps, m);
  TEST_EQUAL(m.size(), 6)
  std::vector<AbsoluteQuantitationStandards::featureConcentration> fc;
  for (Size i = 0; i < 5; ++i)
  {
    fc = m.at(String("component") + i);
    TEST_EQUAL(fc[0].feature.getMetaValue("native_id"), String("component") + i)
    TEST_EQUAL(fc[0].IS_feature.getMetaValue("native_id"), String("IS_component") + i)
  }
  fc = m.at("component10");
  TEST_EQUAL(fc.size(), 1)
  TEST_EQUAL(fc[0].feature.getMetaValue("native_id"), "component10")
  TEST_EQUAL(fc[0].IS_feature.metaValueExists("native_id"), false)
  fc = m.at("component0");
  TEST_EQUAL(fc.size(), 2)
  TEST_EQUAL(fc[1].feature.getMetaValue("native_id"), "component0")
  TEST_EQUAL(fc[1].IS_feature.metaValueExists("native_id"), false)
}
END_SECTION

START_SECTION(void getComponentFeatureConcentrations(
  const std::vector<AbsoluteQuantitationStandards::runConcentration>& run_concentrations,
  const std::vector<FeatureMap>& feature_maps,
  const String& component_name,
  std::vector<AbsoluteQuantitationStandards::featureConcentration>& feature_concentrations
) const)
{
  AbsoluteQuantitationStandards aqs;
  std::vector<AbsoluteQuantitationStandards::featureConcentration> fc;
  aqs.getComponentFeatureConcentrations(runs, fmaps, "component0", fc);
  TEST_EQUAL(fc.size(), 2)
  TEST_EQUAL(fc[0].feature.getMetaValue("native_id"), "component0")
  TEST_EQUAL(fc[0].IS_feature.getMetaValue("native_id"), "IS_component0")
  TEST_EQUAL(fc[1].feature.getMetaValue("native_id"), "component0")
  TEST_EQUAL(fc[1].IS_feature.metaValueExists("native_id"), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
