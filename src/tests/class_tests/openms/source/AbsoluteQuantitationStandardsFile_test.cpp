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
#include <OpenMS/FORMAT/AbsoluteQuantitationStandardsFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(AbsoluteQuantitationStandardsFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

AbsoluteQuantitationStandardsFile* ptr = 0;
AbsoluteQuantitationStandardsFile* null_ptr = 0;

const String csv_path = OPENMS_GET_TEST_DATA_PATH("AbsoluteQuantitationStandardsFile.csv");

START_SECTION(AbsoluteQuantitationStandardsFile())
{
  ptr = new AbsoluteQuantitationStandardsFile();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~AbsoluteQuantitationStandardsFile())
{
  delete ptr;
}
END_SECTION

ptr = new AbsoluteQuantitationStandardsFile();

START_SECTION(void load(
  const String& filename,
  std::vector<AbsoluteQuantitationStandards::runConcentration>& run_concentrations
) const)
{
  std::vector<AbsoluteQuantitationStandards::runConcentration> runs;
  ptr->load(csv_path, runs);
  TEST_EQUAL(runs.size(), 10)
  TEST_EQUAL(runs[0].sample_name, "150516_CM1_Level1")
  TEST_EQUAL(runs[0].component_name, "23dpg.23dpg_1.Light")
  TEST_EQUAL(runs[0].IS_component_name, "23dpg.23dpg_1.Heavy")
  TEST_REAL_SIMILAR(runs[0].actual_concentration, 0)
  TEST_REAL_SIMILAR(runs[0].IS_actual_concentration, 1)
  TEST_EQUAL(runs[0].concentration_units, "uM")
  TEST_REAL_SIMILAR(runs[0].dilution_factor, 1)
  TEST_EQUAL(runs[3].sample_name, "150516_CM1_Level1")
  TEST_EQUAL(runs[3].component_name, "35cgmp.35cgmp_1.Light")
  TEST_EQUAL(runs[3].IS_component_name, "35cgmp.35cgmp_1.Heavy")
  TEST_REAL_SIMILAR(runs[3].actual_concentration, 0)
  TEST_REAL_SIMILAR(runs[3].IS_actual_concentration, 1)
  TEST_EQUAL(runs[3].concentration_units, "uM")
  TEST_REAL_SIMILAR(runs[3].dilution_factor, 1)
  TEST_EQUAL(runs[6].sample_name, "150516_CM1_Level1")
  TEST_EQUAL(runs[6].component_name, "accoa.accoa_1.Light")
  TEST_EQUAL(runs[6].IS_component_name, "accoa.accoa_1.Heavy")
  TEST_REAL_SIMILAR(runs[6].actual_concentration, 0)
  TEST_REAL_SIMILAR(runs[6].IS_actual_concentration, 1)
  TEST_EQUAL(runs[6].concentration_units, "uM")
  TEST_REAL_SIMILAR(runs[6].dilution_factor, 1)
  TEST_EQUAL(runs[9].sample_name, "150516_CM1_Level1")
  TEST_EQUAL(runs[9].component_name, "ade.ade_1.Light")
  TEST_EQUAL(runs[9].IS_component_name, "ade.ade_1.Heavy")
  TEST_REAL_SIMILAR(runs[9].actual_concentration, 40)
  TEST_REAL_SIMILAR(runs[9].IS_actual_concentration, 1)
  TEST_EQUAL(runs[9].concentration_units, "uM")
  TEST_REAL_SIMILAR(runs[9].dilution_factor, 1)
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
