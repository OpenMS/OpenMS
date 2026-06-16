// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/GNPSQuantificationFile.h>
///////////////////////////

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <fstream>
#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(GNPSQuantificationFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static void store(const ConsensusMap& consensus_map, const std::string& output_file)))
{
  // Two input maps and one consensus feature with a sub-feature in each map.
  ConsensusMap cm;
  ConsensusMap::ColumnHeader h0; h0.filename = "/data/runA.mzML"; h0.label = "lf"; h0.size = 10;
  ConsensusMap::ColumnHeader h1; h1.filename = "/data/runB.mzML"; h1.label = "lf"; h1.size = 12;
  cm.getColumnHeaders()[0] = h0;
  cm.getColumnHeaders()[1] = h1;

  ConsensusFeature cf;
  cf.setRT(100.0); cf.setMZ(500.0); cf.setIntensity(2000.0); cf.setCharge(2); cf.setUniqueId(42);
  Peak2D p0; p0.setRT(100.1); p0.setMZ(500.0); p0.setIntensity(1200.0);
  cf.insert(FeatureHandle(0, p0, 111));
  Peak2D p1; p1.setRT(99.9); p1.setMZ(500.0); p1.setIntensity(800.0);
  cf.insert(FeatureHandle(1, p1, 222));
  cm.push_back(cf);

  std::string fn;
  NEW_TMP_FILE(fn)
  GNPSQuantificationFile::store(cm, fn);

  std::ifstream is(fn.c_str());
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(is, line)) lines.push_back(line);

  // #MAP / #CONSENSUS column headers, one MAP row per input map, one CONSENSUS row
  TEST_EQUAL(lines.size(), 5)
  TEST_STRING_EQUAL(lines[0], "#MAP\tid\tfilename\tlabel\tsize")
  TEST_STRING_EQUAL(lines[1], "#CONSENSUS\trt_cf\tmz_cf\tintensity_cf\tcharge_cf\twidth_cf\tquality_cf\trt_0\tmz_0\tintensity_0\tcharge_0\twidth_0\trt_1\tmz_1\tintensity_1\tcharge_1\twidth_1")
  TEST_STRING_EQUAL(lines[2], "MAP\t0\t/data/runA.mzML\tlf\t10")
  TEST_STRING_EQUAL(lines[3], "MAP\t1\t/data/runB.mzML\tlf\t12")

  // consensus row: check the stable (non-RT) numeric fields by column index
  StringList c;
  StringUtils::split(lines[4], '\t', c);
  TEST_STRING_EQUAL(c[0], "CONSENSUS")
  TEST_STRING_EQUAL(c[1], "100.0")   // rt_cf (set exactly)
  TEST_STRING_EQUAL(c[2], "500.0")   // mz_cf
  TEST_STRING_EQUAL(c[3], "2000.0")  // intensity_cf
  TEST_STRING_EQUAL(c[4], "2")       // charge_cf
  TEST_STRING_EQUAL(c[9], "1200.0")  // intensity of the map-0 sub-feature
  TEST_STRING_EQUAL(c[14], "800.0")  // intensity of the map-1 sub-feature
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
