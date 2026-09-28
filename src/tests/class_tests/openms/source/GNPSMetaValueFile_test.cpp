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
#include <OpenMS/FORMAT/GNPSMetaValueFile.h>
///////////////////////////

#include <OpenMS/KERNEL/ConsensusMap.h>

#include <fstream>
#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(GNPSMetaValueFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static void store(const ConsensusMap& consensus_map, const std::string& output_file)))
{
  // a ConsensusMap with two input maps (the GNPS metadata file lists one row
  // per map: index, mzML basename, MAP<index>)
  ConsensusMap cm;
  ConsensusMap::ColumnHeader h0; h0.filename = "/data/runA.mzML"; h0.label = "lf"; h0.size = 10;
  ConsensusMap::ColumnHeader h1; h1.filename = "/data/runB.mzML"; h1.label = "lf"; h1.size = 12;
  cm.getColumnHeaders()[0] = h0;
  cm.getColumnHeaders()[1] = h1;

  std::string fn;
  NEW_TMP_FILE(fn)
  GNPSMetaValueFile::store(cm, fn);

  std::ifstream is(fn.c_str());
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(is, line)) lines.push_back(line);

  TEST_EQUAL(lines.size(), 3) // header + one row per map
  // header
  TEST_STRING_EQUAL(lines[0], "\tfilename\tATTRIBUTE_MAPID")
  // per-map rows: index, basename (path stripped), MAP<index>
  TEST_STRING_EQUAL(lines[1], "0\trunA.mzML\tMAP0")
  TEST_STRING_EQUAL(lines[2], "1\trunB.mzML\tMAP1")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
