// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/IBSpectraFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(IBSpectraFile, "$Id$")

IBSpectraFile* ptr = nullptr;
IBSpectraFile* nullPointer = nullptr;

START_SECTION((IBSpectraFile()))
{
  ptr = new IBSpectraFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((IBSpectraFile(const IBSpectraFile& other)))
{
  IBSpectraFile ibfile(*ptr);
  TEST_NOT_EQUAL(&ibfile, nullPointer)
}
END_SECTION

START_SECTION((IBSpectraFile& operator=(const IBSpectraFile& rhs)))
{
  IBSpectraFile ibfile;
  ibfile = *ptr;
  TEST_NOT_EQUAL(&ibfile, nullPointer)
}
END_SECTION

START_SECTION((void store(const String& filename, const ConsensusMap& cm)))
{
  // test invalid ConsensusMap
  ConsensusMap cm_no_ms2quant;
  cm_no_ms2quant.setExperimentType("labeled_MS1");

  IBSpectraFile ibfile_no_ms2quant;
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ibfile_no_ms2quant.store("not-a-file-name", cm_no_ms2quant), "Given ConsensusMap does not hold any isobaric quantification data.")

  // test wrong channel count
  ConsensusMap cm_wrong_channel_count;
  cm_wrong_channel_count.setExperimentType("labeled_MS2");
  ConsensusMap::ColumnHeader channel1;
  ConsensusMap::ColumnHeader channel2;
  ConsensusMap::ColumnHeader channel3;
  cm_wrong_channel_count.getColumnHeaders()[0] = channel1;
  cm_wrong_channel_count.getColumnHeaders()[1] = channel2;
  cm_wrong_channel_count.getColumnHeaders()[2] = channel3;

  IBSpectraFile ibfile_wrong_channel_count;
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ibfile_wrong_channel_count.store("not-a-file-name", cm_wrong_channel_count), "Could not guess isobaric quantification data from ConsensusMap due to non-matching number of input maps.")

  // test a real example
  ConsensusMap cm;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IBSpectraFile.consensusXML"),cm);

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  IBSpectraFile ibfile;
  ibfile.store(tmp_filename, cm);

  TEST_FILE_SIMILAR(tmp_filename.c_str(), OPENMS_GET_TEST_DATA_PATH("IBSpectraFile.ibspectra.csv"))
}
END_SECTION

delete ptr;

END_TEST
