// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Hendrik Beschorner, Lenny Kovac, Virginia Rossow$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/QC/MQMsmsExporter.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/test_config.h>

///////////////////////////
///////////////////////////

START_TEST(MQMsms, "$ID$")

using namespace OpenMS;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


File::TempDir dir;
const String path = dir.getPath();



START_SECTION(MQMsms())
{
  MQMsms *ptr = nullptr;
  MQMsms *null_ptr = nullptr;
  ptr = new MQMsms(path);
  TEST_NOT_EQUAL(ptr, null_ptr);
  delete ptr;
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
START_SECTION((void exportFeatureMap(
                const OpenMS::FeatureMap& feature_map,
                const OpenMS::ConsensusMap& cmap,
                const OpenMS::MSExperiment& exp,
                const std::map<OpenMS::String,OpenMS::String>& prot_map = {})))

{
  {
    MQMsms msms(path);
    PeakMap exp;
    ConsensusMap cmap_one;
    ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_1.consensusXML"), cmap_one);
    ConsensusMap cmap_two;
    ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_2.consensusXML"), cmap_two);
    FeatureMap fmap_one;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_1.featureXML"), fmap_one);
    msms.exportFeatureMap(fmap_one, cmap_two, exp);
    FeatureMap fmap_two;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_2.featureXML"), fmap_two);
    msms.exportFeatureMap(fmap_two, cmap_two, exp);
    FeatureMap fmap_three;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_3.featureXML"), fmap_three);
    msms.exportFeatureMap(fmap_three, cmap_two, exp);
    FeatureMap fmap_four;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_4.featureXML"), fmap_four);
    msms.exportFeatureMap(fmap_four, cmap_one, exp);
    FeatureMap fmap_five;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_5.featureXML"), fmap_five);
    msms.exportFeatureMap(fmap_five, cmap_one, exp);
    FeatureMap fmap_six;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_6.featureXML"), fmap_six);
    msms.exportFeatureMap(fmap_six, cmap_one, exp);
  }
  String filename = path + "/msms.txt";
  TEST_FILE_SIMILAR(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("MQMsms_result.txt"));
}

END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST