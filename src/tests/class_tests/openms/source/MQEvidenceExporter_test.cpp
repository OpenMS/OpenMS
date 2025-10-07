// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Valentin Noske, Vincent Musch$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/QC/MQEvidenceExporter.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <map>

#include <OpenMS/test_config.h>

///////////////////////////
///////////////////////////

START_TEST(MQEvidence, "$ID$")

using namespace OpenMS;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

File::TempDir dir;
const String path = dir.getPath();

START_SECTION(MQEvidence())
{
  MQEvidence *ptr = nullptr;
  MQEvidence *null_ptr = nullptr;
  ptr = new MQEvidence(path);
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
                    const std::map<String, String>& fasta_map)))

{
  {

    //FASTA-HANDELING
    std::vector<FASTAFile::FASTAEntry> fasta_info;
    std::map<String, String> fasta_map;
    FASTAFile().load(OPENMS_GET_TEST_DATA_PATH("FASTAContainer_test.fasta"), fasta_info);
    //map the identifier to the description so that we can access the description via the cmap-identifier
    for(const auto& entry : fasta_info)
    {
      fasta_map.emplace(entry.identifier, entry.description);
    }

    MQEvidence evd(path);
    PeakMap exp;
    ConsensusMap cmap_one;
    ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_1.consensusXML"), cmap_one);
    ConsensusMap cmap_two;
    ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_2.consensusXML"), cmap_two);
    FeatureMap fmap_one;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_1.featureXML"), fmap_one);
    evd.exportFeatureMap(fmap_one, cmap_two, exp, fasta_map);
    FeatureMap fmap_two;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_2.featureXML"), fmap_two);
    evd.exportFeatureMap(fmap_two, cmap_two, exp, fasta_map);
    FeatureMap fmap_three;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_3.featureXML"), fmap_three);
    evd.exportFeatureMap(fmap_three, cmap_two, exp, fasta_map);
    FeatureMap fmap_four;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_4.featureXML"), fmap_four);
    evd.exportFeatureMap(fmap_four, cmap_one, exp, fasta_map);
    FeatureMap fmap_five;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_5.featureXML"), fmap_five);
    evd.exportFeatureMap(fmap_five, cmap_one, exp, fasta_map);
    FeatureMap fmap_six;
    FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MQEvidence_6.featureXML"), fmap_six);
    evd.exportFeatureMap(fmap_six, cmap_one, exp, fasta_map);
  }
  String filename = path + "/evidence.txt";
  TEST_FILE_SIMILAR(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("MQEvidence_result.txt"));
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST