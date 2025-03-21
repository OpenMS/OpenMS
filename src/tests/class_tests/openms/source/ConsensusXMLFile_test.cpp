// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
///////////////////////////

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;
using namespace std;


DRange<1> makeRange(double a, double b)
{
  DPosition<1> pa(a), pb(b);
  return DRange<1>(pa, pb);
}

START_TEST(ConsensusXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusXMLFile * ptr = nullptr;
ConsensusXMLFile* nullPointer = nullptr;
START_SECTION((ConsensusXMLFile()))
ptr = new ConsensusXMLFile();
TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~ConsensusXMLFile()))
delete ptr;
END_SECTION

TOLERANCE_ABSOLUTE(0.01)

START_SECTION(const PeakFileOptions& getOptions() const)
ConsensusXMLFile file;
TEST_EQUAL(file.getOptions().hasMSLevels(), false)
END_SECTION

START_SECTION(PeakFileOptions& getOptions())
ConsensusXMLFile file;
file.getOptions().addMSLevel(1);
TEST_EQUAL(file.getOptions().hasMSLevels(), true);
END_SECTION

START_SECTION((void load(const String &filename, ConsensusMap & map)))
ConsensusMap map;
ConsensusXMLFile file;
file.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), map);

//test DocumentIdentifier addition
TEST_STRING_EQUAL(map.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"));
TEST_STRING_EQUAL(FileTypes::typeToName(map.getLoadedFileType()), "consensusXML");

//meta data
TEST_EQUAL(map.getIdentifier(), "lsid")
TEST_EQUAL(map.getExperimentType() == "label-free", true)
TEST_EQUAL(map.getMetaValue("name1") == DataValue("value1"), true)
TEST_EQUAL(map.getMetaValue("name2") == DataValue(2), true)
//file descriptions
TEST_EQUAL(map.getColumnHeaders()[0].filename == "data/MapAlignmentFeatureMap1.xml", true)
TEST_EQUAL(map.getColumnHeaders()[0].label, "label")
TEST_EQUAL(map.getColumnHeaders()[0].size, 144)
TEST_EQUAL(map.getColumnHeaders()[0].getMetaValue("name3") == DataValue("value3"), true)
TEST_EQUAL(map.getColumnHeaders()[0].getMetaValue("name4") == DataValue(4), true)
TEST_STRING_EQUAL(map.getColumnHeaders()[1].filename, "data/MapAlignmentFeatureMap2.xml")
TEST_EQUAL(map.getColumnHeaders()[1].label, "")
TEST_EQUAL(map.getColumnHeaders()[1].size, 0)
TEST_EQUAL(map.getColumnHeaders()[1].getMetaValue("name5") == DataValue("value5"), true)
TEST_EQUAL(map.getColumnHeaders()[1].getMetaValue("name6") == DataValue(6.0), true)
//data processing
TEST_EQUAL(map.getDataProcessing().size(), 2)
TEST_STRING_EQUAL(map.getDataProcessing()[0].getSoftware().getName(), "Software1")
TEST_STRING_EQUAL(map.getDataProcessing()[0].getSoftware().getVersion(), "0.91a")
TEST_EQUAL(map.getDataProcessing()[0].getProcessingActions().size(), 1)
TEST_EQUAL(map.getDataProcessing()[0].getProcessingActions().count(DataProcessing::DEISOTOPING), 1)
TEST_STRING_EQUAL(map.getDataProcessing()[0].getMetaValue("name"), "dataProcessing")
TEST_STRING_EQUAL(map.getDataProcessing()[1].getSoftware().getName(), "Software2")
TEST_STRING_EQUAL(map.getDataProcessing()[1].getSoftware().getVersion(), "0.92a")
TEST_EQUAL(map.getDataProcessing()[1].getProcessingActions().size(), 2)
TEST_EQUAL(map.getDataProcessing()[1].getProcessingActions().count(DataProcessing::SMOOTHING), 1)
TEST_EQUAL(map.getDataProcessing()[1].getProcessingActions().count(DataProcessing::BASELINE_REDUCTION), 1)
//protein identifications
TEST_EQUAL(map.getProteinIdentifications().size(), 2)
TEST_EQUAL(map.getProteinIdentifications()[0].getHits().size(), 2)
TEST_EQUAL(map.getProteinIdentifications()[0].getHits()[0].getSequence(), "ABCDEFG")
TEST_EQUAL(map.getProteinIdentifications()[0].getHits()[1].getSequence(), "HIJKLMN")
TEST_EQUAL(map.getProteinIdentifications()[1].getHits().size(), 1)
TEST_EQUAL(map.getProteinIdentifications()[1].getHits()[0].getSequence(), "OPQREST")
//peptide identifications
TEST_EQUAL(map[0].getPeptideIdentifications().size(), 2)
TEST_EQUAL(map[0].getPeptideIdentifications()[0].getHits().size(), 1)
TEST_EQUAL(map[0].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("A"))
TEST_EQUAL(map[0].getPeptideIdentifications()[1].getHits().size(), 2)
TEST_EQUAL(map[0].getPeptideIdentifications()[1].getHits()[0].getSequence(), AASequence::fromString("C"))
TEST_EQUAL(map[0].getPeptideIdentifications()[1].getHits()[1].getSequence(), AASequence::fromString("D"))
TEST_EQUAL(map[1].getPeptideIdentifications().size(), 1)
TEST_EQUAL(map[1].getPeptideIdentifications()[0].getHits().size(), 1)
TEST_EQUAL(map[1].getPeptideIdentifications()[0].getHits()[0].getSequence(),AASequence::fromString( "E"))
//unassigned peptide identifications
TEST_EQUAL(map.getUnassignedPeptideIdentifications().size(), 2)
TEST_EQUAL(map.getUnassignedPeptideIdentifications()[0].getHits().size(), 1)
TEST_EQUAL(map.getUnassignedPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("F"))
TEST_EQUAL(map.getUnassignedPeptideIdentifications()[1].getHits().size(), 2)
TEST_EQUAL(map.getUnassignedPeptideIdentifications()[1].getHits()[0].getSequence(), AASequence::fromString("G"))
TEST_EQUAL(map.getUnassignedPeptideIdentifications()[1].getHits()[1].getSequence(), AASequence::fromString("H"))

//features
TEST_EQUAL(map.size(), 6)
ConsensusFeature cons_feature = map[0];
TEST_REAL_SIMILAR(cons_feature.getRT(), 1273.27)
TEST_REAL_SIMILAR(cons_feature.getMZ(), 904.47)
TEST_REAL_SIMILAR(cons_feature.getIntensity(), 3.12539e+07)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().minPosition()[0], 1273.27)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().maxPosition()[0], 1273.27)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().minPosition()[1], 904.47)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().maxPosition()[1], 904.47)
TEST_REAL_SIMILAR(cons_feature.getIntensityRange().minPosition()[0], 3.12539e+07)
TEST_REAL_SIMILAR(cons_feature.getIntensityRange().maxPosition()[0], 3.12539e+07)
TEST_REAL_SIMILAR(cons_feature.getQuality(), 1.1)
TEST_EQUAL(cons_feature.getMetaValue("peptide_id") == DataValue("RefSeq:NC_1234"), true)
ConsensusFeature::HandleSetType::const_iterator it = cons_feature.begin();
TEST_REAL_SIMILAR(it->getIntensity(), 3.12539e+07)

cons_feature = map[5];
TEST_REAL_SIMILAR(cons_feature.getRT(), 1194.82)
TEST_REAL_SIMILAR(cons_feature.getMZ(), 777.101)
TEST_REAL_SIMILAR(cons_feature.getIntensity(), 1.78215e+07)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().minPosition()[0], 1194.82)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().maxPosition()[0], 1194.82)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().minPosition()[1], 777.101)
TEST_REAL_SIMILAR(cons_feature.getPositionRange().maxPosition()[1], 777.101)
TEST_REAL_SIMILAR(cons_feature.getIntensityRange().minPosition()[0], 1.78215e+07)
TEST_REAL_SIMILAR(cons_feature.getIntensityRange().maxPosition()[0], 1.78215e+07)
TEST_REAL_SIMILAR(cons_feature.getQuality(), 0.0)
it = cons_feature.begin();
TEST_REAL_SIMILAR(it->getIntensity(), 1.78215e+07)
++ it;
TEST_REAL_SIMILAR(it->getIntensity(), 1.78215e+07)

// test meta values:
TEST_EQUAL(map[0].getMetaValue("myIntList") == ListUtils::create<Int>("1,10,12"), true);
TEST_EQUAL(map[0].getMetaValue("myDoubleList") == ListUtils::create<double>("1.111,10.999,12.45"), true);
std::cout << "list: " << map[0].getMetaValue("myStringList") << "\n";
TEST_EQUAL(map[0].getMetaValue("myStringList") == ListUtils::create<String>("myABC1,Stuff,12"), true);
TEST_EQUAL(map[4].getMetaValue("myDoubleList") == ListUtils::create<double>("6.442"), true);

//PeakFileOptions tests

file.getOptions().setRTRange(makeRange(815, 818));
file.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_2_options.consensusXML"), map);
TEST_EQUAL(map.size(), 1)
TEST_REAL_SIMILAR(map[0].getRT(), 817.266)

file.getOptions() = PeakFileOptions();
file.getOptions().setMZRange(makeRange(944, 945));
file.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_2_options.consensusXML"), map);
TEST_EQUAL(map.size(), 1)
TEST_REAL_SIMILAR(map[0].getMZ(), 944.96)

file.getOptions() = PeakFileOptions();
file.getOptions().setIntensityRange(makeRange(15000, 24000));
file.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_2_options.consensusXML"), map);
TEST_EQUAL(map.size(), 1)
TEST_REAL_SIMILAR(map[0].getIntensity(), 23000.238)

END_SECTION

START_SECTION((void store(const String &filename, const ConsensusMap &consensus_map)))
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  ConsensusMap map;
  ConsensusXMLFile f;

  f.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), map);

  // make protIDs non-unique
  map.getProteinIdentifications().push_back(map.getProteinIdentifications()[0]);
  TEST_EXCEPTION(Exception::ParseError, f.store(tmp_filename, map))
  map.getProteinIdentifications().pop_back(); // undo

  f.store(tmp_filename, map);
  WHITELIST("?xml-stylesheet")
  TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), tmp_filename)


END_SECTION

START_SECTION([EXTRA](bool isValid(const String &filename)))
  ConsensusXMLFile f;
  TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), std::cerr), true);
  TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_2_options.consensusXML"), std::cerr), true);

  //test if written empty file
  // - this is invalid, so it is not tested :)

  //test if written full file is valid
  ConsensusMap m;
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  f.load(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile_1.consensusXML"), m);
  f.store(tmp_filename, m);
  TEST_EQUAL(f.isValid(tmp_filename, std::cerr), true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
