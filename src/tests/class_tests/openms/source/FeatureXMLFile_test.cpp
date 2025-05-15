// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
///////////////////////////

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

DRange<1> makeRange(double a, double b)
{
  DPosition<1> pa(a), pb(b);
  return DRange<1>(pa, pb);
}

///////////////////////////

START_TEST(FeatureXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureXMLFile * ptr = nullptr;
FeatureXMLFile* nullPointer = nullptr;
START_SECTION((FeatureXMLFile()))
{
  ptr = new FeatureXMLFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~FeatureXMLFile()))
{
  delete ptr;
}
END_SECTION

START_SECTION((Size loadSize(const String &filename)))
{
  FeatureMap e;
  FeatureXMLFile dfmap_file;
  //test exception
  TEST_EXCEPTION(Exception::FileNotFound, dfmap_file.loadSize("dummy/dummy.MzData"))
  // real test
  Size r = dfmap_file.loadSize(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"));
  TEST_EQUAL(r, 2);
  // again, to test if reset internally works
  r = dfmap_file.loadSize(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"));
  TEST_EQUAL(r, 2);
}
END_SECTION

START_SECTION((void load(const String &filename, FeatureMap&feature_map)))
{
  TOLERANCE_ABSOLUTE(0.01)

  FeatureMap e;
  FeatureXMLFile dfmap_file;

  //test exception
  TEST_EXCEPTION(Exception::FileNotFound, dfmap_file.load("dummy/dummy.MzData", e))

  // real test
  dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"), e);
  TEST_EQUAL(e.getIdentifier(), "lsid");

  //test DocumentIdentifier addition
  TEST_STRING_EQUAL(e.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"));
  TEST_STRING_EQUAL(FileTypes::typeToName(e.getLoadedFileType()), "featureXML");

  TEST_EQUAL(e.size(), 2)
  TEST_REAL_SIMILAR(e[0].getRT(), 25)
  TEST_REAL_SIMILAR(e[0].getMZ(), 0)
  TEST_REAL_SIMILAR(e[0].getIntensity(), 300)
  TEST_EQUAL(e[0].getMetaValue("stringparametername"), "stringparametervalue")
  TEST_EQUAL((UInt)e[0].getMetaValue("intparametername"), 4)
  TEST_REAL_SIMILAR((double)e[0].getMetaValue("floatparametername"), 4.551)
  TEST_REAL_SIMILAR(e[1].getRT(), 0)
  TEST_REAL_SIMILAR(e[1].getMZ(), 35)
  TEST_REAL_SIMILAR(e[1].getIntensity(), 500)
  //data processing
  TEST_EQUAL(e.getDataProcessing().size(), 2)
  TEST_STRING_EQUAL(e.getDataProcessing()[0].getSoftware().getName(), "Software1")
  TEST_STRING_EQUAL(e.getDataProcessing()[0].getSoftware().getVersion(), "0.91a")
  TEST_EQUAL(e.getDataProcessing()[0].getProcessingActions().size(), 1)
  TEST_EQUAL(e.getDataProcessing()[0].getProcessingActions().count(DataProcessing::DEISOTOPING), 1)
  TEST_STRING_EQUAL(e.getDataProcessing()[0].getMetaValue("name"), "dataProcessing")
  TEST_STRING_EQUAL(e.getDataProcessing()[1].getSoftware().getName(), "Software2")
  TEST_STRING_EQUAL(e.getDataProcessing()[1].getSoftware().getVersion(), "0.92a")
  TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().size(), 2)
  TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().count(DataProcessing::SMOOTHING), 1)
  TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().count(DataProcessing::BASELINE_REDUCTION), 1)
  //protein identifications
  TEST_EQUAL(e.getProteinIdentifications().size(), 2)
  TEST_EQUAL(e.getProteinIdentifications()[0].getHits().size(), 2)
  TEST_EQUAL(e.getProteinIdentifications()[0].getHits()[0].getSequence(), "ABCDEFG")
  TEST_EQUAL(e.getProteinIdentifications()[0].getHits()[1].getSequence(), "HIJKLMN")
  TEST_EQUAL(e.getProteinIdentifications()[1].getHits().size(), 1)
  TEST_EQUAL(e.getProteinIdentifications()[1].getHits()[0].getSequence(), "OPQREST")
  //peptide identifications
  TEST_EQUAL(e[0].getPeptideIdentifications().size(), 2)
  TEST_EQUAL(e[0].getPeptideIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(e[0].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("A"))
  TEST_EQUAL(e[0].getPeptideIdentifications()[1].getHits().size(), 2)
  TEST_EQUAL(e[0].getPeptideIdentifications()[1].getHits()[0].getSequence(), AASequence::fromString("C"))
  TEST_EQUAL(e[0].getPeptideIdentifications()[1].getHits()[1].getSequence(), AASequence::fromString("D"))
  TEST_EQUAL(e[1].getPeptideIdentifications().size(), 1)
  TEST_EQUAL(e[1].getPeptideIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(e[1].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("E"))
  //unassigned peptide identifications
  TEST_EQUAL(e.getUnassignedPeptideIdentifications().size(), 2)
  TEST_EQUAL(e.getUnassignedPeptideIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(e.getUnassignedPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("F"))
  TEST_EQUAL(e.getUnassignedPeptideIdentifications()[1].getHits().size(), 2)
  TEST_EQUAL(e.getUnassignedPeptideIdentifications()[1].getHits()[0].getSequence(), AASequence::fromString("G"))
  TEST_EQUAL(e.getUnassignedPeptideIdentifications()[1].getHits()[1].getSequence(), AASequence::fromString("H"))

  // test meta values:
  TEST_EQUAL(e[0].getMetaValue("myIntList") == ListUtils::create<Int>("1,10,12"), true);
  TEST_EQUAL(e[0].getMetaValue("myDoubleList") == ListUtils::create<double>("1.111,10.999,12.45"), true);
  TEST_EQUAL(e[0].getMetaValue("myStringList") == ListUtils::create<String>("myABC1,Stuff,12"), true);
  TEST_EQUAL(e[1].getMetaValue("myDoubleList") == ListUtils::create<double>("6.442"), true);

  //test if loading a second file works (initialization)
  FeatureMap e2;
  dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"), e2);
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  dfmap_file.store(tmp_filename, e2);
  TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"), tmp_filename)

  //test of old file with mzData description (version 1.2)
  //here only the downward-compatibility of the new parser is tested
  //no exception should be thrown
  dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_3_old.featureXML"), e);
  TEST_EQUAL(e.size(), 1)

  //FeatureFileOptions tests
  dfmap_file.getOptions().setRTRange(makeRange(0, 10));
  dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"), e);
  TEST_EQUAL(e.size(), 1)
  TEST_REAL_SIMILAR(e[0].getRT(), 0)
  TEST_REAL_SIMILAR(e[0].getMZ(), 35)
  TEST_REAL_SIMILAR(e[0].getIntensity(), 500)

  dfmap_file.getOptions() = FeatureFileOptions();
  dfmap_file.getOptions().setMZRange(makeRange(10, 50));
  dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"), e);
  TEST_EQUAL(e.size(), 1)
  TEST_REAL_SIMILAR(e[0].getRT(), 0)
  TEST_REAL_SIMILAR(e[0].getMZ(), 35)
  TEST_REAL_SIMILAR(e[0].getIntensity(), 500)

  dfmap_file.getOptions() = FeatureFileOptions();
  dfmap_file.getOptions().setIntensityRange(makeRange(400, 600));
  dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"), e);
  TEST_EQUAL(e.size(), 1)
  TEST_REAL_SIMILAR(e[0].getRT(), 0)
  TEST_REAL_SIMILAR(e[0].getMZ(), 35)
  TEST_REAL_SIMILAR(e[0].getIntensity(), 500)
  {
    // convex hulls:
    dfmap_file.getOptions() = FeatureFileOptions();
    FeatureMap e_full;
    dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"), e_full);
    dfmap_file.getOptions().setLoadConvexHull(false);
    dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"), e);
    // delete CH's manually
    std::vector<ConvexHull2D> empty_hull;
    for (Size ic = 0; ic < e_full.size(); ++ic)
      e_full[ic].setConvexHulls(empty_hull);
    e_full.updateRanges();
    e.updateRanges();
    TEST_EQUAL(e_full, e)
  }

  // subordinates:
  {
    dfmap_file.getOptions() = FeatureFileOptions();
    FeatureMap e_full;
    dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"), e_full);
    dfmap_file.getOptions().setLoadSubordinates(false);
    dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"), e);
    // delete SO's manually
    std::vector<Feature> empty_f;
    for (Size ic = 0; ic < e_full.size(); ++ic)
      e_full[ic].setSubordinates(empty_f);
    TEST_EQUAL(e_full, e)
  }
}
END_SECTION

START_SECTION((void store(const String &filename, const FeatureMap&feature_map)))
{
  FeatureMap map;
  FeatureXMLFile f;

  f.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"), map);
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  f.store(tmp_filename, map);
  WHITELIST("?xml-stylesheet")
  TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"), tmp_filename)
}
END_SECTION

START_SECTION((FeatureFileOptions & getOptions()))
{
  FeatureXMLFile f;
  FeatureMap e;
  f.getOptions().setRTRange(makeRange(1.5, 4.5));
  f.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"), e);
  TEST_EQUAL(e.size(), 5)

  f.getOptions().setMZRange(makeRange(1025.0, 2000.0));
  f.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"), e);
  TEST_EQUAL(e.size(), 3)

  f.getOptions().setIntensityRange(makeRange(290.0, 310.0));
  f.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"), e);
  TEST_EQUAL(e.size(), 1)

  f.getOptions().setMetadataOnly(true);
  f.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"), e);
  TEST_EQUAL(e.getIdentifier(), "lsid2")
  TEST_EQUAL(e.size(), 0)
}
END_SECTION

START_SECTION([EXTRA] static bool isValid(const String& filename))
{
  FeatureXMLFile f;
  TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"), std::cerr), true);
  TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"), std::cerr), true);

  FeatureMap e;
  String filename;

  //test if empty file is valid
  NEW_TMP_FILE(filename)
  f.store(filename, e);
  TEST_EQUAL(f.isValid(filename, std::cerr), true);

  //test if full file is valid
  NEW_TMP_FILE(filename);
  f.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"), e);
  f.store(filename, e);
  TEST_EQUAL(f.isValid(filename, std::cerr), true);
}
END_SECTION

START_SECTION((const FeatureFileOptions &getOptions() const))
{
  FeatureXMLFile f;
  f.getOptions().setRTRange(makeRange(1.5, 4.5));
  f.getOptions().setIntensityRange(makeRange(290.0, 310.0));

  const FeatureFileOptions pfo = f.getOptions();

  TEST_EQUAL(pfo.getRTRange(), makeRange(1.5, 4.5))
  TEST_EQUAL(pfo.getIntensityRange(), makeRange(290.0, 310.0))
}
END_SECTION

START_SECTION(( void setOptions(const FeatureFileOptions &) ))
{
  FeatureXMLFile f;
  FeatureFileOptions pfo = f.getOptions();
  pfo.setMetadataOnly(true);
  pfo.setLoadConvexHull(false);
  pfo.setRTRange(makeRange(1.5, 4.5));
  pfo.setIntensityRange(makeRange(290.0, 310.0));

  f.setOptions(pfo);
  TEST_EQUAL(pfo.getMetadataOnly(), f.getOptions().getMetadataOnly())
  TEST_EQUAL(pfo.getLoadConvexHull(), f.getOptions().getLoadConvexHull())
  TEST_EQUAL(pfo.getRTRange(), f.getOptions().getRTRange())
  TEST_EQUAL(pfo.getIntensityRange(), f.getOptions().getIntensityRange())
}
END_SECTION

START_SECTION([EXTRA])
{
  Feature f1;
  f1.setRT(1001);
  f1.setMZ(1002);
  f1.setCharge(1003);
  Feature f1_cpy(f1);
  Feature f11;
  f11.setRT(1101);
  f11.setMZ(1102);
  Feature f12;
  f12.setRT(1201);
  f12.setMZ(1202);
  Feature f13;
  f13.setRT(1301);
  f13.setMZ(1302);
  TEST_EQUAL(f1.getSubordinates().empty(), true);
  f1.getSubordinates().push_back(f11);
  TEST_EQUAL(f1.getSubordinates().size(), 1);
  f1.getSubordinates().push_back(f12);
  TEST_EQUAL(f1.getSubordinates().size(), 2);
  f1.getSubordinates().push_back(f13);
  TEST_EQUAL(f1.getSubordinates().size(), 3);
  TEST_EQUAL(f1.getRT(), 1001);
  TEST_EQUAL(f1.getSubordinates()[0].getRT(), 1101);
  TEST_EQUAL(f1.getSubordinates()[1].getRT(), 1201);
  TEST_EQUAL(f1.getSubordinates()[2].getRT(), 1301);
  const Feature& f1_cref = f1;
  TEST_EQUAL(f1_cref.getMZ(), 1002);
  TEST_EQUAL(f1_cref.getSubordinates()[0].getMZ(), 1102);
  TEST_EQUAL(f1_cref.getSubordinates()[1].getMZ(), 1202);
  TEST_EQUAL(f1_cref.getSubordinates()[2].getMZ(), 1302);
  TEST_NOT_EQUAL(f1_cref, f1_cpy);
  Feature f1_cpy2(f1);
  TEST_EQUAL(f1_cpy2, f1);
  f1.getSubordinates().clear();
  TEST_EQUAL(f1_cref, f1_cpy);

  Feature f2;
  f2.setRT(1001);
  f2.setMZ(1002);
  f2.setCharge(1003);
  TEST_NOT_EQUAL(f1_cpy2.getSubordinates().empty(), true);
  f2.setSubordinates(f1_cpy2.getSubordinates());
  TEST_EQUAL(f2, f1_cpy2);

  String filename;
  NEW_TMP_FILE(filename);
  FeatureXMLFile f;
  FeatureMap e;
  e.push_back(f1);
  e.push_back(f2);

  // this will print the number of newly assigned unique ids
  STATUS(e.applyMemberFunction(&UniqueIdInterface::ensureUniqueId));

  f.store(filename, e);
  FeatureMap e2;
  f.load(filename, e2);
  e.updateRanges();
  TEST_TRUE(e == e2);
  String filename2;
  NEW_TMP_FILE(filename2);
  f.store(filename2, e2);

}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
