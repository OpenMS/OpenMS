// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/ChromatogramSettings.h>
///////////////////////////

#include <unordered_set>
#include <unordered_map>

using namespace OpenMS;
using namespace std;

START_TEST(ChromatogramSettings, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ChromatogramSettings* ptr = nullptr;
ChromatogramSettings* nullPointer = nullptr;
START_SECTION(ChromatogramSettings())
{
  ptr = new ChromatogramSettings();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~ChromatogramSettings())
{
  delete ptr;
}
END_SECTION

START_SECTION((ChromatogramSettings(const ChromatogramSettings &source)))
{
  ChromatogramSettings tmp;
  tmp.getAcquisitionInfo().setMethodOfCombination("test");
  tmp.getInstrumentSettings().getScanWindows().resize(1);
  tmp.getPrecursor().setMZ(0.11);
  tmp.getProduct().setMZ(0.12);
  tmp.setChromatogramType(ChromatogramSettings::ChromatogramType::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  tmp.setComment("bla");
  tmp.setNativeID("nid");
  tmp.getDataProcessing().resize(1);
  tmp.setMetaValue("bla","bluff");

  ChromatogramSettings tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getComment(), "bla");
  TEST_EQUAL(tmp2.getChromatogramType(), ChromatogramSettings::ChromatogramType::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  TEST_REAL_SIMILAR(tmp2.getPrecursor().getMZ(), 0.11);
  TEST_REAL_SIMILAR(tmp2.getProduct().getMZ(), 0.12);
  TEST_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), false);  
  TEST_EQUAL(tmp2.getAcquisitionInfo()==AcquisitionInfo(), false);
  TEST_EQUAL(tmp2.getAcquisitionInfo().empty(), true);
  TEST_STRING_EQUAL(tmp2.getNativeID(),"nid");
  TEST_EQUAL(tmp2.getDataProcessing().size(),1);
  TEST_STRING_EQUAL(tmp2.getMetaValue("bla"),"bluff");


  tmp2 = ChromatogramSettings();
  TEST_EQUAL(tmp2.getComment(), "");
  TEST_EQUAL(tmp2.getChromatogramType(), ChromatogramSettings::ChromatogramType::MASS_CHROMATOGRAM);
  TEST_REAL_SIMILAR(tmp2.getPrecursor().getMZ(), 0.0);
  TEST_REAL_SIMILAR(tmp2.getProduct().getMZ(), 0.0);
  TEST_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), true);
  TEST_EQUAL(tmp2.getAcquisitionInfo().empty(), true);
  TEST_STRING_EQUAL(tmp2.getNativeID(),"");
  TEST_EQUAL(tmp2.getDataProcessing().size(),0);
  TEST_EQUAL(tmp2.metaValueExists("bla"),false);

}
END_SECTION

START_SECTION((ChromatogramSettings& operator=(const ChromatogramSettings &source)))
{
  ChromatogramSettings tmp;
  tmp.setMetaValue("bla","bluff");
  tmp.getAcquisitionInfo().setMethodOfCombination("test");
  tmp.getInstrumentSettings().getScanWindows().resize(1);
  tmp.getPrecursor().setMZ(0.13);
  tmp.getProduct().setMZ(0.14);
  tmp.setChromatogramType(ChromatogramSettings::ChromatogramType::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  tmp.setComment("bla");
  tmp.setNativeID("nid");
  tmp.getDataProcessing().resize(1);

  ChromatogramSettings tmp2(tmp);
  TEST_EQUAL(tmp2.getComment(), "bla");
  TEST_EQUAL(tmp2.getChromatogramType(), ChromatogramSettings::ChromatogramType::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  TEST_REAL_SIMILAR(tmp2.getPrecursor().getMZ(), 0.13);
  TEST_REAL_SIMILAR(tmp2.getProduct().getMZ(), 0.14);
  TEST_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), false);
  TEST_EQUAL(tmp2.getAcquisitionInfo()==AcquisitionInfo(), false);
  TEST_EQUAL(tmp2.getAcquisitionInfo().empty(), true);
  TEST_STRING_EQUAL(tmp2.getNativeID(),"nid");
  TEST_EQUAL(tmp2.getDataProcessing().size(),1);
  TEST_EQUAL(tmp2.getMetaValue("bla")=="bluff",true);
}
END_SECTION

START_SECTION((bool operator==(const ChromatogramSettings &rhs) const ))
{
  ChromatogramSettings edit, empty;

  TEST_TRUE(edit == empty);

  edit.getAcquisitionInfo().setMethodOfCombination("test");
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.setNativeID("nid");
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.getInstrumentSettings().getScanWindows().resize(1);
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.setComment("comment");
  TEST_EQUAL(edit == empty, false)

  edit = empty;
  edit.getPrecursor().setMZ(0.15);
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.setChromatogramType(ChromatogramSettings::ChromatogramType::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.setComment("bla");
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.getPrecursor().setMZ(0.16);
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.getProduct().setMZ(0.17);
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.getDataProcessing().resize(1);
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.setMetaValue("bla","bluff");
  TEST_EQUAL(edit==empty, false);


}
END_SECTION

START_SECTION((bool operator!=(const ChromatogramSettings &rhs) const ))
{
  ChromatogramSettings edit, empty;

  TEST_EQUAL(edit!=empty, false);

  edit.getAcquisitionInfo().setMethodOfCombination("test");
  TEST_FALSE(edit == empty);

  edit = empty;
  edit.setNativeID("nid");
  TEST_FALSE(edit == empty)

  edit = empty;
  edit.getInstrumentSettings().getScanWindows().resize(1);
  TEST_FALSE(edit == empty);

  edit = empty;
  edit.setComment("comment");
  TEST_FALSE(edit == empty)

  edit = empty;
  edit.setChromatogramType(ChromatogramSettings::ChromatogramType::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  TEST_FALSE(edit == empty);

  edit = empty;
  edit.setComment("bla");
  TEST_FALSE(edit == empty)

  edit = empty;
  Precursor prec;
  prec.setMZ(1.3);
  edit.setPrecursor(prec);
  TEST_FALSE(edit == empty)

  edit = empty;
  Product prod;
  prod.setMZ(1.5);
  edit.setProduct(prod);
  TEST_FALSE(edit == empty)

  edit = empty;
  edit.getDataProcessing().resize(1);
  TEST_FALSE(edit == empty)

  edit = empty;
  edit.setMetaValue("bla","bluff");
  TEST_FALSE(edit == empty)

}
END_SECTION

START_SECTION((const String& getNativeID() const ))
{
  ChromatogramSettings tmp;
  TEST_STRING_EQUAL(tmp.getNativeID(),"")
}
END_SECTION

START_SECTION((void setNativeID(const String &native_id)))
{
  ChromatogramSettings tmp;
  tmp.setNativeID("nid");
  TEST_STRING_EQUAL(tmp.getNativeID(),"nid")
}
END_SECTION

START_SECTION((const String& getComment() const ))
{
  ChromatogramSettings tmp;
  TEST_STRING_EQUAL(tmp.getComment(), "")
}
END_SECTION

START_SECTION((void setComment(const String &comment)))
{
  ChromatogramSettings tmp;
  tmp.setComment("name");
  TEST_STRING_EQUAL(tmp.getComment(), "name")
}
END_SECTION

START_SECTION((const InstrumentSettings& getInstrumentSettings() const ))
{
  ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), true);
}
END_SECTION

START_SECTION((InstrumentSettings& getInstrumentSettings()))
{
  ChromatogramSettings tmp;
  tmp.getInstrumentSettings().getScanWindows().resize(1);
  TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), false);
}
END_SECTION

START_SECTION((void setInstrumentSettings(const InstrumentSettings &instrument_settings)))
{
  ChromatogramSettings tmp;
  InstrumentSettings is;
  is.getScanWindows().resize(1);
  tmp.setInstrumentSettings(is);
  TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), false);
}
END_SECTION

START_SECTION((const AcquisitionInfo& getAcquisitionInfo() const ))
{
  ChromatogramSettings tmp;
  tmp.getAcquisitionInfo().setMethodOfCombination("test");
  TEST_EQUAL(tmp.getAcquisitionInfo()==AcquisitionInfo(), false);  
  TEST_EQUAL(tmp.getAcquisitionInfo().empty(), true);
}
END_SECTION

START_SECTION((AcquisitionInfo& getAcquisitionInfo()))
{
  ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getAcquisitionInfo().empty(), true);
}
END_SECTION

START_SECTION((void setAcquisitionInfo(const AcquisitionInfo &acquisition_info)))
{
  ChromatogramSettings tmp;
  AcquisitionInfo ai;
  ai.setMethodOfCombination("test");
  tmp.setAcquisitionInfo(ai);
  TEST_EQUAL(tmp.getAcquisitionInfo()==AcquisitionInfo(), false);  
  TEST_EQUAL(tmp.getAcquisitionInfo().empty(), true);  
}
END_SECTION

START_SECTION((const SourceFile& getSourceFile() const ))
{
  ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), true);
}
END_SECTION

START_SECTION((SourceFile& getSourceFile()))
{
  ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getSourceFile()==SourceFile(), true);
}
END_SECTION

START_SECTION((void setSourceFile(const SourceFile &source_file)))
{
  ChromatogramSettings tmp;
  SourceFile sf;
  sf.setNameOfFile("test");
  tmp.setSourceFile(sf);
  TEST_EQUAL(tmp.getSourceFile()==SourceFile(), false);
}
END_SECTION

START_SECTION((const Precursor& getPrecursor() const ))
{
  ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getPrecursor() == Precursor(), true)
}
END_SECTION

START_SECTION((Precursor& getPrecursor()))
{
  ChromatogramSettings tmp;
  tmp.getPrecursor().setMZ(0.3);
  TEST_EQUAL(tmp.getPrecursor() == Precursor(), false)
  TEST_REAL_SIMILAR(tmp.getPrecursor().getMZ(), 0.3)
}
END_SECTION

START_SECTION((void setPrecursor(const Precursor &precursor)))
{
  ChromatogramSettings tmp;
  Precursor prec;
  prec.setMZ(0.4);
  tmp.setPrecursor(prec);
  TEST_EQUAL(tmp.getPrecursor() == Precursor(), false)
  TEST_REAL_SIMILAR(tmp.getPrecursor().getMZ(), 0.4)
}
END_SECTION

START_SECTION((const Product& getProduct() const ))
{
  ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getProduct() == Product(), true)
}
END_SECTION

START_SECTION((Product& getProduct()))
{
  ChromatogramSettings tmp;
  tmp.getProduct().setMZ(0.3);
  TEST_EQUAL(tmp.getProduct() == Product(), false)
  TEST_REAL_SIMILAR(tmp.getProduct().getMZ(), 0.3)
}
END_SECTION

START_SECTION((void setProduct(const Product &product)))
{
  ChromatogramSettings tmp;
  Product prod;
  prod.setMZ(0.4);
  tmp.setProduct(prod);
  TEST_EQUAL(tmp.getProduct() == Product(), false)
  TEST_REAL_SIMILAR(tmp.getProduct().getMZ(), 0.4)
}
END_SECTION

START_SECTION((const std::vector<DataProcessing>& getDataProcessing() const ))
{
  ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getDataProcessing().size(),0);
}
END_SECTION

START_SECTION((std::vector<DataProcessing>& getDataProcessing()))
{
  ChromatogramSettings tmp;
  DataProcessingPtr dp = std::shared_ptr<DataProcessing>(new DataProcessing); 
  tmp.getDataProcessing().push_back(dp);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
}
END_SECTION

START_SECTION((void setDataProcessing(const std::vector< DataProcessing > &data_processing)))
{
  ChromatogramSettings tmp;
  std::vector<DataProcessingPtr > dummy;
  dummy.resize(1);
  tmp.setDataProcessing(dummy);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
}
END_SECTION

START_SECTION((ChromatogramType getChromatogramType() const ))
{
  ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getChromatogramType(), ChromatogramSettings::ChromatogramType::MASS_CHROMATOGRAM)
}
END_SECTION

START_SECTION((void setChromatogramType(ChromatogramType type)))
{
  ChromatogramSettings tmp;
  tmp.setChromatogramType(ChromatogramSettings::ChromatogramType::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  TEST_EQUAL(tmp.getChromatogramType(), ChromatogramSettings::ChromatogramType::SELECTED_REACTION_MONITORING_CHROMATOGRAM)
}
END_SECTION


START_SECTION([EXTRA](ENUMs))
{
  // extra stuff tested here:
  TEST_EQUAL(ChromatogramSettings::SIZE_OF_CHROMATOGRAM_TYPE+1, sizeof( ChromatogramSettings::ChromatogramNames ) / sizeof( char* ))
  TEST_EQUAL(String(ChromatogramSettings::ChromatogramNames[ChromatogramSettings::ChromatogramType::MASS_CHROMATOGRAM]), String("mass chromatogram"))
  TEST_EQUAL(String(ChromatogramSettings::ChromatogramNames[ChromatogramSettings::ChromatogramType::EMISSION_CHROMATOGRAM]), String("emission chromatogram"))
  TEST_EQUAL(String(ChromatogramSettings::ChromatogramNames[ChromatogramSettings::SIZE_OF_CHROMATOGRAM_TYPE]), String("unknown chromatogram")) // should be the last entry
}
END_SECTION

START_SECTION([EXTRA] std::hash<ChromatogramSettings>)
{
  std::hash<ChromatogramSettings> hasher;

  // Test that equal objects have equal hashes
  ChromatogramSettings cs1, cs2;
  TEST_EQUAL(cs1 == cs2, true)
  TEST_EQUAL(hasher(cs1), hasher(cs2))

  // Test with populated objects
  cs1.setNativeID("native_id_1");
  cs1.setComment("test comment");
  cs1.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  cs1.getPrecursor().setMZ(500.5);
  cs1.getProduct().setMZ(200.2);
  cs1.getAcquisitionInfo().setMethodOfCombination("sum");
  cs1.getInstrumentSettings().getScanWindows().resize(1);
  cs1.getInstrumentSettings().getScanWindows()[0].begin = 100.0;
  cs1.getInstrumentSettings().getScanWindows()[0].end = 1000.0;

  cs2.setNativeID("native_id_1");
  cs2.setComment("test comment");
  cs2.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  cs2.getPrecursor().setMZ(500.5);
  cs2.getProduct().setMZ(200.2);
  cs2.getAcquisitionInfo().setMethodOfCombination("sum");
  cs2.getInstrumentSettings().getScanWindows().resize(1);
  cs2.getInstrumentSettings().getScanWindows()[0].begin = 100.0;
  cs2.getInstrumentSettings().getScanWindows()[0].end = 1000.0;

  TEST_EQUAL(cs1 == cs2, true)
  TEST_EQUAL(hasher(cs1), hasher(cs2))

  // Test that different objects produce different hashes (not guaranteed but highly likely)
  ChromatogramSettings cs3;
  cs3.setNativeID("different_id");
  TEST_EQUAL(cs1 == cs3, false)
  TEST_NOT_EQUAL(hasher(cs1), hasher(cs3))

  // Test use in unordered_set
  {
    std::unordered_set<ChromatogramSettings> set;
    ChromatogramSettings s1, s2, s3;
    s1.setNativeID("id1");
    s2.setNativeID("id2");
    s3.setNativeID("id1"); // same as s1

    set.insert(s1);
    set.insert(s2);
    TEST_EQUAL(set.size(), 2)

    // s3 is equal to s1, so size should stay 2
    set.insert(s3);
    TEST_EQUAL(set.size(), 2)

    TEST_EQUAL(set.count(s1), 1)
    TEST_EQUAL(set.count(s2), 1)
    TEST_EQUAL(set.count(s3), 1) // s3 == s1
  }

  // Test use in unordered_map
  {
    std::unordered_map<ChromatogramSettings, int> map;
    ChromatogramSettings k1, k2, k3;
    k1.setNativeID("key1");
    k2.setNativeID("key2");
    k3.setNativeID("key1"); // same as k1

    map[k1] = 1;
    map[k2] = 2;
    TEST_EQUAL(map.size(), 2)

    TEST_EQUAL(map[k1], 1)
    TEST_EQUAL(map[k2], 2)
    TEST_EQUAL(map[k3], 1) // k3 == k1

    map[k3] = 3; // should update k1's value
    TEST_EQUAL(map[k1], 3)
  }

  // Test hash consistency (same object hashed multiple times)
  size_t hash1 = hasher(cs1);
  size_t hash2 = hasher(cs1);
  TEST_EQUAL(hash1, hash2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



