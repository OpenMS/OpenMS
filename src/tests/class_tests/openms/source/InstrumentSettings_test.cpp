// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/InstrumentSettings.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(InstrumentSettings, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

InstrumentSettings* ptr = nullptr;
InstrumentSettings* nullPointer = nullptr;
START_SECTION((InstrumentSettings()))
	ptr = new InstrumentSettings();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~InstrumentSettings()))
	delete ptr;
END_SECTION

START_SECTION((IonSource::Polarity getPolarity() const))
	InstrumentSettings tmp;
	TEST_EQUAL(tmp.getPolarity(),IonSource::POLNULL);
END_SECTION

START_SECTION((void setPolarity(IonSource::Polarity polarity)))
	InstrumentSettings tmp;
	tmp.setPolarity(IonSource::NEGATIVE);
	TEST_EQUAL(tmp.getPolarity(),IonSource::NEGATIVE);
END_SECTION

START_SECTION((const std::vector< ScanWindow >&  getScanWindows() const))
	InstrumentSettings tmp;
  TEST_EQUAL(tmp.getScanWindows().size(),0);
END_SECTION

START_SECTION((std::vector< ScanWindow >& getScanWindows()))
	InstrumentSettings tmp;
  tmp.getScanWindows().resize(1);
  TEST_EQUAL(tmp.getScanWindows().size(),1);
END_SECTION

START_SECTION((void setScanWindows(std::vector< ScanWindow >  scan_windows)))
	InstrumentSettings tmp;
  vector<ScanWindow> vec(17);
  tmp.setScanWindows(vec);
  TEST_EQUAL(tmp.getScanWindows().size(),17);
END_SECTION

START_SECTION((ScanMode getScanMode() const))
	InstrumentSettings tmp;
	TEST_EQUAL(tmp.getScanMode(),InstrumentSettings::UNKNOWN);
END_SECTION

START_SECTION((void setScanMode(ScanMode scan_mode)))
	InstrumentSettings tmp;
	tmp.setScanMode(InstrumentSettings::SIM);
	TEST_EQUAL(tmp.getScanMode(),InstrumentSettings::SIM);
END_SECTION

START_SECTION((bool getZoomScan() const))
	InstrumentSettings tmp;
	TEST_EQUAL(tmp.getZoomScan(),false);
END_SECTION

START_SECTION((void setZoomScan(bool zoom_scan)))
	InstrumentSettings tmp;
	tmp.setZoomScan(true);
	TEST_EQUAL(tmp.getZoomScan(),true);
END_SECTION

START_SECTION((InstrumentSettings(const InstrumentSettings& source)))
  InstrumentSettings tmp;
  tmp.setScanMode(InstrumentSettings::SIM);
  tmp.getScanWindows().resize(1);
  tmp.setPolarity(IonSource::NEGATIVE);
  tmp.setMetaValue("label",String("label"));
	tmp.setZoomScan(true);
  
  InstrumentSettings tmp2(tmp);
  TEST_EQUAL(tmp2.getScanMode(),InstrumentSettings::SIM);
  TEST_EQUAL(tmp2.getScanWindows().size(),1);
  TEST_EQUAL(tmp2.getPolarity(),IonSource::NEGATIVE);  
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");  
	TEST_EQUAL(tmp2.getZoomScan(),true);
END_SECTION

START_SECTION((InstrumentSettings& operator= (const InstrumentSettings& source)))
  InstrumentSettings tmp;
  tmp.setScanMode(InstrumentSettings::SIM);
  tmp.getScanWindows().resize(1);
  tmp.setPolarity(IonSource::NEGATIVE);
  tmp.setMetaValue("label",String("label"));
	tmp.setZoomScan(true);
  
  InstrumentSettings tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getScanMode(),InstrumentSettings::SIM);
  TEST_EQUAL(tmp2.getScanWindows().size(),1);
  TEST_EQUAL(tmp2.getPolarity(),IonSource::NEGATIVE);  
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_EQUAL(tmp2.getZoomScan(),true);
  
  tmp2 = InstrumentSettings();
  TEST_EQUAL(tmp2.getScanMode(),InstrumentSettings::UNKNOWN);
  TEST_EQUAL(tmp2.getScanWindows().size(),0);
  TEST_EQUAL(tmp2.getPolarity(),IonSource::POLNULL);  
	TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
	TEST_EQUAL(tmp2.getZoomScan(),false);
END_SECTION

START_SECTION((bool operator== (const InstrumentSettings& rhs) const))
  InstrumentSettings edit, empty;
  
  TEST_EQUAL(edit==empty,true);
  
  edit.setScanMode(InstrumentSettings::SIM);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty; 
  edit.getScanWindows().resize(1);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setPolarity(IonSource::NEGATIVE);
	TEST_EQUAL(edit==empty,false);
	
	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit==empty,false);

	edit = empty;
	edit.setZoomScan(true);
	TEST_EQUAL(edit==empty,false);
END_SECTION

START_SECTION((bool operator!= (const InstrumentSettings& rhs) const))
  InstrumentSettings edit, empty;
  
  TEST_EQUAL(edit!=empty,false);
  
  edit.setScanMode(InstrumentSettings::SIM);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;	
  edit.getScanWindows().resize(1);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setPolarity(IonSource::NEGATIVE);
	TEST_EQUAL(edit!=empty,true);

	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty,true);

	edit = empty;
	edit.setZoomScan(true);
	TEST_EQUAL(edit!=empty,true);
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



