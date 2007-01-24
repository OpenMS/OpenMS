// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/InstrumentSettings.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(InstrumentSettings, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

InstrumentSettings* ptr = 0;
CHECK(InstrumentSettings())
	ptr = new InstrumentSettings();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~InstrumentSettings())
	delete ptr;
RESULT

CHECK(IonSource::Polarity getPolarity() const)
	InstrumentSettings tmp;
	TEST_EQUAL(tmp.getPolarity(),IonSource::POLNULL);
RESULT

CHECK(void setPolarity(IonSource::Polarity polarity))
	InstrumentSettings tmp;
	tmp.setPolarity(IonSource::NEGATIVE);
	TEST_EQUAL(tmp.getPolarity(),IonSource::NEGATIVE);
RESULT

CHECK(float getMzRangeStart() const)
	InstrumentSettings tmp;
	TEST_REAL_EQUAL(tmp.getMzRangeStart(),0.0);
RESULT

CHECK(void setMzRangeStart(float mz_range_start))
	InstrumentSettings tmp;
	tmp.setMzRangeStart(47.11);
	TEST_REAL_EQUAL(tmp.getMzRangeStart(),47.11);
RESULT

CHECK(float getMzRangeStop() const)
	InstrumentSettings tmp;
	TEST_REAL_EQUAL(tmp.getMzRangeStop(),0.0);
RESULT

CHECK(void setMzRangeStop(float mz_range_stop))
	InstrumentSettings tmp;
	tmp.setMzRangeStop(47.11);
	TEST_REAL_EQUAL(tmp.getMzRangeStop(),47.11);
RESULT

CHECK(ScanMode getScanMode() const)
	InstrumentSettings tmp;
	TEST_EQUAL(tmp.getScanMode(),InstrumentSettings::SCANMODENULL);
RESULT

CHECK(void setScanMode(ScanMode scan_mode))
	InstrumentSettings tmp;
	tmp.setScanMode(InstrumentSettings::SELECTEDIONDETECTION);
	TEST_EQUAL(tmp.getScanMode(),InstrumentSettings::SELECTEDIONDETECTION);
RESULT

CHECK(InstrumentSettings(const InstrumentSettings& source))
  InstrumentSettings tmp;
  tmp.setScanMode(InstrumentSettings::SELECTEDIONDETECTION);
  tmp.setMzRangeStart(47.11);
  tmp.setMzRangeStop(47.12);
  tmp.setPolarity(IonSource::NEGATIVE);
  tmp.setMetaValue("label",String("label"));
  
  InstrumentSettings tmp2(tmp);
  TEST_EQUAL(tmp2.getScanMode(),InstrumentSettings::SELECTEDIONDETECTION);
  TEST_REAL_EQUAL(tmp2.getMzRangeStart(),47.11);
  TEST_REAL_EQUAL(tmp2.getMzRangeStop(),47.12);
  TEST_EQUAL(tmp2.getPolarity(),IonSource::NEGATIVE);  
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");  
RESULT

CHECK(InstrumentSettings& operator= (const InstrumentSettings& source))
  InstrumentSettings tmp;
  tmp.setScanMode(InstrumentSettings::SELECTEDIONDETECTION);
  tmp.setMzRangeStart(47.11);
  tmp.setMzRangeStop(47.12);
  tmp.setPolarity(IonSource::NEGATIVE);
  tmp.setMetaValue("label",String("label"));
  
  InstrumentSettings tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getScanMode(),InstrumentSettings::SELECTEDIONDETECTION);
  TEST_REAL_EQUAL(tmp2.getMzRangeStart(),47.11);
  TEST_REAL_EQUAL(tmp2.getMzRangeStop(),47.12);
  TEST_EQUAL(tmp2.getPolarity(),IonSource::NEGATIVE);  
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");  
  
  tmp2 = InstrumentSettings();
  TEST_EQUAL(tmp2.getScanMode(),InstrumentSettings::SCANMODENULL);
  TEST_REAL_EQUAL(tmp2.getMzRangeStart(),0.0);
  TEST_REAL_EQUAL(tmp2.getMzRangeStop(),0.0);
  TEST_EQUAL(tmp2.getPolarity(),IonSource::POLNULL);  
	TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
RESULT

CHECK(bool operator== (const InstrumentSettings& rhs) const)
  InstrumentSettings edit, empty;
  
  TEST_EQUAL(edit==empty,true);
  
  edit.setScanMode(InstrumentSettings::SELECTEDIONDETECTION);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;	
  edit.setMzRangeStart(47.11);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setMzRangeStop(47.12);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setPolarity(IonSource::NEGATIVE);
	TEST_EQUAL(edit==empty,false);
	
	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit==empty,false);
RESULT

CHECK(bool operator!= (const InstrumentSettings& rhs) const)
  InstrumentSettings edit, empty;
  
  TEST_EQUAL(edit!=empty,false);
  
  edit.setScanMode(InstrumentSettings::SELECTEDIONDETECTION);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;	
  edit.setMzRangeStart(47.11);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setMzRangeStop(47.12);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setPolarity(IonSource::NEGATIVE);
	TEST_EQUAL(edit!=empty,true);

	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty,true);
RESULT



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



