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
#include <OpenMS/METADATA/IonDetector.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IonDetector, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IonDetector* ptr = nullptr;
IonDetector* nullPointer = nullptr;
START_SECTION((IonDetector()))
	ptr = new IonDetector();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~IonDetector()))
	delete ptr;
END_SECTION

START_SECTION((Int getOrder() const))
	IonDetector tmp;
	TEST_EQUAL(tmp.getOrder(),0)
END_SECTION

START_SECTION((void setOrder(Int order)))
	IonDetector tmp;
	tmp.setOrder(4711);
	TEST_EQUAL(tmp.getOrder(),4711)
END_SECTION

START_SECTION((Type getType() const))
  IonDetector tmp;
  TEST_EQUAL(tmp.getType(),IonDetector::TYPENULL);
END_SECTION

START_SECTION((void setType(Type type)))
  IonDetector tmp;
  tmp.setType(IonDetector::ELECTRONMULTIPLIER);
  TEST_EQUAL(tmp.getType(),IonDetector::ELECTRONMULTIPLIER);
END_SECTION

START_SECTION((double getADCSamplingFrequency() const ))
  IonDetector tmp;
  TEST_EQUAL(tmp.getADCSamplingFrequency(),0);
END_SECTION

START_SECTION((void setADCSamplingFrequency(double ADC_sampling_frequency)))
  IonDetector tmp;
  tmp.setADCSamplingFrequency(47.11);
  TEST_REAL_SIMILAR(tmp.getADCSamplingFrequency(),47.11);
END_SECTION

START_SECTION((double getResolution() const ))
  IonDetector tmp;
  TEST_EQUAL(tmp.getResolution(),0);
END_SECTION

START_SECTION((void setResolution(double resolution)))
  IonDetector tmp;
  tmp.setResolution(47.11);
  TEST_REAL_SIMILAR(tmp.getResolution(),47.11);
END_SECTION

START_SECTION((AcquisitionMode getAcquisitionMode() const))
  IonDetector tmp;
  TEST_EQUAL(tmp.getAcquisitionMode(),IonDetector::ACQMODENULL);
END_SECTION

START_SECTION((void setAcquisitionMode(AcquisitionMode acquisition_mode)))
  IonDetector tmp;
  tmp.setAcquisitionMode(IonDetector::PULSECOUNTING);
  TEST_EQUAL(tmp.getAcquisitionMode(),IonDetector::PULSECOUNTING);
END_SECTION

START_SECTION((IonDetector(const IonDetector& source)))
 	IonDetector tmp;
  tmp.setResolution(47.11);
  tmp.setADCSamplingFrequency(47.21);
  tmp.setAcquisitionMode(IonDetector::PULSECOUNTING);
  tmp.setType(IonDetector::ELECTRONMULTIPLIER);
  tmp.setMetaValue("label",String("label"));
  tmp.setOrder(45);
  
  IonDetector tmp2(tmp);	
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  TEST_REAL_SIMILAR(tmp2.getResolution(),47.11);
  TEST_REAL_SIMILAR(tmp2.getADCSamplingFrequency(),47.21);
  TEST_EQUAL(tmp2.getAcquisitionMode(),IonDetector::PULSECOUNTING);
  TEST_EQUAL(tmp2.getType(),IonDetector::ELECTRONMULTIPLIER);
	TEST_EQUAL(tmp2.getOrder(),45)
END_SECTION

START_SECTION((IonDetector& operator= (const IonDetector& source)))
 	IonDetector tmp;
  tmp.setResolution(47.11);
  tmp.setADCSamplingFrequency(47.21);
  tmp.setAcquisitionMode(IonDetector::PULSECOUNTING);
  tmp.setType(IonDetector::ELECTRONMULTIPLIER);
  tmp.setMetaValue("label",String("label"));
  tmp.setOrder(45);
  
  IonDetector tmp2;
  tmp2 = tmp;	
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  TEST_REAL_SIMILAR(tmp2.getResolution(),47.11);
  TEST_REAL_SIMILAR(tmp2.getADCSamplingFrequency(),47.21);
  TEST_EQUAL(tmp2.getAcquisitionMode(),IonDetector::PULSECOUNTING);
  TEST_EQUAL(tmp2.getType(),IonDetector::ELECTRONMULTIPLIER);
	TEST_EQUAL(tmp2.getOrder(),45)

  tmp2 = IonDetector();	
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
  TEST_REAL_SIMILAR(tmp2.getResolution(),0.0);
  TEST_REAL_SIMILAR(tmp2.getADCSamplingFrequency(),0.0);
  TEST_EQUAL(tmp2.getAcquisitionMode(),IonDetector::ACQMODENULL);
  TEST_EQUAL(tmp2.getType(),IonDetector::TYPENULL);
	TEST_EQUAL(tmp2.getOrder(),0)
END_SECTION

START_SECTION((bool operator== (const IonDetector& rhs) const))
 	IonDetector edit,empty;
 	
 	TEST_EQUAL(edit==empty,true);
 	
  edit.setResolution(47.11);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setADCSamplingFrequency(47.21);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setAcquisitionMode(IonDetector::PULSECOUNTING);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setType(IonDetector::ELECTRONMULTIPLIER);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit==empty,false);

  edit = empty;
  edit.setOrder(45);
	TEST_EQUAL(edit==empty,false);
END_SECTION

START_SECTION((bool operator!= (const IonDetector& rhs) const))
 	IonDetector edit,empty;
 	
 	TEST_EQUAL(edit!=empty,false);
 	
  edit.setResolution(47.11);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setADCSamplingFrequency(47.21);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setAcquisitionMode(IonDetector::PULSECOUNTING);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setType(IonDetector::ELECTRONMULTIPLIER);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty,true);
	
  edit = empty;
  edit.setOrder(45);
	TEST_EQUAL(edit!=empty,true);
END_SECTION




/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



