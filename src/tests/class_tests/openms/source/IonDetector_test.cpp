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
#include <unordered_set>
#include <unordered_map>

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
  TEST_EQUAL(tmp.getType(),IonDetector::Type::TYPENULL);
END_SECTION

START_SECTION((void setType(Type type)))
  IonDetector tmp;
  tmp.setType(IonDetector::Type::ELECTRONMULTIPLIER);
  TEST_EQUAL(tmp.getType(),IonDetector::Type::ELECTRONMULTIPLIER);
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
  TEST_EQUAL(tmp.getAcquisitionMode(),IonDetector::AcquisitionMode::ACQMODENULL);
END_SECTION

START_SECTION((void setAcquisitionMode(AcquisitionMode acquisition_mode)))
  IonDetector tmp;
  tmp.setAcquisitionMode(IonDetector::AcquisitionMode::PULSECOUNTING);
  TEST_EQUAL(tmp.getAcquisitionMode(),IonDetector::AcquisitionMode::PULSECOUNTING);
END_SECTION

START_SECTION((IonDetector(const IonDetector& source)))
 	IonDetector tmp;
  tmp.setResolution(47.11);
  tmp.setADCSamplingFrequency(47.21);
  tmp.setAcquisitionMode(IonDetector::AcquisitionMode::PULSECOUNTING);
  tmp.setType(IonDetector::Type::ELECTRONMULTIPLIER);
  tmp.setMetaValue("label",String("label"));
  tmp.setOrder(45);
  
  IonDetector tmp2(tmp);	
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  TEST_REAL_SIMILAR(tmp2.getResolution(),47.11);
  TEST_REAL_SIMILAR(tmp2.getADCSamplingFrequency(),47.21);
  TEST_EQUAL(tmp2.getAcquisitionMode(),IonDetector::AcquisitionMode::PULSECOUNTING);
  TEST_EQUAL(tmp2.getType(),IonDetector::Type::ELECTRONMULTIPLIER);
	TEST_EQUAL(tmp2.getOrder(),45)
END_SECTION

START_SECTION((IonDetector& operator= (const IonDetector& source)))
 	IonDetector tmp;
  tmp.setResolution(47.11);
  tmp.setADCSamplingFrequency(47.21);
  tmp.setAcquisitionMode(IonDetector::AcquisitionMode::PULSECOUNTING);
  tmp.setType(IonDetector::Type::ELECTRONMULTIPLIER);
  tmp.setMetaValue("label",String("label"));
  tmp.setOrder(45);
  
  IonDetector tmp2;
  tmp2 = tmp;	
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  TEST_REAL_SIMILAR(tmp2.getResolution(),47.11);
  TEST_REAL_SIMILAR(tmp2.getADCSamplingFrequency(),47.21);
  TEST_EQUAL(tmp2.getAcquisitionMode(),IonDetector::AcquisitionMode::PULSECOUNTING);
  TEST_EQUAL(tmp2.getType(),IonDetector::Type::ELECTRONMULTIPLIER);
	TEST_EQUAL(tmp2.getOrder(),45)

  tmp2 = IonDetector();	
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
  TEST_REAL_SIMILAR(tmp2.getResolution(),0.0);
  TEST_REAL_SIMILAR(tmp2.getADCSamplingFrequency(),0.0);
  TEST_EQUAL(tmp2.getAcquisitionMode(),IonDetector::AcquisitionMode::ACQMODENULL);
  TEST_EQUAL(tmp2.getType(),IonDetector::Type::TYPENULL);
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
  edit.setAcquisitionMode(IonDetector::AcquisitionMode::PULSECOUNTING);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setType(IonDetector::Type::ELECTRONMULTIPLIER);
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
  edit.setAcquisitionMode(IonDetector::AcquisitionMode::PULSECOUNTING);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setType(IonDetector::Type::ELECTRONMULTIPLIER);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty,true);
	
  edit = empty;
  edit.setOrder(45);
	TEST_EQUAL(edit!=empty,true);
END_SECTION




START_SECTION((static StringList getAllNamesOfType()))
  StringList names = IonDetector::getAllNamesOfType();
  TEST_EQUAL(names.size(), IonDetector::SIZE_OF_TYPE);
  TEST_EQUAL(names[IonDetector::Type::ELECTRONMULTIPLIER], "Electron multiplier");
END_SECTION

START_SECTION((static StringList getAllNamesOfAcquisitionMode()))
  StringList names = IonDetector::getAllNamesOfAcquisitionMode();
  TEST_EQUAL(names.size(), IonDetector::SIZE_OF_ACQUISITIONMODE);
  TEST_EQUAL(names[IonDetector::AcquisitionMode::PULSECOUNTING], "Pulse counting");
END_SECTION

START_SECTION(([EXTRA] std::hash<IonDetector>))
{
  // Test that equal detectors have equal hashes
  IonDetector d1, d2;
  d1.setType(IonDetector::ELECTRONMULTIPLIER);
  d1.setAcquisitionMode(IonDetector::PULSECOUNTING);
  d1.setResolution(47.11);
  d1.setADCSamplingFrequency(47.21);
  d1.setOrder(45);

  d2.setType(IonDetector::ELECTRONMULTIPLIER);
  d2.setAcquisitionMode(IonDetector::PULSECOUNTING);
  d2.setResolution(47.11);
  d2.setADCSamplingFrequency(47.21);
  d2.setOrder(45);

  std::hash<IonDetector> hasher;
  TEST_EQUAL(hasher(d1), hasher(d2))

  // Test that hash changes when values change
  IonDetector d3;
  d3.setType(IonDetector::PHOTOMULTIPLIER);
  d3.setAcquisitionMode(IonDetector::PULSECOUNTING);
  d3.setResolution(47.11);
  d3.setADCSamplingFrequency(47.21);
  d3.setOrder(45);
  TEST_NOT_EQUAL(hasher(d1), hasher(d3))

  IonDetector d4;
  d4.setType(IonDetector::ELECTRONMULTIPLIER);
  d4.setAcquisitionMode(IonDetector::ADC);
  d4.setResolution(47.11);
  d4.setADCSamplingFrequency(47.21);
  d4.setOrder(45);
  TEST_NOT_EQUAL(hasher(d1), hasher(d4))

  IonDetector d5;
  d5.setType(IonDetector::ELECTRONMULTIPLIER);
  d5.setAcquisitionMode(IonDetector::PULSECOUNTING);
  d5.setResolution(100.0);
  d5.setADCSamplingFrequency(47.21);
  d5.setOrder(45);
  TEST_NOT_EQUAL(hasher(d1), hasher(d5))

  IonDetector d6;
  d6.setType(IonDetector::ELECTRONMULTIPLIER);
  d6.setAcquisitionMode(IonDetector::PULSECOUNTING);
  d6.setResolution(47.11);
  d6.setADCSamplingFrequency(100.0);
  d6.setOrder(45);
  TEST_NOT_EQUAL(hasher(d1), hasher(d6))

  IonDetector d7;
  d7.setType(IonDetector::ELECTRONMULTIPLIER);
  d7.setAcquisitionMode(IonDetector::PULSECOUNTING);
  d7.setResolution(47.11);
  d7.setADCSamplingFrequency(47.21);
  d7.setOrder(100);
  TEST_NOT_EQUAL(hasher(d1), hasher(d7))

  // Test use in unordered_set
  std::unordered_set<IonDetector> detector_set;
  detector_set.insert(d1);
  TEST_EQUAL(detector_set.size(), 1)
  detector_set.insert(d2); // same as d1
  TEST_EQUAL(detector_set.size(), 1) // should not increase
  detector_set.insert(d3);
  TEST_EQUAL(detector_set.size(), 2)

  // Test use in unordered_map
  std::unordered_map<IonDetector, int> detector_map;
  detector_map[d1] = 42;
  TEST_EQUAL(detector_map[d1], 42)
  TEST_EQUAL(detector_map[d2], 42) // d2 == d1, should return same value
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



