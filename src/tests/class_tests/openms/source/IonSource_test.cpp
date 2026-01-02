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
#include <OpenMS/METADATA/IonSource.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IonSource, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IonSource* ptr = nullptr;
IonSource* nullPointer = nullptr;
START_SECTION((IonSource()))
	ptr = new IonSource();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~IonSource()))
	delete ptr;
END_SECTION

START_SECTION(Int getOrder() const)
	IonSource tmp;
	TEST_EQUAL(tmp.getOrder(),0)
END_SECTION

START_SECTION(void setOrder(Int order))
	IonSource tmp;
	tmp.setOrder(4711);
	TEST_EQUAL(tmp.getOrder(),4711)
END_SECTION

START_SECTION((InletType getInletType() const))
  IonSource tmp;
  TEST_EQUAL(tmp.getInletType(),IonSource::InletType::INLETNULL);
END_SECTION

START_SECTION((void setInletType(InletType inlet_type)))
  IonSource tmp;
  tmp.setInletType(IonSource::InletType::DIRECT);
  TEST_EQUAL(tmp.getInletType(),IonSource::InletType::DIRECT);
END_SECTION

START_SECTION((IonizationMethod getIonizationMethod() const))
  IonSource tmp;
  TEST_EQUAL(tmp.getIonizationMethod(),IonSource::IonizationMethod::IONMETHODNULL);
END_SECTION

START_SECTION((void setIonizationMethod(IonizationMethod ionization_type)))
  IonSource tmp;
  tmp.setIonizationMethod(IonSource::IonizationMethod::ESI);
  TEST_EQUAL(tmp.getIonizationMethod(),IonSource::IonizationMethod::ESI);
END_SECTION

START_SECTION((Polarity getPolarity() const))
  IonSource tmp;
  TEST_EQUAL(tmp.getPolarity(),IonSource::Polarity::POLNULL);
END_SECTION

START_SECTION((void setPolarity(Polarity polarity)))
	IonSource tmp;
  tmp.setPolarity(IonSource::Polarity::POSITIVE);
  TEST_EQUAL(tmp.getPolarity(),IonSource::Polarity::POSITIVE);
END_SECTION

START_SECTION((IonSource(const IonSource& source)))
  IonSource tmp;
  tmp.setInletType(IonSource::InletType::DIRECT);
  tmp.setIonizationMethod(IonSource::IonizationMethod::ESI);
  tmp.setPolarity(IonSource::Polarity::POSITIVE);
  tmp.setMetaValue("label",String("label"));
  tmp.setOrder(45);
  	
  IonSource tmp2(tmp);
  TEST_EQUAL(tmp2.getPolarity(),IonSource::Polarity::POSITIVE);
  TEST_EQUAL(tmp2.getInletType(),IonSource::InletType::DIRECT);
  TEST_EQUAL(tmp2.getIonizationMethod(),IonSource::IonizationMethod::ESI);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_EQUAL(tmp2.getOrder(),45)
END_SECTION

START_SECTION((IonSource& operator= (const IonSource& source)))
  IonSource tmp;
  tmp.setInletType(IonSource::InletType::DIRECT);
  tmp.setIonizationMethod(IonSource::IonizationMethod::ESI);
  tmp.setPolarity(IonSource::Polarity::POSITIVE);
  tmp.setMetaValue("label",String("label"));
  tmp.setOrder(45);
  
  IonSource tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getPolarity(),IonSource::Polarity::POSITIVE);
  TEST_EQUAL(tmp2.getInletType(),IonSource::InletType::DIRECT);
  TEST_EQUAL(tmp2.getIonizationMethod(),IonSource::IonizationMethod::ESI);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_EQUAL(tmp2.getOrder(),45)
  
  tmp2 = IonSource();
  TEST_EQUAL(tmp2.getPolarity(),IonSource::Polarity::POLNULL);
  TEST_EQUAL(tmp2.getInletType(),IonSource::InletType::INLETNULL);
  TEST_EQUAL(tmp2.getIonizationMethod(),IonSource::IonizationMethod::IONMETHODNULL);
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
	TEST_EQUAL(tmp2.getOrder(),0)
END_SECTION

START_SECTION((bool operator== (const IonSource& rhs) const))
  IonSource edit,empty;
  
  TEST_EQUAL(edit==empty,true);
  
  edit = empty;
  edit.setInletType(IonSource::InletType::DIRECT);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setIonizationMethod(IonSource::IonizationMethod::ESI);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setPolarity(IonSource::Polarity::POSITIVE);
	TEST_EQUAL(edit==empty,false);
	
	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit==empty,false);
	
  edit = empty;
  edit.setOrder(45);
	TEST_EQUAL(edit==empty,false);
END_SECTION

START_SECTION((bool operator!= (const IonSource& rhs) const))
  IonSource edit,empty;
  
  TEST_EQUAL(edit!=empty,false);
  
  edit = empty;
  edit.setInletType(IonSource::InletType::DIRECT);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setIonizationMethod(IonSource::IonizationMethod::ESI);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setPolarity(IonSource::Polarity::POSITIVE);
	TEST_EQUAL(edit!=empty,true);

	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty,true);
	
  edit = empty;
  edit.setOrder(45);
	TEST_EQUAL(edit!=empty,true);
END_SECTION

START_SECTION((static StringList getAllNamesOfInletType()))
  StringList names = IonSource::getAllNamesOfInletType();
  TEST_EQUAL(names.size(), IonSource::SIZE_OF_INLETTYPE);
  TEST_EQUAL(names[IonSource::InletType::INLETNULL], "Unknown");
  TEST_EQUAL(names[IonSource::InletType::DIRECT], "Direct");
  TEST_EQUAL(names[IonSource::InletType::NANOSPRAY], "Nanospray inlet");
END_SECTION

START_SECTION((static StringList getAllNamesOfIonizationMethod()))
  StringList names = IonSource::getAllNamesOfIonizationMethod();
  TEST_EQUAL(names.size(), IonSource::SIZE_OF_IONIZATIONMETHOD);
  TEST_EQUAL(names[IonSource::IonizationMethod::IONMETHODNULL], "Unknown");
  TEST_EQUAL(names[IonSource::IonizationMethod::ESI], "Electrospray ionisation");
  TEST_EQUAL(names[IonSource::IonizationMethod::MALDI], "Matrix-assisted laser desorption ionization");
END_SECTION

START_SECTION((static StringList getAllNamesOfPolarity()))
  StringList names = IonSource::getAllNamesOfPolarity();
  TEST_EQUAL(names.size(), IonSource::SIZE_OF_POLARITY);
  TEST_EQUAL(names[IonSource::Polarity::POLNULL], "unknown");
  TEST_EQUAL(names[IonSource::Polarity::POSITIVE], "positive");
  TEST_EQUAL(names[IonSource::Polarity::NEGATIVE], "negative");
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



