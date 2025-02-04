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
  TEST_EQUAL(tmp.getInletType(),IonSource::INLETNULL);
END_SECTION

START_SECTION((void setInletType(InletType inlet_type)))
  IonSource tmp;
  tmp.setInletType(IonSource::DIRECT);
  TEST_EQUAL(tmp.getInletType(),IonSource::DIRECT);
END_SECTION

START_SECTION((IonizationMethod getIonizationMethod() const))
  IonSource tmp;
  TEST_EQUAL(tmp.getIonizationMethod(),IonSource::IONMETHODNULL);
END_SECTION

START_SECTION((void setIonizationMethod(IonizationMethod ionization_type)))
  IonSource tmp;
  tmp.setIonizationMethod(IonSource::ESI);
  TEST_EQUAL(tmp.getIonizationMethod(),IonSource::ESI);
END_SECTION

START_SECTION((Polarity getPolarity() const))
  IonSource tmp;
  TEST_EQUAL(tmp.getPolarity(),IonSource::POLNULL);
END_SECTION

START_SECTION((void setPolarity(Polarity polarity)))
	IonSource tmp;
  tmp.setPolarity(IonSource::POSITIVE);
  TEST_EQUAL(tmp.getPolarity(),IonSource::POSITIVE);
END_SECTION

START_SECTION((IonSource(const IonSource& source)))
  IonSource tmp;
  tmp.setInletType(IonSource::DIRECT);
  tmp.setIonizationMethod(IonSource::ESI);
  tmp.setPolarity(IonSource::POSITIVE);
  tmp.setMetaValue("label",String("label"));
  tmp.setOrder(45);
  	
  IonSource tmp2(tmp);
  TEST_EQUAL(tmp2.getPolarity(),IonSource::POSITIVE);
  TEST_EQUAL(tmp2.getInletType(),IonSource::DIRECT);
  TEST_EQUAL(tmp2.getIonizationMethod(),IonSource::ESI);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_EQUAL(tmp2.getOrder(),45)
END_SECTION

START_SECTION((IonSource& operator= (const IonSource& source)))
  IonSource tmp;
  tmp.setInletType(IonSource::DIRECT);
  tmp.setIonizationMethod(IonSource::ESI);
  tmp.setPolarity(IonSource::POSITIVE);
  tmp.setMetaValue("label",String("label"));
  tmp.setOrder(45);
  
  IonSource tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getPolarity(),IonSource::POSITIVE);
  TEST_EQUAL(tmp2.getInletType(),IonSource::DIRECT);
  TEST_EQUAL(tmp2.getIonizationMethod(),IonSource::ESI);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_EQUAL(tmp2.getOrder(),45)
  
  tmp2 = IonSource();
  TEST_EQUAL(tmp2.getPolarity(),IonSource::POLNULL);
  TEST_EQUAL(tmp2.getInletType(),IonSource::INLETNULL);
  TEST_EQUAL(tmp2.getIonizationMethod(),IonSource::IONMETHODNULL);
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
	TEST_EQUAL(tmp2.getOrder(),0)
END_SECTION

START_SECTION((bool operator== (const IonSource& rhs) const))
  IonSource edit,empty;
  
  TEST_EQUAL(edit==empty,true);
  
  edit = empty;
  edit.setInletType(IonSource::DIRECT);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setIonizationMethod(IonSource::ESI);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setPolarity(IonSource::POSITIVE);
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
  edit.setInletType(IonSource::DIRECT);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setIonizationMethod(IonSource::ESI);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setPolarity(IonSource::POSITIVE);
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



