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

#include <unordered_set>
#include <unordered_map>

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
  TEST_EQUAL(names.size(), static_cast<size_t>(IonSource::InletType::SIZE_OF_INLETTYPE));
  TEST_EQUAL(names[static_cast<size_t>(IonSource::InletType::INLETNULL)], "Unknown");
  TEST_EQUAL(names[static_cast<size_t>(IonSource::InletType::DIRECT)], "Direct");
  TEST_EQUAL(names[static_cast<size_t>(IonSource::InletType::NANOSPRAY)], "Nanospray inlet");
END_SECTION

START_SECTION((static StringList getAllNamesOfIonizationMethod()))
  StringList names = IonSource::getAllNamesOfIonizationMethod();
  TEST_EQUAL(names.size(), static_cast<size_t>(IonSource::IonizationMethod::SIZE_OF_IONIZATIONMETHOD));
  TEST_EQUAL(names[static_cast<size_t>(IonSource::IonizationMethod::IONMETHODNULL)], "Unknown");
  TEST_EQUAL(names[static_cast<size_t>(IonSource::IonizationMethod::ESI)], "Electrospray ionisation");
  TEST_EQUAL(names[static_cast<size_t>(IonSource::IonizationMethod::MALDI)], "Matrix-assisted laser desorption ionization");
END_SECTION

START_SECTION((static StringList getAllNamesOfPolarity()))
  StringList names = IonSource::getAllNamesOfPolarity();
  TEST_EQUAL(names.size(), static_cast<size_t>(IonSource::Polarity::SIZE_OF_POLARITY));
  TEST_EQUAL(names[static_cast<size_t>(IonSource::Polarity::POLNULL)], "unknown");
  TEST_EQUAL(names[static_cast<size_t>(IonSource::Polarity::POSITIVE)], "positive");
  TEST_EQUAL(names[static_cast<size_t>(IonSource::Polarity::NEGATIVE)], "negative");
END_SECTION

START_SECTION(([EXTRA] std::hash<IonSource>))
{
  // Test that equal objects have equal hashes
  IonSource is1, is2;
  is1.setInletType(IonSource::InletType::DIRECT);
  is1.setIonizationMethod(IonSource::IonizationMethod::ESI);
  is1.setPolarity(IonSource::Polarity::POSITIVE);
  is1.setOrder(45);
  is1.setMetaValue("label", String("test"));

  is2.setInletType(IonSource::InletType::DIRECT);
  is2.setIonizationMethod(IonSource::IonizationMethod::ESI);
  is2.setPolarity(IonSource::Polarity::POSITIVE);
  is2.setOrder(45);
  is2.setMetaValue("label", String("test"));

  TEST_EQUAL(is1 == is2, true)
  TEST_EQUAL(std::hash<IonSource>{}(is1), std::hash<IonSource>{}(is2))

  // Test that different objects (likely) have different hashes
  IonSource is3;
  is3.setInletType(IonSource::InletType::BATCH);
  is3.setIonizationMethod(IonSource::IonizationMethod::MALDI);
  is3.setPolarity(IonSource::Polarity::NEGATIVE);
  is3.setOrder(10);

  TEST_EQUAL(is1 == is3, false)
  // Hash values should differ (not guaranteed but highly likely for different inputs)
  TEST_NOT_EQUAL(std::hash<IonSource>{}(is1), std::hash<IonSource>{}(is3))

  // Test use in unordered_set
  std::unordered_set<IonSource> ion_source_set;
  ion_source_set.insert(is1);
  ion_source_set.insert(is2);  // Duplicate, should not increase size
  ion_source_set.insert(is3);
  TEST_EQUAL(ion_source_set.size(), 2)

  // Test use in unordered_map
  std::unordered_map<IonSource, int> ion_source_map;
  ion_source_map[is1] = 1;
  ion_source_map[is2] = 2;  // Should overwrite is1's value
  ion_source_map[is3] = 3;
  TEST_EQUAL(ion_source_map.size(), 2)
  TEST_EQUAL(ion_source_map[is1], 2)
  TEST_EQUAL(ion_source_map[is3], 3)

  // Test that default-constructed objects have equal hashes
  IonSource default1, default2;
  TEST_EQUAL(default1 == default2, true)
  TEST_EQUAL(std::hash<IonSource>{}(default1), std::hash<IonSource>{}(default2))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



