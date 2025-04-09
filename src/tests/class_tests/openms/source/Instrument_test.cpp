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
#include <OpenMS/METADATA/Instrument.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Instrument, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Instrument* ptr = nullptr;
Instrument* nullPointer = nullptr;
START_SECTION(Instrument())
	ptr = new Instrument();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~Instrument())
	delete ptr;
END_SECTION

START_SECTION(const std::vector<IonDetector>& getIonDetectors() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getIonDetectors().size(),0)
END_SECTION

START_SECTION(const std::vector<IonSource>& getIonSources() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getIonSources().size(),0)
END_SECTION

START_SECTION(const String& getCustomizations() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getCustomizations(),"");
END_SECTION

START_SECTION(const String& getModel() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getModel(),"");
END_SECTION

START_SECTION(const String& getName() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getName(),"");
END_SECTION

START_SECTION(const String& getVendor() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getVendor(),"");
END_SECTION

START_SECTION(const std::vector<MassAnalyzer>& getMassAnalyzers() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getMassAnalyzers().size(),0);
END_SECTION

START_SECTION(const Software& getSoftware() const)
	Instrument tmp;
  TEST_STRING_EQUAL(tmp.getSoftware().getName(),"");
END_SECTION

START_SECTION(IonOpticsType getIonOptics() const)
	Instrument tmp;
  TEST_EQUAL(tmp.getIonOptics(),Instrument::UNKNOWN);
END_SECTION

START_SECTION(void setIonOptics(IonOpticsType ion_optics))
	Instrument tmp;
	tmp.setIonOptics(Instrument::REFLECTRON);
  TEST_EQUAL(tmp.getIonOptics(),Instrument::REFLECTRON);
END_SECTION

START_SECTION(void setCustomizations(const String& customizations))
  Instrument tmp;
  tmp.setCustomizations("Customizations");
  TEST_EQUAL(tmp.getCustomizations(),"Customizations");
END_SECTION

START_SECTION(void setIonDetectors(const std::vector<IonDetector>& ion_detectors))
  Instrument tmp;
  std::vector<IonDetector> dummy;
  dummy.resize(1); 
  tmp.setIonDetectors(dummy);
  TEST_EQUAL(tmp.getIonDetectors().size(),1);
END_SECTION

START_SECTION(void setIonSources(const std::vector<IonSource>& ion_sources))
  Instrument tmp;
  std::vector<IonSource> dummy;
  dummy.resize(1);
  tmp.setIonSources(dummy);
  TEST_EQUAL(tmp.getIonSources().size(),1);
END_SECTION

START_SECTION(void setMassAnalyzers(const std::vector<MassAnalyzer>& mass_analyzers))
  Instrument tmp;
  MassAnalyzer dummy;
  dummy.setScanTime(47.11);
  vector<MassAnalyzer> dummy2;
  dummy2.push_back(dummy);
  dummy.setScanTime(47.12);
  dummy2.push_back(dummy);
  tmp.setMassAnalyzers(dummy2);
  TEST_EQUAL(tmp.getMassAnalyzers().size(),2);
  TEST_REAL_SIMILAR(tmp.getMassAnalyzers()[0].getScanTime(),47.11);
  TEST_REAL_SIMILAR(tmp.getMassAnalyzers()[1].getScanTime(),47.12);
END_SECTION

START_SECTION(void setModel(const String& model))
  Instrument tmp;
  tmp.setModel("Model");
  TEST_EQUAL(tmp.getModel(),"Model");
END_SECTION

START_SECTION(void setName(const String& name))
  Instrument tmp;
  tmp.setName("Name");
  TEST_EQUAL(tmp.getName(),"Name");
END_SECTION

START_SECTION(void setVendor(const String& vendor))
  Instrument tmp;
  tmp.setVendor("Vendor");
  TEST_EQUAL(tmp.getVendor(),"Vendor");
END_SECTION

START_SECTION(void setSoftware(const Software& software))
	Instrument tmp;
	Software s;
	s.setName("sn");
	tmp.setSoftware(s);
  TEST_STRING_EQUAL(tmp.getSoftware().getName(),"sn");
END_SECTION

START_SECTION(std::vector<IonDetector>& getIonDetectors())
  Instrument tmp;
  tmp.getIonDetectors().resize(1);
  TEST_EQUAL(tmp.getIonDetectors().size(),1);
END_SECTION

START_SECTION(std::vector<IonSource>& getIonSources())
  Instrument tmp;
  tmp.getIonSources().resize(1);
  TEST_EQUAL(tmp.getIonSources().size(),1);
END_SECTION

START_SECTION(std::vector<MassAnalyzer>& getMassAnalyzers())
  Instrument tmp;
  tmp.getMassAnalyzers().resize(2);
  tmp.getMassAnalyzers()[0].setScanTime(47.11);
  tmp.getMassAnalyzers()[1].setScanTime(47.12);
  TEST_EQUAL(tmp.getMassAnalyzers().size(),2);
  TEST_REAL_SIMILAR(tmp.getMassAnalyzers()[0].getScanTime(),47.11);
  TEST_REAL_SIMILAR(tmp.getMassAnalyzers()[1].getScanTime(),47.12);
END_SECTION

START_SECTION(Software& getSoftware())
	Instrument tmp;
	tmp.getSoftware().setName("sn");
  TEST_STRING_EQUAL(tmp.getSoftware().getName(),"sn");
END_SECTION

START_SECTION(Instrument(const Instrument& source))
  Instrument tmp;
  tmp.getMassAnalyzers().resize(1);
  tmp.getMassAnalyzers()[0].setScanTime(47.11);
	tmp.getIonSources().resize(1);
	tmp.getIonDetectors().resize(1);
  tmp.setModel("Model");
  tmp.setName("Name");
  tmp.setVendor("Vendor");
  tmp.setMetaValue("label",String("label"));
  tmp.getSoftware().setName("sn");
	tmp.setIonOptics(Instrument::REFLECTRON);
  
  Instrument tmp2(tmp);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  TEST_EQUAL(tmp2.getName(),"Name");
  TEST_EQUAL(tmp2.getModel(),"Model");
  TEST_EQUAL(tmp2.getVendor(),"Vendor");
  TEST_EQUAL(tmp2.getIonDetectors().size(),1);
  TEST_EQUAL(tmp2.getIonSources().size(),1);
  TEST_EQUAL(tmp2.getMassAnalyzers().size(),1);
  TEST_REAL_SIMILAR(tmp2.getMassAnalyzers()[0].getScanTime(),47.11);
  TEST_EQUAL(tmp2.getSoftware().getName(),"sn");
  TEST_EQUAL(tmp2.getIonOptics(),Instrument::REFLECTRON);
END_SECTION

START_SECTION(Instrument& operator= (const Instrument& source))
  Instrument tmp;
  tmp.getMassAnalyzers().resize(1);
  tmp.getMassAnalyzers()[0].setScanTime(47.11);
	tmp.getIonSources().resize(1);
	tmp.getIonDetectors().resize(1);
  tmp.setModel("Model");
  tmp.setName("Name");
  tmp.setVendor("Vendor");
  tmp.setMetaValue("label",String("label"));
  tmp.getSoftware().setName("sn");
	tmp.setIonOptics(Instrument::REFLECTRON);
		
  Instrument tmp2;
  tmp2 = tmp;
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  TEST_EQUAL(tmp2.getName(),"Name");
  TEST_EQUAL(tmp2.getModel(),"Model");
  TEST_EQUAL(tmp2.getVendor(),"Vendor");
  TEST_EQUAL(tmp2.getIonDetectors().size(),1);
  TEST_EQUAL(tmp2.getIonSources().size(),1);
  TEST_EQUAL(tmp2.getMassAnalyzers().size(),1);
  TEST_REAL_SIMILAR(tmp2.getMassAnalyzers()[0].getScanTime(),47.11);
  TEST_EQUAL(tmp2.getSoftware().getName(),"sn");
  TEST_EQUAL(tmp2.getIonOptics(),Instrument::REFLECTRON);

  tmp2 = Instrument();
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
  TEST_EQUAL(tmp2.getName(),"");
  TEST_EQUAL(tmp2.getModel(),"");
  TEST_EQUAL(tmp2.getVendor(),"");
  TEST_EQUAL(tmp2.getIonDetectors().size(),0);
  TEST_EQUAL(tmp2.getIonSources().size(),0);
  TEST_EQUAL(tmp2.getMassAnalyzers().size(),0);
  TEST_EQUAL(tmp2.getSoftware().getName(),"");
  TEST_EQUAL(tmp2.getIonOptics(),Instrument::UNKNOWN);
END_SECTION

START_SECTION(bool operator== (const Instrument& rhs) const)
  Instrument edit,empty;
  
  TEST_EQUAL(edit==empty,true);
  
  edit.getMassAnalyzers().resize(1);
  TEST_EQUAL(edit==empty,false);
	
	edit = empty;
	edit.getIonSources().resize(1);
	TEST_EQUAL(edit==empty,false);
	
	edit = empty;
	edit.getIonDetectors().resize(1);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setModel("Model");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setName("Name");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.getSoftware().setName("sn");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setVendor("Vendor");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
	edit.setIonOptics(Instrument::REFLECTRON);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit==empty,false);
END_SECTION

START_SECTION(bool operator!= (const Instrument& rhs) const)
  Instrument edit,empty;
  
  TEST_EQUAL(edit!=empty,false);
  
  edit.getMassAnalyzers().resize(1);
  TEST_EQUAL(edit!=empty,true);
	
	edit = empty;
	edit.getIonSources().resize(1);
	TEST_EQUAL(edit!=empty,true);
	
	edit = empty;
	edit.getIonDetectors().resize(1);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setModel("Model");
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setName("Name");
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setVendor("Vendor");
  TEST_EQUAL(edit!=empty,true);

  edit = empty;
  edit.getSoftware().setName("sn");
  TEST_EQUAL(edit!=empty,true);

  edit = empty;
	edit.setIonOptics(Instrument::REFLECTRON);
  TEST_EQUAL(edit!=empty,true)
  
  edit = empty;
  edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty,true);
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



