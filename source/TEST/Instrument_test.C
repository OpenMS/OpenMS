// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/METADATA/Instrument.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Instrument, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Instrument* ptr = 0;
CHECK(Instrument())
	ptr = new Instrument();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~Instrument())
	delete ptr;
RESULT

CHECK(const std::vector<IonDetector>& getIonDetectors() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getIonDetectors().size(),0)
RESULT

CHECK(const std::vector<IonSource>& getIonSources() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getIonSources().size(),0)
RESULT

CHECK(const String& getCustomizations() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getCustomizations(),"");
RESULT

CHECK(const String& getModel() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getModel(),"");
RESULT

CHECK(const String& getName() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getName(),"");
RESULT

CHECK(const String& getVendor() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getVendor(),"");
RESULT

CHECK(const std::vector<MassAnalyzer>& getMassAnalyzers() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getMassAnalyzers().size(),0);
RESULT

CHECK(const Software& getSoftware() const)
	Instrument tmp;
  TEST_STRING_EQUAL(tmp.getSoftware().getName(),"");
RESULT

CHECK(IonOpticsType getIonOptics() const)
	Instrument tmp;
  TEST_EQUAL(tmp.getIonOptics(),Instrument::UNKNOWN);
RESULT

CHECK(void setIonOptics(IonOpticsType ion_optics))
	Instrument tmp;
	tmp.setIonOptics(Instrument::REFLECTRON);
  TEST_EQUAL(tmp.getIonOptics(),Instrument::REFLECTRON);
RESULT

CHECK(void setCustomizations(const String& customizations))
  Instrument tmp;
  tmp.setCustomizations("Customizations");
  TEST_EQUAL(tmp.getCustomizations(),"Customizations");
RESULT

CHECK(void setIonDetectors(const std::vector<IonDetector>& ion_detectors))
  Instrument tmp;
  std::vector<IonDetector> dummy;
  dummy.resize(1); 
  tmp.setIonDetectors(dummy);
  TEST_EQUAL(tmp.getIonDetectors().size(),1);
RESULT

CHECK(void setIonSources(const std::vector<IonSource>& ion_sources))
  Instrument tmp;
  std::vector<IonSource> dummy;
  dummy.resize(1);
  tmp.setIonSources(dummy);
  TEST_EQUAL(tmp.getIonSources().size(),1);
RESULT

CHECK(void setMassAnalyzers(const std::vector<MassAnalyzer>& mass_analyzers))
  Instrument tmp;
  MassAnalyzer dummy;
  dummy.setScanTime(47.11);
  vector<MassAnalyzer> dummy2;
  dummy2.push_back(dummy);
  dummy.setScanTime(47.12);
  dummy2.push_back(dummy);
  tmp.setMassAnalyzers(dummy2);
  TEST_EQUAL(tmp.getMassAnalyzers().size(),2);
  TEST_REAL_EQUAL(tmp.getMassAnalyzers()[0].getScanTime(),47.11);
  TEST_REAL_EQUAL(tmp.getMassAnalyzers()[1].getScanTime(),47.12);
RESULT

CHECK(void setModel(const String& model))
  Instrument tmp;
  tmp.setModel("Model");
  TEST_EQUAL(tmp.getModel(),"Model");
RESULT

CHECK(void setName(const String& name))
  Instrument tmp;
  tmp.setName("Name");
  TEST_EQUAL(tmp.getName(),"Name");
RESULT

CHECK(void setVendor(const String& vendor))
  Instrument tmp;
  tmp.setVendor("Vendor");
  TEST_EQUAL(tmp.getVendor(),"Vendor");
RESULT

CHECK(void setSoftware(const Software& software))
	Instrument tmp;
	Software s;
	s.setName("sn");
	tmp.setSoftware(s);
  TEST_STRING_EQUAL(tmp.getSoftware().getName(),"sn");
RESULT

CHECK(std::vector<IonDetector>& getIonDetectors())
  Instrument tmp;
  tmp.getIonDetectors().resize(1);
  TEST_REAL_EQUAL(tmp.getIonDetectors().size(),1);
RESULT

CHECK(std::vector<IonSource>& getIonSources())
  Instrument tmp;
  tmp.getIonSources().resize(1);
  TEST_EQUAL(tmp.getIonSources().size(),1);
RESULT

CHECK(std::vector<MassAnalyzer>& getMassAnalyzers())
  Instrument tmp;
  tmp.getMassAnalyzers().resize(2);
  tmp.getMassAnalyzers()[0].setScanTime(47.11);
  tmp.getMassAnalyzers()[1].setScanTime(47.12);
  TEST_EQUAL(tmp.getMassAnalyzers().size(),2);
  TEST_REAL_EQUAL(tmp.getMassAnalyzers()[0].getScanTime(),47.11);
  TEST_REAL_EQUAL(tmp.getMassAnalyzers()[1].getScanTime(),47.12);
RESULT

CHECK(Software& getSoftware())
	Instrument tmp;
	tmp.getSoftware().setName("sn");
  TEST_STRING_EQUAL(tmp.getSoftware().getName(),"sn");
RESULT

CHECK(Instrument(const Instrument& source))
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
  TEST_REAL_EQUAL(tmp2.getMassAnalyzers()[0].getScanTime(),47.11);
  TEST_EQUAL(tmp2.getSoftware().getName(),"sn");
  TEST_EQUAL(tmp2.getIonOptics(),Instrument::REFLECTRON);
RESULT

CHECK(Instrument& operator= (const Instrument& source))
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
  TEST_REAL_EQUAL(tmp2.getMassAnalyzers()[0].getScanTime(),47.11);
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
RESULT

CHECK(bool operator== (const Instrument& rhs) const)
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
RESULT

CHECK(bool operator!= (const Instrument& rhs) const)
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
RESULT



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



