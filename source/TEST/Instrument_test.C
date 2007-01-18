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

CHECK(const IonDetector& getIonDetector() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getIonDetector()==IonDetector(),true);
RESULT

CHECK(const IonSource& getIonSource() const)
  Instrument tmp;
  TEST_EQUAL(tmp.getIonSource()==IonSource(),true);
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

CHECK(void setCustomizations(const String& customizations))
  Instrument tmp;
  tmp.setCustomizations("Customizations");
  TEST_EQUAL(tmp.getCustomizations(),"Customizations");
RESULT

CHECK(void setIonDetector(const IonDetector& ion_detector))
  Instrument tmp;
  IonDetector dummy;
  dummy.setResolution(47.11); 
  tmp.setIonDetector(dummy);
  TEST_REAL_EQUAL(tmp.getIonDetector().getResolution(),47.11);
RESULT

CHECK(void setIonSource(const IonSource& ion_source))
  Instrument tmp;
  IonSource dummy;
  dummy.setPolarity(IonSource::POSITIVE); 
  tmp.setIonSource(dummy);
  TEST_EQUAL(tmp.getIonSource().getPolarity(),IonSource::POSITIVE);
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

CHECK(IonDetector& getIonDetector())
  Instrument tmp;
  tmp.getIonDetector().setResolution(47.11); 
  TEST_REAL_EQUAL(tmp.getIonDetector().getResolution(),47.11);
RESULT

CHECK(IonSource& getIonSource())
  Instrument tmp;
  tmp.getIonSource().setPolarity(IonSource::POSITIVE); 
  TEST_EQUAL(tmp.getIonSource().getPolarity(),IonSource::POSITIVE);
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

CHECK(Instrument(const Instrument& source))
  Instrument tmp;
  tmp.getMassAnalyzers().resize(1);
  tmp.getMassAnalyzers()[0].setScanTime(47.11);
	tmp.getIonSource().setPolarity(IonSource::POSITIVE);
	tmp.getIonDetector().setResolution(47.12); 
  tmp.setModel("Model");
  tmp.setName("Name");
  tmp.setVendor("Vendor");
  tmp.setMetaValue("label",String("label"));
  
  Instrument tmp2(tmp);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  TEST_EQUAL(tmp2.getName(),"Name");
  TEST_EQUAL(tmp2.getModel(),"Model");
  TEST_EQUAL(tmp2.getVendor(),"Vendor");
  TEST_REAL_EQUAL(tmp2.getIonDetector().getResolution(),47.12);
  TEST_EQUAL(tmp2.getIonSource().getPolarity(),IonSource::POSITIVE);
  TEST_EQUAL(tmp2.getMassAnalyzers().size(),1);
  TEST_REAL_EQUAL(tmp2.getMassAnalyzers()[0].getScanTime(),47.11);
RESULT

CHECK(Instrument& operator= (const Instrument& source))
  Instrument tmp;
  tmp.getMassAnalyzers().resize(1);
  tmp.getMassAnalyzers()[0].setScanTime(47.11);
	tmp.getIonSource().setPolarity(IonSource::POSITIVE);
	tmp.getIonDetector().setResolution(47.12); 
  tmp.setModel("Model");
  tmp.setName("Name");
  tmp.setVendor("Vendor");
  tmp.setMetaValue("label",String("label"));
  
  Instrument tmp2;
  tmp2 = tmp;
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  TEST_EQUAL(tmp2.getName(),"Name");
  TEST_EQUAL(tmp2.getModel(),"Model");
  TEST_EQUAL(tmp2.getVendor(),"Vendor");
  TEST_REAL_EQUAL(tmp2.getIonDetector().getResolution(),47.12);
  TEST_EQUAL(tmp2.getIonSource().getPolarity(),IonSource::POSITIVE);
  TEST_EQUAL(tmp2.getMassAnalyzers().size(),1);
  TEST_REAL_EQUAL(tmp2.getMassAnalyzers()[0].getScanTime(),47.11);

  tmp2 = Instrument();
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
  TEST_EQUAL(tmp2.getName(),"");
  TEST_EQUAL(tmp2.getModel(),"");
  TEST_EQUAL(tmp2.getVendor(),"");
  TEST_REAL_EQUAL(tmp2.getIonDetector().getResolution(),0.0);
  TEST_EQUAL(tmp2.getIonSource().getPolarity(),IonSource::POLNULL);
  TEST_EQUAL(tmp2.getMassAnalyzers().size(),0);
RESULT

CHECK(bool operator== (const Instrument& rhs) const)
  Instrument edit,empty;
  
  TEST_EQUAL(edit==empty,true);
  
  edit.getMassAnalyzers().resize(1);
  TEST_EQUAL(edit==empty,false);
	
	edit = empty;
	edit.getIonSource().setPolarity(IonSource::POSITIVE);
	TEST_EQUAL(edit==empty,false);
	
	edit = empty;
	edit.getIonDetector().setResolution(47.12); 
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setModel("Model");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setName("Name");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setVendor("Vendor");
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
	edit.getIonSource().setPolarity(IonSource::POSITIVE);
	TEST_EQUAL(edit!=empty,true);
	
	edit = empty;
	edit.getIonDetector().setResolution(47.12); 
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
  edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty,true);
RESULT



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



