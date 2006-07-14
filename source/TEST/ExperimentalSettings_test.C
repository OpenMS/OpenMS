// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/METADATA/ExperimentalSettings.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ExperimentalSettings, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ExperimentalSettings* ptr = 0;
CHECK(ExperimentalSettings())
	ptr = new ExperimentalSettings();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~ExperimentalSettings())
	delete ptr;
RESULT

CHECK(ExperimentType getType() const)
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getType(),ExperimentalSettings::UNKNOWN);
RESULT

CHECK(const Date& getDate() const)
  ExperimentalSettings tmp;
  String s;
  tmp.getDate().get(s);
  TEST_EQUAL(s,"0000-00-00");
RESULT

CHECK(const HPLC& getHPLC() const)
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getHPLC()==HPLC(),true);
RESULT

CHECK(const Instrument& getInstrument() const)
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getInstrument()==Instrument(),true);
RESULT

CHECK(const ProcessingMethod& getProcessingMethod() const)
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getProcessingMethod()==ProcessingMethod(),true);
RESULT

CHECK(const Sample& getSample() const)
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getSample()==Sample(),true);
RESULT

CHECK(const Software& getSoftware() const)
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getSoftware()==Software(),true);
RESULT

CHECK(const SourceFile& getSourceFile() const)
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getSourceFile()==SourceFile(),true);
RESULT

CHECK(const std::vector<ContactPerson>& getContacts() const)
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getContacts().size(),0);
RESULT

CHECK(void setContacts(const std::vector<ContactPerson>& contacts))
  ExperimentalSettings tmp;
  std::vector<ContactPerson> dummy;
  ContactPerson c;
  c.setName("bla17");
  dummy.push_back(c);
  c.setName("bla18");
  dummy.push_back(c);
  tmp.setContacts(dummy);
  TEST_EQUAL(tmp.getContacts().size(),2);
  TEST_EQUAL(tmp.getContacts()[0].getName(),"bla17");
  TEST_EQUAL(tmp.getContacts()[1].getName(),"bla18");
RESULT

CHECK(void setDate(const Date& date))
  ExperimentalSettings tmp;
  Date dummy;
  String s;
  dummy.set("02/07/2006");
  tmp.setDate(dummy);
  tmp.getDate().get(s);
  TEST_EQUAL(s,"2006-02-07");
RESULT

CHECK(void setHPLC(const HPLC& hplc))
  ExperimentalSettings tmp;
  HPLC dummy;
  dummy.setFlux(5);
  tmp.setHPLC(dummy);
  TEST_EQUAL(tmp.getHPLC().getFlux(),5);
RESULT

CHECK(void setInstrument(const Instrument& instrument))
  ExperimentalSettings tmp;
  Instrument dummy;
  dummy.setName("bla");
  tmp.setInstrument(dummy);
  TEST_EQUAL(tmp.getInstrument().getName(),"bla");
RESULT

CHECK(void setProcessingMethod(const ProcessingMethod& processing_method))
  ExperimentalSettings tmp;
  ProcessingMethod dummy;
  dummy.setDeisotoping(true);
  tmp.setProcessingMethod(dummy);
  TEST_EQUAL(tmp.getProcessingMethod().getDeisotoping(),true);
RESULT

CHECK(void setSample(const Sample& sample))
  ExperimentalSettings tmp;
  Sample dummy;
  dummy.setName("bla2");
  tmp.setSample(dummy);
  TEST_EQUAL(tmp.getSample().getName(),"bla2");
RESULT

CHECK(void setSoftware(const Software& software))
  ExperimentalSettings tmp;
  Software dummy;
  dummy.setName("bla3");
  tmp.setSoftware(dummy);
  TEST_EQUAL(tmp.getSoftware().getName(),"bla3");
RESULT

CHECK(void setSourceFile(const SourceFile& source_file))
  ExperimentalSettings tmp;
  SourceFile dummy;
  dummy.setNameOfFile("bla4");
  tmp.setSourceFile(dummy);
  TEST_EQUAL(tmp.getSourceFile().getNameOfFile(),"bla4");
RESULT

CHECK(void setType(ExperimentType type))
  ExperimentalSettings tmp;
  tmp.setType(ExperimentalSettings::HPLC_MS);
  TEST_EQUAL(tmp.getType(),ExperimentalSettings::HPLC_MS);
RESULT

CHECK(HPLC& getHPLC())
  ExperimentalSettings tmp;
  tmp.getHPLC().setFlux(5);
  TEST_EQUAL(tmp.getHPLC().getFlux(),5);
RESULT

CHECK(Instrument& getInstrument())
  ExperimentalSettings tmp;
  tmp.getInstrument().setName("bla");
  TEST_EQUAL(tmp.getInstrument().getName(),"bla");
RESULT

CHECK(ProcessingMethod& getProcessingMethod())
  ExperimentalSettings tmp;
  tmp.getProcessingMethod().setDeisotoping(true);
  TEST_EQUAL(tmp.getProcessingMethod().getDeisotoping(),true);
RESULT

CHECK(Sample& getSample())
  ExperimentalSettings tmp;
  tmp.getSample().setName("bla2");
  TEST_EQUAL(tmp.getSample().getName(),"bla2");
RESULT

CHECK(Software& getSoftware())
  ExperimentalSettings tmp;
  tmp.getSoftware().setName("bla3");
  TEST_EQUAL(tmp.getSoftware().getName(),"bla3");
RESULT

CHECK(SourceFile& getSourceFile())
  ExperimentalSettings tmp;
  tmp.getSourceFile().setNameOfFile("bla4");
  TEST_EQUAL(tmp.getSourceFile().getNameOfFile(),"bla4");
RESULT

CHECK(std::vector<ContactPerson>& getContacts())
  ExperimentalSettings tmp;
  ContactPerson c;
  c.setName("bla17");
  tmp.getContacts().push_back(c);
  c.setName("bla18");
  tmp.getContacts().push_back(c);
  TEST_EQUAL(tmp.getContacts().size(),2);
  TEST_EQUAL(tmp.getContacts()[0].getName(),"bla17");
  TEST_EQUAL(tmp.getContacts()[1].getName(),"bla18");
RESULT

CHECK(ExperimentalSettings(const ExperimentalSettings& source))
  ExperimentalSettings tmp;
  tmp.getHPLC().setFlux(5);
  tmp.getInstrument().setName("bla");
  tmp.getProcessingMethod().setDeisotoping(true);
  tmp.getSample().setName("bla2");
  tmp.getSoftware().setName("bla3");
  tmp.getSourceFile().setNameOfFile("bla4");
  tmp.getContacts().resize(1);
  tmp.setMetaValue("label",String("label"));
  
  ExperimentalSettings tmp2(tmp);
  TEST_EQUAL(tmp2.getHPLC().getFlux(),5);
  TEST_EQUAL(tmp2.getInstrument().getName(),"bla");
  TEST_EQUAL(tmp2.getProcessingMethod().getDeisotoping(),true);
  TEST_EQUAL(tmp2.getSample().getName(),"bla2");
  TEST_EQUAL(tmp2.getSoftware().getName(),"bla3");
  TEST_EQUAL(tmp2.getSourceFile().getNameOfFile(),"bla4");
  TEST_EQUAL(tmp2.getContacts().size(),1);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");  
RESULT

CHECK(ExperimentalSettings& operator= (const ExperimentalSettings& source))
  ExperimentalSettings tmp;
  tmp.getHPLC().setFlux(5);
  tmp.getInstrument().setName("bla");
  tmp.getProcessingMethod().setDeisotoping(true);
  tmp.getSample().setName("bla2");
  tmp.getSoftware().setName("bla3");
  tmp.getSourceFile().setNameOfFile("bla4");
  tmp.getContacts().resize(1);
  tmp.setMetaValue("label",String("label"));
  
  ExperimentalSettings tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getHPLC().getFlux(),5);
  TEST_EQUAL(tmp2.getInstrument().getName(),"bla");
  TEST_EQUAL(tmp2.getProcessingMethod().getDeisotoping(),true);
  TEST_EQUAL(tmp2.getSample().getName(),"bla2");
  TEST_EQUAL(tmp2.getSoftware().getName(),"bla3");
  TEST_EQUAL(tmp2.getSourceFile().getNameOfFile(),"bla4");
  TEST_EQUAL(tmp2.getContacts().size(),1);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");  

  tmp2 = ExperimentalSettings();
  TEST_EQUAL(tmp2.getHPLC().getFlux(),0);
  TEST_EQUAL(tmp2.getInstrument().getName(),"");
  TEST_EQUAL(tmp2.getProcessingMethod().getDeisotoping(),false);
  TEST_EQUAL(tmp2.getSample().getName(),"");
  TEST_EQUAL(tmp2.getSoftware().getName(),"");
  TEST_EQUAL(tmp2.getSourceFile().getNameOfFile(),"");
  TEST_EQUAL(tmp2.getContacts().size(),0);
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
RESULT

CHECK(bool operator== (const ExperimentalSettings& rhs) const)
  ExperimentalSettings edit, empty;
  
  TEST_EQUAL(edit==empty,true);
  
  edit.getHPLC().setFlux(5);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.getInstrument().setName("bla");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.getProcessingMethod().setDeisotoping(true);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.getSample().setName("bla2");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.getSoftware().setName("bla3");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.getSourceFile().setNameOfFile("bla4");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.getContacts().resize(1);
  TEST_EQUAL(edit==empty,false);

	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit==empty,false);
RESULT

CHECK(bool operator!= (const ExperimentalSettings& rhs) const)
  ExperimentalSettings edit, empty;
  
  TEST_EQUAL(edit!=empty,false);
  
  edit.getHPLC().setFlux(5);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.getInstrument().setName("bla");
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.getProcessingMethod().setDeisotoping(true);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.getSample().setName("bla2");
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.getSoftware().setName("bla3");
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.getSourceFile().setNameOfFile("bla4");
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.getContacts().resize(1);
  TEST_EQUAL(edit!=empty,true);

	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty,true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



