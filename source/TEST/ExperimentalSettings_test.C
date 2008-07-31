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

CHECK((const String& getComment() const))
	ExperimentalSettings tmp;
	TEST_EQUAL(tmp.getComment(), "");
RESULT

CHECK((void setComment(const String& comment)))
	ExperimentalSettings tmp;
	tmp.setComment("bla");
	TEST_EQUAL(tmp.getComment(), "bla");
RESULT

CHECK(const std::vector<ContactPerson>& getContacts() const)
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getContacts().size(),0);
RESULT

CHECK(void setContacts(const std::vector<ContactPerson>& contacts))
  ExperimentalSettings tmp;
  std::vector<ContactPerson> dummy;
  ContactPerson c;
  c.setFirstName("bla17");
  c.setLastName("blubb17");
  dummy.push_back(c);
  c.setFirstName("bla18");
  c.setLastName("blubb18");
  dummy.push_back(c);
  tmp.setContacts(dummy);
  TEST_EQUAL(tmp.getContacts().size(),2);
  TEST_EQUAL(tmp.getContacts()[0].getFirstName(),"bla17");
  TEST_EQUAL(tmp.getContacts()[1].getFirstName(),"bla18");
  TEST_EQUAL(tmp.getContacts()[0].getLastName(),"blubb17");
  TEST_EQUAL(tmp.getContacts()[1].getLastName(),"blubb18");
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
  dummy.setName("bla3");
  tmp.setSample(dummy);
  TEST_EQUAL(tmp.getSample().getName(),"bla3");
RESULT

CHECK(void setSoftware(const Software& software))
  ExperimentalSettings tmp;
  Software dummy;
  dummy.setName("bla33");
  tmp.setSoftware(dummy);
  TEST_EQUAL(tmp.getSoftware().getName(),"bla33");
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
  tmp.getInstrument().setName("bla55");
  TEST_EQUAL(tmp.getInstrument().getName(),"bla55");
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
  c.setFirstName("bla17");
  tmp.getContacts().push_back(c);
  c.setFirstName("bla18");
  tmp.getContacts().push_back(c);
  TEST_EQUAL(tmp.getContacts().size(),2);
  TEST_EQUAL(tmp.getContacts()[0].getFirstName(),"bla17");
  TEST_EQUAL(tmp.getContacts()[1].getFirstName(),"bla18");
RESULT

CHECK(ExperimentalSettings(const ExperimentalSettings& source))
  ExperimentalSettings tmp;
  ProteinIdentification id;
  DateTime date;
	ProteinHit protein_hit;
	Real protein_significance_threshold = 63.2f;

	date.now();

	id.setDateTime(date);
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);
  
  tmp.getHPLC().setFlux(5);
  tmp.setComment("bla");
  tmp.getInstrument().setName("bla");
  tmp.getProcessingMethod().setDeisotoping(true);
  tmp.getSample().setName("bla2");
  tmp.getSoftware().setName("bla3");
  tmp.getSourceFile().setNameOfFile("bla4");
  tmp.getContacts().resize(1);
  tmp.addProteinIdentification(id);
  tmp.setMetaValue("label",String("label"));
  
  ExperimentalSettings tmp2(tmp);
  TEST_EQUAL(tmp2.getComment(),"bla");
  TEST_EQUAL(tmp2.getHPLC().getFlux(),5);
  TEST_EQUAL(tmp2.getInstrument().getName(),"bla");
  TEST_EQUAL(tmp2.getProcessingMethod().getDeisotoping(),true);
  TEST_EQUAL(tmp2.getSample().getName(),"bla2");
  TEST_EQUAL(tmp2.getSoftware().getName(),"bla3");
  TEST_EQUAL(tmp2.getSourceFile().getNameOfFile(),"bla4");
  TEST_EQUAL(tmp2.getContacts().size(),1);
  TEST_EQUAL(id == tmp2.getProteinIdentifications()[0], true);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");  
RESULT

CHECK(ExperimentalSettings& operator= (const ExperimentalSettings& source))
  ExperimentalSettings tmp;
  ProteinIdentification id;
  DateTime date;
	ProteinHit protein_hit;
	Real protein_significance_threshold = 63.2f;

	date.now();

	id.setDateTime(date);
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);

  tmp.getHPLC().setFlux(5);
  tmp.setComment("bla");
  tmp.getInstrument().setName("bla");
  tmp.getProcessingMethod().setDeisotoping(true);
  tmp.getSample().setName("bla2");
  tmp.getSoftware().setName("bla3");
  tmp.getSourceFile().setNameOfFile("bla4");
  tmp.getContacts().resize(1);
	tmp.addProteinIdentification(id);
  tmp.setMetaValue("label",String("label"));
  
  ExperimentalSettings tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getHPLC().getFlux(),5);
  TEST_EQUAL(tmp2.getInstrument().getName(),"bla");
  TEST_EQUAL(tmp2.getComment(),"bla");
  TEST_EQUAL(tmp2.getProcessingMethod().getDeisotoping(),true);
  TEST_EQUAL(tmp2.getSample().getName(),"bla2");
  TEST_EQUAL(tmp2.getSoftware().getName(),"bla3");
  TEST_EQUAL(tmp2.getSourceFile().getNameOfFile(),"bla4");
  TEST_EQUAL(tmp2.getContacts().size(),1);
  TEST_EQUAL(tmp2.getProteinIdentifications().size(), 1);
  TEST_EQUAL(id == tmp2.getProteinIdentifications()[0], true);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");  

  tmp2 = ExperimentalSettings();
  TEST_EQUAL(tmp2.getHPLC().getFlux(),0);
  TEST_EQUAL(tmp2.getInstrument().getName(),"");
  TEST_EQUAL(tmp2.getProcessingMethod().getDeisotoping(),false);
  TEST_EQUAL(tmp2.getSample().getName(),"");
  TEST_EQUAL(tmp2.getSoftware().getName(),"");
  TEST_EQUAL(tmp2.getSourceFile().getNameOfFile(),"");
  TEST_EQUAL(tmp2.getContacts().size(),0);
  TEST_EQUAL(tmp2.getProteinIdentifications().size(), 0);
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
RESULT

CHECK(bool operator== (const ExperimentalSettings& rhs) const)
  ExperimentalSettings edit, empty;
  ProteinIdentification id;
  DateTime date;
	ProteinHit protein_hit;
	Real protein_significance_threshold = 63.2f;

	date.now();

	id.setDateTime(date);
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);
  
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
	edit.addProteinIdentification(id);
  TEST_EQUAL(edit==empty, false);

  edit = empty;
	edit.setComment("bla");
  TEST_EQUAL(edit==empty, false);

	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit==empty,false);
RESULT

CHECK(bool operator!= (const ExperimentalSettings& rhs) const)
  ExperimentalSettings edit, empty;
  ProteinIdentification id;
  DateTime date;
	ProteinHit protein_hit;
	Real protein_significance_threshold = 63.2f;

	date.now();

	id.setDateTime(date);
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);
  
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
	edit.addProteinIdentification(id);
  TEST_EQUAL(edit!=empty, true);

	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty,true);
RESULT

CHECK(const std::vector<ProteinIdentification>& getProteinIdentifications() const)
  ExperimentalSettings settings;
  ProteinIdentification id;
  DateTime date;
	ProteinHit protein_hit;
	Real protein_significance_threshold = 63.2f;

	date.now();

	id.setDateTime(date);
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);
	
	settings.addProteinIdentification(id);
	const ProteinIdentification& test_id = settings.getProteinIdentifications()[0];
	TEST_EQUAL(id == test_id, true)
RESULT

CHECK(std::vector<ProteinIdentification>& getProteinIdentifications())
  ExperimentalSettings settings;
  ProteinIdentification id;
  DateTime date;
	ProteinHit protein_hit;
	Real protein_significance_threshold = 63.2f;

	date.now();

	id.setDateTime(date);
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);
	
	settings.addProteinIdentification(id);
	ProteinIdentification& test_id = settings.getProteinIdentifications()[0];
	TEST_EQUAL(id == test_id, true)
RESULT

CHECK(void addProteinIdentification(ProteinIdentification& protein_identification))
  ExperimentalSettings settings;
  ProteinIdentification id;
  DateTime date;
	ProteinHit protein_hit;
	Real protein_significance_threshold = 63.2f;

	date.now();

	id.setDateTime(date);
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);
	
	settings.addProteinIdentification(id);
	ProteinIdentification& test_id = settings.getProteinIdentifications()[0];
	TEST_EQUAL(id == test_id, true)
RESULT

CHECK(void setProteinIdentifications(const std::vector<ProteinIdentification>& protein_identifications))
  ExperimentalSettings settings;
  ProteinIdentification id;
  DateTime date;
	ProteinHit protein_hit;
	Real protein_significance_threshold = 63.2f;
	vector<ProteinIdentification> ids;

	date.now();

	id.setDateTime(date);
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);
	ids.push_back(id);
	id.setSignificanceThreshold(21.f);
	ids.push_back(id);
	settings.setProteinIdentifications(ids);
	TEST_EQUAL(ids == settings.getProteinIdentifications(), true)
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



