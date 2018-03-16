// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/ExperimentalSettings.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ExperimentalSettings, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ExperimentalSettings* ptr = nullptr;
ExperimentalSettings* nullPointer = nullptr;
START_SECTION((ExperimentalSettings()))
	ptr = new ExperimentalSettings();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~ExperimentalSettings()))
	delete ptr;
END_SECTION

START_SECTION((const DateTime& getDateTime() const))
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getDateTime().get(),"0000-00-00 00:00:00");
END_SECTION

START_SECTION((const HPLC& getHPLC() const))
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getHPLC()==HPLC(),true);
END_SECTION

START_SECTION((const Instrument& getInstrument() const))
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getInstrument()==Instrument(),true);
END_SECTION

START_SECTION((const Sample& getSample() const))
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getSample()==Sample(),true);
END_SECTION

START_SECTION((const std::vector<SourceFile>& getSourceFiles() const))
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getSourceFiles().size(),0);
END_SECTION

START_SECTION((const String& getComment() const))
	ExperimentalSettings tmp;
	TEST_EQUAL(tmp.getComment(), "");
END_SECTION

START_SECTION((void setComment(const String& comment)))
	ExperimentalSettings tmp;
	tmp.setComment("bla");
	TEST_EQUAL(tmp.getComment(), "bla");
END_SECTION

START_SECTION((const String& getFractionIdentifier() const))
	ExperimentalSettings tmp;
	TEST_EQUAL(tmp.getFractionIdentifier(), "");
END_SECTION

START_SECTION((void setFractionIdentifier(const String& fraction_identifier)))
	ExperimentalSettings tmp;
	tmp.setFractionIdentifier("bla");
	TEST_EQUAL(tmp.getFractionIdentifier(), "bla");
END_SECTION


START_SECTION((void setContacts(const std::vector<ContactPerson>& contacts)))
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
END_SECTION

START_SECTION((void setDateTime(const DateTime &date)))
  ExperimentalSettings tmp;
  DateTime dummy;
  dummy.set("02/07/2006 01:02:03");
  tmp.setDateTime(dummy);
  TEST_EQUAL(tmp.getDateTime().get(),"2006-02-07 01:02:03");
END_SECTION

START_SECTION((void setHPLC(const HPLC& hplc)))
  ExperimentalSettings tmp;
  HPLC dummy;
  dummy.setFlux(5);
  tmp.setHPLC(dummy);
  TEST_EQUAL(tmp.getHPLC().getFlux(),5);
END_SECTION

START_SECTION((void setInstrument(const Instrument& instrument)))
  ExperimentalSettings tmp;
  Instrument dummy;
  dummy.setName("bla");
  tmp.setInstrument(dummy);
  TEST_EQUAL(tmp.getInstrument().getName(),"bla");
END_SECTION

START_SECTION((void setSample(const Sample& sample)))
  ExperimentalSettings tmp;
  Sample dummy;
  dummy.setName("bla3");
  tmp.setSample(dummy);
  TEST_EQUAL(tmp.getSample().getName(),"bla3");
END_SECTION

START_SECTION((void setSourceFiles(const std::vector< SourceFile > &source_files)))
  ExperimentalSettings tmp;
  vector<SourceFile> dummy;
  dummy.resize(1);
  tmp.setSourceFiles(dummy);
  TEST_EQUAL(tmp.getSourceFiles().size(),1);
END_SECTION

START_SECTION((HPLC& getHPLC()))
  ExperimentalSettings tmp;
  tmp.getHPLC().setFlux(5);
  TEST_EQUAL(tmp.getHPLC().getFlux(),5);
END_SECTION

START_SECTION((Instrument& getInstrument()))
  ExperimentalSettings tmp;
  tmp.getInstrument().setName("bla55");
  TEST_EQUAL(tmp.getInstrument().getName(),"bla55");
END_SECTION

START_SECTION((Sample& getSample()))
  ExperimentalSettings tmp;
  tmp.getSample().setName("bla2");
  TEST_EQUAL(tmp.getSample().getName(),"bla2");
END_SECTION

START_SECTION((std::vector<SourceFile>& getSourceFiles()))
  ExperimentalSettings tmp;
  tmp.getSourceFiles().resize(1);
  TEST_EQUAL(tmp.getSourceFiles().size(),1)
END_SECTION

START_SECTION((const std::vector<ContactPerson>& getContacts() const))
  ExperimentalSettings tmp;
  TEST_EQUAL(tmp.getContacts().size(),0);
END_SECTION

START_SECTION((std::vector<ContactPerson>& getContacts()))
  ExperimentalSettings tmp;
  ContactPerson c;
  c.setFirstName("bla17");
  tmp.getContacts().push_back(c);
  c.setFirstName("bla18");
  tmp.getContacts().push_back(c);
  TEST_EQUAL(tmp.getContacts().size(),2);
  TEST_EQUAL(tmp.getContacts()[0].getFirstName(),"bla17");
  TEST_EQUAL(tmp.getContacts()[1].getFirstName(),"bla18");
END_SECTION

START_SECTION((ExperimentalSettings(const ExperimentalSettings& source)))
  ExperimentalSettings tmp;
  ProteinIdentification id;
	ProteinHit protein_hit;
	float protein_significance_threshold = 63.2f;

	id.setDateTime(DateTime::now());
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);
  
  tmp.getHPLC().setFlux(5);
  tmp.setComment("bla");
  tmp.setFractionIdentifier("bla2");
  tmp.setIdentifier("lsid");
  tmp.getInstrument().setName("bla");
  tmp.getSample().setName("bla2");
  tmp.getSourceFiles().resize(1);
  tmp.getContacts().resize(1);
  tmp.getProteinIdentifications().push_back(id);
  tmp.setMetaValue("label",String("label"));
  
  ExperimentalSettings tmp2(tmp);
  TEST_EQUAL(tmp2.getComment(),"bla");
  TEST_EQUAL(tmp2.getFractionIdentifier(),"bla2");
  TEST_EQUAL(tmp2.getIdentifier(),"lsid");
  TEST_EQUAL(tmp2.getHPLC().getFlux(),5);
  TEST_EQUAL(tmp2.getInstrument().getName(),"bla");
  TEST_EQUAL(tmp2.getSample().getName(),"bla2");
  TEST_EQUAL(tmp2.getSourceFiles().size(),1);
  TEST_EQUAL(tmp2.getContacts().size(),1);
  TEST_EQUAL(id == tmp2.getProteinIdentifications()[0], true);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");  
END_SECTION

START_SECTION((ExperimentalSettings& operator= (const ExperimentalSettings& source)))
  ExperimentalSettings tmp;
  ProteinIdentification id;
	ProteinHit protein_hit;
	float protein_significance_threshold = 63.2f;

	id.setDateTime(DateTime::now());
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);

  tmp.getHPLC().setFlux(5);
  tmp.setComment("bla");
  tmp.setFractionIdentifier("bla2");
  tmp.setIdentifier("lsid");
  tmp.getInstrument().setName("bla");
  tmp.getSample().setName("bla2");
  tmp.getSourceFiles().resize(1);
  tmp.getContacts().resize(1);
	tmp.getProteinIdentifications().push_back(id);
  tmp.setMetaValue("label",String("label"));
  
  ExperimentalSettings tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getHPLC().getFlux(),5);
  TEST_EQUAL(tmp2.getInstrument().getName(),"bla");
  TEST_EQUAL(tmp2.getComment(),"bla");
  TEST_EQUAL(tmp2.getFractionIdentifier(),"bla2");
  TEST_EQUAL(tmp2.getIdentifier(),"lsid");
  TEST_EQUAL(tmp2.getSample().getName(),"bla2");
  TEST_EQUAL(tmp2.getSourceFiles().size(),1);
  TEST_EQUAL(tmp2.getContacts().size(),1);
  TEST_EQUAL(tmp2.getProteinIdentifications().size(), 1);
  TEST_EQUAL(id == tmp2.getProteinIdentifications()[0], true);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");  

  tmp2 = ExperimentalSettings();
  TEST_EQUAL(tmp2.getHPLC().getFlux(),0);
  TEST_EQUAL(tmp2.getInstrument().getName(),"");
  TEST_EQUAL(tmp2.getComment(),"");
  TEST_EQUAL(tmp2.getFractionIdentifier(),"");
  TEST_EQUAL(tmp2.getIdentifier(),"");
  TEST_EQUAL(tmp2.getSample().getName(),"");
  TEST_EQUAL(tmp2.getSourceFiles().size(),0);
  TEST_EQUAL(tmp2.getContacts().size(),0);
  TEST_EQUAL(tmp2.getProteinIdentifications().size(), 0);
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
END_SECTION

START_SECTION((bool operator== (const ExperimentalSettings& rhs) const))
  ExperimentalSettings edit, empty;
  ProteinIdentification id;
	ProteinHit protein_hit;
	float protein_significance_threshold = 63.2f;

	id.setDateTime(DateTime::now());
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);
  
  TEST_EQUAL(edit==empty,true);
  
  edit.getHPLC().setFlux(5);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.getInstrument().setName("bla");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.getSample().setName("bla2");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.getSourceFiles().resize(1);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.getContacts().resize(1);
  TEST_EQUAL(edit==empty,false);

  edit = empty;
	edit.getProteinIdentifications().push_back(id);
  TEST_EQUAL(edit==empty, false);

  edit = empty;
	edit.setComment("bla");
  TEST_EQUAL(edit==empty, false);

  edit = empty;
	edit.setFractionIdentifier("bla");
  TEST_EQUAL(edit==empty, false);

  edit = empty;
	edit.setIdentifier("bla");
  TEST_EQUAL(edit==empty, false);

	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit==empty,false);
END_SECTION

START_SECTION((bool operator!= (const ExperimentalSettings& rhs) const))
  ExperimentalSettings edit, empty;
  ProteinIdentification id;
	ProteinHit protein_hit;
	float protein_significance_threshold = 63.2f;

	id.setDateTime(DateTime::now());
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);
  
  TEST_EQUAL(edit!=empty,false);
  
  edit.getHPLC().setFlux(5);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.getInstrument().setName("bla");
  TEST_EQUAL(edit!=empty,true);

  edit = empty;
	edit.setComment("bla");
  TEST_EQUAL(edit!=empty, true);

  edit = empty;
	edit.setFractionIdentifier("bla2");
  TEST_EQUAL(edit!=empty, true);

  edit = empty;
  edit.getSample().setName("bla2");
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.getSourceFiles().resize(1);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.getContacts().resize(1);
  TEST_EQUAL(edit!=empty,true);

  edit = empty;
	edit.getProteinIdentifications().push_back(id);
  TEST_EQUAL(edit!=empty, true);

	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty,true);
END_SECTION

START_SECTION((const std::vector<ProteinIdentification>& getProteinIdentifications() const))
  ExperimentalSettings settings;
  ProteinIdentification id;
	ProteinHit protein_hit;
	float protein_significance_threshold = 63.2f;

	id.setDateTime(DateTime::now());
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);
	
	settings.getProteinIdentifications().push_back(id);
	const ProteinIdentification& test_id = settings.getProteinIdentifications()[0];
	TEST_EQUAL(id == test_id, true)
END_SECTION

START_SECTION((std::vector<ProteinIdentification>& getProteinIdentifications()))
  ExperimentalSettings settings;
  ProteinIdentification id;
	ProteinHit protein_hit;
	float protein_significance_threshold = 63.2f;

	id.setDateTime(DateTime::now());
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);
	
	settings.getProteinIdentifications().push_back(id);
	ProteinIdentification& test_id = settings.getProteinIdentifications()[0];
	TEST_EQUAL(id == test_id, true)
END_SECTION

START_SECTION((void setProteinIdentifications(const std::vector<ProteinIdentification>& protein_identifications)))
  ExperimentalSettings settings;
  ProteinIdentification id;
	ProteinHit protein_hit;
	float protein_significance_threshold = 63.2f;
	vector<ProteinIdentification> ids;

	id.setDateTime(DateTime::now());
	id.setSignificanceThreshold(protein_significance_threshold);
	id.insertHit(protein_hit);
	ids.push_back(id);
	id.setSignificanceThreshold(21.f);
	ids.push_back(id);
	settings.setProteinIdentifications(ids);
	TEST_EQUAL(ids == settings.getProteinIdentifications(), true)
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



