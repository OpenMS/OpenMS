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
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/METADATA/Identification.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SpectrumSettings, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumSettings* ptr = 0;

CHECK((SpectrumSettings()))
	ptr = new SpectrumSettings();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~SpectrumSettings()))
	delete ptr;
RESULT

CHECK((AcquisitionInfo& getAcquisitionInfo()))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getAcquisitionInfo()==AcquisitionInfo(), true);
RESULT

CHECK((void setAcquisitionInfo(const AcquisitionInfo& acquisition_info)))
	SpectrumSettings tmp;
	AcquisitionInfo ai;
	ai.setMethodOfCombination("test");
	tmp.setAcquisitionInfo(ai);
	TEST_EQUAL(tmp.getAcquisitionInfo()==AcquisitionInfo(), false);
RESULT

CHECK((const AcquisitionInfo& getAcquisitionInfo() const))
	SpectrumSettings tmp;
	tmp.getAcquisitionInfo().setMethodOfCombination("test");
	TEST_EQUAL(tmp.getAcquisitionInfo()==AcquisitionInfo(), false);
RESULT

CHECK((SourceFile& getSourceFile()))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getSourceFile()==SourceFile(), true);
RESULT

CHECK((void setSourceFile(const SourceFile& source_file)))
	SpectrumSettings tmp;
	SourceFile sf;
	sf.setNameOfFile("test");
	tmp.setSourceFile(sf);
	TEST_EQUAL(tmp.getSourceFile()==SourceFile(), false);
RESULT

CHECK((const SourceFile& getSourceFile() const))
	SpectrumSettings tmp;
	tmp.getSourceFile().setNameOfFile("test");
	TEST_EQUAL(tmp.getSourceFile()==SourceFile(), false);
RESULT

CHECK((const InstrumentSettings& getInstrumentSettings() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), true);	  
RESULT

CHECK((void setInstrumentSettings(const InstrumentSettings& instrument_settings)))
	SpectrumSettings tmp;
	InstrumentSettings is;
	is.setMzRangeStart(47.11);
	tmp.setInstrumentSettings(is);
	TEST_REAL_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), false);
RESULT

CHECK((InstrumentSettings& getInstrumentSettings()))
	SpectrumSettings tmp;
	tmp.getInstrumentSettings().setMzRangeStart(47.11);
	TEST_REAL_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), false);	
RESULT

CHECK((const Precursor& getPrecursor() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getPrecursor()==Precursor(), true);	  
RESULT

CHECK((void setPrecursor(const Precursor& precursor)))
	SpectrumSettings tmp;
	Precursor is;
	is.setActivationEnergy(47.11);
	tmp.setPrecursor(is);
	TEST_REAL_EQUAL(tmp.getPrecursor()==Precursor(), false);
RESULT

CHECK((Precursor& getPrecursor()))
	SpectrumSettings tmp;
	tmp.getPrecursor().setActivationEnergy(47.11);
	TEST_REAL_EQUAL(tmp.getPrecursor()==Precursor(), false);	
RESULT

CHECK((SpectrumType getType() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getType(), SpectrumSettings::UNKNOWN);	  
RESULT

CHECK((void setType(SpectrumType type)))
	SpectrumSettings tmp;
	tmp.setType(SpectrumSettings::PEAKS);
	TEST_EQUAL(tmp.getType(), SpectrumSettings::PEAKS);
RESULT

CHECK((const String& getComment() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getComment(), "");
RESULT

CHECK((void setComment(const String& comment)))
	SpectrumSettings tmp;
	tmp.setComment("bla");
	TEST_EQUAL(tmp.getComment(), "bla");
RESULT

CHECK((const std::vector<Identification>& getIdentifications() const))
	SpectrumSettings tmp;
	vector<Identification> vec(tmp.getIdentifications());
	TEST_EQUAL(vec.size(),0);
RESULT

CHECK((void setIdentifications(const std::vector<Identification>& identifications)))
	SpectrumSettings tmp;
	vector<Identification> vec;
	
	tmp.setIdentifications(vec);
	TEST_EQUAL(tmp.getIdentifications().size(),0);
	
	Identification dbs;
	vec.push_back(dbs);
	tmp.setIdentifications(vec);
	TEST_EQUAL(tmp.getIdentifications().size(),1);
RESULT

CHECK((std::vector<Identification>& getIdentifications()))
	SpectrumSettings tmp;
	vector<Identification> vec;
	
	tmp.getIdentifications().resize(1);
	TEST_EQUAL(tmp.getIdentifications().size(),1);
RESULT

CHECK((const std::map<String,MetaInfoDescription>& getMetaInfoDescriptions() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getMetaInfoDescriptions().size(),0);
RESULT

CHECK((void setMetaInfoDescriptions(const std::map<String, MetaInfoDescription>& meta_info_descriptions)))
	SpectrumSettings tmp;
	std::map<String, MetaInfoDescription> mid;
	mid["key"].setComment("comment");
	tmp.setMetaInfoDescriptions(mid);
	TEST_EQUAL(tmp.getMetaInfoDescriptions().size(),1);
	TEST_EQUAL(tmp.getMetaInfoDescriptions()["key"].getComment(),"comment");	
RESULT

CHECK((std::map<String,MetaInfoDescription>& getMetaInfoDescriptions()))
	SpectrumSettings tmp;
	std::map<String, MetaInfoDescription> mid;
	mid["key"].setComment("comment");
	tmp.setMetaInfoDescriptions(mid);
	tmp.getMetaInfoDescriptions()["key"].setComment("comment2");
	tmp.getMetaInfoDescriptions()["key2"].setComment("comment");
	TEST_EQUAL(tmp.getMetaInfoDescriptions().size(),2);
	TEST_EQUAL(tmp.getMetaInfoDescriptions()["key2"].getComment(),"comment");
RESULT

CHECK((SpectrumSettings& operator= (const SpectrumSettings& source)))
  SpectrumSettings tmp;
	tmp.getAcquisitionInfo().setMethodOfCombination("test");
	tmp.getInstrumentSettings().setMzRangeStart(47.11);
	tmp.getPrecursor().setActivationEnergy(47.11);
	tmp.getIdentifications().resize(1);
	tmp.setType(SpectrumSettings::PEAKS);
	tmp.setComment("bla");
	
	SpectrumSettings tmp2(tmp);
	TEST_EQUAL(tmp2.getComment(), "bla");
	TEST_EQUAL(tmp2.getType(), SpectrumSettings::PEAKS);
	TEST_EQUAL(tmp2.getIdentifications().size(), 1);	
	TEST_REAL_EQUAL(tmp2.getPrecursor()==Precursor(), false);	
	TEST_REAL_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), false);
	TEST_EQUAL(tmp2.getAcquisitionInfo()==AcquisitionInfo(), false);  
RESULT

CHECK((SpectrumSettings(const SpectrumSettings& source)))
  SpectrumSettings tmp;
	tmp.getAcquisitionInfo().setMethodOfCombination("test");
	tmp.getInstrumentSettings().setMzRangeStart(47.11);
	tmp.getPrecursor().setActivationEnergy(47.11);
	tmp.setType(SpectrumSettings::PEAKS);
	tmp.setComment("bla");
	tmp.getIdentifications().resize(1);
	
	SpectrumSettings tmp2;
	tmp2 = tmp;
	TEST_EQUAL(tmp2.getComment(), "bla");
	TEST_EQUAL(tmp2.getType(), SpectrumSettings::PEAKS);
	TEST_REAL_EQUAL(tmp2.getPrecursor()==Precursor(), false);	
	TEST_REAL_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), false);	
	TEST_EQUAL(tmp2.getAcquisitionInfo()==AcquisitionInfo(), false);
	TEST_EQUAL(tmp2.getIdentifications().size(), 1);	

	tmp2 = SpectrumSettings();
	TEST_EQUAL(tmp2.getComment(), "");
	TEST_EQUAL(tmp2.getType(), SpectrumSettings::UNKNOWN);
	TEST_REAL_EQUAL(tmp2.getPrecursor()==Precursor(), true);	
	TEST_REAL_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), true);	
	TEST_EQUAL(tmp2.getAcquisitionInfo()==AcquisitionInfo(), true);
	TEST_EQUAL(tmp2.getIdentifications().size(), 0);	

RESULT

CHECK((bool operator== (const SpectrumSettings& rhs) const))
  SpectrumSettings edit, empty;
  
  TEST_EQUAL(edit==empty, true);
  
	edit.getAcquisitionInfo().setMethodOfCombination("test");
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.getInstrumentSettings().setMzRangeStart(47.11);
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.getPrecursor().setActivationEnergy(47.11);
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.setType(SpectrumSettings::PEAKS);
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.setComment("bla");
	TEST_EQUAL(edit==empty, false);

	edit = empty;
	edit.getIdentifications().resize(1);
	TEST_EQUAL(edit==empty, false);
RESULT

CHECK((bool operator!= (const SpectrumSettings& rhs) const))
  SpectrumSettings edit, empty;
  
  TEST_EQUAL(edit!=empty, false);
  
	edit.getAcquisitionInfo().setMethodOfCombination("test");
	TEST_EQUAL(edit!=empty, true);
	
	edit = empty;
	edit.getInstrumentSettings().setMzRangeStart(47.11);
	TEST_EQUAL(edit!=empty, true);
	
	edit = empty;
	edit.getPrecursor().setActivationEnergy(47.11);
	TEST_EQUAL(edit!=empty, true);
	
	edit = empty;
	edit.setType(SpectrumSettings::PEAKS);
	TEST_EQUAL(edit!=empty, true);
	
	edit = empty;
	edit.setComment("bla");
	TEST_EQUAL(edit!=empty, true);
	
	edit = empty;
	edit.getIdentifications().resize(1);
	TEST_EQUAL(edit!=empty, true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



