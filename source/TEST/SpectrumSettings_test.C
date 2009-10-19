// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SpectrumSettings, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumSettings* ptr = 0;

START_SECTION((SpectrumSettings()))
	ptr = new SpectrumSettings();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~SpectrumSettings()))
	delete ptr;
END_SECTION

START_SECTION((const String& getNativeID() const))
	SpectrumSettings tmp;
	TEST_STRING_EQUAL(tmp.getNativeID(),"");
END_SECTION

START_SECTION((void setNativeID(const String& native_id)))
	SpectrumSettings tmp;
	tmp.setNativeID("nid");
	TEST_STRING_EQUAL(tmp.getNativeID(),"nid");
END_SECTION

START_SECTION((const std::vector<DataProcessing>& getDataProcessing() const))
  SpectrumSettings tmp;
  TEST_EQUAL(tmp.getDataProcessing().size(),0);
END_SECTION

START_SECTION((void setDataProcessing(const std::vector< DataProcessing > &data_processing)))
  SpectrumSettings tmp;
  std::vector<DataProcessing> dummy;
  dummy.resize(1);
  tmp.setDataProcessing(dummy);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
END_SECTION

START_SECTION((std::vector<DataProcessing>& getDataProcessing()))
  SpectrumSettings tmp;
  tmp.getDataProcessing().resize(1);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
END_SECTION

START_SECTION((AcquisitionInfo& getAcquisitionInfo()))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getAcquisitionInfo()==AcquisitionInfo(), true);
END_SECTION

START_SECTION((void setAcquisitionInfo(const AcquisitionInfo& acquisition_info)))
	SpectrumSettings tmp;
	AcquisitionInfo ai;
	ai.setMethodOfCombination("test");
	tmp.setAcquisitionInfo(ai);
	TEST_EQUAL(tmp.getAcquisitionInfo()==AcquisitionInfo(), false);
END_SECTION

START_SECTION((const AcquisitionInfo& getAcquisitionInfo() const))
	SpectrumSettings tmp;
	tmp.getAcquisitionInfo().setMethodOfCombination("test");
	TEST_EQUAL(tmp.getAcquisitionInfo()==AcquisitionInfo(), false);
END_SECTION

START_SECTION((SourceFile& getSourceFile()))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getSourceFile()==SourceFile(), true);
END_SECTION

START_SECTION((void setSourceFile(const SourceFile& source_file)))
	SpectrumSettings tmp;
	SourceFile sf;
	sf.setNameOfFile("test");
	tmp.setSourceFile(sf);
	TEST_EQUAL(tmp.getSourceFile()==SourceFile(), false);
END_SECTION

START_SECTION((const SourceFile& getSourceFile() const))
	SpectrumSettings tmp;
	tmp.getSourceFile().setNameOfFile("test");
	TEST_EQUAL(tmp.getSourceFile()==SourceFile(), false);
END_SECTION

START_SECTION((const InstrumentSettings& getInstrumentSettings() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), true);	  
END_SECTION

START_SECTION((void setInstrumentSettings(const InstrumentSettings& instrument_settings)))
	SpectrumSettings tmp;
	InstrumentSettings is;
	is.getScanWindows().resize(1);
	tmp.setInstrumentSettings(is);
	TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), false);
END_SECTION

START_SECTION((InstrumentSettings& getInstrumentSettings()))
	SpectrumSettings tmp;
	tmp.getInstrumentSettings().getScanWindows().resize(1);
	TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), false);	
END_SECTION

START_SECTION((const std::vector<Precursor>& getPrecursors() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getPrecursors().size(),0);	  
END_SECTION

START_SECTION((void setPrecursors(const std::vector<Precursor>& precursors)))
	SpectrumSettings tmp;
	tmp.setPrecursors(vector<Precursor>(2));
	TEST_EQUAL(tmp.getPrecursors().size(), 2);
END_SECTION

START_SECTION((std::vector<Precursor>& getPrecursors()))
	SpectrumSettings tmp;
	tmp.getPrecursors().resize(4);
	TEST_EQUAL(tmp.getPrecursors().size(), 4);	
END_SECTION

START_SECTION((const std::vector<Product>& getProducts() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getProducts().size(),0);	  
END_SECTION

START_SECTION((void setProducts(const std::vector<Product>& products)))
	SpectrumSettings tmp;
	tmp.setProducts(vector<Product>(2));
	TEST_EQUAL(tmp.getProducts().size(), 2);
END_SECTION

START_SECTION((std::vector<Product>& getProducts()))
	SpectrumSettings tmp;
	tmp.getProducts().resize(4);
	TEST_EQUAL(tmp.getProducts().size(), 4);	
END_SECTION

START_SECTION((SpectrumType getType() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getType(), SpectrumSettings::UNKNOWN);	  
END_SECTION

START_SECTION((void setType(SpectrumType type)))
	SpectrumSettings tmp;
	tmp.setType(SpectrumSettings::PEAKS);
	TEST_EQUAL(tmp.getType(), SpectrumSettings::PEAKS);
END_SECTION

START_SECTION((const String& getComment() const))
	SpectrumSettings tmp;
	TEST_EQUAL(tmp.getComment(), "");
END_SECTION

START_SECTION((void setComment(const String& comment)))
	SpectrumSettings tmp;
	tmp.setComment("bla");
	TEST_EQUAL(tmp.getComment(), "bla");
END_SECTION

START_SECTION((const std::vector<PeptideIdentification>& getPeptideIdentifications() const))
	SpectrumSettings tmp;
	vector<PeptideIdentification> vec(tmp.getPeptideIdentifications());
	TEST_EQUAL(vec.size(),0);
END_SECTION

START_SECTION((void setPeptideIdentifications(const std::vector<PeptideIdentification>& identifications)))
	SpectrumSettings tmp;
	vector<PeptideIdentification> vec;
	
	tmp.setPeptideIdentifications(vec);
	TEST_EQUAL(tmp.getPeptideIdentifications().size(),0);
	
	PeptideIdentification dbs;
	vec.push_back(dbs);
	tmp.setPeptideIdentifications(vec);
	TEST_EQUAL(tmp.getPeptideIdentifications().size(),1);
END_SECTION

START_SECTION((std::vector<PeptideIdentification>& getPeptideIdentifications()))
	SpectrumSettings tmp;
	vector<PeptideIdentification> vec;
	
	tmp.getPeptideIdentifications().resize(1);
	TEST_EQUAL(tmp.getPeptideIdentifications().size(),1);
END_SECTION

START_SECTION((SpectrumSettings& operator= (const SpectrumSettings& source)))
  SpectrumSettings tmp;
  tmp.setMetaValue("bla","bluff");
	tmp.getAcquisitionInfo().setMethodOfCombination("test");
	tmp.getInstrumentSettings().getScanWindows().resize(1);
	tmp.getPrecursors().resize(1);
	tmp.getProducts().resize(1);
	tmp.getPeptideIdentifications().resize(1);
	tmp.setType(SpectrumSettings::PEAKS);
	tmp.setComment("bla");
	tmp.setNativeID("nid");
	tmp.getDataProcessing().resize(1);
	
	SpectrumSettings tmp2(tmp);
	TEST_EQUAL(tmp2.getComment(), "bla");
	TEST_EQUAL(tmp2.getType(), SpectrumSettings::PEAKS);
	TEST_EQUAL(tmp2.getPeptideIdentifications().size(), 1);	
	TEST_EQUAL(tmp2.getPrecursors().size(),1);	
	TEST_EQUAL(tmp2.getProducts().size(),1);	
	TEST_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), false);
	TEST_EQUAL(tmp2.getAcquisitionInfo()==AcquisitionInfo(), false);  
	TEST_STRING_EQUAL(tmp2.getNativeID(),"nid");
	TEST_EQUAL(tmp2.getDataProcessing().size(),1);
	TEST_EQUAL(tmp2.getMetaValue("bla")=="bluff",true);
END_SECTION

START_SECTION((SpectrumSettings(const SpectrumSettings& source)))
  SpectrumSettings tmp;
	tmp.getAcquisitionInfo().setMethodOfCombination("test");
	tmp.getInstrumentSettings().getScanWindows().resize(1);
	tmp.getPrecursors().resize(1);
	tmp.getProducts().resize(1);
	tmp.setType(SpectrumSettings::PEAKS);
	tmp.setComment("bla");
	tmp.getPeptideIdentifications().resize(1);
	tmp.setNativeID("nid");
	tmp.getDataProcessing().resize(1);
	tmp.setMetaValue("bla","bluff");
	
	SpectrumSettings tmp2;
	tmp2 = tmp;
	TEST_EQUAL(tmp2.getComment(), "bla");
	TEST_EQUAL(tmp2.getType(), SpectrumSettings::PEAKS);
	TEST_EQUAL(tmp2.getPrecursors().size(), 1);
	TEST_EQUAL(tmp2.getProducts().size(),1);	
	TEST_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), false);	
	TEST_EQUAL(tmp2.getAcquisitionInfo()==AcquisitionInfo(), false);
	TEST_EQUAL(tmp2.getPeptideIdentifications().size(), 1);	
	TEST_STRING_EQUAL(tmp2.getNativeID(),"nid");
	TEST_EQUAL(tmp2.getDataProcessing().size(),1);
	TEST_STRING_EQUAL(tmp2.getMetaValue("bla"),"bluff");


	tmp2 = SpectrumSettings();
	TEST_EQUAL(tmp2.getComment(), "");
	TEST_EQUAL(tmp2.getType(), SpectrumSettings::UNKNOWN);
	TEST_EQUAL(tmp2.getPrecursors().size(),0);	
	TEST_EQUAL(tmp2.getProducts().size(),0);	
	TEST_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), true);	
	TEST_EQUAL(tmp2.getAcquisitionInfo()==AcquisitionInfo(), true);
	TEST_EQUAL(tmp2.getPeptideIdentifications().size(), 0);	
	TEST_STRING_EQUAL(tmp2.getNativeID(),"");
	TEST_EQUAL(tmp2.getDataProcessing().size(),0);
	TEST_EQUAL(tmp2.metaValueExists("bla"),false);

END_SECTION

START_SECTION((bool operator== (const SpectrumSettings& rhs) const))
  SpectrumSettings edit, empty;
  
  TEST_EQUAL(edit==empty, true);
  
	edit.getAcquisitionInfo().setMethodOfCombination("test");
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.setNativeID("nid");
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.getInstrumentSettings().getScanWindows().resize(1);
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.getPrecursors().resize(1);
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.setType(SpectrumSettings::PEAKS);
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.setComment("bla");
	TEST_EQUAL(edit==empty, false);

	edit = empty;
	edit.getPrecursors().resize(1);
	TEST_EQUAL(edit==empty, false);

	edit = empty;
	edit.getProducts().resize(1);
	TEST_EQUAL(edit==empty, false);
	
	edit = empty;
	edit.getPeptideIdentifications().resize(1);
	TEST_EQUAL(edit==empty, false);

	edit = empty;
	edit.getDataProcessing().resize(1);
	TEST_EQUAL(edit==empty, false);

	edit = empty;
	edit.setMetaValue("bla","bluff");
	TEST_EQUAL(edit==empty, false);

END_SECTION

START_SECTION((bool operator!= (const SpectrumSettings& rhs) const))
  SpectrumSettings edit, empty;
  
  TEST_EQUAL(edit!=empty, false);
  
	edit.getAcquisitionInfo().setMethodOfCombination("test");
	TEST_EQUAL(edit!=empty, true);
	
	edit = empty;
	edit.setNativeID("nid");
	TEST_EQUAL(edit!=empty, true);
	
	edit = empty;
	edit.getInstrumentSettings().getScanWindows().resize(1);
	TEST_EQUAL(edit!=empty, true);
	
	edit = empty;
	edit.getPrecursors().resize(1);
	TEST_EQUAL(edit!=empty, true);
	
	edit = empty;
	edit.setType(SpectrumSettings::PEAKS);
	TEST_EQUAL(edit!=empty, true);
	
	edit = empty;
	edit.setComment("bla");
	TEST_EQUAL(edit!=empty, true);

	edit = empty;
	edit.getPrecursors().resize(1);
	TEST_EQUAL(edit!=empty, true);

	edit = empty;
	edit.getProducts().resize(1);
	TEST_EQUAL(edit!=empty, true);

	edit = empty;
	edit.getPeptideIdentifications().resize(1);
	TEST_EQUAL(edit!=empty, true);

	edit = empty;
	edit.getDataProcessing().resize(1);
	TEST_EQUAL(edit!=empty, true);

	edit = empty;
	edit.setMetaValue("bla","bluff");
	TEST_EQUAL(edit!=empty, true);


END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



