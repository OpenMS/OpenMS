// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/ChromatogramSettings.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ChromatogramSettings, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ChromatogramSettings* ptr = 0;
START_SECTION(ChromatogramSettings())
{
	ptr = new ChromatogramSettings();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(virtual ~ChromatogramSettings())
{
	delete ptr;
}
END_SECTION

START_SECTION((ChromatogramSettings(const ChromatogramSettings &source)))
{
  ChromatogramSettings tmp;
  tmp.getAcquisitionInfo().setMethodOfCombination("test");
  tmp.getInstrumentSettings().getScanWindows().resize(1);
  tmp.getPrecursor().setMZ(0.11);
  tmp.getProduct().setMZ(0.12);
  tmp.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  tmp.setComment("bla");
  tmp.setNativeID("nid");
  tmp.getDataProcessing().resize(1);
  tmp.setMetaValue("bla","bluff");

  ChromatogramSettings tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getComment(), "bla");
  TEST_EQUAL(tmp2.getChromatogramType(), ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  TEST_REAL_SIMILAR(tmp2.getPrecursor().getMZ(), 0.11);
  TEST_REAL_SIMILAR(tmp2.getProduct().getMZ(), 0.12);
  TEST_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), false);
  TEST_EQUAL(tmp2.getAcquisitionInfo()==AcquisitionInfo(), false);
  TEST_STRING_EQUAL(tmp2.getNativeID(),"nid");
  TEST_EQUAL(tmp2.getDataProcessing().size(),1);
  TEST_STRING_EQUAL(tmp2.getMetaValue("bla"),"bluff");


  tmp2 = ChromatogramSettings();
  TEST_EQUAL(tmp2.getComment(), "");
  TEST_EQUAL(tmp2.getChromatogramType(), ChromatogramSettings::MASS_CHROMATOGRAM);
  TEST_REAL_SIMILAR(tmp2.getPrecursor().getMZ(), 0.0);
  TEST_REAL_SIMILAR(tmp2.getProduct().getMZ(), 0.0);
  TEST_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), true);
  TEST_EQUAL(tmp2.getAcquisitionInfo()==AcquisitionInfo(), true);
  TEST_STRING_EQUAL(tmp2.getNativeID(),"");
  TEST_EQUAL(tmp2.getDataProcessing().size(),0);
  TEST_EQUAL(tmp2.metaValueExists("bla"),false);

}
END_SECTION

START_SECTION((ChromatogramSettings& operator=(const ChromatogramSettings &source)))
{
  ChromatogramSettings tmp;
  tmp.setMetaValue("bla","bluff");
  tmp.getAcquisitionInfo().setMethodOfCombination("test");
  tmp.getInstrumentSettings().getScanWindows().resize(1);
  tmp.getPrecursor().setMZ(0.13);
  tmp.getProduct().setMZ(0.14);
  tmp.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  tmp.setComment("bla");
  tmp.setNativeID("nid");
  tmp.getDataProcessing().resize(1);

  ChromatogramSettings tmp2(tmp);
  TEST_EQUAL(tmp2.getComment(), "bla");
  TEST_EQUAL(tmp2.getChromatogramType(), ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  TEST_REAL_SIMILAR(tmp2.getPrecursor().getMZ(), 0.13);
  TEST_REAL_SIMILAR(tmp2.getProduct().getMZ(), 0.14);
  TEST_EQUAL(tmp2.getInstrumentSettings()==InstrumentSettings(), false);
  TEST_EQUAL(tmp2.getAcquisitionInfo()==AcquisitionInfo(), false);
  TEST_STRING_EQUAL(tmp2.getNativeID(),"nid");
  TEST_EQUAL(tmp2.getDataProcessing().size(),1);
  TEST_EQUAL(tmp2.getMetaValue("bla")=="bluff",true);
}
END_SECTION

START_SECTION((bool operator==(const ChromatogramSettings &rhs) const ))
{
  ChromatogramSettings edit, empty;

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
  edit.setComment("comment");
  TEST_EQUAL(edit == empty, false)

  edit = empty;
  edit.getPrecursor().setMZ(0.15);
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.setComment("bla");
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.getPrecursor().setMZ(0.16);
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.getProduct().setMZ(0.17);
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.getDataProcessing().resize(1);
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.setMetaValue("bla","bluff");
  TEST_EQUAL(edit==empty, false);


}
END_SECTION

START_SECTION((bool operator!=(const ChromatogramSettings &rhs) const ))
{
  ChromatogramSettings edit, empty;

  TEST_EQUAL(edit!=empty, false);

  edit.getAcquisitionInfo().setMethodOfCombination("test");
  TEST_EQUAL(edit!=empty, true);

  edit = empty;
  edit.setNativeID("nid");
  TEST_EQUAL(edit!=empty, true)

  edit = empty;
  edit.getInstrumentSettings().getScanWindows().resize(1);
  TEST_EQUAL(edit!=empty, true);

	edit = empty;
	edit.setComment("comment");
	TEST_EQUAL(edit != empty, true)

  edit = empty;
  edit.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
  TEST_EQUAL(edit!=empty, true);

  edit = empty;
  edit.setComment("bla");
  TEST_EQUAL(edit!=empty, true)

  edit = empty;
	Precursor prec;
	prec.setMZ(1.3);
  edit.setPrecursor(prec);
  TEST_EQUAL(edit!=empty, true)

  edit = empty;
	Product prod;
	prod.setMZ(1.5);
  edit.setProduct(prod);
  TEST_EQUAL(edit!=empty, true)

  edit = empty;
  edit.getDataProcessing().resize(1);
  TEST_EQUAL(edit!=empty, true)

  edit = empty;
  edit.setMetaValue("bla","bluff");
  TEST_EQUAL(edit!=empty, true)

}
END_SECTION

START_SECTION((const String& getNativeID() const ))
{
  ChromatogramSettings tmp;
  TEST_STRING_EQUAL(tmp.getNativeID(),"")
}
END_SECTION

START_SECTION((void setNativeID(const String &native_id)))
{
  ChromatogramSettings tmp;
  tmp.setNativeID("nid");
  TEST_STRING_EQUAL(tmp.getNativeID(),"nid")
}
END_SECTION

START_SECTION((const String& getComment() const ))
{
  ChromatogramSettings tmp;
	TEST_STRING_EQUAL(tmp.getComment(), "")
}
END_SECTION

START_SECTION((void setComment(const String &comment)))
{
  ChromatogramSettings tmp;
	tmp.setComment("name");
	TEST_STRING_EQUAL(tmp.getComment(), "name")
}
END_SECTION

START_SECTION((const InstrumentSettings& getInstrumentSettings() const ))
{
  ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), true);
}
END_SECTION

START_SECTION((InstrumentSettings& getInstrumentSettings()))
{
  ChromatogramSettings tmp;
  tmp.getInstrumentSettings().getScanWindows().resize(1);
  TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), false);
}
END_SECTION

START_SECTION((void setInstrumentSettings(const InstrumentSettings &instrument_settings)))
{
  ChromatogramSettings tmp;
  InstrumentSettings is;
  is.getScanWindows().resize(1);
  tmp.setInstrumentSettings(is);
  TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), false);
}
END_SECTION

START_SECTION((const AcquisitionInfo& getAcquisitionInfo() const ))
{
  ChromatogramSettings tmp;
  tmp.getAcquisitionInfo().setMethodOfCombination("test");
  TEST_EQUAL(tmp.getAcquisitionInfo()==AcquisitionInfo(), false);
}
END_SECTION

START_SECTION((AcquisitionInfo& getAcquisitionInfo()))
{
  ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getAcquisitionInfo()==AcquisitionInfo(), true);
}
END_SECTION

START_SECTION((void setAcquisitionInfo(const AcquisitionInfo &acquisition_info)))
{
  ChromatogramSettings tmp;
  AcquisitionInfo ai;
  ai.setMethodOfCombination("test");
  tmp.setAcquisitionInfo(ai);
  TEST_EQUAL(tmp.getAcquisitionInfo()==AcquisitionInfo(), false);
}
END_SECTION

START_SECTION((const SourceFile& getSourceFile() const ))
{
	ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getInstrumentSettings()==InstrumentSettings(), true);
}
END_SECTION

START_SECTION((SourceFile& getSourceFile()))
{
  ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getSourceFile()==SourceFile(), true);
}
END_SECTION

START_SECTION((void setSourceFile(const SourceFile &source_file)))
{
  ChromatogramSettings tmp;
  SourceFile sf;
  sf.setNameOfFile("test");
  tmp.setSourceFile(sf);
  TEST_EQUAL(tmp.getSourceFile()==SourceFile(), false);
}
END_SECTION

START_SECTION((const Precursor& getPrecursor() const ))
{
	ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getPrecursor() == Precursor(), true)
}
END_SECTION

START_SECTION((Precursor& getPrecursor()))
{
  ChromatogramSettings tmp;
	tmp.getPrecursor().setMZ(0.3);
	TEST_EQUAL(tmp.getPrecursor() == Precursor(), false)
	TEST_REAL_SIMILAR(tmp.getPrecursor().getMZ(), 0.3)
}
END_SECTION

START_SECTION((void setPrecursor(const Precursor &precursor)))
{
  ChromatogramSettings tmp;
	Precursor prec;
	prec.setMZ(0.4);
	tmp.setPrecursor(prec);
	TEST_EQUAL(tmp.getPrecursor() == Precursor(), false)
	TEST_REAL_SIMILAR(tmp.getPrecursor().getMZ(), 0.4)
}
END_SECTION

START_SECTION((const Product& getProduct() const ))
{
  ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getProduct() == Product(), true)
}
END_SECTION

START_SECTION((Product& getProduct()))
{
  ChromatogramSettings tmp;
  tmp.getProduct().setMZ(0.3);
  TEST_EQUAL(tmp.getProduct() == Product(), false)
  TEST_REAL_SIMILAR(tmp.getProduct().getMZ(), 0.3)
}
END_SECTION

START_SECTION((void setProduct(const Product &product)))
{
  ChromatogramSettings tmp;
  Product prod;
  prod.setMZ(0.4);
  tmp.setProduct(prod);
  TEST_EQUAL(tmp.getProduct() == Product(), false)
  TEST_REAL_SIMILAR(tmp.getProduct().getMZ(), 0.4)
}
END_SECTION

START_SECTION((const std::vector<DataProcessing>& getDataProcessing() const ))
{
	ChromatogramSettings tmp;
  TEST_EQUAL(tmp.getDataProcessing().size(),0);
}
END_SECTION

START_SECTION((std::vector<DataProcessing>& getDataProcessing()))
{
  ChromatogramSettings tmp;
  tmp.getDataProcessing().resize(1);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
}
END_SECTION

START_SECTION((void setDataProcessing(const std::vector< DataProcessing > &data_processing)))
{
  ChromatogramSettings tmp;
  std::vector<DataProcessing> dummy;
  dummy.resize(1);
  tmp.setDataProcessing(dummy);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
}
END_SECTION

START_SECTION((ChromatogramType getChromatogramType() const ))
{
  ChromatogramSettings tmp;
	TEST_EQUAL(tmp.getChromatogramType(), ChromatogramSettings::MASS_CHROMATOGRAM)
}
END_SECTION

START_SECTION((void setChromatogramType(ChromatogramType type)))
{
  ChromatogramSettings tmp;
	tmp.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
	TEST_EQUAL(tmp.getChromatogramType(), ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



