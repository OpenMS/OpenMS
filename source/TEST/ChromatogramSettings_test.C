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

START_SECTION(~ChromatogramSettings())
{
	delete ptr;
}
END_SECTION

START_SECTION((ChromatogramSettings(const ChromatogramSettings &source)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual ~ChromatogramSettings()))
{
  // TODO
}
END_SECTION

START_SECTION((ChromatogramSettings& operator=(const ChromatogramSettings &source)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const ChromatogramSettings &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const ChromatogramSettings &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getNativeID() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setNativeID(const String &native_id)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getComment() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setComment(const String &comment)))
{
  // TODO
}
END_SECTION

START_SECTION((const InstrumentSettings& getInstrumentSettings() const ))
{
  // TODO
}
END_SECTION

START_SECTION((InstrumentSettings& getInstrumentSettings()))
{
  // TODO
}
END_SECTION

START_SECTION((void setInstrumentSettings(const InstrumentSettings &instrument_settings)))
{
  // TODO
}
END_SECTION

START_SECTION((const AcquisitionInfo& getAcquisitionInfo() const ))
{
  // TODO
}
END_SECTION

START_SECTION((AcquisitionInfo& getAcquisitionInfo()))
{
  // TODO
}
END_SECTION

START_SECTION((void setAcquisitionInfo(const AcquisitionInfo &acquisition_info)))
{
  // TODO
}
END_SECTION

START_SECTION((const SourceFile& getSourceFile() const ))
{
  // TODO
}
END_SECTION

START_SECTION((SourceFile& getSourceFile()))
{
  // TODO
}
END_SECTION

START_SECTION((void setSourceFile(const SourceFile &source_file)))
{
  // TODO
}
END_SECTION

START_SECTION((const Precursor& getPrecursor() const ))
{
  // TODO
}
END_SECTION

START_SECTION((Precursor& getPrecursor()))
{
  // TODO
}
END_SECTION

START_SECTION((void setPrecursor(const Precursor &precursor)))
{
  // TODO
}
END_SECTION

START_SECTION((const Product& getProduct() const ))
{
  // TODO
}
END_SECTION

START_SECTION((Product& getProduct()))
{
  // TODO
}
END_SECTION

START_SECTION((void setProduct(const Product &product)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<DataProcessing>& getDataProcessing() const ))
{
  // TODO
}
END_SECTION

START_SECTION((std::vector<DataProcessing>& getDataProcessing()))
{
  // TODO
}
END_SECTION

START_SECTION((void setDataProcessing(const std::vector< DataProcessing > &data_processing)))
{
  // TODO
}
END_SECTION

START_SECTION((ChromatogramType getChromatogramType() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setChromatogramType(ChromatogramType type)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



