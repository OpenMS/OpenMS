// -*- Mode: C++; tab-width: 2; -*-
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

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MzMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzMLFile* ptr = 0;
CHECK((MzMLFile()))
	ptr = new MzMLFile;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~MzMLFile()))
	delete ptr;
RESULT

CHECK(const PeakFileOptions& getOptions() const)
	MzMLFile file;
	TEST_EQUAL(file.getOptions().hasMSLevels(),false)
RESULT

CHECK(PeakFileOptions& getOptions())
	MzMLFile file;
	file.getOptions().addMSLevel(1);
	TEST_EQUAL(file.getOptions().hasMSLevels(),true);
RESULT

CHECK((template <typename MapType> void load(const String& filename, MapType& map)))
	MzMLFile file;
	MSExperiment<> exp;
	file.load("data/MzMLFile_1.mzML",exp);

	TEST_EQUAL(exp.size(),3)
	
	//-------------------------- spectrum 0 --------------------------
	
	TEST_EQUAL(exp[0].size(),15)
	TEST_EQUAL(exp[0].getMSLevel(),1)
	TEST_EQUAL(exp[0].getInstrumentSettings().getScanMode(),InstrumentSettings::PRODUCT)
	TEST_EQUAL(exp[0].getMetaDataArrays().size(),0)
	TEST_EQUAL(exp[0].getType(),SpectrumSettings::PEAKS)
	TEST_REAL_EQUAL(exp[0].getRT(),5.8905)
	//TODO TEST_EQUAL(exp[0].getInstrumentSettings().getPolarity(),IonSource::POSITIVE)
	TEST_REAL_EQUAL(exp[0].getInstrumentSettings().getMzRangeStart(),400.0)
	TEST_REAL_EQUAL(exp[0].getInstrumentSettings().getMzRangeStop(),1800.0)

	//-------------------------- spectrum 1 --------------------------
	
	TEST_EQUAL(exp[1].size(),10)
	TEST_EQUAL(exp[1].getMSLevel(),2)
	TEST_EQUAL(exp[1].getInstrumentSettings().getScanMode(),InstrumentSettings::PRODUCT)
	TEST_EQUAL(exp[1].getType(),SpectrumSettings::PEAKS)
	TEST_REAL_EQUAL(exp[1].getRT(),5.9905)
	//TODO TEST_EQUAL(exp[1].getInstrumentSettings().getPolarity(),IonSource::POSITIVE)
	TEST_REAL_EQUAL(exp[1].getInstrumentSettings().getMzRangeStart(),110.0)
	TEST_REAL_EQUAL(exp[1].getInstrumentSettings().getMzRangeStop(),905.0)
	
	//meta data arrays
	TEST_EQUAL(exp[1].getMetaDataArrays().size(),2)
	TEST_EQUAL(exp[1].getMetaDataArrays()[0].getName(),"signal to noise")
	TEST_EQUAL(exp[1].getMetaDataArrays()[0].size(),10)
	TEST_EQUAL(exp[1].getMetaDataArrays()[1].getName(),"charge")
	TEST_EQUAL(exp[1].getMetaDataArrays()[1].size(),10)
	
	//-------------------------- spectrum 2 --------------------------
	
	TEST_EQUAL(exp[2].size(),0)
	TEST_EQUAL(exp[2].getMSLevel(),1)
	TEST_EQUAL(exp[2].getInstrumentSettings().getScanMode(),InstrumentSettings::PRODUCT)
	TEST_EQUAL(exp[2].getMetaDataArrays().size(),0)
	TEST_EQUAL(exp[2].getType(),SpectrumSettings::UNKNOWN)
	TEST_REAL_EQUAL(exp[2].getRT(),-1.0)
	TEST_EQUAL(exp[2].getInstrumentSettings().getPolarity(),IonSource::POLNULL)
	TEST_REAL_EQUAL(exp[2].getInstrumentSettings().getMzRangeStart(),0.0)
	TEST_REAL_EQUAL(exp[2].getInstrumentSettings().getMzRangeStop(),0.0)

RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

