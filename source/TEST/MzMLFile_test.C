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

PRECISION(0.01)

CHECK((template <typename MapType> void load(const String& filename, MapType& map)))
	MzMLFile file;
	MSExperiment<> exp;
	file.load("data/MzMLFile_1.mzML",exp);

	//-------------------------- general information --------------------------
  TEST_EQUAL(exp.getMetaValue("flag").valueType(),DataValue::STRING_VALUE)
  TEST_STRING_EQUAL(exp.getMetaValue("flag").toString(),"")
  TEST_EQUAL(exp.getMetaValue("string").valueType(),DataValue::STRING_VALUE)
  TEST_STRING_EQUAL((String)exp.getMetaValue("string"),"bla")
  TEST_EQUAL(exp.getMetaValue("float").valueType(),DataValue::DOUBLE_VALUE)
  TEST_REAL_EQUAL((double)exp.getMetaValue("float"),5.11)
  TEST_EQUAL(exp.getMetaValue("int").valueType(),DataValue::INT_VALUE)
  TEST_EQUAL((Int)exp.getMetaValue("int"),5)
  		
	TEST_EQUAL(exp.size(),3)
	
	//contacts
	TEST_EQUAL(exp.getContacts().size(),1)
	TEST_STRING_EQUAL(exp.getContacts()[0].getFirstName(),"William")
	TEST_STRING_EQUAL(exp.getContacts()[0].getLastName(),"Pennington")
	TEST_STRING_EQUAL(exp.getContacts()[0].getEmail(),"wpennington@higglesworth.edu")
	//source files
	TEST_STRING_EQUAL(exp.getSourceFile().getNameOfFile(),"tiny1.RAW")
	TEST_STRING_EQUAL(exp.getSourceFile().getPathToFile(),"file:///F:/data/Exp01")
	TEST_STRING_EQUAL(exp.getSourceFile().getSha1(),"71be39fb2700ab2f3c8b2234b91274968b6899b1")
	//sample
	TEST_STRING_EQUAL(exp.getSample().getName(),"Sample1")
	TEST_REAL_EQUAL(exp.getSample().getMass(),11.7)
	TEST_STRING_EQUAL(exp.getSample().getNumber(),"5")
	TEST_REAL_EQUAL(exp.getSample().getVolume(),3.1)
	TEST_REAL_EQUAL(exp.getSample().getConcentration(),5.5)
	//instrument (general)
	TEST_STRING_EQUAL(exp.getInstrument().getName(),"LCQ Deca")
	TEST_STRING_EQUAL(exp.getInstrument().getCustomizations(),"Umbau")
	//ion source
	TEST_EQUAL(exp.getInstrument().getIonSource().getInletType(),IonSource::DIRECT)
	TEST_EQUAL(exp.getInstrument().getIonSource().getIonizationMethod(),IonSource::ESI)
	//mass analyzers
	TEST_EQUAL(exp.getInstrument().getMassAnalyzers().size(),2)
	TEST_EQUAL(exp.getInstrument().getMassAnalyzers()[0].getType(),MassAnalyzer::PAULIONTRAP)
	TEST_REAL_EQUAL(exp.getInstrument().getMassAnalyzers()[0].getMagneticFieldStrength(),14.56)
	TEST_EQUAL(exp.getInstrument().getMassAnalyzers()[1].getType(),MassAnalyzer::LIT)
	TEST_REAL_EQUAL(exp.getInstrument().getMassAnalyzers()[1].getMagneticFieldStrength(),1414.14)
	//detector
	TEST_EQUAL(exp.getInstrument().getIonDetector().getType(),IonDetector::ELECTRONMULTIPLIER)
	TEST_EQUAL(exp.getInstrument().getIonDetector().getAcquisitionMode(),IonDetector::TDC)
	TEST_REAL_EQUAL(exp.getInstrument().getIonDetector().getResolution(),5.1)
	TEST_REAL_EQUAL(exp.getInstrument().getIonDetector().getADCSamplingFrequency(),1.1)

	//-------------------------- spectrum 0 --------------------------
	{
		const MSSpectrum<>& spec = exp[0]; 
	
		TEST_EQUAL(spec.size(),15)
		TEST_EQUAL(spec.getMSLevel(),1)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::PRODUCT)
		TEST_EQUAL(spec.getMetaDataArrays().size(),0)
		TEST_EQUAL(spec.getType(),SpectrumSettings::PEAKS)
		TEST_REAL_EQUAL(spec.getRT(),5.8905)
		TEST_EQUAL(spec.getInstrumentSettings().getPolarity(),IonSource::POSITIVE)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getMzRangeStart(),400.0)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getMzRangeStop(),1800.0)
		TEST_STRING_EQUAL(spec.getAcquisitionInfo().getMethodOfCombination(),"median")
		TEST_EQUAL(spec.getAcquisitionInfo().size(),2)
		TEST_EQUAL(spec.getAcquisitionInfo()[0].getNumber(),4711)
		TEST_EQUAL(spec.getAcquisitionInfo()[1].getNumber(),4712)
		TEST_EQUAL(spec.getSourceFile()==SourceFile(),true)
	}
	
	//-------------------------- spectrum 1 --------------------------
	{
		const MSSpectrum<>& spec = exp[1]; 	
		TEST_EQUAL(spec.size(),10)
		TEST_EQUAL(spec.getMSLevel(),2)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::PRODUCT)
		TEST_EQUAL(spec.getType(),SpectrumSettings::PEAKS)
		TEST_REAL_EQUAL(spec.getRT(),5.9905)
		TEST_EQUAL(spec.getInstrumentSettings().getPolarity(),IonSource::POSITIVE)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getMzRangeStart(),110.0)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getMzRangeStop(),905.0)
		TEST_EQUAL(spec.getAcquisitionInfo().getMethodOfCombination(),"")
		TEST_EQUAL(spec.getAcquisitionInfo().size(),0)
		TEST_EQUAL(spec.getSourceFile()==SourceFile(),true)
		
		//meta data arrays
		TEST_EQUAL(spec.getMetaDataArrays().size(),2)
		TEST_STRING_EQUAL(spec.getMetaDataArrays()[0].getName(),"signal to noise")
		TEST_EQUAL(spec.getMetaDataArrays()[0].size(),10)
		TEST_STRING_EQUAL(spec.getMetaDataArrays()[1].getName(),"charge")
		TEST_EQUAL(spec.getMetaDataArrays()[1].size(),10)
		
		//precursor
		TEST_REAL_EQUAL(spec.getPrecursorPeak().getIntensity(),120053)
		TEST_EQUAL(spec.getPrecursorPeak().getCharge(),2)
		TEST_REAL_EQUAL(spec.getPrecursorPeak().getPosition()[0],445.34)
		TEST_EQUAL(spec.getPrecursor().getActivationMethod(),Precursor::CID)
		TEST_REAL_EQUAL(spec.getPrecursor().getActivationEnergy(),35)
	}
	
	//-------------------------- spectrum 2 --------------------------
	{
		const MSSpectrum<>& spec = exp[2]; 
		TEST_EQUAL(spec.size(),0)
		TEST_EQUAL(spec.getMSLevel(),1)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::PRODUCT)
		TEST_EQUAL(spec.getMetaDataArrays().size(),0)
		TEST_EQUAL(spec.getType(),SpectrumSettings::UNKNOWN)
		TEST_REAL_EQUAL(spec.getRT(),-1.0)
		TEST_EQUAL(spec.getInstrumentSettings().getPolarity(),IonSource::POLNULL)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getMzRangeStart(),0.0)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getMzRangeStop(),0.0)
		TEST_STRING_EQUAL(spec.getAcquisitionInfo().getMethodOfCombination(),"")
		TEST_EQUAL(spec.getAcquisitionInfo().size(),0)
		TEST_STRING_EQUAL(spec.getSourceFile().getNameOfFile(),"tiny1.dta")
		TEST_STRING_EQUAL(spec.getSourceFile().getPathToFile(),"file:///F:/data/Exp01")
		TEST_STRING_EQUAL(spec.getSourceFile().getSha1(),"81be39fb2700ab2f3c8b2234b91274968b6899b1")
	}

RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

