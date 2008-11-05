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

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

DRange<1> makeRange(Real a, Real b)
{
	DPosition<1> pa(a), pb(b);
	return DRange<1>(pa, pb);
}

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
  		
	TEST_EQUAL(exp.size(),4)
	//id
	TEST_EQUAL(exp.getIdentifier(),"urn:lsid:psidev.info:mzML.instanceDocuments.tiny.pwiz");
	//contacts
	TEST_EQUAL(exp.getContacts().size(),2)
	TEST_STRING_EQUAL(exp.getContacts()[0].getFirstName(),"William")
	TEST_STRING_EQUAL(exp.getContacts()[0].getLastName(),"Pennington")
	TEST_STRING_EQUAL(exp.getContacts()[0].getEmail(),"wpennington@higglesworth.edu")
	TEST_STRING_EQUAL(exp.getContacts()[0].getURL(),"http://www.higglesworth.edu/")
	TEST_STRING_EQUAL(exp.getContacts()[0].getAddress(),"Higglesworth University, 12 Higglesworth Avenue, 12045, HI, USA")
	TEST_STRING_EQUAL(exp.getContacts()[1].getFirstName(),"")
	TEST_STRING_EQUAL(exp.getContacts()[1].getLastName(),"Drek'Thar")
	TEST_STRING_EQUAL(exp.getContacts()[1].getEmail(),"")
	TEST_STRING_EQUAL(exp.getContacts()[1].getURL(),"")
	TEST_STRING_EQUAL(exp.getContacts()[1].getAddress(),"")
	//source files
	TEST_EQUAL(exp.getSourceFiles().size(),2);
	TEST_STRING_EQUAL(exp.getSourceFiles()[0].getNameOfFile(),"tiny1.RAW")
	TEST_STRING_EQUAL(exp.getSourceFiles()[0].getPathToFile(),"file:///F:/data/Exp01")
	TEST_STRING_EQUAL(exp.getSourceFiles()[0].getChecksum(),"71be39fb2700ab2f3c8b2234b91274968b6899b1")
	TEST_EQUAL(exp.getSourceFiles()[0].getChecksumType(),SourceFile::SHA1)
	TEST_STRING_EQUAL(exp.getSourceFiles()[0].getFileType(),"Xcalibur RAW")
	TEST_STRING_EQUAL(exp.getSourceFiles()[1].getNameOfFile(),"tiny2.RAW")
	TEST_STRING_EQUAL(exp.getSourceFiles()[1].getPathToFile(),"file:///F:/data/Exp02")
	TEST_STRING_EQUAL(exp.getSourceFiles()[1].getChecksum(),"71be39fb2700ab2f3c8b2234b91274968b6899b2")
	TEST_EQUAL(exp.getSourceFiles()[1].getChecksumType(),SourceFile::MD5)
	TEST_STRING_EQUAL(exp.getSourceFiles()[1].getFileType(),"pkl")
	//sample
	TEST_STRING_EQUAL(exp.getSample().getName(),"Sample1")
	TEST_REAL_EQUAL(exp.getSample().getMass(),11.7)
	TEST_STRING_EQUAL(exp.getSample().getNumber(),"5")
	TEST_REAL_EQUAL(exp.getSample().getVolume(),3.1)
	TEST_REAL_EQUAL(exp.getSample().getConcentration(),5.5)
	//instrument (general)
	TEST_STRING_EQUAL(exp.getInstrument().getName(),"LCQ Deca")
	TEST_STRING_EQUAL(exp.getInstrument().getCustomizations(),"Umbau")
	//ion sources
	TEST_EQUAL(exp.getInstrument().getIonSources().size(),2)
	TEST_EQUAL(exp.getInstrument().getIonSources()[0].getOrder(),101)
	TEST_EQUAL(exp.getInstrument().getIonSources()[0].getInletType(),IonSource::DIRECT)
	TEST_EQUAL(exp.getInstrument().getIonSources()[0].getIonizationMethod(),IonSource::ESI)
	TEST_EQUAL(exp.getInstrument().getIonSources()[1].getOrder(),102)
	TEST_EQUAL(exp.getInstrument().getIonSources()[1].getInletType(),IonSource::DIRECT)
	TEST_EQUAL(exp.getInstrument().getIonSources()[1].getIonizationMethod(),IonSource::FAB)
	//mass analyzers
	TEST_EQUAL(exp.getInstrument().getMassAnalyzers().size(),2)
	TEST_EQUAL(exp.getInstrument().getMassAnalyzers()[0].getOrder(),201)
	TEST_EQUAL(exp.getInstrument().getMassAnalyzers()[0].getType(),MassAnalyzer::PAULIONTRAP)
	TEST_REAL_EQUAL(exp.getInstrument().getMassAnalyzers()[0].getMagneticFieldStrength(),14.56)
	TEST_EQUAL(exp.getInstrument().getMassAnalyzers()[1].getOrder(),202)
	TEST_EQUAL(exp.getInstrument().getMassAnalyzers()[1].getType(),MassAnalyzer::LIT)
	TEST_REAL_EQUAL(exp.getInstrument().getMassAnalyzers()[1].getMagneticFieldStrength(),1414.14)
	//detectors
	TEST_EQUAL(exp.getInstrument().getIonDetectors().size(),2)
	TEST_EQUAL(exp.getInstrument().getIonDetectors()[0].getOrder(),301)
	TEST_EQUAL(exp.getInstrument().getIonDetectors()[0].getType(),IonDetector::ELECTRONMULTIPLIER)
	TEST_EQUAL(exp.getInstrument().getIonDetectors()[0].getAcquisitionMode(),IonDetector::TDC)
	TEST_REAL_EQUAL(exp.getInstrument().getIonDetectors()[0].getResolution(),5.1)
	TEST_REAL_EQUAL(exp.getInstrument().getIonDetectors()[0].getADCSamplingFrequency(),1.1)
	TEST_EQUAL(exp.getInstrument().getIonDetectors()[1].getOrder(),302)
	TEST_EQUAL(exp.getInstrument().getIonDetectors()[1].getType(),IonDetector::ELECTRONMULTIPLIER)
	TEST_EQUAL(exp.getInstrument().getIonDetectors()[1].getAcquisitionMode(),IonDetector::TDC)
	TEST_REAL_EQUAL(exp.getInstrument().getIonDetectors()[1].getResolution(),6.1)
	TEST_REAL_EQUAL(exp.getInstrument().getIonDetectors()[1].getADCSamplingFrequency(),1.1)
	//instrument software
	TEST_EQUAL(exp.getInstrument().getSoftware().getName(),"Bioworks")
	TEST_EQUAL(exp.getInstrument().getSoftware().getVersion(),"3.3.1 sp1")
	//processing
	TEST_EQUAL(exp.getDataProcessing().size(),2)
	TEST_EQUAL(exp.getDataProcessing()[1].getSoftware().getName(),"ProteoWizard")
	TEST_EQUAL(exp.getDataProcessing()[1].getSoftware().getVersion(),"1.0")
	TEST_EQUAL(exp.getDataProcessing()[1].getProcessingActions().size(),1)
	TEST_EQUAL(exp.getDataProcessing()[1].getProcessingActions().count(DataProcessing::CONVERSION_MZML),1)
	TEST_EQUAL(exp.getDataProcessing()[1].isMetaEmpty(),true)
	TEST_EQUAL(exp.getDataProcessing()[0].getSoftware().getName(),"Xcalibur")
	TEST_EQUAL(exp.getDataProcessing()[0].getSoftware().getVersion(),"2.0.5")
	TEST_EQUAL(exp.getDataProcessing()[0].getProcessingActions().size(),4)
	TEST_EQUAL(exp.getDataProcessing()[0].getProcessingActions().count(DataProcessing::DEISOTOPING),1)
	TEST_EQUAL(exp.getDataProcessing()[0].getProcessingActions().count(DataProcessing::CHARGE_DECONVOLUTION),1)
	TEST_EQUAL(exp.getDataProcessing()[0].getProcessingActions().count(DataProcessing::LOW_INTENSITY_REMOVAL),1)
	TEST_REAL_EQUAL(DoubleReal(exp.getDataProcessing()[0].getMetaValue("low_intensity_threshold")),5.9)
	TEST_REAL_EQUAL(DoubleReal(exp.getDataProcessing()[0].getMetaValue("high_intensity_threshold")),10.9)
	TEST_EQUAL(exp.getDataProcessing()[0].isMetaEmpty(),false)

	//-------------------------- spectrum 0 --------------------------
	{
		const MSSpectrum<>& spec = exp[0]; 
		//peaks
		TEST_EQUAL(spec.size(),15)
		for (UInt i=0; i<15; ++i)
		{
			TEST_REAL_EQUAL(spec[i].getMZ(),i);
			TEST_REAL_EQUAL(spec[i].getIntensity(),15-i);
		}
		//general info		
		TEST_EQUAL(spec.getMSLevel(),1)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::PRODUCT)
		TEST_EQUAL(spec.getMetaDataArrays().size(),0)
		TEST_EQUAL(spec.getType(),SpectrumSettings::PEAKS)
		TEST_REAL_EQUAL(spec.getRT(),5.1)
		TEST_EQUAL(spec.getInstrumentSettings().getScanWindows().size(),1)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getScanWindows()[0].begin,400.0)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getScanWindows()[0].end,1800.0)
		TEST_STRING_EQUAL(spec.getAcquisitionInfo().getMethodOfCombination(),"median")
		TEST_EQUAL(spec.getAcquisitionInfo().size(),2)
		TEST_EQUAL(spec.getAcquisitionInfo()[0].getNumber(),4711)
		TEST_EQUAL(spec.getAcquisitionInfo()[1].getNumber(),4712)
		TEST_EQUAL(spec.getSourceFile()==SourceFile(),true)
	}
	
	//-------------------------- spectrum 1 --------------------------
	{
		const MSSpectrum<>& spec = exp[1]; 	
		//peaks
		TEST_EQUAL(spec.size(),10)
		for (UInt i=0; i<10; ++i)
		{
			TEST_REAL_EQUAL(spec[i].getMZ(),2.0*i);
			TEST_REAL_EQUAL(spec[i].getIntensity(),20-2.0*i);
		}
		//general info
		TEST_EQUAL(spec.getMSLevel(),2)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::PRODUCT)
		TEST_EQUAL(spec.getType(),SpectrumSettings::PEAKS)
		TEST_REAL_EQUAL(spec.getRT(),5.2)
		TEST_EQUAL(spec.getInstrumentSettings().getPolarity(),IonSource::POSITIVE)
		TEST_EQUAL(spec.getInstrumentSettings().getScanWindows().size(),3)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getScanWindows()[0].begin,100.0)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getScanWindows()[0].end,500.0)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getScanWindows()[1].begin,600.0)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getScanWindows()[1].end,1000.0)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getScanWindows()[2].begin,1100.0)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getScanWindows()[2].end,1500.0)
		TEST_EQUAL(spec.getAcquisitionInfo().getMethodOfCombination(),"")
		TEST_EQUAL(spec.getAcquisitionInfo().size(),0)
		//meta data arrays
		TEST_EQUAL(spec.getMetaDataArrays().size(),2)
		TEST_STRING_EQUAL(spec.getMetaDataArrays()[0].getName(),"signal to noise")
		TEST_EQUAL(spec.getMetaDataArrays()[0].size(),10)
		TEST_STRING_EQUAL(spec.getMetaDataArrays()[1].getName(),"charge")
		TEST_EQUAL(spec.getMetaDataArrays()[1].size(),10)
		//precursor
		TEST_REAL_EQUAL(spec.getPrecursorPeak().getIntensity(),120053)
		TEST_EQUAL(spec.getPrecursorPeak().getCharge(),2)
		TEST_REAL_EQUAL(spec.getPrecursorPeak().getPosition()[0],5.55)
		TEST_EQUAL(spec.getPrecursor().getActivationMethod(),Precursor::CID)
		TEST_REAL_EQUAL(spec.getPrecursor().getActivationEnergy(),35)
		TEST_EQUAL(spec.getPrecursorPeak().getPossibleChargeStates().size(),3)
		TEST_EQUAL(spec.getPrecursorPeak().getPossibleChargeStates()[0],1)
		TEST_EQUAL(spec.getPrecursorPeak().getPossibleChargeStates()[1],3)
		TEST_EQUAL(spec.getPrecursorPeak().getPossibleChargeStates()[2],4)
		//source file
		TEST_STRING_EQUAL(spec.getSourceFile().getNameOfFile(),"tiny1.dta")
		TEST_STRING_EQUAL(spec.getSourceFile().getPathToFile(),"file:///F:/data/Exp01")
		TEST_STRING_EQUAL(spec.getSourceFile().getChecksum(),"81be39fb2700ab2f3c8b2234b91274968b6899b1")
		TEST_EQUAL(spec.getSourceFile().getChecksumType(),SourceFile::SHA1)
		
	}
	
	//-------------------------- spectrum 2 --------------------------
	{
		const MSSpectrum<>& spec = exp[2]; 
		//peaks
		TEST_EQUAL(spec.size(),15)
		for (UInt i=0; i<15; ++i)
		{
			TEST_REAL_EQUAL(spec[i].getMZ(),i);
			TEST_REAL_EQUAL(spec[i].getIntensity(),15-i);
		}
		//general info		
		TEST_EQUAL(spec.getMSLevel(),1)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::PRODUCT)
		TEST_EQUAL(spec.getMetaDataArrays().size(),0)
		TEST_EQUAL(spec.getType(),SpectrumSettings::PEAKS)
		TEST_REAL_EQUAL(spec.getRT(),5.3)
		TEST_EQUAL(spec.getInstrumentSettings().getPolarity(),IonSource::POSITIVE)
		TEST_EQUAL(spec.getInstrumentSettings().getScanWindows().size(),1)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getScanWindows()[0].begin,400.0)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getScanWindows()[0].end,1800.0)
		TEST_STRING_EQUAL(spec.getAcquisitionInfo().getMethodOfCombination(),"median")
		TEST_EQUAL(spec.getAcquisitionInfo().size(),2)
		TEST_EQUAL(spec.getAcquisitionInfo()[0].getNumber(),4711)
		TEST_EQUAL(spec.getAcquisitionInfo()[1].getNumber(),4712)
		TEST_EQUAL(spec.getSourceFile()==SourceFile(),true)
	}

	//-------------------------- spectrum 3 (no peaks) --------------------------
	{
		const MSSpectrum<>& spec = exp[3]; 
		//peaks
		TEST_EQUAL(spec.size(),0)
		//general info
		TEST_EQUAL(spec.getMSLevel(),1)
		TEST_REAL_EQUAL(spec.getRT(),5.4)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::PRODUCT)
		TEST_EQUAL(spec.getMetaDataArrays().size(),0)
		TEST_EQUAL(spec.getType(),SpectrumSettings::RAWDATA)
		TEST_EQUAL(spec.getInstrumentSettings().getScanWindows().size(),1)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getScanWindows()[0].begin,110.0)
		TEST_REAL_EQUAL(spec.getInstrumentSettings().getScanWindows()[0].end,905.0)
		TEST_STRING_EQUAL(spec.getAcquisitionInfo().getMethodOfCombination(),"")
		TEST_EQUAL(spec.getAcquisitionInfo().size(),0)
	}
	
	//-------------------------- userParam --------------------------
	//run
  TEST_EQUAL(exp.getMetaValue("flag").valueType(),DataValue::STRING_VALUE)
  TEST_STRING_EQUAL((String)exp.getMetaValue("flag"),"")
  TEST_EQUAL(exp.getMetaValue("string").valueType(),DataValue::STRING_VALUE)
  TEST_STRING_EQUAL((String)exp.getMetaValue("string"),"bla")
  TEST_EQUAL(exp.getMetaValue("float").valueType(),DataValue::DOUBLE_VALUE)
  TEST_REAL_EQUAL((double)exp.getMetaValue("float"),5.11)
  TEST_EQUAL(exp.getMetaValue("int").valueType(),DataValue::INT_VALUE)
  TEST_EQUAL((Int)exp.getMetaValue("int"),5)
	//instrumentConfiguration
	TEST_STRING_EQUAL((String)exp.getInstrument().getMetaValue("name"),"instrumentConfiguration")
	TEST_STRING_EQUAL((String)exp.getInstrument().getIonSources()[0].getMetaValue("name"),"source1")
	TEST_STRING_EQUAL((String)exp.getInstrument().getIonSources()[1].getMetaValue("name"),"source2")
	TEST_STRING_EQUAL((String)exp.getInstrument().getMassAnalyzers()[0].getMetaValue("name"),"analyzer1")
	TEST_STRING_EQUAL((String)exp.getInstrument().getMassAnalyzers()[1].getMetaValue("name"),"analyzer2")
	TEST_STRING_EQUAL((String)exp.getInstrument().getIonDetectors()[0].getMetaValue("name"),"detector1")
	TEST_STRING_EQUAL((String)exp.getInstrument().getIonDetectors()[1].getMetaValue("name"),"detector2")
	//sample
	TEST_STRING_EQUAL((String)exp.getSample().getMetaValue("name"),"sample")
	//contact
	TEST_STRING_EQUAL((String)exp.getContacts()[0].getMetaValue("name"),"contact1")
	TEST_STRING_EQUAL((String)exp.getContacts()[1].getMetaValue("name"),"contact2")
	//spectrum
	TEST_STRING_EQUAL((String)exp[0].getMetaValue("name"),"spectrum1")
	TEST_STRING_EQUAL((String)exp[1].getMetaValue("name"),"spectrum2")
	TEST_STRING_EQUAL((String)exp[2].getMetaValue("name"),"spectrum3")
	TEST_STRING_EQUAL((String)exp[3].getMetaValue("name"),"spectrum4")
	TEST_STRING_EQUAL((String)exp[0].getMetaValue("sdname"),"spectrumdescription1")
	TEST_STRING_EQUAL((String)exp[1].getMetaValue("sdname"),"spectrumdescription2")
	TEST_STRING_EQUAL((String)exp[2].getMetaValue("sdname"),"spectrumdescription3")
	TEST_STRING_EQUAL((String)exp[3].getMetaValue("sdname"),"spectrumdescription4")
	TEST_STRING_EQUAL((String)exp[0].getMetaValue("mzname"),"mzarray1")
	TEST_STRING_EQUAL((String)exp[0].getMetaValue("itname"),"itarray1")
	TEST_STRING_EQUAL((String)exp[1].getMetaValue("mzname"),"mzarray2")
	TEST_STRING_EQUAL((String)exp[1].getMetaValue("itname"),"itarray2")
	//binaryDataArray
	TEST_STRING_EQUAL((String)exp[1].getMetaDataArrays()[0].getMetaValue("name"),"binaryDataArray_sn")
	TEST_STRING_EQUAL((String)exp[1].getMetaDataArrays()[0].getMetaValue("name2"),"binaryDataArray_sn2")
	TEST_STRING_EQUAL((String)exp[1].getMetaDataArrays()[1].getMetaValue("name"),"binaryDataArray_c")
	TEST_STRING_EQUAL((String)exp[1].getMetaDataArrays()[1].getMetaValue("name2"),"")
	//scan
	TEST_STRING_EQUAL((String)exp[0].getInstrumentSettings().getMetaValue("name"),"scan1")
	TEST_STRING_EQUAL((String)exp[1].getInstrumentSettings().getMetaValue("name"),"scan2")
	TEST_STRING_EQUAL((String)exp[2].getInstrumentSettings().getMetaValue("name"),"scan3")
	TEST_STRING_EQUAL((String)exp[3].getInstrumentSettings().getMetaValue("name"),"")
	//acquisition list
	TEST_STRING_EQUAL((String)exp[0].getAcquisitionInfo().getMetaValue("name"),"acquisition_list")	
	//acquisition
	TEST_STRING_EQUAL((String)exp[0].getAcquisitionInfo()[0].getMetaValue("name"),"acquisition1")
	TEST_STRING_EQUAL((String)exp[0].getAcquisitionInfo()[1].getMetaValue("name"),"acquisition2")
	//source file
	TEST_STRING_EQUAL((String)exp.getSourceFiles()[0].getMetaValue("name"),"sourcefile1")
	TEST_STRING_EQUAL((String)exp.getSourceFiles()[1].getMetaValue("name"),"sourcefile2")
	TEST_STRING_EQUAL((String)exp[1].getSourceFile().getMetaValue("name"),"sourcefile4")
	//data processing
	TEST_STRING_EQUAL(exp.getDataProcessing()[0].getMetaValue("p1").toString(),"value1")
	TEST_STRING_EQUAL(exp.getDataProcessing()[0].getMetaValue("p2").toString(),"value2")
	//precursor
	TEST_STRING_EQUAL(exp[1].getPrecursor().getMetaValue("iwname").toString(),"isolationwindow1")
	TEST_STRING_EQUAL(exp[1].getPrecursor().getMetaValue("siname").toString(),"selectedion1")
	TEST_STRING_EQUAL(exp[1].getPrecursor().getMetaValue("acname").toString(),"activation1")

	//-------------------------- cvParam (but no member => meta data)--------------------------
	//spectrum 1
	TEST_STRING_EQUAL((String)exp[1].getMetaValue("mass resolution"),"4.1")	
	TEST_STRING_EQUAL((String)exp[1].getPrecursor().getMetaValue("collision energy"),"4.2")	
	TEST_STRING_EQUAL((String)exp[0].getAcquisitionInfo()[0].getMetaValue("mass resolution"),"4.3")
	TEST_STRING_EQUAL((String)exp.getSample().getMetaValue("sample batch"),"4.4")
	TEST_STRING_EQUAL((String)exp.getInstrument().getMetaValue("ion optics"),"magnetic deflection")

	/////////////////////// TESTING SPECIAL CASES ///////////////////////
	
	//load a second time to make sure everything is re-initialized correctly
	MSExperiment<> exp2;
	file.load("data/MzMLFile_1.mzML",exp2);
	TEST_EQUAL(exp==exp2,true)
	
	//load minimal file to see if that works
	MSExperiment<> exp3;
	file.load("data/MzMLFile_2_minimal.mzML",exp3);
	TEST_EQUAL(exp3.size(),0)
RESULT


CHECK([EXTRA] load only meta data)
	MzMLFile file;
	file.getOptions().setMetadataOnly(true);
	MSExperiment<> exp;
	file.load("data/MzMLFile_1.mzML",exp);

	TEST_EQUAL(exp.size(),0)
	TEST_EQUAL(exp.getIdentifier(),"urn:lsid:psidev.info:mzML.instanceDocuments.tiny.pwiz");
	TEST_EQUAL(exp.getContacts().size(),2)
	TEST_EQUAL(exp.getSourceFiles().size(),2);
	TEST_EQUAL(exp.getInstrument().getMassAnalyzers().size(),2)
	TEST_EQUAL(exp.getDataProcessing().size(),2)
RESULT


CHECK([EXTRA] load with restricted MS levels)
	MzMLFile file;
	file.getOptions().addMSLevel(1);
	MSExperiment<> exp;
	file.load("data/MzMLFile_1.mzML",exp);

	TEST_EQUAL(exp.size(),3)
	TEST_REAL_EQUAL(exp[0].getRT(),5.1)
	TEST_EQUAL((Int)exp[0].getMetaValue("original_spectrum_number"),0)
	TEST_REAL_EQUAL(exp[1].getRT(),5.3)
	TEST_EQUAL((Int)exp[1].getMetaValue("original_spectrum_number"),2)
	TEST_REAL_EQUAL(exp[2].getRT(),5.4)
	TEST_EQUAL((Int)exp[2].getMetaValue("original_spectrum_number"),3)
RESULT


CHECK([EXTRA] load with restricted RT range)
	MzMLFile file;
	file.getOptions().setRTRange(makeRange(5.15,5.35));
	MSExperiment<> exp;
	file.load("data/MzMLFile_1.mzML",exp);
	TEST_EQUAL(exp.size(),2)
	TEST_REAL_EQUAL(exp[0].getRT(),5.2)
	TEST_REAL_EQUAL(exp[1].getRT(),5.3)
RESULT

CHECK([EXTRA] load with restricted m/z range)
	MzMLFile file;
	file.getOptions().setMZRange(makeRange(6.5,9.5));
	MSExperiment<> exp;
	file.load("data/MzMLFile_1.mzML",exp);
	
	TEST_EQUAL(exp.size(),4)
	TEST_EQUAL(exp[0].size(),3)
	TEST_REAL_EQUAL(exp[0][0].getMZ(),7.0)
	TEST_REAL_EQUAL(exp[0][1].getMZ(),8.0)
	TEST_REAL_EQUAL(exp[0][2].getMZ(),9.0)
	TEST_EQUAL(exp[1].size(),1)
	TEST_REAL_EQUAL(exp[1][0].getMZ(),8.0)
	TEST_EQUAL(exp[2].size(),3)
	TEST_REAL_EQUAL(exp[2][0].getMZ(),7.0)
	TEST_REAL_EQUAL(exp[2][1].getMZ(),8.0)
	TEST_REAL_EQUAL(exp[2][2].getMZ(),9.0)
	TEST_EQUAL(exp[3].size(),0)
RESULT



CHECK([EXTRA] load intensity range)
	MzMLFile file;
	file.getOptions().setIntensityRange(makeRange(6.5,9.5));
	MSExperiment<> exp;
	file.load("data/MzMLFile_1.mzML",exp);

	TEST_EQUAL(exp.size(),4)
	TEST_EQUAL(exp[0].size(),3)
	TEST_REAL_EQUAL(exp[0][0].getIntensity(),9.0)
	TEST_REAL_EQUAL(exp[0][1].getIntensity(),8.0)
	TEST_REAL_EQUAL(exp[0][2].getIntensity(),7.0)
	TEST_EQUAL(exp[1].size(),1)
	TEST_REAL_EQUAL(exp[1][0].getIntensity(),8.0)
	TEST_EQUAL(exp[2].size(),3)
	TEST_REAL_EQUAL(exp[2][0].getIntensity(),9.0)
	TEST_REAL_EQUAL(exp[2][1].getIntensity(),8.0)
	TEST_REAL_EQUAL(exp[2][2].getIntensity(),7.0)
	TEST_EQUAL(exp[3].size(),0)
RESULT


CHECK((template <typename MapType> void store(const String& filename, const MapType& map) const))
	std::string tmp_filename;
 	NEW_TMP_FILE(tmp_filename);
 	
	MzMLFile file;
	MSExperiment<> exp;
	TEST_EXCEPTION(Exception::NotImplemented, file.store(tmp_filename,exp))
RESULT


CHECK([EXTRA] bool isValid(const String& filename))
	MzMLFile file;
	TEST_EQUAL(file.isValid("data/MzMLFile_1.mzML"),true)
	TEST_EQUAL(file.isValid("data/MzMLFile_2_minimal.mzML"),true)
	TEST_EQUAL(file.isValid("data/MzMLFile_4_indexed.mzML"),true)
RESULT

CHECK( bool isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings))
	MzMLFile file;
	StringList errors, warnings;
	
	//valid file
	TEST_EQUAL(file.isSemanticallyValid("data/MzMLFile_1.mzML", errors, warnings),false)
	TEST_EQUAL(errors.size(),4)
	TEST_EQUAL(warnings.size(),2)
	
	//invalid file
	TEST_EQUAL(file.isSemanticallyValid("data/MzMLFile_3_invalid.mzML", errors, warnings),false)
	
	TEST_EQUAL(errors.size(),7)
	TEST_STRING_EQUAL(errors[0],"CV term should not have a value: 'MS:1000580 - MSn spectrum' (value: '4444') at element '/mzML/run/spectrumList/spectrum'")
	TEST_STRING_EQUAL(errors[1],"CV term should have a floating-point value: 'MS:1000528 - lowest m/z value' (value: 'abc') at element '/mzML/run/spectrumList/spectrum/spectrumDescription'")
	TEST_STRING_EQUAL(errors[2],"CV term should have a numerical value: 'MS:1000527 - highest m/z value' (value: '') at element '/mzML/run/spectrumList/spectrum/spectrumDescription'")
	TEST_STRING_EQUAL(errors[3],"CV term used in invalid element: 'MS:1000133 - collision-induced dissociation' at element '/mzML/run/spectrumList/spectrum/spectrumDescription/precursorList/precursor/activation'")
	TEST_STRING_EQUAL(errors[4],"CV term used in invalid element: 'MS:1000509 - activation energy' at element '/mzML/run/spectrumList/spectrum/spectrumDescription/precursorList/precursor/activation'")
	TEST_STRING_EQUAL(errors[5],"CV term used in invalid element: 'MS:1000045 - collision energy' at element '/mzML/run/spectrumList/spectrum/spectrumDescription/precursorList/precursor/activation'")
	TEST_STRING_EQUAL(errors[6],"Violated mapping rule 'R23' at element '/mzML/run/spectrumList/spectrum/spectrumDescription/precursorList/precursor/activation'")
	
	TEST_EQUAL(warnings.size(),2)
	TEST_STRING_EQUAL(warnings[0],"No mapping rule found for element '/mzML/run/chromatogramList/chromatogram'")
	TEST_STRING_EQUAL(warnings[1],"No mapping rule found for element '/mzML/run/chromatogramList/chromatogram'")

	//indexed MzML
	TEST_EQUAL(file.isSemanticallyValid("data/MzMLFile_4_indexed.mzML", errors, warnings),false)
	TEST_EQUAL(errors.size(),7)
	TEST_EQUAL(warnings.size(),2)

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

