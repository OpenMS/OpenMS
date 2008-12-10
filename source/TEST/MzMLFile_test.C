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
START_SECTION((MzMLFile()))
	ptr = new MzMLFile;
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~MzMLFile()))
	delete ptr;
END_SECTION

START_SECTION(const PeakFileOptions& getOptions() const)
	MzMLFile file;
	TEST_EQUAL(file.getOptions().hasMSLevels(),false)
END_SECTION

START_SECTION(PeakFileOptions& getOptions())
	MzMLFile file;
	file.getOptions().addMSLevel(1);
	TEST_EQUAL(file.getOptions().hasMSLevels(),true);
END_SECTION

TOLERANCE_ABSOLUTE(0.01)

START_SECTION((template <typename MapType> void load(const String& filename, MapType& map)))
	MzMLFile file;
	MSExperiment<> exp;
	file.load("data/MzMLFile_1.mzML",exp);

	//-------------------------- general information --------------------------
  		
	TEST_EQUAL(exp.size(),4)
	//run
	TEST_EQUAL(exp.getNativeIDType(),ExperimentalSettings::MULTIPLE_PEAK_LISTS)
	TEST_EQUAL(exp.getIdentifier(),"document_accession")
	TEST_EQUAL(exp.getDateTime().get(),"2007-06-27 15:23:45")
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
	TEST_STRING_EQUAL(exp.getSourceFiles()[0].getFileType(),"Xcalibur RAW file")
	TEST_STRING_EQUAL(exp.getSourceFiles()[1].getNameOfFile(),"tiny2.RAW")
	TEST_STRING_EQUAL(exp.getSourceFiles()[1].getPathToFile(),"file:///F:/data/Exp02")
	TEST_STRING_EQUAL(exp.getSourceFiles()[1].getChecksum(),"71be39fb2700ab2f3c8b2234b91274968b6899b2")
	TEST_EQUAL(exp.getSourceFiles()[1].getChecksumType(),SourceFile::MD5)
	TEST_STRING_EQUAL(exp.getSourceFiles()[1].getFileType(),"pkl file")
	//sample
	TEST_STRING_EQUAL(exp.getSample().getName(),"Sample1")
	TEST_REAL_SIMILAR(exp.getSample().getMass(),11.7)
	TEST_STRING_EQUAL(exp.getSample().getNumber(),"5")
	TEST_REAL_SIMILAR(exp.getSample().getVolume(),3.1)
	TEST_REAL_SIMILAR(exp.getSample().getConcentration(),5.5)
	TEST_EQUAL(exp.getSample().getState(),Sample::SUSPENSION)
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
	TEST_REAL_SIMILAR(exp.getInstrument().getMassAnalyzers()[0].getMagneticFieldStrength(),14.56)
	TEST_EQUAL(exp.getInstrument().getMassAnalyzers()[1].getOrder(),202)
	TEST_EQUAL(exp.getInstrument().getMassAnalyzers()[1].getType(),MassAnalyzer::LIT)
	TEST_REAL_SIMILAR(exp.getInstrument().getMassAnalyzers()[1].getMagneticFieldStrength(),1414.14)
	//detectors
	TEST_EQUAL(exp.getInstrument().getIonDetectors().size(),2)
	TEST_EQUAL(exp.getInstrument().getIonDetectors()[0].getOrder(),301)
	TEST_EQUAL(exp.getInstrument().getIonDetectors()[0].getType(),IonDetector::ELECTRONMULTIPLIER)
	TEST_EQUAL(exp.getInstrument().getIonDetectors()[0].getAcquisitionMode(),IonDetector::TDC)
	TEST_REAL_SIMILAR(exp.getInstrument().getIonDetectors()[0].getResolution(),5.1)
	TEST_REAL_SIMILAR(exp.getInstrument().getIonDetectors()[0].getADCSamplingFrequency(),1.1)
	TEST_EQUAL(exp.getInstrument().getIonDetectors()[1].getOrder(),302)
	TEST_EQUAL(exp.getInstrument().getIonDetectors()[1].getType(),IonDetector::ELECTRONMULTIPLIER)
	TEST_EQUAL(exp.getInstrument().getIonDetectors()[1].getAcquisitionMode(),IonDetector::TDC)
	TEST_REAL_SIMILAR(exp.getInstrument().getIonDetectors()[1].getResolution(),6.1)
	TEST_REAL_SIMILAR(exp.getInstrument().getIonDetectors()[1].getADCSamplingFrequency(),1.1)
	//instrument software
	TEST_EQUAL(exp.getInstrument().getSoftware().getName(),"Bioworks")
	TEST_EQUAL(exp.getInstrument().getSoftware().getVersion(),"3.3.1 sp1")
	//processing
	TEST_EQUAL(exp.getDataProcessing().size(),2)
	TEST_EQUAL(exp.getDataProcessing()[0].getSoftware().getName(),"Xcalibur")
	TEST_EQUAL(exp.getDataProcessing()[0].getSoftware().getVersion(),"2.0.5")
	TEST_EQUAL(exp.getDataProcessing()[0].getProcessingActions().size(),2)
	TEST_EQUAL(exp.getDataProcessing()[0].getProcessingActions().count(DataProcessing::DEISOTOPING),1)
	TEST_EQUAL(exp.getDataProcessing()[0].getProcessingActions().count(DataProcessing::CHARGE_DECONVOLUTION),1)
	TEST_STRING_EQUAL(exp.getDataProcessing()[0].getCompletionTime().get(),"2001-02-03 04:05:00")
	TEST_REAL_SIMILAR(DoubleReal(exp.getDataProcessing()[0].getMetaValue("low_intensity_threshold")),5.9)
	TEST_REAL_SIMILAR(DoubleReal(exp.getDataProcessing()[0].getMetaValue("high_intensity_threshold")),10.9)
	TEST_EQUAL(exp.getDataProcessing()[0].isMetaEmpty(),false)
	TEST_EQUAL(exp.getDataProcessing()[1].getSoftware().getName(),"ProteoWizard")
	TEST_EQUAL(exp.getDataProcessing()[1].getSoftware().getVersion(),"1.0")
	TEST_EQUAL(exp.getDataProcessing()[1].getProcessingActions().size(),1)
	TEST_EQUAL(exp.getDataProcessing()[1].getProcessingActions().count(DataProcessing::CONVERSION_MZML),1)
	TEST_EQUAL(exp.getDataProcessing()[1].isMetaEmpty(),false)

	//-------------------------- spectrum 0 --------------------------
	{
		const MSSpectrum<>& spec = exp[0]; 
		//peaks
		TEST_EQUAL(spec.size(),15)
		for (UInt i=0; i<15; ++i)
		{
			TEST_REAL_SIMILAR(spec[i].getMZ(),i);
			TEST_REAL_SIMILAR(spec[i].getIntensity(),15-i);
		}
		//general info		
		TEST_EQUAL(spec.getMSLevel(),1)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::FULL)
		TEST_EQUAL(spec.getMetaDataArrays().size(),0)
		TEST_EQUAL(spec.getType(),SpectrumSettings::PEAKS)
		TEST_REAL_SIMILAR(spec.getRT(),5.1)
		TEST_EQUAL(spec.getInstrumentSettings().getScanWindows().size(),1)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[0].begin,400.0)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[0].end,1800.0)
		TEST_STRING_EQUAL(spec.getAcquisitionInfo().getMethodOfCombination(),"median of spectra")
		TEST_EQUAL(spec.getAcquisitionInfo().size(),2)
		TEST_EQUAL(spec.getAcquisitionInfo()[0].getNumber(),4711)
		TEST_STRING_EQUAL(spec.getAcquisitionInfo()[0].getMetaValue("source_file_name"),"ac.dta")
		TEST_STRING_EQUAL(spec.getAcquisitionInfo()[0].getMetaValue("source_file_path"),"file:///F:/data/Exp02")
//		TEST_STRING_EQUAL(spec.getAcquisitionInfo()[0].getMetaValue("external_native_id"),"ENI0")
//		TEST_STRING_EQUAL(spec.getAcquisitionInfo()[0].getMetaValue("external_spectrum_id"),"ESI0")
		TEST_EQUAL(spec.getAcquisitionInfo()[1].getNumber(),4712)
		TEST_EQUAL(spec.getSourceFile()==SourceFile(),true)
		//ids
		TEST_STRING_EQUAL(spec.getNativeID(),"index=0")
		TEST_STRING_EQUAL(spec.getMetaValue("maldi_spot_id"),"M0")
	}
	
	//-------------------------- spectrum 1 --------------------------
	{
		const MSSpectrum<>& spec = exp[1]; 	
		//peaks
		TEST_EQUAL(spec.size(),10)
		for (UInt i=0; i<10; ++i)
		{
			TEST_REAL_SIMILAR(spec[i].getMZ(),2.0*i);
			TEST_REAL_SIMILAR(spec[i].getIntensity(),20-2.0*i);
		}
		//general info
		TEST_EQUAL(spec.getMSLevel(),2)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::FULL)
		TEST_EQUAL(spec.getType(),SpectrumSettings::PEAKS)
		TEST_REAL_SIMILAR(spec.getRT(),5.2)
		TEST_EQUAL(spec.getInstrumentSettings().getPolarity(),IonSource::POSITIVE)
		TEST_EQUAL(spec.getInstrumentSettings().getScanWindows().size(),3)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[0].begin,100.0)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[0].end,500.0)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[1].begin,600.0)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[1].end,1000.0)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[2].begin,1100.0)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[2].end,1500.0)
		TEST_EQUAL(spec.getAcquisitionInfo().getMethodOfCombination(),"no combination")
		TEST_EQUAL(spec.getAcquisitionInfo().size(),1)
		TEST_EQUAL(spec.getAcquisitionInfo()[0].getNumber(),0)
		//meta data arrays
		TEST_EQUAL(spec.getMetaDataArrays().size(),2)
		TEST_STRING_EQUAL(spec.getMetaDataArrays()[0].getName(),"signal to noise array")
		TEST_EQUAL(spec.getMetaDataArrays()[0].size(),10)
		TEST_STRING_EQUAL(spec.getMetaDataArrays()[1].getName(),"user-defined name")
		TEST_EQUAL(spec.getMetaDataArrays()[1].size(),10)
		//precursor
		TEST_REAL_SIMILAR(spec.getPrecursorPeak().getIntensity(),120053)
		TEST_EQUAL(spec.getPrecursorPeak().getCharge(),2)
		TEST_REAL_SIMILAR(spec.getPrecursorPeak().getPosition()[0],5.55)
		TEST_EQUAL(spec.getPrecursor().getActivationMethod(),Precursor::CID)
		TEST_REAL_SIMILAR(spec.getPrecursor().getActivationEnergy(),35)
		TEST_EQUAL(spec.getPrecursorPeak().getPossibleChargeStates().size(),3)
		TEST_EQUAL(spec.getPrecursorPeak().getPossibleChargeStates()[0],1)
		TEST_EQUAL(spec.getPrecursorPeak().getPossibleChargeStates()[1],3)
		TEST_EQUAL(spec.getPrecursorPeak().getPossibleChargeStates()[2],4)
		TEST_STRING_EQUAL(spec.getPrecursor().getMetaValue("source_file_name"),"pr.dta")
		TEST_STRING_EQUAL(spec.getPrecursor().getMetaValue("source_file_path"),"file:///F:/data/Exp03")
//		TEST_STRING_EQUAL(spec.getPrecursor().getMetaValue("external_native_id"),"ENI1")
//		TEST_STRING_EQUAL(spec.getPrecursor().getMetaValue("external_spectrum_id"),"ESI1")
		//source file
		TEST_STRING_EQUAL(spec.getSourceFile().getNameOfFile(),"tiny1.dta")
		TEST_STRING_EQUAL(spec.getSourceFile().getPathToFile(),"file:///F:/data/Exp01")
		TEST_STRING_EQUAL(spec.getSourceFile().getChecksum(),"81be39fb2700ab2f3c8b2234b91274968b6899b1")
		TEST_EQUAL(spec.getSourceFile().getChecksumType(),SourceFile::SHA1)
		//ids
		TEST_STRING_EQUAL(spec.getNativeID(),"index=1")
		TEST_STRING_EQUAL(spec.getMetaValue("maldi_spot_id"),"M1")
		
	}
	
	//-------------------------- spectrum 2 --------------------------
	{
		const MSSpectrum<>& spec = exp[2]; 
		//peaks
		TEST_EQUAL(spec.size(),15)
		for (UInt i=0; i<15; ++i)
		{
			TEST_REAL_SIMILAR(spec[i].getMZ(),i);
			TEST_REAL_SIMILAR(spec[i].getIntensity(),15-i);
		}
		//general info		
		TEST_EQUAL(spec.getMSLevel(),1)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::FULL)
		TEST_EQUAL(spec.getMetaDataArrays().size(),0)
		TEST_EQUAL(spec.getType(),SpectrumSettings::PEAKS)
		TEST_REAL_SIMILAR(spec.getRT(),5.3)
		TEST_EQUAL(spec.getInstrumentSettings().getPolarity(),IonSource::POSITIVE)
		TEST_EQUAL(spec.getInstrumentSettings().getScanWindows().size(),1)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[0].begin,400.0)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[0].end,1800.0)
		//acquisition
		TEST_STRING_EQUAL(spec.getAcquisitionInfo().getMethodOfCombination(),"median of spectra")
		TEST_EQUAL(spec.getAcquisitionInfo().size(),2)
		TEST_EQUAL(spec.getAcquisitionInfo()[0].getNumber(),4711)
		TEST_EQUAL(spec.getAcquisitionInfo()[1].getNumber(),4712)
		TEST_EQUAL(spec.getSourceFile()==SourceFile(),true)
		//ids
		TEST_STRING_EQUAL(spec.getNativeID(),"index=2")
		TEST_STRING_EQUAL(spec.getMetaValue("maldi_spot_id"),"M2")
	}

	//-------------------------- spectrum 3 (no peaks) --------------------------
	{
		const MSSpectrum<>& spec = exp[3]; 
		//peaks
		TEST_EQUAL(spec.size(),0)
		//general info
		TEST_EQUAL(spec.getMSLevel(),1)
		TEST_REAL_SIMILAR(spec.getRT(),5.4)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::FULL)
		TEST_EQUAL(spec.getMetaDataArrays().size(),0)
		TEST_EQUAL(spec.getType(),SpectrumSettings::RAWDATA)
		TEST_EQUAL(spec.getInstrumentSettings().getScanWindows().size(),1)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[0].begin,110.0)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[0].end,905.0)
		TEST_STRING_EQUAL(spec.getAcquisitionInfo().getMethodOfCombination(),"no combination")
		TEST_EQUAL(spec.getAcquisitionInfo().size(),1)
		TEST_EQUAL(spec.getAcquisitionInfo()[0].getNumber(),0)
		//ids
		TEST_STRING_EQUAL(spec.getNativeID(),"index=3")
		TEST_EQUAL(spec.metaValueExists("maldi_spot_id"),false)
	}
	
	//-------------------------- userParam --------------------------
	//run
  TEST_STRING_EQUAL(exp.getMetaValue("mzml_id"),"document_id")
  TEST_EQUAL(exp.getMetaValue("flag").valueType(),DataValue::STRING_VALUE)
  TEST_STRING_EQUAL((String)exp.getMetaValue("flag"),"")
  TEST_EQUAL(exp.getMetaValue("string").valueType(),DataValue::STRING_VALUE)
  TEST_STRING_EQUAL((String)exp.getMetaValue("string"),"bla")
  TEST_EQUAL(exp.getMetaValue("float").valueType(),DataValue::DOUBLE_VALUE)
  TEST_REAL_SIMILAR((double)exp.getMetaValue("float"),5.11)
  TEST_EQUAL(exp.getMetaValue("int").valueType(),DataValue::INT_VALUE)
  TEST_EQUAL((Int)exp.getMetaValue("int"),5)
	//instrumentConfiguration
	TEST_EQUAL(exp.getInstrument().getIonOptics(),Instrument::MAGNETIC_DEFLECTION)
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
	TEST_STRING_EQUAL(exp.getDataProcessing()[1].getMetaValue("p2").toString(),"value2")
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

	/////////////////////// TESTING SPECIAL CASES ///////////////////////
	
	//load a second time to make sure everything is re-initialized correctly
	MSExperiment<> exp2;
	file.load("data/MzMLFile_1.mzML",exp2);
	TEST_EQUAL(exp==exp2,true)
	
	//load minimal file
	MSExperiment<> exp3;
	file.load("data/MzMLFile_2_minimal.mzML",exp3);
	TEST_EQUAL(exp3.size(),0)

	//load file with huge CDATA and whitespaces in CDATA
	MSExperiment<> exp4;
	file.load("data/MzMLFile_5_long.mzML",exp4);
	TEST_EQUAL(exp4.size(),1)
	TEST_EQUAL(exp4[0].size(),997530)
	
	//load 32 bit data
	MSExperiment<> exp5;
	file.load("data/MzMLFile_6_32bit.mzML",exp5);
	TEST_EQUAL(exp5.size(),4)
	TEST_EQUAL(exp5[0].size(),15)
	TEST_EQUAL(exp5[1].size(),10)
	TEST_EQUAL(exp5[2].size(),15)
	TEST_EQUAL(exp5[3].size(),0)
	
	//test if it works with different peak types
	MSExperiment<RichPeak1D> e_rich;
  file.load("data/MzMLFile_1.mzML",e_rich);
	
END_SECTION

START_SECTION([EXTRA] load only meta data)
	MzMLFile file;
	file.getOptions().setMetadataOnly(true);
	MSExperiment<> exp;
	file.load("data/MzMLFile_1.mzML",exp);

	TEST_EQUAL(exp.size(),0)
	TEST_EQUAL(exp.getIdentifier(),"document_accession");
	TEST_EQUAL(exp.getContacts().size(),2)
	TEST_EQUAL(exp.getSourceFiles().size(),2);
	TEST_EQUAL(exp.getInstrument().getMassAnalyzers().size(),2)
	TEST_EQUAL(exp.getDataProcessing().size(),2)
END_SECTION


START_SECTION([EXTRA] load with restricted MS levels)
	MzMLFile file;
	file.getOptions().addMSLevel(1);
	MSExperiment<> exp;
	file.load("data/MzMLFile_1.mzML",exp);

	TEST_EQUAL(exp.size(),3)
	TEST_REAL_SIMILAR(exp[0].getRT(),5.1)
	TEST_REAL_SIMILAR(exp[1].getRT(),5.3)
	TEST_REAL_SIMILAR(exp[2].getRT(),5.4)
END_SECTION


START_SECTION([EXTRA] load with restricted RT range)
	MzMLFile file;
	file.getOptions().setRTRange(makeRange(5.15,5.35));
	MSExperiment<> exp;
	file.load("data/MzMLFile_1.mzML",exp);
	TEST_EQUAL(exp.size(),2)
	TEST_REAL_SIMILAR(exp[0].getRT(),5.2)
	TEST_REAL_SIMILAR(exp[1].getRT(),5.3)
END_SECTION

START_SECTION([EXTRA] load with restricted m/z range)
	MzMLFile file;
	file.getOptions().setMZRange(makeRange(6.5,9.5));
	MSExperiment<> exp;
	file.load("data/MzMLFile_1.mzML",exp);
	
	TEST_EQUAL(exp.size(),4)
	TEST_EQUAL(exp[0].size(),3)
	TEST_REAL_SIMILAR(exp[0][0].getMZ(),7.0)
	TEST_REAL_SIMILAR(exp[0][1].getMZ(),8.0)
	TEST_REAL_SIMILAR(exp[0][2].getMZ(),9.0)
	TEST_EQUAL(exp[1].size(),1)
	TEST_REAL_SIMILAR(exp[1][0].getMZ(),8.0)
	TEST_EQUAL(exp[2].size(),3)
	TEST_REAL_SIMILAR(exp[2][0].getMZ(),7.0)
	TEST_REAL_SIMILAR(exp[2][1].getMZ(),8.0)
	TEST_REAL_SIMILAR(exp[2][2].getMZ(),9.0)
	TEST_EQUAL(exp[3].size(),0)
END_SECTION

START_SECTION([EXTRA] load intensity range)
	MzMLFile file;
	file.getOptions().setIntensityRange(makeRange(6.5,9.5));
	MSExperiment<> exp;
	file.load("data/MzMLFile_1.mzML",exp);

	TEST_EQUAL(exp.size(),4)
	TEST_EQUAL(exp[0].size(),3)
	TEST_REAL_SIMILAR(exp[0][0].getIntensity(),9.0)
	TEST_REAL_SIMILAR(exp[0][1].getIntensity(),8.0)
	TEST_REAL_SIMILAR(exp[0][2].getIntensity(),7.0)
	TEST_EQUAL(exp[1].size(),1)
	TEST_REAL_SIMILAR(exp[1][0].getIntensity(),8.0)
	TEST_EQUAL(exp[2].size(),3)
	TEST_REAL_SIMILAR(exp[2][0].getIntensity(),9.0)
	TEST_REAL_SIMILAR(exp[2][1].getIntensity(),8.0)
	TEST_REAL_SIMILAR(exp[2][2].getIntensity(),7.0)
	TEST_EQUAL(exp[3].size(),0)
END_SECTION

START_SECTION((template <typename MapType> void store(const String& filename, const MapType& map) const))
	MzMLFile file;
	
	//load map
	MSExperiment<> exp_original;
	file.load("data/MzMLFile_1.mzML",exp_original);
 	
 	//store map
	std::string tmp_filename;
 	NEW_TMP_FILE(tmp_filename);
	file.store(tmp_filename,exp_original);
	
	//load written map
	MSExperiment<> exp;
	file.load(tmp_filename,exp);
	
	//test if everything worked
	TEST_EQUAL(exp==exp_original,true)
	
END_SECTION

START_SECTION([EXTRA] bool isValid(const String& filename))
	std::string tmp_filename;
  MzMLFile file;
  MSExperiment<> e;
  
  //written empty file
	NEW_TMP_FILE(tmp_filename);
  file.store(tmp_filename,e);
  TEST_EQUAL(file.isValid(tmp_filename),true);

	//written filled file
	NEW_TMP_FILE(tmp_filename);
	file.load("data/MzMLFile_1.mzML",e);
  file.store(tmp_filename,e);
  TEST_EQUAL(file.isValid(tmp_filename),true);
	
	//indexed file
	TEST_EQUAL(file.isValid("data/MzMLFile_4_indexed.mzML"),true)
END_SECTION

START_SECTION( bool isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings))
	std::string tmp_filename;
	MzMLFile file;
	StringList errors, warnings;
  MSExperiment<> e;
	
  //written empty file
	NEW_TMP_FILE(tmp_filename);
  file.store(tmp_filename,e);
  TEST_EQUAL(file.isSemanticallyValid(tmp_filename, errors, warnings),true);
	TEST_EQUAL(errors.size(),0)
	TEST_EQUAL(warnings.size(),0)
	
	//written filled file
	NEW_TMP_FILE(tmp_filename);
	file.load("data/MzMLFile_1.mzML",e);
  file.store(tmp_filename,e);
  TEST_EQUAL(file.isSemanticallyValid(tmp_filename, errors, warnings),true);
	TEST_EQUAL(errors.size(),0)
	TEST_EQUAL(warnings.size(),0)

	//valid file
	TEST_EQUAL(file.isSemanticallyValid("data/MzMLFile_1.mzML", errors, warnings),true)
	TEST_EQUAL(errors.size(),0)
	TEST_EQUAL(warnings.size(),2)

	//indexed MzML
	TEST_EQUAL(file.isSemanticallyValid("data/MzMLFile_4_indexed.mzML", errors, warnings),true)
	TEST_EQUAL(errors.size(),0)
	TEST_EQUAL(warnings.size(),0)
	
	//invalid file
	TEST_EQUAL(file.isSemanticallyValid("data/MzMLFile_3_invalid.mzML", errors, warnings),false)
	TEST_EQUAL(errors.size(),5)
	TEST_EQUAL(warnings.size(),0)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

