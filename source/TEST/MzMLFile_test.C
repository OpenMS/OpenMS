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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

DRange<1> makeRange(DoubleReal a, DoubleReal b)
{
	DPosition<1> pa(a), pb(b);
	return DRange<1>(pa, pb);
}

///////////////////////////

START_TEST(MzMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

//Note: This code generates the test files for meta data arrays of differnet types. Do not delete it!

////template spectrum with 100 peaks
//MSSpectrum<> template_spec;
//for (Size i=0; i<100; ++i)
//{
//	Peak1D p;
//	p.setIntensity(i);
//	p.setMZ(i);
//	template_spec.push_back(p);
//}
//
//MSSpectrum<> spec;
//MSExperiment<> exp;
//Size spectrum_number = 0;
//Size array_number = 1;
//
////spectrum 1 - 3 float arrays of size 50,100,200
//spec = template_spec; ++spectrum_number; array_number = 1;
//spec.setNativeID(String("index=") + spectrum_number);
//spec.setRT(1.0 * spectrum_number);
//spec.setName(String("spectum number=") + spectrum_number);
//Size array_size=50;
//for (Size i=0; i<3; ++i)
//{
//	spec.getFloatDataArrays().resize(i+1);
//	for (Size j=0; j<array_size; ++j)
//	{
//		spec.getFloatDataArrays()[i].push_back(100*(i+1) + j);
//	}
//	spec.getFloatDataArrays()[i].setName(String("array number=") + array_number);
//	array_size *=2;
//	array_number +=1;
//}
//exp.push_back(spec);
//
////spectrum 2 - 3 string arrays of size 50,100,200
//spec = template_spec; ++spectrum_number; array_number = 1;
//spec.setNativeID(String("index=") + spectrum_number);
//spec.setRT(1.0 * spectrum_number);
//spec.setName(String("spectum number=") + spectrum_number);
//array_size=50;
//for (Size i=0; i<3; ++i)
//{
//	spec.getStringDataArrays().resize(i+1);
//	for (Size j=0; j<array_size; ++j)
//	{
//		spec.getStringDataArrays()[i].push_back(String(100*(i+1) + j));
//	}
//	spec.getStringDataArrays()[i].setName(String("array number=") + array_number);
//	array_size *=2;
//	array_number +=1;
//}
//exp.push_back(spec);
//
////spectrum 3 - 3 integer arrays of size 50,100,200
//spec = template_spec; ++spectrum_number; array_number = 1;
//spec.setNativeID(String("index=") + spectrum_number);
//spec.setRT(1.0 * spectrum_number);
//spec.setName(String("spectum number=") + spectrum_number);
//array_size=50;
//for (Size i=0; i<3; ++i)
//{
//	spec.getIntegerDataArrays().resize(i+1);
//	for (Size j=0; j<array_size; ++j)
//	{
//		spec.getIntegerDataArrays()[i].push_back(100*(i+1) + j);
//	}
//	spec.getIntegerDataArrays()[i].setName(String("array number=") + array_number);
//	array_size *=2;
//	array_number +=1;
//}
//exp.push_back(spec);
//
//
////spectrum 4 - 2 float arrays of size 50,100 + 1 string arrays of size 200 + 3 integer arrays of size 50,100,200
//spec = template_spec; ++spectrum_number; array_number = 1;
//spec.setNativeID(String("index=") + spectrum_number);
//spec.setRT(1.0 * spectrum_number);
//spec.setName(String("spectum number=") + spectrum_number);
//array_size=50;
//for (Size i=0; i<2; ++i)
//{
//	spec.getFloatDataArrays().resize(i+1);
//	for (Size j=0; j<array_size; ++j)
//	{
//		spec.getFloatDataArrays()[i].push_back(100*(i+1) + j);
//	}
//	spec.getFloatDataArrays()[i].setName(String("array number=") + array_number);
//	array_size *=2;
//	array_number +=1;
//}
//array_size=200;
//for (Size i=0; i<1; ++i)
//{
//	spec.getStringDataArrays().resize(i+1);
//	for (Size j=0; j<array_size; ++j)
//	{
//		spec.getStringDataArrays()[i].push_back(String(100*(i+1) + j));
//	}
//	spec.getStringDataArrays()[i].setName(String("array number=") + array_number);
//	array_size *=2;
//	array_number +=1;
//}
//array_size=50;
//for (Size i=0; i<3; ++i)
//{
//	spec.getIntegerDataArrays().resize(i+1);
//	for (Size j=0; j<array_size; ++j)
//	{
//		spec.getIntegerDataArrays()[i].push_back(100*(i+1) + j);
//	}
//	spec.getIntegerDataArrays()[i].setName(String("array number=") + array_number);
//	array_size *=2;
//	array_number +=1;
//}
//exp.push_back(spec);
//
//MzMLFile f;
//f.store("data/MzMLFile_6_uncompressed.mzML",exp);
//f.getOptions().setCompression(true);
//f.store("data/MzMLFile_6_compressed.mzML",exp);

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

START_SECTION(([EXTRA] Chromatogram section))
	MzMLFile file;
	PeakMap exp;
	//file.load("/share/usr/bertsch/MRM_example_1.1.0.mzML", exp);
	//TEST_EQUAL(exp.getChromatograms().size(), 2)
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
	file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),exp);

	//test DocumentIdentifier addition
	TEST_STRING_EQUAL(exp.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
	TEST_STRING_EQUAL(FileHandler::typeToName(exp.getLoadedFileType()),"mzML");

	//-------------------------- general information --------------------------

	TEST_EQUAL(exp.size(),4)
	//run
	TEST_EQUAL(exp.getIdentifier(),"document_accession")
	TEST_EQUAL(exp.getFractionIdentifier(),"the_best_fraction_ever")
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
	TEST_EQUAL(exp.getSourceFiles().size(),1);
	TEST_STRING_EQUAL(exp.getSourceFiles()[0].getNameOfFile(),"tiny1.RAW")
	TEST_STRING_EQUAL(exp.getSourceFiles()[0].getPathToFile(),"file:///F:/data/Exp01")
	TEST_STRING_EQUAL(exp.getSourceFiles()[0].getChecksum(),"71be39fb2700ab2f3c8b2234b91274968b6899b1")
	TEST_EQUAL(exp.getSourceFiles()[0].getChecksumType(),SourceFile::SHA1)
	TEST_STRING_EQUAL(exp.getSourceFiles()[0].getFileType(),"Thermo RAW file")
	TEST_STRING_EQUAL(exp.getSourceFiles()[0].getNativeIDType(),"multiple peak list nativeID format")
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
	TEST_REAL_SIMILAR(exp.getInstrument().getMassAnalyzers()[0].getAccuracy(),10.5)
	TEST_REAL_SIMILAR(exp.getInstrument().getMassAnalyzers()[0].getMagneticFieldStrength(),14.56)
	TEST_REAL_SIMILAR(exp.getInstrument().getMassAnalyzers()[0].getTOFTotalPathLength(),11.1)
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
	
	//-------------------------- spectrum 0 --------------------------
	{
		const MSSpectrum<>& spec = exp[0];
		//peaks
		TEST_EQUAL(spec.size(),15)
		for (UInt i=0; i<15; ++i)
		{
			TEST_REAL_SIMILAR(spec[(Size)i].getMZ(),i);
			TEST_REAL_SIMILAR(spec[(Size)i].getIntensity(),15-i);
		}
		//general info
		TEST_EQUAL(spec.getMSLevel(),1)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::MASSSPECTRUM)
		TEST_EQUAL(spec.getFloatDataArrays().size(),0)
		TEST_EQUAL(spec.getType(),SpectrumSettings::PEAKS)
		TEST_REAL_SIMILAR(spec.getRT(),5.1)
		TEST_EQUAL(spec.getInstrumentSettings().getScanWindows().size(),1)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[0].begin,400.0)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[0].end,1800.0)
		TEST_STRING_EQUAL(spec.getAcquisitionInfo().getMethodOfCombination(),"median of spectra")
		TEST_EQUAL(spec.getAcquisitionInfo().size(),2)
		TEST_EQUAL(spec.getAcquisitionInfo()[0].getIdentifier(),"4711")
		TEST_STRING_EQUAL(spec.getAcquisitionInfo()[0].getMetaValue("source_file_name"),"ac.dta")
		TEST_STRING_EQUAL(spec.getAcquisitionInfo()[0].getMetaValue("source_file_path"),"file:///F:/data/Exp02")
		TEST_EQUAL(spec.getAcquisitionInfo()[1].getIdentifier(),"4712")
		TEST_EQUAL(spec.getSourceFile()==SourceFile(),true)
		//ids
		TEST_STRING_EQUAL(spec.getNativeID(),"index=0")
		TEST_STRING_EQUAL(spec.getMetaValue("maldi_spot_id"),"M0")
		//precursors
		TEST_EQUAL(spec.getPrecursors().size(),0)
		TEST_EQUAL(spec.getProducts().size(),0)
		//data processing
		TEST_EQUAL(spec.getDataProcessing().size(),2)
		TEST_EQUAL(spec.getDataProcessing()[0].getSoftware().getName(),"Xcalibur")
		TEST_EQUAL(spec.getDataProcessing()[0].getSoftware().getVersion(),"2.0.5")
		TEST_EQUAL(spec.getDataProcessing()[0].getProcessingActions().size(),2)
		TEST_EQUAL(spec.getDataProcessing()[0].getProcessingActions().count(DataProcessing::DEISOTOPING),1)
		TEST_EQUAL(spec.getDataProcessing()[0].getProcessingActions().count(DataProcessing::CHARGE_DECONVOLUTION),1)
		TEST_STRING_EQUAL(spec.getDataProcessing()[0].getCompletionTime().get(),"2001-02-03 04:05:00")
		TEST_REAL_SIMILAR(DoubleReal(spec.getDataProcessing()[0].getMetaValue("low_intensity_threshold")),5.9)
		TEST_REAL_SIMILAR(DoubleReal(spec.getDataProcessing()[0].getMetaValue("high_intensity_threshold")),10.9)
		TEST_EQUAL(spec.getDataProcessing()[0].isMetaEmpty(),false)
		TEST_EQUAL(spec.getDataProcessing()[1].getSoftware().getName(),"ProteoWizard")
		TEST_EQUAL(spec.getDataProcessing()[1].getSoftware().getVersion(),"1.0")
		TEST_EQUAL(spec.getDataProcessing()[1].getProcessingActions().size(),1)
		TEST_EQUAL(spec.getDataProcessing()[1].getProcessingActions().count(DataProcessing::CONVERSION_MZML),1)
		TEST_EQUAL(spec.getDataProcessing()[1].isMetaEmpty(),false)
}

	//-------------------------- spectrum 1 --------------------------
	{
		const MSSpectrum<>& spec = exp[1];
		//peaks
		TEST_EQUAL(spec.size(),10)
		for (Size i=0; i<10; ++i)
		{
			TEST_REAL_SIMILAR(spec[i].getMZ(),2.0*i);
			TEST_REAL_SIMILAR(spec[i].getIntensity(),20-2.0*i);
		}
		//general info
		TEST_EQUAL(spec.getMSLevel(),2)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::MASSSPECTRUM)
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
		TEST_EQUAL(spec.getAcquisitionInfo()[0].getIdentifier(),"0")
		//meta data arrays
		TEST_EQUAL(spec.getFloatDataArrays().size(),2)
		TEST_STRING_EQUAL(spec.getFloatDataArrays()[0].getName(),"signal to noise array")
		TEST_EQUAL(spec.getFloatDataArrays()[0].size(),10)
		TEST_EQUAL(spec.getFloatDataArrays()[0].getDataProcessing().size(),1)
		TEST_EQUAL(spec.getFloatDataArrays()[0].getDataProcessing()[0].getSoftware().getName(), "FileFilter")
		TEST_EQUAL(spec.getFloatDataArrays()[0].getDataProcessing()[0].getSoftware().getVersion(), "1.6.1")
		TEST_EQUAL(spec.getFloatDataArrays()[0].getDataProcessing()[0].getProcessingActions().size(), 1)
		TEST_EQUAL(spec.getFloatDataArrays()[0].getDataProcessing()[0].getProcessingActions().count(DataProcessing::CHARGE_CALCULATION), 1)
		TEST_STRING_EQUAL(spec.getFloatDataArrays()[0].getDataProcessing()[0].getCompletionTime().get(),"2001-02-03 04:15:00")
		TEST_STRING_EQUAL(spec.getFloatDataArrays()[1].getName(),"user-defined name")
		TEST_EQUAL(spec.getFloatDataArrays()[1].getDataProcessing().size(),0)
		TEST_EQUAL(spec.getFloatDataArrays()[1].size(),10)
		//precursors
		TEST_EQUAL(spec.getPrecursors().size(),2)
		TEST_REAL_SIMILAR(spec.getPrecursors()[0].getIntensity(),120053)
		TEST_EQUAL(spec.getPrecursors()[0].getCharge(),2)
		TEST_REAL_SIMILAR(spec.getPrecursors()[0].getMZ(),5.55)
		TEST_EQUAL(spec.getPrecursors()[0].getActivationMethods().size(),2)
		TEST_EQUAL(spec.getPrecursors()[0].getActivationMethods().count(Precursor::CID),1)
		TEST_EQUAL(spec.getPrecursors()[0].getActivationMethods().count(Precursor::PD),1)
		TEST_REAL_SIMILAR(spec.getPrecursors()[0].getActivationEnergy(),35)
		TEST_REAL_SIMILAR(spec.getPrecursors()[0].getIsolationWindowLowerOffset(),6.66)
		TEST_REAL_SIMILAR(spec.getPrecursors()[0].getIsolationWindowUpperOffset(),7.77)
		TEST_EQUAL(spec.getPrecursors()[0].getPossibleChargeStates().size(),3)
		TEST_EQUAL(spec.getPrecursors()[0].getPossibleChargeStates()[0],1)
		TEST_EQUAL(spec.getPrecursors()[0].getPossibleChargeStates()[1],3)
		TEST_EQUAL(spec.getPrecursors()[0].getPossibleChargeStates()[2],4)
		TEST_REAL_SIMILAR(spec.getPrecursors()[1].getMZ(),15.55)
		TEST_REAL_SIMILAR(spec.getPrecursors()[1].getIsolationWindowLowerOffset(),16.66)
		TEST_REAL_SIMILAR(spec.getPrecursors()[1].getIsolationWindowUpperOffset(),17.77)
		TEST_EQUAL(spec.getPrecursors()[1].getActivationMethods().size(),1)
		TEST_EQUAL(spec.getPrecursors()[1].getActivationMethods().count(Precursor::ETD),1)
		TEST_REAL_SIMILAR(spec.getPrecursors()[1].getActivationEnergy(),36)
		TEST_REAL_SIMILAR(spec.getPrecursors()[1].getIntensity(),0.0f)
		TEST_EQUAL(spec.getPrecursors()[1].getCharge(),0)
		TEST_EQUAL(spec.getPrecursors()[1].getPossibleChargeStates().size(),0)
		//products
		TEST_EQUAL(spec.getProducts().size(),0)
		//source file
		TEST_STRING_EQUAL(spec.getSourceFile().getNameOfFile(),"tiny1.dta")
		TEST_STRING_EQUAL(spec.getSourceFile().getPathToFile(),"file:///F:/data/Exp01")
		TEST_STRING_EQUAL(spec.getSourceFile().getChecksum(),"81be39fb2700ab2f3c8b2234b91274968b6899b1")
		TEST_EQUAL(spec.getSourceFile().getChecksumType(),SourceFile::SHA1)
		//ids
		TEST_STRING_EQUAL(spec.getNativeID(),"index=1")
		TEST_STRING_EQUAL(spec.getMetaValue("maldi_spot_id"),"M1")
		//data processing
		TEST_EQUAL(spec.getDataProcessing().size(),2)
		TEST_EQUAL(spec.getDataProcessing()[0].getSoftware().getName(),"Xcalibur")
		TEST_EQUAL(spec.getDataProcessing()[0].getSoftware().getVersion(),"2.0.5")
		TEST_EQUAL(spec.getDataProcessing()[0].getProcessingActions().size(),2)
		TEST_EQUAL(spec.getDataProcessing()[0].getProcessingActions().count(DataProcessing::DEISOTOPING),1)
		TEST_EQUAL(spec.getDataProcessing()[0].getProcessingActions().count(DataProcessing::CHARGE_DECONVOLUTION),1)
		TEST_STRING_EQUAL(spec.getDataProcessing()[0].getCompletionTime().get(),"2001-02-03 04:05:00")
		TEST_REAL_SIMILAR(DoubleReal(spec.getDataProcessing()[0].getMetaValue("low_intensity_threshold")),5.9)
		TEST_REAL_SIMILAR(DoubleReal(spec.getDataProcessing()[0].getMetaValue("high_intensity_threshold")),10.9)
		TEST_EQUAL(spec.getDataProcessing()[0].isMetaEmpty(),false)
		TEST_EQUAL(spec.getDataProcessing()[1].getSoftware().getName(),"ProteoWizard")
		TEST_EQUAL(spec.getDataProcessing()[1].getSoftware().getVersion(),"1.0")
		TEST_EQUAL(spec.getDataProcessing()[1].getProcessingActions().size(),1)
		TEST_EQUAL(spec.getDataProcessing()[1].getProcessingActions().count(DataProcessing::CONVERSION_MZML),1)
		TEST_EQUAL(spec.getDataProcessing()[1].isMetaEmpty(),false)

	}

	//-------------------------- spectrum 2 --------------------------
	{
		const MSSpectrum<>& spec = exp[2];
		//peaks
		TEST_EQUAL(spec.size(),15)
		for (UInt i=0; i<15; ++i)
		{
			TEST_REAL_SIMILAR(spec[(Size)i].getMZ(),i);
			TEST_REAL_SIMILAR(spec[(Size)i].getIntensity(),15-i);
		}
		//general info
		TEST_EQUAL(spec.getMSLevel(),1)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::MASSSPECTRUM)
		TEST_EQUAL(spec.getFloatDataArrays().size(),0)
		TEST_EQUAL(spec.getType(),SpectrumSettings::PEAKS)
		TEST_REAL_SIMILAR(spec.getRT(),5.3)
		TEST_EQUAL(spec.getInstrumentSettings().getPolarity(),IonSource::POSITIVE)
		TEST_EQUAL(spec.getInstrumentSettings().getScanWindows().size(),1)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[0].begin,400.0)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[0].end,1800.0)
		//acquisition
		TEST_STRING_EQUAL(spec.getAcquisitionInfo().getMethodOfCombination(),"median of spectra")
		TEST_EQUAL(spec.getAcquisitionInfo().size(),2)
		TEST_EQUAL(spec.getAcquisitionInfo()[0].getIdentifier(),"4711")
		TEST_EQUAL(spec.getAcquisitionInfo()[1].getIdentifier(),"4712")
		TEST_EQUAL(spec.getSourceFile()==SourceFile(),true)
		//ids
		TEST_STRING_EQUAL(spec.getNativeID(),"index=2")
		TEST_STRING_EQUAL(spec.getMetaValue("maldi_spot_id"),"M2")
		//precursors
		TEST_EQUAL(spec.getPrecursors().size(),0)
		//products
		TEST_EQUAL(spec.getProducts().size(),2)
		TEST_REAL_SIMILAR(spec.getProducts()[0].getMZ(),18.88)
		TEST_REAL_SIMILAR(spec.getProducts()[0].getIsolationWindowLowerOffset(),1.0)
		TEST_REAL_SIMILAR(spec.getProducts()[0].getIsolationWindowUpperOffset(),2.0)
		TEST_REAL_SIMILAR(spec.getProducts()[1].getMZ(),19.99)
		TEST_REAL_SIMILAR(spec.getProducts()[1].getIsolationWindowLowerOffset(),3.0)
		TEST_REAL_SIMILAR(spec.getProducts()[1].getIsolationWindowUpperOffset(),4.0)
		//data processing
		TEST_EQUAL(spec.getDataProcessing().size(),1)
		TEST_EQUAL(spec.getDataProcessing()[0].getSoftware().getName(),"Xcalibur")
		TEST_EQUAL(spec.getDataProcessing()[0].getSoftware().getVersion(),"2.0.5")
		TEST_EQUAL(spec.getDataProcessing()[0].getProcessingActions().size(),2)
		TEST_EQUAL(spec.getDataProcessing()[0].getProcessingActions().count(DataProcessing::DEISOTOPING),1)
		TEST_EQUAL(spec.getDataProcessing()[0].getProcessingActions().count(DataProcessing::CHARGE_DECONVOLUTION),1)
		TEST_STRING_EQUAL(spec.getDataProcessing()[0].getCompletionTime().get(),"2001-02-03 04:05:00")
		TEST_REAL_SIMILAR(DoubleReal(spec.getDataProcessing()[0].getMetaValue("low_intensity_threshold")),5.9)
		TEST_REAL_SIMILAR(DoubleReal(spec.getDataProcessing()[0].getMetaValue("high_intensity_threshold")),10.9)
		TEST_EQUAL(spec.getDataProcessing()[0].isMetaEmpty(),false)
	}

	//-------------------------- spectrum 3 (no peaks) --------------------------
	{
		const MSSpectrum<>& spec = exp[3];
		//peaks
		TEST_EQUAL(spec.size(),0)
		//general info
		TEST_EQUAL(spec.getMSLevel(),1)
		TEST_REAL_SIMILAR(spec.getRT(),5.4)
		TEST_EQUAL(spec.getInstrumentSettings().getScanMode(),InstrumentSettings::MASSSPECTRUM)
		TEST_EQUAL(spec.getInstrumentSettings().getZoomScan(),true)
		TEST_EQUAL(spec.getFloatDataArrays().size(),0)
		TEST_EQUAL(spec.getType(),SpectrumSettings::RAWDATA)
		TEST_EQUAL(spec.getInstrumentSettings().getScanWindows().size(),1)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[0].begin,110.0)
		TEST_REAL_SIMILAR(spec.getInstrumentSettings().getScanWindows()[0].end,905.0)
		TEST_STRING_EQUAL(spec.getAcquisitionInfo().getMethodOfCombination(),"no combination")
		TEST_EQUAL(spec.getAcquisitionInfo().size(),1)
		TEST_EQUAL(spec.getAcquisitionInfo()[0].getIdentifier(),"0")
		//ids
		TEST_STRING_EQUAL(spec.getNativeID(),"index=3")
		TEST_EQUAL(spec.metaValueExists("maldi_spot_id"),false)
		//precursors
		TEST_EQUAL(spec.getPrecursors().size(),0)
		TEST_EQUAL(spec.getProducts().size(),0)
		//data processing
		TEST_EQUAL(spec.getDataProcessing().size(),1)
		TEST_EQUAL(spec.getDataProcessing()[0].getSoftware().getName(),"ProteoWizard")
		TEST_EQUAL(spec.getDataProcessing()[0].getSoftware().getVersion(),"1.0")
		TEST_EQUAL(spec.getDataProcessing()[0].getProcessingActions().size(),1)
		TEST_EQUAL(spec.getDataProcessing()[0].getProcessingActions().count(DataProcessing::CONVERSION_MZML),1)
		TEST_EQUAL(spec.getDataProcessing()[0].isMetaEmpty(),false)

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
	TEST_STRING_EQUAL((String)exp.getSample().getMetaValue("brenda source tissue"),"cardiac muscle")
	TEST_STRING_EQUAL((String)exp.getSample().getMetaValue("GO cellular component"),"nucleus")
	TEST_STRING_EQUAL((String)exp.getSample().getMetaValue("cellular quality"),"11.11")
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
	TEST_STRING_EQUAL((String)exp[1].getFloatDataArrays()[0].getMetaValue("name"),"binaryDataArray_sn")
	TEST_STRING_EQUAL((String)exp[1].getFloatDataArrays()[0].getMetaValue("name2"),"binaryDataArray_sn2")
	TEST_STRING_EQUAL((String)exp[1].getFloatDataArrays()[1].getMetaValue("name"),"binaryDataArray_c")
	TEST_STRING_EQUAL((String)exp[1].getFloatDataArrays()[1].getMetaValue("name2"),"")
	//acquisition list
	TEST_STRING_EQUAL((String)exp[0].getAcquisitionInfo().getMetaValue("name"),"acquisition_list")
	//acquisition
	TEST_STRING_EQUAL((String)exp[0].getAcquisitionInfo()[0].getMetaValue("name"),"acquisition1")
	TEST_STRING_EQUAL((String)exp[0].getAcquisitionInfo()[1].getMetaValue("name"),"acquisition2")
	//source file
	TEST_STRING_EQUAL((String)exp.getSourceFiles()[0].getMetaValue("name"),"sourcefile1")
	TEST_STRING_EQUAL((String)exp[1].getSourceFile().getMetaValue("name"),"sourcefile4")
	//data processing
	TEST_STRING_EQUAL(exp[0].getDataProcessing()[0].getMetaValue("p1").toString(),"value1")
	TEST_STRING_EQUAL(exp[0].getDataProcessing()[1].getMetaValue("p2").toString(),"value2")
	TEST_STRING_EQUAL(exp[1].getDataProcessing()[0].getMetaValue("p1").toString(),"value1")
	TEST_STRING_EQUAL(exp[1].getDataProcessing()[1].getMetaValue("p2").toString(),"value2")
	TEST_STRING_EQUAL(exp[2].getDataProcessing()[0].getMetaValue("p1").toString(),"value1")
	TEST_STRING_EQUAL(exp[3].getDataProcessing()[0].getMetaValue("p2").toString(),"value2")
	TEST_STRING_EQUAL(exp[1].getFloatDataArrays()[0].getDataProcessing()[0].getMetaValue("p3").toString(),"value3")
	//precursor
	TEST_STRING_EQUAL(exp[1].getPrecursors()[0].getMetaValue("iwname").toString(),"isolationwindow1")
	TEST_STRING_EQUAL(exp[1].getPrecursors()[0].getMetaValue("siname").toString(),"selectedion1")
	TEST_STRING_EQUAL(exp[1].getPrecursors()[0].getMetaValue("acname").toString(),"activation1")
	TEST_STRING_EQUAL(exp[1].getPrecursors()[1].getMetaValue("acname").toString(),"activation2")
	TEST_STRING_EQUAL(exp[1].getPrecursors()[1].getMetaValue("iwname").toString(),"isolationwindow2")
	//product
	TEST_STRING_EQUAL(exp[2].getProducts()[0].getMetaValue("iwname").toString(),"isolationwindow3")
	TEST_STRING_EQUAL(exp[2].getProducts()[1].getMetaValue("iwname").toString(),"isolationwindow4")
	//scan window
	TEST_STRING_EQUAL((String)exp[0].getInstrumentSettings().getScanWindows()[0].getMetaValue("name"),"scanwindow1")
	//-------------------------- cvParam (but no member => meta data)--------------------------
	//general
	TEST_STRING_EQUAL((String)exp.getSample().getMetaValue("sample batch"),"4.4")
	//spectrum 1
	TEST_REAL_SIMILAR((DoubleReal)exp[0].getMetaValue("elution time (seconds)"),55.11)
	TEST_REAL_SIMILAR((DoubleReal)exp[0].getMetaValue("lowest observed m/z"),400.39)
	TEST_REAL_SIMILAR((DoubleReal)exp[0].getMetaValue("highest observed m/z"),1795.56)
	TEST_REAL_SIMILAR((DoubleReal)exp[0].getMetaValue("lowest observed wavelength"),500.39)
	TEST_REAL_SIMILAR((DoubleReal)exp[0].getMetaValue("highest observed wavelength"),795.56)
	TEST_REAL_SIMILAR((DoubleReal)exp[0].getMetaValue("base peak m/z"),445.347)
	TEST_REAL_SIMILAR((DoubleReal)exp[0].getMetaValue("base peak intensity"),120054)
	TEST_REAL_SIMILAR((DoubleReal)exp[0].getMetaValue("total ion current"),16675500)
	TEST_STRING_EQUAL((String)exp[0].getMetaValue("spectrum title"),"title")
	TEST_STRING_EQUAL((String)exp[0].getMetaValue("peak list scans"),"15 scans")
	TEST_STRING_EQUAL((String)exp[0].getMetaValue("peak list raw scans"),"16 scans")

	TEST_STRING_EQUAL((String)exp[0].getMetaValue("mass resolution"),"4.3")
	TEST_REAL_SIMILAR((DoubleReal)exp[0].getMetaValue("analyzer scan offset"),-4.5)
	TEST_REAL_SIMILAR((DoubleReal)exp[0].getMetaValue("dwell time"),123.45)
	TEST_STRING_EQUAL((String)exp[0].getMetaValue("filter string"),"+ c NSI Full ms [ 400.00-1800.00]")
	TEST_STRING_EQUAL((String)exp[0].getMetaValue("preset scan configuration"),"3 abc")
	TEST_REAL_SIMILAR((DoubleReal)exp[0].getMetaValue("scan rate"),17.17)
	//spectrum 2
	TEST_STRING_EQUAL((String)exp[1].getMetaValue("mass resolution"),"4.1")
	TEST_STRING_EQUAL((String)exp[1].getPrecursors()[0].getMetaValue("collision gas"), "Argon")
	TEST_STRING_EQUAL((String)exp[1].getPrecursors()[0].getMetaValue("buffer gas"), "Krypton")
	TEST_STRING_EQUAL((String)exp[1].getPrecursors()[0].getMetaValue("source_file_name"),"pr.dta")
	TEST_STRING_EQUAL((String)exp[1].getPrecursors()[0].getMetaValue("source_file_path"),"file:///F:/data/Exp03")

	/////////////////////// TESTING SPECIAL CASES ///////////////////////

	//load a second time to make sure everything is re-initialized correctly
	MSExperiment<> exp2;
	file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),exp2);
	TEST_EQUAL(exp==exp2,true)

	//load minimal file
	MSExperiment<> exp3;
	file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_2_minimal.mzML"),exp3);
	TEST_EQUAL(exp3.size(),0)

	//load file with huge CDATA and whitespaces in CDATA
	MSExperiment<> exp4;
	file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_5_long.mzML"),exp4);
	TEST_EQUAL(exp4.size(),1)
	TEST_EQUAL(exp4[0].size(),997530)
	
	//test 32/64 bit floats, 32/64 bit integer, null terminated strings, zlib compression
	MSExperiment<> exp_ucomp;
	file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML"),exp_ucomp);
	MSExperiment<> exp_comp;
	file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_compressed.mzML"),exp_comp);
	TEST_EQUAL(exp_ucomp.size(),exp_comp.size())
	for (Size s=0; s< exp_ucomp.size(); ++s)
	{
		//check if the same number of peak and meta data arrays is present
		TEST_EQUAL(exp_ucomp[s].size(),exp_comp[s].size())
		TEST_EQUAL(exp_ucomp[s].getFloatDataArrays().size(),exp_comp[s].getFloatDataArrays().size())
		TEST_EQUAL(exp_ucomp[s].getIntegerDataArrays().size(),exp_comp[s].getIntegerDataArrays().size())
		TEST_EQUAL(exp_ucomp[s].getStringDataArrays().size(),exp_comp[s].getStringDataArrays().size())
		//check content of peak array
		for (Size p=0; p< exp_ucomp[s].size(); ++p)
		{
			TEST_REAL_SIMILAR(exp_ucomp[s][p].getMZ(),exp_comp[s][p].getMZ())
			TEST_REAL_SIMILAR(exp_ucomp[s][p].getIntensity(),exp_comp[s][p].getIntensity())
		}
		//check content of float arrays
		for (Size a=0; a<exp_ucomp[s].getFloatDataArrays().size(); ++a)
		{
			for (Size m=0; m< exp_ucomp[s].getFloatDataArrays()[a].size(); ++m)
			{
				TEST_REAL_SIMILAR(exp_ucomp[s].getFloatDataArrays()[a][m],exp_comp[s].getFloatDataArrays()[a][m])
			}
		}
		//check content of integer arrays
		for (Size a=0; a<exp_ucomp[s].getIntegerDataArrays().size(); ++a)
		{
			for (Size m=0; m< exp_ucomp[s].getIntegerDataArrays()[a].size(); ++m)
			{
				TEST_EQUAL(exp_ucomp[s].getIntegerDataArrays()[a][m],exp_comp[s].getIntegerDataArrays()[a][m])
			}
		}
		//check content of string arrays
		for (Size a=0; a<exp_ucomp[s].getStringDataArrays().size(); ++a)
		{
			for (Size m=0; m< exp_ucomp[s].getStringDataArrays()[a].size(); ++m)
			{
				TEST_STRING_EQUAL(exp_ucomp[s].getStringDataArrays()[a][m],exp_comp[s].getStringDataArrays()[a][m])
			}
		}
	}
	
	//test if it works with different peak types
	MSExperiment<RichPeak1D> e_rich;
  file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),e_rich);

END_SECTION

START_SECTION([EXTRA] load only meta data)
	MzMLFile file;
	file.getOptions().setMetadataOnly(true);
	MSExperiment<> exp;
	file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),exp);

	TEST_EQUAL(exp.size(),0)
	TEST_EQUAL(exp.getIdentifier(),"document_accession");
	TEST_EQUAL(exp.getContacts().size(),2)
	TEST_EQUAL(exp.getSourceFiles().size(),1);
	TEST_EQUAL(exp.getInstrument().getMassAnalyzers().size(),2)
END_SECTION


START_SECTION([EXTRA] load with restricted MS levels)
	MzMLFile file;
	file.getOptions().addMSLevel(1);
	MSExperiment<> exp;
	file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),exp);

	TEST_EQUAL(exp.size(),3)
	TEST_REAL_SIMILAR(exp[0].getRT(),5.1)
	TEST_REAL_SIMILAR(exp[1].getRT(),5.3)
	TEST_REAL_SIMILAR(exp[2].getRT(),5.4)
END_SECTION


START_SECTION([EXTRA] load with restricted RT range)
	MzMLFile file;
	file.getOptions().setRTRange(makeRange(5.15,5.35));
	MSExperiment<> exp;
	file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),exp);
	TEST_EQUAL(exp.size(),2)
	TEST_REAL_SIMILAR(exp[0].getRT(),5.2)
	TEST_REAL_SIMILAR(exp[1].getRT(),5.3)
END_SECTION

START_SECTION([EXTRA] load with restricted m/z range)
	MzMLFile file;
	file.getOptions().setMZRange(makeRange(6.5,9.5));
	MSExperiment<> exp;
	file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),exp);

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
	file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),exp);

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
	
	//test with full file
	{
		//load map
		MSExperiment<> exp_original;
		file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),exp_original);
	 	//store map
		std::string tmp_filename;
	 	NEW_TMP_FILE(tmp_filename);
		file.store(tmp_filename,exp_original);
		//load written map
		MSExperiment<> exp;
		file.load(tmp_filename,exp);
		//test if everything worked
		TEST_EQUAL(exp==exp_original,true)
		//NOTE: If it does not work, use this code to find out where the difference is
		TEST_EQUAL(exp.size()==exp_original.size(),true)
		TEST_EQUAL(exp.ExperimentalSettings::operator==(exp_original),true)
		TEST_EQUAL(exp[0].SpectrumSettings::operator==(exp_original[0]),true)
		TEST_EQUAL(exp[0]==exp_original[0],true)
		TEST_EQUAL(exp[1].SpectrumSettings::operator==(exp_original[1]),true)
		TEST_EQUAL(exp[1]==exp_original[1],true)
		TEST_EQUAL(exp[2].SpectrumSettings::operator==(exp_original[2]),true)
		TEST_EQUAL(exp[2]==exp_original[2],true)
		TEST_EQUAL(exp[3].SpectrumSettings::operator==(exp_original[3]),true)
		TEST_EQUAL(exp[3]==exp_original[3],true)
		TEST_EQUAL(exp.getChromatograms().size(), exp_original.getChromatograms().size());
		TEST_EQUAL(exp.getChromatograms() == exp_original.getChromatograms(), true);
	}
	
	//test with empty map
	{
		
		MSExperiment<> empty, exp;
				
		std::string tmp_filename;
		NEW_TMP_FILE(tmp_filename);
		file.store(tmp_filename,empty);
		file.load(tmp_filename,exp);
		TEST_EQUAL(exp==empty,true)
	}
	
	//test with one empty spectrum
	{
		MSExperiment<> empty, exp;
		empty.resize(1);
		empty[0].setRT(17.1234);
		
		//this will be set when writing (forced by mzML)
		empty[0].setNativeID("spectrum=0");
		empty[0].getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
		empty[0].getDataProcessing().resize(1);
		empty[0].getDataProcessing()[0].getProcessingActions().insert(DataProcessing::CONVERSION_MZML);
		empty[0].getAcquisitionInfo().setMethodOfCombination("no combination");
		empty[0].getAcquisitionInfo().resize(1);
		
		std::string tmp_filename;
		NEW_TMP_FILE(tmp_filename);
		file.store(tmp_filename,empty);
		file.load(tmp_filename,exp);
		TEST_EQUAL(exp==empty,true)
		
		//NOTE: If it does not work, use this code to find out where the difference is
//		TEST_EQUAL(exp.size()==empty.size(),true)
//		TEST_EQUAL(exp.ExperimentalSettings::operator==(empty),true)
//		TEST_EQUAL(exp[0].SpectrumSettings::operator==(empty[0]),true)
//		TEST_EQUAL(exp[0]==empty[0],true);
	}

	//test 32/64 bit floats, 32/64 bit integer, null terminated strings, zlib compression
	{
		//load map
		MSExperiment<> exp_original;
		file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML"),exp_original);
	 	//store map
		std::string tmp_filename;
	 	NEW_TMP_FILE(tmp_filename);
		file.getOptions().setCompression(true);
		file.store(tmp_filename,exp_original);
		//load written map
		MSExperiment<> exp;
		file.load(tmp_filename,exp);
		//test if everything worked
		TEST_EQUAL(exp == exp_original,true)
	}

END_SECTION

START_SECTION(bool isValid(const String& filename, std::ostream& os = std::cerr))
	std::string tmp_filename;
  MzMLFile file;
  MSExperiment<> e;

  //written empty file
	NEW_TMP_FILE(tmp_filename);
  file.store(tmp_filename,e);
  TEST_EQUAL(file.isValid(tmp_filename),true);

	//written filled file
	NEW_TMP_FILE(tmp_filename);
	file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),e);
  file.store(tmp_filename,e);
  TEST_EQUAL(file.isValid(tmp_filename),true);

	//indexed file
	TEST_EQUAL(file.isValid(OPENMS_GET_TEST_DATA_PATH("MzMLFile_4_indexed.mzML")),true)
END_SECTION

START_SECTION(bool isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings))
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
	file.load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),e);
  file.store(tmp_filename,e);
  TEST_EQUAL(file.isSemanticallyValid(tmp_filename, errors, warnings),true);
	TEST_EQUAL(errors.size(),0)
	TEST_EQUAL(warnings.size(),10) // TODO add mappings for chromatogram/precursor/activation and selectedIon

	//valid file
	TEST_EQUAL(file.isSemanticallyValid(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), errors, warnings),true)
	TEST_EQUAL(errors.size(),0)
	TEST_EQUAL(warnings.size(),0)
//	for (Size i=0; i<errors.size(); ++i)
//	{
//		cout << "ERROR: " << errors[i] << endl;
//	}
//	for (Size i=0; i<warnings.size(); ++i)
//	{
//		cout << "WARNING: " << warnings[i] << endl;
//	}

	//indexed MzML
	TEST_EQUAL(file.isSemanticallyValid(OPENMS_GET_TEST_DATA_PATH("MzMLFile_4_indexed.mzML"), errors, warnings),true)
	TEST_EQUAL(errors.size(),0)
	TEST_EQUAL(warnings.size(),0)

	//invalid file
	TEST_EQUAL(file.isSemanticallyValid(OPENMS_GET_TEST_DATA_PATH("MzMLFile_3_invalid.mzML"), errors, warnings),false)
	TEST_EQUAL(errors.size(),8)
	TEST_EQUAL(warnings.size(),1)
//	for (Size i=0; i<errors.size(); ++i)
//	{
//		cout << "ERROR: " << errors[i] << endl;
//	}
//	for (Size i=0; i<warnings.size(); ++i)
//	{
//		cout << "WARNING: " << warnings[i] << endl;
//	}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

