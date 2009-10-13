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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ANDIFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FileHandler.h>

using namespace OpenMS;

///////////////////////////

START_TEST(ANDIFile, "$Id$")

/////////////////////////////////////////////////////////////

ANDIFile* ptr = 0;
START_SECTION((ANDIFile()))
	ptr = new ANDIFile;
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~ANDIFile()))
	delete ptr;
END_SECTION

START_SECTION((template<typename MapType> void load(const String& filename, MapType& map) ))
	TOLERANCE_ABSOLUTE(0.01)

	ANDIFile file;
	MSExperiment<> e;

	//test exception
	TEST_EXCEPTION( Exception::FileNotFound , file.load("dummy/dummy.cdf",e) )

  file.load(OPENMS_GET_TEST_DATA_PATH("ANDIFile_test.cdf"),e);
	//test DocumentIdentifier addition
	TEST_STRING_EQUAL(e.getLoadedFilePath(), File::absolutePath(OPENMS_GET_TEST_DATA_PATH("ANDIFile_test.cdf")));
	TEST_STRING_EQUAL(FileHandler::typeToName(e.getLoadedFileType()),"cdf");

  //---------------------------------------------------------------------------
  // 60 : (120,100)
  // 120: (110,100) (120,200) (130,100)
	// 180: (100,100) (110,200) (120,300) (130,200) (140,100)
	//---------------------------------------------------------------------------
  TEST_EQUAL(e.size(), 3)
	TEST_EQUAL(e[0].getMSLevel(), 1)
	TEST_EQUAL(e[1].getMSLevel(), 1)
	TEST_EQUAL(e[2].getMSLevel(), 1)
	TEST_REAL_SIMILAR(e[0].getRT(), 60)
	TEST_REAL_SIMILAR(e[1].getRT(), 120)
	TEST_REAL_SIMILAR(e[2].getRT(), 180)
	TEST_EQUAL(e[0].size(), 1)
	TEST_EQUAL(e[1].size(), 3)
	TEST_EQUAL(e[2].size(), 5)
	TEST_STRING_EQUAL(e[0].getNativeID(), "index=0")
	TEST_STRING_EQUAL(e[1].getNativeID(), "index=1")
	TEST_STRING_EQUAL(e[2].getNativeID(), "index=2")

	TEST_REAL_SIMILAR(e[0][0].getPosition()[0], 120)
	TEST_REAL_SIMILAR(e[0][0].getIntensity(), 100)
	TEST_REAL_SIMILAR(e[1][0].getPosition()[0], 110)
	TEST_REAL_SIMILAR(e[1][0].getIntensity(), 100)
	TEST_REAL_SIMILAR(e[1][1].getPosition()[0], 120)
	TEST_REAL_SIMILAR(e[1][1].getIntensity(), 200)
	TEST_REAL_SIMILAR(e[1][2].getPosition()[0], 130)
	TEST_REAL_SIMILAR(e[1][2].getIntensity(), 100)
	TEST_REAL_SIMILAR(e[2][0].getPosition()[0], 100)
	TEST_REAL_SIMILAR(e[2][0].getIntensity(), 100)
	TEST_REAL_SIMILAR(e[2][1].getPosition()[0], 110)
	TEST_REAL_SIMILAR(e[2][1].getIntensity(), 200)
	TEST_REAL_SIMILAR(e[2][2].getPosition()[0], 120)
	TEST_REAL_SIMILAR(e[2][2].getIntensity(), 300)
	TEST_REAL_SIMILAR(e[2][3].getPosition()[0], 130)
	TEST_REAL_SIMILAR(e[2][3].getIntensity(), 200)
	TEST_REAL_SIMILAR(e[2][4].getPosition()[0], 140)
	TEST_REAL_SIMILAR(e[2][4].getIntensity(), 100)

  TEST_REAL_SIMILAR(e[0].getRT(), 60)
  TEST_REAL_SIMILAR(e[1].getRT(), 120)
  TEST_REAL_SIMILAR(e[2].getRT(), 180)

	//check data processing (all spectra are assigned the general information)
	for (Size i=0; i< e.size(); ++i)
	{
		TEST_EQUAL(e[i].getDataProcessing().size(),1)
		TEST_EQUAL(e[i].getDataProcessing().back().getSoftware().getName(), "17")
		TEST_EQUAL(e[i].getDataProcessing().back().getCompletionTime().get(), "0000-00-00 00:00:00")
		TEST_REAL_SIMILAR(e[i].getDataProcessing().back().getMetaValue("ProcessingNumer"), 123.0)
		TEST_EQUAL(e[i].getDataProcessing().back().getMetaValue("ErrorLog"), "")
		TEST_EQUAL(e[i].getDataProcessing().back().getMetaValue("CalibrationHistory"), "25262728")
		TEST_REAL_SIMILAR(e[i].getDataProcessing().back().getMetaValue("NumOfCalibrations"), 456.0)
	}

  //---------------------------------------------------------------------------
  // RT = 60
  //---------------------------------------------------------------------------
  TEST_EQUAL(e[0].getType(), SpectrumSettings::UNKNOWN)
	TEST_EQUAL(e[0].getInstrumentSettings().getScanMode(), InstrumentSettings::UNKNOWN)
	TEST_EQUAL(e[0].getInstrumentSettings().getPolarity(), IonSource::POSITIVE)
	TEST_EQUAL(e[0].getInstrumentSettings().getScanWindows().size(), 1)
	TEST_REAL_SIMILAR(e[0].getInstrumentSettings().getScanWindows()[0].begin, 0)
	TEST_REAL_SIMILAR(e[0].getInstrumentSettings().getScanWindows()[0].end, 0)
	TEST_EQUAL(e[0].getPrecursors().size(), 0)

	//---------------------------------------------------------------------------
	// RT = 120
	//---------------------------------------------------------------------------
	TEST_EQUAL(e[1].getType(), SpectrumSettings::UNKNOWN)
	TEST_EQUAL(e[1].getInstrumentSettings().getScanMode(), InstrumentSettings::UNKNOWN)
	TEST_EQUAL(e[1].getInstrumentSettings().getPolarity(), IonSource::POSITIVE)
	TEST_EQUAL(e[1].getInstrumentSettings().getScanWindows().size(), 1)
	TEST_REAL_SIMILAR(e[1].getInstrumentSettings().getScanWindows()[0].begin, 0)
	TEST_REAL_SIMILAR(e[1].getInstrumentSettings().getScanWindows()[0].end, 0)
	TEST_EQUAL(e[1].getPrecursors().size(), 0)

	//---------------------------------------------------------------------------
	// RT = 180
	//---------------------------------------------------------------------------
	TEST_EQUAL(e[2].getType(), SpectrumSettings::UNKNOWN)
	TEST_EQUAL(e[2].getInstrumentSettings().getScanMode(), InstrumentSettings::UNKNOWN)
	TEST_EQUAL(e[2].getInstrumentSettings().getPolarity(), IonSource::POSITIVE)
	TEST_EQUAL(e[2].getInstrumentSettings().getScanWindows().size(), 1)
	TEST_REAL_SIMILAR(e[2].getInstrumentSettings().getScanWindows()[0].begin, 0)
	TEST_REAL_SIMILAR(e[2].getInstrumentSettings().getScanWindows()[0].end, 0)
	TEST_EQUAL(e[2].getPrecursors().size(), 0)

  //---------------------------------------------------------------------------
  // const vector<SourceFile>& getSourceFiles() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getSourceFiles().size(),1)
	TEST_EQUAL(e.getSourceFiles()[0].getNameOfFile(), "18")
  TEST_REAL_SIMILAR(e.getSourceFiles()[0].getFileSize(), 0)
  TEST_EQUAL(e.getSourceFiles()[0].getFileType(), "19")

  //---------------------------------------------------------------------------
  // const std::vector<ContactPerson>& getContacts() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getContacts().size(), 2)
  ABORT_IF(e.getContacts().size() != 2);
  TEST_EQUAL(e.getContacts()[0].getLastName(), "15")
  TEST_EQUAL(e.getContacts()[0].getMetaValue("ContactPosition"), "Operator")
  TEST_EQUAL(e.getContacts()[1].getLastName(), "7")
  TEST_EQUAL(e.getContacts()[1].getContactInfo(), "6")
  TEST_EQUAL(e.getContacts()[1].getMetaValue("ContactPosition"), "Dataset owner")

  //---------------------------------------------------------------------------
  // const Instrument& getInstrument() const;
  //---------------------------------------------------------------------------
	const Instrument& inst = e.getInstrument();
  TEST_EQUAL(inst.getName(), "i1")
	TEST_EQUAL(inst.getVendor(), "i3")
	TEST_EQUAL(inst.getModel(), "i4")
	TEST_EQUAL(inst.getMetaValue("InstSerial"), "i5")
	TEST_EQUAL(inst.getMetaValue("InstComments"), "i10")
	TEST_EQUAL(inst.getMetaValue("InstSoftware"), "i6")
	TEST_EQUAL(inst.getMetaValue("InstFirmware"), "i7")
	TEST_EQUAL(inst.getMetaValue("InstOS"), "i8")
	TEST_EQUAL(inst.getMetaValue("InstID"), "i2")
	TEST_EQUAL(inst.getMetaValue("InstParams"), "50")
	TEST_EQUAL(inst.getIonSources().size(),1)
	TEST_EQUAL(inst.getIonSources()[0].getIonizationMethod(), IonSource::EI)
	TEST_EQUAL(inst.getIonSources()[0].getInletType(), IonSource::MEMBRANESEPARATOR)
	TEST_EQUAL(inst.getIonSources()[0].getPolarity(), IonSource::POLNULL)
	TEST_REAL_SIMILAR(inst.getIonSources()[0].getMetaValue("InletTemp"), 2.7)
	TEST_EQUAL(inst.getIonSources()[0].getMetaValue("IonModeAdd"), "FABType=44 FABMatrix=45 ReagentGas=43 ReagentGasPressure=12.3 ElectronEnergy=23.56 LaserWaveLength=56.23 FilamentCurrent=2.3 EmissionCurrent=3.4 ")
	TEST_REAL_SIMILAR(inst.getIonSources()[0].getMetaValue("SrcTemp"), 1.2)
	TEST_REAL_SIMILAR(inst.getIonSources()[0].getMetaValue("AccPot"), 4.5)
	TEST_EQUAL(inst.getIonDetectors().size(),1)
	TEST_EQUAL(inst.getIonDetectors()[0].getType(), IonDetector::ELECTRONMULTIPLIER)
	TEST_EQUAL(inst.getIonDetectors()[0].getAcquisitionMode(), IonDetector::ACQMODENULL)
	TEST_EQUAL(inst.getIonDetectors()[0].getResolution(), 0)
	TEST_EQUAL(inst.getIonDetectors()[0].getADCSamplingFrequency(), 0)
	TEST_REAL_SIMILAR(inst.getIonDetectors()[0].getMetaValue("DetPot"), 5.6)
	TEST_REAL_SIMILAR(inst.getIonDetectors()[0].getMetaValue("DetEntrPot"), 6.7)
  TEST_EQUAL(inst.getMassAnalyzers().size(), 1)
  ABORT_IF(inst.getMassAnalyzers().size() != 1);
	TEST_EQUAL(inst.getMassAnalyzers()[0].getType(), MassAnalyzer::ANALYZERNULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionMethod(), MassAnalyzer::RESMETHNULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionType(), MassAnalyzer::RESTYPENULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanDirection(), MassAnalyzer::UP)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanLaw(), MassAnalyzer::LINEAR)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getReflectronState(), MassAnalyzer::REFLSTATENULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolution(), 0)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getAccuracy(), 0)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanRate(), 0)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanTime(), 12.2f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getTOFTotalPathLength(), 0)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getIsolationWidth(), 0)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getFinalMSExponent(), 0)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getMagneticFieldStrength(), 0)

	//---------------------------------------------------------------------------
  // const HPLC& getHPLC() const;
  //---------------------------------------------------------------------------
	TEST_EQUAL(e.getHPLC().getTemperature(), 21)
	TEST_EQUAL(e.getHPLC().getPressure(), 0)
	TEST_EQUAL(e.getHPLC().getFlux(), 0)

	/////////////////////// TESTING SPECIAL CASES ///////////////////////

	//load a second time to make sure everything is re-initialized correctly
	MSExperiment<> e2;
	file.load(OPENMS_GET_TEST_DATA_PATH("ANDIFile_test.cdf"),e2);
	TEST_EQUAL(e==e2,true)

	//test if it works with different peak types
	MSExperiment<RichPeak1D> e_rich;
  file.load(OPENMS_GET_TEST_DATA_PATH("ANDIFile_test.cdf"),e_rich);

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
