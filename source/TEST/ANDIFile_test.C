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

#include <OpenMS/FORMAT/ANDIFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;

///////////////////////////

START_TEST(ANDIFile, "$Id$")

/////////////////////////////////////////////////////////////

ANDIFile* ptr = 0;
CHECK((ANDIFile()))
	ptr = new ANDIFile;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~ANDIFile()))
	delete ptr;
RESULT

CHECK((template<typename MapType> void load(const String& filename, MapType& map) ))
	PRECISION(0.01)

	MSExperiment< Peak1D > e2;
	ANDIFile andi;

	//test exception
	TEST_EXCEPTION( Exception::FileNotFound , andi.load("dummy/dummy.cdf",e2) )

	// real test
	andi.load("data/ANDIFile_test.cdf",e2);
  //---------------------------------------------------------------------------
  // 60 : (120,100) 
  // 120: (110,100) (120,200) (130,100)
	// 180: (100,100) (110,200) (120,300) (130,200) (140,100) 
	//---------------------------------------------------------------------------
  TEST_EQUAL(e2.size(), 3)
	TEST_REAL_EQUAL(e2[0].getMSLevel(), 1)
	TEST_REAL_EQUAL(e2[1].getMSLevel(), 1)
	TEST_REAL_EQUAL(e2[2].getMSLevel(), 1)
	TEST_REAL_EQUAL(e2[0].getRT(), 60)
	TEST_REAL_EQUAL(e2[1].getRT(), 120)
	TEST_REAL_EQUAL(e2[2].getRT(), 180)
	TEST_REAL_EQUAL(e2[0].getContainer().size(), 1)
	TEST_REAL_EQUAL(e2[1].getContainer().size(), 3)
	TEST_REAL_EQUAL(e2[2].getContainer().size(), 5)

	TEST_REAL_EQUAL(e2[0].getContainer()[0].getPosition()[0], 120)
	TEST_REAL_EQUAL(e2[0].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e2[1].getContainer()[0].getPosition()[0], 110)
	TEST_REAL_EQUAL(e2[1].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e2[1].getContainer()[1].getPosition()[0], 120)
	TEST_REAL_EQUAL(e2[1].getContainer()[1].getIntensity(), 200)
	TEST_REAL_EQUAL(e2[1].getContainer()[2].getPosition()[0], 130)
	TEST_REAL_EQUAL(e2[1].getContainer()[2].getIntensity(), 100)
	TEST_REAL_EQUAL(e2[2].getContainer()[0].getPosition()[0], 100)
	TEST_REAL_EQUAL(e2[2].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e2[2].getContainer()[1].getPosition()[0], 110)
	TEST_REAL_EQUAL(e2[2].getContainer()[1].getIntensity(), 200)
	TEST_REAL_EQUAL(e2[2].getContainer()[2].getPosition()[0], 120)
	TEST_REAL_EQUAL(e2[2].getContainer()[2].getIntensity(), 300)
	TEST_REAL_EQUAL(e2[2].getContainer()[3].getPosition()[0], 130)
	TEST_REAL_EQUAL(e2[2].getContainer()[3].getIntensity(), 200)
	TEST_REAL_EQUAL(e2[2].getContainer()[4].getPosition()[0], 140)
	TEST_REAL_EQUAL(e2[2].getContainer()[4].getIntensity(), 100)
	
	// meta data:
	
	MSExperiment<> e;
  andi.load("data/ANDIFile_test.cdf",e);
  
  TEST_EQUAL(e.size(), 3)
  
  TEST_REAL_EQUAL(e[0].getRT(), 60)
  TEST_REAL_EQUAL(e[1].getRT(), 120)
  TEST_REAL_EQUAL(e[2].getRT(), 180)
  
  //---------------------------------------------------------------------------
  // RT = 60
  //---------------------------------------------------------------------------
  TEST_EQUAL(e[0].getType(), SpectrumSettings::UNKNOWN)
	TEST_EQUAL(e[0].getInstrumentSettings().getScanMode(), InstrumentSettings::UNKNOWN)
	TEST_EQUAL(e[0].getInstrumentSettings().getPolarity(), IonSource::POSITIVE)
	TEST_REAL_EQUAL(e[0].getInstrumentSettings().getMzRangeStart(), 0)
	TEST_REAL_EQUAL(e[0].getInstrumentSettings().getMzRangeStop(), 0)
	TEST_EQUAL(e[0].getPrecursor().getActivationMethod(), Precursor::ACTMETHNULL)
	TEST_EQUAL(e[0].getPrecursor().getActivationEnergy(), 0)
	TEST_EQUAL(e[0].getPrecursor().getActivationEnergyUnit(), Precursor::UNITSNULL)
	TEST_EQUAL(e[0].getPrecursor().getWindowSize(), 0)
	
	//---------------------------------------------------------------------------
	// RT = 120
	//---------------------------------------------------------------------------
	TEST_EQUAL(e[1].getType(), SpectrumSettings::UNKNOWN)
	TEST_EQUAL(e[1].getInstrumentSettings().getScanMode(), InstrumentSettings::UNKNOWN)
	TEST_EQUAL(e[1].getInstrumentSettings().getPolarity(), IonSource::POSITIVE)
	TEST_REAL_EQUAL(e[1].getInstrumentSettings().getMzRangeStart(), 0)
	TEST_REAL_EQUAL(e[1].getInstrumentSettings().getMzRangeStop(), 0)
	TEST_EQUAL(e[1].getPrecursor().getActivationMethod(), Precursor::ACTMETHNULL)
	TEST_EQUAL(e[1].getPrecursor().getActivationEnergy(), 0)
	TEST_EQUAL(e[1].getPrecursor().getActivationEnergyUnit(), Precursor::UNITSNULL)
	TEST_EQUAL(e[1].getPrecursor().getWindowSize(), 0)
	
	//---------------------------------------------------------------------------
	// RT = 180
	//---------------------------------------------------------------------------
	TEST_EQUAL(e[2].getType(), SpectrumSettings::UNKNOWN)
	TEST_EQUAL(e[2].getInstrumentSettings().getScanMode(), InstrumentSettings::UNKNOWN)
	TEST_EQUAL(e[2].getInstrumentSettings().getPolarity(), IonSource::POSITIVE)
	TEST_REAL_EQUAL(e[2].getInstrumentSettings().getMzRangeStart(), 0)
	TEST_REAL_EQUAL(e[2].getInstrumentSettings().getMzRangeStop(), 0)
	TEST_EQUAL(e[2].getPrecursor().getActivationMethod(), Precursor::ACTMETHNULL)
	TEST_EQUAL(e[2].getPrecursor().getActivationEnergy(), 0)
	TEST_EQUAL(e[2].getPrecursor().getActivationEnergyUnit(), Precursor::UNITSNULL)
	TEST_EQUAL(e[2].getPrecursor().getWindowSize(), 0)
  
  //---------------------------------------------------------------------------
  // const SourceFile& getSourceFile() const;
  //---------------------------------------------------------------------------
	TEST_EQUAL(e.getSourceFile().getNameOfFile(), "18")
  TEST_REAL_EQUAL(e.getSourceFile().getFileSize(), 0)
  TEST_EQUAL(e.getSourceFile().getFileType(), "19")

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
  // const Software& getSoftware() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getSoftware().getName(), "17")
  String tmp;
  e.getSoftware().getCompletionTime().get(tmp);
  TEST_EQUAL(tmp, "0000-00-00 00:00:00")
  
  //---------------------------------------------------------------------------
  // const ProcessingMethod& getProcessingMethod() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getProcessingMethod().getDeisotoping(), false)
  TEST_EQUAL(e.getProcessingMethod().getChargeDeconvolution(), false)
  TEST_EQUAL(e.getProcessingMethod().getSpectrumType(), SpectrumSettings::PEAKS)
  TEST_REAL_EQUAL(e.getProcessingMethod().getIntensityCutoff(), 0)
  TEST_REAL_EQUAL(e.getProcessingMethod().getMetaValue("ProcessingNumer"), 123.0)
  TEST_EQUAL(e.getProcessingMethod().getMetaValue("ErrorLog"), "")
  TEST_EQUAL(e.getProcessingMethod().getMetaValue("CalibrationHistory"), "25262728")
  TEST_REAL_EQUAL(e.getProcessingMethod().getMetaValue("NumOfCalibrations"), 456.0)
  
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
	TEST_EQUAL(inst.getIonSource().getIonizationMethod(), IonSource::EI)
	TEST_EQUAL(inst.getIonSource().getInletType(), IonSource::MEMBRANESEPARATOR)
	TEST_EQUAL(inst.getIonSource().getPolarity(), IonSource::POLNULL)
	TEST_REAL_EQUAL(inst.getIonSource().getMetaValue("InletTemp"), 2.7)
	TEST_EQUAL(inst.getIonSource().getMetaValue("IonModeAdd"), "FABType=44 FABMatrix=45 ReagentGas=43 ReagentGasPressure=12.3 ElectronEnergy=23.56 LaserWaveLength=56.23 FilamentCurrent=2.3 EmissionCurrent=3.4 ")
	TEST_REAL_EQUAL(inst.getIonSource().getMetaValue("SrcTemp"), 1.2)
	TEST_REAL_EQUAL(inst.getIonSource().getMetaValue("AccPot"), 4.5)
	TEST_EQUAL(inst.getIonDetector().getType(), IonDetector::ELECTRONMULTIPLIER)
	TEST_EQUAL(inst.getIonDetector().getAcquisitionMode(), IonDetector::ACQMODENULL)
	TEST_EQUAL(inst.getIonDetector().getResolution(), 0)
	TEST_EQUAL(inst.getIonDetector().getADCSamplingFrequency(), 0)
	TEST_REAL_EQUAL(inst.getIonDetector().getMetaValue("DetPot"), 5.6)
	TEST_REAL_EQUAL(inst.getIonDetector().getMetaValue("DetEntrPot"), 6.7)
  TEST_EQUAL(inst.getMassAnalyzers().size(), 1)
  ABORT_IF(inst.getMassAnalyzers().size() != 1);
	TEST_EQUAL(inst.getMassAnalyzers()[0].getType(), MassAnalyzer::ANALYZERNULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionMethod(), MassAnalyzer::RESMETHNULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionType(), MassAnalyzer::RESTYPENULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanFunction(), MassAnalyzer::MASSSCAN)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanDirection(), MassAnalyzer::UP)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanLaw(), MassAnalyzer::LINEAR)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getTandemScanMethod(), MassAnalyzer::TANDEMNULL)
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
  
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
