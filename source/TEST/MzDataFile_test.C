// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSExperimentExtern.h>

using namespace OpenMS;
using namespace std;

DRange<1> makeRange(float a, float b)
{
	DPosition<1> pa(a), pb(b);
	return DRange<1>(pa, pb);
}

///////////////////////////

START_TEST(MzDataFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzDataFile* ptr = 0;
CHECK((MzDataFile()))
	ptr = new MzDataFile;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~MzDataFile()))
	delete ptr;
RESULT

CHECK(const PeakFileOptions& getOptions() const)
	MzDataFile file;
	TEST_EQUAL(file.getOptions().hasMSLevels(),false)
RESULT

CHECK(PeakFileOptions& getOptions())
	MzDataFile file;
	file.getOptions().addMSLevel(1);
	TEST_EQUAL(file.getOptions().hasMSLevels(),true);
RESULT

CHECK((template<typename MapType> void load(const String& filename, MapType& map) throw(Exception::FileNotFound, Exception::ParseError)))
	PRECISION(0.01)

  //---------------------------------------------------------------------------
  // test with DRawDataPoint (only peak data is tested, no meta data)
  //---------------------------------------------------------------------------

	MSExperiment< DRawDataPoint<1> > e2;
	MzDataFile mzdata;
	
	//test exception
	TEST_EXCEPTION( Exception::FileNotFound , mzdata.load("dummy/dummy.MzData",e2) )
	
	// real test
	mzdata.load("data/MzDataFile_test_1.mzData",e2);
  //---------------------------------------------------------------------------
  // 60 : (120,100)
  // 120: (110,100) (120,200) (130,100)
  // 180: (100,100) (110,200) (120,300) (130,200) (140,100)
	//--------------------------------------------------------------------------- 
  TEST_EQUAL(e2.size(), 3)

	TEST_EQUAL(e2[0].getContainer().size(), 1)
	TEST_EQUAL(e2[1].getContainer().size(), 3)
	TEST_EQUAL(e2[2].getContainer().size(), 5)

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

  //---------------------------------------------------------------------------
  // test with DPickedPeak
  //---------------------------------------------------------------------------
	MSExperiment< DPickedPeak<1> > e;

	// real test
	mzdata.load("data/MzDataFile_test_1.mzData",e);
  //---------------------------------------------------------------------------
  // 60 : (120,100)
  // 120: (110,100) (120,200) (130,100)
  // 180: (100,100) (110,200) (120,300) (130,200) (140,100) 
	//--------------------------------------------------------------------------- 
  TEST_EQUAL(e.size(), 3)
	TEST_REAL_EQUAL(e[0].getMSLevel(), 1)
	TEST_REAL_EQUAL(e[1].getMSLevel(), 2)
	TEST_REAL_EQUAL(e[2].getMSLevel(), 1)
	TEST_REAL_EQUAL(e[0].getRetentionTime(), 60)
	TEST_REAL_EQUAL(e[1].getRetentionTime(), 120)
	TEST_REAL_EQUAL(e[2].getRetentionTime(), 180)
	TEST_EQUAL(e[0].getType(), SpectrumSettings::UNKNOWN)

	TEST_EQUAL(e[0].getMetaInfoDescriptions()["1"].getSourceFile().getNameOfFile(),"area.raw")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["1"].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["1"].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["1"].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["1"].getMetaValue("Comment"),"Area of the peak")

	TEST_EQUAL(e[0].getMetaInfoDescriptions()["2"].getSourceFile().getNameOfFile(),"fwhm.raw")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["2"].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["2"].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["2"].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["2"].getMetaValue("Comment"),"Full width at half max")

	TEST_EQUAL(e[0].getMetaInfoDescriptions()["3"].getSourceFile().getNameOfFile(),"leftWidth.raw")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["3"].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["3"].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["3"].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["3"].getMetaValue("Comment"),"Left width")

	TEST_EQUAL(e[0].getMetaInfoDescriptions()["4"].getSourceFile().getNameOfFile(),"rightWidth.raw")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["4"].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["4"].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["4"].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["4"].getMetaValue("Comment"),"Right width")

	TEST_EQUAL(e[0].getMetaInfoDescriptions()["5"].getSourceFile().getNameOfFile(),"charge.raw")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["5"].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["5"].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["5"].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["5"].getMetaValue("Comment"),"Peak charge")

	TEST_EQUAL(e[0].getMetaInfoDescriptions()["6"].getSourceFile().getNameOfFile(),"signalToNoise.raw")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["6"].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["6"].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["6"].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["6"].getMetaValue("Comment"),"Signal to noise ratio")

	TEST_EQUAL(e[0].getMetaInfoDescriptions()["7"].getSourceFile().getNameOfFile(),"rValue.raw")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["7"].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["7"].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["7"].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["7"].getMetaValue("Comment"),"Correlation value")

	TEST_EQUAL(e[0].getMetaInfoDescriptions()["8"].getSourceFile().getNameOfFile(),"peakShape.raw")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["8"].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["8"].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["8"].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaInfoDescriptions()["8"].getMetaValue("Comment"),"Peak shape")

	TEST_EQUAL(e[1].getPrecursorPeak().getPosition()[0], 1.2f)
	TEST_EQUAL(e[1].getPrecursorPeak().getCharge(), 2)
	TEST_EQUAL(e[1].getPrecursorPeak().getIntensity(), 2.3f)
	TEST_EQUAL(e[1].getPrecursorPeak().getMetaValue("#IntensityUnits"),
																														"NumberOfCounts")
	TEST_EQUAL(e[1].getPrecursorPeak().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(e[1].getPrecursorPeak().getMetaValue("IonSelectionComment"), "selected")
	TEST_EQUAL(e[1].getPrecursor().getActivationMethod(), Precursor::CID)
	TEST_EQUAL(e[1].getPrecursor().getActivationEnergy(), 3.4f)
	TEST_EQUAL(e[1].getPrecursor().getActivationEnergyUnit(), Precursor::PERCENT)
	TEST_EQUAL(e[1].getPrecursor().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(e[1].getPrecursor().getMetaValue("ActivationComment"), "active")

	TEST_EQUAL(e[0].getInstrumentSettings().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(e[1].getInstrumentSettings().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(e[2].getInstrumentSettings().metaValueExists("URL"), false)
	TEST_EQUAL(e[0].getInstrumentSettings().getMetaValue("SpecComment"), "Spectrum 1")
	TEST_EQUAL(e[1].getInstrumentSettings().getMetaValue("SpecComment"), "Spectrum 2")
	TEST_EQUAL(e[2].getInstrumentSettings().metaValueExists("SpecComment"), false)
	TEST_EQUAL(e[0].getInstrumentSettings().getScanMode(), InstrumentSettings::MASSSCAN)
	TEST_EQUAL(e[1].getInstrumentSettings().getScanMode(), InstrumentSettings::MASSSCAN)
	TEST_EQUAL(e[2].getInstrumentSettings().getScanMode(), InstrumentSettings::SELECTEDIONDETECTION)
	TEST_EQUAL(e[0].getInstrumentSettings().getPolarity(), IonSource::POSITIVE)
	TEST_EQUAL(e[1].getInstrumentSettings().getPolarity(), IonSource::POSITIVE)
	TEST_EQUAL(e[2].getInstrumentSettings().getPolarity(), IonSource::NEGATIVE)
	TEST_REAL_EQUAL(e[0].getInstrumentSettings().getMzRangeStart(), 0)
	TEST_REAL_EQUAL(e[1].getInstrumentSettings().getMzRangeStart(), 110)
	TEST_REAL_EQUAL(e[2].getInstrumentSettings().getMzRangeStart(), 100)
	TEST_REAL_EQUAL(e[0].getInstrumentSettings().getMzRangeStop(), 0)
	TEST_REAL_EQUAL(e[1].getInstrumentSettings().getMzRangeStop(), 130)
	TEST_REAL_EQUAL(e[2].getInstrumentSettings().getMzRangeStop(), 140)


	TEST_EQUAL(e[0].getAcquisitionInfo().size(), 0)
  ABORT_IF(e[0].getAcquisitionInfo().size()!=0);
	TEST_EQUAL(e[1].getAcquisitionInfo().size(), 2)

  ABORT_IF(e[1].getAcquisitionInfo().size()!=2);
	TEST_EQUAL(e[1].getType(), SpectrumSettings::RAWDATA)
	TEST_EQUAL(e[1].getAcquisitionInfo().getMethodOfCombination(), "sum")
	TEST_EQUAL(e[1].getAcquisitionInfo()[0].getNumber(), 501)
	TEST_EQUAL(e[1].getAcquisitionInfo()[1].getNumber(), 502)
	TEST_EQUAL(e[1].getAcquisitionInfo()[0].getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(e[1].getAcquisitionInfo()[1].getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(e[1].getAcquisitionInfo()[0].getMetaValue("AcqComment"), "Acquisition 1")
	TEST_EQUAL(e[1].getAcquisitionInfo()[1].getMetaValue("AcqComment"), "Acquisition 2")

	TEST_EQUAL(e[2].getAcquisitionInfo().size(), 1)
  ABORT_IF(e[2].getAcquisitionInfo().size()!=1);
	TEST_EQUAL(e[2].getType(), SpectrumSettings::PEAKS)
	TEST_EQUAL(e[2].getAcquisitionInfo().getMethodOfCombination(), "average")
	TEST_EQUAL(e[2].getAcquisitionInfo()[0].getNumber(), 601)

	TEST_EQUAL(e[0].getContainer().size(), 1)
	TEST_EQUAL(e[1].getContainer().size(), 3)
	TEST_EQUAL(e[2].getContainer().size(), 5)

	TEST_REAL_EQUAL(e[0].getContainer()[0].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getArea(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getFWHM(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getRightWidthParameter(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getLeftWidthParameter(), 100)
	TEST_EQUAL(e[0].getContainer()[0].getCharge(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getRValue(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getSN(), 100)
	TEST_EQUAL(e[0].getContainer()[0].getPeakShape(), 100)

	TEST_REAL_EQUAL(e[1].getContainer()[0].getPosition()[0], 110)
	TEST_REAL_EQUAL(e[1].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[0].getArea(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[0].getFWHM(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[0].getRightWidthParameter(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[0].getLeftWidthParameter(), 100)
	TEST_EQUAL(e[1].getContainer()[0].getCharge(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[0].getRValue(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[0].getSN(), 100)
	TEST_EQUAL(e[1].getContainer()[0].getPeakShape(), 100)

	TEST_REAL_EQUAL(e[1].getContainer()[1].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[1].getContainer()[1].getIntensity(), 200)
	TEST_REAL_EQUAL(e[1].getContainer()[1].getArea(), 200)
	TEST_REAL_EQUAL(e[1].getContainer()[1].getFWHM(), 200)
	TEST_REAL_EQUAL(e[1].getContainer()[1].getRightWidthParameter(), 200)
	TEST_REAL_EQUAL(e[1].getContainer()[1].getLeftWidthParameter(), 200)
	TEST_EQUAL(e[1].getContainer()[1].getCharge(), 200)
	TEST_REAL_EQUAL(e[1].getContainer()[1].getRValue(), 200)
	TEST_REAL_EQUAL(e[1].getContainer()[1].getSN(), 200)
	TEST_EQUAL(e[1].getContainer()[1].getPeakShape(), 200)

	TEST_REAL_EQUAL(e[1].getContainer()[2].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[1].getContainer()[2].getIntensity(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[2].getArea(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[2].getFWHM(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[2].getRightWidthParameter(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[2].getLeftWidthParameter(), 100)
	TEST_EQUAL(e[1].getContainer()[2].getCharge(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[2].getRValue(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[2].getSN(), 100)
	TEST_EQUAL(e[1].getContainer()[2].getPeakShape(), 100)

	TEST_REAL_EQUAL(e[2].getContainer()[0].getPosition()[0], 100)
	TEST_REAL_EQUAL(e[2].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[0].getArea(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[0].getFWHM(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[0].getRightWidthParameter(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[0].getLeftWidthParameter(), 100)
	TEST_EQUAL(e[2].getContainer()[0].getCharge(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[0].getRValue(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[0].getSN(), 100)
	TEST_EQUAL(e[2].getContainer()[0].getPeakShape(), 100)

	TEST_REAL_EQUAL(e[2].getContainer()[1].getPosition()[0], 110)
	TEST_REAL_EQUAL(e[2].getContainer()[1].getIntensity(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[1].getArea(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[1].getFWHM(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[1].getRightWidthParameter(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[1].getLeftWidthParameter(), 200)
	TEST_EQUAL(e[2].getContainer()[1].getCharge(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[1].getRValue(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[1].getSN(), 200)
	TEST_EQUAL(e[2].getContainer()[1].getPeakShape(), 200)

	TEST_REAL_EQUAL(e[2].getContainer()[2].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[2].getContainer()[2].getIntensity(), 300)
	TEST_REAL_EQUAL(e[2].getContainer()[2].getArea(), 300)
	TEST_REAL_EQUAL(e[2].getContainer()[2].getFWHM(), 300)
	TEST_REAL_EQUAL(e[2].getContainer()[2].getRightWidthParameter(), 300)
	TEST_REAL_EQUAL(e[2].getContainer()[2].getLeftWidthParameter(), 300)
	TEST_EQUAL(e[2].getContainer()[2].getCharge(), 300)
	TEST_REAL_EQUAL(e[2].getContainer()[2].getRValue(), 300)
	TEST_REAL_EQUAL(e[2].getContainer()[2].getSN(), 300)
	TEST_EQUAL(e[2].getContainer()[2].getPeakShape(), 300)

	TEST_REAL_EQUAL(e[2].getContainer()[3].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[2].getContainer()[3].getIntensity(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[3].getArea(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[3].getFWHM(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[3].getRightWidthParameter(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[3].getLeftWidthParameter(), 200)
	TEST_EQUAL(e[2].getContainer()[3].getCharge(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[3].getRValue(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[3].getSN(), 200)
	TEST_EQUAL(e[2].getContainer()[3].getPeakShape(), 200)

	TEST_REAL_EQUAL(e[2].getContainer()[4].getPosition()[0], 140)
	TEST_REAL_EQUAL(e[2].getContainer()[4].getIntensity(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[4].getArea(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[4].getFWHM(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[4].getRightWidthParameter(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[4].getLeftWidthParameter(), 100)
	TEST_EQUAL(e[2].getContainer()[4].getCharge(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[4].getRValue(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[4].getSN(), 100)
	TEST_EQUAL(e[2].getContainer()[4].getPeakShape(), 100)

  //---------------------------------------------------------------------------
  // const SourceFile& getSourceFile() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getSourceFile().getNameOfFile(), "MzDataFile_test_1.raw");
  TEST_EQUAL(e.getSourceFile().getPathToFile(), "/share/data/");
  TEST_EQUAL(e.getSourceFile().getFileType(), "MS");

  //---------------------------------------------------------------------------
  // const std::vector<ContactPerson>& getContacts() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getContacts().size(), 2);
  ABORT_IF(e.getContacts().size()!=2);
  TEST_EQUAL(e.getContacts()[0].getFirstName(), "John");
  TEST_EQUAL(e.getContacts()[0].getLastName(), "Doe");
  TEST_EQUAL(e.getContacts()[0].getInstitution(), "department 1");
  TEST_EQUAL(e.getContacts()[0].getContactInfo(), "www.john.doe");
  TEST_EQUAL(e.getContacts()[1].getFirstName(), "Jane");
  TEST_EQUAL(e.getContacts()[1].getLastName(), "Doe");
  TEST_EQUAL(e.getContacts()[1].getInstitution(), "department 2");
  TEST_EQUAL(e.getContacts()[1].getContactInfo(), "www.jane.doe");

  //---------------------------------------------------------------------------
  // const Software& getSoftware() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getSoftware().getName(), "MS-X");
  TEST_EQUAL(e.getSoftware().getVersion(), "1.0");
  TEST_EQUAL(e.getSoftware().getComment(), "none");
	String tmp;
	e.getSoftware().getCompletionTime().get(tmp);
  TEST_EQUAL(tmp, "2001-02-03 04:05:06");

  //---------------------------------------------------------------------------
  // const ProcessingMethod& getProcessingMethod() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getProcessingMethod().getDeisotoping(), false)
	TEST_EQUAL(e.getProcessingMethod().getChargeDeconvolution(), false)
	TEST_EQUAL(e.getProcessingMethod().getSpectrumType(), SpectrumSettings::PEAKS)
	TEST_EQUAL(e.getProcessingMethod().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(e.getProcessingMethod().getMetaValue("ProcessingComment"), "Processed")

  //---------------------------------------------------------------------------
  // const Instrument& getInstrument() const;
  //---------------------------------------------------------------------------
	const Instrument& inst = e.getInstrument();
  TEST_EQUAL(inst.getName(), "MS-Instrument")
	TEST_EQUAL(inst.getVendor(), "MS-Vendor")
	TEST_EQUAL(inst.getModel(), "MS 1")
	TEST_EQUAL(inst.getCustomizations(), "tuned")
	TEST_EQUAL(inst.getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getMetaValue("AdditionalComment"), "Additional")
	TEST_EQUAL(inst.getIonSource().getIonizationMethod(), IonSource::ESI)
	TEST_EQUAL(inst.getIonSource().getInletType(), IonSource::DIRECT)
	TEST_EQUAL(inst.getIonSource().getPolarity(), IonSource::NEGATIVE) 
	TEST_EQUAL(inst.getIonSource().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getIonSource().getMetaValue("SourceComment"), "Source")
	TEST_EQUAL(inst.getIonDetector().getType(), IonDetector::FARADAYCUP)
	TEST_EQUAL(inst.getIonDetector().getAcquisitionMode(), IonDetector::TDC)
	TEST_EQUAL(inst.getIonDetector().getResolution(), 0.815f)
	TEST_EQUAL(inst.getIonDetector().getADCSamplingFrequency(), 11.22f)
	TEST_EQUAL(inst.getIonDetector().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getIonDetector().getMetaValue("DetectorComment"), "Detector")
  TEST_EQUAL(inst.getMassAnalyzers().size(), 2)
  ABORT_IF(inst.getMassAnalyzers().size()!=2);
	TEST_EQUAL(inst.getMassAnalyzers()[0].getType(), MassAnalyzer::PAULIONTRAP)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionMethod(), MassAnalyzer::FWHM)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionType(), MassAnalyzer::CONSTANT)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanFunction(), MassAnalyzer::MASSSCAN)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanDirection(), MassAnalyzer::UP)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanLaw(), MassAnalyzer::LINEAR)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getTandemScanMethod(), MassAnalyzer::PRECURSORIONSCAN)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getReflectronState(), MassAnalyzer::OFF)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolution(), 22.33f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getAccuracy(), 33.44f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanRate(), 44.55f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanTime(), 55.66f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getTOFTotalPathLength(), 66.77f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getIsolationWidth(), 77.88f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getFinalMSExponent(), 2)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getMagneticFieldStrength(), 88.99f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getMassAnalyzers()[0].getMetaValue("AnalyzerComment"), "Analyzer 1")
	TEST_EQUAL(inst.getMassAnalyzers()[1].getType(), MassAnalyzer::QUADRUPOLE)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getResolutionMethod(), MassAnalyzer::BASELINE)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getResolutionType(), MassAnalyzer::PROPORTIONAL)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanFunction(), MassAnalyzer::SELECTEDIONDETECTION)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanDirection(), MassAnalyzer::DOWN)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanLaw(), MassAnalyzer::EXPONENTIAL)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getTandemScanMethod(), MassAnalyzer::PRODUCTIONSCAN)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getReflectronState(), MassAnalyzer::ON)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getResolution(), 12.3f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getAccuracy(), 13.4f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanRate(), 14.5f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanTime(), 15.6f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getTOFTotalPathLength(), 16.7f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getIsolationWidth(), 17.8f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getFinalMSExponent(), -2)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getMagneticFieldStrength(), 18.9f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getMassAnalyzers()[1].getMetaValue("AnalyzerComment"), "Analyzer 2")

  //---------------------------------------------------------------------------
	// const Sample& getSample()
  //---------------------------------------------------------------------------
	TEST_EQUAL(e.getSample().getName(), "MS-Sample")
	TEST_EQUAL(e.getSample().getNumber(), "0-815")
	TEST_EQUAL(e.getSample().getState(), Sample::GAS)
 	TEST_EQUAL(e.getSample().getMass(), 1.01f)
	TEST_EQUAL(e.getSample().getVolume(), 2.02f)
	TEST_EQUAL(e.getSample().getConcentration(), 3.03f)
	TEST_EQUAL(e.getSample().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(e.getSample().getMetaValue("SampleComment"), "Sample")
RESULT

CHECK([EXTRA] load with DRawDataPoint)
	PRECISION(0.01)

  //---------------------------------------------------------------------------
  // test with DRawDataPoint (only peak data is tested, no meta data)
  //---------------------------------------------------------------------------

	MSExperimentExtern< DRawDataPoint<1> > e2;

	// real test
	MzDataFile().load("data/MzDataFile_test_1.mzData",e2);
  //---------------------------------------------------------------------------
  // 60 : (120,100)
  // 120: (110,100) (120,200) (130,100)
  // 180: (100,100) (110,200) (120,300) (130,200) (140,100)
	//--------------------------------------------------------------------------- 
  TEST_EQUAL(e2.size(), 3)

	TEST_EQUAL(e2[0].getContainer().size(), 1)
	TEST_EQUAL(e2[1].getContainer().size(), 3)
	TEST_EQUAL(e2[2].getContainer().size(), 5)

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

RESULT  


CHECK(([EXTRA] load with metadata-only flag))
	PRECISION(0.01)

  //---------------------------------------------------------------------------
  // test with DRawDataPoint (only peak data is tested, no meta data)
  //---------------------------------------------------------------------------

	MzDataFile mzdata;
	mzdata.getOptions().setMetadataOnly(true);
	
	MSExperiment< DPickedPeak<1> > e;

	// real test
	mzdata.load("data/MzDataFile_test_1.mzData",e);
  
  //---------------------------------------------------------------------------
  // const SourceFile& getSourceFile() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getSourceFile().getNameOfFile(), "MzDataFile_test_1.raw");
  TEST_EQUAL(e.getSourceFile().getPathToFile(), "/share/data/");
  TEST_EQUAL(e.getSourceFile().getFileType(), "MS");

  //---------------------------------------------------------------------------
  // const std::vector<ContactPerson>& getContacts() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getContacts().size(), 2);
  ABORT_IF(e.getContacts().size()!=2);
  TEST_EQUAL(e.getContacts()[0].getFirstName(), "John");
  TEST_EQUAL(e.getContacts()[0].getLastName(), "Doe");
  TEST_EQUAL(e.getContacts()[0].getInstitution(), "department 1");
  TEST_EQUAL(e.getContacts()[0].getContactInfo(), "www.john.doe");
  TEST_EQUAL(e.getContacts()[1].getFirstName(), "Jane");
  TEST_EQUAL(e.getContacts()[1].getLastName(), "Doe");
  TEST_EQUAL(e.getContacts()[1].getInstitution(), "department 2");
  TEST_EQUAL(e.getContacts()[1].getContactInfo(), "www.jane.doe");

  //---------------------------------------------------------------------------
  // const Software& getSoftware() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getSoftware().getName(), "MS-X");
  TEST_EQUAL(e.getSoftware().getVersion(), "1.0");
  TEST_EQUAL(e.getSoftware().getComment(), "none");
	String tmp;
	e.getSoftware().getCompletionTime().get(tmp);
  TEST_EQUAL(tmp, "2001-02-03 04:05:06");

  //---------------------------------------------------------------------------
  // const ProcessingMethod& getProcessingMethod() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getProcessingMethod().getDeisotoping(), false)
	TEST_EQUAL(e.getProcessingMethod().getChargeDeconvolution(), false)
	TEST_EQUAL(e.getProcessingMethod().getSpectrumType(), SpectrumSettings::PEAKS)
	TEST_EQUAL(e.getProcessingMethod().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(e.getProcessingMethod().getMetaValue("ProcessingComment"), "Processed")

  //---------------------------------------------------------------------------
  // const Instrument& getInstrument() const;
  //---------------------------------------------------------------------------
	const Instrument& inst = e.getInstrument();
  TEST_EQUAL(inst.getName(), "MS-Instrument")
	TEST_EQUAL(inst.getVendor(), "MS-Vendor")
	TEST_EQUAL(inst.getModel(), "MS 1")
	TEST_EQUAL(inst.getCustomizations(), "tuned")
	TEST_EQUAL(inst.getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getMetaValue("AdditionalComment"), "Additional")
	TEST_EQUAL(inst.getIonSource().getIonizationMethod(), IonSource::ESI)
	TEST_EQUAL(inst.getIonSource().getInletType(), IonSource::DIRECT)
	TEST_EQUAL(inst.getIonSource().getPolarity(), IonSource::NEGATIVE) 
	TEST_EQUAL(inst.getIonSource().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getIonSource().getMetaValue("SourceComment"), "Source")
	TEST_EQUAL(inst.getIonDetector().getType(), IonDetector::FARADAYCUP)
	TEST_EQUAL(inst.getIonDetector().getAcquisitionMode(), IonDetector::TDC)
	TEST_EQUAL(inst.getIonDetector().getResolution(), 0.815f)
	TEST_EQUAL(inst.getIonDetector().getADCSamplingFrequency(), 11.22f)
	TEST_EQUAL(inst.getIonDetector().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getIonDetector().getMetaValue("DetectorComment"), "Detector")
  TEST_EQUAL(inst.getMassAnalyzers().size(), 2)
  ABORT_IF(inst.getMassAnalyzers().size()!=2);
	TEST_EQUAL(inst.getMassAnalyzers()[0].getType(), MassAnalyzer::PAULIONTRAP)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionMethod(), MassAnalyzer::FWHM)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionType(), MassAnalyzer::CONSTANT)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanFunction(), MassAnalyzer::MASSSCAN)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanDirection(), MassAnalyzer::UP)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanLaw(), MassAnalyzer::LINEAR)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getTandemScanMethod(), MassAnalyzer::PRECURSORIONSCAN)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getReflectronState(), MassAnalyzer::OFF)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolution(), 22.33f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getAccuracy(), 33.44f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanRate(), 44.55f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanTime(), 55.66f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getTOFTotalPathLength(), 66.77f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getIsolationWidth(), 77.88f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getFinalMSExponent(), 2)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getMagneticFieldStrength(), 88.99f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getMassAnalyzers()[0].getMetaValue("AnalyzerComment"), "Analyzer 1")
	TEST_EQUAL(inst.getMassAnalyzers()[1].getType(), MassAnalyzer::QUADRUPOLE)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getResolutionMethod(), MassAnalyzer::BASELINE)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getResolutionType(), MassAnalyzer::PROPORTIONAL)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanFunction(), MassAnalyzer::SELECTEDIONDETECTION)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanDirection(), MassAnalyzer::DOWN)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanLaw(), MassAnalyzer::EXPONENTIAL)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getTandemScanMethod(), MassAnalyzer::PRODUCTIONSCAN)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getReflectronState(), MassAnalyzer::ON)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getResolution(), 12.3f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getAccuracy(), 13.4f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanRate(), 14.5f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanTime(), 15.6f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getTOFTotalPathLength(), 16.7f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getIsolationWidth(), 17.8f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getFinalMSExponent(), -2)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getMagneticFieldStrength(), 18.9f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getMassAnalyzers()[1].getMetaValue("AnalyzerComment"), "Analyzer 2")

  //---------------------------------------------------------------------------
	// const Sample& getSample()
  //---------------------------------------------------------------------------
	TEST_EQUAL(e.getSample().getName(), "MS-Sample")
	TEST_EQUAL(e.getSample().getNumber(), "0-815")
	TEST_EQUAL(e.getSample().getState(), Sample::GAS)
 	TEST_EQUAL(e.getSample().getMass(), 1.01f)
	TEST_EQUAL(e.getSample().getVolume(), 2.02f)
	TEST_EQUAL(e.getSample().getConcentration(), 3.03f)
	TEST_EQUAL(e.getSample().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(e.getSample().getMetaValue("SampleComment"), "Sample")
RESULT

CHECK(([EXTRA] load with selected MS levels))
	PRECISION(0.01)

	MSExperiment< DRawDataPoint<1> > e;
	MzDataFile mzdata;
	
	// load only MS level 1
	mzdata.getOptions().addMSLevel(1);
	mzdata.load("data/MzDataFile_test_1.mzData",e);
	TEST_EQUAL(e.size(), 2)
	TEST_EQUAL(e[0].getContainer().size(), 1)
	TEST_EQUAL(e[1].getContainer().size(), 5)
	TEST_REAL_EQUAL(e[0].getMSLevel(), 1)
	TEST_REAL_EQUAL(e[1].getMSLevel(), 1)
	
	// load all MS levels
	mzdata.getOptions().clearMSLevels();
	mzdata.load("data/MzDataFile_test_1.mzData",e);
	TEST_EQUAL(e.size(), 3)
	TEST_EQUAL(e[0].getContainer().size(), 1)
	TEST_EQUAL(e[1].getContainer().size(), 3)
	TEST_EQUAL(e[2].getContainer().size(), 5)
	TEST_REAL_EQUAL(e[0].getMSLevel(), 1)
	TEST_REAL_EQUAL(e[1].getMSLevel(), 2)
	TEST_REAL_EQUAL(e[2].getMSLevel(), 1)
RESULT

CHECK(([EXTRA] load with RT range))
	PRECISION(0.01)

	MSExperiment< DRawDataPoint<1> > e;
	MzDataFile mzdata;
	
	mzdata.getOptions().setRTRange(makeRange(100, 200));
	mzdata.load("data/MzDataFile_test_1.mzData",e);
	//---------------------------------------------------------------------------
	// 60 : (120,100)
	// 120: (110,100) (120,200) (130,100)
	// 180: (100,100) (110,200) (120,300) (130,200) (140,100) 
	//--------------------------------------------------------------------------- 
	TEST_EQUAL(e.size(), 2)
	TEST_REAL_EQUAL(e[0].getMSLevel(), 2)
	TEST_REAL_EQUAL(e[1].getMSLevel(), 1)
	TEST_REAL_EQUAL(e[0].getRetentionTime(), 120)
	TEST_REAL_EQUAL(e[1].getRetentionTime(), 180)
RESULT

CHECK(([EXTRA] load with MZ range))
	PRECISION(0.01)

	MSExperiment< DRawDataPoint<1> > e;
	MzDataFile mzdata;
	
	mzdata.getOptions().setMZRange(makeRange(115, 135));
	mzdata.load("data/MzDataFile_test_1.mzData",e);
	//---------------------------------------------------------------------------
	// 60 : +(120,100)
	// 120: -(110,100) +(120,200) +(130,100)
	// 180: -(100,100) -(110,200) +(120,300) +(130,200) -(140,100)
	//--------------------------------------------------------------------------- 
	TEST_EQUAL(e.size(), 3)

	TEST_EQUAL(e[0].getContainer().size(), 1)
	TEST_EQUAL(e[1].getContainer().size(), 2)
	TEST_EQUAL(e[2].getContainer().size(), 2)

	TEST_REAL_EQUAL(e[0].getContainer()[0].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getIntensity(), 100)

	TEST_REAL_EQUAL(e[1].getContainer()[0].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[1].getContainer()[0].getIntensity(), 200)

	TEST_REAL_EQUAL(e[1].getContainer()[1].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[1].getContainer()[1].getIntensity(), 100)

	TEST_REAL_EQUAL(e[2].getContainer()[0].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[2].getContainer()[0].getIntensity(), 300)

	TEST_REAL_EQUAL(e[2].getContainer()[1].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[2].getContainer()[1].getIntensity(), 200)
RESULT

CHECK(([EXTRA] load with intensity range))
	PRECISION(0.01)

	MSExperiment< DRawDataPoint<1> > e;
	MzDataFile mzdata;
	
	mzdata.getOptions().setIntensityRange(makeRange(150, 350));
	mzdata.load("data/MzDataFile_test_1.mzData",e);
	//---------------------------------------------------------------------------
	// 60 : -(120,100)
	// 120: -(110,100) +(120,200) -(130,100)
	// 180: -(100,100) +(110,200) +(120,300) +(130,200) -(140,100)
	//--------------------------------------------------------------------------- 
	TEST_EQUAL(e.size(), 3)

	TEST_EQUAL(e[0].getContainer().size(), 0)
	TEST_EQUAL(e[1].getContainer().size(), 1)
	TEST_EQUAL(e[2].getContainer().size(), 3)

	TEST_REAL_EQUAL(e[1].getContainer()[0].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[1].getContainer()[0].getIntensity(), 200)

	TEST_REAL_EQUAL(e[2].getContainer()[0].getPosition()[0], 110)
	TEST_REAL_EQUAL(e[2].getContainer()[0].getIntensity(), 200)

	TEST_REAL_EQUAL(e[2].getContainer()[1].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[2].getContainer()[1].getIntensity(), 300)

	TEST_REAL_EQUAL(e[2].getContainer()[2].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[2].getContainer()[2].getIntensity(), 200)
RESULT

CHECK((template<typename MapType> void store(const String& filename, const MapType& map) const throw(Exception::UnableToCreateFile)))
  MSExperiment< DPickedPeak<1> > e1, e2;
  MzDataFile f;
  f.load("data/MzDataFile_test_1.mzData",e1);
	TEST_EQUAL(e1.size(), 3)

	std::string tmp_filename;
 	NEW_TMP_FILE(tmp_filename);
	f.store(tmp_filename,e1);
	f.load(tmp_filename,e2);
	TEST_EQUAL(e1==e2,true);

	MSExperiment< DRawDataPoint<1> > e3, e4;
	NEW_TMP_FILE(tmp_filename);
	f.load("data/MzDataFile_test_2.mzData",e3);
	f.store(tmp_filename,e3);
	f.load(tmp_filename,e4);
	TEST_EQUAL(e3==e4,true);
RESULT

// check load for 64Bit precision and endian conversion
CHECK([EXTRA] load with 64 bit )
	PRECISION(0.01)

	MSExperiment< DPickedPeak<1> > e;
	MzDataFile mzdata;

	// real test
	mzdata.load("data/MzDataFile_test_3.mzData",e);
  //---------------------------------------------------------------------------
  // 120: (110,100) (120,200) (130,100)
	//---------------------------------------------------------------------------
  TEST_EQUAL(e.size(), 1)

	TEST_EQUAL(e[0].getContainer().size(), 3)

	TEST_REAL_EQUAL(e[0].getContainer()[0].getPosition()[0], 110)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getArea(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getFWHM(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getRightWidthParameter(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getLeftWidthParameter(), 100)
	TEST_EQUAL(e[0].getContainer()[0].getCharge(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getRValue(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getSN(), 100)
	TEST_EQUAL(e[0].getContainer()[0].getPeakShape(), 100)

	TEST_REAL_EQUAL(e[0].getContainer()[1].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[0].getContainer()[1].getIntensity(), 200)
	TEST_REAL_EQUAL(e[0].getContainer()[1].getArea(), 200)
	TEST_REAL_EQUAL(e[0].getContainer()[1].getFWHM(), 200)
	TEST_REAL_EQUAL(e[0].getContainer()[1].getRightWidthParameter(), 200)
	TEST_REAL_EQUAL(e[0].getContainer()[1].getLeftWidthParameter(), 200)
	TEST_EQUAL(e[0].getContainer()[1].getCharge(), 200)
	TEST_REAL_EQUAL(e[0].getContainer()[1].getRValue(), 200)
	TEST_REAL_EQUAL(e[0].getContainer()[1].getSN(), 200)
	TEST_EQUAL(e[0].getContainer()[1].getPeakShape(), 200)

	TEST_REAL_EQUAL(e[0].getContainer()[2].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[0].getContainer()[2].getIntensity(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[2].getArea(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[2].getFWHM(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[2].getRightWidthParameter(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[2].getLeftWidthParameter(), 100)
	TEST_EQUAL(e[0].getContainer()[2].getCharge(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[2].getRValue(), 100)
	TEST_REAL_EQUAL(e[0].getContainer()[2].getSN(), 100)
	TEST_EQUAL(e[0].getContainer()[2].getPeakShape(), 100)
RESULT

// check for Float Kernel traits
CHECK(([EXTRA] load/store for Float Kernel Traits))
	std::string tmp_filename;
	NEW_TMP_FILE(tmp_filename);

  	MzDataFile f;
	MSExperiment< DRawDataPoint<1, FloatKernelTraits> > e1, e2;

	f.load("data/MzDataFile_test_2.mzData",e1);
	f.store(tmp_filename,e1);
	f.load(tmp_filename,e2);
	TEST_EQUAL(e1==e2, true);
RESULT

// check for Float Kernel traits
CHECK(([EXTRA]  load/store for empty scans))
	std::string tmp_filename;
	NEW_TMP_FILE(tmp_filename);

  	MzDataFile f;
	MSExperiment<> e2;
	e2.resize(5);

	f.store(tmp_filename,e2);
	f.load(tmp_filename,e2);
	
	e2.updateRanges();
	
	TEST_EQUAL(e2.size(),5);
	TEST_EQUAL(e2.getSize(),0);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
