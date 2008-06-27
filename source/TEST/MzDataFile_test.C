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

#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

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

CHECK((template<typename MapType> void load(const String& filename, MapType& map) ))
	PRECISION(0.01)

  //---------------------------------------------------------------------------
  // test with DPeak (only peak data is tested, no meta data)
  //---------------------------------------------------------------------------

	MSExperiment< Peak1D > e2;
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
  // test with Peak1D
  //---------------------------------------------------------------------------
	MSExperiment<> e;

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
	TEST_REAL_EQUAL(e[0].getRT(), 60)
	TEST_REAL_EQUAL(e[1].getRT(), 120)
	TEST_REAL_EQUAL(e[2].getRT(), 180)
	TEST_EQUAL(e[0].getType(), SpectrumSettings::UNKNOWN)

	TEST_EQUAL(e[0].getMetaDataArrays()[0].getSourceFile().getNameOfFile(),"area.raw")
	TEST_EQUAL(e[0].getMetaDataArrays()[0].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaDataArrays()[0].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaDataArrays()[0].getComment(),"bla|comment|bla")
	TEST_EQUAL(e[0].getMetaDataArrays()[0].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaDataArrays()[0].getMetaValue("Comment"),"Area of the peak")

	TEST_EQUAL(e[0].getMetaDataArrays()[1].getSourceFile().getNameOfFile(),"fwhm.raw")
	TEST_EQUAL(e[0].getMetaDataArrays()[1].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaDataArrays()[1].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaDataArrays()[1].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaDataArrays()[1].getMetaValue("Comment"),"Full width at half max")

	TEST_EQUAL(e[0].getMetaDataArrays()[2].getSourceFile().getNameOfFile(),"leftWidth.raw")
	TEST_EQUAL(e[0].getMetaDataArrays()[2].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaDataArrays()[2].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaDataArrays()[2].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaDataArrays()[2].getMetaValue("Comment"),"Left width")

	TEST_EQUAL(e[0].getMetaDataArrays()[3].getSourceFile().getNameOfFile(),"rightWidth.raw")
	TEST_EQUAL(e[0].getMetaDataArrays()[3].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaDataArrays()[3].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaDataArrays()[3].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaDataArrays()[3].getMetaValue("Comment"),"Right width")

	TEST_EQUAL(e[0].getMetaDataArrays()[4].getSourceFile().getNameOfFile(),"charge.raw")
	TEST_EQUAL(e[0].getMetaDataArrays()[4].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaDataArrays()[4].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaDataArrays()[4].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaDataArrays()[4].getMetaValue("Comment"),"Peak charge")

	TEST_EQUAL(e[0].getMetaDataArrays()[5].getSourceFile().getNameOfFile(),"signalToNoise.raw")
	TEST_EQUAL(e[0].getMetaDataArrays()[5].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaDataArrays()[5].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaDataArrays()[5].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaDataArrays()[5].getMetaValue("Comment"),"Signal to noise ratio")

	TEST_EQUAL(e[0].getMetaDataArrays()[6].getSourceFile().getNameOfFile(),"rValue.raw")
	TEST_EQUAL(e[0].getMetaDataArrays()[6].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaDataArrays()[6].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaDataArrays()[6].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaDataArrays()[6].getMetaValue("Comment"),"Correlation value")

	TEST_EQUAL(e[0].getMetaDataArrays()[7].getSourceFile().getNameOfFile(),"peakShape.raw")
	TEST_EQUAL(e[0].getMetaDataArrays()[7].getSourceFile().getPathToFile(),"/share/data/")
	TEST_EQUAL(e[0].getMetaDataArrays()[7].getSourceFile().getFileType(),"aux")
	TEST_EQUAL(e[0].getMetaDataArrays()[7].getMetaValue("URL"),"www.open-ms.de")
	TEST_EQUAL(e[0].getMetaDataArrays()[7].getMetaValue("Comment"),"Peak shape")

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
	TEST_EQUAL(e[0].getInstrumentSettings().getScanMode(), InstrumentSettings::FULL)
	TEST_EQUAL(e[1].getInstrumentSettings().getScanMode(), InstrumentSettings::FULL)
	TEST_EQUAL(e[2].getInstrumentSettings().getScanMode(), InstrumentSettings::ZOOM)
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
 
//	0) r_value
//  1) area
//  2) FWHM
//  3) left_width
//  4) right_width
//  5) charge
//  5) type
//  6) signal_to_noise

	TEST_REAL_EQUAL(e[0].getContainer()[0].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[1][0], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[2][0], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[4][0], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[3][0], 100)
	TEST_EQUAL(e[0].getMetaDataArrays()[5][0], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[0][0], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[7][0], 100)
	TEST_EQUAL(e[0].getMetaDataArrays()[6][0], 100)

	TEST_REAL_EQUAL(e[1].getContainer()[0].getPosition()[0], 110)
	TEST_REAL_EQUAL(e[1].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[1][0], 100)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[2][0], 100)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[4][0], 100)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[3][0], 100)
	TEST_EQUAL(e[1].getMetaDataArrays()[5][0], 100)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[0][0], 100)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[7][0], 100)
	TEST_EQUAL(e[1].getMetaDataArrays()[6][0], 100)

	TEST_REAL_EQUAL(e[1].getContainer()[1].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[1].getContainer()[1].getIntensity(), 200)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[1][1], 200)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[2][1], 200)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[4][1], 200)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[3][1], 200)
	TEST_EQUAL(e[1].getMetaDataArrays()[5][1], 200)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[0][1], 200)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[7][1], 200)
	TEST_EQUAL(e[1].getMetaDataArrays()[6][1], 200)

	TEST_REAL_EQUAL(e[1].getContainer()[2].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[1].getContainer()[2].getIntensity(), 100)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[1][2], 100)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[2][2], 100)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[4][2], 100)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[3][2], 100)
	TEST_EQUAL(e[1].getMetaDataArrays()[5][2], 100)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[0][2], 100)
	TEST_REAL_EQUAL(e[1].getMetaDataArrays()[7][2], 100)
	TEST_EQUAL(e[1].getMetaDataArrays()[6][2], 100)

	TEST_REAL_EQUAL(e[2].getContainer()[0].getPosition()[0], 100)
	TEST_REAL_EQUAL(e[2].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[1][0], 100)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[2][0], 100)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[4][0], 100)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[3][0], 100)
	TEST_EQUAL(e[2].getMetaDataArrays()[5][0], 100)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[0][0], 100)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[7][0], 100)
	TEST_EQUAL(e[2].getMetaDataArrays()[6][0], 100)

	TEST_REAL_EQUAL(e[2].getContainer()[1].getPosition()[0], 110)
	TEST_REAL_EQUAL(e[2].getContainer()[1].getIntensity(), 200)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[1][1], 200)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[2][1], 200)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[4][1], 200)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[3][1], 200)
	TEST_EQUAL(e[2].getMetaDataArrays()[5][1], 200)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[0][1], 200)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[7][1], 200)
	TEST_EQUAL(e[2].getMetaDataArrays()[6][1], 200)

	TEST_REAL_EQUAL(e[2].getContainer()[2].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[2].getContainer()[2].getIntensity(), 300)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[1][2], 300)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[2][2], 300)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[4][2], 300)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[3][2], 300)
	TEST_EQUAL(e[2].getMetaDataArrays()[5][2], 300)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[0][2], 300)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[7][2], 300)
	TEST_EQUAL(e[2].getMetaDataArrays()[6][2], 300)

	TEST_REAL_EQUAL(e[2].getContainer()[3].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[2].getContainer()[3].getIntensity(), 200)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[1][3], 200)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[2][3], 200)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[4][3], 200)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[3][3], 200)
	TEST_EQUAL(e[2].getMetaDataArrays()[5][3], 200)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[0][3], 200)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[7][3], 200)
	TEST_EQUAL(e[2].getMetaDataArrays()[6][3], 200)

	TEST_REAL_EQUAL(e[2].getContainer()[4].getPosition()[0], 140)
	TEST_REAL_EQUAL(e[2].getContainer()[4].getIntensity(), 100)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[1][4], 100)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[2][4], 100)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[4][4], 100)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[3][4], 100)
	TEST_EQUAL(e[2].getMetaDataArrays()[5][4], 100)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[0][4], 100)
	TEST_REAL_EQUAL(e[2].getMetaDataArrays()[7][4], 100)
	TEST_EQUAL(e[2].getMetaDataArrays()[6][4], 100)

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

CHECK([EXTRA] load with DPeak)
	PRECISION(0.01)

  //---------------------------------------------------------------------------
  // test with DPeak (only peak data is tested, no meta data)
  //---------------------------------------------------------------------------

	MSExperiment< Peak1D > e2;

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
  // test with DPeak (only peak data is tested, no meta data)
  //---------------------------------------------------------------------------

	MzDataFile mzdata;
	mzdata.getOptions().setMetadataOnly(true);
	
	MSExperiment<> e;

	// real test
	mzdata.load("data/MzDataFile_test_1.mzData",e);

	//check number of scans
	TEST_EQUAL(e.size(),0)

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

	MSExperiment< Peak1D > e;
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

	MSExperiment< Peak1D > e;
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
	TEST_REAL_EQUAL(e[0].getRT(), 120)
	TEST_REAL_EQUAL(e[1].getRT(), 180)
RESULT

CHECK(([EXTRA] load with MZ range))
	PRECISION(0.01)

	MSExperiment< Peak1D > e;
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

	MSExperiment< Peak1D > e;
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

CHECK(([EXTRA] load one extremely long spectrum - tests CDATA splitting))
	MSExperiment< Peak1D > e;
	MzDataFile mzdata;
	
	mzdata.load("data/MzDataFile_test_4.mzData",e);
	TEST_EQUAL(e.size(), 1)
	TEST_EQUAL(e[0].size(), 997530)
RESULT

CHECK((template<typename MapType> void store(const String& filename, const MapType& map) const ))
  MSExperiment<> e1, e2;
  MzDataFile f;
  f.load("data/MzDataFile_test_1.mzData",e1);
	TEST_EQUAL(e1.size(), 3)

	std::string tmp_filename;
 	NEW_TMP_FILE(tmp_filename);
	f.store(tmp_filename,e1);
	f.load(tmp_filename,e2);
  TEST_EQUAL(e1==e2,true);
	TEST_EQUAL(e1[0].getMetaDataArrays()==e1[0].getMetaDataArrays(),true);
	TEST_EQUAL(e1[1].getMetaDataArrays()==e1[1].getMetaDataArrays(),true);
	TEST_EQUAL(e1[2].getMetaDataArrays()==e1[2].getMetaDataArrays(),true);
	NEW_TMP_FILE(tmp_filename);
	f.store(tmp_filename,e2);
	
	MSExperiment< Peak1D > e3, e4;
	NEW_TMP_FILE(tmp_filename);
	f.load("data/MzDataFile_test_2.mzData",e3);
	f.store(tmp_filename,e3);
	f.load(tmp_filename,e4);
	TEST_EQUAL(e3==e4,true);
RESULT

CHECK([EXTRA] load with 64 bit precision and endian conversion )
	PRECISION(0.01)

	MSExperiment<> e;
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
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[1][0], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[2][0], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[4][0], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[3][0], 100)
	TEST_EQUAL(e[0].getMetaDataArrays()[5][0], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[0][0], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[7][0], 100)
	TEST_EQUAL(e[0].getMetaDataArrays()[6][0], 100)

	TEST_REAL_EQUAL(e[0].getContainer()[1].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[0].getContainer()[1].getIntensity(), 200)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[1][1], 200)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[2][1], 200)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[4][1], 200)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[3][1], 200)
	TEST_EQUAL(e[0].getMetaDataArrays()[5][1], 200)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[0][1], 200)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[7][1], 200)
	TEST_EQUAL(e[0].getMetaDataArrays()[6][1], 200)

	TEST_REAL_EQUAL(e[0].getContainer()[2].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[0].getContainer()[2].getIntensity(), 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[1][2], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[2][2], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[4][2], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[3][2], 100)
	TEST_EQUAL(e[0].getMetaDataArrays()[5][2], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[0][2], 100)
	TEST_REAL_EQUAL(e[0].getMetaDataArrays()[7][2], 100)
	TEST_EQUAL(e[0].getMetaDataArrays()[6][2], 100)
RESULT

CHECK([EXTRA] static bool isValid(const String& filename))
	std::string tmp_filename;
  MzDataFile f;
  MSExperiment<> e;
  
  //test if empty file is valid
	NEW_TMP_FILE(tmp_filename);
  f.store(tmp_filename,e);
  TEST_EQUAL(f.isValid(tmp_filename),true);

	//test if fill file is valid
	NEW_TMP_FILE(tmp_filename);
	f.load("data/MzDataFile_test_1.mzData",e);
  f.store(tmp_filename,e);
  TEST_EQUAL(f.isValid(tmp_filename),true);
RESULT

CHECK([EXTRA] storing/loading of meta data arrays)
	MzDataFile file;
	//init spectrum/experiment/meta data array
	MSExperiment<> exp;
	MSSpectrum<> spec;
	spec.resize(5);
	spec[0].setIntensity(1.0); spec[0].setMZ(1.0);
	spec[1].setIntensity(2.0); spec[1].setMZ(2.0);
	spec[2].setIntensity(3.0); spec[2].setMZ(3.0);
	spec[3].setIntensity(4.0); spec[3].setMZ(4.0);
	spec[4].setIntensity(5.0); spec[4].setMZ(5.0);
	MSSpectrum<>::MetaDataArray mda1;
	mda1.push_back(1.1);
	mda1.push_back(1.2);
	mda1.push_back(1.3);
	mda1.push_back(1.4);
	mda1.push_back(1.5);
	MSSpectrum<>::MetaDataArray mda2;
	mda2.push_back(-2.1);
	mda2.push_back(-2.2);
	mda2.push_back(-2.3);
	mda2.push_back(-2.4);
	mda2.push_back(-2.5);
	
	//spectrum 1 (one meta data arrays)
	spec.setRT(500.0);
	spec.getMetaDataArrays().push_back(mda1);
	spec.getMetaDataArrays()[0].setName("MDA1");
	exp.push_back(spec);
	
	//spectrum 2 (zero meta data array)
	spec.setRT(600.0);
	spec.getMetaDataArrays().clear();
	exp.push_back(spec);

	//spectrum 3 (two meta data array)
	spec.setRT(700.0);
	spec.getMetaDataArrays().push_back(mda1);
	spec.getMetaDataArrays().push_back(mda2);
	spec.getMetaDataArrays()[0].setName("MDA1");
	spec.getMetaDataArrays()[1].setName("MDA2");
	exp.push_back(spec);
	
	//*******************************************
	//store file
	std::string filename;
 	NEW_TMP_FILE(filename);
 	cout << "Filename: " << filename << std::endl;
	file.store(filename,exp);

	//*******************************************
	//load and check file
	MSExperiment<> exp2;
	file.load(filename,exp2);
	
	TEST_EQUAL(exp2.size(),3)
	TEST_EQUAL(exp2[0].getMetaDataArrays().size(),1)
	TEST_EQUAL(exp2[1].getMetaDataArrays().size(),0)
	TEST_EQUAL(exp2[2].getMetaDataArrays().size(),2)
	
	TEST_EQUAL(exp2[0].getMetaDataArrays()[0].getName(),"MDA1");
	TEST_EQUAL(exp2[2].getMetaDataArrays()[0].getName(),"MDA1");
	TEST_EQUAL(exp2[2].getMetaDataArrays()[1].getName(),"MDA2");
	
	TEST_EQUAL(exp2[0].getMetaDataArrays()[0].size(),5)
	TEST_REAL_EQUAL(exp2[0].getMetaDataArrays()[0][0],1.1);
	TEST_REAL_EQUAL(exp2[0].getMetaDataArrays()[0][1],1.2);
	TEST_REAL_EQUAL(exp2[0].getMetaDataArrays()[0][2],1.3);
	TEST_REAL_EQUAL(exp2[0].getMetaDataArrays()[0][3],1.4);
	TEST_REAL_EQUAL(exp2[0].getMetaDataArrays()[0][4],1.5);
	
	TEST_EQUAL(exp2[2].getMetaDataArrays()[0].size(),5)
	TEST_REAL_EQUAL(exp2[2].getMetaDataArrays()[0][0],1.1);
	TEST_REAL_EQUAL(exp2[2].getMetaDataArrays()[0][1],1.2);
	TEST_REAL_EQUAL(exp2[2].getMetaDataArrays()[0][2],1.3);
	TEST_REAL_EQUAL(exp2[2].getMetaDataArrays()[0][3],1.4);
	TEST_REAL_EQUAL(exp2[2].getMetaDataArrays()[0][4],1.5);
	
	TEST_EQUAL(exp2[2].getMetaDataArrays()[1].size(),5)
	TEST_REAL_EQUAL(exp2[2].getMetaDataArrays()[1][0],-2.1);
	TEST_REAL_EQUAL(exp2[2].getMetaDataArrays()[1][1],-2.2);
	TEST_REAL_EQUAL(exp2[2].getMetaDataArrays()[1][2],-2.3);
	TEST_REAL_EQUAL(exp2[2].getMetaDataArrays()[1][3],-2.4);
	TEST_REAL_EQUAL(exp2[2].getMetaDataArrays()[1][4],-2.5);
	
	//*******************************************	
	//check if filtering of meta data arrays works
	MSExperiment<> exp3;
	file.getOptions().setMZRange(DRange< 1 >(2.5,7.0));
	file.load(filename,exp3);
	
	TEST_EQUAL(exp.size(),3)
	TEST_EQUAL(exp3[0].size(),3)
	TEST_EQUAL(exp3[1].size(),3)
	TEST_EQUAL(exp3[2].size(),3)

	TEST_EQUAL(exp3[0].getMetaDataArrays().size(),1)
	TEST_EQUAL(exp3[1].getMetaDataArrays().size(),0)
	TEST_EQUAL(exp3[2].getMetaDataArrays().size(),2)
	
	TEST_EQUAL(exp3[0].getMetaDataArrays()[0].getName(),"MDA1");
	TEST_EQUAL(exp3[2].getMetaDataArrays()[0].getName(),"MDA1");
	TEST_EQUAL(exp3[2].getMetaDataArrays()[1].getName(),"MDA2");

	TEST_EQUAL(exp3[0].getMetaDataArrays()[0].size(),3)	
	TEST_REAL_EQUAL(exp3[0].getMetaDataArrays()[0][0],1.3);
	TEST_REAL_EQUAL(exp3[0].getMetaDataArrays()[0][1],1.4);
	TEST_REAL_EQUAL(exp3[0].getMetaDataArrays()[0][2],1.5)
	
	TEST_EQUAL(exp3[2].getMetaDataArrays()[0].size(),3)
	TEST_REAL_EQUAL(exp3[2].getMetaDataArrays()[0][0],1.3);
	TEST_REAL_EQUAL(exp3[2].getMetaDataArrays()[0][1],1.4);
	TEST_REAL_EQUAL(exp3[2].getMetaDataArrays()[0][2],1.5);
	
	TEST_EQUAL(exp3[2].getMetaDataArrays()[1].size(),3)
	TEST_REAL_EQUAL(exp3[2].getMetaDataArrays()[1][0],-2.3);
	TEST_REAL_EQUAL(exp3[2].getMetaDataArrays()[1][1],-2.4);
	TEST_REAL_EQUAL(exp3[2].getMetaDataArrays()[1][2],-2.5);

  //*********************************************
  //test if the storing meta data arrays without a name works
  
	exp3[0].getMetaDataArrays()[0].setName("");
	exp3[2].getMetaDataArrays()[0].setName("");
	exp3[2].getMetaDataArrays()[1].setName("");

	MSExperiment<> exp4;
	file.store(filename,exp3);
	file.load(filename,exp4);
	
	TEST_EQUAL(exp.size(),3)
	TEST_EQUAL(exp4[0].size(),3)
	TEST_EQUAL(exp4[1].size(),3)
	TEST_EQUAL(exp4[2].size(),3)

	TEST_EQUAL(exp4[0].getMetaDataArrays().size(),1)
	TEST_EQUAL(exp4[1].getMetaDataArrays().size(),0)
	TEST_EQUAL(exp4[2].getMetaDataArrays().size(),2)
	
	TEST_EQUAL(exp4[0].getMetaDataArrays()[0].getName(),"");
	TEST_EQUAL(exp4[2].getMetaDataArrays()[0].getName(),"");
	TEST_EQUAL(exp4[2].getMetaDataArrays()[1].getName(),"");

	TEST_EQUAL(exp4[0].getMetaDataArrays()[0].size(),3)	
	TEST_REAL_EQUAL(exp4[0].getMetaDataArrays()[0][0],1.3);
	TEST_REAL_EQUAL(exp4[0].getMetaDataArrays()[0][1],1.4);
	TEST_REAL_EQUAL(exp4[0].getMetaDataArrays()[0][2],1.5)
	
	TEST_EQUAL(exp4[2].getMetaDataArrays()[0].size(),3)
	TEST_REAL_EQUAL(exp4[2].getMetaDataArrays()[0][0],1.3);
	TEST_REAL_EQUAL(exp4[2].getMetaDataArrays()[0][1],1.4);
	TEST_REAL_EQUAL(exp4[2].getMetaDataArrays()[0][2],1.5);
	
	TEST_EQUAL(exp4[2].getMetaDataArrays()[1].size(),3)
	TEST_REAL_EQUAL(exp4[2].getMetaDataArrays()[1][0],-2.3);
	TEST_REAL_EQUAL(exp4[2].getMetaDataArrays()[1][1],-2.4);
	TEST_REAL_EQUAL(exp4[2].getMetaDataArrays()[1][2],-2.5);
	

RESULT

CHECK([[EXTRA] loading a minimal file containing one spectrum  - with whitespaces inside the base64 data)
	MSExperiment<> e;
  MzDataFile f;
  f.load("data/MzDataFile_test_5.mzData",e);
  TEST_EQUAL(e.size(),1)
  TEST_EQUAL(e[0].size(),3)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
