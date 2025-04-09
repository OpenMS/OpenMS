// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
///////////////////////////

#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

DRange<1> makeRange(double a, double b)
{
  DPosition<1> pa(a), pb(b);
  return DRange<1>(pa, pb);
}

///////////////////////////

START_TEST(MzDataFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzDataFile * ptr = nullptr;
MzDataFile* nullPointer = nullptr;
START_SECTION((MzDataFile()))
{
ptr = new MzDataFile;
TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~MzDataFile()))
{
delete ptr;
}
END_SECTION

START_SECTION(const PeakFileOptions& getOptions() const)
{
  MzDataFile file;
  TEST_EQUAL(file.getOptions().hasMSLevels(), false)

  const PeakFileOptions o = file.getOptions();
  TEST_EQUAL(o.hasMSLevels(), false)
}
END_SECTION

START_SECTION(setOptions(const PeakFileOptions & options))
{
  MzDataFile file;
  TEST_EQUAL(file.getOptions().hasMSLevels(), false)

  const PeakFileOptions o = file.getOptions();
  PeakFileOptions options = o;
  options.addMSLevel(1);

  file.setOptions(options);
  TEST_EQUAL(file.getOptions().hasMSLevels(), true)
}
END_SECTION

START_SECTION(PeakFileOptions& getOptions())
{
  MzDataFile file;
  file.getOptions().addMSLevel(1);
  TEST_EQUAL(file.getOptions().hasMSLevels(), true);
}
END_SECTION

START_SECTION((template <typename MapType> void load(const String &filename, MapType & map)))
{
  TOLERANCE_ABSOLUTE(0.01)

  MzDataFile file;
  PeakMap e;

  // real test
  file.load(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"), e);

  //test DocumentIdentifier addition
  TEST_STRING_EQUAL(e.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"));
  TEST_STRING_EQUAL(FileTypes::typeToName(e.getLoadedFileType()), "mzData");

  //---------------------------------------------------------------------------
  // ms-level, RT, native ID
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.size(), 3)
  TEST_EQUAL(e[0].getMSLevel(), 1)
  TEST_EQUAL(e[1].getMSLevel(), 2)
  TEST_EQUAL(e[2].getMSLevel(), 1)
  TEST_REAL_SIMILAR(e[0].getRT(), 60)
  TEST_REAL_SIMILAR(e[1].getRT(), 120)
  TEST_REAL_SIMILAR(e[2].getRT(), 180)
  TEST_STRING_EQUAL(e[0].getNativeID(), "spectrum=10")
  TEST_STRING_EQUAL(e[1].getNativeID(), "spectrum=11")
  TEST_STRING_EQUAL(e[2].getNativeID(), "spectrum=12")
  TEST_EQUAL(e[0].getType(), SpectrumSettings::UNKNOWN)

  //---------------------------------------------------------------------------
  //meta data array meta data
  //---------------------------------------------------------------------------
  TEST_EQUAL(e[0].getFloatDataArrays()[0].getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(e[0].getFloatDataArrays()[0].getMetaValue("Comment"), "Area of the peak")
  TEST_EQUAL(e[0].getFloatDataArrays()[0].getMetaValue("comment"), "bla|comment|bla")

  TEST_EQUAL(e[0].getFloatDataArrays()[1].getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(e[0].getFloatDataArrays()[1].getMetaValue("Comment"), "Full width at half max")

  TEST_EQUAL(e[0].getFloatDataArrays()[2].getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(e[0].getFloatDataArrays()[2].getMetaValue("Comment"), "Left width")

  TEST_EQUAL(e[0].getFloatDataArrays()[3].getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(e[0].getFloatDataArrays()[3].getMetaValue("Comment"), "Right width")

  TEST_EQUAL(e[0].getFloatDataArrays()[4].getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(e[0].getFloatDataArrays()[4].getMetaValue("Comment"), "Peak charge")

  TEST_EQUAL(e[0].getFloatDataArrays()[5].getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(e[0].getFloatDataArrays()[5].getMetaValue("Comment"), "Signal to noise ratio")

  TEST_EQUAL(e[0].getFloatDataArrays()[6].getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(e[0].getFloatDataArrays()[6].getMetaValue("Comment"), "Correlation value")

  TEST_EQUAL(e[0].getFloatDataArrays()[7].getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(e[0].getFloatDataArrays()[7].getMetaValue("Comment"), "Peak shape")

  //---------------------------------------------------------------------------
  //precursors
  //---------------------------------------------------------------------------
  TEST_EQUAL(e[0].getPrecursors().size(), 0)
  TEST_EQUAL(e[1].getPrecursors().size(), 2)
  TEST_EQUAL(e[2].getPrecursors().size(), 0)

  TEST_REAL_SIMILAR(e[1].getPrecursors()[0].getMZ(), 1.2)
  TEST_EQUAL(e[1].getPrecursors()[0].getCharge(), 2)
  TEST_REAL_SIMILAR(e[1].getPrecursors()[0].getIntensity(), 2.3f)
  TEST_EQUAL(e[1].getPrecursors()[0].getMetaValue("IonSelectionComment"), "selected")
  TEST_EQUAL(e[1].getPrecursors()[0].getActivationMethods().count(Precursor::CID), 1)
  TEST_REAL_SIMILAR(e[1].getPrecursors()[0].getActivationEnergy(), 3.4)
  TEST_EQUAL(e[1].getPrecursors()[0].getMetaValue("ActivationComment"), "active")

  TEST_REAL_SIMILAR(e[1].getPrecursors()[1].getMZ(), 2.2)
  TEST_EQUAL(e[1].getPrecursors()[1].getCharge(), 3)
  TEST_REAL_SIMILAR(e[1].getPrecursors()[1].getIntensity(), 3.3f)
  TEST_EQUAL(e[1].getPrecursors()[1].getMetaValue("IonSelectionComment"), "selected2")
  TEST_EQUAL(e[1].getPrecursors()[1].getActivationMethods().count(Precursor::SID), 1)
  TEST_REAL_SIMILAR(e[1].getPrecursors()[1].getActivationEnergy(), 4.4)
  TEST_EQUAL(e[1].getPrecursors()[1].getMetaValue("ActivationComment"), "active2")

  //---------------------------------------------------------------------------
  //instrument settings
  //---------------------------------------------------------------------------
  TEST_EQUAL(e[0].getInstrumentSettings().getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(e[1].getInstrumentSettings().getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(e[2].getInstrumentSettings().metaValueExists("URL"), false)
  TEST_EQUAL(e[0].getInstrumentSettings().getMetaValue("SpecComment"), "Spectrum 1")
  TEST_EQUAL(e[1].getInstrumentSettings().getMetaValue("SpecComment"), "Spectrum 2")
  TEST_EQUAL(e[2].getInstrumentSettings().metaValueExists("SpecComment"), false)
  TEST_EQUAL(e[0].getInstrumentSettings().getScanMode(), InstrumentSettings::MASSSPECTRUM)
  TEST_EQUAL(e[1].getInstrumentSettings().getScanMode(), InstrumentSettings::MASSSPECTRUM)
  TEST_EQUAL(e[2].getInstrumentSettings().getScanMode(), InstrumentSettings::SIM)
  TEST_EQUAL(e[0].getInstrumentSettings().getPolarity(), IonSource::POSITIVE)
  TEST_EQUAL(e[1].getInstrumentSettings().getPolarity(), IonSource::POSITIVE)
  TEST_EQUAL(e[2].getInstrumentSettings().getPolarity(), IonSource::NEGATIVE)
  TEST_EQUAL(e[0].getInstrumentSettings().getScanWindows().size(), 0)
  TEST_EQUAL(e[1].getInstrumentSettings().getScanWindows().size(), 1)
  TEST_REAL_SIMILAR(e[1].getInstrumentSettings().getScanWindows()[0].begin, 110)
  TEST_REAL_SIMILAR(e[1].getInstrumentSettings().getScanWindows()[0].end, 0)
  TEST_EQUAL(e[2].getInstrumentSettings().getScanWindows().size(), 1)
  TEST_REAL_SIMILAR(e[2].getInstrumentSettings().getScanWindows()[0].begin, 100)
  TEST_REAL_SIMILAR(e[2].getInstrumentSettings().getScanWindows()[0].end, 140)

  //---------------------------------------------------------------------------
  //acquisition
  //---------------------------------------------------------------------------
  TEST_EQUAL(e[0].getAcquisitionInfo().size(), 0)
  ABORT_IF(!e[0].getAcquisitionInfo().empty());
  TEST_EQUAL(e[1].getAcquisitionInfo().size(), 2)

  ABORT_IF(e[1].getAcquisitionInfo().size() != 2);
  TEST_EQUAL(e[1].getType(), SpectrumSettings::PROFILE)
  TEST_EQUAL(e[1].getAcquisitionInfo().getMethodOfCombination(), "sum")
  TEST_EQUAL(e[1].getAcquisitionInfo()[0].getIdentifier(), "501")
  TEST_EQUAL(e[1].getAcquisitionInfo()[1].getIdentifier(), "502")
  TEST_EQUAL(e[1].getAcquisitionInfo()[0].getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(e[1].getAcquisitionInfo()[1].getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(e[1].getAcquisitionInfo()[0].getMetaValue("AcqComment"), "Acquisition 1")
  TEST_EQUAL(e[1].getAcquisitionInfo()[1].getMetaValue("AcqComment"), "Acquisition 2")

  TEST_EQUAL(e[2].getAcquisitionInfo().size(), 1)
  ABORT_IF(e[2].getAcquisitionInfo().size() != 1);
  TEST_EQUAL(e[2].getType(), SpectrumSettings::CENTROID)
  TEST_EQUAL(e[2].getAcquisitionInfo().getMethodOfCombination(), "average")
  TEST_EQUAL(e[2].getAcquisitionInfo()[0].getIdentifier(), "601")

  //---------------------------------------------------------------------------
  // actual peak data:
  // 60 : (120,100)
  // 120: (110,100) (120,200) (130,100)
  // 180: (100,100) (110,200) (120,300) (130,200) (140,100)
  //
  // meta data array values:
  // 0) r_value
  // 1) area
  // 2) FWHM
  // 3) left_width
  // 4) right_width
  // 5) charge
  // 5) type
  // 6) signal_to_noise
  //---------------------------------------------------------------------------
  TEST_EQUAL(e[0].size(), 1)
  TEST_EQUAL(e[1].size(), 3)
  TEST_EQUAL(e[2].size(), 5)

  TEST_REAL_SIMILAR(e[0][0].getPosition()[0], 120)
  TEST_REAL_SIMILAR(e[0][0].getIntensity(), 100)
  TEST_REAL_SIMILAR(e[0].getFloatDataArrays()[1][0], 100)
  TEST_REAL_SIMILAR(e[0].getFloatDataArrays()[2][0], 100)
  TEST_REAL_SIMILAR(e[0].getFloatDataArrays()[4][0], 100)
  TEST_REAL_SIMILAR(e[0].getFloatDataArrays()[3][0], 100)
  TEST_EQUAL(e[0].getFloatDataArrays()[5][0], 100)
  TEST_REAL_SIMILAR(e[0].getFloatDataArrays()[0][0], 100)
  TEST_REAL_SIMILAR(e[0].getFloatDataArrays()[7][0], 100)
  TEST_EQUAL(e[0].getFloatDataArrays()[6][0], 100)

  TEST_REAL_SIMILAR(e[1][0].getPosition()[0], 110)
  TEST_REAL_SIMILAR(e[1][0].getIntensity(), 100)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[1][0], 100)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[2][0], 100)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[4][0], 100)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[3][0], 100)
  TEST_EQUAL(e[1].getFloatDataArrays()[5][0], 100)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[0][0], 100)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[7][0], 100)
  TEST_EQUAL(e[1].getFloatDataArrays()[6][0], 100)

  TEST_REAL_SIMILAR(e[1][1].getPosition()[0], 120)
  TEST_REAL_SIMILAR(e[1][1].getIntensity(), 200)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[1][1], 200)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[2][1], 200)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[4][1], 200)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[3][1], 200)
  TEST_EQUAL(e[1].getFloatDataArrays()[5][1], 200)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[0][1], 200)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[7][1], 200)
  TEST_EQUAL(e[1].getFloatDataArrays()[6][1], 200)

  TEST_REAL_SIMILAR(e[1][2].getPosition()[0], 130)
  TEST_REAL_SIMILAR(e[1][2].getIntensity(), 100)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[1][2], 100)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[2][2], 100)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[4][2], 100)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[3][2], 100)
  TEST_EQUAL(e[1].getFloatDataArrays()[5][2], 100)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[0][2], 100)
  TEST_REAL_SIMILAR(e[1].getFloatDataArrays()[7][2], 100)
  TEST_EQUAL(e[1].getFloatDataArrays()[6][2], 100)

  TEST_REAL_SIMILAR(e[2][0].getPosition()[0], 100)
  TEST_REAL_SIMILAR(e[2][0].getIntensity(), 100)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[1][0], 100)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[2][0], 100)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[4][0], 100)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[3][0], 100)
  TEST_EQUAL(e[2].getFloatDataArrays()[5][0], 100)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[0][0], 100)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[7][0], 100)
  TEST_EQUAL(e[2].getFloatDataArrays()[6][0], 100)

  TEST_REAL_SIMILAR(e[2][1].getPosition()[0], 110)
  TEST_REAL_SIMILAR(e[2][1].getIntensity(), 200)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[1][1], 200)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[2][1], 200)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[4][1], 200)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[3][1], 200)
  TEST_EQUAL(e[2].getFloatDataArrays()[5][1], 200)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[0][1], 200)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[7][1], 200)
  TEST_EQUAL(e[2].getFloatDataArrays()[6][1], 200)

  TEST_REAL_SIMILAR(e[2][2].getPosition()[0], 120)
  TEST_REAL_SIMILAR(e[2][2].getIntensity(), 300)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[1][2], 300)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[2][2], 300)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[4][2], 300)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[3][2], 300)
  TEST_EQUAL(e[2].getFloatDataArrays()[5][2], 300)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[0][2], 300)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[7][2], 300)
  TEST_EQUAL(e[2].getFloatDataArrays()[6][2], 300)

  TEST_REAL_SIMILAR(e[2][3].getPosition()[0], 130)
  TEST_REAL_SIMILAR(e[2][3].getIntensity(), 200)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[1][3], 200)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[2][3], 200)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[4][3], 200)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[3][3], 200)
  TEST_EQUAL(e[2].getFloatDataArrays()[5][3], 200)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[0][3], 200)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[7][3], 200)
  TEST_EQUAL(e[2].getFloatDataArrays()[6][3], 200)

  TEST_REAL_SIMILAR(e[2][4].getPosition()[0], 140)
  TEST_REAL_SIMILAR(e[2][4].getIntensity(), 100)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[1][4], 100)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[2][4], 100)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[4][4], 100)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[3][4], 100)
  TEST_EQUAL(e[2].getFloatDataArrays()[5][4], 100)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[0][4], 100)
  TEST_REAL_SIMILAR(e[2].getFloatDataArrays()[7][4], 100)
  TEST_EQUAL(e[2].getFloatDataArrays()[6][4], 100)

  //---------------------------------------------------------------------------
  // accession number
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getIdentifier(), "lsid");

  //---------------------------------------------------------------------------
  // source file
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getSourceFiles().size(), 1)
  TEST_STRING_EQUAL(e.getSourceFiles()[0].getNameOfFile(), "MzDataFile_test_1.raw");
  TEST_STRING_EQUAL(e.getSourceFiles()[0].getPathToFile(), "/share/data/");
  TEST_STRING_EQUAL(e.getSourceFiles()[0].getFileType(), "MS");
  TEST_STRING_EQUAL(e.getSourceFiles()[0].getChecksum(), "");
  TEST_EQUAL(e.getSourceFiles()[0].getChecksumType(), SourceFile::UNKNOWN_CHECKSUM);

  //---------------------------------------------------------------------------
  // conteact list
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getContacts().size(), 2);
  ABORT_IF(e.getContacts().size() != 2);
  TEST_EQUAL(e.getContacts()[0].getFirstName(), "John");
  TEST_EQUAL(e.getContacts()[0].getLastName(), "Doe");
  TEST_EQUAL(e.getContacts()[0].getInstitution(), "department 1");
  TEST_EQUAL(e.getContacts()[0].getContactInfo(), "www.john.doe");
  TEST_EQUAL(e.getContacts()[1].getFirstName(), "Jane");
  TEST_EQUAL(e.getContacts()[1].getLastName(), "Doe");
  TEST_EQUAL(e.getContacts()[1].getInstitution(), "department 2");
  TEST_EQUAL(e.getContacts()[1].getContactInfo(), "www.jane.doe");

  //---------------------------------------------------------------------------
  // data processing
  //---------------------------------------------------------------------------
  for (Size i = 0; i < e.size(); ++i)
  {
    TEST_EQUAL(e[i].getDataProcessing().size(), 1)
    TEST_EQUAL(e[i].getDataProcessing()[0]->getMetaValue("URL"), "www.open-ms.de")
    TEST_EQUAL(e[i].getDataProcessing()[0]->getMetaValue("comment"), "ProcessingComment")
    TEST_EQUAL(e[i].getDataProcessing()[0]->getCompletionTime().get(), "2001-02-03 04:05:06");

    TEST_EQUAL(e[i].getDataProcessing()[0]->getSoftware().getName(), "MS-X");
    TEST_EQUAL(e[i].getDataProcessing()[0]->getSoftware().getVersion(), "1.0");
    TEST_EQUAL(e[i].getDataProcessing()[0]->getSoftware().getMetaValue("comment"), "SoftwareComment")
  }
  //---------------------------------------------------------------------------
  // instrument
  //---------------------------------------------------------------------------
  const Instrument& inst = e.getInstrument();
  TEST_EQUAL(inst.getName(), "MS-Instrument")
  TEST_EQUAL(inst.getVendor(), "MS-Vendor")
  TEST_EQUAL(inst.getModel(), "MS 1")
  TEST_EQUAL(inst.getCustomizations(), "tuned")
  TEST_EQUAL(inst.getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(inst.getMetaValue("AdditionalComment"), "Additional")
  TEST_EQUAL(inst.getIonSources().size(), 1)
  TEST_EQUAL(inst.getIonSources()[0].getIonizationMethod(), IonSource::ESI)
  TEST_EQUAL(inst.getIonSources()[0].getInletType(), IonSource::DIRECT)
  TEST_EQUAL(inst.getIonSources()[0].getPolarity(), IonSource::NEGATIVE)
  TEST_EQUAL(inst.getIonSources()[0].getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(inst.getIonSources()[0].getMetaValue("SourceComment"), "Source")
  TEST_EQUAL(inst.getIonDetectors().size(), 1)
  TEST_EQUAL(inst.getIonDetectors()[0].getType(), IonDetector::FARADAYCUP)
  TEST_EQUAL(inst.getIonDetectors()[0].getAcquisitionMode(), IonDetector::TDC)
  TEST_EQUAL(inst.getIonDetectors()[0].getResolution(), 0.815)
  TEST_EQUAL(inst.getIonDetectors()[0].getADCSamplingFrequency(), 11.22)
  TEST_EQUAL(inst.getIonDetectors()[0].getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(inst.getIonDetectors()[0].getMetaValue("DetectorComment"), "Detector")
  TEST_EQUAL(inst.getMassAnalyzers().size(), 2)
  ABORT_IF(inst.getMassAnalyzers().size() != 2);
  TEST_EQUAL(inst.getMassAnalyzers()[0].getType(), MassAnalyzer::PAULIONTRAP)
  TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionMethod(), MassAnalyzer::FWHM)
  TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionType(), MassAnalyzer::CONSTANT)
  TEST_EQUAL(inst.getMassAnalyzers()[0].getScanDirection(), MassAnalyzer::UP)
  TEST_EQUAL(inst.getMassAnalyzers()[0].getScanLaw(), MassAnalyzer::LINEAR)
  TEST_EQUAL(inst.getMassAnalyzers()[0].getReflectronState(), MassAnalyzer::OFF)
  TEST_EQUAL(inst.getMassAnalyzers()[0].getResolution(), 22.33)
  TEST_EQUAL(inst.getMassAnalyzers()[0].getAccuracy(), 33.44)
  TEST_EQUAL(inst.getMassAnalyzers()[0].getScanRate(), 44.55)
  TEST_EQUAL(inst.getMassAnalyzers()[0].getScanTime(), 55.66)
  TEST_EQUAL(inst.getMassAnalyzers()[0].getTOFTotalPathLength(), 66.77)
  TEST_EQUAL(inst.getMassAnalyzers()[0].getIsolationWidth(), 77.88)
  TEST_EQUAL(inst.getMassAnalyzers()[0].getFinalMSExponent(), 2)
  TEST_EQUAL(inst.getMassAnalyzers()[0].getMagneticFieldStrength(), 88.99)
  TEST_EQUAL(inst.getMassAnalyzers()[0].getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(inst.getMassAnalyzers()[0].getMetaValue("AnalyzerComment"), "Analyzer 1")
  TEST_EQUAL(inst.getMassAnalyzers()[1].getType(), MassAnalyzer::QUADRUPOLE)
  TEST_EQUAL(inst.getMassAnalyzers()[1].getResolutionMethod(), MassAnalyzer::BASELINE)
  TEST_EQUAL(inst.getMassAnalyzers()[1].getResolutionType(), MassAnalyzer::PROPORTIONAL)
  TEST_EQUAL(inst.getMassAnalyzers()[1].getScanDirection(), MassAnalyzer::DOWN)
  TEST_EQUAL(inst.getMassAnalyzers()[1].getScanLaw(), MassAnalyzer::EXPONENTIAL)
  TEST_EQUAL(inst.getMassAnalyzers()[1].getReflectronState(), MassAnalyzer::ON)
  TEST_EQUAL(inst.getMassAnalyzers()[1].getResolution(), 12.3)
  TEST_EQUAL(inst.getMassAnalyzers()[1].getAccuracy(), 13.4)
  TEST_EQUAL(inst.getMassAnalyzers()[1].getScanRate(), 14.5)
  TEST_EQUAL(inst.getMassAnalyzers()[1].getScanTime(), 15.6)
  TEST_EQUAL(inst.getMassAnalyzers()[1].getTOFTotalPathLength(), 16.7)
  TEST_EQUAL(inst.getMassAnalyzers()[1].getIsolationWidth(), 17.8)
  TEST_EQUAL(inst.getMassAnalyzers()[1].getFinalMSExponent(), -2)
  TEST_EQUAL(inst.getMassAnalyzers()[1].getMagneticFieldStrength(), 18.9)
  TEST_EQUAL(inst.getMassAnalyzers()[1].getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(inst.getMassAnalyzers()[1].getMetaValue("AnalyzerComment"), "Analyzer 2")

  //---------------------------------------------------------------------------
  // sample
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getSample().getName(), "MS-Sample")
  TEST_EQUAL(e.getSample().getNumber(), "0-815")
  TEST_EQUAL(e.getSample().getState(), Sample::GAS)
  TEST_EQUAL(e.getSample().getMass(), 1.01)
  TEST_EQUAL(e.getSample().getVolume(), 2.02)
  TEST_EQUAL(e.getSample().getConcentration(), 3.03)
  TEST_EQUAL(e.getSample().getMetaValue("URL"), "www.open-ms.de")
  TEST_EQUAL(e.getSample().getMetaValue("SampleComment"), "Sample")

  /////////////////////// TESTING SPECIAL CASES ///////////////////////

  //load a second time to make sure everything is re-initialized correctly
  PeakMap e2;
  file.load(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"), e2);
  TEST_TRUE(e == e2)

  //loading a minimal file containing one spectrum  - with whitespaces inside the base64 data
  PeakMap e3;
  file.load(OPENMS_GET_TEST_DATA_PATH("MzDataFile_3_minimal.mzData"), e3);
  TEST_EQUAL(e3.size(), 1)
  TEST_EQUAL(e3[0].size(), 3)

  //load one extremely long spectrum - tests CDATA splitting
  PeakMap e4;
  file.load(OPENMS_GET_TEST_DATA_PATH("MzDataFile_2_long.mzData"), e4);
  TEST_EQUAL(e4.size(), 1)
  TEST_EQUAL(e4[0].size(), 997530)

  //load with 64 bit precision and endian conversion
  PeakMap e5;
  file.load(OPENMS_GET_TEST_DATA_PATH("MzDataFile_4_64bit.mzData"), e5);
  TEST_EQUAL(e5.getIdentifier(), "");
  TEST_EQUAL(e5.size(), 1)
  TEST_EQUAL(e5[0].size(), 3)
  TEST_REAL_SIMILAR(e5[0][0].getPosition()[0], 110)
  TEST_REAL_SIMILAR(e5[0][0].getIntensity(), 100)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[1][0], 100)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[2][0], 100)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[4][0], 100)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[3][0], 100)
  TEST_EQUAL(e5[0].getFloatDataArrays()[5][0], 100)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[0][0], 100)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[7][0], 100)
  TEST_EQUAL(e5[0].getFloatDataArrays()[6][0], 100)

  TEST_REAL_SIMILAR(e5[0][1].getPosition()[0], 120)
  TEST_REAL_SIMILAR(e5[0][1].getIntensity(), 200)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[1][1], 200)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[2][1], 200)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[4][1], 200)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[3][1], 200)
  TEST_EQUAL(e5[0].getFloatDataArrays()[5][1], 200)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[0][1], 200)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[7][1], 200)
  TEST_EQUAL(e5[0].getFloatDataArrays()[6][1], 200)

  TEST_REAL_SIMILAR(e5[0][2].getPosition()[0], 130)
  TEST_REAL_SIMILAR(e5[0][2].getIntensity(), 100)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[1][2], 100)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[2][2], 100)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[4][2], 100)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[3][2], 100)
  TEST_EQUAL(e5[0].getFloatDataArrays()[5][2], 100)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[0][2], 100)
  TEST_REAL_SIMILAR(e5[0].getFloatDataArrays()[7][2], 100)
  TEST_EQUAL(e5[0].getFloatDataArrays()[6][2], 100)
}
END_SECTION

START_SECTION(([EXTRA] load with metadata - only flag))
{
  TOLERANCE_ABSOLUTE(0.01)

  MzDataFile file;
  file.getOptions().setMetadataOnly(true);

  PeakMap e;

  // real test
  file.load(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"), e);

  //check number of scans
  TEST_EQUAL(e.size(), 0)

  TEST_EQUAL(e.getSourceFiles().size(), 1)
  TEST_STRING_EQUAL(e.getSourceFiles()[0].getNameOfFile(), "MzDataFile_test_1.raw");
  TEST_EQUAL(e.getContacts().size(), 2);
  TEST_EQUAL(e.getContacts()[0].getFirstName(), "John");
  TEST_EQUAL(e.getContacts()[0].getLastName(), "Doe");
  TEST_EQUAL(e.getInstrument().getName(), "MS-Instrument")
  TEST_EQUAL(e.getInstrument().getVendor(), "MS-Vendor")
  TEST_EQUAL(e.getSample().getName(), "MS-Sample")
  TEST_EQUAL(e.getSample().getNumber(), "0-815")
}
END_SECTION

START_SECTION(([EXTRA] load with selected MS levels))
{
  TOLERANCE_ABSOLUTE(0.01)

  PeakMap e;
  MzDataFile file;

  // load only MS level 1
  file.getOptions().addMSLevel(1);
  file.load(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"), e);
  TEST_EQUAL(e.size(), 2)
  TEST_EQUAL(e[0].size(), 1)
  TEST_STRING_EQUAL(e[0].getNativeID(), "spectrum=10")
  TEST_EQUAL(e[1].size(), 5)
  TEST_STRING_EQUAL(e[1].getNativeID(), "spectrum=12")
  TEST_EQUAL(e[0].getMSLevel(), 1)
  TEST_EQUAL(e[1].getMSLevel(), 1)

  // load all MS levels
  file.getOptions().clearMSLevels();
  file.load(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"), e);
  TEST_EQUAL(e.size(), 3)
  TEST_EQUAL(e[0].size(), 1)
  TEST_EQUAL(e[1].size(), 3)
  TEST_EQUAL(e[2].size(), 5)
  TEST_EQUAL(e[0].getMSLevel(), 1)
  TEST_EQUAL(e[1].getMSLevel(), 2)
  TEST_EQUAL(e[2].getMSLevel(), 1)
}
END_SECTION

START_SECTION(([EXTRA] load with RT range))
{
  TOLERANCE_ABSOLUTE(0.01)

  PeakMap e;
  MzDataFile file;

  file.getOptions().setRTRange(makeRange(100, 200));
  file.load(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"), e);
  //---------------------------------------------------------------------------
  // 60 : (120,100)
  // 120: (110,100) (120,200) (130,100)
  // 180: (100,100) (110,200) (120,300) (130,200) (140,100)
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.size(), 2)
  TEST_EQUAL(e[0].getMSLevel(), 2)
  TEST_EQUAL(e[1].getMSLevel(), 1)
  TEST_REAL_SIMILAR(e[0].getRT(), 120)
  TEST_REAL_SIMILAR(e[1].getRT(), 180)
}
END_SECTION

START_SECTION(([EXTRA] load with MZ range))
{
  TOLERANCE_ABSOLUTE(0.01)

  PeakMap e;
  MzDataFile file;

  file.getOptions().setMZRange(makeRange(115, 135));
  file.load(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"), e);
  //---------------------------------------------------------------------------
  // 60 : +(120,100)
  // 120: -(110,100) +(120,200) +(130,100)
  // 180: -(100,100) -(110,200) +(120,300) +(130,200) -(140,100)
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.size(), 3)

  TEST_EQUAL(e[0].size(), 1)
  TEST_EQUAL(e[1].size(), 2)
  TEST_EQUAL(e[2].size(), 2)

  TEST_REAL_SIMILAR(e[0][0].getPosition()[0], 120)
  TEST_REAL_SIMILAR(e[0][0].getIntensity(), 100)

  TEST_REAL_SIMILAR(e[1][0].getPosition()[0], 120)
  TEST_REAL_SIMILAR(e[1][0].getIntensity(), 200)

  TEST_REAL_SIMILAR(e[1][1].getPosition()[0], 130)
  TEST_REAL_SIMILAR(e[1][1].getIntensity(), 100)

  TEST_REAL_SIMILAR(e[2][0].getPosition()[0], 120)
  TEST_REAL_SIMILAR(e[2][0].getIntensity(), 300)

  TEST_REAL_SIMILAR(e[2][1].getPosition()[0], 130)
  TEST_REAL_SIMILAR(e[2][1].getIntensity(), 200)
}
END_SECTION

START_SECTION(([EXTRA] load with intensity range))
{
  TOLERANCE_ABSOLUTE(0.01)

  PeakMap e;
  MzDataFile file;

  file.getOptions().setIntensityRange(makeRange(150, 350));
  file.load(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"), e);
  //---------------------------------------------------------------------------
  // 60 : -(120,100)
  // 120: -(110,100) +(120,200) -(130,100)
  // 180: -(100,100) +(110,200) +(120,300) +(130,200) -(140,100)
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.size(), 3)

  TEST_EQUAL(e[0].size(), 0)
  TEST_EQUAL(e[1].size(), 1)
  TEST_EQUAL(e[2].size(), 3)

  TEST_REAL_SIMILAR(e[1][0].getPosition()[0], 120)
  TEST_REAL_SIMILAR(e[1][0].getIntensity(), 200)

  TEST_REAL_SIMILAR(e[2][0].getPosition()[0], 110)
  TEST_REAL_SIMILAR(e[2][0].getIntensity(), 200)

  TEST_REAL_SIMILAR(e[2][1].getPosition()[0], 120)
  TEST_REAL_SIMILAR(e[2][1].getIntensity(), 300)

  TEST_REAL_SIMILAR(e[2][2].getPosition()[0], 130)
  TEST_REAL_SIMILAR(e[2][2].getIntensity(), 200)
}
END_SECTION

START_SECTION((template <typename MapType> void store(const String &filename, const MapType &map) const))
{
  PeakMap e1, e2;
  MzDataFile f;
  f.load(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"), e1);
  TEST_EQUAL(e1.size(), 3)

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  f.store(tmp_filename, e1);
  f.load(tmp_filename, e2);
  TEST_EQUAL(e2.getIdentifier(), "lsid");
  e2[0].getDataProcessing()[0]->getSoftware().setMetaValue("comment", String("SoftwareComment"));
  e2[1].getDataProcessing()[0]->getSoftware().setMetaValue("comment", String("SoftwareComment"));
  e2[2].getDataProcessing()[0]->getSoftware().setMetaValue("comment", String("SoftwareComment"));
  TEST_TRUE(e1 == e2);
}
END_SECTION

START_SECTION([EXTRA] storing / loading of meta data arrays)
{
  MzDataFile file;
  //init spectrum/experiment/meta data array
  PeakMap exp;
  MSSpectrum spec;
  spec.resize(5);
  spec[0].setIntensity(1.0f); spec[0].setMZ(1.0);
  spec[1].setIntensity(2.0f); spec[1].setMZ(2.0);
  spec[2].setIntensity(3.0f); spec[2].setMZ(3.0);
  spec[3].setIntensity(4.0f); spec[3].setMZ(4.0);
  spec[4].setIntensity(5.0f); spec[4].setMZ(5.0);
  MSSpectrum::FloatDataArray mda1;
  mda1.push_back(1.1f);
  mda1.push_back(1.2f);
  mda1.push_back(1.3f);
  mda1.push_back(1.4f);
  mda1.push_back(1.5f);
  MSSpectrum::FloatDataArray mda2;
  mda2.push_back(-2.1f);
  mda2.push_back(-2.2f);
  mda2.push_back(-2.3f);
  mda2.push_back(-2.4f);
  mda2.push_back(-2.5f);

  //spectrum 1 (one meta data arrays)
  spec.setRT(500.0);
  spec.getFloatDataArrays().push_back(mda1);
  spec.getFloatDataArrays()[0].setName("MDA1");
  exp.addSpectrum(spec);

  //spectrum 2 (zero meta data array)
  spec.setRT(600.0);
  spec.getFloatDataArrays().clear();
  exp.addSpectrum(spec);

  //spectrum 3 (two meta data array)
  spec.setRT(700.0);
  spec.getFloatDataArrays().push_back(mda1);
  spec.getFloatDataArrays().push_back(mda2);
  spec.getFloatDataArrays()[0].setName("MDA1");
  spec.getFloatDataArrays()[1].setName("MDA2");
  exp.addSpectrum(spec);

  //*******************************************
  //store file
  std::string filename;
  NEW_TMP_FILE(filename);
  cout << "Filename: " << filename << std::endl;
  file.store(filename, exp);

  //*******************************************
  //load and check file
  PeakMap exp2;
  file.load(filename, exp2);

  TEST_EQUAL(exp2.size(), 3)
  TEST_EQUAL(exp2[0].getFloatDataArrays().size(), 1)
  TEST_EQUAL(exp2[1].getFloatDataArrays().size(), 0)
  TEST_EQUAL(exp2[2].getFloatDataArrays().size(), 2)

  TEST_EQUAL(exp2[0].getFloatDataArrays()[0].getName(), "MDA1");
  TEST_EQUAL(exp2[2].getFloatDataArrays()[0].getName(), "MDA1");
  TEST_EQUAL(exp2[2].getFloatDataArrays()[1].getName(), "MDA2");

  TEST_EQUAL(exp2[0].getFloatDataArrays()[0].size(), 5)
  TEST_REAL_SIMILAR(exp2[0].getFloatDataArrays()[0][0], 1.1);
  TEST_REAL_SIMILAR(exp2[0].getFloatDataArrays()[0][1], 1.2);
  TEST_REAL_SIMILAR(exp2[0].getFloatDataArrays()[0][2], 1.3);
  TEST_REAL_SIMILAR(exp2[0].getFloatDataArrays()[0][3], 1.4);
  TEST_REAL_SIMILAR(exp2[0].getFloatDataArrays()[0][4], 1.5);

  TEST_EQUAL(exp2[2].getFloatDataArrays()[0].size(), 5)
  TEST_REAL_SIMILAR(exp2[2].getFloatDataArrays()[0][0], 1.1);
  TEST_REAL_SIMILAR(exp2[2].getFloatDataArrays()[0][1], 1.2);
  TEST_REAL_SIMILAR(exp2[2].getFloatDataArrays()[0][2], 1.3);
  TEST_REAL_SIMILAR(exp2[2].getFloatDataArrays()[0][3], 1.4);
  TEST_REAL_SIMILAR(exp2[2].getFloatDataArrays()[0][4], 1.5);

  TEST_EQUAL(exp2[2].getFloatDataArrays()[1].size(), 5)
  TEST_REAL_SIMILAR(exp2[2].getFloatDataArrays()[1][0], -2.1);
  TEST_REAL_SIMILAR(exp2[2].getFloatDataArrays()[1][1], -2.2);
  TEST_REAL_SIMILAR(exp2[2].getFloatDataArrays()[1][2], -2.3);
  TEST_REAL_SIMILAR(exp2[2].getFloatDataArrays()[1][3], -2.4);
  TEST_REAL_SIMILAR(exp2[2].getFloatDataArrays()[1][4], -2.5);

  //*******************************************
  //check if filtering of meta data arrays works
  PeakMap exp3;
  file.getOptions().setMZRange(DRange<1>(2.5, 7.0));
  file.load(filename, exp3);

  TEST_EQUAL(exp.size(), 3)
  TEST_EQUAL(exp3[0].size(), 3)
  TEST_EQUAL(exp3[1].size(), 3)
  TEST_EQUAL(exp3[2].size(), 3)

  TEST_EQUAL(exp3[0].getFloatDataArrays().size(), 1)
  TEST_EQUAL(exp3[1].getFloatDataArrays().size(), 0)
  TEST_EQUAL(exp3[2].getFloatDataArrays().size(), 2)

  TEST_EQUAL(exp3[0].getFloatDataArrays()[0].getName(), "MDA1");
  TEST_EQUAL(exp3[2].getFloatDataArrays()[0].getName(), "MDA1");
  TEST_EQUAL(exp3[2].getFloatDataArrays()[1].getName(), "MDA2");

  TEST_EQUAL(exp3[0].getFloatDataArrays()[0].size(), 3)
  TEST_REAL_SIMILAR(exp3[0].getFloatDataArrays()[0][0], 1.3);
  TEST_REAL_SIMILAR(exp3[0].getFloatDataArrays()[0][1], 1.4);
  TEST_REAL_SIMILAR(exp3[0].getFloatDataArrays()[0][2], 1.5)

  TEST_EQUAL(exp3[2].getFloatDataArrays()[0].size(), 3)
  TEST_REAL_SIMILAR(exp3[2].getFloatDataArrays()[0][0], 1.3);
  TEST_REAL_SIMILAR(exp3[2].getFloatDataArrays()[0][1], 1.4);
  TEST_REAL_SIMILAR(exp3[2].getFloatDataArrays()[0][2], 1.5);

  TEST_EQUAL(exp3[2].getFloatDataArrays()[1].size(), 3)
  TEST_REAL_SIMILAR(exp3[2].getFloatDataArrays()[1][0], -2.3);
  TEST_REAL_SIMILAR(exp3[2].getFloatDataArrays()[1][1], -2.4);
  TEST_REAL_SIMILAR(exp3[2].getFloatDataArrays()[1][2], -2.5);

  //*********************************************
  //test if the storing meta data arrays without a name works

  exp3[0].getFloatDataArrays()[0].setName("");
  exp3[2].getFloatDataArrays()[0].setName("");
  exp3[2].getFloatDataArrays()[1].setName("");

  PeakMap exp4;
  file.store(filename, exp3);
  file.load(filename, exp4);

  TEST_EQUAL(exp.size(), 3)
  TEST_EQUAL(exp4[0].size(), 3)
  TEST_EQUAL(exp4[1].size(), 3)
  TEST_EQUAL(exp4[2].size(), 3)

  TEST_EQUAL(exp4[0].getFloatDataArrays().size(), 1)
  TEST_EQUAL(exp4[1].getFloatDataArrays().size(), 0)
  TEST_EQUAL(exp4[2].getFloatDataArrays().size(), 2)

  TEST_EQUAL(exp4[0].getFloatDataArrays()[0].getName(), "");
  TEST_EQUAL(exp4[2].getFloatDataArrays()[0].getName(), "");
  TEST_EQUAL(exp4[2].getFloatDataArrays()[1].getName(), "");

  TEST_EQUAL(exp4[0].getFloatDataArrays()[0].size(), 3)
  TEST_REAL_SIMILAR(exp4[0].getFloatDataArrays()[0][0], 1.3);
  TEST_REAL_SIMILAR(exp4[0].getFloatDataArrays()[0][1], 1.4);
  TEST_REAL_SIMILAR(exp4[0].getFloatDataArrays()[0][2], 1.5)

  TEST_EQUAL(exp4[2].getFloatDataArrays()[0].size(), 3)
  TEST_REAL_SIMILAR(exp4[2].getFloatDataArrays()[0][0], 1.3);
  TEST_REAL_SIMILAR(exp4[2].getFloatDataArrays()[0][1], 1.4);
  TEST_REAL_SIMILAR(exp4[2].getFloatDataArrays()[0][2], 1.5);

  TEST_EQUAL(exp4[2].getFloatDataArrays()[1].size(), 3)
  TEST_REAL_SIMILAR(exp4[2].getFloatDataArrays()[1][0], -2.3);
  TEST_REAL_SIMILAR(exp4[2].getFloatDataArrays()[1][1], -2.4);
  TEST_REAL_SIMILAR(exp4[2].getFloatDataArrays()[1][2], -2.5);
}
END_SECTION

START_SECTION([EXTRA] static bool isValid(const String& filename))
{
  std::string tmp_filename;
  MzDataFile f;
  PeakMap e;

  //test if empty file is valid
  NEW_TMP_FILE(tmp_filename);
  f.store(tmp_filename, e);
  TEST_EQUAL(f.isValid(tmp_filename, std::cerr), true);

  //test if filled file is valid
  NEW_TMP_FILE(tmp_filename);
  f.load(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"), e);
  f.store(tmp_filename, e);
  TEST_EQUAL(f.isValid(tmp_filename, std::cerr), true);
}
END_SECTION

START_SECTION(bool isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings))
{
  //This is not officially supported - the mapping file was hand-crafted by Marc Sturm
  NOT_TESTABLE
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
