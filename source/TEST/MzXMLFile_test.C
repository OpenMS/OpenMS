// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: MzXMLFile_test.C,v 1.16 2006/05/30 15:46:43 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Jens Joachim $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

///////////////////////////

START_TEST(DTAFile, "$Id: MzXMLFile_test.C,v 1.16 2006/05/30 15:46:43 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

MzXMLFile* ptr = 0;
CHECK(MzXMLFile())
	ptr = new MzXMLFile;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~MzXMLFile())
	delete ptr;
RESULT

CHECK(void load(const String& filename, MSExperiment<>& exp) throw (Exception::FileNotFound))
	PRECISION(0.01)

	MSExperiment< DRawDataPoint<1> > e;
	MzXMLFile mzxml;

	//test exception
	TEST_EXCEPTION( Exception::FileNotFound , mzxml.load("dummy/dummy.mzXML",e) )

	// real test
	mzxml.load("data/MzXMLFile_test_1.mzXML",e);
  //---------------------------------------------------------------------------
  // 60 : (120,100)
  // 120: (110,100) (120,200) (130,100)
  // 180: (100,100) (110,200) (120,300) (130,200) (140,100)
	//--------------------------------------------------------------------------- 
  TEST_EQUAL(e.size(), 4)
	TEST_REAL_EQUAL(e[0].getMSLevel(), 1)
	TEST_REAL_EQUAL(e[1].getMSLevel(), 1)
	TEST_REAL_EQUAL(e[2].getMSLevel(), 1)
	TEST_REAL_EQUAL(e[0].getRetentionTime(), 60)
	TEST_REAL_EQUAL(e[1].getRetentionTime(), 120)
	TEST_REAL_EQUAL(e[2].getRetentionTime(), 180)
	TEST_REAL_EQUAL(e[0].getContainer().size(), 1)
	TEST_REAL_EQUAL(e[1].getContainer().size(), 3)
	TEST_REAL_EQUAL(e[2].getContainer().size(), 5)

	TEST_REAL_EQUAL(e[0].getContainer()[0].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[0].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[0].getPosition()[0], 110)
	TEST_REAL_EQUAL(e[1].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[1].getContainer()[1].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[1].getContainer()[1].getIntensity(), 200)
	TEST_REAL_EQUAL(e[1].getContainer()[2].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[1].getContainer()[2].getIntensity(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[0].getPosition()[0], 100)
	TEST_REAL_EQUAL(e[2].getContainer()[0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[2].getContainer()[1].getPosition()[0], 110)
	TEST_REAL_EQUAL(e[2].getContainer()[1].getIntensity(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[2].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[2].getContainer()[2].getIntensity(), 300)
	TEST_REAL_EQUAL(e[2].getContainer()[3].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[2].getContainer()[3].getIntensity(), 200)
	TEST_REAL_EQUAL(e[2].getContainer()[4].getPosition()[0], 140)
	TEST_REAL_EQUAL(e[2].getContainer()[4].getIntensity(), 100)

	TEST_EQUAL(e[0].getMetaValue("URL1"), "www.open-ms.de")
	TEST_EQUAL(e[0].getMetaValue("URL2"), "www.uni-tuebingen.de")
	TEST_EQUAL(e[0].getComment(), "Scan Comment")
	TEST_EQUAL(e[0].getInstrumentSettings().getPolarity(), IonSource::POLNULL)
	TEST_EQUAL(e[1].getInstrumentSettings().getPolarity(), IonSource::POSITIVE)
	TEST_EQUAL(e[2].getInstrumentSettings().getPolarity(), IonSource::NEGATIVE)

	TEST_REAL_EQUAL(e[0].getInstrumentSettings().getMzRangeStart(), 0)
	TEST_REAL_EQUAL(e[1].getInstrumentSettings().getMzRangeStart(), 110)
	TEST_REAL_EQUAL(e[2].getInstrumentSettings().getMzRangeStart(), 100)
	TEST_REAL_EQUAL(e[0].getInstrumentSettings().getMzRangeStop(), 0)
	TEST_REAL_EQUAL(e[1].getInstrumentSettings().getMzRangeStop(), 130)
	TEST_REAL_EQUAL(e[2].getInstrumentSettings().getMzRangeStop(), 140)

	//---------------------------------------------------------------------------
  // const SourceFile& getSourceFile() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getSourceFile().getNameOfFile(), "File_test_1.raw");
  TEST_EQUAL(e.getSourceFile().getPathToFile(), "");
  TEST_EQUAL(e.getSourceFile().getFileType(), "RAWData");

	//---------------------------------------------------------------------------
  // const Software& getSoftware() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getSoftware().getName(), "MS-X");
  TEST_EQUAL(e.getSoftware().getVersion(), "1.0");
  TEST_EQUAL(e.getSoftware().getComment(), "conversion");
  TEST_EQUAL(e.getSoftware().getCompletionTime(), 4711);

	//---------------------------------------------------------------------------
  // const ProcessingMethod& getProcessingMethod() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getProcessingMethod().getDeisotoping(), true)
	TEST_EQUAL(e.getProcessingMethod().getChargeDeconvolution(), true)
	TEST_EQUAL(e.getProcessingMethod().getSpectrumType(), SpectrumSettings::PEAKS)
	TEST_EQUAL(e.getProcessingMethod().getMetaValue("processing 1#type 1"), "done 1")
	TEST_EQUAL(e.getProcessingMethod().getMetaValue("processing 2#type 2"), "done 2")
	TEST_EQUAL(e.getProcessingMethod().getMetaValue("#Comment"), "Software Comment")

	//---------------------------------------------------------------------------
  // const Instrument& getInstrument() const;
  //---------------------------------------------------------------------------
	const Instrument& inst = e.getInstrument();
	TEST_EQUAL(inst.getVendor(), "MS-Vendor")
	TEST_EQUAL(inst.getModel(), "MS 1")
	TEST_EQUAL(inst.getMetaValue("URL1"), "www.open-ms.de")
	TEST_EQUAL(inst.getMetaValue("URL2"), "www.uni-tuebingen.de")
	TEST_EQUAL(inst.getMetaValue("#Comment"), "Instrument Comment")
	TEST_EQUAL(inst.getIonSource().getIonizationMethod(), IonSource::ESI)
	TEST_EQUAL(inst.getIonDetector().getType(), IonDetector::FARADAYCUP)
	// UNSET:
  TEST_EQUAL(inst.getName(), "")
	TEST_EQUAL(inst.getCustomizations(), "")
	TEST_EQUAL(inst.getIonDetector().getResolution(), 0.0f)
	TEST_EQUAL(inst.getIonDetector().getADCSamplingFrequency(), 0.0f)
	TEST_EQUAL(inst.getIonSource().getInletType(), IonSource::INLETNULL)
	TEST_EQUAL(inst.getIonDetector().getAcquisitionMode(), IonDetector::ACQMODENULL)
	TEST_EQUAL(inst.getIonSource().getPolarity(), IonSource::POLNULL)

	TEST_EQUAL(inst.getMassAnalyzers().size(), 1)
  ABORT_IF(inst.getMassAnalyzers().size()!=1);
	TEST_EQUAL(inst.getMassAnalyzers()[0].getType(), MassAnalyzer::PAULIONTRAP)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionMethod(), MassAnalyzer::FWHM)
	// UNSET:
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionType(), MassAnalyzer::RESTYPENULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanFunction(), MassAnalyzer::SCANFCTNULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanDirection(), MassAnalyzer::SCANDIRNULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanLaw(), MassAnalyzer::SCANLAWNULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getTandemScanMethod(), MassAnalyzer::TANDEMNULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getReflectronState(), MassAnalyzer::REFLSTATENULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolution(), 0.0f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getAccuracy(), 0.0f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanRate(), 0.0f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanTime(), 0.0f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getTOFTotalPathLength(), 0.0f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getIsolationWidth(), 0.0f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getFinalMSExponent(), 0)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getMagneticFieldStrength(), 0.0f)

  //---------------------------------------------------------------------------
	// const Sample& getSample()
  //---------------------------------------------------------------------------
	// UNSET:
	TEST_EQUAL(e.getSample().getName(), "")
	TEST_EQUAL(e.getSample().getNumber(), "")
	TEST_EQUAL(e.getSample().getState(), Sample::SAMPLENULL)
 	TEST_EQUAL(e.getSample().getMass(), 0.0f)
	TEST_EQUAL(e.getSample().getVolume(), 0.0f)
	TEST_EQUAL(e.getSample().getConcentration(), 0.0f)

	// test reading different scheme and 64 bit data

	MSExperiment< DRawDataPoint<1> > e2;
	mzxml.load("data/MzXMLFile_test_3.mzXML",e2);

  TEST_EQUAL(e2.size(), 3)
	TEST_REAL_EQUAL(e2[0].getMSLevel(), 1)
	TEST_REAL_EQUAL(e2[1].getMSLevel(), 1)
	TEST_REAL_EQUAL(e2[2].getMSLevel(), 1)
	TEST_REAL_EQUAL(e2[0].getRetentionTime(), 60)
	TEST_REAL_EQUAL(e2[1].getRetentionTime(), 120)
	TEST_REAL_EQUAL(e2[2].getRetentionTime(), 180)
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


	const Instrument& inst2 = e2.getInstrument();
	TEST_EQUAL(inst2.getVendor(), "ABI")
	TEST_EQUAL(inst2.getModel(), "Neptune")
	TEST_EQUAL(inst2.getMetaValue("URL1"), "www.open-ms.de")
	TEST_EQUAL(inst2.getMetaValue("URL2"), "www.uni-tuebingen.de")
	TEST_EQUAL(inst2.getMetaValue("#Comment"), "Instrument Comment")
	TEST_EQUAL(inst2.getIonSource().getIonizationMethod(), IonSource::ESI)
	TEST_EQUAL(inst2.getIonDetector().getType(), IonDetector::TYPENULL)

	TEST_EQUAL(inst2.getMassAnalyzers().size(), 1)
  ABORT_IF(inst2.getMassAnalyzers().size()!=1);
	TEST_EQUAL(inst2.getMassAnalyzers()[0].getType(), MassAnalyzer::PAULIONTRAP)
	TEST_EQUAL(inst2.getMassAnalyzers()[0].getResolutionMethod(), MassAnalyzer::RESMETHNULL)
RESULT


CHECK(void store(const String& filename, const MSExperiment<>& exp) const throw (Exception::UnableToCreateFile))
	std::string tmp_filename;
  MSExperiment< DRawDataPoint<1> > e;
  MzXMLFile f;

  NEW_TMP_FILE(tmp_filename);
  f.load("data/MzXMLFile_test_1.mzXML",e);
	//test if empty spectra are skipped
	e.push_back(MSSpectrum< DRawDataPoint<1> >());
	e.push_back(MSSpectrum< DRawDataPoint<1> >());
	e.push_back(MSSpectrum< DRawDataPoint<1> >());
	TEST_EQUAL(e.size(), 7)


	f.store(tmp_filename,e);
	TextFile t1(tmp_filename, true);
	TextFile t2("data/MzXMLFile_test_1.mzXML", true);
	TEST_EQUAL(t1==t2,true)

	NEW_TMP_FILE(tmp_filename);
	f.load("data/MzXMLFile_test_2.mzXML",e);
	f.store(tmp_filename,e);
	TextFile t3(tmp_filename, true);
	TextFile t4("data/MzXMLFile_test_2.mzXML", true);
	TEST_EQUAL(t3==t4,true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
