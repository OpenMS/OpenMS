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

#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>

using namespace OpenMS;
using namespace std;

DRange<1> makeRange(float a, float b)
{
	DPosition<1> pa(a), pb(b);
	return DRange<1>(pa, pb);
}

///////////////////////////

START_TEST(MzXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzXMLFile* ptr = 0;
CHECK((MzXMLFile()))
	ptr = new MzXMLFile;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~MzXMLFile()))
	delete ptr;
RESULT

CHECK(const PeakFileOptions& getOptions() const)
	MzXMLFile file;
	TEST_EQUAL(file.getOptions().hasMSLevels(),false)
RESULT

CHECK(PeakFileOptions& getOptions())
	MzXMLFile file;
	file.getOptions().addMSLevel(1);
	TEST_EQUAL(file.getOptions().hasMSLevels(),true);
RESULT

CHECK((template<typename MapType> void load(const String& filename, MapType& map) ))
	PRECISION(0.01)
	
	MzXMLFile mzxml;

	//test exception
	{
		MSExperiment< Peak1D > e;
		TEST_EXCEPTION( Exception::FileNotFound , mzxml.load("dummy/dummy.mzXML",e) )
	}
	
	// first test
	{
		MSExperiment< Peak1D > e;
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
		TEST_REAL_EQUAL(e[0].size(), 1)
		TEST_REAL_EQUAL(e[1].size(), 3)
		TEST_REAL_EQUAL(e[2].size(), 5)
	
		TEST_REAL_EQUAL(e[0][0].getPosition()[0], 120)
		TEST_REAL_EQUAL(e[0][0].getIntensity(), 100)
		TEST_REAL_EQUAL(e[1][0].getPosition()[0], 110)
		TEST_REAL_EQUAL(e[1][0].getIntensity(), 100)
		TEST_REAL_EQUAL(e[1][1].getPosition()[0], 120)
		TEST_REAL_EQUAL(e[1][1].getIntensity(), 200)
		TEST_REAL_EQUAL(e[1][2].getPosition()[0], 130)
		TEST_REAL_EQUAL(e[1][2].getIntensity(), 100)
		TEST_REAL_EQUAL(e[2][0].getPosition()[0], 100)
		TEST_REAL_EQUAL(e[2][0].getIntensity(), 100)
		TEST_REAL_EQUAL(e[2][1].getPosition()[0], 110)
		TEST_REAL_EQUAL(e[2][1].getIntensity(), 200)
		TEST_REAL_EQUAL(e[2][2].getPosition()[0], 120)
		TEST_REAL_EQUAL(e[2][2].getIntensity(), 300)
		TEST_REAL_EQUAL(e[2][3].getPosition()[0], 130)
		TEST_REAL_EQUAL(e[2][3].getIntensity(), 200)
		TEST_REAL_EQUAL(e[2][4].getPosition()[0], 140)
		TEST_REAL_EQUAL(e[2][4].getIntensity(), 100)
	
		TEST_EQUAL(e[0].getMetaValue("URL1"), "www.open-ms.de")
		TEST_EQUAL(e[0].getMetaValue("URL2"), "www.uni-tuebingen.de")
		TEST_EQUAL(e[0].getComment(), "Scan Comment")
	
		//---------------------------------------------------------------------------
	  // const vector<SourceFile>& getSourceFiles() const;
	  //---------------------------------------------------------------------------
	  TEST_EQUAL(e.getSourceFiles().size(),2)
	  TEST_STRING_EQUAL(e.getSourceFiles()[0].getNameOfFile(), "File_test_1.raw");
	  TEST_STRING_EQUAL(e.getSourceFiles()[0].getPathToFile(), "");
	  TEST_STRING_EQUAL(e.getSourceFiles()[0].getFileType(), "RAWData");
	  TEST_STRING_EQUAL(e.getSourceFiles()[0].getChecksum(), "2fd4e1c67a2d28fced849ee1bb76e7391b93eb12");
		TEST_EQUAL(e.getSourceFiles()[0].getChecksumType(),SourceFile::SHA1)
	  TEST_STRING_EQUAL(e.getSourceFiles()[1].getNameOfFile(), "File_test_2.raw");
	  TEST_STRING_EQUAL(e.getSourceFiles()[1].getPathToFile(), "");
	  TEST_STRING_EQUAL(e.getSourceFiles()[1].getFileType(), "processedData");
	  TEST_STRING_EQUAL(e.getSourceFiles()[1].getChecksum(), "2fd4e1c67a2d28fced849ee1bb76e7391b93eb13");
		TEST_EQUAL(e.getSourceFiles()[1].getChecksumType(),SourceFile::SHA1)
	
		//---------------------------------------------------------------------------
	  // const vector<DataProcessing>& getDataProcessing() const;
	  //---------------------------------------------------------------------------
	  TEST_EQUAL(e.getDataProcessing().size(),2)

	  TEST_EQUAL(e.getDataProcessing()[0].getSoftware().getName(), "MS-X");
	  TEST_EQUAL(e.getDataProcessing()[0].getSoftware().getVersion(), "1.0");
	  TEST_STRING_EQUAL(e.getDataProcessing()[0].getMetaValue("#type").toString(), "conversion");
	  TEST_STRING_EQUAL(e.getDataProcessing()[0].getMetaValue("processing 1").toString(), "done 1");
	  TEST_STRING_EQUAL(e.getDataProcessing()[0].getMetaValue("processing 2").toString(), "done 2");
		String tmp;
		e.getDataProcessing()[0].getCompletionTime().get(tmp);
	  TEST_EQUAL(tmp, "2001-02-03 04:05:06");
		TEST_EQUAL(e.getDataProcessing()[0].getProcessingActions().size(),0)

	  TEST_EQUAL(e.getDataProcessing()[1].getSoftware().getName(), "MS-Y");
	  TEST_EQUAL(e.getDataProcessing()[1].getSoftware().getVersion(), "1.1");
	  TEST_STRING_EQUAL(e.getDataProcessing()[1].getMetaValue("#type").toString(), "processing");
	  TEST_REAL_EQUAL((DoubleReal)(e.getDataProcessing()[1].getMetaValue("#intensity_cutoff")), 3.4);
	  TEST_STRING_EQUAL(e.getDataProcessing()[1].getMetaValue("processing 3").toString(), "done 3");
		e.getDataProcessing()[1].getCompletionTime().get(tmp);
	  TEST_EQUAL(tmp, "0000-00-00 00:00:00");
		TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().size(),3)
		TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().count(DataProcessing::DEISOTOPING),1)
		TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().count(DataProcessing::CHARGE_DECONVOLUTION),1)
		TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().count(DataProcessing::PEAK_PICKING),1)

		//---------------------------------------------------------------------------
	  // const Instrument& getInstrument() const;
	  //---------------------------------------------------------------------------
		const Instrument& inst = e.getInstrument();
		TEST_EQUAL(inst.getVendor(), "MS-Vendor")
		TEST_EQUAL(inst.getModel(), "MS 1")
		TEST_EQUAL(inst.getMetaValue("URL1"), "www.open-ms.de")
		TEST_EQUAL(inst.getMetaValue("URL2"), "www.uni-tuebingen.de")
		TEST_EQUAL(inst.getMetaValue("#Comment"), "Instrument Comment")
	  TEST_EQUAL(inst.getName(), "")
		TEST_EQUAL(inst.getCustomizations(), "")
		TEST_EQUAL(inst.getIonSources().size(),1)
		TEST_EQUAL(inst.getIonSources()[0].getIonizationMethod(), IonSource::ESI)
		TEST_EQUAL(inst.getIonSources()[0].getInletType(), IonSource::INLETNULL)
		TEST_EQUAL(inst.getIonSources()[0].getPolarity(), IonSource::POLNULL)
		TEST_EQUAL(inst.getIonDetectors().size(),1)
		TEST_EQUAL(inst.getIonDetectors()[0].getType(), IonDetector::FARADAYCUP)
		TEST_REAL_EQUAL(inst.getIonDetectors()[0].getResolution(), 0.0f)
		TEST_REAL_EQUAL(inst.getIonDetectors()[0].getADCSamplingFrequency(), 0.0f)
		TEST_EQUAL(inst.getIonDetectors()[0].getAcquisitionMode(), IonDetector::ACQMODENULL)
		TEST_EQUAL(inst.getMassAnalyzers().size(), 1)
		TEST_EQUAL(inst.getMassAnalyzers()[0].getType(), MassAnalyzer::PAULIONTRAP)
		TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionMethod(), MassAnalyzer::FWHM)
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
		TEST_EQUAL(inst.getSoftware().getName(),"MS-Z")
		TEST_EQUAL(inst.getSoftware().getVersion(),"3.0")
		
	  //---------------------------------------------------------------------------
		// vector<ContactPerson>& getContacts()
	  //---------------------------------------------------------------------------
		const vector<ContactPerson>& contacts = e.getContacts();
		TEST_EQUAL(contacts.size(),1)
		TEST_STRING_EQUAL(contacts[0].getFirstName(),"FirstName")
		TEST_STRING_EQUAL(contacts[0].getLastName(),"LastName")
		TEST_STRING_EQUAL(contacts[0].getMetaValue("#phone"),"0049")
		TEST_STRING_EQUAL(contacts[0].getEmail(),"a@b.de")
		TEST_STRING_EQUAL(contacts[0].getURL(),"http://bla.de")
		TEST_STRING_EQUAL(contacts[0].getContactInfo(),"")
		
	  //---------------------------------------------------------------------------
		// const Sample& getSample()
	  //---------------------------------------------------------------------------
		TEST_EQUAL(e.getSample().getName(), "")
		TEST_EQUAL(e.getSample().getNumber(), "")
		TEST_EQUAL(e.getSample().getState(), Sample::SAMPLENULL)
	 	TEST_EQUAL(e.getSample().getMass(), 0.0f)
		TEST_EQUAL(e.getSample().getVolume(), 0.0f)
		TEST_EQUAL(e.getSample().getConcentration(), 0.0f)
	}
	
	//second test (to see if everything is initialized properly
	{
		MSExperiment< Peak1D > e;
		MzXMLFile().load("data/MzXMLFile_test_1.mzXML",e);
	  
	  //---------------------------------------------------------------------------
	  // 60 : (120,100)
	  // 120: (110,100) (120,200) (130,100)
	  // 180: (100,100) (110,200) (120,300) (130,200) (140,100)
		//--------------------------------------------------------------------------- 
	  TEST_EQUAL(e.size(), 4)
		TEST_REAL_EQUAL(e[0].getMSLevel(), 1)
		TEST_REAL_EQUAL(e[1].getMSLevel(), 1)
		TEST_REAL_EQUAL(e[2].getMSLevel(), 1)
		TEST_REAL_EQUAL(e[0].size(), 1)
		TEST_REAL_EQUAL(e[1].size(), 3)
		TEST_REAL_EQUAL(e[2].size(), 5)
	
		TEST_REAL_EQUAL(e[0][0].getPosition()[0], 120)
		TEST_REAL_EQUAL(e[0][0].getIntensity(), 100)
		TEST_REAL_EQUAL(e[1][0].getPosition()[0], 110)
		TEST_REAL_EQUAL(e[1][0].getIntensity(), 100)
		TEST_REAL_EQUAL(e[1][1].getPosition()[0], 120)
		TEST_REAL_EQUAL(e[1][1].getIntensity(), 200)
		TEST_REAL_EQUAL(e[1][2].getPosition()[0], 130)
		TEST_REAL_EQUAL(e[1][2].getIntensity(), 100)
		TEST_REAL_EQUAL(e[2][0].getPosition()[0], 100)
		TEST_REAL_EQUAL(e[2][0].getIntensity(), 100)
		TEST_REAL_EQUAL(e[2][1].getPosition()[0], 110)
		TEST_REAL_EQUAL(e[2][1].getIntensity(), 200)
		TEST_REAL_EQUAL(e[2][2].getPosition()[0], 120)
		TEST_REAL_EQUAL(e[2][2].getIntensity(), 300)
		TEST_REAL_EQUAL(e[2][3].getPosition()[0], 130)
		TEST_REAL_EQUAL(e[2][3].getIntensity(), 200)
		TEST_REAL_EQUAL(e[2][4].getPosition()[0], 140)
		TEST_REAL_EQUAL(e[2][4].getIntensity(), 100)
	
		TEST_EQUAL(e[0].getMetaValue("URL1"), "www.open-ms.de")
		TEST_EQUAL(e[0].getMetaValue("URL2"), "www.uni-tuebingen.de")
		TEST_EQUAL(e[0].getComment(), "Scan Comment")
	
		//---------------------------------------------------------------------------
	  // const vector<SourceFile>& getSourceFiles() const;
	  //---------------------------------------------------------------------------
	  TEST_EQUAL(e.getSourceFiles().size(),2)
	  TEST_STRING_EQUAL(e.getSourceFiles()[0].getNameOfFile(), "File_test_1.raw");
	  TEST_STRING_EQUAL(e.getSourceFiles()[0].getPathToFile(), "");
	  TEST_STRING_EQUAL(e.getSourceFiles()[0].getFileType(), "RAWData");
	  TEST_STRING_EQUAL(e.getSourceFiles()[0].getChecksum(), "2fd4e1c67a2d28fced849ee1bb76e7391b93eb12");
		TEST_EQUAL(e.getSourceFiles()[0].getChecksumType(),SourceFile::SHA1)
	  TEST_STRING_EQUAL(e.getSourceFiles()[1].getNameOfFile(), "File_test_2.raw");
	  TEST_STRING_EQUAL(e.getSourceFiles()[1].getPathToFile(), "");
	  TEST_STRING_EQUAL(e.getSourceFiles()[1].getFileType(), "processedData");
	  TEST_STRING_EQUAL(e.getSourceFiles()[1].getChecksum(), "2fd4e1c67a2d28fced849ee1bb76e7391b93eb13");
		TEST_EQUAL(e.getSourceFiles()[1].getChecksumType(),SourceFile::SHA1)
	
		//---------------------------------------------------------------------------
	  // const Software& getSoftware() const;
	  //---------------------------------------------------------------------------
	  TEST_EQUAL(e.getDataProcessing()[0].getSoftware().getName(), "MS-X");
	  TEST_EQUAL(e.getDataProcessing()[0].getSoftware().getVersion(), "1.0");
	  TEST_STRING_EQUAL(e.getDataProcessing()[0].getMetaValue("#type").toString(), "conversion");
	  TEST_STRING_EQUAL(e.getDataProcessing()[0].getMetaValue("processing 1").toString(), "done 1");
	  TEST_STRING_EQUAL(e.getDataProcessing()[0].getMetaValue("processing 2").toString(), "done 2");
		String tmp;
		e.getDataProcessing()[0].getCompletionTime().get(tmp);
	  TEST_EQUAL(tmp, "2001-02-03 04:05:06");
		TEST_EQUAL(e.getDataProcessing()[0].getProcessingActions().size(),0)

	  TEST_EQUAL(e.getDataProcessing()[1].getSoftware().getName(), "MS-Y");
	  TEST_EQUAL(e.getDataProcessing()[1].getSoftware().getVersion(), "1.1");
	  TEST_STRING_EQUAL(e.getDataProcessing()[1].getMetaValue("#type").toString(), "processing");
	  TEST_REAL_EQUAL((DoubleReal)(e.getDataProcessing()[1].getMetaValue("#intensity_cutoff")), 3.4);
	  TEST_STRING_EQUAL(e.getDataProcessing()[1].getMetaValue("processing 3").toString(), "done 3");
		e.getDataProcessing()[1].getCompletionTime().get(tmp);
	  TEST_EQUAL(tmp, "0000-00-00 00:00:00");
		TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().size(),3)
		TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().count(DataProcessing::DEISOTOPING),1)
		TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().count(DataProcessing::CHARGE_DECONVOLUTION),1)
		TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().count(DataProcessing::PEAK_PICKING),1)
	
		//---------------------------------------------------------------------------
	  // const Instrument& getInstrument() const;
	  //---------------------------------------------------------------------------
		const Instrument& inst = e.getInstrument();
		TEST_EQUAL(inst.getVendor(), "MS-Vendor")
		TEST_EQUAL(inst.getModel(), "MS 1")
		TEST_EQUAL(inst.getMetaValue("URL1"), "www.open-ms.de")
		TEST_EQUAL(inst.getMetaValue("URL2"), "www.uni-tuebingen.de")
		TEST_EQUAL(inst.getMetaValue("#Comment"), "Instrument Comment")
	  TEST_EQUAL(inst.getName(), "")
		TEST_EQUAL(inst.getCustomizations(), "")
		TEST_EQUAL(inst.getIonSources().size(),1)
		TEST_EQUAL(inst.getIonSources()[0].getIonizationMethod(), IonSource::ESI)
		TEST_EQUAL(inst.getIonSources()[0].getInletType(), IonSource::INLETNULL)
		TEST_EQUAL(inst.getIonSources()[0].getPolarity(), IonSource::POLNULL)
		TEST_EQUAL(inst.getIonDetectors().size(),1)
		TEST_EQUAL(inst.getIonDetectors()[0].getType(), IonDetector::FARADAYCUP)
		TEST_REAL_EQUAL(inst.getIonDetectors()[0].getResolution(), 0.0f)
		TEST_REAL_EQUAL(inst.getIonDetectors()[0].getADCSamplingFrequency(), 0.0f)
		TEST_EQUAL(inst.getIonDetectors()[0].getAcquisitionMode(), IonDetector::ACQMODENULL)
		TEST_EQUAL(inst.getMassAnalyzers().size(), 1)
		TEST_EQUAL(inst.getMassAnalyzers()[0].getType(), MassAnalyzer::PAULIONTRAP)
		TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionMethod(), MassAnalyzer::FWHM)
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
		TEST_EQUAL(inst.getSoftware().getName(),"MS-Z")
		TEST_EQUAL(inst.getSoftware().getVersion(),"3.0")

	  //---------------------------------------------------------------------------
		// vector<ContactPerson>& getContacts()
	  //---------------------------------------------------------------------------
		const vector<ContactPerson>& contacts = e.getContacts();
		TEST_EQUAL(contacts.size(),1)
		TEST_STRING_EQUAL(contacts[0].getFirstName(),"FirstName")
		TEST_STRING_EQUAL(contacts[0].getLastName(),"LastName")
		TEST_STRING_EQUAL(contacts[0].getMetaValue("#phone"),"0049")
		TEST_STRING_EQUAL(contacts[0].getEmail(),"a@b.de")
		TEST_STRING_EQUAL(contacts[0].getURL(),"http://bla.de")
		TEST_STRING_EQUAL(contacts[0].getContactInfo(),"")
		
	  //---------------------------------------------------------------------------
		// const Sample& getSample()
	  //---------------------------------------------------------------------------
		TEST_EQUAL(e.getSample().getName(), "")
		TEST_EQUAL(e.getSample().getNumber(), "")
		TEST_EQUAL(e.getSample().getState(), Sample::SAMPLENULL)
	 	TEST_EQUAL(e.getSample().getMass(), 0.0f)
		TEST_EQUAL(e.getSample().getVolume(), 0.0f)
		TEST_EQUAL(e.getSample().getConcentration(), 0.0f)
	}

	// test reading 64 bit data

	MSExperiment< Peak1D > e2;
	mzxml.load("data/MzXMLFile_test_3.mzXML",e2);

  TEST_EQUAL(e2.size(), 3)
	TEST_REAL_EQUAL(e2[0].getMSLevel(), 1)
	TEST_REAL_EQUAL(e2[1].getMSLevel(), 1)
	TEST_REAL_EQUAL(e2[2].getMSLevel(), 1)
	TEST_REAL_EQUAL(e2[0].getRT(), 1)
	TEST_REAL_EQUAL(e2[1].getRT(), 121)
	TEST_REAL_EQUAL(e2[2].getRT(), 3661)
	TEST_REAL_EQUAL(e2[0].size(), 1)
	TEST_REAL_EQUAL(e2[1].size(), 3)
	TEST_REAL_EQUAL(e2[2].size(), 5)

	TEST_REAL_EQUAL(e2[0][0].getPosition()[0], 120)
	TEST_REAL_EQUAL(e2[0][0].getIntensity(), 100)
	TEST_REAL_EQUAL(e2[1][0].getPosition()[0], 110)
	TEST_REAL_EQUAL(e2[1][0].getIntensity(), 100)
	TEST_REAL_EQUAL(e2[1][1].getPosition()[0], 120)
	TEST_REAL_EQUAL(e2[1][1].getIntensity(), 200)
	TEST_REAL_EQUAL(e2[1][2].getPosition()[0], 130)
	TEST_REAL_EQUAL(e2[1][2].getIntensity(), 100)
	TEST_REAL_EQUAL(e2[2][0].getPosition()[0], 100)
	TEST_REAL_EQUAL(e2[2][0].getIntensity(), 100)
	TEST_REAL_EQUAL(e2[2][1].getPosition()[0], 110)
	TEST_REAL_EQUAL(e2[2][1].getIntensity(), 200)
	TEST_REAL_EQUAL(e2[2][2].getPosition()[0], 120)
	TEST_REAL_EQUAL(e2[2][2].getIntensity(), 300)
	TEST_REAL_EQUAL(e2[2][3].getPosition()[0], 130)
	TEST_REAL_EQUAL(e2[2][3].getIntensity(), 200)
	TEST_REAL_EQUAL(e2[2][4].getPosition()[0], 140)
	TEST_REAL_EQUAL(e2[2][4].getIntensity(), 100)

RESULT

CHECK((template<typename MapType> void store(const String& filename, const MapType& map) const ))
	std::string tmp_filename;
  MSExperiment< Peak1D > e1, e2;
  MzXMLFile f;

  NEW_TMP_FILE(tmp_filename);
 	 f.load("data/MzXMLFile_test_1.mzXML",e1);
	TEST_EQUAL(e1.size(), 4)

	f.store(tmp_filename,e1);
	f.load(tmp_filename,e2);
	TEST_EQUAL(e1==e2, true);

	NEW_TMP_FILE(tmp_filename);
	f.load("data/MzXMLFile_test_2.mzXML",e1);
	f.store(tmp_filename,e1);
	f.load(tmp_filename,e2);
	TEST_EQUAL(e1==e2,true);
RESULT

CHECK(([EXTRA] load/store for Float Kernel Traits))
	std::string tmp_filename;
	NEW_TMP_FILE(tmp_filename);
	
  MSExperiment< Peak1D > e1, e2;
  MzXMLFile f;
	
	f.load("data/MzXMLFile_test_2.mzXML",e1);
	f.store(tmp_filename,e1);
	f.load(tmp_filename,e2);
	TEST_EQUAL(e1==e2,true);
RESULT

// check for Float Kernel traits
CHECK(([EXTRA] load/store for empty scans))
	std::string tmp_filename;
	NEW_TMP_FILE(tmp_filename);

  MzXMLFile f;
	MSExperiment<> e2;
	e2.resize(5);
	
	f.store(tmp_filename,e2);
	f.load(tmp_filename,e2);
	
	e2.updateRanges();
	
	TEST_EQUAL(e2.size(),5);
	TEST_EQUAL(e2.getSize(),0);
RESULT

CHECK(([EXTRA] load with optional attributes))
	PRECISION(0.01)

	MSExperiment< Peak1D > e;
	MzXMLFile mzxml;
	
	// real test
	mzxml.load("data/MzXMLFile_test_4.mzXML",e);
  //---------------------------------------------------------------------------
  // 60 : (120,100)
  // 120: (110,100) (120,200) (130,100)
  // 180: (100,100) (110,200) (120,300) (130,200) (140,100)
	//--------------------------------------------------------------------------- 
  TEST_EQUAL(e.size(), 4)
	TEST_REAL_EQUAL(e[0].getRT(), 60)
	TEST_REAL_EQUAL(e[1].getRT(), 120)
	TEST_REAL_EQUAL(e[2].getRT(), 180)

	TEST_EQUAL(e[0].getInstrumentSettings().getPolarity(), IonSource::POLNULL)
	TEST_EQUAL(e[1].getInstrumentSettings().getPolarity(), IonSource::POSITIVE)
	TEST_EQUAL(e[2].getInstrumentSettings().getPolarity(), IonSource::NEGATIVE)

	TEST_EQUAL(e[0].getInstrumentSettings().getScanWindows().size(), 0)
	TEST_EQUAL(e[1].getInstrumentSettings().getScanWindows().size(), 1)
	TEST_REAL_EQUAL(e[1].getInstrumentSettings().getScanWindows()[0].begin, 110)
	TEST_REAL_EQUAL(e[1].getInstrumentSettings().getScanWindows()[0].end, 0)
	TEST_EQUAL(e[2].getInstrumentSettings().getScanWindows().size(), 1)
	TEST_REAL_EQUAL(e[2].getInstrumentSettings().getScanWindows()[0].begin, 100)
	TEST_REAL_EQUAL(e[2].getInstrumentSettings().getScanWindows()[0].end, 140)
	TEST_EQUAL(e[3].getInstrumentSettings().getScanWindows().size(), 0)

	//---------------------------------------------------------------------------
  // const vector<DataProcessing>& getDataProcessing() const;
  //---------------------------------------------------------------------------
	TEST_EQUAL(e.getDataProcessing().size(),1)
	String tmp;
	e.getDataProcessing().back().getCompletionTime().get(tmp);
  TEST_EQUAL(tmp, "2001-02-03 04:05:06");

RESULT

CHECK(([EXTRA] load with metadata only flag))
	PRECISION(0.01)

	MSExperiment< Peak1D > e;
	MzXMLFile mzxml;
	mzxml.getOptions().setMetadataOnly(true);

	// real test
	mzxml.load("data/MzXMLFile_test_1.mzXML",e);
	
	//check number of scans
	TEST_EQUAL(e.size(),0)
	
	//---------------------------------------------------------------------------
  // const vector<SourceFile>& getSourceFiles() const;
  //---------------------------------------------------------------------------
	  TEST_EQUAL(e.getSourceFiles().size(),2)
	  TEST_STRING_EQUAL(e.getSourceFiles()[0].getNameOfFile(), "File_test_1.raw");
	  TEST_STRING_EQUAL(e.getSourceFiles()[0].getPathToFile(), "");
	  TEST_STRING_EQUAL(e.getSourceFiles()[0].getFileType(), "RAWData");
	  TEST_STRING_EQUAL(e.getSourceFiles()[0].getChecksum(), "2fd4e1c67a2d28fced849ee1bb76e7391b93eb12");
		TEST_EQUAL(e.getSourceFiles()[0].getChecksumType(),SourceFile::SHA1)
	  TEST_STRING_EQUAL(e.getSourceFiles()[1].getNameOfFile(), "File_test_2.raw");
	  TEST_STRING_EQUAL(e.getSourceFiles()[1].getPathToFile(), "");
	  TEST_STRING_EQUAL(e.getSourceFiles()[1].getFileType(), "processedData");
	  TEST_STRING_EQUAL(e.getSourceFiles()[1].getChecksum(), "2fd4e1c67a2d28fced849ee1bb76e7391b93eb13");
		TEST_EQUAL(e.getSourceFiles()[1].getChecksumType(),SourceFile::SHA1)

	//---------------------------------------------------------------------------
  // const Software& getSoftware() const;
  //---------------------------------------------------------------------------
	TEST_EQUAL(e.getDataProcessing().size(),2)
  TEST_EQUAL(e.getDataProcessing().back().getSoftware().getName(), "MS-Y");
  TEST_EQUAL(e.getDataProcessing().back().getSoftware().getVersion(), "1.1");

	//---------------------------------------------------------------------------
  // const Instrument& getInstrument() const;
  //---------------------------------------------------------------------------
	const Instrument& inst = e.getInstrument();
	TEST_EQUAL(inst.getVendor(), "MS-Vendor")
	TEST_EQUAL(inst.getModel(), "MS 1")
  TEST_EQUAL(inst.getName(), "")
	TEST_EQUAL(inst.getCustomizations(), "")
	TEST_EQUAL(inst.getMetaValue("URL1"), "www.open-ms.de")
	TEST_EQUAL(inst.getMetaValue("URL2"), "www.uni-tuebingen.de")
	TEST_EQUAL(inst.getMetaValue("#Comment"), "Instrument Comment")
	TEST_EQUAL(inst.getIonSources().size(),1)
	TEST_EQUAL(inst.getIonSources()[0].getIonizationMethod(), IonSource::ESI)
	TEST_EQUAL(inst.getIonSources()[0].getInletType(), IonSource::INLETNULL)
	TEST_EQUAL(inst.getIonSources()[0].getPolarity(), IonSource::POLNULL)
	TEST_EQUAL(inst.getIonDetectors().size(),1)
	TEST_EQUAL(inst.getIonDetectors()[0].getType(), IonDetector::FARADAYCUP)
	TEST_REAL_EQUAL(inst.getIonDetectors()[0].getResolution(), 0.0f)
	TEST_REAL_EQUAL(inst.getIonDetectors()[0].getADCSamplingFrequency(), 0.0f)
	TEST_EQUAL(inst.getIonDetectors()[0].getAcquisitionMode(), IonDetector::ACQMODENULL)
	TEST_EQUAL(inst.getMassAnalyzers().size(), 1)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getType(), MassAnalyzer::PAULIONTRAP)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionMethod(), MassAnalyzer::FWHM)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionType(), MassAnalyzer::RESTYPENULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanFunction(), MassAnalyzer::SCANFCTNULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanDirection(), MassAnalyzer::SCANDIRNULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanLaw(), MassAnalyzer::SCANLAWNULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getTandemScanMethod(), MassAnalyzer::TANDEMNULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getReflectronState(), MassAnalyzer::REFLSTATENULL)
	TEST_REAL_EQUAL(inst.getMassAnalyzers()[0].getResolution(), 0.0f)
	TEST_REAL_EQUAL(inst.getMassAnalyzers()[0].getAccuracy(), 0.0f)
	TEST_REAL_EQUAL(inst.getMassAnalyzers()[0].getScanRate(), 0.0f)
	TEST_REAL_EQUAL(inst.getMassAnalyzers()[0].getScanTime(), 0.0f)
	TEST_REAL_EQUAL(inst.getMassAnalyzers()[0].getTOFTotalPathLength(), 0.0f)
	TEST_REAL_EQUAL(inst.getMassAnalyzers()[0].getIsolationWidth(), 0.0f)
	TEST_REAL_EQUAL(inst.getMassAnalyzers()[0].getFinalMSExponent(), 0)
	TEST_REAL_EQUAL(inst.getMassAnalyzers()[0].getMagneticFieldStrength(), 0.0f)
	TEST_EQUAL(inst.getSoftware().getName(),"MS-Z")
	TEST_EQUAL(inst.getSoftware().getVersion(),"3.0")

  //---------------------------------------------------------------------------
	// vector<ContactPerson>& getContacts()
  //---------------------------------------------------------------------------
	const vector<ContactPerson>& contacts = e.getContacts();
	TEST_EQUAL(contacts.size(),1)
	TEST_STRING_EQUAL(contacts[0].getFirstName(),"FirstName")
	TEST_STRING_EQUAL(contacts[0].getLastName(),"LastName")
	TEST_STRING_EQUAL(contacts[0].getMetaValue("#phone"),"0049")
	TEST_STRING_EQUAL(contacts[0].getEmail(),"a@b.de")
	TEST_STRING_EQUAL(contacts[0].getURL(),"http://bla.de")
		TEST_STRING_EQUAL(contacts[0].getContactInfo(),"")
		

  //---------------------------------------------------------------------------
	// const Sample& getSample()
  //---------------------------------------------------------------------------
	TEST_EQUAL(e.getSample().getName(), "")
	TEST_EQUAL(e.getSample().getNumber(), "")
	TEST_EQUAL(e.getSample().getState(), Sample::SAMPLENULL)
 	TEST_EQUAL(e.getSample().getMass(), 0.0f)
	TEST_EQUAL(e.getSample().getVolume(), 0.0f)
	TEST_EQUAL(e.getSample().getConcentration(), 0.0f)
RESULT

CHECK(([EXTRA] load with selected MS levels))
	PRECISION(0.01)

	MSExperiment< Peak1D > e;
	MzXMLFile mzxml;
	
	// load only MS level 1
	mzxml.getOptions().addMSLevel(1);
	mzxml.load("data/MzXMLFile_test_6.mzXML",e);

	TEST_EQUAL(e.size(), 3)
	TEST_REAL_EQUAL(e[0].getMSLevel(), 1)
	TEST_EQUAL((UInt)(e[0].getMetaValue("original_spectrum_number")),0)
	TEST_REAL_EQUAL(e[1].getMSLevel(), 1)
	TEST_EQUAL((UInt)(e[1].getMetaValue("original_spectrum_number")),2)
	TEST_REAL_EQUAL(e[2].getMSLevel(), 1)
	TEST_EQUAL((UInt)(e[2].getMetaValue("original_spectrum_number")),4)
	TEST_REAL_EQUAL(e[0].size(), 1)
	TEST_REAL_EQUAL(e[1].size(), 3)
	TEST_REAL_EQUAL(e[2].size(), 5)
	
	// load all levels
	mzxml.getOptions().clearMSLevels();
	mzxml.load("data/MzXMLFile_test_1.mzXML",e);

	TEST_EQUAL(e.size(), 4)
RESULT

CHECK(([EXTRA] load with selected MZ range))
	PRECISION(0.01)

	MSExperiment< Peak1D > e;
	MzXMLFile mzxml;

	mzxml.getOptions().setMZRange(makeRange(115,135));
	mzxml.load("data/MzXMLFile_test_1.mzXML",e);
	//---------------------------------------------------------------------------
	// 60 : +(120,100)
	// 120: -(110,100) +(120,200) +(130,100)
	// 180: -(100,100) -(110,200) +(120,300) +(130,200) -(140,100)
	//--------------------------------------------------------------------------- 
	
	TEST_REAL_EQUAL(e[0].size(), 1)
	TEST_REAL_EQUAL(e[1].size(), 2)
	TEST_REAL_EQUAL(e[2].size(), 2)

	TEST_REAL_EQUAL(e[0][0].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[0][0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[1][0].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[1][0].getIntensity(), 200)
	TEST_REAL_EQUAL(e[1][1].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[1][1].getIntensity(), 100)
	TEST_REAL_EQUAL(e[2][0].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[2][0].getIntensity(), 300)
	TEST_REAL_EQUAL(e[2][1].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[2][1].getIntensity(), 200)
RESULT

CHECK(([EXTRA] load with RT range))
	PRECISION(0.01)

	MSExperiment< Peak1D > e;
	MzXMLFile mzxml;
	mzxml.getOptions().setRTRange(makeRange(100, 200));
	mzxml.load("data/MzXMLFile_test_2.mzXML",e);
	//---------------------------------------------------------------------------
	// 120: (110,100) (120,200) (130,100)
	// 180: (100,100) (110,200) (120,300) (130,200) (140,100)
	//--------------------------------------------------------------------------- 
 	TEST_REAL_EQUAL(e.size(), 2)
 	TEST_REAL_EQUAL(e[0].size(), 3)
 	TEST_REAL_EQUAL(e[1].size(), 5)

	TEST_REAL_EQUAL(e[0][0].getPosition()[0], 110)
	TEST_REAL_EQUAL(e[0][0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[0][1].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[0][1].getIntensity(), 200)
	TEST_REAL_EQUAL(e[0][2].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[0][2].getIntensity(), 100)
	TEST_REAL_EQUAL(e[1][0].getPosition()[0], 100)
	TEST_REAL_EQUAL(e[1][0].getIntensity(), 100)
	TEST_REAL_EQUAL(e[1][1].getPosition()[0], 110)
	TEST_REAL_EQUAL(e[1][1].getIntensity(), 200)
	TEST_REAL_EQUAL(e[1][2].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[1][2].getIntensity(), 300)
	TEST_REAL_EQUAL(e[1][3].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[1][3].getIntensity(), 200)
	TEST_REAL_EQUAL(e[1][4].getPosition()[0], 140)
	TEST_REAL_EQUAL(e[1][4].getIntensity(), 100)
RESULT

CHECK(([EXTRA] load with intensity range))
	PRECISION(0.01)

	MSExperiment< Peak1D > e;
	MzXMLFile mzxml;
	mzxml.getOptions().setIntensityRange(makeRange(150, 350));
	mzxml.load("data/MzXMLFile_test_1.mzXML",e);
	//---------------------------------------------------------------------------
	// 60 : -(120,100)
	// 120: -(110,100) +(120,200) -(130,100)
	// 180: -(100,100) +(110,200) +(120,300) +(130,200) -(140,100)
	//--------------------------------------------------------------------------- 
	TEST_REAL_EQUAL(e[0].size(), 0)
	TEST_REAL_EQUAL(e[1].size(), 1)
	TEST_REAL_EQUAL(e[2].size(), 3)

	TEST_REAL_EQUAL(e[1][0].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[1][0].getIntensity(), 200)
	TEST_REAL_EQUAL(e[2][0].getPosition()[0], 110)
	TEST_REAL_EQUAL(e[2][0].getIntensity(), 200)
	TEST_REAL_EQUAL(e[2][1].getPosition()[0], 120)
	TEST_REAL_EQUAL(e[2][1].getIntensity(), 300)
	TEST_REAL_EQUAL(e[2][2].getPosition()[0], 130)
	TEST_REAL_EQUAL(e[2][2].getIntensity(), 200)
RESULT

CHECK(([EXTRA] load/store for nested scans))
	std::string tmp_filename;
	NEW_TMP_FILE(tmp_filename);
  MzXMLFile f;
	MSExperiment<> e2;
	e2.resize(5);
	
	//alternating
	e2[0].setMSLevel(1);
	e2[1].setMSLevel(2);
	e2[2].setMSLevel(1);
	e2[3].setMSLevel(2);
	e2[4].setMSLevel(1);
	f.store(tmp_filename,e2);
	f.load(tmp_filename,e2);
	TEST_EQUAL(e2.size(),5);

	//ending with ms level 2
	e2[0].setMSLevel(1);
	e2[1].setMSLevel(2);
	e2[2].setMSLevel(1);
	e2[3].setMSLevel(2);
	e2[4].setMSLevel(2);
	f.store(tmp_filename,e2);
	f.load(tmp_filename,e2);
	TEST_EQUAL(e2.size(),5);

	//MS level 1-3
	e2[0].setMSLevel(1);
	e2[1].setMSLevel(2);
	e2[2].setMSLevel(3);
	e2[3].setMSLevel(2);
	e2[4].setMSLevel(3);
	f.store(tmp_filename,e2);
	f.load(tmp_filename,e2);
	TEST_EQUAL(e2.size(),5);

	//MS level 2
	e2[0].setMSLevel(2);
	e2[1].setMSLevel(2);
	e2[2].setMSLevel(2);
	e2[3].setMSLevel(2);
	e2[4].setMSLevel(2);
	f.store(tmp_filename,e2);
	f.load(tmp_filename,e2);
	TEST_EQUAL(e2.size(),5);

	//MS level 2-3
	e2[0].setMSLevel(2);
	e2[1].setMSLevel(2);
	e2[2].setMSLevel(3);
	e2[3].setMSLevel(2);
	e2[4].setMSLevel(3);
	f.store(tmp_filename,e2);
	f.load(tmp_filename,e2);
	TEST_EQUAL(e2.size(),5);

	//MS level 1-3 (not starting with 1)
	e2[0].setMSLevel(2);
	e2[1].setMSLevel(1);
	e2[2].setMSLevel(2);
	e2[3].setMSLevel(3);
	e2[4].setMSLevel(1);
	f.store(tmp_filename,e2);
	f.load(tmp_filename,e2);
	TEST_EQUAL(e2.size(),5);
RESULT

CHECK([EXTRA] static bool isValid(const String& filename))
	std::string tmp_filename;
  MzXMLFile f;
  MSExperiment<> e;
  
  //Note: empty mzXML files are not valid, thus this test is omitted

	//test if fill file is valid
	NEW_TMP_FILE(tmp_filename);
	f.load("data/MzXMLFile_test_1.mzXML",e);
  f.store(tmp_filename,e);
  TEST_EQUAL(f.isValid(tmp_filename),true);
RESULT

CHECK([EXTRA] loading a minimal file containing one spectrum  - with whitespaces inside the base64 data)
	MSExperiment<> e;
  MzXMLFile f;
  f.load("data/MzXMLFile_test_minimal.mzXML",e);
  TEST_EQUAL(e.size(),1)
  TEST_EQUAL(e[0].size(),1)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
