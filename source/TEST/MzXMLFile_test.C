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

#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>

using namespace OpenMS;
using namespace std;

DRange<1> makeRange(DoubleReal a, DoubleReal b)
{
	DPosition<1> pa(a), pb(b);
	return DRange<1>(pa, pb);
}

///////////////////////////

START_TEST(MzXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzXMLFile* ptr = 0;
START_SECTION((MzXMLFile()))
	ptr = new MzXMLFile;
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~MzXMLFile()))
	delete ptr;
END_SECTION

START_SECTION(const PeakFileOptions& getOptions() const)
	MzXMLFile file;
	TEST_EQUAL(file.getOptions().hasMSLevels(),false)
END_SECTION

START_SECTION(PeakFileOptions& getOptions())
	MzXMLFile file;
	file.getOptions().addMSLevel(1);
	TEST_EQUAL(file.getOptions().hasMSLevels(),true);
END_SECTION

START_SECTION((template<typename MapType> void load(const String& filename, MapType& map) ))
	TOLERANCE_ABSOLUTE(0.01)

	MzXMLFile file;

	//exception
	MSExperiment<> e;
	TEST_EXCEPTION( Exception::FileNotFound , file.load("dummy/dummy.mzXML",e) )

	//real test
	file.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"),e);

	//test DocumentIdentifier addition
	TEST_STRING_EQUAL(e.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"));
	TEST_STRING_EQUAL(FileHandler::typeToName(e.getLoadedFileType()),"mzXML");

  //---------------------------------------------------------------------------
	// actual peak data
  // 60 : (120,100)
  // 120: (110,100) (120,200) (130,100)
  // 180: (100,100) (110,200) (120,300) (130,200) (140,100)
	//---------------------------------------------------------------------------
  TEST_EQUAL(e.size(), 4)
	TEST_EQUAL(e[0].getMSLevel(), 1)
	TEST_EQUAL(e[1].getMSLevel(), 1)
	TEST_EQUAL(e[2].getMSLevel(), 1)
	TEST_EQUAL(e[3].getMSLevel(), 2)
	TEST_EQUAL(e[0].size(), 1)
	TEST_EQUAL(e[1].size(), 3)
	TEST_EQUAL(e[2].size(), 5)
	TEST_EQUAL(e[3].size(), 5)
	TEST_STRING_EQUAL(e[0].getNativeID(),"scan=10")
	TEST_STRING_EQUAL(e[1].getNativeID(),"scan=11")
	TEST_STRING_EQUAL(e[2].getNativeID(),"scan=12")
	TEST_STRING_EQUAL(e[3].getNativeID(),"scan=13")

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

	TEST_EQUAL(e[0].getMetaValue("URL1"), "www.open-ms.de")
	TEST_EQUAL(e[0].getMetaValue("URL2"), "www.uni-tuebingen.de")
	TEST_EQUAL(e[0].getComment(), "Scan Comment")

	//---------------------------------------------------------------------------
  // source file
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
  // data processing (assigned to each spectrum)
  //---------------------------------------------------------------------------
	for (Size i=0; i<e.size(); ++i)
	{
		TEST_EQUAL(e[i].getDataProcessing().size(),2)

		TEST_EQUAL(e[i].getDataProcessing()[0].getSoftware().getName(), "MS-X");
		TEST_EQUAL(e[i].getDataProcessing()[0].getSoftware().getVersion(), "1.0");
		TEST_STRING_EQUAL(e[i].getDataProcessing()[0].getMetaValue("#type").toString(), "conversion");
		TEST_STRING_EQUAL(e[i].getDataProcessing()[0].getMetaValue("processing 1").toString(), "done 1");
		TEST_STRING_EQUAL(e[i].getDataProcessing()[0].getMetaValue("processing 2").toString(), "done 2");
		TEST_EQUAL(e[i].getDataProcessing()[0].getCompletionTime().get(), "2001-02-03 04:05:06");
		TEST_EQUAL(e[i].getDataProcessing()[0].getProcessingActions().size(),0)

		TEST_EQUAL(e[i].getDataProcessing()[1].getSoftware().getName(), "MS-Y");
		TEST_EQUAL(e[i].getDataProcessing()[1].getSoftware().getVersion(), "1.1");
		TEST_STRING_EQUAL(e[i].getDataProcessing()[1].getMetaValue("#type").toString(), "processing");
		TEST_REAL_SIMILAR((DoubleReal)(e[i].getDataProcessing()[1].getMetaValue("#intensity_cutoff")), 3.4);
		TEST_STRING_EQUAL(e[i].getDataProcessing()[1].getMetaValue("processing 3").toString(), "done 3");
		TEST_EQUAL(e[i].getDataProcessing()[1].getCompletionTime().get(), "0000-00-00 00:00:00");
		TEST_EQUAL(e[i].getDataProcessing()[1].getProcessingActions().size(),3)
		TEST_EQUAL(e[i].getDataProcessing()[1].getProcessingActions().count(DataProcessing::DEISOTOPING),1)
		TEST_EQUAL(e[i].getDataProcessing()[1].getProcessingActions().count(DataProcessing::CHARGE_DECONVOLUTION),1)
		TEST_EQUAL(e[i].getDataProcessing()[1].getProcessingActions().count(DataProcessing::PEAK_PICKING),1)
	}
	
	//---------------------------------------------------------------------------
  // instrument
  //---------------------------------------------------------------------------
	const Instrument& inst = e.getInstrument();
	TEST_EQUAL(inst.getVendor(), "MS-Vendor")
	TEST_EQUAL(inst.getModel(), "MS 1")
	TEST_EQUAL(inst.getMetaValue("URL1"), "www.open-ms.de")
	TEST_EQUAL(inst.getMetaValue("URL2"), "www.uni-tuebingen.de")
	TEST_EQUAL(inst.getMetaValue("#comment"), "Instrument Comment")
  TEST_EQUAL(inst.getName(), "")
	TEST_EQUAL(inst.getCustomizations(), "")
	TEST_EQUAL(inst.getIonSources().size(),1)
	TEST_EQUAL(inst.getIonSources()[0].getIonizationMethod(), IonSource::ESI)
	TEST_EQUAL(inst.getIonSources()[0].getInletType(), IonSource::INLETNULL)
	TEST_EQUAL(inst.getIonSources()[0].getPolarity(), IonSource::POLNULL)
	TEST_EQUAL(inst.getIonDetectors().size(),1)
	TEST_EQUAL(inst.getIonDetectors()[0].getType(), IonDetector::FARADAYCUP)
	TEST_REAL_SIMILAR(inst.getIonDetectors()[0].getResolution(), 0.0f)
	TEST_REAL_SIMILAR(inst.getIonDetectors()[0].getADCSamplingFrequency(), 0.0f)
	TEST_EQUAL(inst.getIonDetectors()[0].getAcquisitionMode(), IonDetector::ACQMODENULL)
	TEST_EQUAL(inst.getMassAnalyzers().size(), 1)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getType(), MassAnalyzer::PAULIONTRAP)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionMethod(), MassAnalyzer::FWHM)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionType(), MassAnalyzer::RESTYPENULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanDirection(), MassAnalyzer::SCANDIRNULL)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanLaw(), MassAnalyzer::SCANLAWNULL)
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
	// contact persons
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
	// sample
  //---------------------------------------------------------------------------
	TEST_EQUAL(e.getSample().getName(), "")
	TEST_EQUAL(e.getSample().getNumber(), "")
	TEST_EQUAL(e.getSample().getState(), Sample::SAMPLENULL)
 	TEST_EQUAL(e.getSample().getMass(), 0.0f)
	TEST_EQUAL(e.getSample().getVolume(), 0.0f)
	TEST_EQUAL(e.getSample().getConcentration(), 0.0f)
	
	//---------------------------------------------------------------------------
	// precursors
	//---------------------------------------------------------------------------
	TEST_EQUAL(e[0].getPrecursors().size(),0)
	TEST_EQUAL(e[1].getPrecursors().size(),0)
	TEST_EQUAL(e[2].getPrecursors().size(),0)
	TEST_EQUAL(e[3].getPrecursors().size(),3)
	
	TEST_REAL_SIMILAR(e[3].getPrecursors()[0].getMZ(),101.0)
	TEST_REAL_SIMILAR(e[3].getPrecursors()[0].getIntensity(),100.0)
	TEST_REAL_SIMILAR(e[3].getPrecursors()[0].getIsolationWindowLowerOffset(),96.0)
	TEST_REAL_SIMILAR(e[3].getPrecursors()[0].getIsolationWindowUpperOffset(),106.0)
	TEST_EQUAL(e[3].getPrecursors()[0].getCharge(),1)
	
	TEST_REAL_SIMILAR(e[3].getPrecursors()[1].getMZ(),201.0)
	TEST_REAL_SIMILAR(e[3].getPrecursors()[1].getIntensity(),200.0)
	TEST_REAL_SIMILAR(e[3].getPrecursors()[1].getIsolationWindowLowerOffset(),191.0)
	TEST_REAL_SIMILAR(e[3].getPrecursors()[1].getIsolationWindowUpperOffset(),211.0)
	TEST_EQUAL(e[3].getPrecursors()[1].getCharge(),2)

	TEST_REAL_SIMILAR(e[3].getPrecursors()[2].getMZ(),301.0)
	TEST_REAL_SIMILAR(e[3].getPrecursors()[2].getIntensity(),300.0)
	TEST_REAL_SIMILAR(e[3].getPrecursors()[2].getIsolationWindowLowerOffset(),286.0)
	TEST_REAL_SIMILAR(e[3].getPrecursors()[2].getIsolationWindowUpperOffset(),316.0)
	TEST_EQUAL(e[3].getPrecursors()[2].getCharge(),3)
	
	/////////////////////// TESTING SPECIAL CASES ///////////////////////

	//load a second time to make sure everything is re-initialized correctly
	MSExperiment<> e2;
	file.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"),e2);
	TEST_EQUAL(e==e2,true)

	//test reading 64 bit data
	MSExperiment<> e3;
	file.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_3_64bit.mzXML"),e3);

  TEST_EQUAL(e3.size(), 3)
	TEST_EQUAL(e3[0].getMSLevel(), 1)
	TEST_EQUAL(e3[1].getMSLevel(), 1)
	TEST_EQUAL(e3[2].getMSLevel(), 1)
	TEST_REAL_SIMILAR(e3[0].getRT(), 1)
	TEST_REAL_SIMILAR(e3[1].getRT(), 121)
	TEST_REAL_SIMILAR(e3[2].getRT(), 3661)
	TEST_EQUAL(e3[0].size(), 1)
	TEST_EQUAL(e3[1].size(), 3)
	TEST_EQUAL(e3[2].size(), 5)

	TEST_REAL_SIMILAR(e3[0][0].getPosition()[0], 120)
	TEST_REAL_SIMILAR(e3[0][0].getIntensity(), 100)
	TEST_REAL_SIMILAR(e3[1][0].getPosition()[0], 110)
	TEST_REAL_SIMILAR(e3[1][0].getIntensity(), 100)
	TEST_REAL_SIMILAR(e3[1][1].getPosition()[0], 120)
	TEST_REAL_SIMILAR(e3[1][1].getIntensity(), 200)
	TEST_REAL_SIMILAR(e3[1][2].getPosition()[0], 130)
	TEST_REAL_SIMILAR(e3[1][2].getIntensity(), 100)
	TEST_REAL_SIMILAR(e3[2][0].getPosition()[0], 100)
	TEST_REAL_SIMILAR(e3[2][0].getIntensity(), 100)
	TEST_REAL_SIMILAR(e3[2][1].getPosition()[0], 110)
	TEST_REAL_SIMILAR(e3[2][1].getIntensity(), 200)
	TEST_REAL_SIMILAR(e3[2][2].getPosition()[0], 120)
	TEST_REAL_SIMILAR(e3[2][2].getIntensity(), 300)
	TEST_REAL_SIMILAR(e3[2][3].getPosition()[0], 130)
	TEST_REAL_SIMILAR(e3[2][3].getIntensity(), 200)
	TEST_REAL_SIMILAR(e3[2][4].getPosition()[0], 140)
	TEST_REAL_SIMILAR(e3[2][4].getIntensity(), 100)

	//loading a minimal file containing one spectrum - with whitespaces inside the base64 data
	MSExperiment<> e4;
  file.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_2_minimal.mzXML"),e4);
  TEST_EQUAL(e4.size(),1)
  TEST_EQUAL(e4[0].size(),1)

	//load one extremely long spectrum - tests CDATA splitting
	MSExperiment<> e5;
	file.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_4_long.mzXML"),e5);
	TEST_EQUAL(e5.size(), 1)
	TEST_EQUAL(e5[0].size(), 997530)

	//test if it works with different peak types
	MSExperiment<RichPeak1D> e_rich;
  file.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"),e_rich);

	//zllib funcionality
	MSExperiment<> zlib;
	MSExperiment<> none;
	file.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"),none);
	file.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1_compressed.mzXML"),zlib);
	TEST_EQUAL(zlib==none,true)
END_SECTION

START_SECTION(([EXTRA] load with metadata only flag))
	TOLERANCE_ABSOLUTE(0.01)

	MSExperiment<> e;
	MzXMLFile file;
	file.getOptions().setMetadataOnly(true);

	// real test
	file.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"),e);

	TEST_EQUAL(e.size(),0)
  TEST_EQUAL(e.getSourceFiles().size(),2)
  TEST_STRING_EQUAL(e.getSourceFiles()[0].getNameOfFile(), "File_test_1.raw");
  TEST_STRING_EQUAL(e.getSourceFiles()[0].getPathToFile(), "");
	TEST_EQUAL(e.getContacts().size(),1)
	TEST_STRING_EQUAL(e.getContacts()[0].getFirstName(),"FirstName")
	TEST_STRING_EQUAL( e.getContacts()[0].getLastName(),"LastName")
	TEST_STRING_EQUAL(e.getSample().getName(), "")
	TEST_STRING_EQUAL(e.getSample().getNumber(), "")
END_SECTION

START_SECTION(([EXTRA] load with selected MS levels))
	TOLERANCE_ABSOLUTE(0.01)

	MSExperiment<> e;
	MzXMLFile file;

	// load only MS level 1
	file.getOptions().addMSLevel(1);
	file.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"),e);

	TEST_EQUAL(e.size(), 3)
	TEST_EQUAL(e[0].getMSLevel(), 1)
	TEST_EQUAL(e[1].getMSLevel(), 1)
	TEST_EQUAL(e[2].getMSLevel(), 1)
	TEST_EQUAL(e[0].size(), 1)
	TEST_EQUAL(e[1].size(), 3)
	TEST_EQUAL(e[2].size(), 5)
	TEST_STRING_EQUAL(e[0].getNativeID(),"scan=10")
	TEST_STRING_EQUAL(e[1].getNativeID(),"scan=11")
	TEST_STRING_EQUAL(e[2].getNativeID(),"scan=12")

	// load all levels
	file.getOptions().clearMSLevels();
	file.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"),e);

	TEST_EQUAL(e.size(), 4)
END_SECTION

START_SECTION(([EXTRA] load with selected MZ range))
	TOLERANCE_ABSOLUTE(0.01)

	MSExperiment<> e;
	MzXMLFile file;

	file.getOptions().setMZRange(makeRange(115,135));
	file.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"),e);
	//---------------------------------------------------------------------------
	// 60 : +(120,100)
	// 120: -(110,100) +(120,200) +(130,100)
	// 180: -(100,100) -(110,200) +(120,300) +(130,200) -(140,100)
	//---------------------------------------------------------------------------

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
END_SECTION

START_SECTION(([EXTRA] load with RT range))
	TOLERANCE_ABSOLUTE(0.01)

	MSExperiment<> e;
	MzXMLFile file;
	file.getOptions().setRTRange(makeRange(100, 200));
	file.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"),e);
	//---------------------------------------------------------------------------
	// 120: (110,100) (120,200) (130,100)
	// 180: (100,100) (110,200) (120,300) (130,200) (140,100)
	//---------------------------------------------------------------------------
 	TEST_EQUAL(e.size(), 2)
 	TEST_EQUAL(e[0].size(), 3)
 	TEST_EQUAL(e[1].size(), 5)

	TEST_REAL_SIMILAR(e[0][0].getPosition()[0], 110)
	TEST_REAL_SIMILAR(e[0][0].getIntensity(), 100)
	TEST_REAL_SIMILAR(e[0][1].getPosition()[0], 120)
	TEST_REAL_SIMILAR(e[0][1].getIntensity(), 200)
	TEST_REAL_SIMILAR(e[0][2].getPosition()[0], 130)
	TEST_REAL_SIMILAR(e[0][2].getIntensity(), 100)
	TEST_REAL_SIMILAR(e[1][0].getPosition()[0], 100)
	TEST_REAL_SIMILAR(e[1][0].getIntensity(), 100)
	TEST_REAL_SIMILAR(e[1][1].getPosition()[0], 110)
	TEST_REAL_SIMILAR(e[1][1].getIntensity(), 200)
	TEST_REAL_SIMILAR(e[1][2].getPosition()[0], 120)
	TEST_REAL_SIMILAR(e[1][2].getIntensity(), 300)
	TEST_REAL_SIMILAR(e[1][3].getPosition()[0], 130)
	TEST_REAL_SIMILAR(e[1][3].getIntensity(), 200)
	TEST_REAL_SIMILAR(e[1][4].getPosition()[0], 140)
	TEST_REAL_SIMILAR(e[1][4].getIntensity(), 100)
END_SECTION

START_SECTION(([EXTRA] load with intensity range))
	TOLERANCE_ABSOLUTE(0.01)

	MSExperiment<> e;
	MzXMLFile file;
	file.getOptions().setIntensityRange(makeRange(150, 350));
	file.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"),e);
	//---------------------------------------------------------------------------
	// 60 : -(120,100)
	// 120: -(110,100) +(120,200) -(130,100)
	// 180: -(100,100) +(110,200) +(120,300) +(130,200) -(140,100)
	//---------------------------------------------------------------------------
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
	
	

	
END_SECTION

START_SECTION(([EXTRA] load/store for nested scans))
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
END_SECTION

START_SECTION((template<typename MapType> void store(const String& filename, const MapType& map) const ))
	std::string tmp_filename;
  MSExperiment<> e1, e2;
  MzXMLFile f;

  NEW_TMP_FILE(tmp_filename);
 	 f.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"),e1);
	TEST_EQUAL(e1.size(), 4)

	f.store(tmp_filename,e1);
	f.load(tmp_filename,e2);
	TEST_EQUAL(e1==e2, true);
	
END_SECTION

START_SECTION([EXTRA] static bool isValid(const String& filename))
	std::string tmp_filename;
  MzXMLFile f;
  MSExperiment<> e;

  //Note: empty mzXML files are not valid, thus this test is omitted

	//test if full file is valid
	NEW_TMP_FILE(tmp_filename);
	f.load(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"),e);
  f.store(tmp_filename,e);
  TEST_EQUAL(f.isValid(tmp_filename),true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
