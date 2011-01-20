// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

///////////////////////////

START_TEST(FileHandler, "FileHandler")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_SECTION((static String typeToName(FileTypes::Type type)))
	FileHandler tmp;
	TEST_EQUAL(tmp.typeToName(FileTypes::UNKNOWN),"Unknown");
	TEST_EQUAL(tmp.typeToName(FileTypes::DTA),"DTA");
	TEST_EQUAL(tmp.typeToName(FileTypes::DTA2D),"DTA2D");
	TEST_EQUAL(tmp.typeToName(FileTypes::MZDATA),"mzData");
	TEST_EQUAL(tmp.typeToName(FileTypes::MZXML),"mzXML");
	TEST_EQUAL(tmp.typeToName(FileTypes::MZML),"mzML");
	TEST_EQUAL(tmp.typeToName(FileTypes::FEATUREXML),"FeatureXML");
	TEST_EQUAL(tmp.typeToName(FileTypes::ANDIMS),"cdf");
	TEST_EQUAL(tmp.typeToName(FileTypes::IDXML),"IdXML");
	TEST_EQUAL(tmp.typeToName(FileTypes::CONSENSUSXML),"ConsensusXML");
	TEST_EQUAL(tmp.typeToName(FileTypes::TRANSFORMATIONXML),"TrafoXML");
	TEST_EQUAL(tmp.typeToName(FileTypes::INI),"ini");
	TEST_EQUAL(tmp.typeToName(FileTypes::PNG),"PNG");
END_SECTION

START_SECTION((static FileTypes::Type nameToType(const String &name)))
	FileHandler tmp;
	TEST_EQUAL(FileTypes::UNKNOWN, tmp.nameToType("Unknown"));
	TEST_EQUAL(FileTypes::DTA, tmp.nameToType("DTA"));
	TEST_EQUAL(FileTypes::DTA2D, tmp.nameToType("DTA2D"));
	TEST_EQUAL(FileTypes::MZDATA, tmp.nameToType("mzData"));
	TEST_EQUAL(FileTypes::MZXML, tmp.nameToType("mzXML"));
	TEST_EQUAL(FileTypes::FEATUREXML, tmp.nameToType("FeatureXML"));
	TEST_EQUAL(FileTypes::ANDIMS, tmp.nameToType("cdf"));
	TEST_EQUAL(FileTypes::ANDIMS, tmp.nameToType("CdF"));
	TEST_EQUAL(FileTypes::IDXML, tmp.nameToType("IdXmL"));
	TEST_EQUAL(FileTypes::CONSENSUSXML, tmp.nameToType("ConsensusXMl"));
  TEST_EQUAL(FileTypes::MGF, tmp.nameToType("mgf"));
  TEST_EQUAL(FileTypes::INI, tmp.nameToType("ini"));
	TEST_EQUAL(FileTypes::TRANSFORMATIONXML, tmp.nameToType("TrafoXML"));
	TEST_EQUAL(FileTypes::MZML, tmp.nameToType("mzML"));
  TEST_EQUAL(FileTypes::MS2, tmp.nameToType("ms2"));
	TEST_EQUAL(FileTypes::PEPXML, tmp.nameToType("pepXML"));
	TEST_EQUAL(FileTypes::PROTXML, tmp.nameToType("protXML"));
	TEST_EQUAL(FileTypes::MZIDENTML, tmp.nameToType("mzIdentML"));
	TEST_EQUAL(FileTypes::GELML, tmp.nameToType("GelML"));
	TEST_EQUAL(FileTypes::TRAML, tmp.nameToType("TraML"));
	TEST_EQUAL(FileTypes::MSP, tmp.nameToType("MSP"));
	TEST_EQUAL(FileTypes::OMSSAXML, tmp.nameToType("OMSSAXML"));
  TEST_EQUAL(FileTypes::PNG, tmp.nameToType("PNG"));
  TEST_EQUAL(FileTypes::XMASS, tmp.nameToType("fid"));
	TEST_EQUAL(FileTypes::TSV, tmp.nameToType("tsv"));
	TEST_EQUAL(FileTypes::PEPLIST, tmp.nameToType("pepList"));
	TEST_EQUAL(FileTypes::HARDKLOER, tmp.nameToType("hardkloer"));
	TEST_EQUAL(FileTypes::KROENIK, tmp.nameToType("kroenik"));
	TEST_EQUAL(FileTypes::FASTA, tmp.nameToType("fasta"));
	TEST_EQUAL(FileTypes::EDTA, tmp.nameToType("edta"));
END_SECTION


START_SECTION((static FileTypes::Type getTypeByFileName(const String &filename)))
	FileHandler tmp;
	TEST_EQUAL(tmp.getTypeByFileName("test.bla"), FileTypes::UNKNOWN)
	TEST_EQUAL(tmp.getTypeByFileName("test.dta"), FileTypes::DTA)
	TEST_EQUAL(tmp.getTypeByFileName("test.DTA2D"), FileTypes::DTA2D)
	TEST_EQUAL(tmp.getTypeByFileName("test.MzData"), FileTypes::MZDATA)
	TEST_EQUAL(tmp.getTypeByFileName("test.MZXML"), FileTypes::MZXML)
	TEST_EQUAL(tmp.getTypeByFileName("test.featureXML"), FileTypes::FEATUREXML)
	TEST_EQUAL(tmp.getTypeByFileName("test.cDf"), FileTypes::ANDIMS)
	TEST_EQUAL(tmp.getTypeByFileName("test.idXML"), FileTypes::IDXML)
	TEST_EQUAL(tmp.getTypeByFileName("test.consensusXML"), FileTypes::CONSENSUSXML)
	TEST_EQUAL(tmp.getTypeByFileName("test.mGf"), FileTypes::MGF)
	TEST_EQUAL(tmp.getTypeByFileName("test.ini"), FileTypes::INI)
  TEST_EQUAL(tmp.getTypeByFileName("test.TraFoXML"), FileTypes::TRANSFORMATIONXML)
	TEST_EQUAL(tmp.getTypeByFileName("test.MzML"), FileTypes::MZML)
	TEST_EQUAL(tmp.getTypeByFileName(OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.bz2")), FileTypes::MZML)
	TEST_EQUAL(tmp.getTypeByFileName(OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.gz")), FileTypes::MZML)
	TEST_EQUAL(tmp.getTypeByFileName("test.mS2"), FileTypes::MS2)
	TEST_EQUAL(tmp.getTypeByFileName("test.pepXML"), FileTypes::PEPXML)
 	TEST_EQUAL(tmp.getTypeByFileName("test.pep.xml"), FileTypes::PEPXML)
	TEST_EQUAL(tmp.getTypeByFileName("test.protXML"), FileTypes::PROTXML)
	TEST_EQUAL(tmp.getTypeByFileName("test.prot.xml"), FileTypes::PROTXML)
	TEST_EQUAL(tmp.getTypeByFileName("test.mzIdentML"), FileTypes::MZIDENTML)
	TEST_EQUAL(tmp.getTypeByFileName("test.GELML"), FileTypes::GELML)
	TEST_EQUAL(tmp.getTypeByFileName("test.TRAML"), FileTypes::TRAML)
	TEST_EQUAL(tmp.getTypeByFileName("test.MSP"), FileTypes::MSP)
	TEST_EQUAL(tmp.getTypeByFileName("test.OMSSAXML"), FileTypes::OMSSAXML)
  TEST_EQUAL(tmp.getTypeByFileName("test.png"), FileTypes::PNG)
	TEST_EQUAL(tmp.getTypeByFileName("./foo.bar/XMass/fid"), FileTypes::XMASS);
	TEST_EQUAL(tmp.getTypeByFileName("test.TSV"), FileTypes::TSV)
	TEST_EQUAL(tmp.getTypeByFileName("test.PEPLIST"), FileTypes::PEPLIST)
	TEST_EQUAL(tmp.getTypeByFileName("test.HARDKLOER"), FileTypes::HARDKLOER)
	TEST_EQUAL(tmp.getTypeByFileName("test.fasta"), FileTypes::FASTA)
	TEST_EQUAL(tmp.getTypeByFileName("test.EDTA"), FileTypes::EDTA)
END_SECTION

START_SECTION((static FileTypes::Type getTypeByContent(const String &filename)))
	FileHandler tmp;
	TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData")), FileTypes::MZDATA)
	TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML")), FileTypes::FEATUREXML)
	TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML")), FileTypes::MZXML)
	TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML")), FileTypes::MZML)
	TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.bz2")), FileTypes::MZML)
	TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.gz")), FileTypes::MZML)
	TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("DTAFile_test.dta")), FileTypes::DTA)
	TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d")), FileTypes::DTA2D)
	TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_2.dta2d")), FileTypes::DTA2D)
	TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("ANDIFile_test.cdf")), FileTypes::ANDIMS)
	TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("class_test_infile.txt")), FileTypes::UNKNOWN)
	TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML")), FileTypes::IDXML)
	TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_1.trafoXML")), FileTypes::TRANSFORMATIONXML)
	TEST_EQUAL(tmp.getTypeByContent(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta")), FileTypes::FASTA)

	TEST_EXCEPTION(Exception::FileNotFound,tmp.getTypeByContent("/bli/bla/bluff"))
END_SECTION

START_SECTION((static FileTypes::Type getType(const String &filename)))
	FileHandler tmp;
	TEST_EQUAL(tmp.getType(OPENMS_GET_TEST_DATA_PATH("class_test_infile.txt")), FileTypes::UNKNOWN)
	TEST_EQUAL(tmp.getType(OPENMS_GET_TEST_DATA_PATH("IdXMLFile_whole.idXML")), FileTypes::IDXML)
	TEST_EQUAL(tmp.getType(OPENMS_GET_TEST_DATA_PATH("ConsensusXMLFile.consensusXML")), FileTypes::CONSENSUSXML)
	TEST_EQUAL(tmp.getType(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_1.trafoXML")), FileTypes::TRANSFORMATIONXML)

	TEST_EXCEPTION(Exception::FileNotFound,tmp.getType("/bli/bla/bluff"))
END_SECTION

START_SECTION((template < class PeakType > bool loadExperiment(const String &filename, MSExperiment< PeakType > &exp, FileTypes::Type force_type=FileTypes::UNKNOWN, ProgressLogger::LogType log=ProgressLogger::NONE)))
	FileHandler tmp;
	MSExperiment<> exp;
	TEST_EQUAL(tmp.loadExperiment("test.bla",exp), false)
	TEST_EQUAL(tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("DTAFile_test.dta"),exp), true)

	TEST_EQUAL(tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"),exp), true)
	TEST_REAL_SIMILAR(exp[1][0].getPosition()[0], 110)
	TEST_REAL_SIMILAR(exp[1][1].getPosition()[0], 120)
	TEST_REAL_SIMILAR(exp[1][2].getPosition()[0], 130)

  // starts with 110, so this one should skip the first
  tmp.getOptions().setMZRange(DRange<1> (115, 1000));
	TEST_EQUAL(tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("MzDataFile_1.mzData"),exp), true)
	TEST_REAL_SIMILAR(exp[1][0].getPosition()[0], 120)
	TEST_REAL_SIMILAR(exp[1][1].getPosition()[0], 130)

  tmp.getOptions() = PeakFileOptions();
  TEST_EQUAL(tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"),exp), true)
	TEST_REAL_SIMILAR(exp[2][0].getPosition()[0], 100)
	TEST_REAL_SIMILAR(exp[2][1].getPosition()[0], 110)
	TEST_REAL_SIMILAR(exp[2][2].getPosition()[0], 120)

  tmp.getOptions().setMZRange(DRange<1> (115, 1000));
  TEST_EQUAL(tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("MzXMLFile_1.mzXML"),exp), true)
	TEST_REAL_SIMILAR(exp[2][0].getPosition()[0], 120)
	TEST_REAL_SIMILAR(exp[2][1].getPosition()[0], 130)
	TEST_REAL_SIMILAR(exp[2][2].getPosition()[0], 140)

  tmp.getOptions() = PeakFileOptions();
  TEST_EQUAL(tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),exp), true)
	TEST_EQUAL(exp.size(),4)


#ifdef USE_ANDIMS
  TEST_EQUAL(tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("ANDIFile_test.cdf"),exp), true)
#else
	TEST_EQUAL(tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("ANDIFile_test.cdf"),exp), false)
#endif

  tmp.getOptions() = PeakFileOptions();
  TEST_EQUAL(tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"),exp), true)
	TEST_REAL_SIMILAR(exp[0][0].getPosition()[0], 230.02)
	TEST_REAL_SIMILAR(exp[0][1].getPosition()[0], 430.02)
	TEST_REAL_SIMILAR(exp[0][2].getPosition()[0], 630.02)

  tmp.getOptions().setMZRange(DRange<1> (300, 1000));
  TEST_EQUAL(tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("DTA2DFile_test_1.dta2d"),exp), true)
	TEST_REAL_SIMILAR(exp[0][0].getPosition()[0], 430.02)
	TEST_REAL_SIMILAR(exp[0][1].getPosition()[0], 630.02)
	
	TEST_EQUAL(tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("XMassFile_test/fid"),exp), true)

	TEST_EXCEPTION(Exception::ParseError,tmp.loadExperiment(OPENMS_GET_TEST_DATA_PATH("DTAFile_test.dta"),exp, FileTypes::DTA2D))
END_SECTION

START_SECTION((static bool isSupported(FileTypes::Type type)))
	FileHandler tmp;
	TEST_EQUAL(false, tmp.isSupported(FileTypes::UNKNOWN));
	TEST_EQUAL(true, tmp.isSupported(FileTypes::DTA));
	TEST_EQUAL(true, tmp.isSupported(FileTypes::DTA2D));
	TEST_EQUAL(true, tmp.isSupported(FileTypes::MZDATA));
	TEST_EQUAL(true, tmp.isSupported(FileTypes::MZML));
  TEST_EQUAL(true, tmp.isSupported(FileTypes::MZXML));
  TEST_EQUAL(true, tmp.isSupported(FileTypes::XMASS));
	TEST_EQUAL(true, tmp.isSupported(FileTypes::FEATUREXML));
#ifdef USE_ANDIMS
  TEST_EQUAL(true, tmp.isSupported(FileTypes::ANDIMS));
#else
	TEST_EQUAL(false, tmp.isSupported(FileTypes::ANDIMS));
#endif
END_SECTION

START_SECTION((const PeakFileOptions& getOptions() const))
	FileHandler a;
	TEST_EQUAL(a.getOptions().hasMSLevels(),false)
END_SECTION

START_SECTION((PeakFileOptions& getOptions()))
	FileHandler a;
	a.getOptions().addMSLevel(1);
	TEST_EQUAL(a.getOptions().hasMSLevels(),true);
END_SECTION

START_SECTION((template < class FeatureType > bool loadFeatures(const String &filename, FeatureMap< FeatureType > &map, FileTypes::Type force_type=FileTypes::UNKNOWN)))
  FileHandler tmp;
	FeatureMap<> map;
	TEST_EQUAL(tmp.loadFeatures("test.bla",map), false)
	TEST_EQUAL(tmp.loadFeatures(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"),map), true)
	TEST_EQUAL(map.size(),7);
	TEST_EQUAL(tmp.loadFeatures(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"),map), true)
	TEST_EQUAL(map.size(),7);
END_SECTION

START_SECTION((template <class PeakType> void storeExperiment(const String& filename, const MSExperiment<PeakType>& exp, ProgressLogger::LogType log = ProgressLogger::NONE)))
	FileHandler fh;
	MSExperiment<> exp;
	fh.loadExperiment(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"),exp);
	
	//test mzML
	String filename;
	NEW_TMP_FILE(filename)
	fh.storeExperiment(filename,exp);
	TEST_EQUAL(fh.getTypeByContent(filename),FileTypes::MZML)
		
	//other types cannot be tested, because the NEW_TMP_FILE template does not support file extensions...
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
