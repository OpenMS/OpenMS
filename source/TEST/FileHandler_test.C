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

#include <OpenMS/FORMAT/FileHandler.h>

///////////////////////////

START_TEST(FileHandler, "FileHandler")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

CHECK(String typeToName(Type type))
	FileHandler tmp;
	TEST_EQUAL("Unknown", tmp.typeToName(FileHandler::UNKNOWN));
	TEST_EQUAL("DTA", tmp.typeToName(FileHandler::DTA));
	TEST_EQUAL("DTA2D", tmp.typeToName(FileHandler::DTA2D));
	TEST_EQUAL("mzData", tmp.typeToName(FileHandler::MZDATA));
	TEST_EQUAL("mzXML", tmp.typeToName(FileHandler::MZXML));
	TEST_EQUAL("FeatureFile", tmp.typeToName(FileHandler::FEATURE));
	TEST_EQUAL("ANDIMS", tmp.typeToName(FileHandler::ANDIMS));
	TEST_EQUAL("FeaturePairs", tmp.typeToName(FileHandler::FEATURE_PAIRS));
RESULT

CHECK(Type nameToType(const String& name))
	FileHandler tmp;
	TEST_EQUAL(FileHandler::UNKNOWN, tmp.nameToType("Unknown"));
	TEST_EQUAL(FileHandler::DTA, tmp.nameToType("DTA"));
	TEST_EQUAL(FileHandler::DTA2D, tmp.nameToType("DTA2D"));
	TEST_EQUAL(FileHandler::MZDATA, tmp.nameToType("mzData"));
	TEST_EQUAL(FileHandler::MZXML, tmp.nameToType("mzXML"));
	TEST_EQUAL(FileHandler::FEATURE, tmp.nameToType("FeatureFile"));
	TEST_EQUAL(FileHandler::FEATURE_PAIRS, tmp.nameToType("FeaturePairs"));
	TEST_EQUAL(FileHandler::ANDIMS, tmp.nameToType("ANDIMS"));
	TEST_EQUAL(FileHandler::ANDIMS, tmp.nameToType("aNdIMs"));
RESULT

CHECK(Type getTypeByFileName(const String& filename))
	FileHandler tmp;
	TEST_EQUAL(tmp.getTypeByFileName("test.bla"), FileHandler::UNKNOWN)
	TEST_EQUAL(tmp.getTypeByFileName("test.dta"), FileHandler::DTA)
	TEST_EQUAL(tmp.getTypeByFileName("test.MzData"), FileHandler::MZDATA)
	TEST_EQUAL(tmp.getTypeByFileName("test.DTA2D"), FileHandler::DTA2D)
	TEST_EQUAL(tmp.getTypeByFileName("test.feat"), FileHandler::FEATURE)
	TEST_EQUAL(tmp.getTypeByFileName("test.pairs"), FileHandler::FEATURE_PAIRS)
	TEST_EQUAL(tmp.getTypeByFileName("test.MZXML"), FileHandler::MZXML)
	TEST_EQUAL(tmp.getTypeByFileName("test.cdf"), FileHandler::ANDIMS)
RESULT

CHECK(Type getTypeByContent(const String& filename) throw(Exception::FileNotFound))
	FileHandler tmp;
	TEST_EQUAL(tmp.getTypeByContent("data/MzDataFile_test_1.mzData"), FileHandler::MZDATA)
	TEST_EQUAL(tmp.getTypeByContent("data/DFeatureMapFile.xml"), FileHandler::FEATURE)
	TEST_EQUAL(tmp.getTypeByContent("data/DFeaturePairsFile.xml"), FileHandler::FEATURE_PAIRS)
	TEST_EQUAL(tmp.getTypeByContent("data/MzXMLFile_test_1.mzXML"), FileHandler::MZXML)
	TEST_EQUAL(tmp.getTypeByContent("data/DTAFile_test.dta"), FileHandler::DTA)
	TEST_EQUAL(tmp.getTypeByContent("data/DTA2DFile_test_1.dta2d"), FileHandler::DTA2D)
	TEST_EQUAL(tmp.getTypeByContent("data/DTA2DFile_test_2.dta2d"), FileHandler::DTA2D)
	TEST_EQUAL(tmp.getTypeByContent("data/ANDIFile_test.cdf"), FileHandler::ANDIMS)
	TEST_EQUAL(tmp.getTypeByContent("data/class_test_infile.txt"), FileHandler::UNKNOWN)
RESULT

CHECK(template<class PeakType> bool loadExperiment(const String& filename, MSExperiment<PeakType>& exp, Type force_type = UNKNOWN))
	FileHandler tmp;
	MSExperiment<> exp;
	TEST_EQUAL(tmp.loadExperiment("test.bla",exp), false)	
	TEST_EQUAL(tmp.loadExperiment("data/DTAFile_test.dta",exp), true)
	TEST_EQUAL(tmp.loadExperiment("data/MzDataFile_test_1.mzData",exp), true)	
  TEST_EQUAL(tmp.loadExperiment("data/MzXMLFile_test_1.mzXML",exp), true)	
#ifdef ANDIMS_DEF
  TEST_EQUAL(tmp.loadExperiment("data/ANDIFile_test.cdf",exp), true)
#else
	TEST_EQUAL(tmp.loadExperiment("data/ANDIFile_test.cdf",exp), false)
#endif
  TEST_EQUAL(tmp.loadExperiment("data/DTA2DFile_test_1.dta2d",exp), true)	
	TEST_EXCEPTION(Exception::ParseError,tmp.loadExperiment("data/DTAFile_test.dta",exp, FileHandler::DTA2D))
RESULT

CHECK(bool isSupported(Type type))
	FileHandler tmp;
	TEST_EQUAL(false, tmp.isSupported(FileHandler::UNKNOWN));
	TEST_EQUAL(true, tmp.isSupported(FileHandler::DTA));
	TEST_EQUAL(true, tmp.isSupported(FileHandler::DTA2D));
	TEST_EQUAL(true, tmp.isSupported(FileHandler::MZDATA));
  TEST_EQUAL(true, tmp.isSupported(FileHandler::MZXML));
	TEST_EQUAL(true, tmp.isSupported(FileHandler::FEATURE));
	TEST_EQUAL(true, tmp.isSupported(FileHandler::FEATURE_PAIRS));
#ifdef ANDIMS_DEF
  TEST_EQUAL(true, tmp.isSupported(FileHandler::ANDIMS));
#else
	TEST_EQUAL(false, tmp.isSupported(FileHandler::ANDIMS));
#endif
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
