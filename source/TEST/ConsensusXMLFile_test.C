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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusXMLFile* ptr = 0;
CHECK((ConsensusXMLFile()))
	ptr = new ConsensusXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~ConsensusXMLFile()))
	delete ptr;
RESULT

CHECK((template<typename AlignmentT> void store(const String& filename, const AlignmentT& alignment) const throw(Exception::UnableToCreateFile)))
  std::string tmp_filename;
  ConsensusMap cons_map;
  ConsensusXMLFile cons_file;
  
  cons_file.load("data/ConsensusXMLFile.xml",cons_map);
  std::vector<LinearMapping> mapping(2);
  mapping[0].setSlope(0.5);
  mapping[0].setIntercept(-5.99959);
  mapping[1].setSlope(0.999999);
  mapping[1].setIntercept(-0.0990517);

  NEW_TMP_FILE(tmp_filename);
  cons_file.store(tmp_filename,cons_map);
  PRECISION(0.01);
	FuzzyStringComparator fsc;
	fsc.setVerboseLevel(0);
	fsc.setAcceptableRelative(1.0);
	fsc.setAcceptableAbsolute(0.0);
	bool file_is_okay = fsc.compare_files(tmp_filename.c_str(),"data/ConsensusXMLFile2.xml");
	TEST_EQUAL(file_is_okay,true);
  TEST_EQUAL(cons_file.isValid(tmp_filename),true);
RESULT

CHECK((void load(const String &filename, ConsensusMap&map, bool load_element_maps=true) throw (Exception::FileNotFound, Exception::ParseError)))
  ConsensusMap cons_map;
  ConsensusXMLFile cons_file;
  cons_file.load("data/ConsensusXMLFile.xml", cons_map);
  TEST_EQUAL(cons_map.getFileNames()[0] == "data/MapAlignmentFeatureMap1.xml", true)
  TEST_EQUAL(cons_map.getFileNames()[1] == "data/MapAlignmentFeatureMap2.xml", true)

  ConsensusFeature cons_feature = cons_map[0];
  TEST_REAL_EQUAL(cons_feature.getPosition()[0],1273.27)  
  TEST_REAL_EQUAL(cons_feature.getPosition()[1],904.47)
  TEST_REAL_EQUAL(cons_feature.getIntensity(),3.12539e+07)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[0],1273.27)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[0],1273.27)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[1],904.47)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[1],904.47)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().min()[0],3.12539e+07)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().max()[0],3.12539e+07)
  Group::const_iterator it = cons_feature.begin();
//  TEST_REAL_EQUAL(it->getElement().getPosition()[0],1273.27)  
//  TEST_REAL_EQUAL(it->getElement().getPosition()[1],904.47)
  TEST_REAL_EQUAL(it->getIntensity(),3.12539e+07)
    
  cons_feature = cons_map[5];
  TEST_REAL_EQUAL(cons_feature.getPosition()[0],1194.82)  
  TEST_REAL_EQUAL(cons_feature.getPosition()[1],777.101)
  TEST_REAL_EQUAL(cons_feature.getIntensity(),1.78215e+07)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[0],1194.82)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[0],1194.82)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[1],777.101)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[1],777.101)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().min()[0],1.78215e+07)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().max()[0],1.78215e+07)
  it = cons_feature.begin();
//  TEST_REAL_EQUAL(it->getElement().getPosition()[0],1194.82)  
//  TEST_REAL_EQUAL(it->getElement().getPosition()[1],777.101)
  TEST_REAL_EQUAL(it->getIntensity(),1.78215e+07)
  ++it;
//  TEST_REAL_EQUAL(it->getElement().getPosition()[0],2401.64)  
//  TEST_REAL_EQUAL(it->getElement().getPosition()[1],777.201)
  TEST_REAL_EQUAL(it->getIntensity(),1.78215e+07)
RESULT

CHECK(static bool isValid(const String& filename))
	//tested above
	NOT_TESTABLE;
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



