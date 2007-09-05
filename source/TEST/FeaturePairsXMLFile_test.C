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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
///////////////////////////

#include <OpenMS/FORMAT/FeaturePairsXMLFile.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ElementPair.h>
///////////////////////////

START_TEST(FeaturePairsXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

FeaturePairsXMLFile* ptr = 0;
CHECK((FeaturePairsXMLFile()))
	ptr = new FeaturePairsXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~FeaturePairsXMLFile()))
	delete ptr;
RESULT

CHECK((void load(String filename, std::vector< ElementPair< Feature > > &pairs) throw (Exception::FileNotFound, Exception::ParseError)))
	std::vector< ElementPair < Feature > >  pvector;
	FeaturePairsXMLFile pfile;
		   
  pfile.load("data/FeaturePairsXMLFile.xml",pvector);
  ElementPair< Feature > pair = pvector.back();
  
  Feature first  = pair.getFirst();
  Feature second = pair.getSecond();
	
	TEST_REAL_EQUAL(first.getIntensity(),5);
  TEST_REAL_EQUAL(first.getPosition()[0],0);
  TEST_REAL_EQUAL(first.getPosition()[1],0);
	
  TEST_REAL_EQUAL(second.getIntensity(),0);
  TEST_REAL_EQUAL(second.getPosition()[0],1.4);
  TEST_REAL_EQUAL(second.getPosition()[1],2.5);
	
RESULT

CHECK((void store(String filename, const std::vector< ElementPair< Feature > > &pairs) const throw (Exception::UnableToCreateFile)))
	std::string tmp_filename;
  std::vector< ElementPair < Feature > >  pvector;
	FeaturePairsXMLFile pfile;
  
  NEW_TMP_FILE(tmp_filename);
  pfile.load("data/FeaturePairsXMLFile.xml",pvector);
	pfile.store(tmp_filename,pvector);
	
	TEST_FILE(tmp_filename.c_str(), "data/FeaturePairsXMLFile.xml");
RESULT

CHECK((static void pairsToFeatures(const std::vector< ElementPair< Feature > > &pairs, FeatureMap<> &map)))
  FeatureMap<> map;
  Feature feat_1, feat_2;
  feat_1.setIntensity(10);
  feat_2.setIntensity(20);
  ElementPair< Feature > pair(feat_1,feat_2);
  std::vector< ElementPair< Feature > > elem_vec;
  elem_vec.push_back(pair);
    
  FeaturePairsXMLFile fp;
  fp.pairsToFeatures(elem_vec,map);
  
  TEST_EQUAL(map[0] == feat_1,true) 
  TEST_EQUAL(map[1] == feat_2,true)
RESULT

CHECK(static bool isValid(const String& filename))
	std::vector< ElementPair< Feature > > e;
	String filename;
	FeaturePairsXMLFile f;

  //test if empty file is valid
	NEW_TMP_FILE(filename)
	f.store(filename, e);	
  TEST_EQUAL(f.isValid(filename),true);	
	
	//test if full file is valid
	NEW_TMP_FILE(filename);
	f.load("data/FeaturePairsXMLFile.xml",e);
	f.store(filename, e);	
  TEST_EQUAL(f.isValid(filename),true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
