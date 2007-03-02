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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
///////////////////////////

#include <OpenMS/FORMAT/FeaturePairsFile.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/ElementPair.h>
///////////////////////////

START_TEST(FeaturePairsFile, "$Id: FeaturePairsFile_test.C 1586 2007-03-01 17:59:10Z elange $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

FeaturePairsFile* ptr = 0;
CHECK((FeaturePairsFile()))
	ptr = new FeaturePairsFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~FeaturePairsFile()))
	delete ptr;
RESULT

CHECK((void load(String filename, DFeaturePairVector<D>& pairs) throw(Exception::FileNotFound, Exception::ParseError)))
	std::vector< ElementPair < Feature > >  pvector;
	FeaturePairsFile pfile;
		   
  pfile.load("data/FeaturePairsFile.xml",pvector);
  ElementPair< Feature > pair = pvector.back();
  
  Feature first  = pair.getFirst();
  Feature second = pair.getSecond();
	
	TEST_REAL_EQUAL(first.getIntensity(),5);
  TEST_REAL_EQUAL(first.getPos()[0],0);
  TEST_REAL_EQUAL(first.getPos()[1],0);
	
  TEST_REAL_EQUAL(second.getIntensity(),0);
  TEST_REAL_EQUAL(second.getPos()[0],1.4);
  TEST_REAL_EQUAL(second.getPos()[1],2.5);
	
RESULT

CHECK((void store(String filename, const DFeaturePairVector<D>& pairs) const throw(Exception::UnableToCreateFile)))
	std::string tmp_filename;
  std::vector< ElementPair < Feature > >  pvector;
	FeaturePairsFile pfile;
  
  NEW_TMP_FILE(tmp_filename);
  pfile.load("data/FeaturePairsFile.xml",pvector);
	pfile.store(tmp_filename,pvector);
	
	TEST_FILE(tmp_filename.c_str(), "data/FeaturePairsFile.xml");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
