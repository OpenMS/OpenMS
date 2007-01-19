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

#include <OpenMS/FORMAT/DFeaturePairsFile.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
///////////////////////////

START_TEST(DFeaturePairsFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

DFeaturePairsFile* ptr = 0;
CHECK((DFeaturePairsFile()))
	ptr = new DFeaturePairsFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~DFeaturePairsFile()))
	delete ptr;
RESULT

CHECK((template<Size D> void load(String filename, DFeaturePairVector<D>& pairs) throw(Exception::FileNotFound, Exception::ParseError)))
	PRECISION(0.01)
	
	DFeaturePairVector<2> pvector;
	DFeaturePairsFile pfile;
		   
  pfile.load("data/DFeaturePairsFile.xml",pvector);
  DFeaturePair<2> pair = pvector.back();
  
  DFeature<2> first  = pair.getFirst();
  DFeature<2> second = pair.getSecond();
	
	TEST_EQUAL(first.getIntensity(),5);
	TEST_EQUAL(first.getPosition()[0],0);
	TEST_EQUAL(first.getPosition()[1],0);
	
	TEST_EQUAL(second.getIntensity(),0);
	TEST_EQUAL(second.getPosition()[0],1.4);
	TEST_EQUAL(second.getPosition()[1],2.5);
	
RESULT

CHECK((template<Size D> void store(String filename, const DFeaturePairVector<D>& pairs) const throw(Exception::UnableToCreateFile)))
	
	std::string tmp_filename;
  DFeaturePairVector<2> pvector;
	DFeaturePairsFile pfile;
  
  NEW_TMP_FILE(tmp_filename);
  pfile.load("data/DFeaturePairsFile.xml",pvector);
	pfile.store(tmp_filename,pvector);
	
	TEST_FILE(tmp_filename.c_str(), "data/DFeaturePairsFile.xml");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
