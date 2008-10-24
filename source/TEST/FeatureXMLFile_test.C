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
// $Maintainer: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
///////////////////////////

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>

using namespace OpenMS;
using namespace std;

DRange<1> makeRange(float a, float b)
{
  DPosition<1> pa(a), pb(b);
  return DRange<1>(pa, pb);
}

///////////////////////////

START_TEST(FeatureXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureXMLFile* ptr = 0;
CHECK((FeatureXMLFile()))
	ptr = new FeatureXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~FeatureXMLFile()))
	delete ptr;
RESULT
 
CHECK((void load(String filename, FeatureMap<>& feature_map)))
	PRECISION(0.01)
	
	FeatureMap<> e;
	FeatureXMLFile dfmap_file;
	
	//test exception
	TEST_EXCEPTION( Exception::FileNotFound , dfmap_file.load("dummy/dummy.MzData",e) )
	
	// real test
	dfmap_file.load("data/FeatureXMLFile.xml",e);

  // id
	TEST_EQUAL(e.getIdentifier(),"lsid");
	
	TEST_EQUAL(e.size(),2)
	TEST_REAL_EQUAL(e[0].getRT(), 25)
	TEST_REAL_EQUAL(e[0].getMZ(), 0)
	TEST_REAL_EQUAL(e[0].getIntensity(), 300)
	TEST_EQUAL(e[0].getMetaValue("stringparametername"),"stringparametervalue")
	TEST_EQUAL((UInt)e[0].getMetaValue("intparametername"),4)
	TEST_REAL_EQUAL((DoubleReal)e[0].getMetaValue("floatparametername"),4.551)

	TEST_REAL_EQUAL(e[1].getRT(), 0)
	TEST_REAL_EQUAL(e[1].getMZ(), 35)
	TEST_REAL_EQUAL(e[1].getIntensity(), 500)

	//PeakFileOptions tests
	dfmap_file.getOptions().setRTRange(makeRange(0, 10));
	dfmap_file.load("data/FeatureXMLFile.xml",e);
	TEST_EQUAL(e.size(),1)
	TEST_REAL_EQUAL(e[0].getRT(), 0)
	TEST_REAL_EQUAL(e[0].getMZ(), 35)
	TEST_REAL_EQUAL(e[0].getIntensity(), 500)

	dfmap_file.getOptions() = PeakFileOptions();
	dfmap_file.getOptions().setMZRange(makeRange(10, 50));
	dfmap_file.load("data/FeatureXMLFile.xml",e);
	TEST_EQUAL(e.size(),1)
	TEST_REAL_EQUAL(e[0].getRT(), 0)
	TEST_REAL_EQUAL(e[0].getMZ(), 35)
	TEST_REAL_EQUAL(e[0].getIntensity(), 500)

	dfmap_file.getOptions() = PeakFileOptions();
	dfmap_file.getOptions().setIntensityRange(makeRange(400, 600));
	dfmap_file.load("data/FeatureXMLFile.xml",e);
	TEST_EQUAL(e.size(),1)
	TEST_REAL_EQUAL(e[0].getRT(), 0)
	TEST_REAL_EQUAL(e[0].getMZ(), 35)
	TEST_REAL_EQUAL(e[0].getIntensity(), 500)
RESULT

CHECK((void store(String filename, const FeatureMap<>& feature_map) const))
  
  std::string tmp_filename;
  FeatureMap<> e;
  FeatureXMLFile f;
  
  NEW_TMP_FILE(tmp_filename);
  f.load("data/FeatureXMLFile.xml",e);
  f.store(tmp_filename,e);
  TEST_FILE(tmp_filename.c_str(),"data/FeatureXMLFile.xml");

RESULT

CHECK( PeakFileOptions& getOptions() )
	FeatureXMLFile f;
  FeatureMap<> e;
	f.getOptions().setRTRange(makeRange(1.5, 4.5));
	f.load("data/FeatureXMLFile2.xml",e);
	TEST_EQUAL(e.size(), 5)

	f.getOptions().setMZRange(makeRange(1025.0, 2000.0));
	f.load("data/FeatureXMLFile2.xml",e);
	TEST_EQUAL(e.size(), 3)

	f.getOptions().setIntensityRange(makeRange(290.0, 310.0));
	f.load("data/FeatureXMLFile2.xml",e);
	TEST_EQUAL(e.size(), 1)
	
	f.getOptions().setMetadataOnly(true);
	f.load("data/FeatureXMLFile2.xml",e);
	TEST_EQUAL(e.getIdentifier(), "lsid2")
	TEST_EQUAL(e.size(), 0)
RESULT

CHECK([EXTRA] static bool isValid(const String& filename))
	FeatureMap<> e;
  FeatureXMLFile f;
	String filename;
	
  //test if empty file is valid
	NEW_TMP_FILE(filename)
	f.store(filename,e);	
  TEST_EQUAL(f.isValid(filename),true);	
	
	//test if full file is valid
	NEW_TMP_FILE(filename);
	f.load("data/FeatureXMLFile.xml",e);
	f.store(filename, e);	
  TEST_EQUAL(f.isValid(filename),true);
RESULT

CHECK( const PeakFileOptions& getOptions() const )
	FeatureXMLFile f;
 	FeatureMap<> e;
	f.getOptions().setRTRange(makeRange(1.5, 4.5));
	f.getOptions().setIntensityRange(makeRange(290.0, 310.0));
	
	const PeakFileOptions pfo = f.getOptions();
	
	TEST_EQUAL(pfo.getRTRange(),makeRange(1.5, 4.5))
	TEST_EQUAL(pfo.getIntensityRange(),makeRange(290.0, 310.0))
	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
 
