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
// $Maintainer: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

DRange<1> makeRange(float a, float b)
{
  DPosition<1> pa(a), pb(b);
  return DRange<1>(pa, pb);
}

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

PRECISION(0.01)

CHECK((void load(const String &filename, ConsensusMap &map) throw (Exception::FileNotFound, Exception::ParseError)))
  ConsensusMap cons_map;
  ConsensusXMLFile cons_file;
  cons_file.load("data/ConsensusXMLFile.xml", cons_map);
  TEST_EQUAL(cons_map.getFileDescriptions()[0].filename == "data/MapAlignmentFeatureMap1.xml", true)
  TEST_EQUAL(cons_map.getFileDescriptions()[0].label,"label")
  TEST_EQUAL(cons_map.getFileDescriptions()[0].size, 144)
  TEST_STRING_EQUAL(cons_map.getFileDescriptions()[1].filename,"data/MapAlignmentFeatureMap2.xml")
  TEST_EQUAL(cons_map.getFileDescriptions()[1].label,"")
  TEST_EQUAL(cons_map.getFileDescriptions()[1].size, 0)

  ConsensusFeature cons_feature = cons_map[0];
  TEST_REAL_EQUAL(cons_feature.getRT(),1273.27)  
  TEST_REAL_EQUAL(cons_feature.getMZ(),904.47)
  TEST_REAL_EQUAL(cons_feature.getIntensity(),3.12539e+07)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[0],1273.27)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[0],1273.27)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[1],904.47)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[1],904.47)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().min()[0],3.12539e+07)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().max()[0],3.12539e+07)
  TEST_REAL_EQUAL(cons_feature.getQuality(),1.1)
  ConsensusFeature::HandleSetType::const_iterator it = cons_feature.begin();
  TEST_REAL_EQUAL(it->getIntensity(),3.12539e+07)
  
  cons_feature = cons_map[5];
  TEST_REAL_EQUAL(cons_feature.getRT(),1194.82)  
  TEST_REAL_EQUAL(cons_feature.getMZ(),777.101)
  TEST_REAL_EQUAL(cons_feature.getIntensity(),1.78215e+07)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[0],1194.82)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[0],1194.82)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[1],777.101)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[1],777.101)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().min()[0],1.78215e+07)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().max()[0],1.78215e+07)
  TEST_REAL_EQUAL(cons_feature.getQuality(),0.0)
  it = cons_feature.begin();
  TEST_REAL_EQUAL(it->getIntensity(),1.78215e+07)
  ++it;
  TEST_REAL_EQUAL(it->getIntensity(),1.78215e+07)


	//PeakFileOptions tests
	
	cons_file.getOptions().setRTRange(makeRange(815, 818));
	cons_file.load("data/ConsensusXMLFile3.xml",cons_map);
	TEST_EQUAL(cons_map.size(),1)
	TEST_REAL_EQUAL(cons_map[0].getRT(), 817.266)

	cons_file.getOptions() = PeakFileOptions();
	cons_file.getOptions().setMZRange(makeRange(944, 945));
	cons_file.load("data/ConsensusXMLFile3.xml",cons_map);
	TEST_EQUAL(cons_map.size(),1)
	TEST_REAL_EQUAL(cons_map[0].getMZ(), 944.96)

	cons_file.getOptions() = PeakFileOptions();
	cons_file.getOptions().setIntensityRange(makeRange(15000,24000));
	cons_file.load("data/ConsensusXMLFile3.xml",cons_map);
	TEST_EQUAL(cons_map.size(),1)
	TEST_REAL_EQUAL(cons_map[0].getIntensity(),23000.238)

RESULT

CHECK((void store(const String &filename, const ConsensusMap &map)))
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  
  ConsensusMap cons_map, cons_map2;
  ConsensusXMLFile cons_file;
  
  cons_file.load("data/ConsensusXMLFile2.xml",cons_map);  
  cons_file.store(tmp_filename,cons_map);
  cons_file.load(tmp_filename,cons_map2);  
  TEST_EQUAL(cons_map.size(),cons_map2.size())
  TEST_EQUAL(cons_map[0]==cons_map2[0],true)
  TEST_EQUAL(cons_map[1]==cons_map2[1],true)
RESULT

CHECK([EXTRA] (bool isValid(const String& filename)))
  ConsensusXMLFile cons_file;
  TEST_EQUAL(cons_file.isValid("data/ConsensusXMLFile.xml"),true);
  TEST_EQUAL(cons_file.isValid("data/ConsensusXMLFile2.xml"),true);
  TEST_EQUAL(cons_file.isValid("data/ConsensusXMLFile3.xml"),true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



