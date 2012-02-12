// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/SpecArrayFile.h>
///////////////////////////

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

using namespace OpenMS;
using namespace std;

START_TEST(SpecArrayFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpecArrayFile* ptr = 0;
SpecArrayFile* null_ptr = 0;
START_SECTION(SpecArrayFile())
{
	ptr = new SpecArrayFile();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~SpecArrayFile())
{
	delete ptr;
}
END_SECTION


START_SECTION((template < typename FeatureMapType > void load(const String &filename, FeatureMapType &feature_map)))
{
  SpecArrayFile f;
  FeatureMap<> fm;
  f.load(OPENMS_GET_TEST_DATA_PATH("SpecArrayFile_test_1.pepList"), fm);
  TEST_EQUAL(fm.size(),2)
  ABORT_IF(fm.size()!=2)
  TEST_EQUAL(fm[0].getRT(), 60.1*60)
  TEST_REAL_SIMILAR(fm[0].getMZ(), 500.1)
  TEST_EQUAL(fm[0].getIntensity(), 4343534)
  TEST_EQUAL(fm[0].getCharge(), 5)
  TEST_EQUAL(double(fm[0].getMetaValue("s/n")), 3.2)
  TEST_EQUAL(fm[1].getRT(),  40.1*60)
  TEST_REAL_SIMILAR(fm[1].getMZ(), 700.1	)
  TEST_EQUAL(fm[1].getIntensity(), 222432)
  TEST_EQUAL(fm[1].getCharge(), 3)
  TEST_EQUAL(double(fm[1].getMetaValue("s/n")), 2.2)

  

  TEST_EXCEPTION(Exception::ParseError, f.load(OPENMS_GET_TEST_DATA_PATH("SpecArrayFile_test_2.pepList"), fm));
  
  TEST_EXCEPTION(Exception::FileNotFound, f.load(OPENMS_GET_TEST_DATA_PATH("SpecArrayFile_test_2_doesnotexist.pepList"), fm));
}
END_SECTION

START_SECTION((template < typename SpectrumType > void store(const String &filename, const SpectrumType &spectrum) const ))
{
  SpecArrayFile f;
  MSSpectrum<> spec;
  TEST_EXCEPTION(Exception::NotImplemented, f.store("bla", spec))
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



