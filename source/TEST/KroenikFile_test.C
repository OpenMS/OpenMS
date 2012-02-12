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
#include <OpenMS/FORMAT/KroenikFile.h>
///////////////////////////

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

using namespace OpenMS;
using namespace std;

START_TEST(KroenikFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

KroenikFile* ptr = 0;
KroenikFile* null_ptr = 0;
START_SECTION(KroenikFile())
{
	ptr = new KroenikFile();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~KroenikFile())
{
	delete ptr;
}
END_SECTION


START_SECTION((template < typename FeatureMapType > void load(const String &filename, FeatureMapType &feature_map)))
{
  KroenikFile f;
  FeatureMap<> fm;
  f.load(OPENMS_GET_TEST_DATA_PATH("KroenikFile_test_1.krf"), fm);
  TEST_EQUAL(fm.size(),3)
  ABORT_IF(fm.size()!=3)
  TEST_EQUAL(fm[0].getRT(), 63.2)
  TEST_REAL_SIMILAR(fm[0].getMZ(), 1002.11)
  TEST_EQUAL(fm[0].getIntensity(), 999999)
  TEST_EQUAL(fm[0].getCharge(), 1)
  TEST_EQUAL(String(fm[0].getMetaValue("AveragineModifications")), String("Carbamido"))
  TEST_EQUAL(fm[1].getRT(),  62.2)
  TEST_REAL_SIMILAR(fm[1].getMZ(), 252.057	)
  TEST_EQUAL(fm[1].getIntensity(), 9999)
  TEST_EQUAL(fm[1].getCharge(), 2)
  TEST_EQUAL(String(fm[1].getMetaValue("AveragineModifications")), String("Carbamido2"))

  

  TEST_EXCEPTION(Exception::ParseError, f.load(OPENMS_GET_TEST_DATA_PATH("KroenikFile_test_2.krf"), fm));
  
  TEST_EXCEPTION(Exception::FileNotFound, f.load(OPENMS_GET_TEST_DATA_PATH("KroenikFile_test_2_doesnotexist.edta"), fm));
}
END_SECTION

START_SECTION((template < typename SpectrumType > void store(const String &filename, const SpectrumType &spectrum) const ))
{
  KroenikFile f;
  MSSpectrum<> spec;
  TEST_EXCEPTION(Exception::NotImplemented, f.store("bla", spec))
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



