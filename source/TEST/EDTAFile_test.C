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
#include <OpenMS/FORMAT/EDTAFile.h>
///////////////////////////

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

using namespace OpenMS;
using namespace std;

START_TEST(EDTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EDTAFile* ptr = 0;
EDTAFile* null_ptr = 0;
START_SECTION(EDTAFile())
{
	ptr = new EDTAFile();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~EDTAFile())
{
	delete ptr;
}
END_SECTION


START_SECTION(void load(const String &filename, ConsensusMap &consensus_map))
{
  EDTAFile f;
  ConsensusMap fm;
  f.load(OPENMS_GET_TEST_DATA_PATH("EDTAFile_test_1.edta"), fm);
  TEST_EQUAL(fm.size(),2)
  ABORT_IF(fm.size()!=2)
  TEST_EQUAL(fm[0].getRT(), 321)
  TEST_EQUAL(fm[0].getMZ(), 405.233)
  TEST_EQUAL(fm[0].getIntensity(), 24543534)
  TEST_EQUAL(fm[0].getCharge(), 2)
  TEST_EQUAL(String(fm[0].getMetaValue("mymeta")), String("lala"))
  TEST_EQUAL(fm[1].getRT(), 322)
  TEST_EQUAL(fm[1].getMZ(), 406.207)
  TEST_EQUAL(fm[1].getIntensity(), 4343344)
  TEST_EQUAL(fm[1].getCharge(), 3)
  TEST_EQUAL(String(fm[1].getMetaValue("mymeta")), String("blubb"))

  
  f.load(OPENMS_GET_TEST_DATA_PATH("EDTAFile_test_3.edta"), fm);
  TEST_EQUAL(fm.size(),3)

  TEST_EXCEPTION(Exception::ParseError, f.load(OPENMS_GET_TEST_DATA_PATH("EDTAFile_test_2.edta"), fm));
  
  TEST_EXCEPTION(Exception::FileNotFound, f.load(OPENMS_GET_TEST_DATA_PATH("EDTAFile_test_3_doesnotexist.edta"), fm));
      
}
END_SECTION

START_SECTION((void store(const String& filename, const ConsensusMap& map) const))
{
  EDTAFile f;
  ConsensusMap fm;
  TEST_EXCEPTION(Exception::NotImplemented, f.store("bla", fm))
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



