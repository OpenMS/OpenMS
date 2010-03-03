// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Guillaume Belz $
// $Authors: Guillaune Belz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/XMassFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include "data/XMassFile_test.h"

///////////////////////////

START_TEST(XMassFile, "$Id: XMassFile_test.C 2010-01-18 guillaume_belz $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

XMassFile* ptr = 0;
START_SECTION(XMassFile())
	ptr = new XMassFile;
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~XMassFile())
	delete ptr;
END_SECTION

START_SECTION(template<typename SpectrumType> void load(const String& filename, SpectrumType& spectrum) )
	TOLERANCE_ABSOLUTE(0.001)
	MSSpectrum<> s;
	MSSpectrum<>::ConstIterator it;
	XMassFile f;
	Size index;
	
	TEST_EXCEPTION(Exception::FileNotFound, f.load("data_Idontexist", s);)

	f.load(OPENMS_GET_TEST_DATA_PATH("XMassFile_test/fid"),s);
	
	TEST_EQUAL(s.size(), 80478)
	ABORT_IF(s.size() != 80478)

	for(it=s.begin(), index=0; it!=s.end(); it++, index++)
	{
	  TEST_REAL_SIMILAR(it->getPosition()[0], XMassFile_test_data[2*index])
	  TEST_REAL_SIMILAR(it->getIntensity(), XMassFile_test_data[2*index+1])
	}

END_SECTION

START_SECTION(template<typename SpectrumType> void store(const String& filename, const SpectrumType& spectrum) const)
  // not implemented
	TEST_EXCEPTION(Exception::NotImplemented, XMassFile().store(String(), MSSpectrum<>()))
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

