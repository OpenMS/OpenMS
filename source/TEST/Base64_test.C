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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/Base64.h>

///////////////////////////

START_TEST(Base64, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::string;

// default ctor
Base64* ptr = 0;

CHECK((Base64()))
	ptr = new Base64;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK((~Base64()))
	delete ptr;
RESULT

CHECK((void encode(std::vector<Real>& in, ByteOrder to_byte_order, std::string& out)))
  PRECISION(0.001)

	Base64 b64;
  std::vector<Real> data;
  std::vector<Real> res;
  string dest;

  data.push_back(300.15);
  data.push_back(303.998);
  data.push_back(304.6);
	b64.encode(data, Base64::LITTLEENDIAN, dest);
  b64.decode(dest.c_str(), Base64::LITTLEENDIAN, res);
  TEST_REAL_EQUAL(res[0], 300.15)
  TEST_REAL_EQUAL(res[1], 303.998)
  TEST_REAL_EQUAL(res[2], 304.6)
  
	b64.encode(data, Base64::BIGENDIAN, dest);
  b64.decode(dest.c_str(), Base64::BIGENDIAN, res);
  TEST_REAL_EQUAL(res[0], 300.15)
  TEST_REAL_EQUAL(res[1], 303.998)
  TEST_REAL_EQUAL(res[2], 304.6)

  data = std::vector<Real>();
  data.push_back(4711.08);
  b64.encode(data, Base64::LITTLEENDIAN, dest);
	TEST_EQUAL(dest, "pDiTRQ==")
RESULT

CHECK((void decode(const std::string& in, ByteOrder from_byte_order, std::vector<Real>& out)))
  PRECISION(0.001)

	Base64 b64;
	string src;
	std::vector<Real> res;

	src = "QvAAAELIAA==";
  b64.decode(src.c_str(), Base64::BIGENDIAN, res);
	TEST_REAL_EQUAL(res[0], 120)
	TEST_REAL_EQUAL(res[1], 100)

	src = "Q+vIuEec9YBD7TgoR/HTgEPt23hHA8UA";
  b64.decode(src.c_str(), Base64::BIGENDIAN, res);
	TEST_REAL_EQUAL(res[0], 471.568)
	TEST_REAL_EQUAL(res[1], 80363)
	TEST_REAL_EQUAL(res[2], 474.439)
	TEST_REAL_EQUAL(res[3], 123815)
	TEST_REAL_EQUAL(res[4], 475.715)
	TEST_REAL_EQUAL(res[5], 33733)

	src = "JhOWQ8b/l0PMTJhD";
  b64.decode(src.c_str(), Base64::LITTLEENDIAN, res);
	TEST_REAL_EQUAL(res[0], 300.15)
	TEST_REAL_EQUAL(res[1], 303.998)
	TEST_REAL_EQUAL(res[2], 304.6)

	src = "QGYTSADLaUgAAABA";
  b64.decode(src.c_str(), Base64::LITTLEENDIAN, res);
	TEST_REAL_EQUAL(res[0], 150937)
	TEST_REAL_EQUAL(res[1], 239404)
	TEST_REAL_EQUAL(res[2], 2)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
