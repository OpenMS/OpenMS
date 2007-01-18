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

CHECK((Size getOutputBufferSize()))
	Base64 b64;
  string src = "SGFsbG8gV29ybA==";  //"Hallo Worl"
	b64.decode(src.c_str(), src.size());
  TEST_EQUAL(b64.getOutputBufferSize(), 12)

  src = "SGFsbG8gV29ybGQ=";  //"Hallo World"
	b64.decode(src.c_str(), src.size());
  TEST_EQUAL(b64.getOutputBufferSize(), 12)

	src = "SGFsbG8gV29ybGQh";  //"Hallo World!"
	b64.decode(src.c_str(), src.size());
  TEST_EQUAL(b64.getOutputBufferSize(), 12)

	src = "SGFsbG8gV29ybGQhPw==";  //"Hallo World!?"
	b64.decode(src.c_str(), src.size());
  TEST_EQUAL(b64.getOutputBufferSize(), 15)

	src = "VGhpcyBpcyBvbmUgdGVzdA==";  //"This is one test"
	b64.decode(src.c_str(), src.size());
  TEST_EQUAL(b64.getOutputBufferSize(), 18)

	src = "VGhlIEVuZA==";  //"The End"
  b64.decode(src.c_str(), src.size());
  TEST_EQUAL(b64.getOutputBufferSize(), 18)
RESULT

CHECK((void setOutputBufferSize(Size s)))
	Base64 b64;
  b64.setOutputBufferSize(22);
  string src = "SGFsbG8gV29ybGQ=";  //"Hallo World"
	b64.decode(src.c_str(), src.size());
  TEST_EQUAL(b64.getOutputBufferSize(), 24)
RESULT

CHECK((char* decode(const char* src, Size size)))
	Base64 b64;
	string src = "SGFsbG8gV29ybGQ=";
	string dest = b64.decode(src.c_str(), src.size());
  TEST_EQUAL(dest, string("Hallo World"))

	src = "SGFsbG8gV29ybGQh";
	dest = b64.decode(src.c_str() ,src.size());
  TEST_EQUAL(dest, "Hallo World!")
	src = "VGhpcyBpcyBvbmUgdGVzdA==";
	dest = b64.decode(src.c_str() ,src.size());
  TEST_EQUAL(dest, "This is one test")
RESULT

CHECK((float* decodeFloat(const char* src, Size size)))
  PRECISION(0.001)
	Base64 b64;
	string src = "JhOWQ8b/l0PMTJhD";
  float* res = b64.decodeFloat(src.c_str() ,src.size());
	TEST_REAL_EQUAL(res[0], 300.15)
	TEST_REAL_EQUAL(res[1], 303.998)
	TEST_REAL_EQUAL(res[2], 304.6)

	src = "QGYTSADLaUgAAABA";
  res = b64.decodeFloat(src.c_str(), src.size());
	TEST_REAL_EQUAL(res[0], 150937)
	TEST_REAL_EQUAL(res[1], 239404)
	TEST_REAL_EQUAL(res[2], 2)
RESULT


CHECK((float* decodeFloatCorrected(const char* src, Size size)))
  PRECISION(0.001)
	Base64 b64;
	string src = "Q+vIuEec9YBD7TgoR/HTgEPt23hHA8UA";
  float* res = b64.decodeFloatCorrected(src.c_str() ,src.size());
	TEST_REAL_EQUAL(res[0], 471.568)
	TEST_REAL_EQUAL(res[1], 80363)
	TEST_REAL_EQUAL(res[2], 474.439)
	TEST_REAL_EQUAL(res[3], 123815)
	TEST_REAL_EQUAL(res[4], 475.715)
	TEST_REAL_EQUAL(res[5], 33733)
	
	src = "QvAAAELIAA==";
  res = b64.decodeFloatCorrected(src.c_str() ,src.size());
	TEST_REAL_EQUAL(res[0], 120)
	TEST_REAL_EQUAL(res[1], 100)
RESULT

CHECK((char* encode(const char* src, Size size)))
	Base64 b64;
	string src = "Hallo World";
	string dest = b64.encode(src.c_str(), src.size());
  TEST_EQUAL(dest, string("SGFsbG8gV29ybGQ="))

	src = "Hallo World!";
	dest = b64.encode(src.c_str() ,src.size());
  TEST_EQUAL(dest, "SGFsbG8gV29ybGQh")

	src = "This is one test";
	dest = b64.encode(src.c_str() ,src.size());
  TEST_EQUAL(dest, "VGhpcyBpcyBvbmUgdGVzdA==")

	src = "All around";
  dest = b64.encode(src.c_str() ,src.size());
  dest = b64.decode(dest.c_str(), dest.size());
  TEST_EQUAL(dest, src)

	src = "";
	dest = b64.encode(src.c_str() ,src.size());
	dest = b64.encode(dest.c_str(), dest.size());
	TEST_EQUAL(dest, src)
RESULT

CHECK((float* getFloatBuffer(Size size)))
	Base64 b64;
  float* data = b64.getFloatBuffer(3);
	TEST_NOT_EQUAL(data,0);
RESULT

CHECK((char* encodeFloat()))
  PRECISION(0.001)

	Base64 b64;
  float* data = b64.getFloatBuffer(3);
  data[0] = 300.15;
  data[1] = 303.998;
  data[2] = 304.6;
  string print1 = (char*) data;

  string dest = b64.encodeFloat();
  string print2 = b64.decode(dest.c_str(), dest.size());
  TEST_EQUAL(print1, print2)
  float* res = b64.decodeFloat(dest.c_str() ,dest.size());
	TEST_REAL_EQUAL(res[0], 300.15)
	TEST_REAL_EQUAL(res[1], 303.998)
	TEST_REAL_EQUAL(res[2], 304.6)

  data = b64.getFloatBuffer(1);
  data[0] = 4711.08;
  dest = b64.encodeFloat();
  res = b64.decodeFloat(dest.c_str() ,dest.size());
	TEST_EQUAL(dest, "pDiTRQ==")
	TEST_REAL_EQUAL(res[0], 4711.08)
RESULT

CHECK((char* encodeFloatCorrected()))
  PRECISION(0.001)

	Base64 b64;
  float* data = b64.getFloatBuffer(3);
  data[0] = 300.15;
  data[1] = 303.998;
  data[2] = 304.6;

  string dest = b64.encodeFloatCorrected();
  float* res = b64.decodeFloatCorrected(dest.c_str() ,dest.size());
	TEST_REAL_EQUAL(res[0], 300.15)
	TEST_REAL_EQUAL(res[1], 303.998)
	TEST_REAL_EQUAL(res[2], 304.6)
RESULT

CHECK((double* getDoubleBuffer(Size size)))
	Base64 b64;
  double* data = b64.getDoubleBuffer(3);
	TEST_NOT_EQUAL(data,0);
RESULT

CHECK((char* encodeDouble()))
  PRECISION(0.001)

	Base64 b64;
  double* data = b64.getDoubleBuffer(3);
  data[0] = 300.15;
  data[1] = 303.998;
  data[2] = 304.6;
  string print1 = (char*) data;

  string dest = b64.encodeDouble();
  string print2 = b64.decode(dest.c_str(), dest.size());
  TEST_EQUAL(print1, print2)
  double* res = b64.decodeDouble(dest.c_str() ,dest.size());
	TEST_REAL_EQUAL(res[0], 300.15)
	TEST_REAL_EQUAL(res[1], 303.998)
	TEST_REAL_EQUAL(res[2], 304.6)

  data = b64.getDoubleBuffer(1);
  data[0] = 4711.08;
  dest = b64.encodeDouble();
  res = b64.decodeDouble(dest.c_str() ,dest.size());
	TEST_REAL_EQUAL(res[0], 4711.08)
RESULT

CHECK((char* encodeDoubleCorrected()))
  PRECISION(0.001)

	Base64 b64;
  double* data = b64.getDoubleBuffer(3);
  data[0] = 300.15;
  data[1] = 303.998;
  data[2] = 304.6;

  string dest = b64.encodeDoubleCorrected();
  double* res = b64.decodeDoubleCorrected(dest.c_str() ,dest.size());
	TEST_REAL_EQUAL(res[0], 300.15)
	TEST_REAL_EQUAL(res[1], 303.998)
	TEST_REAL_EQUAL(res[2], 304.6)

  data = b64.getDoubleBuffer(1);
  data[0] = 4711.08;
  dest = b64.encodeDoubleCorrected();
  res = b64.decodeDoubleCorrected(dest.c_str() ,dest.size());
	TEST_REAL_EQUAL(res[0], 4711.08)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
