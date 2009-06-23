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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/Base64.h>

///////////////////////////

#include <fstream>
#include <string>

#include <QtCore/QString>
#include <OpenMS/CONCEPT/Types.h>

using namespace std;






START_TEST(Base64, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;


// default ctor
Base64* ptr = 0;

START_SECTION((Base64()))
	ptr = new Base64;
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

// destructor
START_SECTION((virtual ~Base64()))
	delete ptr;
END_SECTION

START_SECTION((template < typename FromType > void encode(std::vector< FromType > &in, ByteOrder to_byte_order, std::string &out, bool zlib_compression=false)))
  TOLERANCE_ABSOLUTE(0.001)

	Base64 b64;
  std::vector<Real> data;
  std::vector<Real> res;
  string dest;

	b64.encode(data, Base64::BYTEORDER_LITTLEENDIAN, dest);
	TEST_EQUAL(dest, "");

  data.push_back(300.15f);
  data.push_back(303.998f);
  data.push_back(304.6f);
	b64.encode(data, Base64::BYTEORDER_LITTLEENDIAN, dest);
	TEST_EQUAL(dest, "MxOWQ77/l0PNTJhD");
	// please remember that it is possible that two different strings can
	// decode to the "same" floating point number (considering such a low
	// precision like 0.001).
	
  data = std::vector<Real>();
  data.push_back(4711.08f);
  b64.encode(data, Base64::BYTEORDER_LITTLEENDIAN, dest);
	TEST_EQUAL(dest, "pDiTRQ==")

	// testing the encoding of double vectors
  std::vector<DoubleReal> data_double;
  std::vector<DoubleReal> res_double;
  data_double.push_back(300.15);
  data_double.push_back(303.998);
  data_double.push_back(304.6);
	b64.encode(data_double, Base64::BYTEORDER_BIGENDIAN, dest);
	TEST_EQUAL(dest, "QHLCZmZmZmZAcv/3ztkWh0BzCZmZmZma");
	
	
	
END_SECTION

START_SECTION((template < typename ToType > void decode(const std::string &in, ByteOrder from_byte_order, std::vector< ToType > &out, bool zlib_compression=false)))
  TOLERANCE_ABSOLUTE(0.001)

	Base64 b64;
	string src;
	std::vector<Real> res;
	std::vector<DoubleReal> res_double;

  b64.decode(src, Base64::BYTEORDER_BIGENDIAN, res);
	TEST_EQUAL(res.size(), 0)

	src = "QvAAAELIAA==";
  b64.decode(src, Base64::BYTEORDER_BIGENDIAN, res);
	TEST_REAL_SIMILAR(res[0], 120)
	TEST_REAL_SIMILAR(res[1], 100)

	src = "Q+vIuEec9YBD7TgoR/HTgEPt23hHA8UA";
  b64.decode(src, Base64::BYTEORDER_BIGENDIAN, res);
	TEST_REAL_SIMILAR(res[0], 471.568)
	TEST_REAL_SIMILAR(res[1], 80363)
	TEST_REAL_SIMILAR(res[2], 474.439)
	TEST_REAL_SIMILAR(res[3], 123815)
	TEST_REAL_SIMILAR(res[4], 475.715)
	TEST_REAL_SIMILAR(res[5], 33733)

	src = "JhOWQ8b/l0PMTJhD";
  b64.decode(src, Base64::BYTEORDER_LITTLEENDIAN, res);
	TEST_REAL_SIMILAR(res[0], 300.15)
	TEST_REAL_SIMILAR(res[1], 303.998)
	TEST_REAL_SIMILAR(res[2], 304.6)

	src = "QGYTSADLaUgAAABA";
  b64.decode(src, Base64::BYTEORDER_LITTLEENDIAN, res);
	TEST_REAL_SIMILAR(res[0], 150937)
	TEST_REAL_SIMILAR(res[1], 239404)
	TEST_REAL_SIMILAR(res[2], 2)

	src = "QHLCZmZmZmZAcv/3ztkWh0BzCZmZmZma";
  b64.decode(src, Base64::BYTEORDER_BIGENDIAN, res_double);
	TEST_REAL_SIMILAR(res_double[0], 300.15)
	TEST_REAL_SIMILAR(res_double[1], 303.998)
	TEST_REAL_SIMILAR(res_double[2], 304.6)
END_SECTION

START_SECTION([EXTRA] zlib functionality)
	TOLERANCE_ABSOLUTE(0.001)
	Base64 b64;
	string str,src;
	std::vector<Real> data,res;
	std::vector<DoubleReal> data_double,res_double;
	
	//double_real .- big endian
	data_double.push_back(300.15);
	data_double.push_back(303.998);
	data_double.push_back(304.6);
	b64.encode(data_double,Base64::BYTEORDER_BIGENDIAN,str,true);
	b64.decode(str,Base64::BYTEORDER_BIGENDIAN, res_double,true);
	TEST_REAL_SIMILAR(res_double[0],300.15);
	TEST_REAL_SIMILAR(res_double[1],303.998);
	TEST_REAL_SIMILAR(res_double[2],304.6);
	
	data.clear();
	data.push_back(120);
	data.push_back(100);
	b64.encode(data,Base64::BYTEORDER_BIGENDIAN,str,true);
	b64.decode(str,Base64::BYTEORDER_BIGENDIAN,res,true);

	TEST_REAL_SIMILAR(res[0], 120)
	TEST_REAL_SIMILAR(res[1], 100)	
	//real -big endian
	data.clear();
	data.push_back(471.568);
	data.push_back(80363);
	data.push_back(474.439);
	data.push_back(123815);
	data.push_back(475.715);
	data.push_back(33733);

	b64.encode(data,Base64::BYTEORDER_BIGENDIAN,str,true);
	b64.decode(str,Base64::BYTEORDER_BIGENDIAN, res,true);
	
	TEST_REAL_SIMILAR(res[0], 471.568)
	TEST_REAL_SIMILAR(res[1], 80363)
	TEST_REAL_SIMILAR(res[2], 474.439)
	TEST_REAL_SIMILAR(res[3], 123815)
	TEST_REAL_SIMILAR(res[4], 475.715)
	TEST_REAL_SIMILAR(res[5], 33733)
	
	//double_real - little endian
	data.clear();
	data.push_back(300.15f);
	data.push_back(303.998f);
	data.push_back(304.61);
	
	b64.encode(data,Base64::BYTEORDER_BIGENDIAN,str,true);
	b64.decode(str,Base64::BYTEORDER_BIGENDIAN, res,true);

	TEST_REAL_SIMILAR(res[0], 300.151)
	TEST_REAL_SIMILAR(res[1],  303.9981)
	TEST_REAL_SIMILAR(res[2], 304.61)
	
	src = "JhOWQ8b/l0PMTJhD";
	b64.decode(src, Base64::BYTEORDER_LITTLEENDIAN, res);
	b64.encode(res,Base64::BYTEORDER_LITTLEENDIAN,str,true);
	b64.decode(str,Base64::BYTEORDER_LITTLEENDIAN,data,true);

	TEST_REAL_SIMILAR(data[0], 300.15)
	TEST_REAL_SIMILAR(data[1], 303.998)
	TEST_REAL_SIMILAR(data[2], 304.6)	
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
