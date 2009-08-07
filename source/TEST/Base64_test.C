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
#include <QtCore/QTime>
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
	data.push_back(120.0f);
	data.push_back(100.0f);
	b64.encode(data,Base64::BYTEORDER_BIGENDIAN,str,true);
	b64.decode(str,Base64::BYTEORDER_BIGENDIAN,res,true);

	TEST_REAL_SIMILAR(res[0], 120)
	TEST_REAL_SIMILAR(res[1], 100)	
	//real -big endian
	data.clear();
	data.push_back(471.568f);
	data.push_back(80363.0f);
	data.push_back(474.439f);
	data.push_back(123815.0f);
	data.push_back(475.715f);
	data.push_back(33733.0f);

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
	data.push_back(304.61f);
	
	b64.encode(data,Base64::BYTEORDER_BIGENDIAN,str,true);
	b64.decode(str,Base64::BYTEORDER_BIGENDIAN, res,true);

	TEST_REAL_SIMILAR(res[0], 300.151)
	TEST_REAL_SIMILAR(res[1],  303.9981)
	TEST_REAL_SIMILAR(res[2], 304.61)
	
	src = "JhOWQ8b/l0PMTJhD";
	b64.decode(src, Base64::BYTEORDER_LITTLEENDIAN, res);
	b64.encode(res,Base64::BYTEORDER_LITTLEENDIAN,str,true);
	b64.decode(str,Base64::BYTEORDER_LITTLEENDIAN,data,true);

	TEST_REAL_SIMILAR(data[0], 300.15f)
	TEST_REAL_SIMILAR(data[1], 303.998f)
	TEST_REAL_SIMILAR(data[2], 304.6f)	
	
END_SECTION

START_SECTION([EXTRA] integers and strings)
	Base64 b64;
	String src,str;
	vector<Real> res;
	vector<DoubleReal> double_res;
	//with zlib compression
	src="eJwNw4c2QgEAANAniezMIrKyUrKyMooIIdki4/8/wr3n3CAIgjZDthu2w4iddhm12x577bPfAQeNOeSwI4465rhxE044adIpp00546xzzrtg2kWXXHbFVTOumTXnunk33HTLbXcsuOue+x54aNEjjz3x1JJlzzy34oWXVr3y2htr3nrnvXUfbPjok8+++Oqb737Y9NMvW377469//gPgoxL0";

	b64.decode(src, Base64::BYTEORDER_LITTLEENDIAN,res,true,Base64::INTEGER);
	
	for(Size i = 0 ; i < res.size();++i)
	{
		TEST_REAL_SIMILAR(res[i], (Real) i)
	}
	
	src="eJwtxdciAgAAAMDMZBWyiUrZLdlkZJRC9l79/0f04O7lAoF/bW53hzvd5W4H3eOQe93nfg940GFHPORhjzjqUY953BOe9JSnPeNZxzznecedcNILTjntRS952Ste9ZrXnXHWOedd8IaL3vSWt73jXe953wc+dMlHPvaJT132mc994UtXXPWVa6772je+dcN3vveDH/3kZ7/41W9+94c//eVv//jXf266BcFVEvQ=";
	b64.decode(src,Base64::BYTEORDER_LITTLEENDIAN,double_res,true,Base64::INTEGER);
	
	for(Size i = 0 ; i < double_res.size();++i)
	{
			TEST_REAL_SIMILAR(double_res[i], (DoubleReal) i)
	}
	
	src="eJxjZGBgYAJiZiAGAAA0AAc=";
	b64.decode(src,Base64::BYTEORDER_BIGENDIAN,res,true,Base64::INTEGER);
	TEST_REAL_SIMILAR(res[0],16777215)
	TEST_REAL_SIMILAR(res[1],33554432 )
	TEST_REAL_SIMILAR(res[2],50331648 )
	
	//without zlib compression 32bit
	src = "AAAAAQAAAAUAAAAGAAAABwAAAAgAAAAJAAACCg==";
	
	b64.decode(src, Base64::BYTEORDER_BIGENDIAN,res,false,Base64::INTEGER);
	
	TEST_REAL_SIMILAR(res[0],1)
	TEST_REAL_SIMILAR(res[1],5)
	TEST_REAL_SIMILAR(res[2],6)
	TEST_REAL_SIMILAR(res[3],7)
	TEST_REAL_SIMILAR(res[4],8)
	TEST_REAL_SIMILAR(res[5],9)
	TEST_REAL_SIMILAR(res[6],522)
	//64bit
	src = "AAAAAAAAAAUAAAAAAAAAAwAAAAAAAAAJ";	
	b64.decode(src, Base64::BYTEORDER_BIGENDIAN,double_res,false,Base64::INTEGER);	
	TEST_REAL_SIMILAR(double_res[0],5)
	TEST_REAL_SIMILAR(double_res[1],3)
	TEST_REAL_SIMILAR(double_res[2],9)	

	//64bit
	src = "BQAAAAAAAAADAAAAAAAAAAkAAAAAAAAA";	
	b64.decode(src, Base64::BYTEORDER_LITTLEENDIAN,double_res,false,Base64::INTEGER);	
	TEST_REAL_SIMILAR(double_res[0],5)
	TEST_REAL_SIMILAR(double_res[1],3)
	TEST_REAL_SIMILAR(double_res[2],9)	
	//32bit
	src ="AQAAAAUAAAAGAAAABwAAAAgAAAAJAAAACgIAAA==";
	b64.decode(src, Base64::BYTEORDER_LITTLEENDIAN,res,false,Base64::INTEGER);
	
	TEST_REAL_SIMILAR(res[0],1)
	TEST_REAL_SIMILAR(res[1],5)
	TEST_REAL_SIMILAR(res[2],6)
	TEST_REAL_SIMILAR(res[3],7)
	TEST_REAL_SIMILAR(res[4],8)
	TEST_REAL_SIMILAR(res[5],9)
	TEST_REAL_SIMILAR(res[6],522)
	//encode and decode of strings
	//without zlib compression
	src="ZGFzAGlzdABlaW4AdGVzdAAxMjM0";
	vector<String> strings;
	b64.decodeStrings(src,strings,false);
	TEST_EQUAL(strings.size() == 5,true 	)
	TEST_EQUAL(strings[0],"das")
	TEST_EQUAL(strings[1],"ist")
	TEST_EQUAL(strings[2],"ein")
	TEST_EQUAL(strings[3],"test")
	TEST_EQUAL(strings[4],"1234")
	//same as above but this time the hole string is null-terminated as well
	src="ZGFzAGlzdABlaW4AdGVzdAAxMjM0AA==";
	b64.decodeStrings(src,strings,false);
	TEST_EQUAL(strings.size() == 5,true 	)
	TEST_EQUAL(strings[0],"das")
	TEST_EQUAL(strings[1],"ist")
	TEST_EQUAL(strings[2],"ein")
	TEST_EQUAL(strings[3],"test")
	TEST_EQUAL(strings[4],"1234")
	
	//zlib compressed			
	src = "eJxLSSxmyCwuYUjNzGMoSQUyDI2MTRgAUX4GTw==";
	b64.decodeStrings(src,strings,true);
	TEST_EQUAL(strings.size() == 5,true )
	TEST_EQUAL(strings[0],"das")
	TEST_EQUAL(strings[1],"ist")
	TEST_EQUAL(strings[2],"ein")
	TEST_EQUAL(strings[3],"test")
	TEST_EQUAL(strings[4],"1234")
	
	//without zlib compression
	b64.encodeStrings(strings,str,false);
	b64.decodeStrings(str,strings,false);
	TEST_EQUAL(strings[0],"das")
	TEST_EQUAL(strings[1],"ist")
	TEST_EQUAL(strings[2],"ein")
	TEST_EQUAL(strings[3],"test")
	TEST_EQUAL(strings[4],"1234")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
