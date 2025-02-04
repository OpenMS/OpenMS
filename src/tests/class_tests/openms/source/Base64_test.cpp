// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

using namespace std;

START_TEST(Base64, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;


// default ctor
Base64* ptr = nullptr;
Base64* nullPointer = nullptr;

START_SECTION((Base64()))
  ptr = new Base64;
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~Base64()))
  delete ptr;
END_SECTION

/*

Python 

# Little Endian floats
>>> import base64
>>> import struct
>>> mynr = base64.standard_b64decode("pDiTRQ==")
>>> [struct.unpack('<f', mynr[i:i+4]) for i in range(0, len(mynr), 4) ]
[(4711.080078125,)]
 
# Big Endian doubles
>>> import base64
>>> import struct
>>> mynr = base64.standard_b64decode("QHLCZmZmZmZAcv/3ztkWh0BzCZmZmZma")
>>> [struct.unpack('>d', mynr[i:i+8]) for i in range(0, len(mynr), 8) ]
[(300.15,), (303.998,), (304.6,)]

*/

START_SECTION((template < typename FromType > void encode(std::vector< FromType > &in, ByteOrder to_byte_order, String &out, bool zlib_compression=false)))
  TOLERANCE_ABSOLUTE(0.001)
{
  Base64 b64;
  std::vector<float> data;
  std::vector<float> res;
  String dest;

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
  
  data = std::vector<float>();
  data.push_back(4711.08f);
  b64.encode(data, Base64::BYTEORDER_LITTLEENDIAN, dest);
  TEST_EQUAL(dest, "pDiTRQ==")

  // testing the encoding of double vectors
  std::vector<double> data_double;
  std::vector<double> res_double;
  data_double.push_back(300.15);
  data_double.push_back(303.998);
  data_double.push_back(304.6);
  b64.encode(data_double, Base64::BYTEORDER_BIGENDIAN, dest);
  TEST_EQUAL(dest, "QHLCZmZmZmZAcv/3ztkWh0BzCZmZmZma");
  b64.decode(dest,Base64::BYTEORDER_BIGENDIAN,res_double);
}
END_SECTION

START_SECTION((template < typename ToType > void decode(const String &in, ByteOrder from_byte_order, std::vector< ToType > &out, bool zlib_compression=false)))
  TOLERANCE_ABSOLUTE(0.001)
{
  Base64 b64;
  String src;
  std::vector<float> res;
  std::vector<double> res_double;

  b64.decode(src, Base64::BYTEORDER_BIGENDIAN, res);
  TEST_EQUAL(res.size(), 0)

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

  // test some corrupted strings
  src = "==";
  b64.decode(src, Base64::BYTEORDER_BIGENDIAN, res);
  TEST_EQUAL(res.size(), 0)

  src = "Q==";
  b64.decode(src, Base64::BYTEORDER_BIGENDIAN, res);
  TEST_EQUAL(res.size(), 0)

  src = "====";
  b64.decode(src, Base64::BYTEORDER_BIGENDIAN, res);
  TEST_EQUAL(res.size(), 0)

  // corrupted data
  src = "whoPutMeHere:somecrazyperson,obviously!WhatifIcontaininvalidcharacterslikethese";
  TEST_EXCEPTION(Exception::ConversionError, b64.decode(src, Base64::BYTEORDER_BIGENDIAN, res) );

  // TODO : some error checking and handling
  // currently there is no "safe" Base64 decoding that checks that all
  // characters are actually valid and the string is actually encoding to
  // floats.
  // 
  // src = "Q A..A=="; // spaces and dots are not allowed
  // b64.decode(src, Base64::BYTEORDER_BIGENDIAN, res);
  // TEST_EQUAL(res.size(), 0)
}
END_SECTION

START_SECTION([EXTRA] zlib functionality)
{
  TOLERANCE_ABSOLUTE(0.001)
  Base64 b64;
  String str,src;
  std::vector<float> data, res;
  std::vector<double> data_double, res_double;
  
  // double data - big endian
  data_double.push_back(300.15);
  data_double.push_back(15.124);
  data_double.push_back(304.2);
  b64.encode(data_double,Base64::BYTEORDER_BIGENDIAN,str,true);
  b64.decode(str,Base64::BYTEORDER_BIGENDIAN, res_double,true);
  TEST_REAL_SIMILAR(res_double[0],300.15);
  TEST_REAL_SIMILAR(res_double[1],15.124);
  TEST_REAL_SIMILAR(res_double[2],304.2);
  
  data.clear();
  data.push_back(120.0f);
  data.push_back(100.0f);
  b64.encode(data,Base64::BYTEORDER_BIGENDIAN,str,true);
  b64.decode(str,Base64::BYTEORDER_BIGENDIAN,res,true);

  TEST_REAL_SIMILAR(res[0], 120)
  TEST_REAL_SIMILAR(res[1], 100)  
  // float data -big endian
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
  
  // double data - little endian
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
}
END_SECTION

START_SECTION(( void encodeStrings(const std::vector<String> & in, String & out, bool zlib_compression = false, bool append_zero_byte = true)))
{
  Base64 b64;
  String src,str;
  
  //without zlib compression
  src="ZGFzAGlzdABlaW4AdGVzdAAxMjM0";
  vector<String> strings;
  b64.decodeStrings(src,strings,false);
  TEST_EQUAL(strings.size() == 5,true   )
  TEST_EQUAL(strings[0],"das")
  TEST_EQUAL(strings[1],"ist")
  TEST_EQUAL(strings[2],"ein")
  TEST_EQUAL(strings[3],"test")
  TEST_EQUAL(strings[4],"1234")

  //same as above but this time the whole String is null-terminated as well
  src="ZGFzAGlzdABlaW4AdGVzdAAxMjM0AA==";
  b64.decodeStrings(src,strings,false);
  TEST_EQUAL(strings.size() == 5,true   )
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
  TEST_EQUAL(strings.size() == 5,true )
  TEST_EQUAL(strings[0],"das")
  TEST_EQUAL(strings[1],"ist")
  TEST_EQUAL(strings[2],"ein")
  TEST_EQUAL(strings[3],"test")
  TEST_EQUAL(strings[4],"1234")

  // test some corrupted strings
  src = "==";
  b64.decodeStrings(src, strings, false);
  TEST_EQUAL(strings.size(), 0)

  src = "Q==";
  b64.decodeStrings(src, strings, false);
  TEST_EQUAL(strings.size(), 0)

  src = "====";
  b64.decodeStrings(src, strings, false);
  TEST_EQUAL(strings.size(), 0)

  src = "Q A..A=="; // spaces and dots are not allowed
  b64.decodeStrings(src, strings, false);
  // TODO : some error checking and handling
  // TEST_EQUAL(strings.size(), 0)
}
END_SECTION
  
START_SECTION((void decodeStrings(const String& in, std::vector<String>& out, bool zlib_compression = false)))
  //this functionality is tested in the encodeString test
  NOT_TESTABLE
END_SECTION

START_SECTION((void decodeSingleString(const String & in, QByteArray & base64_uncompressed, bool zlib_compression)))
  //this functionality is tested in the decodeStrings test
  NOT_TESTABLE
END_SECTION

START_SECTION((template < typename ToType > void decodeIntegers(const String &in, ByteOrder from_byte_order, std::vector< ToType > &out, bool zlib_compression=false)))
{
  Base64 b64;
  String src,str;
  vector<Int32> res;
  vector<Int64> double_res;
  //with zlib compression
  src="eJwNw4c2QgEAANAniezMIrKyUrKyMooIIdki4/8/wr3n3CAIgjZDthu2w4iddhm12x577bPfAQeNOeSwI4465rhxE044adIpp00546xzzrtg2kWXXHbFVTOumTXnunk33HTLbXcsuOue+x54aNEjjz3x1JJlzzy34oWXVr3y2htr3nrnvXUfbPjok8+++Oqb737Y9NMvW377469//gPgoxL0";

  b64.decodeIntegers(src, Base64::BYTEORDER_LITTLEENDIAN,res,true);
  
  for(Size i = 0 ; i < res.size();++i)
  {
    TEST_EQUAL(res[i], i)
  }
  
  src="eJwtxdciAgAAAMDMZBWyiUrZLdlkZJRC9l79/0f04O7lAoF/bW53hzvd5W4H3eOQe93nfg940GFHPORhjzjqUY953BOe9JSnPeNZxzznecedcNILTjntRS952Ste9ZrXnXHWOedd8IaL3vSWt73jXe953wc+dMlHPvaJT132mc994UtXXPWVa6772je+dcN3vveDH/3kZ7/41W9+94c//eVv//jXf266BcFVEvQ=";
  b64.decodeIntegers(src,Base64::BYTEORDER_LITTLEENDIAN,double_res,true);
  
  for(Size i = 0 ; i < double_res.size();++i)
  {
      TEST_EQUAL(double_res[i], i)
  }
  
  src="eJxjZGBgYAJiZiAGAAA0AAc=";
  b64.decodeIntegers(src,Base64::BYTEORDER_BIGENDIAN,res,true);
  TEST_EQUAL(res[0],16777216)
  TEST_EQUAL(res[1],33554432)
  TEST_EQUAL(res[2],50331648)
  
  //without zlib compression 32bit
  src = "AAAAAQAAAAUAAAAGAAAABwAAAAgAAAAJAAACCg==";
  
  b64.decodeIntegers(src, Base64::BYTEORDER_BIGENDIAN,res,false);
  
  TEST_EQUAL(res[0],1)
  TEST_EQUAL(res[1],5)
  TEST_EQUAL(res[2],6)
  TEST_EQUAL(res[3],7)
  TEST_EQUAL(res[4],8)
  TEST_EQUAL(res[5],9)
  TEST_EQUAL(res[6],522)
  //64bit
  src = "AAAAAAAAAAUAAAAAAAAAAwAAAAAAAAAJ";  
  b64.decodeIntegers(src, Base64::BYTEORDER_BIGENDIAN,double_res,false);  
  TEST_EQUAL(double_res[0],5)
  TEST_EQUAL(double_res[1],3)
  TEST_EQUAL(double_res[2],9)  

  //64bit
  src = "BQAAAAAAAAADAAAAAAAAAAkAAAAAAAAA";  
  b64.decodeIntegers(src, Base64::BYTEORDER_LITTLEENDIAN,double_res,false);  
  TEST_EQUAL(double_res[0],5)
  TEST_EQUAL(double_res[1],3)
  TEST_EQUAL(double_res[2],9)  
  //32bit
  src ="AQAAAAUAAAAGAAAABwAAAAgAAAAJAAAACgIAAA==";
  b64.decodeIntegers(src, Base64::BYTEORDER_LITTLEENDIAN,res,false);
  
  TEST_EQUAL(res[0],1)
  TEST_EQUAL(res[1],5)
  TEST_EQUAL(res[2],6)
  TEST_EQUAL(res[3],7)
  TEST_EQUAL(res[4],8)
  TEST_EQUAL(res[5],9)
  TEST_EQUAL(res[6],522)

  // test some corrupted strings
  src = "==";
  b64.decodeIntegers(src, Base64::BYTEORDER_BIGENDIAN,res,false);
  TEST_EQUAL(res.size(), 0)

  src = "Q==";
  b64.decodeIntegers(src, Base64::BYTEORDER_BIGENDIAN,res,false);
  TEST_EQUAL(res.size(), 0)

  src = "====";
  b64.decodeIntegers(src, Base64::BYTEORDER_BIGENDIAN,res,false);
  TEST_EQUAL(res.size(), 0)

  // src = "Q A..A=="; // spaces and dots are not allowed
  // b64.decodeIntegers(src, Base64::BYTEORDER_BIGENDIAN,res,false);
  // TODO : some error checking and handling
  // TEST_EQUAL(res.size(), 0)
}
END_SECTION

START_SECTION((template <typename FromType> void encodeIntegers(std::vector<FromType>& in, ByteOrder to_byte_order, String& out, bool zlib_compression=false)))
{
  Base64 b64;
  String tmp;
  
  //64 bit tests
  vector<Int64> vec64, vec64_in, vec64_out;
  vec64.push_back(0);
  vec64.push_back(1);
  vec64.push_back(2);
  vec64.push_back(3);
  vec64.push_back(4);
  vec64.push_back(5);
  
  //test with little endian and without compression
  tmp="";
  vec64_in = vec64;
  vec64_out.clear();
  b64.encodeIntegers(vec64_in, Base64::BYTEORDER_LITTLEENDIAN, tmp, false);
  b64.decodeIntegers(tmp, Base64::BYTEORDER_LITTLEENDIAN, vec64_out, false);
  TEST_EQUAL(vec64.size(),vec64_out.size())
  for (Size i=0; i<vec64.size(); ++i)
  {
    TEST_EQUAL(vec64[i],vec64_out[i])
  }

  //test with big endian and compression
  vec64.push_back(999999);
  tmp = "";
  vec64_in = vec64;
  vec64_out.clear();
  b64.encodeIntegers(vec64_in, Base64::BYTEORDER_BIGENDIAN, tmp, true);
  b64.decodeIntegers(tmp, Base64::BYTEORDER_BIGENDIAN, vec64_out, true);
  TEST_EQUAL(vec64.size(),vec64_out.size())
  for (Size i=0; i<vec64.size(); ++i)
  {
    TEST_EQUAL(vec64[i],vec64_out[i])
  }

  //32 bit tests  
  vector<Int32> vec32, vec32_in, vec32_out;
  vec32.push_back(0);
  vec32.push_back(5);
  vec32.push_back(10);
  vec32.push_back(15);
  vec32.push_back(20);
  vec32.push_back(25);

  //test with little endian and without compression
  tmp = "";
  vec32_in = vec32;
  vec32_out.clear();
  b64.encodeIntegers(vec32_in, Base64::BYTEORDER_LITTLEENDIAN, tmp, false);
  b64.decodeIntegers(tmp, Base64::BYTEORDER_LITTLEENDIAN, vec32_out, false);
  TEST_EQUAL(vec32.size(),vec32_out.size())
  for (Size i=0; i<vec32.size(); ++i)
  {
    TEST_EQUAL(vec32[i],vec32_out[i])
  }

  //test with big endian and compression
  vec32.push_back(999999);
  tmp = "";
  vec32_in = vec32;
  vec32_out.clear();
  b64.encodeIntegers(vec32_in, Base64::BYTEORDER_BIGENDIAN, tmp, true);
  b64.decodeIntegers(tmp, Base64::BYTEORDER_BIGENDIAN, vec32_out, true);
  TEST_EQUAL(vec32.size(),vec32_out.size())
  for (Size i=0; i<vec32.size(); ++i)
  {
    TEST_EQUAL(vec32[i],vec32_out[i])
  }
}
END_SECTION

ptr = new Base64;

START_SECTION(inline UInt32 endianize32(const UInt32& n))
  TEST_EQUAL(0, endianize32(0))  // swapping 0 should do nothing
  TEST_EQUAL(std::numeric_limits<UInt32>::max(), endianize32(std::numeric_limits<UInt32>::max()))  // swapping MAX should do nothing
  TEST_EQUAL(0x000000FF, endianize32(0xFF000000)) 
  TEST_EQUAL(0x0000FF00, endianize32(0x00FF0000)) 
  TEST_EQUAL(0x00FF0000, endianize32(0x0000FF00)) 
  TEST_EQUAL(0xFF000000, endianize32(0x000000FF))
  // random value should stay the same upon double call
  UInt32 r = (UInt32)UniqueIdGenerator::getUniqueId();
  TEST_EQUAL(r, endianize32(endianize32(r)))
END_SECTION

START_SECTION(inline UInt64 endianize64(const UInt64& n))
  TEST_EQUAL(0, endianize64(0))  // swapping 0 should do nothing
  TEST_EQUAL(std::numeric_limits<UInt64>::max(), endianize64(std::numeric_limits<UInt64>::max()))  // swapping MAX should do nothing
  TEST_EQUAL(0x00000000000000FF, endianize64(0xFF00000000000000)) 
  TEST_EQUAL(0x000000000000FF00, endianize64(0x00FF000000000000)) 
  TEST_EQUAL(0x0000000000FF0000, endianize64(0x0000FF0000000000)) 
  TEST_EQUAL(0x00000000FF000000, endianize64(0x000000FF00000000)) 
  TEST_EQUAL(0x000000FF00000000, endianize64(0x00000000FF000000)) 
  TEST_EQUAL(0x0000FF0000000000, endianize64(0x0000000000FF0000)) 
  TEST_EQUAL(0x00FF000000000000, endianize64(0x000000000000FF00)) 
  TEST_EQUAL(0xFF00000000000000, endianize64(0x00000000000000FF))
  // random value should stay the same upon double call
  UInt64 r = UniqueIdGenerator::getUniqueId();
  TEST_EQUAL(r, endianize64(endianize64(r)))
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
