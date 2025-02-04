// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/MSNumpressCoder.h>

///////////////////////////

#include <OpenMS/CONCEPT/Types.h>
#include <cmath>       /* pow */

using namespace std;


std::vector< double > setup_test_vec1()
{
  std::vector< double > in;
  in.push_back(100.0);
  in.push_back(200.0);
  in.push_back(300.00005);
  in.push_back(400.00010);
  return in;
}

std::vector< double > setup_test_vec2()
{
  // Compute a series of values which adds small values up to 1e-9 to an
  // integer value, giving raise to differences to the integer value of up to
  // 2e-12 - for example:
  //
  // 400                   10^-1
  // 401.01                10^-2
  // 402.002               10^-3
  // 403.0003              10^-4
  // 404.00004             10^-5
  // 405.000005            10^-6
  // 406.0000006           10^-7
  // 407.00000007          10^-8
  // 408.000000008         10^-9
  // 409.0000000009        10^-10
  // 411                   10^-1
  // 411.11                10^-2
  // 412.012               10^-3
  // 413.0013              10^-4
  // 414.00014             10^-5
  // 415.000015            10^-6
  // 416.0000016           10^-7
  // 417.00000017          10^-8
  // 418.000000018         10^-9
  // 419.0000000019        10^-10
  //
  // [ ... ]             
  //
  // 499                   10^-1
  // 491.91                10^-2
  // 492.092               10^-3
  // 493.0093              10^-4
  // 494.00094             10^-5
  // 495.000095            10^-6
  // 496.0000096           10^-7
  // 497.00000097          10^-8
  // 498.000000098         10^-9
  // 499.0000000099        10^-10
  //
  std::vector< double > in;
  for (int i = 0; i < 100; i++)
  {
    // compute a value 100 + i + i * exp(10, -i%10 -1)
    double val = 400 + i + i * pow(10.0, -i%10 -1);
    in.push_back(val);
    // std::cout << val << "\t\t\t10^" << -i%10 -1 <<std::endl;
  }
  return in;
}

bool check_vec2_abs(std::vector<double> vec, double eps)
{
  if (vec.size() != 100) return false;
  for (int i = 0; i < 100; i++)
  {
    double val = 400 + i + i * pow(10.0, -i%10 -1);
    //std::cout << val << "\t\t\tvs " << vec[i] << " diff " << fabs(val-vec[i])<<std::endl;
    if (fabs(val - vec[i]) > eps) return false;
  }
  return true;
}

bool check_vec2_rel(std::vector<double> vec, double eps)
{
  if (vec.size() != 100) return false;
  for (int i = 0; i < 100; i++)
  {
    double val = 400 + i + i * pow(10.0, -i%10 -1);
    double ratio = val / vec[i];
    if (ratio < 1.0) ratio = vec[i] / val;
    //std::cout << val << "\t\t\tvs " << vec[i] << " ratio " << val/vec[i]<<std::endl;
    if (fabs(ratio - 1) > eps) return false;

  }
  return true;
}

START_TEST(MSNumpressCoder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

// default ctor
MSNumpressCoder* ptr = nullptr;
MSNumpressCoder* nullPointer = nullptr;

START_SECTION((MSNumpressCoder()))
	ptr = new MSNumpressCoder;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~MSNumpressCoder()))
	delete ptr;
END_SECTION

START_SECTION(( void encodeNP(const std::vector<double> & in, String & result, bool zlib_compression, const NumpressConfig & config)))
{
  std::vector< double > in = setup_test_vec1();
  String out;
  MSNumpressCoder::NumpressConfig config;
  config.np_compression = MSNumpressCoder::PIC;

  bool zlib_compression = false;
  MSNumpressCoder().encodeNP(in, out, zlib_compression, config);

  TEST_EQUAL(out.size(), 12)
}
END_SECTION

START_SECTION(( void encodeNP(const std::vector<float> & in, String & result, bool zlib_compression, const NumpressConfig & config)))
{
  // tested using the encodeNP double function
  NOT_TESTABLE
}
END_SECTION

START_SECTION(( void decodeNP(const String & in, std::vector<double> & out, bool zlib_compression, const NumpressConfig & config) ))
{
  String in = "ZGaMXCFQkQ==";

  MSNumpressCoder::NumpressConfig config;
  config.np_compression = MSNumpressCoder::PIC;

  std::vector<double> out;

  bool zlib_compression = false;
  MSNumpressCoder().decodeNP(in, out, zlib_compression, config);

  TEST_EQUAL(out.size(), 4)

  TOLERANCE_ABSOLUTE(0.001)
  TEST_REAL_SIMILAR(out[0], 100.0)
  TEST_REAL_SIMILAR(out[0], 100.0)
  TEST_REAL_SIMILAR(out[1], 200.0)
  TEST_REAL_SIMILAR(out[2], 300.00005)
  TEST_REAL_SIMILAR(out[3], 400.00010)
}
END_SECTION

START_SECTION(([MSNumpressCoder::NumpressConfig] NumpressConfig()))
{
  MSNumpressCoder::NumpressConfig * config = new MSNumpressCoder::NumpressConfig();
  MSNumpressCoder::NumpressConfig * nullConfigPtr = nullptr;
	TEST_NOT_EQUAL(config, nullConfigPtr)
  delete config;
}
END_SECTION

///////////////////////////////////////////////////////////////////////////
// Encode / Decode a small vector
///////////////////////////////////////////////////////////////////////////

START_SECTION([EXTRA] encodeNP_LINEAR)
{
  std::vector< double > in = setup_test_vec1();
  String out;

  MSNumpressCoder::NumpressConfig config;
  config.np_compression = MSNumpressCoder::LINEAR;
  config.estimate_fixed_point = true; // critical

  bool zlib_compression = false;
  MSNumpressCoder().encodeNP(in, out, zlib_compression, config);

  TEST_EQUAL(out.size(), 28)
  TEST_EQUAL(out, "QWR64UAAAADo//8/0P//f1kSgA==")
}
END_SECTION

START_SECTION([EXTRA] decodeNP_LINEAR)
{
  String in = "QWR64UAAAADo//8/0P//f1kSgA==";

  MSNumpressCoder::NumpressConfig config;
  config.np_compression = MSNumpressCoder::LINEAR;

  std::vector<double> out;

  bool zlib_compression = false;
  MSNumpressCoder().decodeNP(in, out, zlib_compression, config);

  TEST_EQUAL(out.size(), 4)

  TOLERANCE_ABSOLUTE(0.001)
  TEST_REAL_SIMILAR(out[0], 100.0)
  TEST_REAL_SIMILAR(out[0], 100.0)
  TEST_REAL_SIMILAR(out[1], 200.0)
  TEST_REAL_SIMILAR(out[2], 300.00005)
  TEST_REAL_SIMILAR(out[3], 400.00010)
}
END_SECTION

START_SECTION([EXTRA] encodeNP_PIC)
{
  std::vector< double > in = setup_test_vec1();
  String out;

  MSNumpressCoder::NumpressConfig config;
  config.np_compression = MSNumpressCoder::PIC;

  bool zlib_compression = false;
  MSNumpressCoder().encodeNP(in, out, zlib_compression, config);

  TEST_EQUAL(out.size(), 12)
  TEST_EQUAL(out, "ZGaMXCFQkQ==")
}
END_SECTION

START_SECTION([EXTRA] decodeNP_PIC)
{
  String in = "ZGaMXCFQkQ==";

  MSNumpressCoder::NumpressConfig config;
  config.np_compression = MSNumpressCoder::PIC;

  std::vector<double> out;

  bool zlib_compression = false;
  MSNumpressCoder().decodeNP(in, out, zlib_compression, config);

  TEST_EQUAL(out.size(), 4)

  TOLERANCE_ABSOLUTE(0.001)
  TEST_REAL_SIMILAR(out[0], 100.0)
  TEST_REAL_SIMILAR(out[0], 100.0)
  TEST_REAL_SIMILAR(out[1], 200.0)
  TEST_REAL_SIMILAR(out[2], 300.00005)
  TEST_REAL_SIMILAR(out[3], 400.00010)
}
END_SECTION

START_SECTION([EXTRA] encodeNP_SLOF)
{
  std::vector< double > in = setup_test_vec1();
  String out;

  MSNumpressCoder::NumpressConfig config;
  config.np_compression = MSNumpressCoder::SLOF;
  config.estimate_fixed_point = true; // critical

  bool zlib_compression = false;
  MSNumpressCoder().encodeNP(in, out, zlib_compression, config);

  TEST_EQUAL(out.size(), 24)
  TEST_EQUAL(out, "QMVagAAAAAAZxX3ivPP8/w==")
}
END_SECTION

START_SECTION([EXTRA] decodeNP_SLOF)
{
  String in = "QMVagAAAAAAZxX3ivPP8/w==";

  MSNumpressCoder::NumpressConfig config;
  config.np_compression = MSNumpressCoder::SLOF;

  std::vector<double> out;

  bool zlib_compression = false;
  MSNumpressCoder().decodeNP(in, out, zlib_compression, config);

  TEST_EQUAL(out.size(), 4)

  TOLERANCE_RELATIVE(1+ 1e-4)
  TEST_REAL_SIMILAR(out[0], 100.0)
  TEST_REAL_SIMILAR(out[0], 100.0)
  TEST_REAL_SIMILAR(out[1], 200.0)
  TEST_REAL_SIMILAR(out[2], 300.00005)
  TEST_REAL_SIMILAR(out[3], 400.00010)

  //MSNumpressCoder().encodeNP(in, out, zlib_compression, config);
}
END_SECTION

///////////////////////////////////////////////////////////////////////////
// Large test
///////////////////////////////////////////////////////////////////////////

START_SECTION([EXTRA] test_large_LINEAR)
{
  std::vector< double > in = setup_test_vec2();
  String base64_string;
  std::vector<double> result;

  MSNumpressCoder::NumpressConfig config;
  config.np_compression = MSNumpressCoder::LINEAR;
  config.estimate_fixed_point = true; // critical

  bool zlib_compression = false;
  MSNumpressCoder().encodeNP(in, base64_string, zlib_compression, config);
  TEST_EQUAL(base64_string.size(), 360)
  MSNumpressCoder().decodeNP(base64_string, result, zlib_compression, config);
  TEST_EQUAL(result.size(), 100)

  TEST_EQUAL(check_vec2_abs(result, 1e-5), true)
  TEST_EQUAL(check_vec2_rel(result, 0.1e-6), true) // accurate to 0.1 ppm
}
END_SECTION

START_SECTION([EXTRA] test_large_PIC)
{
  std::vector< double > in = setup_test_vec2();
  String base64_string;
  std::vector<double> result;

  MSNumpressCoder::NumpressConfig config;
  config.np_compression = MSNumpressCoder::PIC;

  bool zlib_compression = false;
  MSNumpressCoder().encodeNP(in, base64_string, zlib_compression, config);
  TEST_EQUAL(base64_string.size(), 268)
  MSNumpressCoder().decodeNP(base64_string, result, zlib_compression, config);
  TEST_EQUAL(result.size(), 100)

  TEST_EQUAL(check_vec2_abs(result, 0.99), true)
  TEST_EQUAL(check_vec2_rel(result, 10000e-6), true) // accurate to 10 000 ppm
}
END_SECTION

START_SECTION([EXTRA] test_large_SLOF)
{
  std::vector< double > in = setup_test_vec2();
  String base64_string;
  std::vector<double> result;

  MSNumpressCoder::NumpressConfig config;
  config.np_compression = MSNumpressCoder::SLOF;
  config.estimate_fixed_point = true; // critical

  bool zlib_compression = false;
  MSNumpressCoder().encodeNP(in, base64_string, zlib_compression, config);
  TEST_EQUAL(base64_string.size(), 280)
  MSNumpressCoder().decodeNP(base64_string, result, zlib_compression, config);
  TEST_EQUAL(result.size(), 100)

  TEST_EQUAL(check_vec2_abs(result, 0.05), true)
  TEST_EQUAL(check_vec2_rel(result, 100e-6), true) // accurate to 100 ppm
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
