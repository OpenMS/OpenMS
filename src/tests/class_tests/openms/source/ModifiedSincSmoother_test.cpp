// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/PROCESSING/SMOOTHING/ModifiedSincSmoother.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <vector>
#include <numeric>

///////////////////////////

START_TEST(ModifiedSincSmoother<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

ModifiedSincSmoother* ms_ptr = nullptr;
ModifiedSincSmoother* ms_nullPointer = nullptr;
START_SECTION((ModifiedSincSmoother(bool, int, int)))
  ms_ptr = new ModifiedSincSmoother(false, 6, 7);
  TEST_NOT_EQUAL(ms_ptr, ms_nullPointer)
END_SECTION

START_SECTION((virtual ~ModifiedSincSmoother()))
  delete ms_ptr;
END_SECTION

START_SECTION((std::vector<double> smooth(const std::vector<double>&)))
{
  int degree = 6;
  int m = 7;
  std::vector<double> data = {0, 1, -2, 3, -4, 5, -6, 7, -8, 9, 10, 6, 3, 1, 0};
  std::vector<double> expected = {
    0.15835885, 0.11657466, -0.09224721, 0.03165689, -0.05481473,
   -0.05436219, 0.51054827, -0.59067866, -1.21928695, 5.28610520,
   10.46161952, 6.82674246, 2.49236743, 1.04220381, 0.03264660
  };

  ModifiedSincSmoother smoother(false, degree, m);
  std::vector<double> output = smoother.smooth(data);

  for (size_t i = 0; i < data.size(); ++i)
  {
    TEST_REAL_SIMILAR(output[i], expected[i])
  }
}
END_SECTION

START_SECTION((smooth MS1 variant))
{
  int degree = 6;
  int m = 7;
  std::vector<double> data = {0, 1, -2, 3, -4, 5, -6, 7, -8, 9, 10, 6, 3, 1, 0};

  ModifiedSincSmoother smoother(true, degree, m); // MS1 variant
  std::vector<double> output = smoother.smooth(data);

  double sum_input = std::accumulate(data.begin(), data.end(), 0.0);
  double sum_output = std::accumulate(output.begin(), output.end(), 0.0);
  TEST_REAL_SIMILAR(sum_input, sum_output)

  TEST_EQUAL(data.size(), output.size())
}
END_SECTION

START_SECTION((short input, expect safe behavior))
{
  int degree = 4;
  int m = 3;
  std::vector<double> data = {1.0, 2.0, 3.0};

  ModifiedSincSmoother smoother(false, degree, m);
  std::vector<double> output = smoother.smooth(data);

  TEST_EQUAL(output.size(), data.size())
}
END_SECTION

START_SECTION((constant input stays constant))
{
  int degree = 6;
  int m = 7;
  std::vector<double> data(21, 5.0);

  ModifiedSincSmoother smoother(false, degree, m);
  std::vector<double> output = smoother.smooth(data);

  for (double val : output)
  {
    TEST_REAL_SIMILAR(val, 5.0)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST