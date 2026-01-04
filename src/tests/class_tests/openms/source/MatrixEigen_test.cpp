// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
// This tests the internal MatrixEigen utilities.
// NOTE: MatrixEigen.h is an internal header not part of the public API.
#include <OpenMS/DATASTRUCTURES/MatrixEigen.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MatrixEigen, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((eigenView - mutable Matrix))
{
  // Create a 3x4 matrix with distinct values
  Matrix<double> mat(3, 4, 0.0);
  for (Size i = 0; i < 3; ++i)
  {
    for (Size j = 0; j < 4; ++j)
    {
      mat(i, j) = static_cast<double>(i * 10 + j);
    }
  }

  // Get mutable Eigen view
  auto view = eigenView(mat);

  // Test dimensions
  TEST_EQUAL(view.rows(), 3)
  TEST_EQUAL(view.cols(), 4)

  // Test that view sees correct values
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      TEST_REAL_SIMILAR(view(i, j), static_cast<double>(i * 10 + j))
    }
  }

  // Test zero-copy: modify through view, check original
  view(1, 2) = 999.0;
  TEST_REAL_SIMILAR(mat(1, 2), 999.0)

  // Test zero-copy: modify original, check view
  mat(0, 3) = 888.0;
  TEST_REAL_SIMILAR(view(0, 3), 888.0)
}
END_SECTION

START_SECTION((eigenView - const Matrix))
{
  Matrix<double> mat(2, 3, 0.0);
  mat(0, 0) = 1.0;
  mat(0, 1) = 2.0;
  mat(0, 2) = 3.0;
  mat(1, 0) = 4.0;
  mat(1, 1) = 5.0;
  mat(1, 2) = 6.0;

  const Matrix<double>& const_mat = mat;
  auto view = eigenView(const_mat);

  // Test dimensions
  TEST_EQUAL(view.rows(), 2)
  TEST_EQUAL(view.cols(), 3)

  // Test values
  TEST_REAL_SIMILAR(view(0, 0), 1.0)
  TEST_REAL_SIMILAR(view(0, 1), 2.0)
  TEST_REAL_SIMILAR(view(0, 2), 3.0)
  TEST_REAL_SIMILAR(view(1, 0), 4.0)
  TEST_REAL_SIMILAR(view(1, 1), 5.0)
  TEST_REAL_SIMILAR(view(1, 2), 6.0)
}
END_SECTION

START_SECTION((eigenVectorView - mutable std::vector))
{
  std::vector<double> vec = {1.0, 2.0, 3.0, 4.0, 5.0};

  auto view = eigenVectorView(vec);

  // Test size
  TEST_EQUAL(view.size(), 5)

  // Test values
  for (int i = 0; i < 5; ++i)
  {
    TEST_REAL_SIMILAR(view(i), static_cast<double>(i + 1))
  }

  // Test zero-copy: modify through view
  view(2) = 99.0;
  TEST_REAL_SIMILAR(vec[2], 99.0)

  // Test zero-copy: modify original
  vec[4] = 88.0;
  TEST_REAL_SIMILAR(view(4), 88.0)
}
END_SECTION

START_SECTION((eigenVectorView - const std::vector))
{
  const std::vector<double> vec = {10.0, 20.0, 30.0};

  auto view = eigenVectorView(vec);

  TEST_EQUAL(view.size(), 3)
  TEST_REAL_SIMILAR(view(0), 10.0)
  TEST_REAL_SIMILAR(view(1), 20.0)
  TEST_REAL_SIMILAR(view(2), 30.0)
}
END_SECTION

START_SECTION((eigenVectorView - raw pointer mutable))
{
  double data[4] = {1.5, 2.5, 3.5, 4.5};

  auto view = eigenVectorView(data, 4);

  TEST_EQUAL(view.size(), 4)
  TEST_REAL_SIMILAR(view(0), 1.5)
  TEST_REAL_SIMILAR(view(3), 4.5)

  // Zero-copy test
  view(1) = 100.0;
  TEST_REAL_SIMILAR(data[1], 100.0)
}
END_SECTION

START_SECTION((eigenVectorView - raw pointer const))
{
  const double data[3] = {7.0, 8.0, 9.0};

  auto view = eigenVectorView(data, 3);

  TEST_EQUAL(view.size(), 3)
  TEST_REAL_SIMILAR(view(0), 7.0)
  TEST_REAL_SIMILAR(view(1), 8.0)
  TEST_REAL_SIMILAR(view(2), 9.0)
}
END_SECTION

START_SECTION((eigenMatrixView - raw pointer mutable))
{
  // Column-major storage: 2x3 matrix
  // Layout: [col0_row0, col0_row1, col1_row0, col1_row1, col2_row0, col2_row1]
  double data[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

  auto view = eigenMatrixView(data, 2, 3);

  TEST_EQUAL(view.rows(), 2)
  TEST_EQUAL(view.cols(), 3)

  // Column-major: data[0]=mat(0,0), data[1]=mat(1,0), data[2]=mat(0,1), etc.
  TEST_REAL_SIMILAR(view(0, 0), 1.0)
  TEST_REAL_SIMILAR(view(1, 0), 2.0)
  TEST_REAL_SIMILAR(view(0, 1), 3.0)
  TEST_REAL_SIMILAR(view(1, 1), 4.0)
  TEST_REAL_SIMILAR(view(0, 2), 5.0)
  TEST_REAL_SIMILAR(view(1, 2), 6.0)

  // Zero-copy test
  view(1, 1) = 999.0;
  TEST_REAL_SIMILAR(data[3], 999.0)  // col1_row1 is at index 3
}
END_SECTION

START_SECTION((eigenMatrixView - raw pointer const))
{
  const double data[4] = {10.0, 20.0, 30.0, 40.0};

  auto view = eigenMatrixView(data, 2, 2);

  TEST_EQUAL(view.rows(), 2)
  TEST_EQUAL(view.cols(), 2)
  TEST_REAL_SIMILAR(view(0, 0), 10.0)
  TEST_REAL_SIMILAR(view(1, 0), 20.0)
  TEST_REAL_SIMILAR(view(0, 1), 30.0)
  TEST_REAL_SIMILAR(view(1, 1), 40.0)
}
END_SECTION

START_SECTION((Eigen operations on views))
{
  // Test that Eigen operations work correctly on views
  Matrix<double> mat(3, 3, 0.0);
  mat(0, 0) = 1.0; mat(0, 1) = 0.0; mat(0, 2) = 0.0;
  mat(1, 0) = 0.0; mat(1, 1) = 2.0; mat(1, 2) = 0.0;
  mat(2, 0) = 0.0; mat(2, 1) = 0.0; mat(2, 2) = 3.0;

  auto view = eigenView(mat);

  // Test isIdentity (should be false, it's diagonal but not identity)
  TEST_EQUAL(view.isIdentity(0.0), false)

  // Test isDiagonal
  TEST_EQUAL(view.isDiagonal(0.0), true)

  // Test sum
  TEST_REAL_SIMILAR(view.sum(), 6.0)

  // Test trace
  TEST_REAL_SIMILAR(view.trace(), 6.0)

  // Test matrix-vector multiplication
  std::vector<double> vec = {1.0, 1.0, 1.0};
  auto vec_view = eigenVectorView(vec);
  Eigen::VectorXd result = view * vec_view;
  TEST_REAL_SIMILAR(result(0), 1.0)
  TEST_REAL_SIMILAR(result(1), 2.0)
  TEST_REAL_SIMILAR(result(2), 3.0)
}
END_SECTION

START_SECTION((Matrix and view consistency with column-major storage))
{
  // Verify that OpenMS Matrix and Eigen view agree on element positions
  Matrix<double> mat(3, 4, 0.0);

  // Fill via Matrix interface
  for (Size i = 0; i < 3; ++i)
  {
    for (Size j = 0; j < 4; ++j)
    {
      mat(i, j) = static_cast<double>(i * 100 + j);
    }
  }

  auto view = eigenView(mat);

  // Both should agree on all elements
  for (Size i = 0; i < 3; ++i)
  {
    for (Size j = 0; j < 4; ++j)
    {
      double expected = static_cast<double>(i * 100 + j);
      TEST_REAL_SIMILAR(mat(i, j), expected)
      TEST_REAL_SIMILAR(view(i, j), expected)
    }
  }

  // Modify via Eigen, verify via Matrix
  view(2, 3) = 12345.0;
  TEST_REAL_SIMILAR(mat(2, 3), 12345.0)

  // Modify via Matrix, verify via Eigen
  mat(0, 0) = 54321.0;
  TEST_REAL_SIMILAR(view(0, 0), 54321.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
