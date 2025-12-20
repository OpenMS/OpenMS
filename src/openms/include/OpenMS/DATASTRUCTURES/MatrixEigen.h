// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// NOTE: This is an INTERNAL header file and should NOT be included in public headers.
// It provides Eigen-based utilities for internal OpenMS algorithms.
// Users of the public API should use Matrix.h directly.

#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <Eigen/Core>

namespace OpenMS
{
  /// Type alias for Eigen dynamic matrix (internal use only)
  template <typename Value>
  using EigenMatrixType = Eigen::Matrix<Value, Eigen::Dynamic, Eigen::Dynamic>;

  /// Type alias for Eigen dynamic column vector (internal use only)
  template <typename Value>
  using EigenVectorType = Eigen::Matrix<Value, Eigen::Dynamic, 1>;

  /// Type alias for Eigen dynamic row vector (internal use only)
  template <typename Value>
  using EigenRowVectorType = Eigen::Matrix<Value, 1, Eigen::Dynamic>;

  /**
   * @brief Create a mutable Eigen::Map view of an OpenMS Matrix.
   *
   * Zero-copy: the returned Map operates directly on the Matrix's buffer.
   * The Matrix must outlive the returned Map.
   *
   * @param mat The matrix to view
   * @return Eigen::Map providing full Eigen matrix operations
   *
   * @note Internal use only. This function requires Eigen headers.
   */
  template <typename Value>
  Eigen::Map<EigenMatrixType<Value>> eigenView(Matrix<Value>& mat)
  {
    return Eigen::Map<EigenMatrixType<Value>>(mat.data(), mat.rows(), mat.cols());
  }

  /**
   * @brief Create a const Eigen::Map view of an OpenMS Matrix.
   *
   * Zero-copy: the returned Map operates directly on the Matrix's buffer.
   * The Matrix must outlive the returned Map.
   *
   * @param mat The matrix to view
   * @return Const Eigen::Map for read-only Eigen operations
   */
  template <typename Value>
  Eigen::Map<const EigenMatrixType<Value>> eigenView(const Matrix<Value>& mat)
  {
    return Eigen::Map<const EigenMatrixType<Value>>(mat.data(), mat.rows(), mat.cols());
  }

  /**
   * @brief Create an Eigen::Map view of a std::vector as a column vector.
   *
   * Zero-copy: the returned Map operates directly on the vector's buffer.
   *
   * @param vec The vector to view
   * @return Eigen::Map treating the vector as a column vector
   */
  template <typename Value>
  Eigen::Map<EigenVectorType<Value>> eigenVectorView(std::vector<Value>& vec)
  {
    return Eigen::Map<EigenVectorType<Value>>(vec.data(), vec.size());
  }

  /**
   * @brief Create a const Eigen::Map view of a std::vector as a column vector.
   */
  template <typename Value>
  Eigen::Map<const EigenVectorType<Value>> eigenVectorView(const std::vector<Value>& vec)
  {
    return Eigen::Map<const EigenVectorType<Value>>(vec.data(), vec.size());
  }

  /**
   * @brief Create an Eigen::Map view of a raw pointer + size as column vector.
   *
   * @param data Pointer to the first element
   * @param size Number of elements
   * @return Eigen::Map treating the data as a column vector
   */
  template <typename Value>
  Eigen::Map<EigenVectorType<Value>> eigenVectorView(Value* data, size_t size)
  {
    return Eigen::Map<EigenVectorType<Value>>(data, size);
  }

  /**
   * @brief Create a const Eigen::Map view of a raw pointer + size as column vector.
   */
  template <typename Value>
  Eigen::Map<const EigenVectorType<Value>> eigenVectorView(const Value* data, size_t size)
  {
    return Eigen::Map<const EigenVectorType<Value>>(data, size);
  }

  /**
   * @brief Create an Eigen::Map view of a raw pointer + dimensions as a matrix.
   *
   * @param data Pointer to the first element (column-major storage)
   * @param rows Number of rows
   * @param cols Number of columns
   * @return Eigen::Map treating the data as a matrix
   */
  template <typename Value>
  Eigen::Map<EigenMatrixType<Value>> eigenMatrixView(Value* data, size_t rows, size_t cols)
  {
    return Eigen::Map<EigenMatrixType<Value>>(data, rows, cols);
  }

  /**
   * @brief Create a const Eigen::Map view of a raw pointer + dimensions as a matrix.
   */
  template <typename Value>
  Eigen::Map<const EigenMatrixType<Value>> eigenMatrixView(const Value* data, size_t rows, size_t cols)
  {
    return Eigen::Map<const EigenMatrixType<Value>>(data, rows, cols);
  }

} // namespace OpenMS
