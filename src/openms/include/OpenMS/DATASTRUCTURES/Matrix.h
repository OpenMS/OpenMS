// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Macros.h>

#include <Eigen/Dense>

#include <algorithm>
#include <iostream>
#include <iomanip>

namespace OpenMS
{
    /**
   * @brief A class representing a thin wrapper around an Eigen matrix.
   * 
   * The Matrix class provides functionality for creating, manipulating, and accessing matrices.
   * It is implemented using the Eigen library and supports various operations such as resizing, clearing,
   * accessing elements, setting values, and comparing matrices.
   * 
   * @ingroup Datastructures
   */
  template <typename Value>
  class Matrix : public Eigen::Matrix<Value, Eigen::Dynamic, Eigen::Dynamic>
  {
  public:
    /**
     * @brief Eigen matrix type.
     */
    using EigenMatrixType = Eigen::Matrix<Value, Eigen::Dynamic, Eigen::Dynamic>;
    using EigenMatrixType::fill;
    using EigenMatrixType::innerStride;
    using EigenMatrixType::outerStride;

    // Default constructor. Creates the "null" matrix.
    Matrix() = default;

    // Destructor
    ~Matrix() = default;

    // Copy constructor
    Matrix(const Matrix& other) = default;

    // Copy assignment operator
    Matrix& operator=(const Matrix& other) = default;

    // Move constructor
    Matrix(Matrix&& other) noexcept = default;

    // Move assignment operator
    Matrix& operator=(Matrix&& other) noexcept = default;

    /**
     * @brief Constructor to create a matrix with specified dimensions and fill value.
     * 
     * @param rows Number of rows in the matrix.
     * @param cols Number of columns in the matrix.
     * @param value Initial value to fill the matrix.
     */
    Matrix(Size rows, Size cols, Value value = Value()) : EigenMatrixType(rows, cols)
    {
      fill(value);
    }

    /*
      @brief get matrix entry
      Note: pyOpenMS can't easily wrap operator() so we provide additional getter / setter.
    */
    const Value& getValue(size_t const i, size_t const j) const
    {
      return *this(i, j);
    }

    /*
      @brief get matrix entry
      Note: pyOpenMS can't easily wrap operator() so we provide additional getter / setter.
    */
    Value& getValue(size_t const i, size_t const j)
    {
      return this->operator()(i, j);
    }    

    /*
      @brief set matrix entry
      Note: pyOpenMS can't easily wrap operator() so we provide additional getter / setter.
    */
    void setValue(size_t const i, size_t const j, const Value& value)
    {
      this->operator()(i, j) = value;
    }

    // apparently needed for cython
    void resize(size_t rows, size_t cols)
    {
      EigenMatrixType::resize(rows, cols);
    }

    int innerStride() const
    {
      return EigenMatrixType::innerStride();
    }

    int outerStride() const
    {
      return EigenMatrixType::outerStride();
    }

    bool rowMajor() const
    {
      return EigenMatrixType::IsRowMajor;
    }

    /**
     * @brief Sets the matrix values using a 2D array.
     * 
     * This function resizes the matrix to the specified number of rows and columns,
     * and then assigns the values from the 2D array to the corresponding elements
     * in the matrix.
     * 
     * @tparam T The type of the matrix elements.
     * @tparam ROWS The number of rows in the matrix.
     * @tparam COLS The number of columns in the matrix.
     * @param array The 2D array containing the values to be assigned to the matrix.
     */
    template <typename T, long int ROWS, long int COLS>
    void setMatrix(T const (&array)[ROWS][COLS]) 
    {
      resize(ROWS, COLS);
      for (int i = 0; i < ROWS; ++i) 
      {
        for (int j = 0; j < COLS; ++j) 
        {
          this->operator()(i, j) = array[i][j];
        }
      }
    }
 
    /**
     * @brief Equality operator. Compares two matrices for equality.
     * 
     * @param rhs The matrix to be compared.
     * @return True if matrices are equal, false otherwise.
     */
    bool operator==(const Matrix& rhs) const { 
      return EigenMatrixType::operator==(rhs);
    }

    /**
     * @brief Friend function to output the matrix to an output stream.
     * 
     * @param os Output stream.
     * @param matrix Matrix to be output.
     * @return Reference to the output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const Matrix<Value>& matrix)
    {
      for (long int i = 0; i < matrix.rows(); ++i)
      {
        for (long int j = 0; j < matrix.cols(); ++j)
        {
          os << std::setprecision(6) << std::setw(6) << matrix(i, j) << ' ';
        }
        os << '\n';
      }
      return os;
    }

    /**
     * @brief Get a const reference to the underlying Eigen matrix.
     * 
     * @return Const reference to the Eigen matrix.
     */
    const EigenMatrixType& getEigenMatrix() const
    {
      return *this;
    }

    /**
     * @brief Get a reference to the underlying Eigen matrix.
     * 
     * @return reference to the Eigen matrix.
     */
    EigenMatrixType& getEigenMatrix()
    {
      return *this;
    }
  };
} // namespace OpenMS
