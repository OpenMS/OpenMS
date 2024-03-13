// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
  class Matrix
  {
  public:
    /**
     * @brief Eigen matrix type.
     */
    using EigenMatrixType = Eigen::Matrix<Value, Eigen::Dynamic, Eigen::Dynamic>;

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
    Matrix(Size rows, Size cols, Value value = Value()) : data_(rows, cols)
    {
      data_.fill(value);
    }
    
    /**
     * @brief Accessor to get the value at the specified position in the matrix.
     * 
     * @param i Row index.
     * @param j Column index.
     * @return reference to the value at the specified position.
     */
    const Value& operator()(int i, int j) const
    {
      return data_(i, j);
    }

    /**
     * @brief Accessor to get the value at the specified position in the matrix.
     * 
     * @param i Row index.
     * @param j Column index.
     * @return const reference to the value at the specified position.
     */
    Value& operator()(int i, int j)
    {
      return data_(i, j);
    }

    /*
      @brief get matrix entry
      Note: pyOpenMS can't easily wrap operator() so we provide additional getter / setter.
    */
    const Value& getValue(size_t const i, size_t const j) const
    {
      return data_(i, j);
    }

    /*
      @brief get matrix entry
      Note: pyOpenMS can't easily wrap operator() so we provide additional getter / setter.
    */
    Value& getValue(size_t const i, size_t const j)
    {
      return data_(i, j);
    }    

    /*
      @brief set matrix entry
      Note: pyOpenMS can't easily wrap operator() so we provide additional getter / setter.
    */
    void setValue(size_t const i, size_t const j, const Value& value)
    {
      data_(i, j) = value;
    }

    /**
     * @brief Number of rows
     */
    size_t rows() const 
    { 
      return data_.rows(); 
    }

    /**
     * @brief Number of columns
     */
    size_t cols() const 
    { 
      return data_.cols(); 
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
    template <typename T, int ROWS, int COLS>
    void setMatrix(T const (&array)[ROWS][COLS]) 
    {
      data_.resize(ROWS, COLS);
      for (int i = 0; i < ROWS; ++i) 
      {
        for (int j = 0; j < COLS; ++j) 
        {
          data_(i, j) = array[i][j];
        }
      }
    }
 
    /**
     * @brief Equality operator. Compares two matrices for equality.
     * 
     * @param rhs The matrix to be compared.
     * @return True if matrices are equal, false otherwise.
     */
    bool operator==(const Matrix& rhs) const
    {
      return data_ == rhs.data_;
    }

    /**
     * @brief Get the total number of elements in the matrix. Useful for checking if the matrix is empty
     * or iterating over raw data.
     * 
     * @return The total number of elements.
     */
    size_t size() const
    {
      return data_.size();
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
      for (size_t i = 0; i < matrix.rows(); ++i)
      {
        for (size_t j = 0; j < matrix.cols(); ++j)
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
      return data_;
    }

    /**
     * @brief Get a reference to the underlying Eigen matrix.
     * 
     * @return reference to the Eigen matrix.
     */
    EigenMatrixType& getEigenMatrix()
    {
      return data_;
    }
  private:
    EigenMatrixType data_; ///< Eigen matrix storing the actual data.
  };
} // namespace OpenMS
