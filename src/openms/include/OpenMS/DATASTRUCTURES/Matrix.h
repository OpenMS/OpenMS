// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Macros.h>

  /**
    @brief A two-dimensional matrix based on Eigen matrices.

    @ingroup Datastructures
  */

#include <Eigen/Dense>

#include <algorithm>
#include <iostream>
#include <iomanip>

namespace OpenMS
{
  template <typename Value>
  class Matrix
  /**
   * @class Matrix
   * @brief A class representing a matrix using Eigen library.
   * 
   * The Matrix class provides functionality for creating, manipulating, and accessing matrices.
   * It is implemented using the Eigen library and supports various operations such as resizing, clearing,
   * accessing elements, setting values, and comparing matrices.
   */
  {
  public:
    typedef Eigen::Matrix<Value, Eigen::Dynamic, Eigen::Dynamic> EigenMatrixType;

    // Constructors
    Matrix() : data_(0, 0) 
    {      
    }

    Matrix(Size rows, Size cols, Value value = Value()) : data_(rows, cols)
    {
      data_.fill(value);
    }
    
    // Copy constructor
    Matrix(const Matrix& source) : data_(source.data_) {}

    // Assignment operator
    Matrix& operator=(const Matrix& rhs)
    {
      if (this != &rhs)
      {
        data_ = rhs.data_;
      }
      return *this;
    }

    // Accessors
    Value operator()(int i, int j) const
    {
      return data_(i, j);
    }

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
     * @brief Resizes the matrix to the specified number of rows and columns.
     * 
     * @param rows The number of rows in the resized matrix.
     * @param cols The number of columns in the resized matrix.
     * @param value The value to fill the resized matrix with (default: Value()).
     */
    void resize(size_t rows, size_t cols, Value value = Value())
    {
      EigenMatrixType newData(rows, cols);
      newData.fill(value);
      data_.swap(newData);
    }

    /**
     * @brief Clears the matrix by resizing it to 0x0.
     */
    void clear()
    {
      data_.resize(0, 0);
    }

    // Rows and Columns
    size_t rows() const 
    { 
      return data_.rows(); 
    }
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
      resize(ROWS, COLS);
      for (int i = 0; i < ROWS; ++i) 
      {
        for (int j = 0; j < COLS; ++j) 
        {
          data_(i, j) = array[i][j];
        }
      }
    }
 
    // Equality and Comparison
    bool operator==(const Matrix& rhs) const
    {
      return data_ == rhs.data_;
    }

    size_t size() const
    {
      return data_.size();
    }

    bool empty() const
    {
      return data_.size() == 0;
    }

    /**
     * Iterator to the beginning of the data array.
     * Warning: Order depends on the storage order of the Eigen matrix. Use with care.
     * @return Iterator to the beginning of the data array.
     */
    Value* begin()
    {
      return data_.data();
    }

    Value* end()
    {
      return data_.data() + data_.size();
    }

    const Value* begin() const
    {
      return data_.data();
    }

    const Value* end() const
    {
      return data_.data() + data_.size();
    }

    // Print
    friend std::ostream& operator<<(std::ostream& os, const Matrix<Value>& matrix)
    {
      for (size_t i = 0; i < matrix.rows(); ++i)
      {
        for (size_t j = 0; j < matrix.cols(); ++j)
        {
          os << std::setprecision(6) << std::setw(6) << matrix(i, j) << ' ';
        }
        os << std::endl;
      }
      return os;
    }

    const EigenMatrixType& getEigenMatrix() const
    {
      return data_;
    }
  private:
    EigenMatrixType data_;
  };
} // namespace OpenMS
