// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
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
  {
  public:
    typedef Eigen::Matrix<Value, Eigen::Dynamic, Eigen::Dynamic> EigenMatrixType;

    // Constructors
    Matrix() : data_(0, 0) {}

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

    // Resize
    void resize(size_t rows, size_t cols, Value value = Value())
    {
      EigenMatrixType newData(rows, cols);
      newData.fill(value);
      data_.swap(newData);
    }

    // Clear
    void clear()
    {
      data_.resize(0, 0);
    }

    // Rows and Columns
    size_t rows() const { return data_.rows(); }
    size_t cols() const { return data_.cols(); }

    // Row and Column Access
    std::vector<Value> row(size_t i) const
    {
      return std::vector<Value>(data_.row(i).data(), data_.row(i).data() + data_.cols());
    }

    std::vector<Value> col(size_t j) const
    {
      return std::vector<Value>(data_.col(j).data(), data_.col(j).data() + data_.rows());
    }

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

    bool operator<(const Matrix& rhs) const
    {
      for (int i = 0; i < rows(); ++i)
      {
        for (int j = 0; j < cols(); ++j)
        {
          if ((*this)(i, j) < rhs(i, j)) return true;
          if ((*this)(i, j) > rhs(i, j)) return false;
        }
      }
      return false;
    }

    size_t size() const
    {
      return data_.size();
    }

    bool empty() const
    {
      return data_.size() == 0;
    }

    // Iterator access
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
