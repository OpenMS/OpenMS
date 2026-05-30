// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/Types.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdexcept>

namespace OpenMS
{
  /**
   * @brief A 2D matrix class with efficient buffer access for NumPy interoperability.
   *
   * Data is stored in column-major order (Fortran style) for compatibility with
   * Eigen and most numerical libraries. Internal OpenMS code can obtain zero-copy
   * Eigen views via the eigenView() function in MatrixEigen.h (internal header).
   *
   * @tparam Value The element type (typically double or float)
   *
   * @ingroup Datastructures
   */
  template <typename Value>
  class Matrix
  {
  public:
    ///@name Type definitions
    ///@{
    using value_type = Value;
    using iterator = typename std::vector<Value>::iterator;
    using const_iterator = typename std::vector<Value>::const_iterator;
    ///@}

    ///@name Constructors and assignment
    ///@{

    /// Default constructor creates empty matrix
    Matrix() : rows_(0), cols_(0) {}

    /**
     * @brief Constructor to create a matrix with specified dimensions and fill value.
     *
     * @param rows Number of rows in the matrix.
     * @param cols Number of columns in the matrix.
     * @param value Initial value to fill the matrix.
     */
    Matrix(Size rows, Size cols, Value value = Value())
      : data_(rows * cols, value), rows_(rows), cols_(cols) {}

    /// Copy constructor
    Matrix(const Matrix&) = default;

    /// Move constructor
    Matrix(Matrix&&) noexcept = default;

    /// Copy assignment operator
    Matrix& operator=(const Matrix&) = default;

    /// Move assignment operator
    Matrix& operator=(Matrix&&) noexcept = default;

    /// Destructor
    ~Matrix() = default;

    ///@}

    ///@name Dimension accessors
    ///@{

    /// Number of rows
    Size rows() const { return rows_; }

    /// Number of columns
    Size cols() const { return cols_; }

    /// Total number of elements
    Size size() const { return data_.size(); }

    /// Check if matrix is empty
    bool empty() const { return data_.empty(); }

    ///@}

    ///@name Element access
    ///@{

    /**
     * @brief Get element at (row, col)
     *
     * Note: pyOpenMS can't easily wrap operator() so we provide additional getter.
     *
     * @param[in] row Zero-based row index.
     * @param[in] col Zero-based column index.
     * @return Const reference to the value at the given position.
     */
    const Value& getValue(size_t const row, size_t const col) const
    {
      return data_[col * rows_ + row];  // Column-major
    }

    /**
     * @brief Get mutable element at (row, col)
     *
     * Note: pyOpenMS can't easily wrap operator() so we provide additional getter.
     *
     * @param[in] row Zero-based row index.
     * @param[in] col Zero-based column index.
     * @return Mutable reference to the value at the given position.
     */
    Value& getValue(size_t const row, size_t const col)
    {
      return data_[col * rows_ + row];  // Column-major
    }

    /**
     * @brief Set element at (row, col)
     *
     * Note: pyOpenMS can't easily wrap operator() so we provide additional setter.
     *
     * @param[in] row Zero-based row index.
     * @param[in] col Zero-based column index.
     * @param[in] value Value to set at the given position.
     */
    void setValue(size_t const row, size_t const col, const Value& value)
    {
      data_[col * rows_ + row] = value;  // Column-major
    }

    /// Unchecked element access (row, col)
    Value& operator()(size_t row, size_t col)
    {
      return data_[col * rows_ + row];
    }

    /// Unchecked const element access (row, col)
    const Value& operator()(size_t row, size_t col) const
    {
      return data_[col * rows_ + row];
    }

    ///@}

    ///@name Buffer access (for NumPy interoperability)
    ///@{

    /// Pointer to raw data buffer (column-major storage)
    Value* data() { return data_.data(); }

    /// Const pointer to raw data buffer
    const Value* data() const { return data_.data(); }

    /// Stride between consecutive elements in same column (always 1 for column-major)
    int innerStride() const { return 1; }

    /// Stride between consecutive columns
    int outerStride() const { return static_cast<int>(rows_); }

    /// Returns false (column-major storage, not row-major)
    static constexpr bool rowMajor() { return false; }

    ///@}

    ///@name Modifiers
    ///@{

    /// Resize matrix (contents become undefined after resize)
    void resize(size_t rows, size_t cols)
    {
      rows_ = rows;
      cols_ = cols;
      data_.resize(rows * cols);
    }

    /// Fill all elements with value
    void fill(Value value)
    {
      std::fill(data_.begin(), data_.end(), value);
    }

    /// Clear matrix (set to 0x0)
    void clear()
    {
      data_.clear();
      rows_ = 0;
      cols_ = 0;
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
     * @param[in] array The 2D array containing the values to be assigned to the matrix.
     */
    template <typename T, long int ROWS, long int COLS>
    void setMatrix(T const (&array)[ROWS][COLS])
    {
      resize(ROWS, COLS);
      for (long int i = 0; i < ROWS; ++i)
      {
        for (long int j = 0; j < COLS; ++j)
        {
          (*this)(i, j) = array[i][j];
        }
      }
    }

    ///@}

    ///@name Iterators
    ///@{
    iterator begin() { return data_.begin(); }
    iterator end() { return data_.end(); }
    const_iterator begin() const { return data_.begin(); }
    const_iterator end() const { return data_.end(); }
    const_iterator cbegin() const { return data_.cbegin(); }
    const_iterator cend() const { return data_.cend(); }
    ///@}

    ///@name Reduction operations
    ///@{

    /**
     * @brief Returns the maximum value in the matrix.
     *
     * @return The maximum value. For empty matrices, returns Value() (default-constructed value).
     */
    Value maxValue() const
    {
      if (data_.empty()) return Value();
      return *std::max_element(data_.begin(), data_.end());
    }

    /**
     * @brief Returns the minimum value in the matrix.
     *
     * @return The minimum value. For empty matrices, returns Value() (default-constructed value).
     */
    Value minValue() const
    {
      if (data_.empty()) return Value();
      return *std::min_element(data_.begin(), data_.end());
    }

    ///@}

    ///@name Comparison
    ///@{

    /**
     * @brief Equality operator. Compares two matrices for equality.
     *
     * @param[in] rhs The matrix to be compared.
     * @return True if matrices are equal, false otherwise.
     *
     * @throw Exception::Precondition if matrices have different dimensions (Debug mode only)
     */
    bool operator==(const Matrix& rhs) const
    {
      OPENMS_PRECONDITION(rows_ == rhs.rows_ && cols_ == rhs.cols_,
                          "Matrices must have the same dimensions for comparison.");
      return data_ == rhs.data_;
    }

    bool operator!=(const Matrix& rhs) const
    {
      return !(*this == rhs);
    }

    ///@}

    /**
     * @brief Friend function to output the matrix to an output stream.
     *
     * @param[in,out] os Output stream.
     * @param[in] matrix Matrix to be output.
     * @return Reference to the output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const Matrix<Value>& matrix)
    {
      for (Size i = 0; i < matrix.rows(); ++i)
      {
        for (Size j = 0; j < matrix.cols(); ++j)
        {
          os << std::setprecision(6) << std::setw(6) << matrix(i, j) << ' ';
        }
        os << '\n';
      }
      return os;
    }

  private:
    std::vector<Value> data_;  ///< Column-major storage
    Size rows_;                ///< Number of rows
    Size cols_;                ///< Number of columns
  };

} // namespace OpenMS
