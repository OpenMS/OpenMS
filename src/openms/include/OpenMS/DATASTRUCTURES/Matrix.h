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
     * @param[in] rows Number of rows in the matrix.
     * @param[in] cols Number of columns in the matrix.
     * @param[in] value Initial value to fill the matrix.
     */
    Matrix(Size rows, Size cols, Value value = Value()) : EigenMatrixType(rows, cols)
    {
      fill(value);
    }

    /**
     * @brief Unchecked element access operator (const version).
     *
     * Provides direct access to the element at the specified row and column
     * without bounds checking.
     *
     * @param row Index of the row.
     * @param col Index of the column.
     * @return Const reference to the matrix element at (row, col).
     *
     * @note This operator does not perform bounds checking. Accessing elements
     *       outside the matrix dimensions results in undefined behavior.
     */
    const Value& operator()(size_t row, size_t col) const
    {
      return EigenMatrixType::operator()(row, col);
    }

    /**
     * @brief Unchecked element access operator.
     *
     * Provides direct access to the element at the specified row and column
     * without bounds checking.
     *
     * @param row Index of the row.
     * @param col Index of the column.
     * @return Reference to the matrix element at (row, col).
     *
     * @note This operator does not perform bounds checking. Accessing elements
     *       outside the matrix dimensions results in undefined behavior.
     */
    Value& operator()(size_t row, size_t col)
    {
      return EigenMatrixType::operator()(row, col);
    }

    /**
     * @brief Get matrix entry (const version).
     *
     * @param i Row index.
     * @param j Column index.
     * @return Const reference to the matrix element at (i, j).
     *
     * @note pyOpenMS can't easily wrap operator() so we provide additional getter/setter.
     */
    const Value& getValue(size_t const i, size_t const j) const
    {
      return this->operator()(i, j);
    }

    /**
     * @brief Get matrix entry.
     *
     * @param i Row index.
     * @param j Column index.
     * @return Reference to the matrix element at (i, j).
     *
     * @note pyOpenMS can't easily wrap operator() so we provide additional getter/setter.
     */
    Value& getValue(size_t const i, size_t const j)
    {
      return this->operator()(i, j);
    }

    /**
     * @brief Set matrix entry.
     *
     * @param i Row index.
     * @param j Column index.
     * @param value The value to set at position (i, j).
     *
     * @note pyOpenMS can't easily wrap operator() so we provide additional getter/setter.
     */
    void setValue(size_t const i, size_t const j, const Value& value)
    {
      this->operator()(i, j) = value;
    }

    /**
     * @brief Set all matrix elements to the given value.
     *
     * @param value The value to fill the matrix with.
     *
     * @note Affects all elements in the matrix while maintaining the current
     *       rows and columns dimensions.
     */
    void fill(const Value& value)
    {
      EigenMatrixType::fill(value);
    }

    /**
     * @brief Reset matrix to empty 0x0 state.
     *
     * Clears the internal storage and sets both row and column dimensions to 0.
     */
    void clear()
    {
      EigenMatrixType::resize(0, 0);
    }

    /**
     * @brief Resize the matrix to the specified dimensions.
     *
     * Resizes the matrix to the given number of rows and columns.
     * After resizing, the contents of the matrix are undefined and
     * the underlying storage may be reallocated/invalidated, so any
     * existing pointers or iterators into the matrix become invalid.
     *
     * @param rows New number of rows.
     * @param cols New number of columns.
     *
     * @note This method is needed for Cython/pyOpenMS compatibility.
     *       Wrappers relying on stable memory or contiguous layout must
     *       handle potential reallocation after calling this method.
     */
    void resize(size_t rows, size_t cols)
    {
      EigenMatrixType::resize(rows, cols);
    }

    /// @name Buffer access for NumPy/Cython interoperability
    /// @{

    /**
     * @brief Pointer to raw data buffer in column-major order.
     *
     * @return Pointer to the first element of the underlying data storage.
     */
    Value* data()
    {
      return EigenMatrixType::data();
    }

    /**
     * @brief Pointer to raw data buffer in column-major order (const version).
     *
     * @return Const pointer to the first element of the underlying data storage.
     */
    const Value* data() const
    {
      return EigenMatrixType::data();
    }

    /**
     * @brief Stride between consecutive elements in the same column.
     *
     * @return The inner stride (typically 1 for column-major storage).
     */
    int innerStride() const
    {
      return EigenMatrixType::innerStride();
    }

    /**
     * @brief Stride between consecutive columns.
     *
     * @return The outer stride (typically equals the number of rows for column-major storage).
     */
    int outerStride() const
    {
      return EigenMatrixType::outerStride();
    }

    /**
     * @brief Indicates storage order of the matrix.
     *
     * @return False for column-major storage (default), true for row-major.
     *
     * @note This matrix uses column-major storage, so this method always returns false.
     */
    bool rowMajor() const
    {
      return EigenMatrixType::IsRowMajor;
    }

    /// @}

    /// @name Iterator access
    /// @{
    /// Iterators traverse the matrix in column-major order: all elements of column 0,
    /// then column 1, etc.

    /**
     * @brief Returns a mutable iterator to the beginning of the data.
     *
     * Iteration proceeds in column-major order (all elements of column 0,
     * then column 1, etc.).
     *
     * @return Iterator to the first element of the underlying data container.
     */
    Value* begin()
    {
      return EigenMatrixType::data();
    }

    /**
     * @brief Returns a mutable iterator to the end of the data.
     *
     * Iteration proceeds in column-major order (all elements of column 0,
     * then column 1, etc.).
     *
     * @return Iterator to one past the last element of the underlying data container.
     */
    Value* end()
    {
      return EigenMatrixType::data() + EigenMatrixType::size();
    }

    /**
     * @brief Returns a const iterator to the beginning of the data.
     *
     * Iteration proceeds in column-major order (all elements of column 0,
     * then column 1, etc.).
     *
     * @return Const iterator to the first element of the underlying data container.
     */
    const Value* begin() const
    {
      return EigenMatrixType::data();
    }

    /**
     * @brief Returns a const iterator to the end of the data.
     *
     * Iteration proceeds in column-major order (all elements of column 0,
     * then column 1, etc.).
     *
     * @return Const iterator to one past the last element of the underlying data container.
     */
    const Value* end() const
    {
      return EigenMatrixType::data() + EigenMatrixType::size();
    }

    /**
     * @brief Returns a const iterator to the beginning of the data.
     *
     * Iteration proceeds in column-major order (all elements of column 0,
     * then column 1, etc.).
     *
     * @return Const iterator to the first element of the underlying data container.
     */
    const Value* cbegin() const
    {
      return EigenMatrixType::data();
    }

    /**
     * @brief Returns a const iterator to the end of the data.
     *
     * Iteration proceeds in column-major order (all elements of column 0,
     * then column 1, etc.).
     *
     * @return Const iterator to one past the last element of the underlying data container.
     */
    const Value* cend() const
    {
      return EigenMatrixType::data() + EigenMatrixType::size();
    }

    /// @}

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
     * @param[in] rhs The matrix to be compared.
     * @return True if matrices are equal, false otherwise.
     *
     * @throw Exception::Precondition if matrices have different dimensions (Debug mode only)
     */
    bool operator==(const Matrix& rhs) const {
      OPENMS_PRECONDITION(this->rows() == rhs.rows() && this->cols() == rhs.cols(), "Matrices must have the same dimensions for comparison.");
      return EigenMatrixType::operator==(rhs);
    }

    /**
     * @brief Inequality operator. Returns true if two matrices are not equal.
     *
     * @param rhs Matrix to compare against.
     * @return True if matrices differ in dimensions or element values, false otherwise.
     *
     * @see operator==
     */
    bool operator!=(const Matrix& rhs) const
    {
      return !(*this == rhs);
    }

    /**
     * @brief Friend function to output the matrix to an output stream.
     *
     * @param[in,out] os Output stream.
     * @param[in] matrix Matrix to be output.
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
