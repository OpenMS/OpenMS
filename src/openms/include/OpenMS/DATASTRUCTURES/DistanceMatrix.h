// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Colorizer.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/config.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace OpenMS
{

/**
  @brief A two-dimensional distance matrix, similar to OpenMS::Matrix

  Similar to OpenMS::Matrix, but contains only elements above the main
  diagonal, hence translating access with operator(,) for elements of above
  the main diagonal to corresponding elements below the main diagonal and
  returning 0 for requested elements in the main diagonal, since
  self-distance is assumed to be 0. Keeps track of the minimal element in the
  Matrix with OpenMS::DistanceMatrix::min_element_ if only for setting a
  value OpenMS::DistanceMatrix::setValue is used. Other
  OpenMS::DistanceMatrix altering methods may require a manual update by call
  of OpenMS::DistanceMatrix::updateMinElement, see the respective methods
  documentation.

  @ingroup Datastructures
*/
template<typename Value>
class DistanceMatrix
{
public:
  ///@name STL compliance type definitions
  //@{
  typedef Value value_type;
  //@}

  ///@name OpenMS compliance type definitions
  //@{
  typedef Size SizeType;
  typedef value_type ValueType;
  //@}

  /** @brief default constructor

  */
  DistanceMatrix(): matrix_(nullptr), init_size_(0), dimensionsize_(0), min_element_(0, 0)
  {
  }

  /**
    @brief detailed constructor

    @param dimensionsize the number of rows (and therewith cols)
    @param value DistanceMatrix will be filled with this element (main diagonal will still "hold" only zeros)
    @throw Exception::OutOfMemory if requested dimensionsize is to big to fit into memory
  */
  DistanceMatrix(SizeType dimensionsize, Value value = Value()):
      matrix_(new ValueType*[dimensionsize]),
      init_size_(dimensionsize),
      dimensionsize_(dimensionsize),
      min_element_(0, 0)
  {
    matrix_[0] = NULL;
    SizeType i = 1;
    for (i = 1; i < dimensionsize; ++i)
    {
      matrix_[i] = new ValueType[i];
      if (matrix_[i] == NULL)
      {
        SizeType j = i;
        for (i = 1; i < j; i++)
        {
          delete[] matrix_[i];
        }
        delete[] matrix_;
        matrix_ = NULL;
        dimensionsize_ = 0;
        init_size_ = 0;
        throw Exception::OutOfMemory(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     (UInt)((((dimensionsize - 2) * (dimensionsize - 1)) / 2) * sizeof(ValueType)));
      }
    }
    if (matrix_ != NULL)
    {
      for (i = 1; i < dimensionsize; ++i)
      {
        for (SizeType j = 0; j < i; ++j)
        {
          matrix_[i][j] = value;
        }
      }
      min_element_ = std::make_pair(1, 0);
    }
  }

  /**
    @brief copy constructor

    @param source  this DistanceMatrix will be copied
    @throw Exception::OutOfMemory if requested dimensionsize is to big to fit into memory
  */
  DistanceMatrix(const DistanceMatrix& source):
      matrix_(new ValueType*[source.dimensionsize_]),
      init_size_(source.dimensionsize_),
      dimensionsize_(source.dimensionsize_),
      min_element_(source.min_element_)
  {
    matrix_[0] = NULL;
    SizeType i = 1;
    for (i = 1; i < dimensionsize_; ++i)
    {
      matrix_[i] = new ValueType[i];
      if (matrix_[i] == NULL)
      {
        SizeType j = i;
        for (i = 1; i < j; i++)
        {
          delete[] matrix_[i];
        }
        delete[] matrix_;
        matrix_ = NULL;
        dimensionsize_ = 0;
        init_size_ = 0;
        min_element_ = std::make_pair(0, 0);
        throw Exception::OutOfMemory(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     (UInt)((((dimensionsize_ - 2) * (dimensionsize_ - 1)) / 2) * sizeof(ValueType)));
      }
    }
    if (matrix_ != NULL)
    {
      for (i = 1; i < dimensionsize_; ++i)
      {
        std::copy(source.matrix_[i], source.matrix_[i] + i, matrix_[i]);
      }
    }
  }

  /// destructor
  ~DistanceMatrix()
  {
    for (SizeType i = 1; i < init_size_; i++)
    {
      delete[] matrix_[i];
    }
    delete[] matrix_;
  }

  /**
    @brief gets a value at a given position (read only):

    @param i the i-th row
    @param j the j-th col
  */
  const ValueType operator()(SizeType i, SizeType j) const
  {
    return getValue(i, j);
  }

  /**
    @brief gets a value at a given position (read only):

    @param i the i-th row
    @param j the j-th col
  */
  ValueType operator()(SizeType i, SizeType j)
  {
    return getValue(i, j);
  }

  /**
    @brief gets a value at a given position:

    @param i the i-th row
    @param j the j-th col
    @throw Exception::OutOfRange if given coordinates are out of range
  */
  const ValueType getValue(SizeType i, SizeType j) const
  {
    if (i >= dimensionsize_ || j >= dimensionsize_) { throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
    // elements on main diagonal are not stored and assumed to be 0
    if (i == j) { return 0; }
    if (i < j) { std::swap(i, j); }
    return (const ValueType)(matrix_[i][j]);
  }

  /**
    @brief gets a value at a given position:

    @param i the i-th row
    @param j the j-th col
    @throw Exception::OutOfRange if given coordinates are out of range
  */
  ValueType getValue(SizeType i, SizeType j)
  {
    if (i >= dimensionsize_ || j >= dimensionsize_) { throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
    // elements on main diagonal are not stored and assumed to be 0
    if (i == j) { return 0; }
    if (i < j) { std::swap(i, j); }
    return matrix_[i][j];
  }

  /**
    @brief sets a value at a given position:

    @param i the i-th row
    @param j the j-th col
    @param value the set-value
    @throw Exception::OutOfRange if given coordinates are out of range
  */
  void setValue(SizeType i, SizeType j, ValueType value)
  {
    if (i >= dimensionsize_ || j >= dimensionsize_) { throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
    // elements on main diagonal are not stored and assumed to be 0
    if (i != j)
    {
      if (i < j) { std::swap(i, j); }
      if (i != min_element_.first && j != min_element_.second)
      {
        matrix_[i][j] = value;
        if (value < matrix_[min_element_.first][min_element_.second]) // keep min_element_ up-to-date
        {
          min_element_ = std::make_pair(i, j);
        }
      }
      else
      {
        if (value <= matrix_[min_element_.first][min_element_.second]) { matrix_[i][j] = value; }
        else
        {
          matrix_[i][j] = value;
          updateMinElement();
        }
      }
    }
  }

  /**
    @brief sets a value at a given position:

    @param i the i-th row
    @param j the j-th col
    @param value the set-value
    @throw Exception::OutOfRange if given coordinates are out of range

    possible invalidation of min_element_ - make sure to update before further usage of matrix
  */
  void setValueQuick(SizeType i, SizeType j, ValueType value)
  {
    if (i >= dimensionsize_ || j >= dimensionsize_) { throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
    // elements on main diagonal are not stored and assumed to be 0
    if (i != j)
    {
      if (i < j) { std::swap(i, j); }
      matrix_[i][j] = value;
    }
  }

  /// reset all
  void clear()
  {
    for (SizeType i = 1; i < init_size_; i++)
    {
      delete[] matrix_[i];
    }
    delete[] matrix_;
    matrix_ = nullptr;
    min_element_ = std::make_pair(0, 0);
    dimensionsize_ = 0;
    init_size_ = 0;
  }

  /**
    @brief resizing the container

    @param dimensionsize the desired number of rows (and therewith cols)
    @param value which the matrix will be filled with
    @throw Exception::OutOfMemory thrown if size of DistanceMatrix requested does not fit into memory

    invalidates all content
  */
  void resize(SizeType dimensionsize, Value value = Value())
  {
    for (SizeType j = 1; j < init_size_; j++)
    {
      delete[] matrix_[j];
    }
    delete[] matrix_;
    dimensionsize_ = dimensionsize;
    init_size_ = dimensionsize;
    min_element_ = std::make_pair(0, 0);
    matrix_ = new ValueType*[dimensionsize_];
    for (SizeType j = 1; j < dimensionsize_; ++j)
    {
      matrix_[j] = new ValueType[j];
      if (matrix_[j] == nullptr)
      {
        for (SizeType k = 1; k < j; ++k)
        {
          delete[] matrix_[k];
        }
        delete[] matrix_;
        matrix_ = nullptr;
        dimensionsize_ = 0;
        init_size_ = 0;
        throw Exception::OutOfMemory(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     (UInt)((((dimensionsize_ - 2) * (dimensionsize_ - 1)) / 2) * sizeof(Value)));
      }
    }
    if (matrix_ != nullptr)
    {
      for (SizeType j = 0; j < dimensionsize; ++j)
      {
        for (SizeType k = 0; k < j; ++k)
        {
          matrix_[j][k] = value;
        }
      }
      min_element_ = std::make_pair(1, 0);
    }
  }

  /**
    @brief reduces DistanceMatrix by one dimension. first the jth row, then jth column

    @param j the jth row (and therewith also jth col) to be removed
    @throw Exception::OutOfRange if @p j is grater than the greatest row number

    May invalidates min_element_, make sure to update min_element_ if necessary before used
  */
  void reduce(SizeType j)
  {
    if (j >= dimensionsize_) { throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
    // delete row j and therefor overwrite with row j+1 and iterate like this to last row
    SizeType i = j + 1;
    while (i < dimensionsize_ && matrix_[i] != nullptr)
    {
      // left out in the copy is each rows jth element, pointer working here as iterators just fine
      std::copy(matrix_[i] + j + 1, matrix_[i] + i, std::copy(matrix_[i], matrix_[i] + j, matrix_[i - 1]));
      ++i;
    }
    // last row is freed and the pointer set to NULL (outer array's size is not changed)
    delete[] matrix_[i - 1];
    matrix_[i - 1] = nullptr;
    --dimensionsize_;
  }

  /// gives the number of rows (i.e. number of columns)
  SizeType dimensionsize() const
  {
    return dimensionsize_;
  }

  /**
    @brief keep track of the actual minimum element after altering the matrix

    @throw Exception::OutOfRange thrown if there is no element to access
  */
  void updateMinElement()
  {
    min_element_ = std::make_pair(1, 0);
    // error if dimensionsize_<1, return if dimensionsize_ == 1, else
    if (dimensionsize_ < 1) { throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
    if (dimensionsize_ != 1) // else matrix has one element: (1,0)
    {
      ValueType* row_min_;
      for (SizeType r = 2; r < dimensionsize_ && matrix_[r] != nullptr; ++r)
      {
        row_min_ = std::min_element(matrix_[r], matrix_[r] + r);
        if (*row_min_ < matrix_[min_element_.first][min_element_.second]) { min_element_ = std::make_pair(r, row_min_ - matrix_[r]); }
      }
    }
  }

  /**
    @brief Equality comparator.

    @throw Exception::Precondition thrown if given DistanceMatrix is not compatible in size
  */
  bool operator==(DistanceMatrix<ValueType> const& rhs) const
  {
    OPENMS_PRECONDITION(dimensionsize_ == rhs.dimensionsize_, "DistanceMatrices have different sizes.");
    for (Size i = 1; i < rhs.dimensionsize(); ++i)
    {
      for (Size j = 0; j < i; ++j)
      {
        if (matrix_[i][j] != rhs.matrix_[i][j]) { return false; }
      }
    }
    return true;
  }

  /**
    @brief Indexpair of minimal element

    @throw Exception::OutOfRange thrown if there is no element to access
  */
  std::pair<SizeType, SizeType> getMinElementCoordinates() const
  {
    if (dimensionsize_ == 0) { throw Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
    return min_element_;
  }

protected:
  /// sparse element not to be included in base container
  ValueType** matrix_;
  /// number of actually stored rows
  SizeType init_size_; // actual size of outer array
  /// number of accessibly stored rows (i.e. number of columns)
  SizeType dimensionsize_; // number of virtual elements: ((dimensionsize-1)*(dimensionsize))/2
  /// index of minimal element(i.e. number in underlying SparseVector)
  std::pair<SizeType, SizeType> min_element_;

private:
  /// assignment operator (unsafe)
  DistanceMatrix& operator=(const DistanceMatrix& rhs)
  {
    matrix_ = rhs.matrix_;
    init_size_ = rhs.init_size_;
    dimensionsize_ = rhs.dimensionsize_;
    min_element_ = rhs.min_element_;

    return *this;
  }

}; // class DistanceMatrix

/**
  @brief Print the contents to a stream (and colors the diagonal, if the stream is cout/cerr)

  @relatesalso DistanceMatrix
*/
template<typename Value>
std::ostream& operator<<(std::ostream& os, const DistanceMatrix<Value>& matrix)
{
  using SizeType = typename DistanceMatrix<Value>::SizeType;

  // we need to print a square matrix. So we set the width
  //std::ios_base::fmtflags flag_backup = os.setf(std::ios::scientific); // 'scientific' messes with the width; don't do it
  std::streamsize precision_backup = os.precision(6); // we could go with `writtenDigits<Value>(Value())`, but it becomes unreadable...
  auto width_backup = os.width(8);

  for (SizeType i = 0; i < matrix.dimensionsize(); ++i)
  {
    for (SizeType j = 0; j < matrix.dimensionsize(); ++j)
    {
      if (i == j)
      { // color the diagonal in red (conditional, see Colorizer)
        os << red(matrix(i, j)) << '\t';
      }
      else 
      {
        os << matrix(i, j) << '\t';
      }
    }
    os << '\n';
  }
  //os.flags(flag_backup);
  os.precision(precision_backup);
  os.width(width_backup);
  return os;
}

} // namespace OpenMS
