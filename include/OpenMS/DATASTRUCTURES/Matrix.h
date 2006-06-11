// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Id: Matrix.h,v 1.5 2006/03/08 17:48:20 groepl Exp $
// $Author: groepl $
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATSTRUCTURES_MATRIX_H
#define OPENMS_DATSTRUCTURES_MATRIX_H

#include <OpenMS/CONCEPT/Macros.h>

#include <iomanip>
#include <vector>

namespace OpenMS
{
  
  /**
     @brief A two-dimensional matrix.  Similar to std::vector, but uses
     operator() for element access.
      
     This is not intended to be used for linear algebra.  Think of it as a
     random access container.

		 (Yes, one could overload operator[] to make things like mat[i][j] work.
		 But mat(i,j) isn't so bad, either.)

  */
  template <typename Value>
	class Matrix : protected std::vector < Value >
	{
	 protected:
		typedef std::vector < Value > Base;
		
	 public:

		///@brief STL compliance type definitions
		//@{
		typedef Base container_type;

		typedef typename Base::difference_type difference_type;
		typedef typename Base::size_type size_type;

		typedef typename Base::const_iterator const_iterator;
		typedef typename Base::const_reverse_iterator const_reverse_iterator;
		typedef typename Base::iterator iterator;
		typedef typename Base::reverse_iterator reverse_iterator;

		typedef typename Base::const_reference const_reference;
		typedef typename Base::pointer pointer;
		typedef typename Base::reference reference;
		typedef typename Base::value_type value_type;

		typedef typename Base::allocator_type allocator_type;
		//@}

		///@brief OpenMS compliance type definitions
		//@{
		typedef Base ContainerType;
		typedef difference_type DifferenceType;
		typedef size_type SizeType;

		typedef const_iterator ConstIterator;
		typedef const_reverse_iterator ConstReverseIterator;
		typedef iterator Iterator;
		typedef reverse_iterator ReverseIterator;

		typedef const_reference ConstReference;
		typedef pointer Pointer;
		typedef reference Reference;
		typedef value_type ValueType;

		typedef allocator_type AllocatorType;
		//@}

		///@name Constructors, assignment, and destructor
		//@{
		Matrix ()
			: Base(),
				rows_(0),
				cols_(0)
		{}

		Matrix (SizeType rows, SizeType cols, ValueType value = ValueType())
			: Base(rows*cols,value),
				rows_(rows),
				cols_(cols)
		{}

		Matrix (const Matrix & source)
			: Base(source),
				rows_(source.rows_),
				cols_(source.cols_)
		{}

		Matrix & operator = (const Matrix & rhs)
		{
			Base::operator=(rhs);
			rows_ = rhs.rows_;
			cols_ = rhs.cols_;
			return *this;
		}

		~Matrix() {}
		//@}

		///@name Accessors
		//@{
		const_reference operator() (size_type const i, size_type const j) const
		{
			return getValue(i,j);
		}

		reference operator() (size_type const i, size_type const j)
		{
			return getValue(i,j);
		}

    const_reference getValue(size_type const i, size_type const j) const
		{
			return Base::operator[](index(i,j));
		}

    reference getValue(size_type const i, size_type const j)
		{
			return Base::operator[](index(i,j));
		}

    void setValue(size_type const i, size_type const j, value_type value)
		{
			Base::operator[](index(i,j)) = value;
		}

		//@}


		/**
			 @name Pure access declarations
			 
			 These make begin(), end() etc. from container_type accessible.

			 \todo Fix check_test so that it will parse Matrix.h !  Workaround: To
			 see if Matrix_test.C is really complete, comment these access
			 declarations out.
		*/
		//@{
	 public:

		Base::begin;
		Base::end;
		Base::rbegin;
		Base::rend;

		Base::front;
		Base::back;
	
		Base::empty;
		Base::size;

		Base::capacity;
		Base::max_size;

		//@}

		void clear()
		{
			Base::clear();
			rows_ = 0;
			cols_ = 0;
		}

		void resize(size_type i, size_type j, value_type value = value_type())
		{
			rows_ = i;
			cols_ = j;
			Base::resize(rows_*cols_, value);
		}

		void resize(std::pair<Size,Size> const & size_pair, value_type value = value_type())
		{
			rows_ = size_pair.first;
			cols_ = size_pair.second;
			Base::resize(rows_*cols_, value);
		}

		/// Number of rows
		SizeType rows() const throw()
		{
			return rows_;
		}

		/// Number of columns
		SizeType cols() const throw()
		{
			return cols_;
		}

		std::pair<Size,Size> sizePair() const
		{
			return std::pair<Size,Size>(rows_,cols_);
		}

		/**@brief Calculate the index into the underlying vector from row and
			 column.  Note that Matrix uses the (row,columnm) lexicographic ordering
			 for indexing.
		*/
		SizeType const index(SizeType row, SizeType col) const
		{
#ifdef OPENMS_DEBUG
			if ( row >= rows_ ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,row, rows_);
			if ( col >= cols_ ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,col, cols_);
#endif
			return row * cols_ + col;
		}

		/**@brief Calculate the row and column from an index into the underlying
			 vector.  Note that Matrix uses the (row,columnm) lexicographic ordering
			 for indexing.
		*/
		std::pair<Size,Size> const indexPair(Size index) const
		{
#ifdef OPENMS_DEBUG
			if ( index >= size() ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,index,size()-1);
#endif
			return std::pair<SizeType,SizeType>(index/cols_,index%cols_);
		}

		/**@brief Calculate the column from an index into the underlying vector.
			 Note that Matrix uses the (row,columnm) lexicographic ordering for
			 indexing.
		*/
		SizeType colIndex(SizeType index) const
		{
#ifdef OPENMS_DEBUG
			if ( index >= size() ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,index,size()-1);
#endif
			return index%cols_;
		}

		/**@brief Calculate the row from an index into the underlying vector.
			 Note that Matrix uses the (row,columnm) lexicographic ordering for
			 indexing.
		*/
		SizeType rowIndex(SizeType index) const
		{
#ifdef OPENMS_DEBUG
			if ( index >= size() ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,index,size()-1);
#endif
			return index/cols_;
		}


	 protected:

		///@name Data members
		//@{
		/// Number of rows (height of a column)
		SizeType rows_; 
		/// Number of columns (width of a row)
		SizeType cols_;
		//@}

	}; // class Matrix

	template <typename Value>
	bool operator == (Matrix<Value> const & left, Matrix<Value> const & right)
		throw (Exception::Precondition)
	{
		OPENMS_PRECONDITION(left.cols_ == right.cols_,
												"Matrices have different row sizes.");
		OPENMS_PRECONDITION(left.rows_ == right.rows_,
												"Matrices have different column sizes.");
		return left.operator==(right);
	}
	
	template <typename Value>
	bool operator < (Matrix<Value> const & left, Matrix<Value> const & right)
		throw (Exception::Precondition)
	{
		OPENMS_PRECONDITION(left.cols_ == right.cols_,
												"Matrices have different row sizes.");
		OPENMS_PRECONDITION(left.rows_ == right.rows_,
												"Matrices have different column sizes.");
		return left.operator<(right);
	}
	
	///Print the contents to a stream.
	template <typename Value>
	std::ostream& operator << (std::ostream& os, const Matrix<Value>& matrix)
	{
		typedef typename Matrix<Value>::size_type size_type;
	  for ( size_type i = 0; i < matrix.rows(); ++i ) {
	    for ( size_type j = 0; j < matrix.cols(); ++j ) {
				os << std::setprecision(6) << std::setw(6) << matrix(i,j) << ' ';
	    }
	    os << std::endl;
	  }
		return os;
	}


} // namespace OpenMS

#endif // OPENMS_DATSTRUCTURES_MATRIX_H

