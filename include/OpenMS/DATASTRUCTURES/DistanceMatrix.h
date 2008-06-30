// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DISTANCEMATRIX_H
#define OPENMS_DATASTRUCTURES_DISTANCEMATRIX_H

#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/SparseVector.h>

#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iostream>


namespace OpenMS
{

  /**
	@brief A two-dimensional distance matrix, similar to OpenMS::Matrix
	
	similar to OpenMS::Matrix, but contains only elements above the main diagonal, hence translating access with operator(,) 
	for elements of lower triangular matrix to corresponing elements in upper triangular matrix and returning 0 for requested 
	elements in the main diagonal, since selfdistance is assumed to be 0. Inherits OpenMS::SparseVector and is therefor optimal 
	for data with lots of redundant values. Keeps track of the minimal element in the Matrix with OpenMS::DistanceMatrix::min_element_
	if only for setting a value OpenMS::DistanceMatrix::setValue is used. Other Matrix altering functions may require a maual update 
	by call of OpenMS::DistanceMatrix::updateMinElement
    
	@ingroup Datastructures
  */
  template <typename Value>
	class DistanceMatrix : SparseVector<Value>
	{
	 protected:
		typedef SparseVector<Value> Base;

	 public:

		///@name STL compliance type definitions
		//@{
		typedef Base container_type;
		typedef Value value_type;

		typedef typename Base::difference_type difference_type;
		typedef typename Base::size_type size_type;
		typedef typename Base::const_reference const_reference;
		typedef typename Base::reference reference;
		typedef typename Base::pointer pointer;
		typedef typename Base::allocator_type allocator_type;

		typedef typename Base::const_iterator const_iterator;
		typedef typename Base::const_reverse_iterator const_reverse_iterator;
		typedef typename Base::iterator iterator;
		typedef typename Base::reverse_iterator reverse_iterator;
		//@}

		///@name OpenMS compliance type definitions
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

		/// constructor @param se the sparse element
		DistanceMatrix (Value se=1)
			: Base(se),
				dimensionsize_(0), min_element_(0)
		{
		}

		/** @brief detailed constructor, discouraged unless made safe that filling element is same as sparse element
		
			@param dimensionsize the number of rows (and therewith cols)
			@param value Matrix will be filled with this element (main diagonal will still "hold" only zeros)
			@param se the sparse element (@see OpenMS::SparseVector)
		*/
		DistanceMatrix (SizeType dimensionsize, Value value = Value(), Value se=1)
			: Base(((dimensionsize-1)*(dimensionsize))/2,value,se),
				dimensionsize_(dimensionsize-1), min_element_(0)
		{
		}

		/// copy constructor
		DistanceMatrix (const DistanceMatrix& source)
			: Base(source),
				dimensionsize_(source.dimensionsize_), min_element_(source.min_element_)
		{
		}

		/// assignment operator
		DistanceMatrix& operator= (const DistanceMatrix& rhs)
		{
			Base::operator= (rhs);
			dimensionsize_ = rhs.dimensionsize_;
			min_element_ = rhs.min_element_;
			return *this;
		}
		
		/// destructor
		~DistanceMatrix() 
		{
		}

		/** @brief gets a value at a given position:
			
			@param i the i-th row 
			@param j the j-th col
		*/
		const value_type operator() (size_type const i, size_type const j) const
		{
			return /*const_cast<value_type>*/(getValue(i,j));
		}

		/** @brief gets a value at a given position: 
			
			@param i the i-th row
			@param j the j-th col
		*/
		value_type operator() (size_type const i, size_type const j)
		{
			return getValue(i,j);
		}

		/** @brief gets a value at a given position: 
			
			@param i the i-th row
			@param j the j-th col
		*/
		const value_type getValue(size_type const i, size_type const j) const
		{
			// elements on main diagonal are not stored and assumed to be 0
			if(i==j)
			{
				return 0;
			}
			return Base::at(index(i,j));
		}

		/** @brief gets a value at a given position: 
			
			@param i the i-th row
			@param j the j-th col
		*/
		value_type getValue(size_type const i, size_type const j)
		{
			// elements on main diagonal are not stored and assumed to be 0
			if(i==j)
			{
				return 0;
			}
			return Base::at(index(i,j));
		}
		
		/** @brief sets a value at a given position: 
			
			@param i the i-th row
			@param j the j-th col
			@param value the set-value
		*/
		void setValue(size_type const i, size_type const j, value_type value)
		{
			// elements on main diagonal are not stored and assumed to be 0
			if(i!=j)
			{
				UInt pos = index(i,j); 
				
				//this is for keeping min_element_ position in underlying SparseVector up to date
				if(value < Base::at(min_element_))
				{
					Base::operator[](pos) = value; 
					min_element_ = pos;
				}
				else // value >=
				{
					//same as min_element, but maybe at a earlier position
					if(value == Base::at(min_element_))
					{
						Base::operator[](pos) = value; 
						min_element_ = min(min_element_,pos);
					}
					else // value >
					{
						Base::operator[](pos) = value;
						//overwriting min_element_
						if (pos == min_element_)
						{
							updateMinElement();	
						}
					}
				}
			}
		}

		/** @brief sets a value at a given position: 
			
			@param i the i-th row
			@param j the j-th col
			@param value the set-value
			
			possible invalidation of min_element_ - make sure to update before further usage of matrix
		*/
		void setValueQuick(size_type const i, size_type const j, value_type value)
		{
			// elements on main diagonal are not stored and assumed to be 0
			if(i!=j)
			{
				UInt pos = index(i,j); 
				Base::operator[](pos) = value;
			}
		}

		/**
			 @name Pure access declarations

			 These make begin(), end() etc. from container_type accessible.
		*/
		//@{
		Base::begin;
		Base::end;
		Base::rbegin;
		Base::rend;
		//@}

		/// reset all (except allocated mem.)
		void clear()
		{
			Base::clear();
			dimensionsize_ = 0;
			min_element_ = 0;
		}

		/// resizing the container (invalidates content)
		void resize(size_type i) throw (Exception::OutOfRange)
		{
			if (i <= 1) 
			{
				throw Exception::OutOfRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);
			}
			dimensionsize_ = i-1;
			Base::resize((dimensionsize_*(dimensionsize_+1))/2);
			if(min_element_ >= (dimensionsize_*(dimensionsize_+1))/2)
			{
				updateMinElement();
			}
		}
		
		/// reduces triangular matrix by one dimension. first the jth row, then jth collumn - invalidates min_element_ - make sure to update before used @param j j-th row/col to be reduced
		void reduce(size_type j)
		{
			//behind last element in row j 
			UInt x = index(j,dimensionsize_);
		 	if(j!=dimensionsize_)
		  	{
				//delete row j
				Base::erase(Base::begin()+x+1-(dimensionsize_-j),Base::begin()+x+1);
		  	}
		  	if(j!=0)
		  	{
				//delete col j
				iterator it = Base::begin()+x-(dimensionsize_-j);
				for(UInt c=0; c<(j); ++c)
				{
					it = Base::erase(it-(dimensionsize_-j)-c);	
				}
	  		}	
			--dimensionsize_;
			
			//updateMinElement();
		}
			
		/// gives the number of rows (i.e. number of collumns)
		SizeType dimensionsize() const
		{
			return dimensionsize_+1;
		}
		
		/// keep track of the actual minimum element after altering the matrix
		void updateMinElement() throw (Exception::OutOfRange)
		{
	    		iterator pos = Base::getMinElement();
	    		if(pos==Base::end())
	    		{
			    	throw Exception::OutOfRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);
	    		}
	    		min_element_ = pos.position();
		}	
				
		/// Indexpair of minimal element
		std::pair<UInt,UInt> getMinElementCoordinates() const throw(Exception::IndexOverflow)
		{
			if ( Base::size() == 0 ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__);
			return indexPair(min_element_);
		}

		/// Calculate the index into the underlying vector from row and column.  Note that Matrix uses the (row,column) lexicographic ordering for indexing.
		UInt index(UInt row, UInt col) const
		{
		#ifdef DISTANCEMATRIX_DEBUG
			if ( row > dimensionsize_ ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,row, dimensionsize_);
			if ( col > dimensionsize_ ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,col, dimensionsize_);
		#endif
			// (i,j) -> (i,j-1)
			if(row<=col)
			{
				return (((row+1)*dimensionsize_)-(dimensionsize_ - (col))-(row*(row+1)/2)-1);
			}
			else //row>col
			{
				return (((col+1)*dimensionsize_)-(dimensionsize_ - (row))-(col*(col+1)/2)-1);
			}
		}

	
		/// Calculate the row and column from an index into the underlying vector. Note that Matrix uses the (row,column) lexicographic ordering for indexing.
		std::pair<UInt,UInt> const indexPair(UInt index) const
		{
		#ifdef DISTANCEMATRIX_DEBUG
			if ( index >= dimensionsize_ ) throw Exception::IndexOverflow(__FILE__,__LINE__,__PRETTY_FUNCTION__,index,size()-1);
		#endif
			//index counting from 0, formula from 1!
			//rowindex: ceiling[dimensionsize_ + 0,5 -sqrt(dimensionsize_² + dimensionsize_ + 0,25-2*index+1)]
			//colindex: index - (rowindex-1 * (dimensionsize_ -(rowindex-1)/2))
			double row_index = ceil( 
						dimensionsize_ + 0.5 - sqrt( 
										(dimensionsize_*dimensionsize_) + dimensionsize_ + 0.25 - (2 * double(index+1) ) 
										) 
						);
			
			double col_index = (index+1) - ((row_index-1) * (dimensionsize_ - (row_index)/2 ) );

			//(i,j) -> (i,j+1)
			return std::pair<UInt,UInt>(UInt(row_index-1),UInt(col_index));
		}


		/**@brief Equality comparator.
		
			If matrices have different row or column numbers, throws a precondition exception.
		*/
		bool operator== ( DistanceMatrix<Value> const & rhs ) const
			throw (Exception::Precondition)
		{
			OPENMS_PRECONDITION(dimensionsize_ == rhs.dimensionsize_,"DistanceMatrices have different sizes.");
 			return static_cast < typename DistanceMatrix<Value>::Base const &>(*this) == static_cast < typename DistanceMatrix<Value>::Base const &>(rhs);
		}

		/**@brief less-than comparator.  Comparison is done lexicographically: first by row, then by column.

			If matrices have different row or column numbers, throws a precondition exception.
		*/
		bool operator< (DistanceMatrix<Value> const & rhs) const
			throw (Exception::Precondition)
		{
			OPENMS_PRECONDITION(dimensionsize_ == rhs.dimensionsize_,"DistanceMatrices have different sizes.");
			return static_cast < typename DistanceMatrix<Value>::Base const &>(*this) < static_cast < typename DistanceMatrix<Value>::Base const &>(rhs);
		}

	 protected:

		/// number of actually stored rows (i.e. number of columns)
		SizeType dimensionsize_;
		/// index of minimal element(i.e. number in underlying SparseVector)
		UInt min_element_;


	}; // class DistanceMatrix

	/**@brief Print the contents to a stream.

		@relatesalso DistanceMatrix
	*/
	template <typename Value>
	std::ostream& operator<< (std::ostream& os, const DistanceMatrix<Value>& matrix)
	{
		typedef typename DistanceMatrix<Value>::size_type size_type;
		//evtl. color lower triangular matrix o.s.l.t.
		for ( size_type i = 0; i < matrix.dimensionsize(); ++i ) {
			for ( size_type j = 0; j < matrix.dimensionsize(); ++j ) {
				os << std::setprecision(6) << std::setw(6) << matrix(i,j) << ' ';
			}
			os << std::endl;
		}
		return os;
	}


} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_DISTANCEMATRIX_H
