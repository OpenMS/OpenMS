// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_KERNEL_COMPARATORUTILS_H
#define OPENMS_KERNEL_COMPARATORUTILS_H

#include <functional>

/**
	@defgroup ComparatorUtils ComparatorUtils

	@brief A collection of utilities for comparators.

	@ingroup Kernel

	@code
	#include <OpenMS/KERNEL/ComparatorUtils.h>
	@endcode

	This file contains some lightweight class templates
	which simplify the (re-)usage and composition of <b>comparator classes</b>:
	- ReverseComparator (reverse the direction of comparison)
	- LexicographicComparator (combine comparators lexicographically)
	- PointerComparator (compare pointers like the type they point to)
	.
	We provide corresponding "make-functions" so that you will not need to write
	out the type names in the template instantiation.

	We explain this with a simple example.  First a few prerequisites.

	@dontinclude Tutorial_ComparatorUtils.C
	@until String;

	The class IntRealString has three data members.  The class
	IntRealStringVector is a vector of IntRealString.
	We add some print() methods for convenience.
	@until };
	@until };

	Now we will exercise various ways of sorting such a vector.

	Of course, we could use the std::sort algorithm with a comparison function
	like the following.
	@until }

	This is straightforward but does not generalize well.  Instead we introduce
	three <b>comparator classes</b>:
	@until };
	@until };
	@until };
	Note how the std::binary_function template provides the necessary type
	information (sometimes this is called "reflection").

	Now we show various uses of the reverseComparator and lexicographicComparator
	function templates.
	@until vec.print();
	@until vec.print();
	@until vec.print();
	@until vec.print();
	@until vec.print();
	@until vec.print();

	And here is an application of the pointerComparator function template:
	@until main

	The output of the example program is:
	@code
	After initialization:
	(1, 4.5, paul)
	(2, 4.5, josie)
	(1, 4.5, john)
	(2, 3.9, kim)

	Sorted using lessByInt function:
	(1, 4.5, paul)
	(1, 4.5, john)
	(2, 4.5, josie)
	(2, 3.9, kim)

	Sorted using LessByInt comparator class:
	(1, 4.5, paul)
	(1, 4.5, john)
	(2, 4.5, josie)
	(2, 3.9, kim)

	Sorted using reversed LessByInt comparator class:
	(2, 4.5, josie)
	(2, 3.9, kim)
	(1, 4.5, paul)
	(1, 4.5, john)

	Sorted using lexicographic order: 1. LessByInt, 2. LessByReal
	(1, 4.5, paul)
	(1, 4.5, john)
	(2, 3.9, kim)
	(2, 4.5, josie)

	Sorted using lexicographic order: 1. reversed LessByInt, 2. LessByReal, 3. LessByString
	(2, 3.9, kim)
	(2, 4.5, josie)
	(1, 4.5, john)
	(1, 4.5, paul)

	ptr_vec before sorting
	(2, 3.9, kim)
	(2, 4.5, josie)
	(1, 4.5, john)
	(1, 4.5, paul)

	ptr_vec after sorting with pointerComparator(LessByString())
	(1, 4.5, john)
	(2, 4.5, josie)
	(2, 3.9, kim)
	(1, 4.5, paul)

	@endcode

	Note that pointerComparator can also be used with full-blown iterator classes.
	(It should work with everything that provides an operator*(), but the typedefs will always be pointers.)

	Note that these templates can also be used with different types for "left" and "right".
*/

namespace OpenMS
{

  /**
  	@brief Wrapper that takes a comparator for `something' and makes a
		comparator for <i>pointers</i> to `something' out of it.  Normally you should use
		the make-function pointerComparator() because then you do not need to
		specify the template arguments.

		This works by dereferencing the arguments (unary operator*) before
   	comparing them.
		<br> E.g. you can use
   	<code>PointerComparator<Peak1D::IntensityLess></code> to compare
   	<code>Peak1D*</code> in the same way as
   	<code>Peak1D::IntensityLess</code> works for <code>Peak1D</code> .
  */
  template < class Cmp >
	struct PointerComparator
    : public std::binary_function<typename Cmp::first_argument_type *, typename Cmp::second_argument_type *, typename Cmp::result_type>
	{
		PointerComparator(PointerComparator const& pCmp) : cmp_ ( pCmp.cmp_ ) {}
 		PointerComparator(Cmp const& cmp = Cmp() ) : cmp_(cmp) {}

		template < typename T1, typename T2 >
		typename Cmp::result_type
		operator () ( T1 left, T2 right) const
		{
			return cmp_ ( *left, *right ); // T must have operator* defined
		}

	 protected:
		Cmp const & cmp_;
	};

	/**
		@brief  Make-function to create a PointerComparator from another comparator without the need to specify the template arguments.

		@relatesalso PointerComparator

		For example,
	  <pre>
	  int i = 88, j = 99;
	  if ( pointerComparator(std::less<int>())(&i,&j) )
	  {
	    //&nbsp;   // yes, 88 < 99.
	  }
	  </pre>
	*/
  template < class Cmp >
	PointerComparator < Cmp> pointerComparator ( Cmp const & cmp )
	{ return PointerComparator < Cmp > ( cmp ); }



	//======================================================================



  /**
  	@brief Wrapper that reverses (exchanges) the two arguments of a comparator.
		Normally you should use the make-function reverseComparator()
		because then you do not need to specify the template arguments.


		For example, <code>ReverseComparator<less<T> ></code>  works like  <code>greater<T></code> .
  */
  template < class Cmp >
  struct ReverseComparator
    : std::binary_function<typename Cmp::second_argument_type, typename Cmp::first_argument_type, typename Cmp::result_type>
	// (Note that here we must reverse the order of template args!)
  {
		ReverseComparator(ReverseComparator const& cmp) : cmp_(cmp.cmp_) {}

		ReverseComparator(Cmp const& cmp = Cmp()) : cmp_(cmp) {}

		template < typename T1, typename T2 >
		typename Cmp::result_type
		operator () ( T1 left, T2 right) const
		{
			return cmp_ ( right, left ); // the other way round
		}

	 protected:
		Cmp const & cmp_;
	};

	/**
		@brief  Make-function to create a ReverseComparator from another comparator without the need to specify the template arguments.

		@relatesalso ReverseComparator

		For example,
	  <pre>
	  int i = 88, j = 99;
	  if ( reverseComparator(std::less<int>())(j,i) )
	  {
	    //&nbsp;   // yes, 99 > 88.
	  }
	  </pre>
	*/
  template < class Cmp >
	ReverseComparator < Cmp> reverseComparator ( Cmp const & cmp )
	{ return ReverseComparator < Cmp > ( cmp ); }



	//======================================================================


	/**
		@brief A wrapper class that combines two comparators lexicographically.
		Normally you should use the make-function lexicographicComparator()
		because then you do not need to specify the template arguments.

		Both comparators should of course have the same argument types.  The
		result_type is bool, that is, we perform a two-way comparison like
		<code>less<></code> and its relatives.
	*/
	template < typename Cmp1, typename Cmp2 >
	struct LexicographicComparator
		: std::binary_function<typename Cmp1::first_argument_type, typename Cmp1::second_argument_type, bool >
	{
		LexicographicComparator(Cmp1 const& cmp1 = Cmp1(), Cmp2 const& cmp2 = Cmp2()) : cmp1_(cmp1), cmp2_(cmp2) {}

 		template < typename T1, typename T2 >
		bool
		operator () ( T1 left, T2 right) const
		{
			if (cmp1_(left, right))
			{
				return true;
			}
			else
			{
				if (cmp1_(right, left))
				{
					return false;
				}
				else
				{
					return cmp2_(left, right);
				}
			}
		}

	 protected:
		Cmp1 const & cmp1_;
		Cmp2 const & cmp2_;
	};

	/**
		@brief  Make-function to create a LexicographicComparator from two other comparators without the need to specify the template arguments.

		@relatesalso LexicographicComparator

		The usage is similar to pointerComparator() or reverseComparator(), which see.
	*/
	template < typename Cmp1, typename Cmp2 >
	LexicographicComparator < Cmp1, Cmp2 > lexicographicComparator ( Cmp1 const & cmp1, Cmp2 const & cmp2 )
	{ return LexicographicComparator < Cmp1, Cmp2 > ( cmp1, cmp2 ); }

	//======================================================================



	/**
		@brief Class for comparison of std::pair using first ONLY e.g. for use with std::sort
	*/
	template<typename PairType>
	struct PairComparatorFirstElement
	: std::binary_function<PairType, PairType, bool >
	{
			bool operator() (const PairType& left, const PairType& right)const
			{
					return ( left.first < right.first ) ;
			}

	};

	/**
		@brief Class for comparison of std::pair using second ONLY e.g. for use with std::sort
	*/
	template<typename PairType>
	struct PairComparatorSecondElement
	: std::binary_function<PairType, PairType, bool >
	{
			bool operator() (const PairType& left, const PairType& right)const
			{
					return ( left.second < right.second ) ;
			}

	};

	/**
		@brief Class for comparison of std::pair using first ONLY e.g. for use with std::sort
	*/
	template<typename PairType>
	struct PairComparatorFirstElementMore
	: std::binary_function<PairType, PairType, bool >
	{
			bool operator() (const PairType& left, const PairType& right)const
			{
					return ( left.first > right.first ) ;
			}

	};

	/**
		@brief Class for comparison of std::pair using second ONLY e.g. for use with std::sort
	*/
	template<typename PairType>
	struct PairComparatorSecondElementMore
	: std::binary_function<PairType, PairType, bool >
	{
			bool operator() (const PairType& left, const PairType& right)const
			{
					return ( left.second > right.second ) ;
			}

	};

	//======================================================================

	/**
		@brief Class for comparison of std::pair using first ONLY e.g. for use with std::sort
	*/
	template<typename PairType>
	struct PairMatcherFirstElement
	: std::binary_function<PairType, PairType, bool >
	{
			bool operator() (const PairType& left, const PairType& right)const
			{
					return ( left.first == right.first ) ;
			}
	};

	/**
		@brief Struct for comparison of std::pair using second ONLY e.g. for use with std::sort
	*/
	template<typename PairType>
	struct PairMatcherSecondElement
	: std::binary_function<PairType, PairType, bool >
	{
			bool operator() (const PairType& left, const PairType& right)const
			{
					return ( left.second == right.second ) ;
			}
	};

	/**
		@brief Struct for binary predicate to consider equality with a certain tolerance

		@param i first value
		@param i second value
		@return bool if the two parameters are in tolerance close to each other
	*/
	template<typename CompareType>
	struct EqualInTolerance
	: public std::binary_function<CompareType, CompareType, bool>
	{
		CompareType& tolerance;
		EqualInTolerance( CompareType& c ) : tolerance(c) {}
		bool operator()(CompareType i, CompareType j)
		{
			CompareType diff = fabs(i-j);
			return (diff <= tolerance);
		}
	};
}

#endif // OPENMS_KERNEL_COMPARATORUTILS_H
