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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------


#ifndef OPENMS_KERNEL_COMPARATORUTILS_H
#define OPENMS_KERNEL_COMPARATORUTILS_H

#include <functional>

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
   	<code>PointerComparator<DPeak<1>::IntensityLess></code> to compare
   	<code>DPeak<1>*</code> in the same way as
   	<code>DPeak<1>::IntensityLess</code> works for <code>DPeak<1></code> .
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

	/**@brief  Make-function to create a PointerComparator from another comparator without the need to specify the template arguments.
		 \relatesalso PointerComparator

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

	/**@brief  Make-function to create a ReverseComparator from another comparator without the need to specify the template arguments.
		 \relatesalso ReverseComparator

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
	
	/**@brief  Make-function to create a LexicographicComparator from two other comparators without the need to specify the template arguments.
		 \relatesalso LexicographicComparator

	The usage is similar to pointerComparator() or reverseComparator(), which
	see.
	*/
	template < typename Cmp1, typename Cmp2 >
	LexicographicComparator < Cmp1, Cmp2 > lexicographicComparator ( Cmp1 const & cmp1, Cmp2 const & cmp2 )
	{ return LexicographicComparator < Cmp1, Cmp2 > ( cmp1, cmp2 ); }
	

}

#endif // OPENMS_KERNEL_COMPARATORUTILS_H
