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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_DPEAKARRAY_H
#define OPENMS_KERNEL_DPEAKARRAY_H

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/FORMAT/PersistentObject.h>

#include <vector>
#include <algorithm>

namespace OpenMS
{

	/**	
		@brief Peak container implemented as an array.
		
		This class represents an array of D-dimensional peaks.  The peak type must
		provide an enum DIMENSION (values 1, 2, and 3 are supported).  The
		container is based on the STL vector class, but provides more a more
		convenient interface to manipulate these vectors, sort with respect to
		specific dimensions or intensity and a convenient interface to the other
		OpenMS classes.
		
		Note that this is a non-polymorphic container, i.e. you cannot store
		objects of different types in it.

		@todo Implement clearChildIds_ (Marc)
		
		@ingroup Kernel
	*/
	template <typename PeakT>
	class DPeakArray
		:	public std::vector<PeakT>, 
			public PersistentObject
	{
		public:
		
		/// Peak type
		typedef PeakT PeakType;

		/// Dimensionality of the peaks
		/** Values 1, 2, and 3 are supported. */
		enum { DIMENSION = PeakType::DIMENSION };

		/// Base class type
		typedef std::vector<PeakType> Base;

		/// Mutable iterator		
		typedef typename std::vector<PeakType>::iterator Iterator;
		/// Non-mutable iterator
		typedef typename std::vector<PeakType>::const_iterator ConstIterator;
		/// Mutable reverse iterator
		typedef typename std::vector<PeakType>::reverse_iterator ReverseIterator;
		/// Non-mutable reverse iterator
		typedef typename std::vector<PeakType>::const_reverse_iterator ConstReverseIterator;

		/**	@name Constructors and Destructor
		*/
		//@{

		/// See std::vector documentation.
		inline DPeakArray() : Base(), PersistentObject() {} 

		/// See std::vector documentation.
		inline DPeakArray(const DPeakArray& p) : Base(p), PersistentObject(p) {}

		/// See std::vector documentation.
		inline ~DPeakArray() {}

		/// See std::vector documentation.
		DPeakArray(typename std::vector<PeakType>::size_type n) : Base(n), PersistentObject() {} 

		/// See std::vector documentation.
		DPeakArray(typename std::vector<PeakType>::size_type n, const PeakType& peak) : Base(n, peak), PersistentObject() {} 

		/// See std::vector documentation.
		template <class InputIterator>
		DPeakArray(InputIterator f, InputIterator l) : Base(f,l), PersistentObject() {}

		//@}
		
		/// See std::vector documentation.
		DPeakArray& operator = (const DPeakArray& rhs) 
		{ 
			if (this==&rhs) return *this;
			Base::operator=(rhs);
			// don't return Base immediately to avoid a cast
			return *this;
		}

		/**	
			@name Sorting.
			These simplified sorting methods are supported in addition to	
			the standard sorting methods of std::vector.
		*/
		//@{
		/// Sorts the peaks according to ascending intensity
		void sortByIntensity(bool reverse=false) 
		{ 
			if (reverse)
			{
				std::sort(Base::begin(), Base::end(), reverseComparator(typename PeakType::IntensityLess())); 
			}
			else
			{
				std::sort(Base::begin(), Base::end(), typename PeakType::IntensityLess()); 
			}
		}
		
		/// Lexicographically sorts the peaks by their position.
		void sortByPosition() 
		{ 
			std::sort(Base::begin(), Base::end(), typename PeakType::PositionLess());
		}

		/**
			@brief Sorts the peaks by one dimension of their position.
			
			It is only sorted according to dimentsion @p i . 
		*/
		void sortByNthPosition(UInt i) throw (Exception::NotImplemented);
		
		/** 
			@name Generic sorting function templates.
			Any peak comparator can begiven as template argument.
			
			<p> Thus your can e.g. write <code>peaks.sortByComparator <
			Peak1D::IntensityLess > ()</code>, if peaks has type
			<code>DPeakArray < Peak1D ></code>.
		*/
		//@{
		template < typename ComparatorType >
		void sortByComparator ( ComparatorType const & comparator )
		{ 
			std::sort(Base::begin(), Base::end(), ComparatorType( comparator ) ); 
		}
		template < typename ComparatorType >
		void sortByComparator ()
		{ 
			std::sort(Base::begin(), Base::end(), ComparatorType() ); 
		}
		//@}

		/// See std::vector documentation.
		bool operator == (const DPeakArray& array) const
		{
			return std::operator==(*this,array);
		}
		
		/// See std::vector documentation.
		bool operator !=(const DPeakArray& array) const
		{
			return !(operator==(array));
		}

		/// Comparison of container sizes
		bool operator < (const DPeakArray& array) const
		{
			return Base::size() < array.size();
		}

		/// Comparison of container sizes
		bool operator > (const DPeakArray& array) const
		{
			return Base::size() > array.size();
		}
		
		/// Comparison of container sizes
		bool operator <= (const DPeakArray& array) const
		{
			return operator<(array) || operator==(array);
		}

		/// Comparison of container sizes
		bool operator >= (const DPeakArray& array) const
		{
			return operator>(array) || operator==(array);
		}

		protected:
			// Docu in base class
	    virtual void clearChildIds_()
	    {
	    	//TODO Persistence
	    };
	};

	/// Print the contents to a stream.
	template <typename PeakT>
	std::ostream& operator << (std::ostream& os, const DPeakArray<PeakT>& array)
	{
		os << "-- DPEAKARRAY BEGIN --"<<std::endl;
		for (typename DPeakArray<PeakT>::const_iterator it = array.begin(); it!=array.end(); ++it)
		{
			os << *it << std::endl;
		}
		os << "-- DPEAKARRAY END --"<<std::endl;
		return os;
	}

//---------------------------------------------------------------
//  Implementation of the inline / template functions
//---------------------------------------------------------------

	template <typename PeakT > 
	void DPeakArray<PeakT>::sortByNthPosition(UInt i) throw (Exception::NotImplemented)
	{ 
		OPENMS_PRECONDITION(i < UInt(DIMENSION), "illegal dimension")
		if (i==0)
		{
			std::sort(Base::begin(), Base::end(), typename PeakType::template NthPositionLess<0>() );
		}
		else if (i==1)
		{
			std::sort(Base::begin(), Base::end(), typename PeakType::template NthPositionLess<1>() );
		}
		else if (i==2)
		{
			std::sort(Base::begin(), Base::end(), typename PeakType::template NthPositionLess<2>() );
		}
		else
		{
			throw Exception::NotImplemented(__FILE__,__LINE__,__FUNCTION__);
		}
	}

} // namespace OpenMS

#endif // OPENMS_KERNEL_DPEAKARRAY_H
