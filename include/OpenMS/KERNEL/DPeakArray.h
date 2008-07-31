// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

		@ingroup Kernel
	*/
	template <typename PeakT, typename AllocT = std::allocator<PeakT> >
	class DPeakArray
		:	public std::vector<PeakT, AllocT>, 
			public PersistentObject
	{
		public:
		
		/// Peak type
		typedef PeakT PeakType;

    /// Allocator type
    typedef AllocT AllocType;
    
		/// Dimensionality of the peaks
		/** Values 1, 2, and 3 are supported. */
		enum { DIMENSION = PeakType::DIMENSION };

		/// Base class type
		typedef std::vector<PeakType, AllocT> Base;

		/// Mutable iterator		
		typedef typename std::vector<PeakType, AllocT>::iterator Iterator;
		/// Non-mutable iterator
		typedef typename std::vector<PeakType, AllocT>::const_iterator ConstIterator;
		/// Mutable reverse iterator
		typedef typename std::vector<PeakType, AllocT>::reverse_iterator ReverseIterator;
		/// Non-mutable reverse iterator
		typedef typename std::vector<PeakType, AllocT>::const_reverse_iterator ConstReverseIterator;

    // allow c'tors of other template instances to access private members
    template < typename ContainerT_, typename AllocT_ > friend class DPeakArray;    
    
		/**	@name Constructors and Destructor
		*/
		//@{

		/// See std::vector documentation.
		inline DPeakArray() : Base(), PersistentObject() {} 

		/// See std::vector documentation.
		inline DPeakArray(const DPeakArray& p) : Base(p), PersistentObject(p) {}

    /// See std::vector documentation. (different allocator)
    template <typename AllocT2>
    inline DPeakArray(const DPeakArray<PeakT, AllocT2>& p) : Base(p.begin(),p.end()), PersistentObject(p) {}

    /// Constructor, using a given DPeakArray (with other alloc?) @p p and the allocator instance that shall be used
    template <typename AllocT2>
    inline DPeakArray(const DPeakArray<PeakT, AllocT2>& p, const AllocT& alloc) 
    : 
      Base(alloc),
      PersistentObject(p)
    {
      // copy data manually into new vector
      //TODO #ifdef OPENMS_64BIT_ARCHITECTURE .. use large capacity for ExternalAlloc() to avoid reallocation?! Check cache locality (benchmark)!
      this->reserve(p.capacity());
      //std::copy(p.begin(), p.end(), std::back_inserter(this->back()));        
      for(typename std::vector<PeakT, AllocT2>::const_iterator it = p.begin(); it != p.end(); ++it)
      {
        this->push_back(*it);
      }
      
      if (p.size() != this->size())
      {
        std::cout << "ERROR: given size: " << p.size() << " pushed size: " << this->size() << std::endl;
      }
    }

    
		/// See std::vector documentation.
		inline ~DPeakArray() {}

		/// See std::vector documentation.
		DPeakArray(typename std::vector<PeakType, AllocT>::size_type n) : Base(n), PersistentObject() {} 

		/// See std::vector documentation.
		DPeakArray(typename std::vector<PeakType, AllocT>::size_type n, const PeakType& peak) : Base(n, peak), PersistentObject() {} 


    /// See std::vector documentation.
    DPeakArray(const AllocT& a) : Base(a), PersistentObject() {} 

    /// See std::vector documentation.
    DPeakArray(typename std::vector<PeakType>::size_type n, const PeakType& peak, const AllocT& a) : Base(n, peak, a), PersistentObject() {} 
    
        
    
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
			
			It is only sorted according to dimentsion @p i.
			
			@exception Exception::NotImplemented is thrown if sorting for dimension @p i is not implemented
		*/
		void sortByNthPosition(UInt i);
		
		/** 
			@name Generic sorting function templates.
			Any peak comparator can be given as template argument.
			
			<p> Thus your can e.g. write <code>peaks.sortByComparator<Peak1D::IntensityLess>()</code>, if peaks has type
			<code>DPeakArray <></code>.
		*/
		//@{
		template < typename ComparatorType >
		void sortByComparator ( ComparatorType const & comparator  = ComparatorType() )
		{ 
			std::sort(Base::begin(), Base::end(), ComparatorType( comparator ) ); 
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
		template <typename AllocT2>
		bool operator < (const DPeakArray<PeakT, AllocT2>& array) const
		{
			return Base::size() < array.size();
		}

		/// Comparison of container sizes
		template <typename AllocT2>
		bool operator > (const DPeakArray<PeakT, AllocT2>& array) const
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
				//nothing to do here
	    };
	};

	/// Print the contents to a stream.
	template <typename PeakT, typename AllocT>
	std::ostream& operator << (std::ostream& os, const DPeakArray<PeakT, AllocT>& array)
	{
		os << "-- DPEAKARRAY BEGIN --"<<std::endl;
		for (typename DPeakArray<PeakT, AllocT>::const_iterator it = array.begin(); it!=array.end(); ++it)
		{
			os << *it << std::endl;
		}
		os << "-- DPEAKARRAY END --"<<std::endl;
		return os;
	}

//---------------------------------------------------------------
//  Implementation of the inline / template functions
//---------------------------------------------------------------

	template <typename PeakT, typename Alloc > 
	void DPeakArray<PeakT,Alloc>::sortByNthPosition(UInt i)
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
