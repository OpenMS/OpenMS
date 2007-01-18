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
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/FORMAT/PersistentObject.h>

#include <vector>
#include <algorithm>

namespace OpenMS
{

	/**	
		@brief Non-polymorphic peak container implemented as an array
		
		This class represents an array of D-dimensional peaks.
		It is based on the STL vector class, but provides more a more 
		convenient interface to manipulate these vectors, sort with
		respect to specific dimensions or intensity and 
		a convenient interface to the other OpenMS classes.
		
		This non-polymorphic peak array should not be used with objects 
		that are derived from DPeak and DPeaks at a time. See DPeakArray
		for a container that can handle mixed object types.
	
		@ingroup Kernel, Serialization
	*/
	template <Size D, typename PeakT = DPeak<D> >
	class DPeakArray
		:	public std::vector<PeakT>, 
			public PersistentObject
	{
		public:
		
		/// Peak type
		typedef PeakT PeakType;
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
		inline DPeakArray() : PersistentObject() {} 
		/// See std::vector documentation.
		inline DPeakArray(const DPeakArray& p) : Base(p), PersistentObject(p) {}
		/// See std::vector documentation.
		inline ~DPeakArray() {}
		/// See std::vector documentation.
		DPeakArray(typename std::vector<PeakType>::size_type n) : PersistentObject()
		{
			resize(n);
		} 
		/// See std::vector documentation.
		DPeakArray(typename std::vector<PeakType>::size_type n, const PeakType& peak) : PersistentObject()
		{
			reserve(n);
			for (typename std::vector<PeakType>::size_type i=0;i<n;i++)
			{
				push_back(peak);
			}
		} 
		/// See std::vector documentation.
		template <class InputIterator>
		DPeakArray(InputIterator f, InputIterator l) : PersistentObject()
		{
			for (InputIterator it=f;it!=l;++it)
			{
				push_back(*it);
			}
		}
		//@}
		
		/// See std::vector documentation.
		DPeakArray& operator = (const DPeakArray& rhs) 
		{ 
			if (this==&rhs) return *this;
								
			Base::resize(rhs.size());
			std::copy(rhs.begin(), rhs.end(), Base::begin());
			
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
		void sortByNthPosition(UnsignedInt i) throw (Exception::NotImplemented);
		
		/** 
			@name Generic sorting function templates.
			Any peak comparator can begiven as template argument.
			
			<p> Thus your can e.g. write <code>peaks.sortByComparator <
			DPeak<1>::IntensityLess > ()</code>, if peaks has type
			<code>DPeakArray < 1, DPeak <1> ></code>.
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

		///@name Serialization
		//@{
	 private:
		/// Serialization interface
		template<class Archive>
		void serialize(Archive & ar, const unsigned int /* version */ )
		{
			ar & boost::serialization::make_nvp("vector",boost::serialization::base_object<Base>(*this));
		}
		//@}

		/// Serialization
		friend class boost::serialization::access;
			
	};

	///Print the contents to a stream.
	template <Size D, typename Peak>
	std::ostream& operator << (std::ostream& os, const DPeakArray<D, Peak>& array)
	{
		os << "-- DPEAKARRAY-NONPOLYMORPHIC BEGIN --"<<std::endl;
		for (typename DPeakArray<D, Peak>::const_iterator it = array.begin(); it!=array.end(); ++it)
		{
			os << *it << std::endl;
		}
		os << "-- DPEAKARRAY-NONPOLYMORPHIC END --"<<std::endl;
		return os;
	}

//---------------------------------------------------------------
//  Implementation of the inline / template functions
//---------------------------------------------------------------

	template <Size D, typename PeakT > 
	void DPeakArray<D,PeakT>::sortByNthPosition(UnsignedInt i) throw (Exception::NotImplemented)
	{ 
		OPENMS_PRECONDITION(i < Index(D), "illegal dimension")
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
