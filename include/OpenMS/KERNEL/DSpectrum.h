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

#ifndef OPENMS_KERNEL_DSPECTRUM_H
#define OPENMS_KERNEL_DSPECTRUM_H

#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/KERNEL/RangeManager.h>

#include <gsl/gsl_randist.h>


#include <list>

namespace OpenMS
{

	class Peak1D;

	namespace Internal
	{
		/** @brief Internal class used to store some information about
				precursor ions.

				This class is designed for limited use cases, such as storing
				precursor information from DTA files.  No data processing!
				In particular does not interact well with classes like Peak1D or Peak2D.

				@internal If you ever think about using it for more than
				the most trivial tasks, please contact the maintainer!
				We could easily replace DPeak with a better class, but
				at the moment this does not seem to pay off the effort.
				The class has been pulled out of the scope of the DSpectrum class
				because it does not depend on the container type, only on the dimension.
				Thus we can avoid unnecessary code duplication and incompatible types,
				e.g. when raw data and picked data is present during peak picking.

		*/
		template < UInt D >
		class PrecursorPeak : public DPeak<D>
		{

			/// Base class (do not even think of using this outside the scope of this class)
			typedef DPeak<D> Base;

		 public:

			/// Charge Type
			typedef Int ChargeType;

			/// Dimensionality
			enum
				{
					DIMENSION = D
				};

			/// Default constructor
			PrecursorPeak()
				: Base(),
					charge_(0)
			{
			}

			/// Copy constructor
			PrecursorPeak(const PrecursorPeak& rhs)
				: Base(rhs),
					charge_(rhs.charge_)
			{
			}

			/// Assignment operator
			PrecursorPeak & operator=(const PrecursorPeak& rhs)
			{
				Base::operator=(rhs);
				charge_=rhs.charge_;
				return *this;
			}

			/// Destructor
			~PrecursorPeak()
			{
			}

			/// Non-mutable access to the charge
			ChargeType const & getCharge() const
			{
				return charge_;
			}

			/// Mutable access to the charge
			void setCharge( ChargeType charge )
			{
				charge_ = charge;
				return;
			}

		 protected:

			ChargeType charge_;

		};

	} // namespace Internal

	/**
		 @brief Representation of a D-dimensional spectrum.

		 The peak data itself is stored in a container class, which can be a DPeakArray
		 or a STL container like std::list or std::vector.

		 Some meta information about the spectrum (ms-level, precursor peak, ...) is
		 also stored. If you want to store more meta information
		 see the MSSpectrum and MSExperiment classes.

		 The interface to the container is wrapped for convenience. Only  members
		 and types contained in both std::list and std::vector are available.

		 Additionally an interface for the minimum and maximum position, and the minimum and maximum
		 intensity of the peaks is provided by RangeManager.

		 @ingroup Kernel
	*/
	template < typename ContainerT = DPeakArray<Peak1D> >
	class DSpectrum	:
		public MetaInfoInterface,
		public RangeManager<ContainerT::value_type::DIMENSION>
	{
	 public:

		/**	@name	Type definitions */
		//@{

		/// Peak container type
		typedef ContainerT ContainerType;

		/// Peak type
		typedef typename ContainerType::value_type PeakType;

		/// Dimensionality of the peaks
		enum
			{
				DIMENSION = PeakType::DIMENSION
			};

		/// Coordinate type
		typedef typename PeakType::CoordinateType CoordinateType;

		/// Precursor peak type
		typedef Internal::PrecursorPeak<DIMENSION> PrecursorPeakType;

		/// Rangemanger type
		typedef RangeManager<DIMENSION> RangeManagerType;

		//@}

		/**	@name	STL-compliance type definitions of the container interface*/
		//@{
		typedef typename ContainerType::iterator iterator;
		typedef typename ContainerType::const_iterator const_iterator;
		typedef typename ContainerType::reverse_iterator reverse_iterator;
		typedef typename ContainerType::const_reverse_iterator const_reverse_iterator;
		typedef typename ContainerType::value_type value_type;
		typedef typename ContainerType::reference reference;
		typedef typename ContainerType::const_reference const_reference;
		typedef typename ContainerType::pointer pointer;
		typedef typename ContainerType::difference_type difference_type;
		typedef typename ContainerType::size_type size_type;
		//@}

		/**	@name	Type definitions of the container interface*/
		//@{
		/// Mutable iterator
		typedef typename ContainerType::iterator Iterator;
		/// Non-mutable iterator
		typedef typename ContainerType::const_iterator ConstIterator;
		/// Mutable reverse iterator
		typedef typename ContainerType::reverse_iterator ReverseIterator;
		/// Non-mutable reverse iterator
		typedef typename ContainerType::const_reverse_iterator ConstReverseIterator;
		//@}

		/**	@name Constructors and Destructor */
		//@{

		/// Default constructor
		DSpectrum()
			:	MetaInfoInterface(),
				RangeManagerType(),
				container_(),
				precursor_peak_(),
				retention_time_(-1), // warning: don't change this !! Otherwise MSExperimentExtern might not behave as expected !!
				ms_level_(1),
				name_()
		{
		}

		/// Copy constructor
		DSpectrum(const DSpectrum& rhs)
			: MetaInfoInterface(rhs),
				RangeManagerType(rhs),
				container_(rhs.container_),
				precursor_peak_(rhs.precursor_peak_),
				retention_time_(rhs.retention_time_),
				ms_level_(rhs.ms_level_),
				name_(rhs.name_)
		{
		}

		/// Destructor
		inline ~DSpectrum()
		{
		}
		//@}

		/// Assignment operator
		DSpectrum& operator = (const DSpectrum& rhs)
		{
			if (this==&rhs) return *this;

			MetaInfoInterface::operator=(rhs);
			RangeManagerType::operator=(rhs);
			container_ = rhs.container_;
			precursor_peak_ = rhs.precursor_peak_;
			retention_time_ = rhs.retention_time_;
			ms_level_ = rhs.ms_level_;
			name_ = rhs.name_;
			return *this;
		}

		/// Equality operator
		bool operator == (const DSpectrum& rhs) const
		{
			return
				MetaInfoInterface::operator==(rhs) &&
				RangeManagerType::operator==(rhs) &&
				container_ == rhs.container_ &&
				precursor_peak_ == rhs.precursor_peak_ &&
				retention_time_ == rhs.retention_time_ &&
				ms_level_ == rhs.ms_level_
				;
			//name_ == rhs.name_  // the name can differ => do not test it
		}

		/// Equality operator
		bool operator != (const DSpectrum& rhs) const
		{
			return !(operator==(rhs));
		}

		/**	@name	Wrappers of container accessors */
		//@{
		/// Non-mutable access to the peak container
		inline const ContainerType& getContainer() const
		{
			return container_;
		}
		/// Mutable access to the peak container.
		inline ContainerType& getContainer()
		{
			return container_;
		}
		/// Mutable access to the peak container.
		inline void setContainer(const ContainerType& container)
		{
			container_ = container;
		}

		/// Returns the const begin iterator of the container
		inline ConstIterator begin() const
		{
			return container_.begin();
		}
		/// Returns the const end iterator of the container
		inline ConstIterator end() const
		{
			return container_.end();
		}

		/// Returns the begin iterator of the container
		inline Iterator begin()
		{
			return container_.begin();
		}
		/// Returns the end iterator of the container
		inline Iterator end()
		{
			return container_.end();
		}

		/// returns the element with index n
		reference operator[] (size_type n)
		{
			return container_[n];
		}

		/// returns the element with index n
		const_reference operator[] (size_type n) const
		{
			return container_[n];
		}

		/// returns the maxium size possbile (the number of peaks)
		inline size_type max_size() const
		{
			return container_.max_size();
		}
		/// returns the size (the number of peaks)
		inline UInt size() const
		{
			return container_.size();
		}
		/// Returns if the container is empty
		inline bool empty() const
		{
			return container_.empty();
		}

		/// Swaps two containers
		inline void swap(ContainerType& rhs)
		{
			container_.swap(rhs);
		}

		/// Comparison of container sizes
		inline bool operator<(const DSpectrum& rhs)
		{
			return container_<rhs.getContainer();
		}

		/// Comparison of container sizes
		inline bool operator>(const DSpectrum& rhs)
		{
			return container_>rhs.getContainer();
		}

		/// Comparison of container sizes
		inline bool operator<=(const DSpectrum& rhs)
		{
			return container_<=rhs.getContainer();
		}

		/// Comparison of container sizes
		inline bool operator>=(const DSpectrum& rhs)
		{
			return container_>=rhs.getContainer();
		}

		/// See STL documentation
		inline ReverseIterator rbegin()
		{
			return container_.rbegin();
		}

		/// See STL documentation
		inline ConstReverseIterator rbegin() const
		{
			return container_.rbegin();
		}

		/// See STL documentation
		inline ReverseIterator rend()
		{
			return container_.rend();
		}

		/// See STL documentation
		inline ConstReverseIterator rend() const
		{
			return container_.rend();
		}

		/// Inserts an element
		inline Iterator insert( Iterator loc, const value_type& val )
		{
			return container_.insert(loc, val);
		}

		/// Inserts an element several times
		inline void insert( iterator loc, size_type num, const value_type& val )
		{
			container_.insert(loc, num, val);
		}

		/// Inserts a range of elements
		template<class InputIterator> void insert( iterator loc, InputIterator start, InputIterator end )
		{
			container_.insert(loc, start, end);
		}

		/// Erases an element
		inline Iterator erase( iterator loc )
		{
			return container_.erase(loc);
		}

		/// Erases a range of elements
		inline Iterator erase( iterator start, iterator end )
		{
			return container_.erase(start, end);
		}

		/// Returns the first element
		inline value_type& front()
		{
			return container_.front();
		}

		/// Returns the first element
		inline const value_type& front() const
		{
			return container_.front();
		}

		/// Returns the last element
		inline value_type& back()
		{
			return container_.back();
		}

		/// Returns the last element
		inline const value_type& back() const
		{
			return container_.back();
		}

		/// Removes the last element
		inline void pop_back()
		{
			container_.pop_back();
		}

		/// Inserts an element at the end
		inline void push_back( const value_type& val )
		{
			container_.push_back(val);
		}

		/// Fills the container with serval copies of a value
		inline void assign( size_type num, const value_type& val )
		{
			container_.assign(num, val);
		}

		/// Fills the container with a range of values
		template<class InputIterator> void assign( InputIterator start, InputIterator end )
		{
			container_.assign(start, end);
		}

		/// Removes all elements
		inline void clear()
		{
			container_.clear();
		}

		/// Resizes the container to size @p num. Uses @p val to fill up if it is shorter than @p num.
		inline void resize( size_type num, const value_type& val = value_type() )
		{
			container_.resize(num, val);
		}

		//@}

		// Docu in base class
		virtual void updateRanges()
		{
			this->clearRanges();
			updateRanges_(container_.begin(), container_.end());
		}

		/**	@name Accessors for meta information*/
		//@{
		/// const accessor for the precorsor peak
		const PrecursorPeakType& getPrecursorPeak() const
		{
			return precursor_peak_;
		}

		/// accessor for the precorsor peak
		PrecursorPeakType& getPrecursorPeak()
		{
			return precursor_peak_;
		}

		/// sets the precursor peak
		void setPrecursorPeak(const PrecursorPeakType& peak)
		{
			precursor_peak_ = peak;
		}

		/// returns the absolute retention time (unit is seconds)
		CoordinateType getRT() const
		{
			return retention_time_;
		}


		/**
			 Sets the retention time and the start/stop time of the gradient.
			 The latter two are needed for calculating the normalized retention time
		*/
		void setRT(CoordinateType rt)
		{
			retention_time_= rt;
		}

		/**
			 @brief Returns the MS level.

			 For survey scans this is 1, for MS/MS scans 2, ...
		*/
		UInt getMSLevel() const
		{
			return ms_level_;
		}

		///Sets the MS level.
		void setMSLevel(UInt ms_level)
		{
			ms_level_ = ms_level;
		}

		///Returns the name
		String getName() const
		{
			return name_;
		}

		///Sets the name
		void setName(const String& name)
		{
			name_ = name;
		}

		//@}

		/**
			 @brief Fast search for peak range begin

			 @note Make sure the spectrum is sorted with respect to m/z ratio! Otherwise the result is undefined.
		*/
		Iterator MZBegin(double mz)
		{
			PeakType p;
			p.setPosition(mz);
			return lower_bound(container_.begin(), container_.end(), p, typename PeakType::PositionLess());
		}

		/**
			 @brief Fast search for peak range end (returns the past-the-end iterator)

			 @note Make sure the spectrum is sorted with respect to m/z ratio. Otherwise the result is undefined.
		*/
		Iterator MZEnd(double mz)
		{
			PeakType p;
			p.setPosition(mz);
			return upper_bound(container_.begin(), container_.end(), p, typename PeakType::PositionLess());
		}

		/**
			 @brief Fast search for peak range begin

			 @note Make sure the spectrum is sorted with respect to m/z ratio! Otherwise the result is undefined.
		*/
		ConstIterator MZBegin(double mz) const
		{
			PeakType p;
			p.setPosition(mz);
			return lower_bound(container_.begin(), container_.end(), p, typename PeakType::PositionLess());
		}

		/**
			 @brief Fast search for peak range end (returns the past-the-end iterator)

			 @note Make sure the spectrum is sorted with respect to m/z ratio. Otherwise the result is undefined.
		*/
		ConstIterator MZEnd(double mz) const
		{
			PeakType p;
			p.setPosition(mz);
			return upper_bound(container_.begin(), container_.end(), p, typename PeakType::PositionLess());
		}

	 protected:

		/// The container with all the peak data
		ContainerType		container_;

		/// Precursor information
		PrecursorPeakType precursor_peak_;

		/// retention time
		CoordinateType retention_time_;

		/// MS level
		UInt ms_level_;

		/// Name
		String name_;
	};

	///Print the contents to a stream.
	template <typename Container>
	std::ostream& operator << (std::ostream& os, const DSpectrum<Container>& rhs)
	{
		os << "-- DSpectrum BEGIN --"<<std::endl;
		os << "MS-LEVEL:" <<rhs.getMSLevel() << std::endl;
		os << "RT:" <<rhs.getRT() << std::endl;
		os << "NAME:" <<rhs.getName() << std::endl;
		os << "-- DSpectrum END --"<<std::endl;

		return os;
	}

} // namespace OpenMS



#endif // OPENMS_KERNEL_SPECTRUM_H
