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

#ifndef OPENMS_KERNEL_DSPECTRUM_H
#define OPENMS_KERNEL_DSPECTRUM_H

#include <OpenMS/KERNEL/DRichPeak.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>

#include <gsl/gsl_randist.h>

#include <cmath>

namespace OpenMS
{

	class Peak1D;

	namespace Internal
	{
		/**
			@brief Internal class used to store some information about
			precursor ions.

			This class is designed for limited use cases, such as storing
			precursor information from DTA files. No data processing!
		*/
		template < UInt D >
		class OPENMS_DLLAPI PrecursorPeak
			: public DRichPeak<D>::Type
		{

			/// Base class (do not even think of using this outside the scope of this class)
			typedef typename DRichPeak<D>::Type Base;

		 public:

			/// Dimensionality
			enum
				{
					DIMENSION = D
				};

			/// Default constructor
			PrecursorPeak()
				: Base(),
					charge_(0),
					possible_charge_states_()
			{
			}

			/// Copy constructor
			PrecursorPeak(const PrecursorPeak& rhs)
				: Base(rhs),
					charge_(rhs.charge_),
					possible_charge_states_(rhs.possible_charge_states_)
			{
			}

			/// Assignment operator
			PrecursorPeak & operator=(const PrecursorPeak& rhs)
			{
				Base::operator=(rhs);

				charge_=rhs.charge_;
				possible_charge_states_ = rhs.possible_charge_states_;

				return *this;
			}

			/// Destructor
			~PrecursorPeak()
			{
			}

			/// Non-mutable access to the charge
			Int const & getCharge() const
			{
				return charge_;
			}

			/// Mutable access to the charge
			void setCharge( Int charge )
			{
				charge_ = charge;
				return;
			}

			std::vector<Int>& getPossibleChargeStates()
			{
				return possible_charge_states_;
			}

			const std::vector<Int>& getPossibleChargeStates() const
			{
				return possible_charge_states_;
			}

			void setPossibleChargeStates(const std::vector<Int>& possible_charge_states)
			{
				possible_charge_states_ = possible_charge_states;
			}

		 protected:

			Int charge_;
			std::vector<Int> possible_charge_states_;

		};

	} // namespace Internal

	/**
		@brief Representation of a D-dimensional spectrum.

		Some meta information about the spectrum (ms-level, precursor peak, ...) is
		also stored. If you want to store more meta information
		see the MSSpectrum and MSExperiment classes.

		Additionally an interface for the minimum and maximum position, and the minimum and maximum
		intensity of the peaks is provided by RangeManager.

		@ingroup Kernel
	*/
	template < typename PeakT = Peak1D, typename AllocT = std::allocator<PeakT> >
	class OPENMS_DLLAPI DSpectrum
		: public std::vector<PeakT, AllocT>,
			public MetaInfoInterface,
			public RangeManager<PeakT::DIMENSION>
	{
	 public:

	  ///Meta data array struct containing meta information and a name
		class OPENMS_DLLAPI MetaDataArray
	    : public MetaInfoDescription,
	    	public std::vector<Real>
	  {
	  };


		/**	@name	Type definitions */
		//@{
		/// Peak type
		typedef PeakT PeakType;
		/// Peak container type
		typedef std::vector<PeakType, AllocT> ContainerType;
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
		/// MetaDataArrays type
		typedef std::vector<MetaDataArray> MetaDataArrays;
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
			:	ContainerType(),
      	MetaInfoInterface(),
				RangeManagerType(),
				precursor_peak_(),
				retention_time_(-1), // warning: don't change this !! Otherwise MSExperimentExtern might not behave as expected !!
				ms_level_(1),
				name_(),
				meta_data_arrays_()
		{
		}

    /// constructor with custom allocator
    DSpectrum(const AllocT& alloc)
      : ContainerType(alloc),
      	MetaInfoInterface(),
        RangeManagerType(),
        precursor_peak_(),
        retention_time_(-1), // warning: don't change this !! Otherwise MSExperimentExtern might not behave as expected !!
        ms_level_(1),
        name_(),
        meta_data_arrays_()
    {
    }

		/// Copy constructor
		DSpectrum(const DSpectrum& rhs)
			: ContainerType(rhs),
      	MetaInfoInterface(rhs),
				RangeManagerType(rhs),
				precursor_peak_(rhs.precursor_peak_),
				retention_time_(rhs.retention_time_),
				ms_level_(rhs.ms_level_),
				name_(rhs.name_),
				meta_data_arrays_(rhs.meta_data_arrays_)
		{
		}

    /// Copy constructor for different allocator
    template < typename AllocT2>
    DSpectrum(const DSpectrum<PeakType,AllocT2>& rhs)
      : ContainerType(rhs),
      	MetaInfoInterface(rhs),
        RangeManagerType(rhs),
        precursor_peak_(rhs.precursor_peak_),
        retention_time_(rhs.retention_time_),
        ms_level_(rhs.ms_level_),
        name_(rhs.name_),
        meta_data_arrays_(rhs.meta_data_arrays_)
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

			ContainerType::operator=(rhs);
			MetaInfoInterface::operator=(rhs);
			RangeManagerType::operator=(rhs);
			precursor_peak_ = rhs.precursor_peak_;
			retention_time_ = rhs.retention_time_;
			ms_level_ = rhs.ms_level_;
			name_ = rhs.name_;
			meta_data_arrays_ = rhs.meta_data_arrays_;

			return *this;
		}

    /// Assignment operator for different allocator
    template < typename AllocT2>
    DSpectrum& operator = (const DSpectrum< PeakType, AllocT2 >& rhs)
    {
      if (this==&rhs) return *this;

			ContainerType::operator=(rhs);
      MetaInfoInterface::operator=(rhs);
      RangeManagerType::operator=(rhs);
      precursor_peak_ = rhs.precursor_peak_;
      retention_time_ = rhs.retention_time_;
      ms_level_ = rhs.ms_level_;
      name_ = rhs.name_;
      meta_data_arrays_ = rhs.meta_data_arrays_;

      return *this;
    }

		/// Equality operator
		bool operator == (const DSpectrum& rhs) const
		{
			return
				std::operator==(*this, rhs) &&
				MetaInfoInterface::operator==(rhs) &&
				RangeManagerType::operator==(rhs) &&
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

		// Docu in base class
		virtual void updateRanges()
		{
			this->clearRanges();
			updateRanges_(ContainerType::begin(), ContainerType::end());
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

		///Sets the absolute retention time (unit is seconds)
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
			@brief Lexicographically sorts the peaks by their intensity.

			Sorts the peaks according to ascending intensity. Meta data arrays will be sorted accordingly.
		*/
		void sortByIntensity(bool reverse=false)
		{
			if(meta_data_arrays_.size() == 0)
			{
				if (reverse)
				{
					std::sort(ContainerType::begin(), ContainerType::end(), reverseComparator(typename PeakType::IntensityLess()));
				}
				else
				{
					std::sort(ContainerType::begin(), ContainerType::end(), typename PeakType::IntensityLess());
				}
			}
			else
			{
				//sort index list
				std::vector< std::pair<typename PeakType::IntensityType,UInt> > sorted_indices;
				sorted_indices.reserve(ContainerType::size());
				for(UInt i(0); i < ContainerType::size(); ++i)
				{
					sorted_indices.push_back(std::make_pair(ContainerType::operator[](i).getIntensity(),i));
				}

				if (reverse)
				{
					std::sort(sorted_indices.begin(), sorted_indices.end(), reverseComparator(PairComparatorFirstElement< std::pair<typename PeakType::IntensityType,UInt> >()));
				}
				else
				{
					std::sort(sorted_indices.begin(), sorted_indices.end(), PairComparatorFirstElement< std::pair<typename PeakType::IntensityType,UInt> >());
				}

				//apply sorting to ContainerType and to metadataarrays
				ContainerType tmp;
				for(UInt i(0); i < sorted_indices.size(); ++i)
				{
					tmp.push_back(*(ContainerType::begin()+(sorted_indices[i].second)));
				}
				ContainerType::swap(tmp);

				for(UInt i(0); i < meta_data_arrays_.size(); ++i)
				{
					std::vector<Real> mda_tmp;
					for(UInt j(0); j < meta_data_arrays_[i].size(); ++j)
					{
						mda_tmp.push_back(*(meta_data_arrays_[i].begin()+(sorted_indices[j].second)));
					}
					meta_data_arrays_[i].swap(mda_tmp);
				}
			}

		}

		/**
			@brief Lexicographically sorts the peaks by their position.

			The spectrum is sorted with respect to position. Meta data arrays will be sorted accordingly.
		*/
		void sortByPosition()
		{
			if(meta_data_arrays_.size() == 0)
			{
				std::sort(ContainerType::begin(), ContainerType::end(), typename PeakType::PositionLess());
			}
			else
			{
				//sort index list
				std::vector< std::pair<typename PeakType::PositionType,UInt> > sorted_indices;
				sorted_indices.reserve(ContainerType::size());
				for(UInt i(0); i < ContainerType::size(); ++i)
				{
					sorted_indices.push_back(std::make_pair(ContainerType::operator[](i).getPosition(),i));
				}
				std::sort(sorted_indices.begin(), sorted_indices.end(), PairComparatorFirstElement< std::pair<typename PeakType::PositionType,UInt> >());

				//apply sorting to ContainerType and to metadataarrays
				ContainerType tmp;
				for(UInt i(0); i < sorted_indices.size(); ++i)
				{
					tmp.push_back(*(ContainerType::begin()+(sorted_indices[i].second)));
				}
				ContainerType::swap(tmp);

				//what if metadataarray is unequal in size to the spectrum?? throw error?
				for(UInt i(0); i < meta_data_arrays_.size(); ++i)
				{
					std::vector<Real> mda_tmp;
					for(UInt j(0); j < meta_data_arrays_[i].size(); ++j)
					{
						mda_tmp.push_back(*(meta_data_arrays_[i].begin()+(sorted_indices[j].second)));
					}
					std::swap(meta_data_arrays_[i],mda_tmp);
				}
			}

		}


		///@name Searching a peak or peak range
		//@{
		/**
			@brief Binary search for the peak nearest to a specific m/z

			@param mz The searched for mass-to-charge ratio searched
			@return Returns the index of the peak.

			@note Make sure the spectrum is sorted with respect to m/z ratio! Otherwise the result is undefined.

			@exception Exception::Precondition is thrown if the spectrum is empty (not only in debug mode)
		*/
		UInt findNearest(CoordinateType mz) const
		{
			//no peak => no search
			if (ContainerType::size()==0) throw Exception::Precondition(__FILE__,__LINE__,__PRETTY_FUNCTION__,"There must be at least one peak to determine the nearest peak!");

			//searh for position for inserting
			ConstIterator it = MZBegin(mz);
			//border cases
			if (it==ContainerType::begin()) return 0;
			if (it==ContainerType::end()) return ContainerType::size()-1;
			//the peak before or the current peak are closest
			ConstIterator it2 = it;
			--it2;
			if (std::fabs(it->getMZ()-mz)<std::fabs(it2->getMZ()-mz))
			{
				return it - ContainerType::begin();
			}
			else
			{
				return it2 - ContainerType::begin();
			}
		}
		/**
			 @brief Binary search for peak range begin

			 @note Make sure the spectrum is sorted with respect to m/z ratio! Otherwise the result is undefined.
		*/
		Iterator MZBegin(CoordinateType mz)
		{
			PeakType p;
			p.setPosition(mz);
			return lower_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
		}
		/**
			 @brief Binary search for peak range end (returns the past-the-end iterator)

			 @note Make sure the spectrum is sorted with respect to m/z ratio. Otherwise the result is undefined.
		*/
		Iterator MZEnd(CoordinateType mz)
		{
			PeakType p;
			p.setPosition(mz);
			return upper_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
		}
		/**
			 @brief Binary search for peak range begin

			 @note Make sure the spectrum is sorted with respect to m/z ratio! Otherwise the result is undefined.
		*/
		ConstIterator MZBegin(CoordinateType mz) const
		{
			PeakType p;
			p.setPosition(mz);
			return lower_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
		}
		/**
			 @brief Binary search for peak range end (returns the past-the-end iterator)

			 @note Make sure the spectrum is sorted with respect to m/z ratio. Otherwise the result is undefined.
		*/
		ConstIterator MZEnd(CoordinateType mz) const
		{
			PeakType p;
			p.setPosition(mz);
			return upper_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
		}
		//@}


		/**
			@name Peak meta data array methods

			These methods are used to annotate each peak in a spectrum with meta information.
			It is an intermediate way between storing the information in the peak's MetaInfoInterface
			and deriving a new peak type with members for this information.

		  These statements should help you chose which approach to use
		  - Access to meta info arrays is slower than to a member variable
		  - Access to meta info arrays is faster than to a %MetaInfoInterface
		  - Meta info arrays are stored when using mzData or mzML format for storing
		*/
		//@{
		/// Returns a const reference to the integer meta arrays
		inline const MetaDataArrays& getMetaDataArrays() const
		{
			return meta_data_arrays_;
		}
		/// Returns a mutable reference to the integer meta arrays
		inline MetaDataArrays& getMetaDataArrays()
		{
			return meta_data_arrays_;
		}
		//@}

	protected:

		/// Precursor information
		PrecursorPeakType precursor_peak_;

		/// retention time
		CoordinateType retention_time_;

		/// MS level
		UInt ms_level_;

		/// Name
		String name_;

		///Meta info arrays
		MetaDataArrays meta_data_arrays_;
	};

	///Print the contents to a stream.
	template <typename PeakT, typename AllocT>
	std::ostream& operator << (std::ostream& os, const DSpectrum<PeakT,AllocT>& rhs)
	{
		os << "-- DSpectrum BEGIN --"<<std::endl;
		os << "MS-LEVEL:" <<rhs.getMSLevel() << std::endl;
		os << "RT:" <<rhs.getRT() << std::endl;
		os << "NAME:" <<rhs.getName() << std::endl;
		for (typename DSpectrum<PeakT, AllocT>::const_iterator it = rhs.begin(); it!=rhs.end(); ++it)
		{
			os << *it << std::endl;
		}
		os << "-- DSpectrum END --"<<std::endl;

		return os;
	}

} // namespace OpenMS

#endif // OPENMS_KERNEL_DSPECTRUM_H
