// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//									 OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//	Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//	This library is free software; you can redistribute it and/or
//	modify it under the terms of the GNU Lesser General Public
//	License as published by the Free Software Foundation; either
//	version 2.1 of the License, or (at your option) any later version.
//
//	This library is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
//	Lesser General Public License for more details.
//
//	You should have received a copy of the GNU Lesser General Public
//	License along with this library; if not, write to the Free Software
//	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA	02111-1307	USA
//
// --------------------------------------------------------------------------
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MSSPECTRUM_H
#define OPENMS_KERNEL_MSSPECTRUM_H

#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <OpenMS/FORMAT/DB/PersistentObject.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>

namespace OpenMS
{
	class Peak1D;
	
	/**
		@brief The representation of a 1D spectrum.
		
		It contains peak data and metadata about specific instrument settings,
		aquisition settings, description of the meta values used in the peaks and precursor info (SpectrumSettings).
		
		Several MSSpectrum instances are contained in a peak map (MSExperiment), which is essentially
		a vector of spectra with additional information about the experiment.
		
		Precursor info from SpectrumSettings should only be used if this spectrum is a tandem-MS spectrum.
		The precursor spectrum is the first spectrum in MSExperiment, that has a lower MS-level than the current spectrum.
		
		@note For range operations, see \ref RangeUtils "RangeUtils module"!
		
		@ingroup Kernel
	*/
	template <typename PeakT = Peak1D>
	class MSSpectrum
		: public std::vector<PeakT>,
			public RangeManager<1>,
			public SpectrumSettings,
			public PersistentObject
	{

		public:

		  ///Meta data array class
			class OPENMS_DLLAPI MetaDataArray
		    : public MetaInfoDescription,
		    	public std::vector<Real>
		  {
		  };

			///Comparator for the retention time.
			struct RTLess
				: public std::binary_function <MSSpectrum, MSSpectrum, bool>
			{
				inline bool operator () (const MSSpectrum& a, const MSSpectrum& b) const
				{
					return (a.getRT() < b.getRT());
				}
			};
			
			///@name Base type definitions
			//@{
			/// Peak type
			typedef PeakT PeakType;
			/// Coordinate (m/z) type
			typedef typename PeakType::CoordinateType CoordinateType;
			/// Spectrum base type
			typedef std::vector<PeakType> ContainerType;
			/// Metadata array vector type
			typedef std::vector<MetaDataArray> MetaDataArrays;
			//@}

			///@name Peak container iterator type definitions
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


			/// Constructor
			MSSpectrum():
				ContainerType(),
				RangeManager<1>(),
				SpectrumSettings(),
				PersistentObject(),
				retention_time_(-1), // warning: don't change this !! Otherwise MSExperimentExtern might not behave as expected !!
				ms_level_(1),
				name_(),
				meta_data_arrays_()
			{
			}
			
	    /// Copy constructor
			MSSpectrum(const MSSpectrum& source):
				ContainerType(source),
				RangeManager<1>(source),
				SpectrumSettings(source),
				PersistentObject(source),
				retention_time_(source.retention_time_),
				ms_level_(source.ms_level_),
				name_(source.name_),
				meta_data_arrays_(source.meta_data_arrays_)
			{
			}
	    
	 		/// Destructor
			~MSSpectrum()
			{
			}
	
			/// Assignment operator
			MSSpectrum& operator= (const MSSpectrum& source)
			{
				if (&source == this) return *this;
	
				ContainerType::operator=(source);
				RangeManager<1>::operator=(source);
				SpectrumSettings::operator=(source);
				PersistentObject::operator=(source);
				
				retention_time_ = source.retention_time_;
				ms_level_ = source.ms_level_;
				name_ = source.name_;
				meta_data_arrays_ = source.meta_data_arrays_;
				
				return *this;
			}
	
	    
			/// Equality operator
			bool operator== (const MSSpectrum& rhs) const
			{
				return
					std::operator==(*this, rhs) &&
					RangeManager<1>::operator==(rhs) &&
					SpectrumSettings::operator==(rhs)  &&
					
					retention_time_ == rhs.retention_time_ &&
					ms_level_ == rhs.ms_level_ &&
					meta_data_arrays_ == rhs.meta_data_arrays_
					;
					//name_ can differ => it is not checked
			}
			/// Equality operator
			bool operator!= (const MSSpectrum& rhs) const
			{
				return !(operator==(rhs));
			}
			
			// Docu in base class (RangeManager)
			virtual void updateRanges()
			{
				this->clearRanges();
				updateRanges_(ContainerType::begin(), ContainerType::end());
			}
	
			/**	@name Accessors for meta information*/
			//@{
			/// returns the absolute retention time (is seconds)
			inline DoubleReal getRT() const
			{
				return retention_time_;
			}
			///Sets the absolute retention time (is seconds)
			inline void setRT(DoubleReal rt)
			{
				retention_time_= rt;
			}
			/**
				 @brief Returns the MS level.
	
				 For survey scans this is 1, for MS/MS scans 2, ...
			*/
			inline UInt getMSLevel() const
			{
				return ms_level_;
			}
			///Sets the MS level.
			inline void setMSLevel(UInt ms_level)
			{
				ms_level_ = ms_level;
			}
			///Returns the name
			inline const String& getName() const
			{
				return name_;
			}
			///Sets the name
			inline void setName(const String& name)
			{
				name_ = name;
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

			///@name Sorting peaks
			//@{
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
					std::vector< std::pair<typename PeakType::IntensityType,Size> > sorted_indices;
					sorted_indices.reserve(ContainerType::size());
					for (Size i(0); i < ContainerType::size(); ++i)
					{
						sorted_indices.push_back(std::make_pair(ContainerType::operator[](i).getIntensity(),i));
					}
	
					if (reverse)
					{
						std::sort(sorted_indices.begin(), sorted_indices.end(), reverseComparator(PairComparatorFirstElement< std::pair<typename PeakType::IntensityType,Size> >()));
					}
					else
					{
						std::sort(sorted_indices.begin(), sorted_indices.end(), PairComparatorFirstElement< std::pair<typename PeakType::IntensityType,Size> >());
					}
	
					//apply sorting to ContainerType and to meta data arrays
					ContainerType tmp;
					for (Size i(0); i < sorted_indices.size(); ++i)
					{
						tmp.push_back(*(ContainerType::begin()+(sorted_indices[i].second)));
					}
					ContainerType::swap(tmp);
	
					for (Size i(0); i < meta_data_arrays_.size(); ++i)
					{
						std::vector<Real> mda_tmp;
						for (Size j(0); j < meta_data_arrays_[i].size(); ++j)
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
					std::vector< std::pair<typename PeakType::PositionType,Size> > sorted_indices;
					sorted_indices.reserve(ContainerType::size());
					for (Size i(0); i < ContainerType::size(); ++i)
					{
						sorted_indices.push_back(std::make_pair(ContainerType::operator[](i).getPosition(),i));
					}
					std::sort(sorted_indices.begin(), sorted_indices.end(), PairComparatorFirstElement< std::pair<typename PeakType::PositionType,Size> >());
	
					//apply sorting to ContainerType and to metadataarrays
					ContainerType tmp;
					for (Size i(0); i < sorted_indices.size(); ++i)
					{
						tmp.push_back(*(ContainerType::begin()+(sorted_indices[i].second)));
					}
					ContainerType::swap(tmp);
	
					for (Size i(0); i < meta_data_arrays_.size(); ++i)
					{
						std::vector<Real> mda_tmp;
						for (Size j(0); j < meta_data_arrays_[i].size(); ++j)
						{
							mda_tmp.push_back(*(meta_data_arrays_[i].begin()+(sorted_indices[j].second)));
						}
						std::swap(meta_data_arrays_[i],mda_tmp);
					}
				}
			}
			///Checks if all peaks are sorted with respect to ascending m/z
			bool isSorted() const
			{
				for (Size i=1; i<this->size(); ++i)
				{
					if (this->operator[](i-1).getMZ()>this->operator[](i).getMZ()) return false;
				}
				return true;
			}
			//@}

			///@name Searching a peak or peak range
			//@{
			/**
				@brief Binary search for the peak nearest to a specific m/z
	
				@param mz The searched for mass-to-charge ratio searched
				@return Returns the index of the peak.
	
				@note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
	
				@exception Exception::Precondition is thrown if the spectrum is empty (not only in debug mode)
			*/
			Size findNearest(CoordinateType mz) const
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
					return Size(it - ContainerType::begin());
				}
				else
				{
					return Size(it2 - ContainerType::begin());
				}
			}
			/**
				 @brief Binary search for peak range begin
	
				 @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
			*/
			Iterator MZBegin(CoordinateType mz)
			{
				PeakType p;
				p.setPosition(mz);
				return lower_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
			}
			/**
				 @brief Binary search for peak range begin
	
				 @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
			*/
			Iterator MZBegin(Iterator begin, CoordinateType mz, Iterator end)
			{
				PeakType p;
				p.setPosition(mz);
				return lower_bound(begin, end, p, typename PeakType::PositionLess());
			}
			/**
				 @brief Binary search for peak range end (returns the past-the-end iterator)
	
				 @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
			*/
			Iterator MZEnd(CoordinateType mz)
			{
				PeakType p;
				p.setPosition(mz);
				return upper_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
			}		
			/**
				 @brief Binary search for peak range end (returns the past-the-end iterator)
	
				 @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
			*/
			Iterator MZEnd(Iterator begin, CoordinateType mz, Iterator end)
			{
				PeakType p;
				p.setPosition(mz);
				return upper_bound(begin, end, p, typename PeakType::PositionLess());
			}
	
			/**
				 @brief Binary search for peak range begin
	
				 @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
			*/
			ConstIterator MZBegin(CoordinateType mz) const
			{
				PeakType p;
				p.setPosition(mz);
				return lower_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
			}
			/**
				 @brief Binary search for peak range begin
	
				 @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
			*/
			ConstIterator MZBegin(ConstIterator begin, CoordinateType mz, ConstIterator end) const
			{
				PeakType p;
				p.setPosition(mz);
				return lower_bound(begin, end, p, typename PeakType::PositionLess());
			}
			/**
				 @brief Binary search for peak range end (returns the past-the-end iterator)
	
				 @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
			*/
			ConstIterator MZEnd(CoordinateType mz) const
			{
				PeakType p;
				p.setPosition(mz);
				return upper_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
			}		
			/**
				 @brief Binary search for peak range end (returns the past-the-end iterator)
	
				 @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
			*/
			ConstIterator MZEnd(ConstIterator begin, CoordinateType mz, ConstIterator end) const
			{
				PeakType p;
				p.setPosition(mz);
				return upper_bound(begin, end, p, typename PeakType::PositionLess());
			}
	
			//@}

		 protected:
		 	
			// Docu in base class
			virtual void clearChildIds_()
			{
			}
	
			/// Retention time
			DoubleReal retention_time_;
	
			/// MS level
			UInt ms_level_;
	
			/// Name
			String name_;
			
			///Meta info arrays
			MetaDataArrays meta_data_arrays_;
	};

	///Print the contents to a stream.
	template <typename PeakT>
	std::ostream& operator << (std::ostream& os, const MSSpectrum<PeakT>& spec)
	{
		os << "-- MSSPECTRUM BEGIN --"<<std::endl;

		//spectrum settings
		os << static_cast<const SpectrumSettings&>(spec);

		//peaklist
		os << static_cast<const typename MSSpectrum<PeakT>::ContainerType&>(spec);

		os << "-- MSSPECTRUM END --"<<std::endl;

		return os;
	}

} // namespace OpenMS

#endif // OPENMS_KERNEL_MSSPECTRUM_H
