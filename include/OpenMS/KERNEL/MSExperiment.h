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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MSEXPERIMENT_H
#define OPENMS_KERNEL_MSEXPERIMENT_H

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/FORMAT/DB/PersistentObject.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/AreaIterator.h>

#include <vector>
#include <algorithm>
#include <limits>

namespace OpenMS
{
	class Peak1D;
	
	/**
		@brief Representation of a mass spectrometry experiment.
		
		Contains the data and metadata of an experiment performed with an MS (or HPLC and MS).
		
		Be carefull when changing the order of contained MSSpectrum instances, if tandem-MS data is
		stored in this class. The only way to find a precursor spectrum of MSSpectrum x is to 
		search for the first spectrum before x that has a lower MS-level!
		
		@note For range operations, see \ref RangeUtils "RangeUtils module"!

		@ingroup Kernel
	*/
	template <typename PeakT = Peak1D, typename ChromatogramPeakT = ChromatogramPeak>
	class MSExperiment
		:	public std::vector<MSSpectrum<PeakT> >,
			public RangeManager<2>,
			public ExperimentalSettings,
			public PersistentObject
	{
		public:
			/// @name Base type definitions
			//@{
			/// Peak type
			typedef PeakT PeakType;
			/// Chromatogram peak type
			typedef ChromatogramPeakT ChromatogramPeakType;
			/// Area type
			typedef DRange<2> AreaType;
			/// Coordinate type of peak positions
			typedef typename PeakType::CoordinateType CoordinateType;
			/// Intenstiy type of peaks
			typedef typename PeakType::IntensityType IntensityType;
			/// RangeManager type
			typedef RangeManager<2> RangeManagerType;
			/// Spectrum Type
			typedef MSSpectrum<PeakType> SpectrumType;
			/// Chromatogram type
			typedef MSChromatogram<ChromatogramPeakType> ChromatogramType;
			/// STL base class type
			typedef std::vector<SpectrumType> Base;
			//@}

			/// @name Iterator type definitions
			//@{
			/// Mutable iterator
			typedef typename std::vector<SpectrumType>::iterator Iterator;
			/// Non-mutable iterator
			typedef typename std::vector<SpectrumType>::const_iterator ConstIterator;
			/// Mutable area iterator type (for traversal of a rectangular subset of the peaks)
			typedef Internal::AreaIterator<PeakT, PeakT&, PeakT*, Iterator, typename SpectrumType::Iterator> AreaIterator;
			/// Immutable area iterator type (for traversal of a rectangular subset of the peaks)
			typedef Internal::AreaIterator<const PeakT, const PeakT&, const PeakT*, ConstIterator, typename SpectrumType::ConstIterator> ConstAreaIterator;
			//@}
      
			/// Constructor
			MSExperiment() :
			 	Base(),
				RangeManagerType(),
				ExperimentalSettings(),
				PersistentObject(),
				ms_levels_(),
				total_size_(0)
			{
			}

			/// Copy constructor
			MSExperiment(const MSExperiment& source) :
				std::vector<MSSpectrum<PeakT> >(source),
				RangeManagerType(source),
				ExperimentalSettings(source),
				PersistentObject(source),
				ms_levels_(source.ms_levels_),
				total_size_(source.total_size_),
				chromatograms_(source.chromatograms_)
			{
			}

			/// Assignment operator
			MSExperiment& operator= (const MSExperiment& source)
			{
				if (&source == this) return *this;

				Base::operator=(source);
				RangeManagerType::operator=(source);
				ExperimentalSettings::operator=(source);
				PersistentObject::operator=(source);

				ms_levels_           = source.ms_levels_;
				total_size_					 = source.total_size_;
				chromatograms_       = source.chromatograms_;
        
        //no need to copy the alloc?!
        //alloc_
        
				return *this;
			}

      /// Assignment operator
			MSExperiment& operator= (const ExperimentalSettings& source)
			{
				ExperimentalSettings::operator=(source);

				return *this;
			}

			/// Equality operator
			bool operator== (const MSExperiment& rhs) const
			{
				return ExperimentalSettings::operator==(rhs) && std::operator==(rhs,*this) && chromatograms_ == rhs.chromatograms_;
			}
			/// Equality operator
			bool operator!= (const MSExperiment& rhs) const
			{
				return !(operator==(rhs));
			}
			
			///@name Conversion to/from 2D data
			//@{
			/**
				@brief Reads out a 2D Spectrum

				Container can be a PeakArray or an STL container of peaks which
				supports push_back(), end() and back()
			*/
			template <class Container>
			void get2DData(Container& cont) const
			{
				for (typename Base::const_iterator spec = Base::begin(); spec != Base::end(); ++spec)
				{
					if (spec->getMSLevel()!=1)
					{
						continue;
					}
					for (typename SpectrumType::const_iterator it = spec-> begin(); it!=spec->end(); ++it)
					{
						cont.push_back(typename Container::value_type());
						cont.back().setRT(spec->getRT());
						cont.back().setMZ(it->getMZ());
						cont.back().setIntensity(it->getIntensity());
					}
				}
			}

			/**
				@brief Assignment of a 2D spectrum to MSExperiment

				Container can be a PeakArray or an STL container of peaks.

				@exception Exception::Precondition is thrown if the container is not sorted according to retention time (not only in debug mode)
			*/
			template <class Container>
			void set2DData(const Container& cont)
			{
				SpectrumType* spectrum = 0;
				// If the container is empty, nothing will happen
				if (cont.size() == 0) return;

				typename PeakType::CoordinateType current_rt = - (std::numeric_limits<typename PeakType::CoordinateType>::max)();

				for (typename Container::const_iterator iter = cont.begin(); iter != cont.end(); ++iter)
				{
					// check if the retention time time has changed
					if (current_rt != iter->getRT() || spectrum == 0)
					{
						if (current_rt > iter->getRT())
						{
							throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Input container is not sorted!");
						}
						current_rt =  iter->getRT();
						Base::insert(Base::end(),SpectrumType());
						spectrum = &(Base::back());
						spectrum->setRT(current_rt);
						spectrum->setMSLevel(1);
					}

					// create temporary peak and insert it into spectrum
					spectrum->insert(spectrum->end(), PeakType());
					spectrum->back().setIntensity(iter->getIntensity());
					spectrum->back().setPosition(iter->getMZ());
				}
			}
			//@}

			///@name Iterating ranges and areas
			//@{
			///Returns an area iterator for @p area
			AreaIterator areaBegin(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz)
			{
				OPENMS_PRECONDITION(min_rt<=max_rt, "Swapped RT range boundaries!")
					OPENMS_PRECONDITION(min_mz<=max_mz, "Swapped MZ range boundaries!")
					//std::cout << "areaBegin: " << min_rt << " " << max_rt << " " << min_mz << " " << max_mz << std::endl;
					return AreaIterator(this->begin(),RTBegin(min_rt), RTEnd(max_rt), min_mz, max_mz);
			}

			/// Returns an invalid area iterator marking the end of an area
			AreaIterator areaEnd()
			{
				return AreaIterator();
			}

			///Returns a non-mutable area iterator for @p area
			ConstAreaIterator areaBeginConst(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz) const
			{
				OPENMS_PRECONDITION(min_rt<=max_rt, "Swapped RT range boundaries!")
					OPENMS_PRECONDITION(min_mz<=max_mz, "Swapped MZ range boundaries!") 
					//std::cout << "areaBeginConst: " << min_rt << " " << max_rt << " " << min_mz << " " << max_mz << std::endl;
					return ConstAreaIterator(this->begin(),RTBegin(min_rt), RTEnd(max_rt), min_mz, max_mz);
			}

			/// Returns an non-mutable invalid area iterator marking the end of an area
			ConstAreaIterator areaEndConst() const
			{
				return ConstAreaIterator();
			}

			/**
				@brief Fast search for spectrum range begin

        Returns the first scan which has equal or higher (>=) RT than @p rt.

				@note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
			*/
			ConstIterator RTBegin(CoordinateType rt) const
			{
				SpectrumType s;
				s.setRT(rt);
				return lower_bound(Base::begin(), Base::end(), s, typename SpectrumType::RTLess());
			}

			/**
				@brief Fast search for spectrum range end (returns the past-the-end iterator)
        
        Returns the first scan which has higher (>) RT than @p rt.

				@note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
			*/
			ConstIterator RTEnd(CoordinateType rt) const
			{
				SpectrumType s;
				s.setRT(rt);
				return upper_bound(Base::begin(),Base::end(), s, typename SpectrumType::RTLess());
			}

			/**
				@brief Fast search for spectrum range begin

				@note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
			*/
			Iterator RTBegin(CoordinateType rt)
			{
				SpectrumType s;
				s.setRT(rt);
				return lower_bound(Base::begin(), Base::end(), s, typename SpectrumType::RTLess());
			}

			/**
				@brief Fast search for spectrum range end (returns the past-the-end iterator)

				@note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
			*/
			Iterator RTEnd(CoordinateType rt)
			{
				SpectrumType s;
				s.setRT(rt);
				return upper_bound(Base::begin(),Base::end(), s, typename SpectrumType::RTLess());
			}
			//@}

			/**
				@name Range methods  

				@note The range values (min, max, etc.) are not updated automatically. Call updateRanges() to update the values!
			*/
			//@{
			// Docu in base class
			virtual void updateRanges()
			{
				updateRanges(-1);
			}

			/**
				@brief Updates the m/z, intensity, retention time and MS level ranges of all spectra with a certain ms level

				@param ms_level MS level to consider for m/z range , RT range and intensity range (All MS levels if negative)
			*/
			void updateRanges(Int ms_level)
			{
				//clear MS levels
				ms_levels_.clear();

				//reset mz/rt/int range
				this->clearRanges();
				//reset point count
				total_size_ = 0;

				//empty
				if (this->size()==0 && chromatograms_.size() == 0)
				{
					return;
				}

				//update
				for (typename Base::iterator it = this-> begin(); it!=this->end(); ++it)
				{
					if (ms_level < Int(0) || Int(it->getMSLevel())==ms_level)
					{
						//ms levels
						if (std::find(ms_levels_.begin(),ms_levels_.end(),it->getMSLevel())==ms_levels_.end())
						{
							ms_levels_.push_back(it->getMSLevel());
						}

						// calculate size
						total_size_ += it->size();

						//rt
						if (it->getRT() < RangeManagerType::pos_range_.minX()) RangeManagerType::pos_range_.setMinX(it->getRT());
						if (it->getRT() > RangeManagerType::pos_range_.maxX()) RangeManagerType::pos_range_.setMaxX(it->getRT());

						//do not update mz and int when the spectrum is empty
						if (it->size()==0) continue;

						it->updateRanges();

						//mz
						if (it->getMin()[0] < RangeManagerType::pos_range_.minY()) RangeManagerType::pos_range_.setMinY(it->getMin()[0]);
						if (it->getMax()[0] > RangeManagerType::pos_range_.maxY()) RangeManagerType::pos_range_.setMaxY(it->getMax()[0]);

						//int
						if (it->getMinInt() < RangeManagerType::int_range_.minX()) RangeManagerType::int_range_.setMinX(it->getMinInt());
						if (it->getMaxInt() > RangeManagerType::int_range_.maxX()) RangeManagerType::int_range_.setMaxX(it->getMaxInt());

					}
				}
				std::sort(ms_levels_.begin(), ms_levels_.end());

				if (this->chromatograms_.size() == 0)
				{
					return;
				}				

				//TODO CHROM update intensity, m/z and RT according to chromatograms as well! (done????)

				for (typename std::vector<ChromatogramType>::iterator it = chromatograms_.begin(); it != chromatograms_.end(); ++it)
				{

          // ignore TICs and ECs (as these are usually positioned at 0 and therefor lead to a large white margin in plots if included)
          if (it->getChromatogramType() == ChromatogramSettings::TOTAL_ION_CURRENT_CHROMATOGRAM ||
              it->getChromatogramType() == ChromatogramSettings::EMISSION_CHROMATOGRAM)
          {
            continue;
          }

					// update MZ
					if (it->getMZ() < RangeManagerType::pos_range_.minY()) RangeManagerType::pos_range_.setMinY(it->getMZ());
					if (it->getMZ() > RangeManagerType::pos_range_.maxY()) RangeManagerType::pos_range_.setMaxY(it->getMZ());

					// do not update RT and in if the specturm is empty
					if (it->size() == 0) continue;

					total_size_ += it->size();

					it->updateRanges();

					// RT
					if (it->getMin()[0] < RangeManagerType::pos_range_.minX()) RangeManagerType::pos_range_.setMinX(it->getMin()[0]);
					if (it->getMax()[0] > RangeManagerType::pos_range_.maxX()) RangeManagerType::pos_range_.setMaxX(it->getMax()[0]);

					// int
					if (it->getMinInt() < RangeManagerType::int_range_.minX()) RangeManagerType::int_range_.setMinX(it->getMinInt());
          if (it->getMaxInt() > RangeManagerType::int_range_.maxX()) RangeManagerType::int_range_.setMaxX(it->getMaxInt());				
				}
			}

			/// returns the minimal m/z value
			CoordinateType getMinMZ() const
			{
				return RangeManagerType::pos_range_.minPosition()[1];
			}

			/// returns the maximal m/z value
			CoordinateType getMaxMZ() const
			{
				return RangeManagerType::pos_range_.maxPosition()[1];
			}

			/// returns the minimal retention time value
			CoordinateType getMinRT() const
			{
				return RangeManagerType::pos_range_.minPosition()[0];
			}

			/// returns the maximal retention time value
			CoordinateType getMaxRT() const
			{
				return RangeManagerType::pos_range_.maxPosition()[0];
			}

			/**
				@brief Returns RT and m/z range the data lies in.

				RT is dimension 0, m/z is dimension 1
			*/
			const AreaType& getDataRange() const
			{
				return RangeManagerType::pos_range_;
			}

			/// returns the total number of peaks
			UInt64 getSize() const
		 	{
			 	return total_size_;
		 	}

			/// returns an array of MS levels
			const std::vector<UInt>& getMSLevels() const
		 	{
			 	return ms_levels_;
		 	}
			//@}
			
			///@name Sorting spectra and peaks
			//@{
			/**
				@brief Sorts the data points by retention time
				
				@param sort_mz if @em true, spectra are sorted by m/z position as well
			*/
			void sortSpectra(bool sort_mz = true)
			{
				std::sort(this->begin(),this->end(),typename SpectrumType::RTLess());

				if (sort_mz)
				{
					// sort each spectrum by m/z
					for (Iterator iter = this->begin(); iter != this->end(); ++iter)
					{
						iter->sortByPosition();
					}
				}
			}

			/** 
				@brief Sorts the data points of the chromatograms by m/z

				@param sort_rt if @em true, chromatograms are sorted by rt position as well
			*/
			void sortChromatograms(bool sort_rt = true)
			{
				// sort the chromatograms according to their product m/z
				std::sort(chromatograms_.begin(), chromatograms_.end(), typename ChromatogramType::MZLess());

				if (sort_rt)
				{
					for (typename std::vector<ChromatogramType>::iterator it = chromatograms_.begin(); it != chromatograms_.end(); ++it)
					{
						it->sortByPosition();
					}
				}
			}
			
			/**
				@brief Checks if all spectra are sorted with respect to ascending RT
				
				@param check_mz if @em true, checks if all peaks are sorted with respect to ascending m/z
			*/
			bool isSorted(bool check_mz = true ) const
			{
				//check RT positions
				for (Size i=1; i<this->size(); ++i)
				{
					if (this->operator[](i-1).getRT()>this->operator[](i).getRT()) return false;
				}
				//check spectra
				if (check_mz)
				{
					for (Size i=0; i<this->size(); ++i)
					{
						if (!this->operator[](i).isSorted()) return false;
					}
				}
				// TODO CHROM
				return true;
			}
			//@}
			
			/// Resets all internal values
			void reset()
			{
				Base::clear(); //remove data
				RangeManagerType::clearRanges(); //reset range manager
				ExperimentalSettings::operator=(ExperimentalSettings()); //reset meta info
			}
			
			/**
				@brief Clears the meta data arrays of all contained spectra (float, integer and string arrays)
				
				@return @em true if meta data arrays were present and removed. @em false otherwise.
			*/
			bool clearMetaDataArrays()
			{
				bool meta_present = false;
				for (Size i=0; i< this->size(); ++i)
				{
					if (this->operator[](i).getFloatDataArrays().size()!=0 || this->operator[](i).getIntegerDataArrays().size()!=0 || this->operator[](i).getStringDataArrays().size()!=0)
					{
						meta_present = true;
					}
					this->operator[](i).getStringDataArrays().clear();
					this->operator[](i).getIntegerDataArrays().clear();
					this->operator[](i).getFloatDataArrays().clear();
				}
				return meta_present;
			}
			
			/// returns the meta information of this experiment (const access)
			const ExperimentalSettings& getExperimentalSettings() const
			{ 
				return *this; 
			}
			/// returns the meta information of this experiment (mutable access)
			ExperimentalSettings& getExperimentalSettings()
			{ 
				return *this; 
			}

			/**
				@brief Returns the precursor spectrum of the scan pointed to by @p iterator

				If there is no precursor scan the past-the-end iterator is returned.
			*/
			ConstIterator getPrecursorSpectrum(ConstIterator iterator) const
			{
				if (iterator==this->end() || iterator==this->begin())
				{
					return this->end();
				}
				UInt ms_level = iterator->getMSLevel();
				do
				{
					--iterator;
					if (iterator->getMSLevel() < ms_level)
					{
						return iterator;
					}
				}
				while (iterator!=this->begin());

				return this->end();
			}
			
			/// Swaps the content of this map with the content of @p from
			void swap(MSExperiment& from)
			{
				MSExperiment tmp;
				
				//swap range information
				tmp.RangeManagerType::operator=(*this);
				this->RangeManagerType::operator=(from);
				from.RangeManagerType::operator=(tmp);
				
				//swap experimental settings
				tmp.ExperimentalSettings::operator=(*this);
				this->ExperimentalSettings::operator=(from);
				from.ExperimentalSettings::operator=(tmp);

				//swap persistent object
				tmp.PersistentObject::operator=(*this);
				this->PersistentObject::operator=(from);
				from.PersistentObject::operator=(tmp);
			
				// swap chromatograms
				std::swap(chromatograms_, from.chromatograms_);

				//swap peaks
				Base::swap(from);

				//swap remaining members
				ms_levels_.swap(from.ms_levels_);
				std::swap(total_size_,from.total_size_);
			}

			/// sets the chromatogram list
			void setChromatograms(const std::vector<MSChromatogram<ChromatogramPeakType> >& chromatograms)
			{
				chromatograms_ = chromatograms;
			}

			/// adds a chromatogram to the list
			void addChromatogram(const MSChromatogram<ChromatogramPeakType> & chromatogram)
			{
				chromatograms_.push_back(chromatogram);
			}

      /// returns the chromatogram list
      const std::vector<MSChromatogram<ChromatogramPeakType> >& getChromatograms() const
      {
        return chromatograms_;
      }
    
      /// returns the total ion chromatogram (TIC)
      const MSChromatogram<ChromatogramPeakType> getTIC() const
      {
        // The TIC is (re)calculated from the MS1 spectra. Even if MSExperiment does not contain a TIC chromatogram explicitly, it can be reported.
        MSChromatogram<ChromatogramPeakType> TIC;
        for (typename Base::const_iterator spec_it = this->begin(); spec_it != this->end(); ++spec_it)
				{
          if (spec_it->getMSLevel() == 1)
          {
            DoubleReal totalIntensity = 0;
            // sum intensities of a spectrum
            for (typename SpectrumType::const_iterator peak_it = spec_it->begin(); peak_it != spec_it->end(); ++peak_it)
            {
              totalIntensity += peak_it->getIntensity();
            }
            // fill chromatogram
            ChromatogramPeakType peak;
            peak.setRT(spec_it->getRT());
            peak.setIntensity(totalIntensity);
            TIC.push_back(peak);
          }
				}
        return TIC;
      }
    
      /**
				@brief Clears all data and meta data
				
				@param clear_meta_data If @em true, all meta data is cleared in addition to the data.
			*/ 
			void clear(bool clear_meta_data)
			{
				Base::clear();
				
				if (clear_meta_data)
				{
					clearRanges();
					clearId();
					this->ExperimentalSettings::operator=(ExperimentalSettings()); // no "clear" method
					chromatograms_.clear();
					ms_levels_.clear();
					total_size_ = 0;
				}
			}
			
		protected:
	
	    // Docu in base class
	    virtual void clearChildIds_()
	    {
	    	for (Size i=0; i<this->size(); ++i)
				{
					this->operator[](i).clearId(true);
				}
	    }

	    /// MS levels of the data
	    std::vector<UInt> ms_levels_;
	    /// Number of all data points
	    UInt64 total_size_;

			/// chromatograms 
			std::vector<MSChromatogram<ChromatogramPeakType> > chromatograms_;
	};
	
	///Print the contents to a stream.
	template <typename PeakT, typename ChromatogramPeakT>
	std::ostream& operator << (std::ostream& os, const MSExperiment<PeakT, ChromatogramPeakT>& exp)
	{
	    os << "-- MSEXPERIMENT BEGIN --"<<std::endl;
	
	    //experimental settings
	    os <<static_cast<const ExperimentalSettings&>(exp);
	
	    //spectra
	    for (typename MSExperiment<PeakT>::const_iterator it=exp.begin(); it!=exp.end(); ++it)
	    {
	      os << *it;
	    }

			//chromatograms
			for (typename std::vector<MSChromatogram<ChromatogramPeakT> >::const_iterator it = exp.getChromatograms().begin(); it != exp.getChromatograms().end(); ++it)
			{
				os << *it;
			}

	    os << "-- MSEXPERIMENT END --"<<std::endl;
	
	    return os;
	}

} // namespace OpenMS

#endif // OPENMS_KERNEL_MSEXPERIMENT_H
