// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_KERNEL_MSEXPERIMENT_H
#define OPENMS_KERNEL_MSEXPERIMENT_H

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/FORMAT/PersistentObject.h>
#include <OpenMS/CONCEPT/Exception.h>

#include<vector>
#include<algorithm>
#include<limits>

namespace OpenMS
{
/**
	@brief Representation of a mass spectrometry experiment.
	
	Contains the data and metadata of an experiment performed with an MS (or HPLC and MS).
	
	Be carefull when changing the order of contained MSSpectrum instances, if tandem-MS data is
	stored in this class. The only way to find a precursor spectrum of MSSpectrum x is to 
	search for the first spectrum before x that has a lower MS-level!
	
	@note For range operations, see \ref RangeUtils "RangeUtils module"!
	
	@note To iterate over the peaks in all spectra use PeakIterator
			
	@ingroup Kernel
*/
template <typename PeakT = DPeak<1> >
class MSExperiment
            : public std::vector<MSSpectrum<PeakT> >,
            public RangeManager<2, typename PeakT::TraitsType>,
            public ExperimentalSettings,
            public PersistentObject
{
public:

 /**
  	@brief Adaptor class for bidirectional iterator on objects of DPeak<1>  
  	
  	This iterator allows us to move through the data structure in a linear
  	manner i.e. we don't need to jump to the next spectrum manually.
  	
  	The class has a member  DPeakArray<>::iterator pointing
  	to the current peak. The class also remembers the retention time of the current 
  	scan.			
  */
	    template <class IteratorPeakT >
		class PeakIterator : public std::iterator<std::bidirectional_iterator_tag,  IteratorPeakT>
    {
		        
    public:
				
				typedef double CoordinateType;
        typedef IteratorPeakT IteratorPeakType;
				typedef unsigned int difference_type;

        /// Default constructor
        PeakIterator()
                : peak_index_(), rt_(), scan_index_(), exp_()
        {}

        /// Constructor
        PeakIterator(UnsignedInt pind, CoordinateType & co, UnsignedInt sind, MSExperiment<IteratorPeakType>& exp)
                : peak_index_(pind), rt_(co), scan_index_(sind), exp_(&exp)
        {}

        /// Destructor
        ~PeakIterator()
        {}

        /// Copy constructor
        PeakIterator(const PeakIterator& rhs)
                : peak_index_(rhs.peak_index_), rt_(rhs.rt_),
                scan_index_(rhs.scan_index_), exp_(rhs.exp_)
        { }

        /// Assignment operator
        PeakIterator& operator=(const PeakIterator& rhs)
        {
            if (&rhs == this) return *this;

            peak_index_    = rhs.peak_index_;
            rt_                   = rhs.rt_;
            scan_index_    = rhs.scan_index_;
						exp_               = rhs.exp_;

            return (*this);
        }

        /// Test for equality
        bool operator==(const PeakIterator& rhs)
        {
           return ( peak_index_     == rhs.peak_index_ &&
                     								rt_       == rhs.rt_       &&
													scan_index_  == rhs.scan_index_ );
        }

        /// Test for inequality
        bool operator!=(const PeakIterator& rhs)
        {
            return !(*this  == rhs);
        }

        /// Step forward by one (prefix operator)
        PeakIterator& operator++()
        {
            ++peak_index_;
            // test whether we arrived at the end of the current scan
            if ( peak_index_ >=   (*exp_)[scan_index_].size() && scan_index_ !=  ( (*exp_).size() - 1) )
            {
                // we are at the end of a scan, but this scan is not the very last one
                // so we can jump into the next scan
                peak_index_ = 0;
                ++scan_index_;
                rt_   = (*exp_)[scan_index_].getRetentionTime();
            }
            return (*this);
        }

        /// Step backward by one (prefix operator)
        PeakIterator& operator--()
        {
            // test whether we are at the start of a scan
            if (peak_index_  == 0)
            {
                // update scan index and move to end of previous scan
                if (scan_index_ == 0)
                {
                    std::cout << "PeakIterator: In first scan and moving backwards ! " << std::endl;
                    return (*this);
                }
                --scan_index_;
                peak_index_  = ( (*exp_)[scan_index_].size() -1) ;
                rt_                = (*exp_)[scan_index_].getRetentionTime();
            }
            else
            {
                // simply one step backwards
                --peak_index_;
            }
            return (*this);
        }				

        /// Step forward by one (postfix operator)
        PeakIterator operator++(int)
        {
            PeakIterator tmp(*this);
            ++(*this);
            return tmp;
        }

        /// Step backward by one (postfix operator)
        PeakIterator operator--(int)
        {
            PeakIterator tmp(*this);
            --(*this);
            return tmp;
        }

        /// Dereferencing of this pointer yields the underlying peak
        PeakT& operator * ()
        {
            return (*exp_)[scan_index_][peak_index_ ];
        }

        /// Dereferencing of this pointer yields the underlying peak
        PeakT* operator-> ()
        {
            return &((*exp_)[scan_index_][peak_index_ ]);
        }
				
        /** @name Accesssors
        */
        //@{
        /// Returns the current retention time (mutable)
        CoordinateType& getRt() { return rt_; }
        /// Returns the current retention time (not mutable)
        const CoordinateType& getRt() const { return rt_; }
				/// Returns the index of the peak this iterator points to 
				/// NOTE: Call updateRanges() before using this function
				UnsignedInt getPeakNumber()  
				{ 
					if (scan_index_ > 0)
						return (exp_->spectra_lengths_[ (scan_index_-1) ] + peak_index_);
					else
						return peak_index_;
				}
        //@}

    private:
        /// Points to the current peak
        UnsignedInt peak_index_;
        /// Retention time of the current spectrum
        CoordinateType rt_;
        /// Index of the current spectrum
        UnsignedInt scan_index_;
				/// Pointer to the experiment
        MSExperiment<IteratorPeakType> * exp_;
    }
    ; // end of inner class PeakIterator


    /// Spectrum Type
    typedef MSSpectrum<PeakT> SpectrumType;
    /// STL base class type
    typedef std::vector<SpectrumType> Base;
    /// Mutable iterator
    typedef typename std::vector<SpectrumType>::iterator Iterator;
    /// Non-mutable iterator
    typedef typename std::vector<SpectrumType>::const_iterator ConstIterator;
		/// Peak iterator type (for a linear traversal of the data structure)
		typedef PeakIterator< PeakT> PIterator;
	
    /// Peak type
    typedef PeakT PeakType;
    /// Traits types
    typedef typename PeakType::TraitsType TraitsType;
    /// Area type
    typedef DRange<2, TraitsType> AreaType;
    /// Coordinate type of peak positions
    typedef typename TraitsType::CoordinateType CoordinateType;
    /// Intenstiy type of peaks
    typedef typename TraitsType::IntensityType IntensityType;
    /// RangeManager type
    typedef RangeManager<2, TraitsType> RangeManagerType;
		/// const peak reference type
		typedef typename SpectrumType::const_reference ConstPeakReference;
		/// peak reference type
		typedef typename SpectrumType::reference PeakReference;
			 
    /// Constructor
    MSExperiment()
            : std::vector<MSSpectrum<PeakT> >(),
            RangeManagerType(),
            ExperimentalSettings(),
            PersistentObject(),
            ms_levels_(),
            nr_dpoints_(0),
            spectra_lengths_(),
            name_(),
						last_scan_index_(0)
    {
    }

    /// Copy constructor
    MSExperiment(const MSExperiment& source):
            std::vector<MSSpectrum<PeakT> >(source),
            RangeManagerType(source),
            ExperimentalSettings(source),
            PersistentObject(source),
            ms_levels_(source.ms_levels_),
            nr_dpoints_(source.nr_dpoints_),
            spectra_lengths_(source.spectra_lengths_),
            name_(source.name_),
						last_scan_index_(source.last_scan_index_)
    {
    }

    /// Destructor
    ~MSExperiment()
    {
    }

    /// Assignment operator
    MSExperiment& operator= (const MSExperiment& source)
    {
        if (&source == this)
            return *this;

        std::vector<MSSpectrum<PeakT> >::operator=(source);
        ExperimentalSettings::operator=(source);
        PersistentObject::operator=(source);
				
				ms_levels_           = source.ms_levels_;
        nr_dpoints_					 = source.nr_dpoints_;
				spectra_lengths_	 = source.spectra_lengths_;
				name_                 = source.name_;
				last_scan_index_ = source.last_scan_index_;
				
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
        return ExperimentalSettings::operator==(rhs) && std::operator==(rhs,*this);
    }
    /// Equality operator
    bool operator!= (const MSExperiment& rhs) const
    {
        return !(operator==(rhs));
    }

    /**
    	@brief Reads out a 2D Spectrum
      	
      	Container is a DPeakArray<2> or a STL container of DPeak<2> 
      	or DRawDataPoint<2> which supports insert(), end() and back()
    */
    template <class Container>
    void get2DData(Container& cont) const
    {
        const int MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ;
        const int RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT;

        for (typename Base_::const_iterator spec = Base_::begin(); spec != Base_::end(); ++spec)
        {
            if (spec->getMSLevel()!=1)
            {
                continue;
            }
            for (typename MSSpectrum<PeakT>::const_iterator it = spec->
              begin();
              it!=spec->end();
              ++it)
            {
              cont.insert(cont.end(), typename Container::value_type());
              cont.back().getPosition()[RT] = spec->getRetentionTime();
              cont.back().setIntensity(it->getIntensity());
              cont.back().getPosition()[MZ] = it->getPosition()[0];
            }
        }
    }

    /**
    	@brief Assignment of a 2D spectrum to MSExperiment
    	  	
    	Container is a DPeakArray<2> or a STL container of DPeak<2> or DRawDataPoint<2>
    	
    	@note The container has to be sorted according to retention time. Otherwise a Precondition exception is thrown.
    */
    template <class Container>
    void set2DData(const Container& cont) throw (Exception::Precondition)
    {
			SpectrumType* spectrum = 0;
			// If the container is empty, nothing will happen
			if (cont.size() == 0) return;

			const int MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ;
			const int RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT;

			typename PeakType::CoordinateType current_rt = - std::numeric_limits<typename PeakType::CoordinateType>::max();

			for (typename Container::const_iterator iter = cont.begin(); iter != cont.end(); ++iter)
			{
				// check if the retention time time has changed
				if (current_rt != iter->getPosition()[RT] || spectrum == 0)
				{
					if (current_rt > iter->getPosition()[RT])
					{
						throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Input container is not sorted!");
					}
					current_rt =  iter->getPosition()[RT];
					Base_::insert(Base_::end(),SpectrumType());
					spectrum = &(Base_::back());
					spectrum->setRetentionTime(current_rt);
					spectrum->setMSLevel(1);
				}

				// create temporary peak and insert it into spectrum
				spectrum->insert(spectrum->end(), PeakType());
				spectrum->back().setIntensity(iter->getIntensity());
				spectrum->back().getPosition()[0] = iter->getPosition()[MZ];
			}
    }

    /// Returns the name
    const String& getName() const
    {
      return name_;
    }

    /// Sets the name
    void setName(const String& name)
    {
      name_ = name;
    }

   /// Returns an iterator pointing at the first peak
    PIterator peakBegin()
    {
      return PIterator( (unsigned int) 0 , this->at(0).getRetentionTime(), (unsigned int) 0 , *this) ;
    }

    /// Returns an iterator pointing at the last peak
    PIterator peakEnd()
    {
			unsigned int sz = (this->size() - 1);
      return(PIterator( (unsigned int) ( (*this)[sz].size()), (*this)[ sz ].getRetentionTime(), (unsigned int) (sz), *this ) );
    }

    /**
    	@brief Fast search for spectrum range begin
    	
    	@note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    ConstIterator RTBegin(double rt) const
    {
      SpectrumType s;
      s.setRetentionTime(rt);
      return lower_bound(Base_::begin(), Base_::end(), s, typename SpectrumType::RTLess());
    }

    /**
    	@brief Fast search for spectrum range end (returns the past-the-end iterator)
    	
    	@note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    ConstIterator RTEnd(double rt) const
    {
      SpectrumType s;
      s.setRetentionTime(rt);
      return upper_bound(Base_::begin(),Base_::end(), s, typename SpectrumType::RTLess());
    }

    /**
    	@brief Fast search for spectrum range begin
    	
    	@note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    Iterator RTBegin(double rt)
    {
      SpectrumType s;
      s.setRetentionTime(rt);
      return lower_bound(Base_::begin(), Base_::end(), s, typename SpectrumType::RTLess());
    }

    /**
    	@brief Fast search for spectrum range end (returns the past-the-end iterator)
    	
    	@note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
    */
    Iterator RTEnd(double rt)
    {
      SpectrumType s;
      s.setRetentionTime(rt);
      return upper_bound(Base_::begin(),Base_::end(), s, typename SpectrumType::RTLess());
    }


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
		  @brief Returns an array with the accumulated spectra lengths

      The first entry contains the number of peaks in the first spectrum.
			The second entry contains the number of peaks in the first and second spectrum.
			And so on!
		*/
    const std::vector<UnsignedInt>& getSpectraLengths() { return spectra_lengths_; }

		/**
    	@brief Updates the m/z, intensity, retention time and MS level ranges of all spectra with a certain ms level
    	
    	@param ms_level MS level to consider for m/z range , RT range and intensity range (All MS levels if negative)
    */
    void updateRanges(SignedInt ms_level)
    {
      //clear MS levels
      ms_levels_.clear();
      // clear spectra lengths
      spectra_lengths_.clear();
      spectra_lengths_.reserve(this->size());

      //reset mz/rt/int range
      this->clearRanges();
      //reset point count
      nr_dpoints_ = 0;
	
      //empty
      if (this->size()==0)
      {
      	return;
      }

      //update
      for (typename std::vector<MSSpectrum<PeakT> >::iterator it = this->
        begin();
        it!=this->end();
        ++it)
      {
        if (ms_level < SignedInt(0) || SignedInt(it->getMSLevel())==ms_level)
        {  
	        //ms levels
	        if (std::find(ms_levels_.begin(),ms_levels_.end(),it->getMSLevel())==ms_levels_.end())
	        {
	        	ms_levels_.push_back(it->getMSLevel());
	        }
		
					// calculate size
	        nr_dpoints_ += it->size();
		
	        //spectrum lengths
	        spectra_lengths_.push_back( nr_dpoints_ );
              
          //rt
          if (it->getRetentionTime() < RangeManagerType::pos_range_.minX()) RangeManagerType::pos_range_.setMinX(it->getRetentionTime());
          if (it->getRetentionTime() > RangeManagerType::pos_range_.maxX()) RangeManagerType::pos_range_.setMaxX(it->getRetentionTime());
					
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
    }

    /// returns the minimal m/z value
    CoordinateType getMinMZ() const
    {
      return RangeManagerType::pos_range_.min()[1];
    }

    /// returns the maximal m/z value
    CoordinateType getMaxMZ() const
    {
      return RangeManagerType::pos_range_.max()[1];
    }

    /// returns the minimal retention time value
    CoordinateType getMinRT() const
    {
      return RangeManagerType::pos_range_.min()[0];
    }

    /// returns the maximal retention time value
    CoordinateType getMaxRT() const
    {
      return RangeManagerType::pos_range_.max()[0];
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
    UnsignedInt getSize() const { return nr_dpoints_; }

    /// returns an array of the spectrum lengths (all ms level)
    const std::vector<UnsignedInt>& getSpectraLengths() const { return spectra_lengths_; }

    /// returns an array of MS levels
    const std::vector<UnsignedInt>& getMSLevels() const { return ms_levels_; }
    //@}


    /// Mutable access to peak with index @p
    DRawDataPoint<2> getPeak(const UnsignedInt index) throw (Exception::IndexOverflow)
    {
      if (index > nr_dpoints_) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, nr_dpoints_);
				
				UnsignedInt scan_index = 0;
				UnsignedInt peak_index = 0;
				
				// test if requested peak is in the same scan as the last one requested
				UnsignedInt test_offset = 0;
				if (last_scan_index_ == 0)
					test_offset = index;
				else
					test_offset = (index - spectra_lengths_[ (last_scan_index_ - 1) ]);
				
				if (test_offset < (*this)[last_scan_index_].size() && test_offset > 0)
				{
						// good, no binary search necessary
						scan_index = last_scan_index_;
						peak_index = test_offset;
				}
				else
				{
					// bad luck, perform binary search
        	std::vector<UnsignedInt>::iterator it = std::upper_bound(spectra_lengths_.begin(),spectra_lengths_.end(),index);

        	// index of scan is simply the distance to the begin() iterator
        	scan_index =  (it - spectra_lengths_.begin() );        
					last_scan_index_ = scan_index;
					
        	// determine index of peak
        	if (scan_index == 0)
        	{
					  peak_index          =  index;
        	}
        	else
        	{
          	// upper_bound gives last iterator (if several equal values occur),
            // so we have to walk back one step.
            --it;
            peak_index = (index - *it);
        	}
// 					last_scan_index_ = (it - spectra_lengths_.begin() );    
				}
				
				// all information was  collected, compile peak and continue
				DRawDataPoint<2> rp;
				rp.getPosition()[0] = ((*this)[scan_index]).getRetentionTime();
				rp.getPosition()[1] = (*this)[scan_index][peak_index].getPosition()[0];
				rp.getIntensity()    = (*this)[scan_index][peak_index].getIntensity();
				
        return rp;
    }

    /// const access to peak with index @p (call updateRanges() before using this method)
    const DRawDataPoint<2> getPeak(const UnsignedInt index) const throw (Exception::IndexOverflow)
    {
        if (index > nr_dpoints_) throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, nr_dpoints_);
				
				UnsignedInt scan_index = 0;
				UnsignedInt peak_index = 0;
				
				// test if requested peak is in the same scan as the last one requested
				UnsignedInt test_offset = 0;
				if (last_scan_index_ == 0)
					test_offset = index;
				else
					test_offset = (index - spectra_lengths_[ (last_scan_index_ - 1) ]);
				
				if (test_offset < (*this)[last_scan_index_].size() && test_offset > 0)
				{
						// good, no binary search necessary
						scan_index = last_scan_index_;
						peak_index = test_offset;
				}
				else
				{
					// bad luck, perform binary search
        	std::vector<UnsignedInt>::const_iterator it = std::upper_bound(spectra_lengths_.begin(),spectra_lengths_.end(),index);

        	// index of scan is simply the distance to the begin() iterator
        	scan_index =  (it - spectra_lengths_.begin() );        
					last_scan_index_ = scan_index;
					
        	// determine index of peak
        	if (scan_index == 0)
        	{
						last_scan_index_ = 0;
            peak_index          =  index;
        	}
        	else
        	{
          	// upper_bound gives last iterator (if several equal values occur),
            // so we have to walk back one step.
            --it;
            peak_index = (index - *it);
        	}
// 					last_scan_index_ = (it - spectra_lengths_.begin() );    
				}
				
				// all information was  collected, compile peak and continue
				DRawDataPoint<2> rp;
				rp.getPosition()[0] = ((*this)[scan_index]).getRetentionTime();
				rp.getPosition()[1] = (*this)[scan_index][peak_index].getPosition()[0];
				rp.getIntensity()    = (*this)[scan_index][peak_index].getIntensity();
				
        return rp;
    }

	/// Mutable access to retention time of peak with index @p index
// 	CoordinateType& getPeakRt(UnsignedInt index) throw (Exception::IndexOverflow) 
// 	{
// 		if (index > nr_dpoints_)
//             throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, nr_dpoints_);
// 			
// 		std::vector<unsigned int>::iterator it = upper_bound(spectra_lengths_.begin(),spectra_lengths_.end(),index);
// 		unsigned int scan_index =  (it - spectra_lengths_.begin() );
// 		return ((*this)[scan_index]).getRetentionTime();
// 		
// 	}
	
	/// Const access to retention time of peak with index @p index
//	const CoordinateType& getPeakRt(UnsignedInt index) const throw (Exception::IndexOverflow) 
// 	{
// 		if (index > nr_dpoints_)
//             throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, nr_dpoints_);
// 			
// 		std::vector<unsigned int>::iterator it = upper_bound(spectra_lengths_.begin(),spectra_lengths_.end(),index);
// 		unsigned int scan_index =  (it - spectra_lengths_.begin() );
// 		return (*this)[scan_index].getRetentionTime();
// 		
// 	}


    /**
    	@brief Sorts the data points by retention time
    					
    	@p sort_mz : indicates whether the points should be sorted by m/z as well			
    */
    void sortSpectra(bool sort_mz = true)
    {

        std::sort(this->begin(),this->end(),typename SpectrumType::RTLess());

	      if (sort_mz)
      	{
            // sort each spectrum by m/z
            for (Iterator iter = this->begin(); iter != this->end(); ++iter)
            {
							iter->getContainer().sortByPosition();
            }
        }
      
    }
		
	/// Resets all internal values
	void reset()
	{
		Base::clear(); //remove data
		RangeManagerType::clearRanges(); //reset range manager
		ExperimentalSettings::operator=(ExperimentalSettings()); //reset meta info
	}
	
	/// returns the meta information of this experiment
	const ExperimentalSettings& getExperimentalSettings() const { return *this; }
		
protected:

    // Docu in base class
    virtual void clearChildIds_()
    {
        //TODO Persistence
    }
    ;

    ///Base class typedef
    typedef typename std::vector<MSSpectrum<PeakT> > Base_;

    /// MS levels of the data
    std::vector<UnsignedInt> ms_levels_;
    /// Number of all data points
    UnsignedInt nr_dpoints_;
    /// Sums of consecutive spectrum lengths
    std::vector<UnsignedInt> spectra_lengths_;
    /// Name string
    String name_;
		/// Index of last scan retrieved
		mutable UnsignedInt last_scan_index_;
};

///Print the contents to a stream.
template <typename PeakT>
std::ostream& operator << (std::ostream& os, const MSExperiment<PeakT>& exp)
{
    os << "-- MSEXPERIMENT BEGIN --"<<std::endl;

    //experimental settings
    os <<static_cast<const ExperimentalSettings&>(exp);

    //spectra
    for (typename MSExperiment<PeakT>::const_iterator it=exp.begin(); it!=exp.end(); ++it)
    {
        os << *it;
    }

    os << "-- MSEXPERIMENT END --"<<std::endl;

    return os;
}

} // namespace OpenMS

#endif // OPENMS_KERNEL_MSEXPERIMENT_H
