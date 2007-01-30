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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MSEXPERIMENTEXTERN_H
#define OPENMS_KERNEL_MSEXPERIMENTEXTERN_H

#define _FILE_OFFSET_BITS 64

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include<vector>
#include<algorithm>
#include<limits>

#include <iostream>
#include <fstream>

#include <stdio.h>
// #include <sys/mman.h>
// #include <sys/types.h>
#include <sys/errno.h>

namespace OpenMS
{
	/**
		@brief Representation of a mass spectrometry experiment using an external datastructure to store large data sets.
		
		This data structures has the same interface as MSExperiment but uses a ring buffer and stores only a subset of
		all scans in the RAM. Scans are dynamically written to the hard disk and re-loaded if needed.
		
		@note If your LC-MS map is really large, you might want to compile this class with LFS (large file support) such
		that Linux / C can access files > 2 GB. In this case, you will need to comple OpenMS with -D_FILE_OFFSET_BITS = 64.
		
		@note This container works only with DPeak's. Other point types are not supported.
					
		@ingroup Kernel
	*/
	template < typename PeakT = DPeak<1> >
	class MSExperimentExtern
		: public RangeManager<2, typename PeakT::TraitsType>
	{  
		public:
	
	    /// Iterator class
	    template <class IteratorPeakT, class ValueT, class ReferenceT, class PointerT>
	    class MSExperimentExternIterator : public std::iterator<std::bidirectional_iterator_tag,  ValueT>
	    {
				friend class MSExperimentExtern;
	
	    public:
	        /**	@name Type definitions */
	        //@{
	        typedef IteratorPeakT IteratorPeakType;
	        typedef ValueT value_type;
	        typedef MSExperimentExtern<IteratorPeakType> ExperimentType;
	        typedef ReferenceT reference;
	        typedef PointerT pointer;
	        typedef std::random_access_iterator_tag iterator_category;
	        typedef unsigned int difference_type;
	        //@}
	
	        MSExperimentExternIterator()
	                : exp_(0), position_(0)
	        {}
	
	        MSExperimentExternIterator(const ExperimentType * exp , unsigned int pos)
	        {
	            exp_       =  (ExperimentType* ) exp;
	            position_ = pos;
	        }
	
	        MSExperimentExternIterator(ExperimentType * exp , unsigned int pos)
	                : exp_(exp), position_(pos)
	        {}
	
	        MSExperimentExternIterator(const ExperimentType& source)
	                : exp_(source.exp_), position_(source.position_)
	        {}
	
	        ~MSExperimentExternIterator()
	        {}
	
	        MSExperimentExternIterator& operator = (const MSExperimentExternIterator& rhs)
	        {
	            if (this==&rhs)
	                return *this;
	
	            exp_      = rhs.exp_;
	            position_= rhs.position_;
	
	            return *this;
	        }
	
	        bool operator < (const MSExperimentExternIterator& it) const
	        {
	            return position_ < it.position_;
	        }
	
	        bool operator > (const MSExperimentExternIterator& it) const
	        {
	            return position_ > it.position_;
	        }
	
	        bool operator <= (const MSExperimentExternIterator& it) const
	        {
	            return (position_ < it.position_ || position_ == it.position_);
	        }
	
	        bool operator >= (const MSExperimentExternIterator& it) const
	        {
	            return (position_ > it.position_ || position_ == it.position_);
	        }
	
	        bool operator == (const MSExperimentExternIterator& it) const
	        {
	            return position_ == it.position_ && exp_ == it.exp_;
	        }
	
	        bool operator != (const MSExperimentExternIterator& it) const
	        {
	            return position_ != it.position_ || exp_ != it.exp_;
	        }
	
	        MSExperimentExternIterator& operator ++ ()
	        {
	            position_ += 1;
	            return *this;
	        }
	
	        MSExperimentExternIterator operator ++ (int)
	        {
	            MSExperimentExternIterator tmp(*this);
	            ++(*this);
	            return tmp;
	        }
	
	        MSExperimentExternIterator& operator -- ()
	        {
	            position_ -= 1;
	            return *this;
	        }
	
	        MSExperimentExternIterator operator -- (int)
	        {
	            MSExperimentExternIterator tmp(*this);
	            --(*this);
	            return tmp;
	        }
	
	        MSExperimentExternIterator operator - (difference_type n) const
	        {
	            MSExperimentExternIterator tmp(*this);
	            tmp.position_ -= n;
	            return tmp;
	        }
	
	        MSExperimentExternIterator operator + (difference_type n) const
	        {
	            MSExperimentExternIterator tmp(*this);
	            tmp.position_ += n;
	            return tmp;
	        }
	
	        MSExperimentExternIterator& operator += (difference_type n)
	        {
	            position_ += n;
	            return *this;
	        }
	
	        MSExperimentExternIterator& operator -= (difference_type n)
	        {
	            position_ -= n;
	            return *this;
	        }
	
	        friend difference_type operator - ( const MSExperimentExternIterator& i1, const MSExperimentExternIterator& i2 )
	        {
	            return (i1.position_ - i2.position_);
	        }
	
	        friend MSExperimentExternIterator operator + ( difference_type n, const MSExperimentExternIterator& i )
	        {
	            MSExperimentExternIterator tmp(i);
	            tmp.position_ += n;
	            return tmp;
	        }
	
	        reference operator * ()
	        {
	            return (*exp_)[position_];
	        }
	
	        pointer operator -> ()
	        {
	            return &((*exp_)[position_]);
	        }
	
	        pointer operator -> () const
	        {
	            return &((*exp_)[position_]);
	        }
	
	        reference operator [] (difference_type n)
	        {
	            return (*this)+n;
	        }
	
	        friend void swap(MSExperimentExternIterator& i1, MSExperimentExternIterator& i2)
	        {
	            unsigned int tmp = i1.position_;
	            i1.position_ = i2.position_;
	            i2.position_ = tmp;
	        }
	
	    protected:
	
	        ExperimentType * exp_;
	        unsigned int position_;
	    };
	    // end of class MSExperimentExternIterator
	
	    typedef PeakT PeakType;
	    typedef typename PeakT::IntensityType IntensityType;
	    typedef typename PeakT::PositionType PositionType;
	    typedef MSSpectrum<PeakT> SpectrumType;
	    typedef typename SpectrumType::ContainerType ContainerType;
	    typedef typename PeakType::CoordinateType CoordinateType;
	    typedef MSExperiment<PeakType> ExperimentType;
	
	    typedef typename PeakType::TraitsType TraitsType;
	    /// Area type
	    typedef DRange<2, TraitsType> AreaType;
	    typedef RangeManager<2, TraitsType> RangeManagerType;
	
	    typedef MSExperimentExternIterator<PeakType, SpectrumType, SpectrumType&, SpectrumType*> Iterator;
	    typedef MSExperimentExternIterator<PeakType, SpectrumType, const SpectrumType&, const SpectrumType*> ConstIterator;
	
			/// Mutable area iterator type (for traversal of a rectangular subset of the peaks)
			typedef Internal::AreaIterator<PeakType, PeakType&, PeakType*, Iterator, typename SpectrumType::Iterator> AIterator;
			/// Immutable area iterator type (for traversal of a rectangular subset of the peaks)
			typedef Internal::AreaIterator<PeakType, const PeakType&, const PeakType*, ConstIterator, typename SpectrumType::ConstIterator> AConstIterator;
	    typedef std::reverse_iterator<Iterator> ReverseIterator;
	    typedef std::reverse_iterator<ConstIterator> ConstReverseIterator;
	
	    typedef typename ExperimentType::value_type value_type;
	    typedef typename ExperimentType::size_type size_type;
	    typedef typename ExperimentType::difference_type difference_type;
	    typedef typename ExperimentType::reference reference;
	    typedef typename ExperimentType::const_reference const_reference;
	    typedef typename ExperimentType::pointer pointer;
	
	    /// const reference to peak type
	    typedef typename SpectrumType::const_reference const_preference;
	    /// mutable reference to peak type
	    typedef typename SpectrumType::reference preference;
	
	    typedef Iterator iterator;
	    typedef ConstIterator const_iterator;
	    typedef ReverseIterator reverse_iterator;
	    typedef ConstReverseIterator const_reverse_iterator;
	
	    /// Standard constructor, allocates a buffer of size 100
	    MSExperimentExtern()
	            : buffer_size_(100),
	            scan_location_(), current_scan_(0),
	            buffer_index_(0), scan2buffer_(),
	            buffer2scan_(), exp_(), pFile_(0),
	            nr_dpoints_(0), ms_levels_()
	    { 
	        file_name_ = String("msexp_") + std::rand();
	        exp_.resize(buffer_size_);
	        buffer2scan_.resize(buffer_size_);
	    }
	
	    /// Copy constructor: copies the content of the temporary file as well (slow !)
	    MSExperimentExtern(const MSExperimentExtern& source)
	            : RangeManagerType(source), buffer_size_(source.buffer_size_),
	            scan_location_(source.scan_location_), current_scan_(source.current_scan_),
	            buffer_index_(source.buffer_index_),  scan2buffer_(source.scan2buffer_),
	            buffer2scan_(source.buffer2scan_), scan_sizes_(source.scan_sizes_), 
							exp_(source.exp_), nr_dpoints_(source.nr_dpoints_), 
							ms_levels_(source.ms_levels_)
	    {
	        // genarete new temp file and copy the old one
					std::remove( file_name_ .c_str());
	        file_name_ = String("msexp_") + std::rand();
	        copyTmpFile_(source.file_name_);
	    }
	
	    /// Destructor
	    virtual ~MSExperimentExtern()
	    {
	        // delete temporary file
	        std::remove( file_name_ .c_str());
	    }
	
	    /// Assignment operator
	    MSExperimentExtern & operator= (const MSExperimentExtern& source)
	    {
	        if (&source == this) return *this;
	        
					RangeManagerType::operator=(source);
	        buffer_size_         = source.buffer_size_;
	        scan_location_     = source.scan_location_;
	        current_scan_      = source.current_scan_;
	        buffer_index_        = source.buffer_index_;
					scan2buffer_        = source.scan2buffer_;
	        buffer2scan_        = source.buffer2scan_;
	        scan_sizes_        = source.scan_sizes_;
	        exp_				            = source.exp_;
	        nr_dpoints_         = source.nr_dpoints_;
					ms_levels_         = source.ms_levels_;
	
	        // generate new name for temp file
	        std::remove( file_name_ .c_str());
	        file_name_ = String("msexp_") + std::rand();
	        // and copy the old one
	        copyTmpFile_(source.file_name_);
					
	        return *this;
	    }
	
	    /// Equality operator
	    bool operator== (const MSExperimentExtern& rhs) const
	    {
	        return (buffer_size_      == rhs.buffer_size_ &&
	                scan_location_     == rhs.scan_location_ &&
	                buffer_index_        == rhs.buffer_index_ &&
	                scan2buffer_         == rhs.scan2buffer_ &&
	                buffer2scan_         == rhs.buffer2scan_ &&
									scan_sizes_         == rhs.scan_sizes_ &&
	                exp_				             == rhs.exp_             &&
									nr_dpoints_					 == rhs.nr_dpoints_ &&
	                current_scan_       == rhs.current_scan_ &&
									ms_levels_           == rhs.ms_levels_);
	    }
	
	    /// Equality operator
	    bool operator!= (const MSExperimentExtern& rhs) const
	    {
	        return !(operator==(rhs));
	    }
	
	    template <class Container>
	    void get2DData(Container& cont) const
	    {
	       	const int MZ = DimensionDescription < LCMS_Tag >::MZ;
	       	const int RT = DimensionDescription < LCMS_Tag >::RT;
	
					SpectrumType spec;
					
					for (UnsignedInt i=0; i<scan2buffer_.size(); ++i)
	        {
							spec = (*this)[i];
	            if (spec.getMSLevel()!=1)
	            {
	                continue;
	            }
	            for (typename MSSpectrum<PeakT>::const_iterator it = spec.begin(); it!=spec.end(); ++it)
	            {
	              cont.insert(cont.end(), typename Container::value_type());
	              cont.back().getPosition()[RT] = spec.getRetentionTime();
	              cont.back().setIntensity(it->getIntensity());
	              cont.back().getPosition()[MZ] = it->getPosition()[0];
	            }
	        }
	    }
	
	    template <class Container>
	    void set2DData(Container& cont)
	    {
	        SpectrumType* spectrum = 0;
	        /// If the container is emptry, nothing will happen
	        if (cont.size() == 0) return;
	
	        const int MZ = DimensionDescription < LCMS_Tag >::MZ;
	        const int RT = DimensionDescription < LCMS_Tag >::RT;
	
	        typename PeakType::CoordinateType current_rt = -1.0*std::numeric_limits<typename PeakType::CoordinateType>::max();
	
	        for (typename Container::const_iterator iter = cont.begin(); iter != cont.end(); ++iter)
	        {
	          // check if retention time has changed
	          if (current_rt != iter->getPosition()[RT] || spectrum == 0)
	          {
	            if (current_rt > iter->getPosition()[RT])
	            {
	              throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Input container is not sorted!");
	            }
	            current_rt =  iter->getPosition()[RT];
	            this->push_back(SpectrumType());
							spectrum = &(this->back());
	            spectrum->setRetentionTime(current_rt);
	            spectrum->setMSLevel(1);
	          }
	
	          // create temporary peak and insert it into spectrum
	          spectrum->insert(spectrum->end(), PeakType());
	          spectrum->back().setIntensity(iter->getIntensity());
	          spectrum->back().getPosition()[0] = iter->getPosition()[MZ];
	        }
	    }
	
	    /// Update the range informations
			virtual void updateRanges()
	    {
	      updateRanges(-1);
	    }
			
			/**
	    	@brief Updates the m/z, intensity, retention time and MS level ranges of all spectra with a certain ms level
	    	
	    	@param ms_level MS level to consider for m/z range , RT range and intensity range (All MS levels if negative)
	    */
	    void updateRanges(SignedInt ms_level)
	    {
				//clear MS levels
	      ms_levels_.clear();
	
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
	      for (UnsignedInt i=0;i<scan2buffer_.size();++i)
	      {
					SpectrumType spec_temp = (*this)[i];
					
	        if (ms_level < SignedInt(0) || SignedInt(spec_temp.getMSLevel())==ms_level)
	        {  
		        //ms levels
		        if (std::find(ms_levels_.begin(),ms_levels_.end(),spec_temp.getMSLevel())==ms_levels_.end())
		        {
		        	ms_levels_.push_back(spec_temp.getMSLevel());
		        }
			
						// calculate size
		        nr_dpoints_ += spec_temp.size();
	              
	          //rt
	          if (spec_temp.getRetentionTime() < RangeManagerType::pos_range_.minX()) RangeManagerType::pos_range_.setMinX(spec_temp.getRetentionTime());
	          if (spec_temp.getRetentionTime() > RangeManagerType::pos_range_.maxX()) RangeManagerType::pos_range_.setMaxX(spec_temp.getRetentionTime());
						
						//do not update mz and intensities if spectrum is empty
						if (spec_temp.size()==0) continue;
						
	          spec_temp.updateRanges();
	          
						//mz
	          if (spec_temp.getMin()[0] < RangeManagerType::pos_range_.minY()) RangeManagerType::pos_range_.setMinY(spec_temp.getMin()[0]);
	          if (spec_temp.getMax()[0] > RangeManagerType::pos_range_.maxY()) RangeManagerType::pos_range_.setMaxY(spec_temp.getMax()[0]);
	
	          //int
	          if (spec_temp.getMinInt() < RangeManagerType::int_range_.minX()) RangeManagerType::int_range_.setMinX(spec_temp.getMinInt());
	          if (spec_temp.getMaxInt() > RangeManagerType::int_range_.maxX()) RangeManagerType::int_range_.setMaxX(spec_temp.getMaxInt());
	
	        }
	      }
	      std::sort(ms_levels_.begin(), ms_levels_.end());
	    
	    }
	
	    /// Sorts the spectra (if @p sort_mz is set to true, the data points are sorted by m/z as well)
	    void sortSpectra(bool sort_mz = true)
	    {
	        std::sort(this->begin(),this->end(),typename SpectrumType::RTLess());
	
	        if (sort_mz)
	        {
	            // sort each spectrum by m/z
	            for (Iterator iter = this->begin();
	                    iter != this->end();
	                    ++iter)
	            {
	                iter->getContainer().sortByPosition();
	            }
	        }
	    }
	
	    /**
	    	@brief Fast search for spectrum range begin
	    	
	    	@note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
	    */
	    ConstIterator RTBegin(double rt) const
	    {
	        SpectrumType s;
	        s.setRetentionTime(rt);
	        return lower_bound(begin(), end(), s, typename SpectrumType::RTLess());
	    }
	
	    /**
	    	@brief Fast search for spectrum range end (returns the past-the-end iterator)
	    	
	    	@note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
	    */
	    ConstIterator RTEnd(double rt) const
	    {
	        SpectrumType s;
	        s.setRetentionTime(rt);
	        return upper_bound(begin(),end(), s, typename SpectrumType::RTLess());
	    }
	
	    /**
	    	@brief Binary search for rt range begin
	    	
	    	@note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
	    */
	    Iterator RTBegin(double rt)
	    {
	        SpectrumType s;
	        s.setRetentionTime(rt);
	        return lower_bound(begin(),end(), s, typename SpectrumType::RTLess());
	    }
	
	    /**
	       	@brief Binary search for rt range end (returns the past-the-end iterator)
	       	
	       	@note Make sure the spectra are sorted with respect to retention time! Otherwise the result is undefined.
	       */
	    Iterator RTEnd(double rt)
	    {
	        SpectrumType s;
	        s.setRetentionTime(rt);
	        return upper_bound(begin(),end(), s, typename SpectrumType::RTLess());
	    }
	
	    /// See std::vector documentation.
	    Iterator begin() { return Iterator(this,(unsigned int)0); }
	
	    /// See std::vector documentation.
	    Iterator end() { return Iterator(this,this->size()); }
	
	    /// See std::vector documentation.
	    ConstIterator begin() const { return ConstIterator(this,(unsigned int)0); }
	
	    /// See std::vector documentation.
	    ConstIterator end() const { return ConstIterator(this,this->size()); }
	
	    /// See std::vector documentation.
	    ReverseIterator rbegin() { return ReverseIterator(end()); }
	
	    /// See std::vector documentation.
	    ReverseIterator rend() { return ReverseIterator(begin()); }
	
	    /// See std::vector documentation.
	    ConstReverseIterator rbegin() const { return ConstReverseIterator(end()); }
	
	    /// See std::vector documentation.
	    ConstReverseIterator rend() const { return ConstReverseIterator(begin()); }
	
	    /// See std::vector documentation.
	    reference back() { return  this->at( (scan2buffer_.size() -1) ) ; }
	
	    void push_back(const SpectrumType& spec)
	    {
	//          std::cout << "Inserting scan " << current_scan_ << std::endl;
	//          std::cout << "buffer capacity: " << buffer_size_ << " buffer index: " << buffer_index_ << " buffer size: " << exp_.size() << std::endl;
	        if (buffer_index_ < buffer_size_)
	        {
	//      				std::cout << "Writing in buffer at pos: " << buffer_index_ << std::endl;
	            // test if we already wrote at this buffer position
	            if (current_scan_ >= buffer_size_)
	            {
	                // yes => store scan at current buffer position and then overwrite
	// 								std::cout << "Position in buffer occupied. Writing to hard disk and then overwrite." << std::endl;
	                writeScan_(buffer2scan_[buffer_index_], exp_[buffer_index_] );
	            }
	
	            exp_[buffer_index_] = spec;
	            scan2buffer_.push_back(buffer_index_);
	            buffer2scan_[buffer_index_++] = current_scan_++;
	        }
	        else
	        {
	// 	         	std::cout << "Buffer full. Overwriting buffer at pos 0."   << std::endl;
	            buffer_index_ = 0; 																		// reset buffer index
	            writeScan_(buffer2scan_[buffer_index_],  exp_[buffer_index_] ); 		// write content of buffer to temp. file
	            exp_[buffer_index_] = spec;														// store new spectrum
	
	            scan2buffer_.push_back(buffer_index_);
	            buffer2scan_[buffer_index_++] = current_scan_++;
	        }
	//         std::cout << "scan2buffer : " << scan2buffer_[ (current_scan_-1)] << std::endl;
	//         std::cout << "buffer2scan: " << buffer2scan_[ ( scan2buffer_[ (current_scan_-1)]  )] << std::endl;
	    }
	
	    /// see std::vector (additionally test if scan is in buffer or needs to be read from temp file)
	    reference operator[] (size_type n)
	    {
	//         std::cout << "operator[" << n << "]" << std::endl;
	        // test if current scan is in buffer
	        UnsignedInt b = scan2buffer_[n];
	        if (buffer2scan_[b] != n)
					{	
	// 						std::cout << "scan not in buffer." << std::endl;
	            storeInBuffer_(n);	// scan is not in buffer, needs to be read from file
					 }
	        b = scan2buffer_[n];
	        return exp_[b];
	    }
	
	    /// see std::vector (additionally test if scan is in buffer or needs to be read from temp file)
	    const_reference operator[] (size_type n) const
	    {
	// 				std::cout << "operator[" << n << "] const" << std::endl;
	        // test if current scan is in buffer
	        UnsignedInt b = scan2buffer_[n];
	        if (buffer2scan_[b] != n)
	        {
	//             std::cout << "scan not in buffer." << std::endl;
	            storeInBuffer_(n);	// scan is not in buffer, needs to be read from file
	        }
	        b = scan2buffer_[n];
	        return exp_[b];
	    }
	
	    /// see std::vector (additionally test if scan is in buffer or needs to be read from temp file)
	    reference at(size_type n)
	    {
	//         std::cout << "at(" << n << ")" << std::endl;
	        // test if current scan is in buffer
	        UnsignedInt b = scan2buffer_[n];
	        if (buffer2scan_[b] != n) 
					{
	// 					std::cout << "scan not in buffer." << std::endl;
	           storeInBuffer_(n);	// scan is not in buffer, needs to be read from file
					}
	        b = scan2buffer_[n];
	        return exp_.at(b);
	    }
	
	    /// see std::vector (additionally test if scan is in buffer or needs to be read from temp file)
	    const_reference at(size_type n) const
	    {
	//         std::cout << "at(" << n << ") const" << std::endl;
	        // test if current scan is in buffer
	        UnsignedInt b = scan2buffer_[n];
	        if (buffer2scan_[b] != n)
	        {    
	// 					std::cout << "scan not in buffer." << std::endl;
						storeInBuffer_(n);	// scan is not in buffer, needs to be read from file
					}
	        b = scan2buffer_[n];
	        return exp_.at(b);
	    }
	
	    /// Sets buffer size
	    void setBufferSize(UnsignedInt sz) {  buffer_size_ = sz; }
	    /// Returns the buffer size (mutable)
	    UnsignedInt& getBufferSize() { return buffer_size_; }
	    /// Returns the buffer size (not mutable)
	    const UnsignedInt& getBufferSize() const { return buffer_size_; }
	
	    /// Changes the size of the internal buffer
	    void updateBuffer()
	    {
	        exp_.resize(buffer_size_);
	        buffer2scan_.resize(buffer_size_);
	    }
	
	    /// See std::vector documentation.
	    /// Note: the internal vector is of size 100 by default but this
	    /// function returns the actual number of scans stored so far
	    size_type size() const { return scan2buffer_.size(); }
			
	    /// empties buffer and removes temporary file
	    void clear()
	    {
	
	        scan_location_.clear();
	        buffer_index_    = 0;
	        scan2buffer_.clear();
	        buffer2scan_.clear();
	        exp_.clear();
	        exp_.resize(buffer_size_);
	
	        // generate new name for temp file
	        std::remove( file_name_ .c_str());
	        file_name_ = String("msexp_") + std::rand();
	
	    }
	
	    /// Returns the number of data points in the buffer (not scans)
	    UnsignedInt getSize() const { return nr_dpoints_; }
	
	    /// Same effect as updateBuffer()
	    void resize(UnsignedInt new_size) 
			{ 
					exp_.resize(new_size);
	        buffer2scan_.resize(buffer_size_);
	    }
	
	    /// See std::vector documentation.
	    void reserve(size_type n) {	exp_.reserve(n); }
	
	    /// returns a const reference to the sample description
	    const Sample& getSample() const { return exp_.getSample(); }
	    /// returns a mutable reference to the sample description
	    Sample& getSample() { return exp_.getSample();  }
	    /// sets the sample description 
			void setSample(const Sample& sample) { exp_.setSample(sample); }
			
	    /// returns a const reference to the source date file
	    const SourceFile& getSourceFile() const { return exp_.getSourceFile(); }
	    /// returns a mutable reference to the source date file
	    SourceFile& getSourceFile() { return exp_.getSourceFile(); }
	    /// sets the source date file 
			void setSourceFile(const SourceFile& source_file) { exp_.setSourceFile(source_file); }
	
	    /// returns a const reference to the list of contact persons
	    const std::vector<ContactPerson>& getContacts() const { return exp_.getContacts(); }
	    /// returns a mutable reference to the list of contact persons
	    std::vector<ContactPerson>& getContacts() { return exp_.getContacts(); }
	    /// sets the list of contact persons
	    void setContacts(const std::vector<ContactPerson>& contacts) { return exp_.setContacts(contacts); }
	
	    /// returns a const reference to the MS instrument description
	    const Instrument& getInstrument() const { return exp_.getInstrument(); }
	    /// returns a mutable reference to the MS instrument description
	    Instrument& getInstrument() { return exp_.getInstrument(); }
	    /// sets the MS instrument description
	    void setInstrument(const Instrument& instrument) { exp_.setInstrument(instrument); }
	
	    /// returns a const reference to the software used for processing
	    const Software& getSoftware() const { return exp_.getSoftware(); }
	    /// returns a mutable reference to the software used for processing
	    Software& getSoftware() { return exp_.getSoftware(); }
	    /// sets the software used for processing
	    void setSoftware(const Software& software) { exp_.setSoftware(software); }
	
	    /// returns a const reference to the description of the applied processing
	    const ProcessingMethod& getProcessingMethod() const { return exp_.getProcessingMethod(); }
	    /// returns a mutable reference to the description of the applied processing
	    ProcessingMethod& getProcessingMethod() { return exp_.getProcessingMethod(); }
	    /// sets the description of the applied processing
	    void setProcessingMethod(const ProcessingMethod& processing_method) { exp_.setProcessingMethod(processing_method); }
	
	    /// returns a const reference to the description of the HPLC run
	    const HPLC& getHPLC() const { return exp_.getHPLC(); }
	    /// returns a mutable reference to the description of the HPLC run
	    HPLC& getHPLC() { return exp_.getHPLC(); }
	    /// sets the description of the HPLC run
	    void setHPLC(const HPLC& hplc) { exp_.setHPLC(hplc); }
	
	    /// returns the experiment type
	    ExperimentType& getType() const { return exp_.getType(); }
	    /// sets the experiment type
	    void setType(ExperimentType type) { exp_.setType(type); }
	
	    /// returns the date the experiment was performed
	    const Date& getDate() const { return exp_.getDate(); }
	    /// sets the date the experiment was performed
	    void setDate(const Date& date) {	exp_.setDate(date); }
	
	    /// deletes the temporary file (or what did you expect ?)
	    void deleteTempFile_()	{ std::remove( file_name_ .c_str()); }
	
	    /// resets the internal data
	    void reset()
	    {
	        clear(); //remove data
	        exp_ = ExperimentalSettings(); //remove meta data
	        RangeManagerType::clearRanges(); // clear RangeManager
	    }
			
	    /// returns the meta information of this experiment (mutable access)
	    ExperimentalSettings& getExperimentalSettings() { return ( exp_); }
			/// returns the meta information of this experiment (const access)
			const ExperimentalSettings& getExperimentalSettings() const { return  exp_; }
	
			/// Returns an area iterator for @p area
			AIterator areaBegin(const AreaType& area)
			{
				return AIterator(RTBegin(area.minX()), RTEnd(area.maxX()), area.minY(), area.maxY());
			}
			
			/// Returns an invalid area iterator marking as the end iterator
			AIterator areaEnd()
			{
				return AIterator();
			}
			
			/// Returns an immutable area iterator for @p area
			AConstIterator areaBegin(const AreaType& area) const
			{
				return AConstIterator(RTBegin(area.minX()), RTEnd(area.maxX()), area.minY(), area.maxY());
			}
			
			/// Returns an immutable invalid area iterator marking as the end iterator
			AConstIterator areaEnd() const
			{
				return AConstIterator();
			}
	
	protected:
	    /// size of the internal buffer
	    UnsignedInt buffer_size_;
	
	    /// stores the offset of each scan on the hard disk
	    mutable std::vector<off_t> scan_location_;
	
	    /// number of scans added so far
	    mutable UnsignedInt current_scan_;
	
	    /// index in buffer
	    mutable UnsignedInt buffer_index_;
	
	    /// Maps scan index to index in buffer
	    mutable std::vector<size_type> scan2buffer_;
	
	    /// Maps buffer index to scan number
	    mutable std::vector<size_type> buffer2scan_;
	
	    /// Size of all scans
	    mutable std::vector<size_type> scan_sizes_;
	
	    /// Name of the temporary file to store the peak data
	    String file_name_;
	
	    /// The internal ms experiment instance.
	    mutable ExperimentType exp_;
	
	    /// File descriptor for temporary file
	    mutable FILE * pFile_;
	
	    /// The number of data points (peaks) in spectra of all MS levels (!)
	    UnsignedInt nr_dpoints_;
			
			/// MS levels of the data
	    std::vector<UnsignedInt> ms_levels_;
			
	    /// reads a scan from the temp file and stores it in the buffer
	    void storeInBuffer_(const size_type& n)
	    {
	        // check if buffer is full
	        if (buffer_index_ < buffer_size_)
	        {
	            //             std::cout << "buffer is not full, inserting scan at " << buffer_index_ << std::endl;
	            // 			std::cout << scan2buffer_.size() << "  " << buffer2scan_.size() << std::endl;
	
	            // test if we already wrote at this buffer position
	            if (current_scan_ > buffer_size_)
	            {
	                // yes => store scan at current buffer position and then overwrite
	                writeScan_(buffer2scan_[buffer_index_], exp_[buffer_index_] );
	            }
	
	            readScan_(n,exp_[buffer_index_]);
	            scan2buffer_[n] = buffer_index_;
	            buffer2scan_[buffer_index_++] = n;
	        }
	        else
	        {
	            // 			std::cout << "buffer is full, inserting scan at first position " << std::endl;
	            // buffer is full, therefore we overwrite the first entry
	            buffer_index_ = 0;
	
	            // check if size of buffer is set to zero
	            if (buffer_size_ > 0 )
	            {
	                writeScan_(buffer2scan_[buffer_index_], exp_[buffer_index_] );
	                readScan_(n,exp_[buffer_index_]);
	                scan2buffer_[n]                      = buffer_index_;
	                buffer2scan_[buffer_index_++] = n;
	            }
	            else // buffer size is set to zero
	            {
	                throw Exception::OutOfRange(__FILE__, __LINE__,"MSExperimentExtern::storeInBuffer_()");
	            }
	        }
	    }
			
			/// open the temporary file in mode @p mode
			off_t openFile_(const  off_t& pos, const char* mode) const
			{
					pFile_ = fopen(file_name_.c_str(),mode);
					
					if ( ftello(pFile_) < 0)
	        {
	            std::cout << "MSExperimentExtern:: Error determining writing position!" << std::endl;
	            std::cout << "Error code: " << errno << std::endl;
	            if (errno == EOVERFLOW)
	            {
	                std::cout << "An overflow of the position index was encountered." << std::endl;
	                std::cout << "Try re-compiling this class using -D_FILE_OFFSET_BITS=64"  << std::endl;
	                std::cout << "e.g. you might need to enable large file support for OpenMS since the temporary" << std::endl;
	                std::cout << "file became too large." << std::endl;
	                throw Exception::IndexOverflow(__FILE__, __LINE__,"MSExperimentExtern::writeScan_()",pos,sizeof(off_t));
	            }
	
	        }
					
					if (pos > 0)
					{
						//std::cout << "Setting stream to position : " << pos << std::endl;
						//fseeko(pFile_, pos,SEEK_SET);
						if ( fseeko(pFile_,pos,SEEK_SET) != 0)
	        	{
	            std::cout << "MSExperimentExtern:: Error determining reading position!" << std::endl;
	            std::cout << "Error code: " << errno << std::endl;
	            if (errno == EOVERFLOW)
	            {
	                std::cout << "An overflow of the position index was encountered." << std::endl;
	                std::cout << "Try re-compiling this class using -D_FILE_OFFSET_BITS=64"  << std::endl;
	                std::cout << "e.g. you might need to enable large file support for OpenMS since the temporary" << std::endl;
	                std::cout << "file became too large." << std::endl;
	                throw Exception::IndexOverflow(__FILE__, __LINE__,"MSExperimentExtern::readScan_()",pos,sizeof(off_t));
	            }
	        	}
					}
					
					// return current stream position (should be identical to pos by now) 
					return ftello(pFile_);		
			}
	
			/// Stores spectrum with number @p n in buffer
	    void storeInBuffer_(const size_type& n) const
	    {
	        // 		std::cout << "storeInBuffer_ :: spectra at " << n << " is not in buffer. " << std::endl;
	
	        // check if buffer is full
	        if (buffer_index_ < buffer_size_)
	        {
	            //             std::cout << "buffer is not full, inserting scan at " << buffer_index_ << std::endl;
	            // 			std::cout << scan2buffer_.size() << "  " << buffer2scan_.size() << std::endl;
	
	            // test if we already wrote at this buffer position
	            if (current_scan_ > buffer_size_)
	            {
	                // yes => store scan at current buffer position and then overwrite
	                writeScan_(buffer2scan_[buffer_index_], exp_[buffer_index_] );
	            }
	
	            readScan_(n,exp_[buffer_index_]);
	            scan2buffer_[n] = buffer_index_;
	            buffer2scan_[buffer_index_++] = n;
	        }
	        else
	        {
	            // 			std::cout << "buffer is full, inserting scan at first position " << std::endl;
	            // buffer is full, therefore we overwrite the first entry
	            buffer_index_ = 0;
	
	            // check if size of buffer is set to zero
	            if (buffer_size_ > 0 )
	            {
	                writeScan_(buffer2scan_[buffer_index_], exp_[buffer_index_] );
	                readScan_(n,exp_[buffer_index_]);
	                scan2buffer_[n]                      = buffer_index_;
	                buffer2scan_[buffer_index_++] = n;
	            }
	            else // buffer size is set to zero
	            {
	                throw Exception::OutOfRange(__FILE__, __LINE__,"MSExperimentExtern::storeInBuffer_()");
	            }
	        }
	    }
	
	    /// write spectrum to file
	    void writeScan_(const size_type& index, const SpectrumType& spec) const
	    {
					//std::cout << "Writing scan " << index << std::endl;
	// 				pFile_ = fopen(file_name_.c_str(),"a");
	        CoordinateType rt  = spec.getRetentionTime();
					UnsignedInt  mslvl = spec.getMSLevel();
					
	//         off_t pos = 0.0;
	        // determine position in file and store it
	//         if ( ( pos = ftello(pFile_) ) < 0)
	//         {
	//             std::cout << "MSExperimentExtern:: Error determining writing position!" << std::endl;
	//             std::cout << "Error code: " << errno << std::endl;
	//             if (errno == EOVERFLOW)
	//             {
	//                 std::cout << "An overflow of the position index was encountered." << std::endl;
	//                 std::cout << "Try re-compiling this class using -D_FILE_OFFSET_BITS=64"  << std::endl;
	//                 std::cout << "e.g. you might need to enable large file support for OpenMS since the temporary" << std::endl;
	//                 std::cout << "file became too large." << std::endl;
	//                 throw Exception::IndexOverflow(__FILE__, __LINE__,"MSExperimentExtern::writeScan_()",pos,sizeof(off_t));
	//             }
	// 
	//         }
	
	        // test if this scan was already written and store its offset
	        if (index >= scan_sizes_.size() )
	        {
							// scan was not written yet => append writing position
	            scan_location_.push_back( openFile_(0,"ab") );			
							scan_sizes_.push_back(spec.getContainer().size());			
	        }
	        else
	        {
	            // scan has already been written, check if size has changed
	            if (scan_sizes_[index] == spec.size() )
	            {
	                // write at old position
	               // std::cout << "Size has not changed. Overwriting" << std::endl;
									off_t pos = scan_location_[index];
									openFile_(pos,"rb+");
	                //fseeko(pFile_,pos,SEEK_SET);
	            }
	            else
	            {
	                // size has changed, forget old position and append
								//	std::cout << "Size has changed. Old "  << scan_sizes_[index] << " new " << spec.size() << std::endl;
									scan_location_[index] = openFile_(0,"a");
									scan_sizes_[index] = spec.getContainer().size();	// update scan size
	//                 scan_location_[index] = pos;
	            }
	
	        }
	
	       // std::cout << "writeScan_: writing scan " << index << " at " << ftello(pFile_) << std::endl;
				//	std::cout << "Writing rt: " << rt << " and mslvl " << mslvl << std::endl;
					// store retention time and ms level first
	//         fwrite(&rt,sizeof(CoordinateType),1,pFile_);
					if ( fwrite(&rt,sizeof(rt),1,pFile_) == 0 )
					{
					  std::cout << "Error writing retention time" << std::endl;
						throw Exception::IndexOverflow(__FILE__, __LINE__,"MSExperimentExtern::readScan_()",1,sizeof(off_t));
					}
					//fwrite(&mslvl,sizeof(UnsignedInt),1,pFile_);
					if (fwrite(&mslvl,sizeof(mslvl),1,pFile_) == 0)
					{
						std::cout << "Error writing ms level" << std::endl;
						throw Exception::IndexOverflow(__FILE__, __LINE__,"MSExperimentExtern::readScan_()",1,sizeof(off_t));
					}			
					
	//         size_t sizeof_peak =  sizeof(PeakType);
	//				std::cout << "Writing " << spec.getContainer().size() << " peaks." << std::endl;
	// 				std::cout << "Peak size: " << sizeof_peak << std::endl;
	//         for (typename ContainerType::const_iterator cit = spec.getContainer().begin();
	//                 cit != spec.getContainer().end(); ++cit)
	//         {
	// // 						std::cout << "Writing peak: "  << *cit << std::endl;
	//             if ( fwrite(&(*cit),sizeof_peak,1,pFile_) != sizeof_peak )
	// 						{
	// 							std::cout << "Error writing peak data" << std::endl;
	// 							throw Exception::IndexOverflow(__FILE__, __LINE__,"MSExperimentExtern::readScan_()",1,sizeof(off_t));
	// 						}
	//         }
					
					for (UnsignedInt i = 0; i < spec.getContainer().size();++i)
	        {
							PeakType p = spec.getContainer()[i];
	// 						std::cout << "Writing peak: "  << *cit << std::endl;
	//             if ( fwrite(&(spec.getContainer()[i]),sizeof(spec.getContainer()[i]),1,pFile_) != sizeof(spec.getContainer()[i]) )
							if ( fwrite(&p,sizeof(p),1,pFile_) != 1 )
							{
								std::cout << "Error writing peak data" << std::endl;
								throw Exception::IndexOverflow(__FILE__, __LINE__,"MSExperimentExtern::readScan_()",1,sizeof(off_t));
							}
	        }
	
	        fclose(pFile_);
	
	        // test if this scan was already written and store its size
	//         if (index >= scan_sizes_.size() )
	//         {
	//             scan_sizes_.push_back(spec.getContainer().size());
	//         }
	//         else
	//         {
	//             scan_sizes_[index] = spec.getContainer().size();
	//         }
	
	    } // end of write(spectrum)
	
	    /// Reads a spectrum from the temporary file
	    void readScan_(const size_type& index, SpectrumType& spec)  const
	    {
	        //pFile_ = fopen(file_name_.c_str(),"r");
					        
					// set stream to starting point of last writing action
	// 				std::cout << "Reading from " << file_name_ << std::endl;
	// 				std::cout << "Retrieving reading offset: " << scan_location_.size() << " " << index << std::endl;
	        off_t pos = scan_location_[index];
					
					/*off_t p2 =*/ openFile_(pos,"rb");
					
	//				std::cout << "readScan_: reading scan " << index << " from " << p2 << std::endl;
	//         if ( fseeko(pFile_,pos,SEEK_SET) != 0)
	//         {
	//             std::cout << "MSExperimentExtern:: Error determining reading position!" << std::endl;
	//             std::cout << "Error code: " << errno << std::endl;
	//             if (errno == EOVERFLOW)
	//             {
	//                 std::cout << "An overflow of the position index was encountered." << std::endl;
	//                 std::cout << "Try re-compiling this class using -D_FILE_OFFSET_BITS=64"  << std::endl;
	//                 std::cout << "e.g. you might need to enable large file support for OpenMS since the temporary" << std::endl;
	//                 std::cout << "file became too large." << std::endl;
	//                 throw Exception::IndexOverflow(__FILE__, __LINE__,"MSExperimentExtern::readScan_()",pos,sizeof(off_t));
	//             }
	//         }
	
	        
	        // read retention time
	        CoordinateType rt = 0;
					UnsignedInt mslvl = 0;
	        fread(&rt,sizeof(CoordinateType),1,pFile_);
		//			std::cout << "Reading rt: " << rt << std::endl;
					fread(&mslvl,sizeof(UnsignedInt),1,pFile_);
		//			std::cout << "Reading ms level: " << mslvl << std::endl;
		
	        spec.setRetentionTime(rt);
					spec.setMSLevel(mslvl);
					unsigned int nr_peaks = scan_sizes_[index];
	   //     std::cout << "Reading peaks: " << nr_peaks << std::endl;
		//			std::cout << "Stream position: " << p2 << std::endl;
		//			std::cout << "Desired position: " << pos << std::endl;
	        spec.getContainer().clear();
	        spec.resize(nr_peaks);
	
					size_t sizeof_peak =  sizeof(PeakType);
	
	        //read coordinates of each peak
	        for (typename SpectrumType::Iterator piter = spec.begin(); piter != spec.end(); ++piter)
	        {
	            if (fread(&(*piter),sizeof_peak,1,pFile_) == 0) 
							{
	                std::cout << "Error reading peak data" << std::endl;
									if ( feof(pFile_) ) std::cout << "End of file was reached. " << std::endl;
									
									throw Exception::IndexOverflow(__FILE__, __LINE__,"MSExperimentExtern::readScan_()",pos,sizeof(off_t));
							}
	        }
	        // 		std::cout << "Done."<< std::endl;
	        fclose(pFile_);
	
	    }	// end of read const
	
	    /// copies the content of the tempory file
	    void copyTmpFile_(String source)
	    {
					std::ifstream input(source.c_str());
	   			std::ofstream output(file_name_.c_str()); 
					
					output << input.rdbuf();
			
					input.close();
					output.close();
			}
	
	
	};  // end of class MSExperimentExtern

} // namespace OpenMS

#endif // OPENMS_KERNEL_MSEXPERIMENTEXTERN_H
