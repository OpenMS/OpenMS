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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MSEXPERIMENTEXTERN_H
#define OPENMS_KERNEL_MSEXPERIMENTEXTERN_H

#define _FILE_OFFSET_BITS 64

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/SYSTEM/StopWatch.h>

#include<vector>
#include<algorithm>
#include<limits>

#include <iostream>

#include <stdio.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/errno.h>

namespace OpenMS
{
/**
	@brief Representation of a mass spectrometry experiment using an external datastructure to store large data sets.
	
	This data structures has the same interface as MSExperiment but uses a ring buffer and stores only a subset of
	all scans in the RAM. Scans are dynamically written to the hard disk and re-loaded if needed.
	
	NOTE: If your LC-MS map is really large, you might want to compile this class with LFS (large file support) such
	that Linux / C can access files > 2 GB. In this case, you will need to comple OpenMS with -D_FILE_OFFSET_BITS = 64.
		
	TODO: Implement range informations. (ole)
	
	TODO: Remember last scan accessed in getPeak(unsigned int) and getPeakRt(unsigned int) (ole)
			
	@ingroup Kernel
**/
template < typename PeakT = DPeak<1> >
class MSExperimentExtern
	: public RangeManager<2, typename PeakT::TraitsType>
	{
	public:

    ///ConstIterator
    template <class IteratorPeakT>
    class MSExperimentExternConstIterator
    {

        friend class MSExperimentExtern;

    public:
        /**	@name Type definitions */
        //@{
        typedef IteratorPeakT IteratorPeakType;
        typedef MSSpectrum<IteratorPeakType> value_type;
        typedef MSExperimentExtern<IteratorPeakType> ExperimentType;
        typedef const value_type& reference;
        typedef const value_type* pointer;
        typedef std::random_access_iterator_tag iterator_category;
        typedef unsigned int difference_type;
        //@}

        MSExperimentExternConstIterator()
                : exp_(0), position_(0)
        {}

        MSExperimentExternConstIterator(const ExperimentType * exp , unsigned int pos)
        {
            exp_       =  (ExperimentType* ) exp;
            position_ = pos;
        }

        MSExperimentExternConstIterator(ExperimentType * exp , unsigned int pos)
                : exp_(exp), position_(pos)
        {}

        MSExperimentExternConstIterator(const ExperimentType& source)
                : exp_(source.exp_), position_(source.position_)
        {}

        ~MSExperimentExternConstIterator()
        {}

        MSExperimentExternConstIterator& operator = (const MSExperimentExternConstIterator& rhs)
        {
            if (this==&rhs)
                return *this;

            exp_      = rhs.exp_;
            position_= rhs.position_;

            return *this;
        }

        bool operator < (const MSExperimentExternConstIterator& it) const
        {
            return position_ < it.position_;
        }

        bool operator > (const MSExperimentExternConstIterator& it) const
        {
            return position_ > it.position_;
        }

        bool operator <= (const MSExperimentExternConstIterator& it) const
        {
            return (position_ < it.position_ || position_ == it.position_);
        }

        bool operator >= (const MSExperimentExternConstIterator& it) const
        {
            return (position_ > it.position_ || position_ == it.position_);
        }

        bool operator == (const MSExperimentExternConstIterator& it) const
        {
            return position_ == it.position_ && exp_ == it.exp_;
        }

        bool operator != (const MSExperimentExternConstIterator& it) const
        {
            return position_ != it.position_ || exp_ != it.exp_;
        }

        MSExperimentExternConstIterator& operator ++ ()
        {
            position_ += 1;
            return *this;
        }

        MSExperimentExternConstIterator operator ++ (int)
        {
            MSExperimentExternConstIterator tmp(*this);
            ++(*this);
            return tmp;
        }

        MSExperimentExternConstIterator& operator -- ()
        {
            position_ -= 1;
            return *this;
        }

        MSExperimentExternConstIterator operator -- (int)
        {
            MSExperimentExternConstIterator tmp(*this);
            --(*this);
            return tmp;
        }

        MSExperimentExternConstIterator operator - (difference_type n) const
        {
            MSExperimentExternConstIterator tmp(*this);
            tmp.position_ -= n;
            return tmp;
        }

        MSExperimentExternConstIterator operator + (difference_type n) const
        {
            MSExperimentExternConstIterator tmp(*this);
            tmp.position_ += n;
            return tmp;
        }

        MSExperimentExternConstIterator& operator += (difference_type n)
        {
            position_ += n;
            return *this;
        }

        MSExperimentExternConstIterator& operator -= (difference_type n)
        {
            position_ -= n;
            return *this;
        }

        friend difference_type operator - ( const MSExperimentExternConstIterator& i1, const MSExperimentExternConstIterator& i2 )
        {
            return (i1.position_ - i2.position_);
        }

        friend MSExperimentExternConstIterator operator + ( difference_type n, const MSExperimentExternConstIterator& i )
        {
            MSExperimentExternConstIterator tmp(i);
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

    protected:

        ExperimentType * exp_;
        unsigned int position_;
    };


    /// Mutable iterator
    template <class IteratorPeakT>
	class MSExperimentExternIterator : public MSExperimentExternConstIterator<IteratorPeakT>
    {
        friend class MSExperimentExtern;

    public:
        typedef IteratorPeakT IteratorPeakType;
        typedef typename MSExperimentExternConstIterator<IteratorPeakType>::value_type& reference;
        typedef typename MSExperimentExternConstIterator<IteratorPeakType>::value_type* pointer;
        typedef typename MSExperimentExternConstIterator<IteratorPeakType>::ExperimentType ExperimentType;
        typedef typename MSExperimentExternConstIterator<IteratorPeakType>::difference_type difference_type;

        using MSExperimentExternConstIterator<IteratorPeakType>::exp_;
        using MSExperimentExternConstIterator<IteratorPeakType>::position_;


        MSExperimentExternIterator()
                : MSExperimentExternConstIterator<IteratorPeakType>()
        {}

        MSExperimentExternIterator(ExperimentType * exp, unsigned int position)
                : MSExperimentExternConstIterator<IteratorPeakType>(exp,position)
        {}

        MSExperimentExternIterator(const MSExperimentExternIterator<IteratorPeakType>& it)
                : MSExperimentExternConstIterator<IteratorPeakType>(it)
        {}

        ~MSExperimentExternIterator()
        {}

        reference operator * ()
        {
            return (*exp_)[position_];
        }

        pointer operator -> ()
        {
            return &((*exp_)[position_]);
        }

        const pointer operator -> () const
        {
            return &((*exp_)[position_]);
        }

        typename MSExperimentExternIterator::reference operator [] (difference_type n)
        {
            return ((*this)+n);
        }

        MSExperimentExternIterator& operator ++ ()
        {
            MSExperimentExternConstIterator<IteratorPeakType>::position_+=1;
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
            MSExperimentExternConstIterator<IteratorPeakType>::position_-=1;
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

        friend MSExperimentExternIterator operator + (difference_type n, const MSExperimentExternIterator& i )
        {
            MSExperimentExternIterator tmp(i);
            tmp.position_ += n;
            return tmp;
        }

        MSExperimentExternIterator& operator += (difference_type n)
        {
            MSExperimentExternConstIterator<IteratorPeakType>::position_ += n;
            return *this;
        }

        MSExperimentExternIterator& operator -= (difference_type n)
        {
            MSExperimentExternConstIterator<IteratorPeakType>::position_ -= n;
            return *this;
        }

        friend void swap(MSExperimentExternIterator& i1, MSExperimentExternIterator& i2)
        {
            unsigned int tmp = i1.position_;
            i1.position_ = i2.position_;
            i2.position_ = tmp;
        }

    protected:

    } ;	// end of class MSExperimentExternIterator
	
	    /**
    	@brief Adaptor class for linear iterator on objects of DPeak<1>  
    	
    	This iterator allows us to move through the data structure in a linear
    	manner i.e. we don't need to jump to the next spectrum manually.
    	
    */
	template <class IteratorPeakT >
 	class PeakIterator : public std::iterator<std::bidirectional_iterator_tag,  IteratorPeakT>
    {
        typedef double CoordinateType;
        typedef IteratorPeakT IteratorPeakType;
		
    public:

        /// Default constructor
        PeakIterator()
                : peak_index_(), rt_(), scan_index_(), exp_()
        {
        }

        /// Constructor
        PeakIterator(unsigned int pind, CoordinateType & co, unsigned int sind, MSExperimentExtern<IteratorPeakType>& exp)
                : peak_index_(pind), rt_(co), scan_index_(sind), exp_(&exp)
        {
        }

        /// Destructor
        ~PeakIterator()
        {
        }

        /// Copy constructor
        PeakIterator(const PeakIterator& rhs)
                : peak_index_(rhs.peak_index_), rt_(rhs.rt_),
                scan_index_(rhs.scan_index_), exp_(rhs.exp_)
        { }

        /// Assignment operator
        PeakIterator& operator=(const PeakIterator& rhs)
        {
            if (&rhs == this)
                return *this;

            peak_index_    = rhs.peak_index_;
            rt_                   = rhs.rt_;
            scan_index_    = rhs.scan_index_;
            exp_               = rhs.exp_;

            return (*this);
        }

        /// Test for equality
        friend bool operator==(const PeakIterator& lhs, const PeakIterator& rhs)
        {
            return ( lhs.peak_index_     == rhs.peak_index_ &&
                     				 lhs.rt_       == rhs.rt_       &&
                     	    lhs.scan_index_  == rhs.scan_index_ );
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
            if ( peak_index_ ==   (*exp_)[scan_index_].size() && scan_index_ !=  ( (*exp_).size() - 1) )
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
					std::cout << "In first scan and moving backwards ! " << std::endl;
					return (*this);
				}
				--scan_index_;
                peak_index_  =  ( (*exp_)[scan_index_].size() -1) ;
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
        CoordinateType& getRt()  { return rt_; }
        /// Returns the current retention time (not mutable)
        const CoordinateType& getRt() const  {  return rt_; }
        ///Returns the current retention time (mutable)
        unsigned int& getPeakIndex() { return peak_index_; }
        /// Returns the current retention time (not mutable)
        const unsigned int& getPeakIndex() const { return peak_index_; }
		/// Mutable access to current scan
		unsigned int& getScanIndex() { return scan_index_;}
		/// Const access to  current scan
		const unsigned int& getScanIndex() const { return scan_index_; } 
        //@}

    private:
        /// Points to the current peak
        unsigned int peak_index_;
        /// Retention time of the current spectrum
        CoordinateType rt_;
        /// Index of the current spectrum
        unsigned int scan_index_;
        /// Pointer to the experiment
        MSExperimentExtern<IteratorPeakType> * exp_;

    }
    ; // end of inner class PeakIterator
    
	typedef PeakT PeakType;
    typedef typename PeakT::IntensityType IntensityType;
	typedef typename PeakT::PositionType PositionType;
	typedef MSSpectrum<PeakT> SpectrumType;
    typedef typename SpectrumType::ContainerType ContainerType;
    typedef typename PeakType::CoordinateType CoordinateType;
    typedef MSExperiment<PeakType> ExperimentType;

    typedef typename PeakType::TraitsType TraitsType;
    typedef RangeManager<2, TraitsType> RangeManagerType;	

    typedef MSExperimentExternIterator<PeakType> Iterator;
    typedef MSExperimentExternConstIterator<PeakType> ConstIterator;
	typedef PeakIterator<PeakType> PIterator;
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
			nr_dpoints_(0)
    {
        file_name_ = "msexp_" + String(std::rand());
		exp_.resize(buffer_size_);
		buffer2scan_.resize(buffer_size_);
    }

    /// Copy constructor: copies the content of the temporary file as well (slow !)
    MSExperimentExtern(const MSExperimentExtern& source)
            : buffer_size_(source.buffer_size_),
            scan_location_(source.scan_location_), current_scan_(source.current_scan_),
            buffer_index_(source.buffer_index_),  scan2buffer_(source.scan2buffer_),
            buffer2scan_(source.buffer2scan_), exp_(source.exp_),
			scan_sizes_(source.scan_sizes_), nr_dpoints_(source.nr_dpoints_)
	 {
	 	// genarete new temp file and copy the old one
        file_name_ = "msexp_" + String(std::rand());
		copyTmpFile__(source.file_name_);
    }

    /// Destructor
    ~MSExperimentExtern()
    {
        // delete temporary file
        std::remove ( file_name_ .c_str());
	}

    /// Assignment operator
    MSExperimentExtern & operator= (const MSExperimentExtern& source)
    {
        if (&source == this)
            return *this;

        buffer_size_      = source.buffer_size_;
        scan_location_ = source.scan_location_;
        buffer_index_    = source.buffer_index_;
		current_scan_  = source.current_scan_;
        scan2buffer_    = source.scan2buffer_;
        buffer2scan_    = source.buffer2scan_;
		scan_sizes_   = source.scan_sizes_;
        exp_				  = source.exp_;
		nr_dpoints_    = source.nr_dpoints_;

        // generate new name for temp file
		std::remove ( file_name_ .c_str());
        file_name_ = "msexp_" + String(std::rand());
		// and copy the old one
        copyTmpFile__(source.file_name_);

        return *this;
    }

    /// Equality operator
    bool operator== (const MSExperimentExtern& rhs) const
    {
        return (buffer_size_   == rhs.buffer_size_ &&
                scan_location_ == rhs.scan_location_ &&
                buffer_index_    == rhs.buffer_index_ &&
                scan2buffer_    == rhs.scan2buffer_ &&
                buffer2scan_    == rhs.buffer2scan_ &&
                exp_				  == rhs.exp_             &&
				scan_sizes_   == rhs.scan_sizes_ &&
				pFile_			  == rhs.pFile_);
    }

    /// Equality operator
    bool operator!= (const MSExperimentExtern& rhs) const
    {
        return !(operator==(rhs));
    }

    template <class Container>
    void get2DData(Container& cont) const
    {
        exp_.get2DData(cont);
    }

    template <class Container>
    void set2DData(Container& cont)
    {
        exp_.set2DData(cont);
    }

	/// Updates the range information 
    void updateRanges()
    {
       // exp_.updateRanges();
	   	   
	   nr_dpoints_ = 0;
	   spectra_lengths_.clear();
       	   
	   for (unsigned int i =0; i< scan2buffer_.size(); ++i)
	   {
	   		nr_dpoints_ += (*this)[i].size();
// 	   		std::cout << "Accumulated size: " << nr_dpoints_ << std::endl;
			spectra_lengths_.push_back(nr_dpoints_);
	   }
	   
    }
	
	/// Returns the minimum position
	const PositionType& getMin() const	
	{ 
		return exp_.getMin(); 
	}
	
	/// Returns the maximum position
	const PositionType& getMax() const 
	{ 
		return exp_.getMax(); 
	}
	
	/// Returns the minimum intensity
	const IntensityType getMinInt() const	
	{ 
		return exp_.getMin()[0]; 
	}
	
	/// Returns the maximum intensity
	const IntensityType getMaxInt() const 
	{ 
	 	return exp_.getMax()[0]; 
	}

    void sortSpectra(bool sort_mz)
    {
        exp_.sortSpectra(sort_mz);
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
    Iterator begin()
    {
        return Iterator(this,(unsigned int)0);
    }

    /// See std::vector documentation.
    Iterator end()
    {
        return Iterator(this,this->size());
    }

    /// See std::vector documentation.
    ConstIterator begin() const
    {
        return ConstIterator(this,(unsigned int)0);
    }

    /// See std::vector documentation.
    ConstIterator end() const
    {
        return ConstIterator(this,this->size());
    }

    /// See std::vector documentation.
    ReverseIterator rbegin()
    {
        return ReverseIterator(end());
       
    }

    /// See std::vector documentation.
    ReverseIterator rend()
    {
        return ReverseIterator(begin());
       
    }

    /// See std::vector documentation.
    ConstReverseIterator rbegin() const
    {
        return ConstReverseIterator(end());
    }

    /// See std::vector documentation.
    ConstReverseIterator rend() const
    {
        return ConstReverseIterator(begin());
    }
   	
	/// See std::vector documentation.
	reference back()
	{
		return  this->at( (scan2buffer_.size() -1) ) ;
	}

    void push_back(const SpectrumType& spec)
    {
        std::cout << "Inserting scan " << current_scan_ << std::endl;
        std::cout << "buffer capacity: " << buffer_size_ << " buffer index: " << buffer_index_ << " buffer size: " << exp_.size() << std::endl;
		
        if (buffer_index_ < buffer_size_)
        {
			std::cout << "Writing in buffer at pos: " << buffer_index_ << std::endl;
		
			// test if we already wrote at this buffer position
			if (current_scan_ > buffer_size_)
			{
				// yes => store scan at current buffer position and then overwrite
				writeScan(current_scan_, exp_[buffer_index_] );
			}
			
			exp_[buffer_index_] = spec;
           	scan2buffer_.push_back(buffer_index_);
            buffer2scan_[buffer_index_++] = current_scan_++;			
        }
        else
        {
//             std::cout << "Buffer full. Overwriting buffer at pos 0."   << std::endl;
			
			buffer_index_ = 0; 																		// reset buffer index
			writeScan(buffer2scan_[buffer_index_],  exp_[buffer_index_] ); 		// write content of buffer to temp. file
			exp_[buffer_index_] = spec;														// store new spectrum 
			
            scan2buffer_.push_back(buffer_index_);
            buffer2scan_[buffer_index_++] = current_scan_++;
		}
		
// 		std::cout << "scan2buffer : " << scan2buffer_[ (current_scan_-1)] << std::endl;
// 		std::cout << "buffer2scan: " << buffer2scan_[ ( scan2buffer_[ (current_scan_-1)]  )] << std::endl;
		
        //writeScan(spec);
    }

    /// see std::vector (additionally test if scan is in buffer or needs to be read from temp file)
    reference operator[] (size_type n)
    {

//         std::cout << "operator[" << n << "]" << std::endl;
        // test if current scan is in buffer
        UnsignedInt b = scan2buffer_[n];
        if (buffer2scan_[b] != n)
            storeInBuffer(n);	// scan is not in buffer, needs to be read from file

        b = scan2buffer_[n];
        return exp_[b];
    }

    /// see std::vector (additionally test if scan is in buffer or needs to be read from temp file)
    const_reference operator[] (size_type n) const
    {
//         std::cout << "operator[" << n << "] const" << std::endl;
        // test if current scan is in buffer
        UnsignedInt b = scan2buffer_[n];
        if (buffer2scan_[b] != n)
		{
//             std::cout << "scan not in buffer." << std::endl;
			storeInBuffer(n);	// scan is not in buffer, needs to be read from file
		}
        b = scan2buffer_[n];
        return exp_[b];
    }

    /// see std::vector (additionally test if scan is in buffer or needs to be read from temp file)
    reference at(size_type n)
    {
        //std::cout << "at(" << n << ")" << std::endl;
        // test if current scan is in buffer
        UnsignedInt b = scan2buffer_[n];
        if (buffer2scan_[b] != n)
            storeInBuffer(n);	// scan is not in buffer, needs to be read from file

        b = scan2buffer_[n];
        return exp_.at(b);
    }

    /// see std::vector (additionally test if scan is in buffer or needs to be read from temp file)
    const_reference at(size_type n) const
    {
        //std::cout << "at(" << n << ") const" << std::endl;
        // test if current scan is in buffer
        UnsignedInt b = scan2buffer_[n];
        if (buffer2scan_[b] != n)
            storeInBuffer(n);	// scan is not in buffer, needs to be read from file

        b = scan2buffer_[n];
        return exp_.at(b);
    }

    /// Sets buffer size
    void setBufferSize(UnsignedInt sz)
    {
        buffer_size_ = sz;
    }
    /// Returns the buffer size (mutable)
    UnsignedInt& getBufferSize()
    {
        return buffer_size_;
    }
    /// Returns the buffer size (not mutable)
    const UnsignedInt getBufferSize() const
    {
        return buffer_size_;
    }

    /// Changes the size of the internal buffer
    void updateBuffer()
    {
        exp_.resize(buffer_size_);
		buffer2scan_.resize(buffer_size_);
    }

    /// See std::vector documentation.
    /// Note: the internal vector is of size 100 by default but this
    /// function returns the number of scans stored so far
    size_type size() const
    {
        return scan2buffer_.size();
    }

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
		std::remove ( file_name_ .c_str());
        file_name_ = "msexp_" + String(std::rand());		

    }

    /// See MSExperiment documentation.
    UnsignedInt getSize() const
    {
        return exp_.getSize();
    }
	
	/// Same effect as updateBuffer()
	void resize(UnsignedInt new_size) 
	{
		exp_.resize(new_size);
		buffer2scan_.resize(buffer_size_);
	}

    /// See std::vector documentation.
    void reserve(size_type n)
    {
        exp_.reserve(n);
    }

    /// returns a const reference to the sample description
    const Sample& getSample() const
    {
        return exp_.getSample();
    }
    /// returns a mutable reference to the sample description
    Sample& getSample()
    {
        return exp_.getSample();
    }
    /// sets the sample description
    void setSample(const Sample& sample)
    {
        exp_.setSample(sample);
    }

    /// returns a const reference to the source date file
    const SourceFile& getSourceFile() const
    {
        return exp_.getSourceFile();
    }
    /// returns a mutable reference to the source date file
    SourceFile& getSourceFile()
    {
        return exp_.getSourceFile();
    }
    /// sets the source date file
    void setSourceFile(const SourceFile& source_file)
    {
        exp_.setSourceFile(source_file);
    }

    /// returns a const reference to the list of contact persons
    const std::vector<ContactPerson>& getContacts() const
    {
        return exp_.getContacts();
    }
    /// returns a mutable reference to the list of contact persons
    std::vector<ContactPerson>& getContacts()
    {
        return exp_.getContacts();
    }
    /// sets the list of contact persons
    void setContacts(const std::vector<ContactPerson>& contacts)
    {
        return exp_.setContacts(contacts);
    }

    /// returns a const reference to the MS instrument description
    const Instrument& getInstrument() const
    {
        return exp_.getInstrument();
    }
    /// returns a mutable reference to the MS instrument description
    Instrument& getInstrument()
    {
        return exp_.getInstrument();
    }
    /// sets the MS instrument description
    void setInstrument(const Instrument& instrument)
    {
        exp_.setInstrument(instrument);
    }

    /// returns a const reference to the software used for processing
    const Software& getSoftware() const
    {
        return exp_.getSoftware();
    }
    /// returns a mutable reference to the software used for processing
    Software& getSoftware()
    {
        return exp_.getSoftware();
    }
    /// sets the software used for processing
    void setSoftware(const Software& software)
    {
        exp_.setSoftware(software);
    }

    /// returns a const reference to the description of the applied processing
    const ProcessingMethod& getProcessingMethod() const
    {
        return exp_.getProcessingMethod();
    }
    /// returns a mutable reference to the description of the applied processing
    ProcessingMethod& getProcessingMethod()
    {
        return exp_.getProcessingMethod();
    }
    /// sets the description of the applied processing
    void setProcessingMethod(const ProcessingMethod& processing_method)
    {
        exp_.setProcessingMethod(processing_method);
    }

    /// returns a const reference to the description of the HPLC run
    const HPLC& getHPLC() const
    {
        return exp_.getHPLC();
    }
    /// returns a mutable reference to the description of the HPLC run
    HPLC& getHPLC()
    {
        return exp_.getHPLC();
    }
    /// sets the description of the HPLC run
    void setHPLC(const HPLC& hplc)
    {
        exp_.setHPLC(hplc);
    }

    /// returns the experiment type
    ExperimentType& getType() const
    {
        return exp_.getType();
    }
    /// sets the experiment type
    void setType(ExperimentType type)
    {
        exp_.setType(type);
    }

    /// returns the date the experiment was performed
    const Date& getDate() const
    {
        return exp_.getDate();
    }
    /// sets the date the experiment was performed
    void setDate(const Date& date)
    {
        exp_.setDate(date);
    }

	/// deletes the temporary file (or what did you expect ?)
    void deleteTempFile_()
    {
        std::remove( file_name_ .c_str());
    }
		
	/// resets the internal data
    void reset()
    {
    	clear(); //remove data
    	exp_ = ExperimentalSettings(); //remove meta data
      	RangeManagerType::clearRanges(); // clear RangeManager
    }
	
	/// returns the meta information of this experiment
	const ExperimentalSettings getExperimentalSettings() const 
	{
		return ((ExperimentalSettings) exp_);
	}
	
	/// Mutable access to peak with index @p  
   	preference getPeak(UnsignedInt index) throw (Exception::IndexOverflow)
    {
        if (index > nr_dpoints_)
            throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, nr_dpoints_);

		std::vector<UnsignedInt>::iterator it = std::upper_bound(spectra_lengths_.begin(),spectra_lengths_.end(),index);
	
		unsigned int scan_index =  (it - spectra_lengths_.begin() );
		unsigned int peak_index = 0;
	
		if (scan_index == 0)
		{
			peak_index  =  index;
		}
		else
		{
			--it;
			peak_index = (index - *it);
		}		
		
		return (*this)[scan_index][peak_index];
    }
	
	/// const access to peak with index @p  
   	const_preference getPeak(UnsignedInt index) const throw (Exception::IndexOverflow) 
    {
        if (index > nr_dpoints_)
            throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, nr_dpoints_);

		std::vector<UnsignedInt>::iterator it = std::upper_bound(spectra_lengths_.begin(),spectra_lengths_.end(),index);
	
		unsigned int scan_index =  (it - spectra_lengths_.begin() );
		unsigned int peak_index = 0;
	
		if (scan_index == 0)
		{
			peak_index  =  index;
		}
		else
		{
			--it;
			peak_index = (index - *it);
		}		
		
		return (*this)[scan_index][peak_index];
    }
	
	/// Mutable access to retention time of peak @p
	CoordinateType& getPeakRt(UnsignedInt index) throw (Exception::IndexOverflow) 
	{
		if (index > nr_dpoints_)
            throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, nr_dpoints_);
			
		std::vector<unsigned int>::iterator it = std::upper_bound(spectra_lengths_.begin(),spectra_lengths_.end(),index);
		unsigned int scan_index =  (it - spectra_lengths_.begin() );
		return ((*this)[scan_index]).getRetentionTime();
		
	}
	
	/// Const access to retention time of peak @p
	const CoordinateType& getPeakRt(UnsignedInt index) const throw (Exception::IndexOverflow) 
	{
		if (index > nr_dpoints_)
            throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, nr_dpoints_);
			
		std::vector<unsigned int>::iterator it = std::upper_bound(spectra_lengths_.begin(),spectra_lengths_.end(),index);
		unsigned int scan_index =  (it - spectra_lengths_.begin() );
		return (*this)[scan_index].getRetentionTime();
		
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
	
	/// Sums of consecutive spectrum lengths 
	std::vector<UnsignedInt> spectra_lengths_; 
	
    /// reads a scan from the temp file and stores it in the buffer
    void storeInBuffer(const size_type& n)
    {
//         std::cout << "storeInBuffer :: spectra at " << n << " is not in buffer. " << std::endl;
        
		// check if buffer is full
        if (buffer_index_ < buffer_size_)
        {
//             std::cout << "buffer is not full, inserting scan at " << buffer_index_ << std::endl;
// 			std::cout << scan2buffer_.size() << "  " << buffer2scan_.size() << std::endl;
			
			// test if we already wrote at this buffer position
			if (current_scan_ > buffer_size_)
			{
				// yes => store scan at current buffer position and then overwrite
				writeScan(buffer2scan_[buffer_index_], exp_[buffer_index_] );
			}
			
			readScan(n,exp_[buffer_index_]);         	
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
				writeScan(buffer2scan_[buffer_index_], exp_[buffer_index_] );
               	readScan(n,exp_[buffer_index_]);
                scan2buffer_[n]                      = buffer_index_;
                buffer2scan_[buffer_index_++] = n;
            }
            else // buffer size is set to zero
            {
                throw Exception::OutOfRange(__FILE__, __LINE__,"MSExperimentExtern::storeInBuffer()");
            }
        }
    }

    void storeInBuffer(const size_type& n) const
    {
// 		std::cout << "storeInBuffer :: spectra at " << n << " is not in buffer. " << std::endl;
        
		// check if buffer is full
        if (buffer_index_ < buffer_size_)
        {
//             std::cout << "buffer is not full, inserting scan at " << buffer_index_ << std::endl;
// 			std::cout << scan2buffer_.size() << "  " << buffer2scan_.size() << std::endl;
			
			// test if we already wrote at this buffer position
			if (current_scan_ > buffer_size_)
			{
				// yes => store scan at current buffer position and then overwrite
				writeScan(buffer2scan_[buffer_index_], exp_[buffer_index_] );
			}
			
			readScan(n,exp_[buffer_index_]);         	
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
				writeScan(buffer2scan_[buffer_index_], exp_[buffer_index_] );
               	readScan(n,exp_[buffer_index_]);
                scan2buffer_[n]                      = buffer_index_;
                buffer2scan_[buffer_index_++] = n;
            }
            else // buffer size is set to zero
            {
                throw Exception::OutOfRange(__FILE__, __LINE__,"MSExperimentExtern::storeInBuffer()");
            }
        }
    }

    /// write spectrum to file
    void writeScan(const size_type& index, const SpectrumType& spec) const
    {
       	pFile_ = fopen(file_name_.c_str(),"a");
        CoordinateType rt = spec.getRetentionTime();
				
		off_t pos;
		// determine position in file and store it
		if ( ( pos = ftello(pFile_) ) < 0) 
		{
			std::cout << "MSExperimentExtern:: Error determining writing position!" << std::endl;
			std::cout << "Error code: " << errno << std::endl;	
			if (errno == EOVERFLOW)
			{
				std::cout << "An overflow of the position index was encountered." << std::endl;
				std::cout << "Try re-compiling this class using -D_FILE_OFFSET_BITS=64"  << std::endl;
				std::cout << "e.g. you might need to enable large file support for OpenMS since the temporary" << std::endl;
				std::cout << "file became too large." << std::endl;
				throw Exception::IndexOverflow(__FILE__, __LINE__,"MSExperimentExtern::writeScan()",pos,sizeof(off_t)); 
			}
			
		}
				
		// test if this scan was already written and store its offset
		if (index >= scan_sizes_.size() )
		{
        	scan_location_.push_back( pos );
		}
		else
		{
			// scan has already been written, check if size has changed			
			if (scan_sizes_[index] == spec.size() )
			{
				// write at old position
				pos = scan_location_[index];
				fseeko(pFile_,pos,SEEK_SET);		
			}
			else
			{
				// size has changed, forget old position and append
				scan_location_[index] = pos;
			}
				
		}
				
// 		std::cout << "writeScan: writing scan " << index << " at " << ftello(pFile_) << std::endl;
		
        fwrite(&rt,sizeof(CoordinateType),1,pFile_);
		size_t sizeof_peak =  sizeof(PeakType);
		
        for (typename ContainerType::const_iterator cit = spec.getContainer().begin();
                cit != spec.getContainer().end(); ++cit)
        {
            fwrite(&(*cit),sizeof_peak,1,pFile_);
			
// 			CoordinateType it    = cit->getIntensity();
// 			CoordinateType mz = cit->getPosition()[0];
// 			fwrite(&it,sizeof(CoordinateType),1,pFile_);
// 			fwrite(&mz,sizeof(CoordinateType),1,pFile_);
			
        }
        fclose(pFile_);

		// test if this scan was already written and store its size
		if (index >= scan_sizes_.size() )
		{
        	scan_sizes_.push_back(spec.getContainer().size());
		}
		else
		{
			scan_sizes_[index] = spec.getContainer().size();
		}
		
    } // end of write(spectrum)

    /// Reads a spectrum from a file
    void readScan(const size_type& index, SpectrumType& spec)  const
	{
        pFile_ = fopen(file_name_.c_str(),"r");
		
		// set stream to starting point of last writing action
        off_t pos = scan_location_.at(index);
// 		std::cout << " readScan: reading scan " << index << " from " << pos << std::endl;
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
				throw Exception::IndexOverflow(__FILE__, __LINE__,"MSExperimentExtern::readScan()",pos,sizeof(off_t)); 
			}
		}

// 		std::cout << "Reading rt. " << std::endl;
        // read retention time
        CoordinateType rt = 0;
 		fread(&rt,sizeof(CoordinateType),1,pFile_);
        spec.setRetentionTime(rt);

		unsigned int nr_peaks = scan_sizes_.at(index);
// 		std::cout << "Reading peaks: " << nr_peaks << std::endl;
		
		spec.getContainer().clear();
		spec.resize(nr_peaks);
		
		size_t sizeof_peak =  sizeof(PeakType);
				
        //read coordinates of each peak 
		for (typename SpectrumType::Iterator piter = spec.begin(); piter != spec.end(); ++piter)
        {    
			if (fread(&(*piter),sizeof_peak,1,pFile_) == 0)
				std::cout << "Error reading peak data" << std::endl;		

// 			spec[i] = peak;
							
			//++peak_iter;
            //spec.getContainer().push_back(point);
        }
// 		std::cout << "Done."<< std::endl;
        fclose(pFile_);

    }	// end of read const

    /// copies the content of the tempory file
    void copyTmpFile__(String source)
    {
        //std::cout << "Copying temporary file: " << file << std::endl;
        FILE * outFile = fopen(file_name_.c_str(),"w");
        FILE * inFile   = fopen(source.c_str(),"r");

        // check if source file exists
        if (inFile == NULL)
            return;

        // copy file
        while( !feof(inFile) && !ferror(inFile))
            putc( fgetc(inFile), outFile);

        fclose(inFile);
        fclose(outFile);
    }


}
;  // end of class MSExperimentExtern

} // namespace OpenMS

#endif // OPENMS_KERNEL_MSEXPERIMENTEXTERN_H
