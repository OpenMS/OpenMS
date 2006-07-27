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

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include<vector>
#include<algorithm>
#include<limits>

#include <cstdlib>
#include <iostream>

namespace OpenMS
{
/**
	@brief Representation of a mass spectrometry experiment using a temporary file  to store very large data sets.
			
	@ingroup Kernel
**/
template < typename PeakT = DPeak<1> >
class MSExperimentExtern
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
        typedef MSExperimentExtern<IteratorPeakType> IterExperimentType;
        typedef const value_type& reference;
        typedef const value_type* pointer;
        typedef std::random_access_iterator_tag iterator_category;
        typedef unsigned int difference_type;
        //@}

        MSExperimentExternConstIterator()
                : exp_(0), position_(0)
        {}

        MSExperimentExternConstIterator(const IterExperimentType * exp , unsigned int pos)
        {
            exp_       =  (IterExperimentType* ) exp;
            position_ = pos;
        }

        MSExperimentExternConstIterator(IterExperimentType * exp , unsigned int pos)
                : exp_(exp), position_(pos)
        {}

        MSExperimentExternConstIterator(const IterExperimentType& source)
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

        IterExperimentType * exp_;
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
        typedef typename MSExperimentExternConstIterator<IteratorPeakType>::IterExperimentType IterExperimentType;
        typedef typename MSExperimentExternConstIterator<IteratorPeakType>::difference_type difference_type;

        using MSExperimentExternConstIterator<IteratorPeakType>::exp_;
        using MSExperimentExternConstIterator<IteratorPeakType>::position_;


        MSExperimentExternIterator()
                : MSExperimentExternConstIterator<IteratorPeakType>()
        {}

        MSExperimentExternIterator(IterExperimentType * exp, unsigned int position)
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

    }
    ;	// end of class MSExperimentExternIterator


    typedef PeakT PeakType;
    typedef MSSpectrum<PeakT> SpectrumType;
    typedef typename SpectrumType::ContainerType ContainerType;
    typedef typename PeakType::CoordinateType CoordinateType;
    typedef MSExperiment<PeakType> ExperimentType;

    typedef MSExperimentExternIterator<PeakType> Iterator;
    typedef MSExperimentExternConstIterator<PeakType> ConstIterator;
    typedef std::reverse_iterator<Iterator> ReverseIterator;
    typedef std::reverse_iterator<ConstIterator> ConstReverseIterator;

    typedef typename ExperimentType::value_type value_type;
    typedef typename ExperimentType::size_type size_type;
    typedef typename ExperimentType::difference_type difference_type;
    typedef typename ExperimentType::reference reference;
    typedef typename ExperimentType::const_reference const_reference;
    typedef typename ExperimentType::pointer pointer;

    typedef Iterator iterator;
    typedef ConstIterator const_iterator;
    typedef ReverseIterator reverse_iterator;
    typedef ConstReverseIterator const_reverse_iterator;

    /// Standard constructor, allocates a buffer of size 1500
    MSExperimentExtern()
            : buffer_size_(1500),
            scan_location_(), current_scan_(0),
            buffer_index_(0), scan2buffer_(),
            buffer2scan_(), exp_()
    {
        file_name_ = "msexp_" + String(std::rand());
		exp_.resize(buffer_size_);
    }

    /// Copy constructor: copies the content of the temporary file as well (slow !)
    MSExperimentExtern(const MSExperimentExtern& source)
            : buffer_size_(source.buffer_size_),
            scan_location_(source.scan_location_), current_scan_(source.current_scan_),
            buffer_index_(source.buffer_index_),  scan2buffer_(source.scan2buffer_),
            buffer2scan_(source.buffer2scan_), exp_(source.exp_)
    {
        file_name_ = "msexp_" + String(std::rand());
		copyTmpFile_(source.file_name_);
        exp_.resize(buffer_size_);
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
        scan2buffer_    = source.scan2buffer_;
        buffer2scan_    = source.buffer2scan_;
        exp_				  = source.exp_;

        // generate new name for temp file
		std::remove ( file_name_ .c_str());
        file_name_ = "msexp_" + String(std::rand());
		// and copy the old one
        copyTmpFile_(source.file_name_);

        exp_.resize(buffer_size_);

        return *this;
    }

    /// Equality operator
    bool operator== (const MSExperimentExtern& rhs) const
    {
        return (buffer_size_     == rhs.buffer_size_ &&
                scan_location_ == rhs.scan_location_ &&
                buffer_index_    == rhs.buffer_index_ &&
                scan2buffer_    == rhs.scan2buffer_ &&
                buffer2scan_    == rhs.buffer2scan_ &&
                exp_				  == rhs.exp_);
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

    void updateRanges()
    {
        exp_.updateRanges();
    }

    void sortSpectra(bool sort_mz)
    {
        exp_.sortSpectra(sort_mz);
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
        ;
    }

    /// See std::vector documentation.
    ReverseIterator rend()
    {
        return ReverseIterator(begin());
        ;
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
    Iterator erase(Iterator pos)
    {
        delete(&(*pos));
        exp_.erase(exp_.begin()+pos.position_);
        return pos;
    }

    /// See std::vector documentation.
    Iterator erase(Iterator first,Iterator last)
    {
        for (Iterator it=first;it!=last;++it)
        {
            delete(&(*it));
        }
        exp_.erase(exp_.begin()+first.position_,exp_.begin()+last.position_);
        return first;
    }
	
	/// See std::vector documentation.
	reference back()
	{
		return *(exp_.end()-1);
	}
	
	/// See std::vector documentation.
	Iterator insert(Iterator pos, const SpectrumType& spec)
	{
		return pos;
	}

    void push_back(const SpectrumType& spec)
    {
        //std::cout << "Inserting scan " << current_scan_ << std::endl;
        //std::cout << "size: " << buffer_size_ << " index: " << buffer_index_ << std::endl;
        if (buffer_index_ < buffer_size_)
        {
            //std::cout << "Writing in buffer at pos: " << buffer_index_ << std::endl;
            exp_[buffer_index_] = spec;
            scan2buffer_.push_back(buffer_index_++);
            buffer2scan_.push_back(current_scan_++);
        }
        else
        {
            //std::cout << "Inserting dummy values: " << buffer_index_ << std::endl;
            scan2buffer_.push_back(buffer_index_);
            buffer2scan_.push_back(++buffer_index_);
        }
        writeScan(spec);
    }

    /// see std::vector (additionally test if scan is in buffer or needs to be read from temp file)
    reference operator[] (size_type n)
    {

        //std::cout << "operator[" << n << "]" << std::endl;
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
        //std::cout << "operator[" << n << "] const" << std::endl;
        // test if current scan is in buffer
        UnsignedInt b = scan2buffer_[n];
        if (buffer2scan_[b] != n)
            storeInBuffer(n);	// scan is not in buffer, needs to be read from file

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
    }

    /// See std::vector documentation.
    /// Note: the internal vector is of size 1000 by default but this
    /// function returns the number of scans stored so far
    size_type size() const
    {
        return scan2buffer_.size();
    }

    void clear() const
    {
        exp_.clear();
    }

    /// See MSExperiment documentation.
    UnsignedInt getSize() const
    {
        return exp_.getSize();
    }
	
	void resize(UnsignedInt new_size) 
	{
		exp_.resize(new_size);
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

    void deleteTempFile()
    {
        std::remove( file_name_ .c_str());
    }


protected:
    /// size of the internal buffer
    UnsignedInt buffer_size_;

    /// stores the location of the scan in the file
    mutable std::vector<size_type> scan_location_;

    /// number of scans added so far
    mutable UnsignedInt current_scan_;

    /// index in buffer
    mutable UnsignedInt buffer_index_;

    /// The index of each scan in the buffer
    mutable std::vector<size_type> scan2buffer_;

    /// The index of each scan in the buffer
    mutable std::vector<size_type> buffer2scan_;

    /// The index of each scan in the buffer
    mutable std::vector<size_type> scan_sizes_;

    /// Name of the temporary file to store the peak data
    String file_name_;

	/// The internal ms experiment instance.
    mutable ExperimentType exp_;

    /// reads a scan from the temp file and stores it in the buffer
    void storeInBuffer(const size_type& n)
    {
        //std::cout << "spectra at: " << n << " is not in buffer. " << std::endl;
        SpectrumType spec = readScan(n);

        // check if buffer is full
        if (buffer_index_ < buffer_size_)
        {
            //std::cout << "buffer is not full, inserting scan at " << buffer_index_ << std::endl;
            // not full, add scan to buffer
            exp_[buffer_index_] = spec;
            scan2buffer_[n] = buffer_index_;
            buffer2scan_[buffer_index_++] = n;
        }
        else
        {
            // buffer is full, therefore we overwrite the first entry
            buffer_index_ = 0;

            // check if size of buffer is set to zero
            if (buffer_size_ > 0 )
            {
                // if this entry is already occupied by a scan, save it before overwriting
                if (exp_[buffer_index_].getRetentionTime() > 0)
                {
                    writeScan(exp_[buffer_index_]);
                }
                exp_[buffer_index_] = spec;
                scan2buffer_[n] = buffer_index_;
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
        SpectrumType spec = readScan(n);

        // check if buffer is full
        if (buffer_index_ < buffer_size_)
        {
            // not full, add scan to buffer
            exp_[buffer_index_] = spec;
            scan2buffer_[n] = buffer_index_;
            buffer2scan_[buffer_index_++] = n;
        }
        else
        {
            // buffer is full and we overwrite the first entry
            buffer_index_ = 0;

            // check if size of buffer is set to zero
            if (buffer_size_ > 0 )
            {
                exp_[buffer_index_] = spec;
                scan2buffer_[n] = buffer_index_;
                buffer2scan_[buffer_index_++] = n;
            }
            else
            {
                // if we are working on a zero-sized buffer, we need
                // append the current scan
                exp_.push_back(spec);
                scan2buffer_[n] = buffer_index_;
                buffer2scan_[buffer_index_++] = n;
            }
        }
    }

    /// write spectrum to file
    /// format: rt_of_scan first_mz first_intensity second_mz ... last_intensity n
    void writeScan(const SpectrumType& spec)
    {
        FILE * pFile = fopen (file_name_.c_str(),"a");
        float rt = spec.getRetentionTime();
        scan_location_.push_back( ftell(pFile) );


        fwrite(&rt,sizeof(rt),1,pFile);

        for (typename ContainerType::const_iterator cit = spec.getContainer().begin();
                cit != spec.getContainer().end(); ++cit)
        {
            fwrite(&(*cit),sizeof(*cit),1,pFile);
        }
        fclose (pFile);

        scan_sizes_.push_back(spec.getContainer().size());

    } // end of write(spectrum)

    /// Reads a spectrum from a file
    SpectrumType readScan(const size_type& index) const
    {
        SpectrumType spec;

        FILE * pFile = fopen (file_name_.c_str(),"r");

        // set stream to starting point of last writing action
        size_type pos = scan_location_.at(index);
        //std::cout << "Reading from " << pos << std::endl;
        fseek(pFile,pos,SEEK_SET);

        // read retention time
        float rt = 0;
        fread(&rt,sizeof(rt),1,pFile);
        spec.setRetentionTime(rt);

        // read coordinates of each peak and skip
        // spaces in between
        for (unsigned int i=0; i<scan_sizes_.at(index); ++i)
        {
            PeakType point;
            fread(&point,sizeof(point),1,pFile);
            spec.getContainer().push_back(point);
        }
        fclose (pFile);

        return spec;
    }	// end of read const

    /// copies the content of the tempory file
    void copyTmpFile_(String source)
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
