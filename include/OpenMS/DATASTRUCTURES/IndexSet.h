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

#ifndef OPENMS_DATASTRUCTURES_INDEXSET_H
#define OPENMS_DATASTRUCTURES_INDEXSET_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <vector>
#include <iostream>

namespace OpenMS
{
/**
	@brief Sorted and compressed set of indices (UnsignedInt)

	Successive indices are compressed to intervals and thereby a lot of space
	can be saved. 
	
	@ingroup Datastructures
*/
class IndexSet
{

public:
    /**
    	@brief ConstIterator for IndexSet

    	
    */
    class IndexSetConstIterator
    {
        friend class IndexSet;

    public:
        typedef std::vector<UnsignedInt>::size_type VecIndex;

        IndexSetConstIterator(): index_(), pos_(), ref_(0)
        {}

        IndexSetConstIterator(const UnsignedInt index, const VecIndex pos, const IndexSet* ref)
                : index_(index), pos_(pos), ref_(ref)
        {}

        IndexSetConstIterator(const IndexSetConstIterator& it)
                : index_(it.index_), pos_(it.pos_), ref_(it.ref_)
        {}

        ~IndexSetConstIterator()
        {}

        IndexSetConstIterator& operator = (const IndexSetConstIterator& rhs)
        {
            index_=rhs.index_;
            pos_=rhs.pos_;
            ref_=rhs.ref_;
            return *this;
        }

        bool operator < (const IndexSetConstIterator& it) const
        {
            return index_ < it.index_;
        }

        bool operator > (const IndexSetConstIterator& it) const
        {
            return index_ > it.index_;
        }

        bool operator <= (const IndexSetConstIterator& it) const
        {
            return (index_ < it.index_ || index_ == it.index_);
        }

        bool operator >= (const IndexSetConstIterator& it) const
        {
            return (index_ > it.index_ || index_ == it.index_);
        }

        bool operator == (const IndexSetConstIterator& it) const
        {
            return index_ == it.index_ && ref_ == it.ref_;
        }

        bool operator != (const IndexSetConstIterator& it) const
        {
            return !(*this==it);
        }

        const UnsignedInt& operator * () const
        {
            return index_;
        }

        IndexSetConstIterator& operator ++ ()
        {
            // not inside a range or at end or next index equals range end
            if ( ref_->set_[pos_]!=DOTS || pos_+1 >= ref_->set_.size() || index_+1==ref_->set_[pos_+1])
                pos_++;

            if (pos_ >= ref_->set_.size()) //at end
            {
                index_ = END;
                pos_ = END;
            }
            else if(ref_->set_[pos_]==DOTS) // inside a range
                index_++;
            else
                index_ = ref_->set_[pos_];
            return *this;
        }

        IndexSetConstIterator operator ++ (int)
        {
            IndexSetConstIterator tmp(*this);
            ++(*this);
            return tmp;
        }

        IndexSetConstIterator& operator -- ()
        {
            // at end or not inside a range or prev index equals range begin
            if (pos_==END)
                pos_ = ref_->set_.size()-1;
            else if ( ref_->set_[pos_]!=DOTS || index_-1==ref_->set_[pos_-1])
                pos_--;

            if (ref_->set_[pos_]==DOTS) // inside a range
                index_--;
            else
                index_ = ref_->set_[pos_];

            return *this;
        }

        IndexSetConstIterator operator -- (int)
        {
            IndexSetConstIterator tmp(*this);
            --(*this);
            return tmp;
        }

    protected:
        UnsignedInt index_;    // actual index
        VecIndex pos_;  // Position in index vector
        const IndexSet* ref_;
    };

    //========================================================================================

public:
    /// Const iterator
    typedef IndexSetConstIterator ConstIterator;
    /// Const iterator definition for STL compliance
    typedef ConstIterator const_iterator;

    /// standard constructor
    IndexSet();
    /// copy constructor
    IndexSet(const IndexSet& source);
    /// destructor
    virtual ~IndexSet();
    /// constructor of index @p index
    IndexSet(UnsignedInt index);
    /// constructor of interval [ @p index_from , ... , @p index_to ]
    IndexSet(UnsignedInt index_from, UnsignedInt index_to);

    /// assignment operator
    IndexSet& operator = (const IndexSet& source);

    /// tests if the set is empty
    bool isEmpty();

    /// set is empty afterwards
    void clear();

    /// Equality operator
    bool operator == (const IndexSet& rhs) const;
    /// Equality operator
    bool operator != (const IndexSet& rhs) const;

    /// iterator pointing to the first index
    ConstIterator begin() const throw(Exception::Precondition);

    /// iterator pointing to index @p index
    ConstIterator begin(UnsignedInt index) const throw(Exception::Precondition);

    /// iterator pointing behind the last index
    ConstIterator end() const throw(Exception::Precondition);

    /**
    	@brief append index to set (lazy append) 
      
      \return enables concatenation of multiple calls to add
    */
    IndexSet& add
        (UnsignedInt index);


    /**
    	@brief append indices [ @p index_from , ... , @p index_to ] to set (lazy append)

    	\param index_from begin of set
    	\param index_to end of set (omitted if @p to < 0)
    	\return enables concatenation of multiple calls to add
    */
    IndexSet& add
        (UnsignedInt index_from, UnsignedInt index_to);

    /**
    	@brief remove index from set (always starts sorting)
    	
    	\return enables concatenation of multiple calls to remove
    */
    IndexSet& remove
        (UnsignedInt index);


    /**
    	@brief remove indices [ @p index_from , ... , @p index_to ] from set (always starts sorting)

    	\param index_from begin of set
    	\param index_to end of set (omitted if @p to < 0)
    	\return enables concatenation of multiple calls to remove
    */
    IndexSet& remove
        (UnsignedInt index_from, UnsignedInt index_to);

    /**
    	@brief sort and compress indices

    Call before accessing set.
    */
    void sort();

    /// Returns the number of indices in set
    Size size() const;

    friend std::ostream& operator << (std::ostream& os, const IndexSet& set );

protected:
    /**
    @brief sort and compress indices

    indices [ @p skip_from , ... , @p skip_to ] are removed
    */
    void sort_(UnsignedInt skip_from=END, UnsignedInt skip_to=END);

		/// internal value for "..." to define a range of indices
    static const UnsignedInt DOTS;
		/// internal value for end iterator
    static const UnsignedInt END;

		/**
			@brief set of indices

			facilitates compression by containining enum-value DOTS to indicate index ranges
		*/
    std::vector<UnsignedInt> set_;

		/**
			@brief is index set sorted?

			required because of lazy add-function
		*/
    bool is_sorted_;

		/**
			@brief number of indices in set

			number of indices in set independent from compression <br>
			updated by every call to add/remove
		*/
    Size size_;

};

///Print the contents of a IndexSet to a stream.
std::ostream& operator << (std::ostream& os, const IndexSet& set);

}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_INDEXSET_H
