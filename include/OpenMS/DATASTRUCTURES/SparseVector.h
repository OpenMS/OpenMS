// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_DATASTRUCTURES_SPARSEVECTOR_H
#define OPENMS_DATASTRUCTURES_SPARSEVECTOR_H

#include <map>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <sstream>
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>

namespace OpenMS
{
	/**
		@brief SparseVector implementation. The container will not actually store a specified type of element - the sparse element, e.g. zero (by default)

		SparseVector for allround usage, will work with Int, UInt, DoubleReal, Real. This should use less space than a normal vector
		(if more than half of the elements are sparse elements, since the underlying structure is a map) and functions can just
		ignore sparse elements (hop(), @see SparseVectorIterator) for faster look over the elements of the container

		@ingroup Datastructures
	*/
	template <typename Value>
	class SparseVector
	{

		public:

			//forward declarations
			class SparseVectorConstIterator;
			class SparseVectorIterator;
			class SparseVectorReverseIterator;
			class SparseVectorConstReverseIterator;
			class ValueProxy;

			//made available from this classes
			typedef SparseVectorConstIterator const_iterator;
			typedef SparseVectorConstReverseIterator const_reverse_iterator;
			typedef SparseVectorIterator iterator;
			typedef SparseVectorReverseIterator reverse_iterator;

			//remapping
			typedef typename std::map<size_t,Value>::difference_type difference_type; //needed?
			typedef typename std::map<size_t,Value>::size_type size_type;
			typedef typename std::map<size_t,Value>::allocator_type allocator_type; //needed?
			typedef Value value_type;
			typedef Value* pointer; //needed?
			typedef ValueProxy& reference;
			typedef const ValueProxy& const_reference;

			//internal use
			typedef typename std::map<size_t,Value>::const_iterator map_const_iterator;
			typedef typename std::map<size_t,Value>::iterator map_iterator;
			typedef typename std::map<size_t,Value>::const_reverse_iterator reverse_map_const_iterator;
			typedef typename std::map<size_t,Value>::reverse_iterator reverse_map_iterator;

			typedef SparseVectorConstIterator ConstIterator;
			typedef SparseVectorConstReverseIterator ConstReverseIterator;
			typedef SparseVectorIterator Iterator;
			typedef SparseVectorReverseIterator ReverseIterator;


			void print() const
			{
				std::cout << std::endl;
				for (map_const_iterator it=values_.begin(); it!=values_.end(); ++it)
				{
					std::cout << it->first << ": " << it->second << std::endl;
				}
			}

		/// default constructor
		SparseVector():values_(),size_(0),sparse_element_(0)
		{
		}

		/// constructor with chosen sparse element
		SparseVector(Value se):values_(),size_(0),sparse_element_(se)
		{
		}

		 /// detailed constructor, use with filling element value is discouraged unless it is the same as sparse element se
		 SparseVector(size_type size, Value value, Value se=0):values_(),size_(size),sparse_element_(se)
		{
				if(value != sparse_element_) //change, if sparse element is another
				{
					map_iterator i = values_.begin();
					for(size_type s=0; s<size; ++s)
					{
						//makes each insertion in amortized constant time inserted direct after last one
						i = values_.insert(i,make_pair(s,value));
					}
				}
			}

		/// copy constructor
		SparseVector(const SparseVector& source):values_(source.values_),size_(source.size_),sparse_element_(source.sparse_element_)
		{
		}

			/// assignment operator
		SparseVector& operator= (const SparseVector& source)
		{
				if (this != &source)
			{
				values_ = source.values_;
				size_ = source.size_;
				sparse_element_ = source.sparse_element_;
			}
				return *this;
			}

			/// destructor
			~SparseVector()
			{
			}

			/// equality operator
	  	bool operator== (const SparseVector& rhs) const
	  	{
				return ((values_ == rhs.values_)&&(size_ == rhs.size_)&&(sparse_element_ == rhs.sparse_element_));
	  	}

			/// less than operator
			bool operator< (const SparseVector& rhs) const
			{
				return (values_ < rhs.values_);
			}

			/// number of nonzero elements, i.e. the space actually used
			size_type nonzero_size() const
			{
				return values_.size();
			}

			/// size of the represented vector
			size_type size() const
			{
				return size_;
			}

			/// push_back (see stl vector docs)
			void push_back(Value value)
			{
				operator[](size_++) = value;
			}

			/** at (see stl vector docs)

					@param pos index at which the desired element stays
					@throw OutOfRange is thrown if the index is greater or equal than the size of the vector
			*/
			Value at(size_type pos) const
			{
				if (pos >= size_)
				{
					throw Exception::OutOfRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				}
				else
				{
					return operator[](pos);
				}
			}

			/// ValueProxy handles the conversion to int and ,the writing ( if != sparseElement )
			const Value/*Proxy*/ operator[] (size_type pos) const
			{
				assert(pos < size_);
				return (Value)ValueProxy(const_cast<SparseVector&>(*this),pos);
			}

			/// ValueProxy handles the conversion and the writing ( if != sparseElement )
			ValueProxy operator[] (size_type pos)
			{
				assert(pos < size_);
				return ValueProxy(*this,pos);
			}

			/// removes all elements
			void clear()
			{
				values_.clear();
				size_ = 0;
			}

			/// resizes the the vector to param newsize
			void resize(size_type newsize)
			{
				// if the vector is to be smaller
				// delete all invalid entries
				if (newsize < size_)
				{
					for (map_iterator mit = values_.begin(); mit != values_.end();)
					{
						if (mit->first >= newsize)
						{
						  size_type nextvalue = (++mit)->first;
						  values_.erase(--mit);
						  mit = values_.find(nextvalue);
						}
						else
						{
							++mit;
						}
					}
				}
				size_ = newsize;
			}

			/** erase indicated element(iterator) and imediately update indices in map

					@param it parameter which specifies the element which should be deleted
					@throw OutOfRange is thrown if the iterator does not point to an element
			*/
			void erase(SparseVectorIterator it)
			{
				if (it.position() >= size_)
				{
					throw Exception::OutOfRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				}
				//store pointer to the element after the current element
				//erase element
				bool update = false;
				map_iterator mit = values_.find(it.position());
				map_iterator mit_next;
				if ( mit != values_.end() ) //element exists => erase it and update indices of elements after it
				{
					mit_next = mit;
					++mit_next;
					values_.erase(mit);
					update = true;
				}
				else  //element does not exists => update indices of elements after it
				{
					mit_next = values_.lower_bound(it.position());
					update = true;
				}

				//update indices if necessary
				if (update) update_(mit_next,1);

				--size_;
			}

			/** erase indicated element(halfopen iterator-range) and imediately update indices in map
					@param first iterator that points to the begin of the range which should be erased
					@param last iterator that points one position behind the last position which should be erased
			*/
			void erase(SparseVectorIterator first, SparseVectorIterator last)
			{
				if (first.position() >= size_ || last.position() > size_ || last.position() < first.position())
				{
					throw Exception::OutOfRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				}

				size_type amount_deleted = last.position() - first.position();
				map_iterator mfirst = values_.lower_bound(first.position());
				map_iterator mlast = values_.lower_bound(last.position());

				if (mfirst==values_.begin())
				{
					values_.erase(mfirst,mlast);
					update_(values_.begin(),amount_deleted);
				}
				else
				{
					map_iterator start_it = mfirst;
					--start_it;
					values_.erase(mfirst,mlast);
					++start_it;
					update_(start_it,amount_deleted);
				}

				size_ -= amount_deleted;
			}

			///gets an Iterator to the element (including sparseElements) with the minimal value
			SparseVectorIterator getMinElement()
			{
				switch(size_)
				{
					case 0:
						break;
					case 1:
						return begin();
						break;
					default:
						if(values_.size()==0)
						{
							//only sparse elements left
							return begin();
						}
						bool first_sparse_found = false;
						size_type pos = 0;
						map_iterator lowest = values_.begin();
						map_iterator second = values_.begin();
						map_iterator first = second++;
						map_iterator last = values_.end();

						if(lowest->first>0)
						{
							first_sparse_found = true;
						}

						while (second!=last)
						{
							if (second->second<lowest->second) //the first element is covered by initial lowest == frst
							{
								lowest=second;
							}
							if (size_ > values_.size() && !first_sparse_found )
							{
								if((second->first)-(first->first) > 1)
								{
									pos = first->first+1;
									first_sparse_found = true;
								}
							}
							++first; ++second;
						}

						if(size_ == values_.size() || lowest->second < SparseVector::sparse_element_)
						{
							return SparseVectorIterator(*this,lowest->first);
						}
						else //lowest->second >(=) sparseElement
						{
							if(!first_sparse_found)
							{
								return SparseVectorIterator(*this,first->first+1);
							}
							return SparseVectorIterator(*this, pos);
						}
						break;
				}
				return end();
				//map_iterator pos = min_element(values_.begin(), values_.end()); //sorts by map.key :(
			}

			/// begin iterator
			iterator begin()
			{
				return SparseVectorIterator(*this,0);
			}

			/// end iterator
			iterator end()
			{
				return SparseVectorIterator(*this,this->size());
			}

			/// rbegin iterator
			reverse_iterator rbegin()
			{
				return SparseVectorReverseIterator(*this,this->size());
			}

			/// rend iterator
			reverse_iterator rend()
			{
				return SparseVectorReverseIterator(*this,0);
			}

			/// const begin iterator
			const_iterator begin() const
			{
				return SparseVectorConstIterator(*this,0);
			}

			/// const end iterator
			const_iterator end() const
			{
				return SparseVectorConstIterator(*this,this->size());
			}

			/// const begin reverse_iterator
			const_reverse_iterator rbegin() const
			{
				return SparseVectorConstIterator(*this,this->size());
			}

			/// const end reverse_iterator
			const_reverse_iterator rend() const
			{
				return SparseVectorConstIterator(*this,0);
			}

		private:
			/// underlying map
			std::map<size_type, Value> values_;

			/// size including sparse elements
			size_type size_;

		protected:

			/// sparse element
			Value sparse_element_;

			///Updates position of @p it and all larger elements
			void update_(map_iterator it, Size amount_deleted)
			{
				while(it != values_.end())
				{
					size_type tmp_index = it->first;
					Value tmp_value = it->second;
					if (it!=values_.begin())
					{
						//makes insertion in amortized constant time if really inserted directly after mit
						map_iterator tmp_it = it;
						--tmp_it;
						values_.erase(it);
						it = values_.insert(tmp_it,make_pair(tmp_index-amount_deleted,tmp_value));
					}
					else
					{
						//simply insert, as we have no element to insert after
						values_.erase(it);
						it = values_.insert(make_pair(tmp_index-amount_deleted,tmp_value)).first;
					}
					++it;
				}
			}

		public:

		/**
			@brief class ValueProxy allows the SparseVector to differentiate between writing and reading, so zeros can be ignored
			See "more effective c++" section 30
		*/
		class ValueProxy
		{

			public:

				/// public constructor
				ValueProxy(SparseVector& vec,size_type index): vec_(vec), index_(index)
				{
				}

				// if there is a entry in the map from SparseVector, return that
				// if not it is a zero, so return sparseElement
				/// cast operator for implicit casting in case of reading in the vector
				operator double() const
				{
					double value = vec_.sparse_element_;
					map_const_iterator cmit = vec_.values_.find(index_);
					if ( cmit != vec_.values_.end() )
					{
						value = cmit->second;
					}
					return value;
				}

				/// cast operator for implicit casting in case of reading in the vector
				operator int() const
					{
					int value = vec_.sparse_element_;
					map_const_iterator cmit = vec_.values_.find(index_);
					if ( cmit != vec_.values_.end() )
					{
						value = cmit->second;
					}
				return value;
				}

				/// cast operator for implicit casting in case of reading in the vector
				operator float() const
				{
					float value = vec_.sparse_element_;
					map_const_iterator cmit = vec_.values_.find(index_);
					if ( cmit != vec_.values_.end() )
					{
						value = cmit->second;
					}
					return value;
				}

				// maybe more cast-operators for other types

				/// assignment operator, ditches the sparse elements
				ValueProxy& operator= (const ValueProxy& rhs)
				{
					if ((this != &rhs) && (vec_ == rhs.vec_))
					{
						//if rhs' value != sparseElement, cmit!=rhs.vec_.values_.end()
						map_const_iterator cmit = rhs.vec_.values_.find(rhs.index_);
						if (cmit != rhs.vec_.values_.end())
						{
							vec_.values_[rhs.index_] = cmit->second;
						}
						//instead of setting value to zero erase it
						else
						{
							map_iterator mit = vec_.values_.find(rhs.index_);
							if (mit != vec_.values_.end())
					  		{
					   			vec_.values_.erase(mit);
					  		}
						}
						index_ = rhs.index_;
					}
				  	return *this;
				}

				/// assignment operator, ditches the sparse elements
				ValueProxy& operator= (Value val)
				{
					if (val != vec_.sparse_element_) //if (fabs(val) > 1e-8)
					{
						vec_.values_[index_] = val;
					}
					else
					{
						map_iterator mit = vec_.values_.find(index_);
						if (mit != vec_.values_.end())
						{
							vec_.values_.erase(mit);
						}
					}
					return *this;
				}

				/// inequality operator
				bool operator!= (const ValueProxy& other)
				{

					return ((index_ != other.index_) || (&vec_ != &other.vec_));
				}

				/// equality operator
				bool operator== (const ValueProxy& other)
				{
					return !(this != other);
				}

				/// less than operator
				bool operator< (const ValueProxy& other)
				{
					return ((Value)*this < (Value)other);
				}

				/// greater than operator
				bool operator> (const ValueProxy& other)
				{
					return ((Value)*this > (Value)other);
				}

				/// less or equal than operator
				bool operator<= (const ValueProxy& other)
				{
					return ((Value)*this <= (Value)other);
				}

				/// greater or equal than operator
				bool operator>= (const ValueProxy& other)
				{
					return ((Value)*this >= (Value)other);
				}

			private:

				/// the referring SparseVector
				SparseVector& vec_;

				/// the reference into the SparseVector
				size_type index_;

		};  //end of class ValueProxy

		/**
			@brief random access iterator for SparseVector
			including the hop() function to jump to the next non-sparse element
		*/
		class SparseVectorIterator
		{
			friend class SparseVector<Value>;
			friend class SparseVectorConstIterator;

			public:

				/// copy constructor
				SparseVectorIterator(const SparseVectorIterator& source)
				: position_(source.position_),
					vector_(source.vector_),
					valit_(source.valit_)
				{
				}

				/// destructor
				virtual ~SparseVectorIterator()
				{
				}

				/// assignment operator
				SparseVectorIterator& operator= (const SparseVectorIterator& source)
				{
					if (this != &source)
					{
						position_ = source.position_;
						vector_ = source.vector_;
						valit_ = source.valit_;
					}
					return *this;
				}

				///prefix increment
				SparseVectorIterator& operator++()
				{
					++position_;
					return *this;
				}

				///postfix increment
				SparseVectorIterator operator++(int)
				{
					SparseVectorIterator tmp(*this);
					++position_;
					return tmp;
				}

				///prefix decrement
				SparseVectorIterator& operator--()
				{
					--position_;
					return *this;
				}

				///postfix decrement
				SparseVectorIterator operator--(int)
				{
					SparseVectorIterator tmp(*this);
					--position_;
					return tmp;
				}

				/// dereference operator
				ValueProxy operator*()
				{
					assert(position_ < vector_.size_);
					return ValueProxy(this->vector_,position_);
				}

				/// const dereference operator
				const Value operator*() const
				{
					assert(position_ < vector_.size_);
					return (Value)ValueProxy(this->vector_,position_);
				}

				/// indexing
				ValueProxy operator[](size_type n)
				{
					position_ += n;
					assert (position_ < vector_.size_);
					return ValueProxy(this->vector_,position_);
				}

				/// compound assignment +
				SparseVectorIterator& operator+= (const size_type rhs)
				{
					position_ += rhs;
					return *this;
				}

				/// compound assignment -
				SparseVectorIterator& operator-= (const size_type rhs)
				{
					position_ -= rhs;
					return *this;
				}

				/// binary arithmetic +
				SparseVectorIterator operator+ (const size_type rhs) const
				{
					return SparseVectorIterator(vector_, position_+rhs);
				}

				/// binary arithmetic +
				difference_type operator+ (const SparseVectorIterator rhs) const
				{
					return (position_+rhs.position());
				}

				/// binary arithmetic -
				SparseVectorIterator operator- (const size_type rhs) const
				{
					return SparseVectorIterator(vector_, position_-rhs);
				}

				/// binary arithmetic -
				difference_type operator- (const SparseVectorIterator rhs) const
				{
					return (position_-rhs.position());
				}

				/// inequality operator
				bool operator!= (const SparseVectorIterator& other)
				{
					return (position_ != other.position_ || &vector_ != &other.vector_);
				}

				/// equality operator
				bool operator== (const SparseVectorIterator& other)
				{
					return !(*this != other);
				}

				/// less than operator
				bool operator< (const SparseVectorIterator& other)
				{
					return (position_ < other.position());
				}

				/// greater than operator
				bool operator> (const SparseVectorIterator& other)
				{
					return (position_ > other.position());
				}

				/// less or equal than operator
				bool operator<= (const SparseVectorIterator& other)
				{
					return (position_ <= other.position());
				}

				/// greater or equal than operator
				bool operator>= (const SparseVectorIterator& other)
				{
					return (position_ >= other.position());
				}

				/// go to the next nonempty position
				SparseVectorIterator& hop()
				{
					//assert(valit_ != vector_.values_.end() );
					//look for first entry if this is the first call. Go one step otherwise
					if (position_ != valit_->first) //first call
					{
					  valit_ = vector_.values_.upper_bound(position_);
					}
					else
					{
					  ++valit_;
					}
					//check if we are at the end
				  if ( valit_ == vector_.values_.end() )
				  {
				  	position_ = vector_.size_;
				  }
				  else
				  {
				  	position_ = valit_->first;
					}
					return *this;
				}

				/// find out at what position the iterator is; useful in combination with hop()
				size_type position() const
				{
					return position_;
				}

			protected:

				///
				SparseVectorIterator(SparseVector& vector, size_type position)
					: position_(position),
						vector_(vector),
						valit_(vector.values_.begin())
				{
				}

				/// the position in the referred SparseVector
				size_type position_;

				/// the referred SparseVector
				SparseVector& vector_;

				/// the position in the underlying map of SparseVector
				map_const_iterator valit_;

			private:

				/// Not implemented => private
				SparseVectorIterator();

		};//end of class SparseVectorIterator

		/**
			@brief random access reverse iterator for SparseVector
			including the hop() function to jump to the next non-sparse element
		*/
		class SparseVectorReverseIterator
		{
			friend class SparseVector<Value>;
			friend class SparseVectorConstReverseIterator;

			public:

				/// copy constructor
				SparseVectorReverseIterator(const SparseVectorReverseIterator& source)
				: position_(source.position_),
					vector_(source.vector_),
					valrit_(source.valrit_)
				{
				}

				/// destructor
				virtual ~SparseVectorReverseIterator()
				{
				}

				/// assignment operator
				SparseVectorReverseIterator& operator = (const SparseVectorReverseIterator& source)
				{
					if (this != &source)
					{
						position_ = source.position_;
						vector_ = source.vector_;
						valrit_ = source.valrit_;
					}
					return *this;
				}

				/// prefix increment
				SparseVectorReverseIterator& operator++ ()
				{
					--position_;
					return *this;
				}

				/// postfix increment
				SparseVectorReverseIterator operator++ (int)
				{
					SparseVectorReverseIterator tmp(*this);
					--position_;
					return tmp;
				}

				/// prefix decrement
				SparseVectorReverseIterator& operator-- ()
				{
					++position_;
					return *this;
				}

				/// postfix decrement
				SparseVectorReverseIterator operator-- (int)
				{
					SparseVectorReverseIterator tmp(*this);
					++position_;
					return tmp;
				}

				/// dereference operator
				Value operator* ()
				{
					assert(position_ <= vector_.size_);
					assert(position_ != 0);
					return ValueProxy(this->vector_,position_-1);
				}

				/// indexing
				ValueProxy operator[](size_type n)
				{
					position_ -= n;
					assert (position_ < vector_.size_);
					return ValueProxy(this->vector_,position_);
				}

				/// compound assignment +
				SparseVectorReverseIterator& operator+= (const size_type rhs)
				{
					position_ -= rhs;
					return *this;
				}

				/// compound assignment -
				SparseVectorReverseIterator& operator-= (const size_type rhs)
				{
					position_ += rhs;
					return *this;
				}

				/// binary arithmetic +
				SparseVectorReverseIterator operator+ (const size_type rhs) const
				{
					return SparseVectorReverseIterator(vector_, position_-rhs);
				}

				/// binary arithmetic +
				difference_type operator+ (const SparseVectorReverseIterator rhs) const
				{
					return (position_+rhs.position());
				}


				/// binary arithmetic -
				SparseVectorReverseIterator operator- (const size_type rhs) const
				{
					return SparseVectorReverseIterator(vector_, position_+rhs);
				}

				/// binary arithmetic -
				difference_type operator- (const SparseVectorReverseIterator rhs) const
				{
					//what about negatives?
					return -1*(position_-rhs.position());
				}

				/// inequality operator
				bool operator!=(const SparseVectorReverseIterator& other)
				{
					return (position_ != other.position_ || &vector_ != &other.vector_);
				}

				/// equality operator
				bool operator== (const SparseVectorReverseIterator& other)
				{
					return !(*this != other);
				}

				/// less than operator
				bool operator< (const SparseVectorReverseIterator& other)
				{
					return !(*this.position < other.position());
				}

				/// greater than operator
				bool operator> (const SparseVectorReverseIterator& other)
				{
					return !(*this.position > other.position());
				}

				/// less or equal than operator
				bool operator<= (const SparseVectorReverseIterator& other)
				{
					return !(*this.position <= other.position());
				}

				/// greater or equal than operator
				bool operator>= (const SparseVectorReverseIterator& other)
				{
					return !(*this.position >= other.position());
				}

				/// go to the next nonempty position
				SparseVectorReverseIterator& rhop()
				{
					assert(valrit_ != reverse_map_const_iterator( vector_.values_.rend() ) );
					//look for first entry if this is the first call. Go one step otherwise
					if (position_-1 != valrit_->first )
					{
						valrit_ = reverse_map_const_iterator(--(vector_.values_.find(position_-1)));
					}
					else
					{
						++valrit_;
					}
					//check if we are at the end(begin)
					if ( valrit_ == reverse_map_const_iterator( vector_.values_.rend() ) )
					{
						position_ = 0;
					}
					else
					{
						position_ = valrit_->first+1;
					}
					return *this;
				}

				/// find out at what position the iterator is; useful in combination with hop()
				size_type position() const
				{
					return position_;
				}

			public:

				/// detailed constructor
				SparseVectorReverseIterator(SparseVector& vector, size_type position)
					: position_(position),
						vector_(vector),
						valrit_(vector.values_.rbegin())
				{
				}

			protected:

				/// the position in the referred SparseVector
				size_type position_;

			private:
		 		/// reffered sparseVector
				SparseVector& vector_;

				/// the position in the underlying map of SparseVector
				reverse_map_const_iterator valrit_;

				/// Not implemented => private
				SparseVectorReverseIterator();



		};//end of class SparseVectorReverseIterator

		/// const_iterator for SparseVector
		class SparseVectorConstIterator
		{
				friend class SparseVector<Value>;
				friend class SparseVectorIterator;

			public:

				/// copy constructor
				SparseVectorConstIterator(const SparseVectorConstIterator& source)
				: position_(source.position_),
					vector_(source.vector_),
					valit_(source.valit_)
				{
				}

				/// copy constructor from SparseVector::SparseVectorIterator
				SparseVectorConstIterator(const SparseVectorIterator& source)
				: position_(source.position_),
					vector_(source.vector_),
					valit_(source.valit_)
				{
				}

				/// destructor
				virtual ~SparseVectorConstIterator()
				{
				}

				/// assignment operator
		  	SparseVectorConstIterator& operator= (const SparseVectorConstIterator& source)
		  	{
					if (this != &source)
					{
						position_ = source.position_;
						const_cast<SparseVector&>(this->vector_) = source.vector_;
						valit_ = source.valit_;
					}
					return *this;
				}

				/// postincrement operator
				SparseVectorConstIterator& operator++ ()
				{
					assert(position_ <= vector_.size_);
					++position_;
					return *this;
				}

				/// immidiate increment operator
				SparseVectorConstIterator operator++ (int)
				{
					SparseVectorConstIterator tmp(*this);
					++position_;
					assert(position_ <= vector_.size_);
					return tmp;
				}

				/// postincrement operator
				SparseVectorConstIterator& operator-- ()
				{
					assert(position_ <= vector_.size_);
					--position_;
					return *this;
				}

				/// immidiate increment operator
				SparseVectorConstIterator operator-- (int)
				{
					SparseVectorConstIterator tmp(*this);
					--position_;
					assert(position_ <= vector_.size_);
					return tmp;
				}


				/// derefence operator
				const Value operator* () const
				{
					assert(position_ < vector_.size_);
					return (Value)ValueProxy(const_cast<SparseVector&>(this->vector_),position_);
				}

				// indexing
				const ValueProxy operator[](size_type n) const
				{
					position_ += n;
					assert (position_ < vector_.size_);
					return ValueProxy(const_cast<SparseVector&>(this->vector_),position_);
				}

				/// compound assignment +
				SparseVectorConstIterator& operator+= (const size_type rhs)
				{
					position_ += rhs;
					return *this;
				}

				/// compound assignment -
				SparseVectorConstIterator& operator-= (const size_type rhs)
				{
					position_ -= rhs;
					return *this;
				}

				/// binary arithmetic +
				SparseVectorConstIterator operator+ (const size_type rhs) const
				{
					return SparseVectorConstIterator(const_cast<SparseVector&>(this->vector_), position_+rhs);
				}

				/// binary arithmetic -
				SparseVectorConstIterator operator- (const size_type rhs) const
				{
					return SparseVectorConstIterator(const_cast<SparseVector&>(this->vector_), position_-rhs);
				}

				/// inequality operator
				bool operator!= (const SparseVectorConstIterator& other)
				{
					return (position_ != other.position_ || &vector_ != &other.vector_);
				}

				/// equality operator
				bool operator== (const SparseVectorConstIterator& other)
				{
					return !(*this != other);
				}

				/// less than operator
				bool operator< (const SparseVectorConstIterator& other)
				{
					return (*this.position < other.position());
				}

				/// greater than operator
				bool operator> (const SparseVectorConstIterator& other)
				{
					return (*this.position > other.position());
				}

				/// less or equal than operator
				bool operator<= (const SparseVectorConstIterator& other)
				{
					return (*this.position <= other.position());
				}

				/// greater or equal than operator
				bool operator>= (const SparseVectorConstIterator& other)
				{
					return (*this.position >= other.position());
				}

				/// go to the next nonempty position
				SparseVectorConstIterator& hop()
				{
					assert(valit_ != vector_.values_.end() );
					//look for first entry if this is the first call. Go one step otherwise
					if (position_ != valit_->first) //first call
					{
					  valit_ = vector_.values_.upper_bound(position_);
					}
					else
					{
						++valit_;
					}
					//check if we are at the end
					if ( valit_ == vector_.values_.end() )
					{
						position_ = vector_.size_;
					}
					else
					{
						position_ = valit_->first;
					}
					return *this;
				}

				/// find out at what position the iterator is, useful in combination with hop()
				size_type position() const
				{
					return position_;
				}

			protected:
				/// default constructor
				SparseVectorConstIterator();

				/// detailed constructor
				SparseVectorConstIterator(const SparseVector& vector, size_type position)
					: position_(position),
						vector_(vector),
						valit_(vector.values_.begin())
				{
				}

			private:
				/// position in reffered SparseVector
				mutable size_type position_;

				/// referring to this SparseVector
				const SparseVector& vector_;

				/// the position in the underlying map of SparseVector
				map_const_iterator valit_;

		};	//end of class	SparseVectorConstIterator

		/// const_reverse_iterator for SparseVector
		class SparseVectorConstReverseIterator
		{
			friend class SparseVector<Value>;

			public:

				/// copy constructor
				SparseVectorConstReverseIterator(const SparseVectorConstIterator& source)
					: position_(source.position_),
						vector_(source.vector_),
						valrit_(source.valrit_)
				{
				}

				/// copy constructor from SparseVector::SparseVectorIterator
				SparseVectorConstReverseIterator(const SparseVectorReverseIterator& source)
					: position_(source.position_),
						vector_(source.vector_),
						valrit_(source.valrit_)
				{
				}

				/// destructor
				virtual ~SparseVectorConstReverseIterator()
				{
				}

				/// assignment operator
				SparseVectorConstReverseIterator& operator = (const SparseVectorConstReverseIterator& source)
				{
					if (this != &source)
					{
						position_ = source.position_;
						const_cast<SparseVector&>(this->vector_) = source.vector_;
						valrit_ = source.valrit_;
					}
					return *this;
				}

				/// postincrement operator
				SparseVectorConstReverseIterator& operator ++ ()
				{
					//assert(position_ < 0);
					--position_;
					return *this;
				}

				/// immidiate increment operator
				SparseVectorConstReverseIterator operator++ (int)
				{
					SparseVectorConstIterator tmp(*this);
					--position_;
					//assert(position_ < 0);
					return tmp;
				}
				/// postdecrement operator
				SparseVectorConstReverseIterator& operator-- ()
				{
					//assert(position_ < 0);
					++position_;
					return *this;
				}

				/// immidiate decrement operator
				SparseVectorConstReverseIterator operator-- (int)
				{
					SparseVectorConstIterator tmp(*this);
					++position_;
					//assert(position_ < 0);
					return tmp;
				}

				/// derefence operator
				ValueProxy operator* ()
				{
					assert(position_ <= vector_.size_);
					assert(position_ != 0);
					return ValueProxy(const_cast<SparseVector&>(this->vector_),position_-1);
				}

				/// go to the next nonempty position
				SparseVectorConstReverseIterator& rhop()
				{
					assert(valrit_ != vector_.values_.rend() );
					//look for first entry if this is the first call. Go one step otherwise
					if (position_-1 != valrit_->first )
					{
						valrit_ = reverse_map_const_iterator(--(vector_.values_.find(position_-1)));
					}
					else
					{
						++valrit_;
					}
					//check if we are at the end(begin)
					if ( valrit_ == reverse_map_const_iterator( vector_.values_.rend() ) )
					{
						position_ = 0;
					}
					else
					{
						position_ = valrit_->first+1;
					}
					return *this;
				}

				/// find out at what position the iterator is, useful in combination with hop()
				size_type position() const
				{
					return position_;
				}

				/// inequality operator
				bool operator!= (const SparseVectorConstReverseIterator& other)
				{
					return (position_ != other.position_ || &vector_ != &other.vector_);
				}

	 		protected:

				/// default constructor
				SparseVectorConstReverseIterator();

				/// detailed constructor
					SparseVectorConstReverseIterator(const SparseVector& vector, size_type position)
					:position_(position),vector_(vector),valrit_(vector.values_.rbegin())
				{
				}

			private:

				// the position in SparseVector
				mutable size_type position_;

				/// referenc to the vector operating on
				const SparseVector& vector_;

				// the position in the underlying map of SparseVector
				reverse_map_const_iterator valrit_;

		};	//end of class	SparseVectorConstReverseIterator


	};//end of class SparseVector

}
#endif //OPENMS_DATASTRUCTURES_SPARSEVECTOR_H

