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
// $Maintainer: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_HASHMAP_H
#define OPENMS_DATASTRUCTURES_HASHMAP_H

#include <OpenMS/CONCEPT/Exception.h>
#include<OpenMS/CONCEPT/Types.h>
#include<OpenMS/CONCEPT/HashFunction.h>

#include <utility>
#include <algorithm>
#include <iterator>
#include <vector>

namespace OpenMS
{
	/**	
		@brief Generic Hash Map Class.
		
		This class implements a simple hash map using hashing by chaining.	
		
		@ingroup Datastructures
	*/
	template <class Key, class T>
	class HashMap
	{
		public:

		// Convenience typedefs
		///
		typedef std::pair<Key, T> ValueType;
		///
		typedef Key KeyType;
		///
		typedef std::pair<Key, T>* PointerType;
			

		/**
			@brief HashMap node (internal use only)
			
			
		*/
		struct Node
		{
			Node*			next;
			ValueType	value;

			Node(const ValueType& my_value, const Node* my_next)
				throw()
				: next(const_cast<Node*>(my_next)),
					value(const_cast<ValueType&>(my_value))
			{
			}
		};
		
		class ConstIterator;

		/**
			@brief HashMap iterator class
			
			
			
		*/
		class Iterator
		{
			
			friend class HashMap;
			friend class ConstIterator;
			
			public:
			typedef ValueType value_type;
			typedef Index difference_type;
			typedef std::forward_iterator_tag iterator_category;
			typedef value_type& reference;
			typedef value_type* pointer;
			
			Iterator() throw() {}

			Iterator(const Iterator& it) throw()
				: position_(it.position_),
					bucket_(it.bucket_),
					bound_(it.bound_)
			{
			}

			~Iterator() throw() {}

			Iterator& operator = (const Iterator& it) throw()
			{
				position_ = it.position_;
				bucket_ = it.bucket_;
				bound_ = it.bound_;
				return *this;
			}

			bool operator == (const Iterator& it) const throw()
			{
				return (position_ == it.position_);
			}

			bool operator != (const Iterator& it) const throw()
			{
				return (position_ != it.position_);
			}
			
			reference operator * () throw()
			{
				return **position_;
			}

			const reference operator * () const throw()
			{
				return *position_;
			}

			pointer operator -> () throw()
			{
				return &position_->value;
			}

			const pointer operator -> () const
			{
				return &position_->value;
			}

			Iterator& operator ++ ()
			{
        position_ = position_->next;

        if (position_ == 0)
        {
					for (++bucket_;  bucket_ < (Position)bound_->bucket_.size();  ++bucket_)
					{
						position_ = bound_->bucket_[bucket_];

						if (position_ != 0)
						{
							return *this;
						}
					}
				}

				return *this;
			} 

			Iterator operator ++ (int)
			{
	      Iterator tmp(*this);
	      ++(*this);//????
	      return tmp;
			} 

			static Iterator end(const HashMap& hm)
			{
				Iterator it;
				it.position_ = 0;
				it.bound_ = const_cast<HashMap*>(&hm);
				it.bucket_ = 0;
				return it;
			}

			static Iterator begin(const HashMap& hm)
			{
				Iterator it;
				it.bound_ = const_cast<HashMap*>(&hm);
        for (it.bucket_ = 0;  it.bucket_ < hm.bucket_.size();  ++it.bucket_)
        {
          it.position_ = hm.bucket_[it.bucket_];

          if (it.position_ != 0)
          {
            return it;
					}
				}
				return it;
			}

			protected:
			Node* position_;
			Position bucket_;
			HashMap* bound_;
		};	

		/**
			@brief HashMap const_iterator class
			
			
			
		*/
		class ConstIterator
		{			
			friend class HashMap;
			
			public:
			typedef ValueType value_type;
			typedef Index difference_type;
			typedef const ValueType& reference;
			typedef const ValueType* pointer;
			typedef std::forward_iterator_tag iterator_category;
			
			ConstIterator() throw() 
				:	position_(0),
					bucket_(0),
					bound_(0)
			{
			}

			ConstIterator(const ConstIterator& it) throw()
				:	position_(it.position_),
					bucket_(it.bucket_),
					bound_(it.bound_)
			{
			}

			ConstIterator(const Iterator& it)
				:	position_(it.position_),
					bucket_(it.bucket_),
					bound_(it.bound_)
			{
			}

			~ConstIterator() throw() {}

			ConstIterator& operator = (const ConstIterator& it)
			{
				position_ = it.position_;
				bucket_ = it.bucket_;
				bound_ = it.bound_;
				return *this;
			}

			ConstIterator& operator = (const Iterator& it)
			{
				position_ = it.position_;
				bucket_ = it.bucket_;
				bound_ = it.bound_;
				return *this;
			}

			bool operator == (const Iterator& it) const
			{
				return position_ == it.position_;
			}

			bool operator != (const Iterator& it) const
			{
				return position_ != it.position_;
			}
			
			bool operator == (const ConstIterator& it) const
			{
				return position_ == it.position_;
			}

			bool operator != (const ConstIterator& it) const
			{
				return position_ != it.position_;
			}
			
			ConstIterator& operator ++ ()
			{
        position_ = position_->next;

        if (position_ == 0)
        {
					for (++bucket_;  bucket_ < (Position)bound_->bucket_.size();  ++bucket_)
					{
						position_ = bound_->bucket_[bucket_];

						if (position_ != 0)
						{
							return *this;
						}
					}
				}

				return *this;
			} 

			ConstIterator operator ++ (int)
			{
	      ConstIterator tmp(*this);
	      ++(*this);
	      return tmp;
			} 

			ConstIterator& operator -- ()
			{
				--position_;
				return *this;
			} 

			ConstIterator operator -- (int)
			{
	      ConstIterator tmp(*this);
	      --(*this);
	      return tmp;
			} 

			reference operator * () const
			{
				return position_->value;
			}

			pointer operator -> () const
			{
				return &position_->value;
			}

			static ConstIterator end(const HashMap& hm)
			{
				ConstIterator it;
				it.position_ = 0;
				it.bucket_ = 0;
				it.bound_ = &hm;
				return it;
			}

			static ConstIterator begin(const HashMap& hm)
			{
				ConstIterator it;
				it.bound_ = &hm;
        for (it.bucket_ = 0;  it.bucket_ < hm.bucket_.size();  ++it.bucket_)
        {
          it.position_ = hm.bucket_[it.bucket_];

          if (it.position_ != 0)
          {
            return it;
					}
				}
				return it;
			}

			protected:
			Node* position_;
			Position bucket_;
			const HashMap* bound_;
		};


		typedef std::reverse_iterator<Iterator> ReverseIterator;
		typedef std::reverse_iterator<ConstIterator> ConstReverseIterator;

		/** @name STL-compliance type definitions
		*/
		//@{
		///
		typedef ValueType value_type;
		typedef Size size_type;
		typedef Index difference_type;
		///
		typedef ValueType& reference;
		typedef const ValueType& const_reference;
		///
		typedef ValueType& pointer;
		///
		typedef Iterator iterator;
		typedef ConstIterator const_iterator;
		typedef ReverseIterator reverse_iterator;
		typedef ConstReverseIterator const_reverse_iterator;
		//@}

		/**	@name	 Enums and Constants
		*/
		//@{
		enum
		{
			/// Initial capacity of the empty hash map
			INITIAL_CAPACITY          = 100,

			/// Initial number of buckets of the empty hash map
			INITIAL_NUMBER_OF_BUCKETS = 50
		};
		//@}

		/**	@name	Exceptions
		*/
		//@{
			
		/**	
			@brief HashMap illegal key exception
			
			@ingroup Exceptions
		*/
		class IllegalKey
			:	public Exception::Base
		{
			public:
			IllegalKey(const char* file, int line, const char* function)
				:	Exception::Base(file, line, function) {}
		};

		//@}

		/**	@name Constructors and Destructors 
		*/
		//@{

		/**	Default constructor.
				Create a new and empty hash map.
				@param initial_capacity the capacity of the hash map
				@param number_of_buckets the number of buckets to create
		*/
		HashMap(Size initial_capacity = INITIAL_CAPACITY, Size number_of_buckets = INITIAL_NUMBER_OF_BUCKETS)
			throw();
			
		/**	Copy Constructor.
		*/
		HashMap(const HashMap& hash_map) throw();

		/**	Destructor.
		*/
		inline
		virtual ~HashMap() throw()
		{
			destroy();
			deleteBuckets_();
		}

		/**	Clear the hash map.
				Remove all nodes from all buckets.
				The capacity and the number of buckets remain unchanged.
		*/
		virtual void clear() throw();
	
		/**	Clear the hash map.
				Remove all nodes from all buckets.
				The capacity and the number of buckets remain unchanged.
				Simply calls clear.
		*/
		void destroy() throw();

		//@}
		/**	@name Assignment 
		*/
		//@{

		/**	Assignment from another hash map.
				@param hash_map the hash map to assign from
		*/
		void set(const HashMap& hash_map) throw();

		/**	Assignment operator.
				Assign the contents of a hash map to another.
				@param hash_map the hash map to assign from
		*/
		const HashMap& operator = (const HashMap& hash_map) throw();

		/**	Assign the contents of this hash map to another map.
		*/
		void get(HashMap& hash_map) const throw();

		/**	Swap the contents of two hash maps.
		*/
		void swap(HashMap& hash_map) throw();

		//@}
		/**	@name	Accessors
		*/
		//@{

		/**	Return the number of buckets
		*/
		Size getBucketSize() const throw();

		/** Return the capcacity of the hash map.
		*/
		Size getCapacity() const throw();

		/**	Return the number of entries in the map.
		*/
		Size getSize() const throw();
			
		/**	Return the number of entries in the map.
		*/
		Size size() const throw();

		/** Find the element whose key is <tt>key</tt>.
		*/
		Iterator find(const Key& key) throw()
		{
			Iterator it = end();
			HashIndex bucket = hash_(key);
			Node*	node_ptr = bucket_[hash_(key)];
			
			for (; node_ptr != 0; node_ptr = node_ptr->next)
			{
				if (node_ptr->value.first == key)
				{
					it.position_ = node_ptr;
					it.bucket_ = bucket;
					break;
				}
			} 

			return it;
		}
	
		/** Find the element whose key is <tt>key</tt>.
		*/
		ConstIterator find(const Key& key) const throw()
		{
			ConstIterator it = end();
			HashIndex bucket = hash_(key);
			Node*	node_ptr = bucket_[hash_(key)];
			
			for (; node_ptr != 0; node_ptr = node_ptr->next)
			{
				if (node_ptr->value.first == key)
				{
					it.position_ = node_ptr;
					it.bucket_ = bucket;
					break;
				}
			} 

			return it;
		}

		/**	Return a mutable reference to the element whose key is <tt>key</tt>.
				If an element with the key <tt>key</tt> does not exist, it is inserted.
				@param	key the key
		*/
		T& operator [] (const Key& key) throw();

		/**	Return a constant reference to the element whose key is <tt>key</tt>.
				@exception IllegalKey if the given key does not exist
				@param	key the key
		*/
		const T& operator [] (const Key& key) const 
			throw(typename HashMap<Key, T>::IllegalKey);

		/**	Insert a new entry into the hash map.
		*/
		std::pair<Iterator, bool> insert(const ValueType& entry) throw();

		/**	Insert a new entry into the hash map.
				For STL compatibility. The value of <tt>pos</tt> is ignored.
		*/
		Iterator insert(Iterator pos, const ValueType& entry) throw();

		/**	Erase element with key <tt>key</tt>.
				@return Size the number of elements erased (0 or 1)
		*/
		Size erase(const Key& key) throw();

		/**	Erase element at a given position.
				@param pos an iterator pointing to the element to delete
		*/
		void erase(Iterator pos) throw();

		/**	Erase a range of elements.
				Erase all elements in the range <tt>first - last</tt>.
		*/
		void erase(Iterator first, Iterator last) throw();

		//@}


		/**	@name	Predicates
		*/
		//@{

		/**	Test whether the map contains the given key.
		*/
		bool has(const Key& key) const throw();

		/**	Test whether the map is empty.
		*/
		bool isEmpty() const throw();

		/**	Compare two hash maps.
		*/
		bool operator == (const HashMap& hash_map) const throw();

		/**	Compare two hash maps.
		*/
		bool operator != (const HashMap& hash_map) const throw();

		//@}

		inline Iterator begin()	throw()
		{
			return Iterator::begin(*this);
		}

		inline Iterator end() throw()
		{
			return Iterator::end(*this);
		}

		inline ConstIterator begin() const throw()
		{
			return ConstIterator::begin(*this);
		}

		inline ConstIterator end() const throw()
		{
			return ConstIterator::end(*this);
		}


		protected:

		virtual Node* newNode_(const ValueType& value, Node* next) const throw();

		virtual void deleteNode_(Node* node) const throw();
	
		virtual HashIndex hash(const Key& key) const throw();

		virtual bool needRehashing_() const throw();

		virtual void rehash() throw();

		PointerType find_(const Key& key, HashIndex& index) throw();
			
		PointerType find_(const Key& key, HashIndex& index) const throw();

		void deleteBuckets_() throw();

		HashIndex hash_(const Key& key) const throw();

		void rehash_() throw();

		/**	@name Attributes
		*/
		//@{

		/**	The number of entries in the map
		*/
		Size size_;
		/**	The maximum number of entries before a resize operation is required
		*/
		Size capacity_;

		/**	Buckets are stored as a vector of linked lists of Nodes 
		*/
		std::vector<Node*> bucket_;

		//@}
	};

	template <class Key, class T>
	inline
	HashMap<Key, T>::HashMap(Size initial_capacity, Size number_of_buckets)
		throw()
		:	size_(0),
			capacity_(initial_capacity),
			bucket_(number_of_buckets)
	{
		for (Position bucket = 0; bucket < (Position)bucket_.size(); ++bucket)
		{
			bucket_[bucket] = 0;
		}
	}

	template <class Key, class T>
	inline
	HashMap<Key, T>::HashMap(const HashMap& hash_map)
		throw()
		:	size_(hash_map.size_),
			capacity_(hash_map.capacity_),
			bucket_((Size)hash_map.bucket_.size())
	{
		Node* node = 0;
		
		for (Position bucket = 0; bucket < (Position)bucket_.size(); ++bucket)
		{
			bucket_[bucket] = 0;

			for (node = hash_map.bucket_[bucket]; node != 0; node = node->next)
			{
				bucket_[bucket] = newNode_(node->value, bucket_[bucket]);
			}
		}
	}

	template <class Key, class T>
	void HashMap<Key, T>::clear()
		throw()
	{
		Node* node = 0;
		Node* next_node = 0;
		
		for (Position bucket = 0; bucket < (Position)bucket_.size(); ++bucket)
		{
			for (node = bucket_[bucket]; node != 0; node = next_node)
			{
				next_node = node->next;
				deleteNode_(node);
			}
			
			bucket_[bucket] = 0;
		}

		size_ = 0;
	}

	template <class Key, class T>
	inline 
	void HashMap<Key, T>::destroy()
		throw()
	{
		clear();
	}

	template <class Key, class T>
	void HashMap<Key, T>::set(const HashMap& hash_map)
		throw()
	{
		if (&hash_map == this)
		{
			return;
		}

		destroy();
		deleteBuckets_();

		size_ = hash_map.size_;
		capacity_ = hash_map.capacity_;
		bucket_.resize(hash_map.bucket_.size());

		Node* node = 0;
		
		for (Position bucket = 0; bucket < (Position)bucket_.size(); ++bucket)
		{
			bucket_[bucket] = 0;

			for (node = hash_map.bucket_[bucket]; node != 0; node = node->next)
			{
				bucket_[bucket] = newNode_(node->value, bucket_[bucket]);
			}
		}
	}

	template <class Key, class T>
	inline 
	const HashMap<Key, T>& HashMap<Key, T>::operator = (const HashMap& hash_map)
		throw()
	{
		set(hash_map);
		return *this;
	}

	template <class Key, class T>
	inline 
	void HashMap<Key, T>::get(HashMap& hash_map) const
		throw()
	{
		hash_map.set(*this);
	}

	template <class Key, class T>
	inline 
	void HashMap<Key, T>::swap(HashMap& hash_map)
		throw()
	{
		std::swap(size_, hash_map.size_);
		std::swap(capacity_, hash_map.capacity_);
		std::swap(bucket_, hash_map.bucket_);
	}

	template <class Key, class T>
	inline 
	Size HashMap<Key, T>::getBucketSize() const
		throw()
	{
		return (Size)bucket_.size();
	}

	template <class Key, class T>
	inline 
	Size HashMap<Key, T>::getCapacity() const
		throw()
	{
		return capacity_;
	}

	template <class Key, class T>
	inline 
	Size HashMap<Key, T>::getSize() const
		throw()
	{
		return size_;
	}

	template <class Key, class T>
	inline 
	Size HashMap<Key, T>::size() const
		throw()
	{
		return size_;
	}

	template <class Key, class T>
	inline 
	T& HashMap<Key, T>::operator [] (const Key& key)
		throw()
	{
		Iterator it = find(key);
		if (it == end())
		{
			T value;
			std::pair<Iterator, bool> result = insert(ValueType(key, value));
			it = result.first;
		} 
		
		return it->second;
	}

	template <class Key, class T>
	inline 
	const T& HashMap<Key, T>::operator [] (const Key& key) const
		throw(typename HashMap<Key, T>::IllegalKey)
	{
		ConstIterator it = find(key);
		if (it == end())
		{
			throw IllegalKey(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		} 
		else 
		{
			return it->second;
		}
	}

	template <class Key, class T>
	::std::pair<typename HashMap<Key, T>::Iterator, bool> HashMap<Key, T>::insert
		(const ValueType& item)	throw()
	{
		Iterator it = find(item.first);
		if (it == end())
		{
			if (needRehashing_() == true)
			{
				rehash_();
			}
			
			HashIndex bucket = hash_(item.first);
			
			Node* node_ptr = bucket_[bucket];
			bucket_[bucket] = newNode_(item, node_ptr);
			
			++size_;
			it.position_	= bucket_[bucket];
			it.bucket_		= bucket;

			return std::pair<Iterator, bool>(it, true);
		} 
		else 
		{
			// replace the existing value
			it->second = item.second;

			return ::std::pair<Iterator, bool>(it, false);
		}
	}

	template <class Key, class T>
	inline
	typename HashMap<Key, T>::Iterator HashMap<Key, T>::insert
		(typename HashMap<Key, T>::Iterator /* pos */, const ValueType& entry)	throw()
	{
		return insert(entry).first;
	}

	template <class Key, class T>
	Size HashMap<Key, T>::erase(const Key& key)
		throw()
	{
		Node*	previous = 0;
		HashIndex bucket = hash_(key);
		Node*	node_ptr = bucket_[bucket];

		while (node_ptr != 0 && node_ptr->value.first != key)
		{
			previous = node_ptr;
			node_ptr = node_ptr->next;
		}

		if (node_ptr == 0)
		{
			return false;
		}

		if (node_ptr == bucket_[bucket])
		{
			bucket_[bucket] = node_ptr->next;
		} 
		else 
		{
			previous->next = node_ptr->next;
		}

		deleteNode_(node_ptr);
		--size_;

		return true;
	}
		
	template <class Key, class T>
	void HashMap<Key, T>::erase(Iterator pos)
		throw()
	{
		if ((pos == end()) || (size_ == 0))
		{
			return;
		}
				
		if (pos.position_ == bucket_[pos.bucket_])
		{
			bucket_[pos.bucket_] = pos.position_->next;
		} 
		else 
		{
			// walk over all nodes in this bucket and identify the predecessor
			// of the node refered to by the iterator pos
			Node* prev = bucket_[pos.bucket_];
			for (; (prev != 0) && (prev->next != pos.position_); prev = prev->next);
			if (prev != 0)
			{
				// remove the node and reconnect the list
				prev->next = pos.position_->next;
			}
			else 
			{
				throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}
		}

		// delete the node and decrement the set size
		deleteNode_(pos.position_);
		--size_;
	}

	template <class Key, class T>
	void HashMap<Key, T>::erase(Iterator f, Iterator l)
		throw()
	{
		if (f == end())
		{
			return;
		}

		Position last_bucket = l.bucket_;
		if (l == end())
		{
			last_bucket = (Position)bucket_.size() - 1;
		}

		if (f.bucket_ > last_bucket)
		{
			// empty range - l < f
			return;
		}

		// count the deleted entries to correct the set size
		Size no_deletions = 0;

		Position bucket = f.bucket_;
		for (; bucket <= last_bucket; bucket++)
		{
			if (bucket_[bucket] == 0)
			{
				// skip all empty buckets
				continue;
			}

			if ((bucket == f.bucket_) && (bucket_[bucket] != f.position_))
			{
				// find the predecessor of f
				Node* n = bucket_[bucket];
				Node* next;
				for (; (n->next != f.position_) && (n->next != 0); n = n->next);
				
				if (bucket == last_bucket)
				{
					// delete everything from f to l in this bucket

					next = n->next;
					n->next = l.position_;
					for (n = next; (n != 0) && (n != l.position_); n = next)
					{
						next = n->next;
						deleteNode_(n);
						no_deletions++;
					}
				}
				else
				{
					// delete everything from f to the end in this bucket

					if (n != 0)
					{
						// mark the end of the list
						next = n->next;
						n->next = 0;

						// delete all remaining nodes
						for (n = next; n != 0; n = next)
						{
							next = n->next;
							deleteNode_(n);
							no_deletions++;
						}
					}
				}
			} 
			// if the current bucket lies between the first and the last bucket...
			else if (bucket < last_bucket)
			{
				// ...delete the whole bucket
				Node* next;
				for (Node* n = bucket_[bucket]; n != 0; n = next)
				{
					next = n->next;
					deleteNode_(n);
					no_deletions++;
				}
				bucket_[bucket] = 0;
			}
			else if (bucket == last_bucket)
			{
				// we delete everything in this bucket up to the iterator l

				// find the predecessor of l
				Node* n = bucket_[bucket];
				Node* next;
				for (; (n != 0) && (n != l.position_); n = next)
				{
					next = n->next;
					deleteNode_(n);
					no_deletions++;
				}

				bucket_[bucket] = l.position_;
			}
		}

		// correct the set size
		size_ -= no_deletions;
	}

	template <class Key, class T>
	inline 
	bool HashMap<Key, T>::has(const Key& key) const
		throw()
	{
		return (find(key) != end());
	}

	template <class Key, class T>
	inline 
	bool HashMap<Key, T>::isEmpty() const
		throw()
	{
		return (size_ == 0);
	}

	template <class Key, class T>
	bool HashMap<Key, T>::operator == (const HashMap& hash_map) const
		throw()
	{
		if (size_ != hash_map.size_) 
		{
			return false;
		}
		
		for (ConstIterator it(begin()); it != end(); ++it)
		{
			ConstIterator hash_map_it(hash_map.find(it->first));
			if ((hash_map_it == hash_map.end()) || (hash_map_it->second != it->second))
			{
				return false;
			}
		}
		
		return true;
	}

	template <class Key, class T>
	inline
	bool HashMap<Key, T>::operator != (const HashMap& hash_map) const
		throw()
	{
		return !(*this == hash_map);
	}

	template <class Key, class T>
	inline HashIndex HashMap<Key, T>::hash(const Key& key) const
		throw()
	{
		return Hash(key);
	}

	template <class Key, class T>
	inline void HashMap<Key, T>::rehash()
		throw()
	{
		capacity_ = (Size)getNextPrime((Size)bucket_.size() * 2);
	}


	template <class Key, class T>
	void HashMap<Key, T>::deleteBuckets_()
		throw()
	{
		Node*	node = 0;
		Node*	next_node = 0;
		for (Position i = 0; i < (Position)bucket_.size(); i++)
		{
			node = bucket_[i];
			while (node != 0)
			{
				next_node = node->next;
				deleteNode_(node);
				node = next_node;
			}
			bucket_[i] = 0;
		}
	}

	template <class Key, class T>
	inline typename HashMap<Key, T>::Node* HashMap<Key, T>::newNode_
		(const ValueType& value, typename HashMap<Key, T>::Node* next) const
		throw()
	{
		return new Node(value, next);
	}

	template <class Key, class T>
	inline void HashMap<Key, T>::deleteNode_(typename HashMap<Key, T>::Node* node) const
		throw()
	{
		delete node;
	}

	template <class Key, class T>
	inline bool HashMap<Key, T>::needRehashing_() const
		throw()
	{
		return (size_ >= capacity_);
	}


	template <class Key, class T>
	inline HashIndex HashMap<Key, T>::hash_(const Key& key) const
		throw()
	{
		return (HashIndex)(hash(key) % bucket_.size());
	}
 
	template <class Key, class T>
	void HashMap<Key, T>::rehash_()
		throw()
	{
		// calculate the new number of buckets (in capacity_)
		rehash();
		
		// save the old contents
		std::vector<Node*> old_buckets(bucket_);

		// resize the bucket vector and initialize it with zero
		bucket_.clear();
		bucket_.resize(capacity_);
		Position i;
		for (i = 0; i < capacity_; i++)
		{
			bucket_[i] = 0;
		}

		// rehash the old contents into the new buckets
		Node*	node;
		Node* next_node;
		for (Position i = 0; i < (Position)old_buckets.size(); ++i)
		{
			for (node = old_buckets[i]; node != 0; node = next_node)
			{
				next_node = node->next;
				Position new_bucket = (Position)hash_(node->value.first);
				node->next = bucket_[new_bucket];
				bucket_[new_bucket] = node; 
			}
		}		
	}
} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_HASHMAP_H
