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
// $Id: DPeakArray.sip,v 1.1.1.1 2005/04/04 13:03:04 marc_sturm Exp $
// $Author: marc_sturm $
// --------------------------------------------------------------------------
//

//	template <Size D, typename Peak = DPeak<D> >
class DPeakArray1
	: PersistenceManagerInterface
{
%TypeHeaderCode
	#include <OpenMS/KERNEL/DPeakArray.h>
	using namespace OpenMS;
	typedef DPeakArray<1> DPeakArray1;
%End

	public:
	
	/*
	//ConstIterator for the DPeakArray
	template <class IteratorPeakType>
	class DPeakArrayConstIterator
	{
		
		friend class DPeakArray;
		
		public:
		typedef IteratorPeakType value_type;
		typedef typename std::vector<IteratorPeakType*>::difference_type difference_type;
		typedef const value_type& reference;
		typedef const value_type* pointer;
		typedef std::random_access_iterator_tag iterator_category;
		
		DPeakArrayConstIterator()
		{
		}
		
		DPeakArrayConstIterator(const typename std::vector<IteratorPeakType*>* vec , unsigned int position)
		{
			vector_ = (typename std::vector<IteratorPeakType*>*)vec;
			position_ = position;
		}

		DPeakArrayConstIterator(typename std::vector<IteratorPeakType*>* vec , unsigned int position)
		{
			vector_ = vec;
			position_ = position;
		}

		DPeakArrayConstIterator(const DPeakArrayConstIterator& it)
		{
			vector_=it.vector_;
			position_=it.position_;
		}

		~DPeakArrayConstIterator()
		{
			
		}

		DPeakArrayConstIterator& operator = (const DPeakArrayConstIterator& rhs)
		{
			vector_=rhs.vector_;
			position_=rhs.position_;
			return *this;
		}
	
		bool operator < (const DPeakArrayConstIterator& it) const
		{
			return position_ < it.position_;
		}

		bool operator > (const DPeakArrayConstIterator& it) const
		{
			return position_ > it.position_;
		}

		bool operator <= (const DPeakArrayConstIterator& it) const
		{
			return (position_ < it.position_ || position_ == it.position_);
		}

		bool operator >= (const DPeakArrayConstIterator& it) const
		{
			return (position_ > it.position_ || position_ == it.position_);
		}
		
		bool operator == (const DPeakArrayConstIterator& it) const
		{
			return position_ == it.position_ && vector_ == it.vector_;
		}

		bool operator != (const DPeakArrayConstIterator& it) const
		{
			return position_ != it.position_ || vector_ != it.vector_;
		}
		
		DPeakArrayConstIterator& operator ++ ()
		{
			position_ += 1;
			return *this;
		} 

		DPeakArrayConstIterator operator ++ (int)
		{
			DPeakArrayConstIterator tmp(*this);
			++(*this);
			return tmp;
		} 

		DPeakArrayConstIterator& operator -- ()
		{
			position_-=1;
			return *this;
		} 

		DPeakArrayConstIterator operator -- (int)
		{
			DPeakArrayConstIterator tmp(*this);
			--(*this);
			return tmp;
		} 

		DPeakArrayConstIterator operator - (difference_type n) const
		{
			DPeakArrayConstIterator tmp(*this);
			tmp.position_ -= n;
			return tmp;
		} 

		DPeakArrayConstIterator operator + (difference_type n) const
		{
			DPeakArrayConstIterator tmp(*this);
			tmp.position_ += n;
			return tmp;
		} 

		DPeakArrayConstIterator& operator += (difference_type n)
		{
			position_ += n;
			return *this;
		} 

		DPeakArrayConstIterator& operator -= (difference_type n)
		{
			position_ -= n;
			return *this;
		}

		friend difference_type operator - ( const DPeakArrayConstIterator& i1, const DPeakArrayConstIterator& i2 )
		{
			return (i1.position_ - i2.position_);
		} 

		friend DPeakArrayConstIterator operator + ( difference_type n, const DPeakArrayConstIterator& i )
		{
			DPeakArrayConstIterator tmp(i);
			tmp.position_ += n;
			return tmp;
		} 
		
		reference operator * ()
		{
			return *((*vector_)[position_]);
		}
		
		pointer operator -> ()
		{
			return (*vector_)[position_];
		}

		pointer operator -> () const
		{
			return (*vector_)[position_];
		}
		
		reference operator [] (difference_type n)
		{
			return *((*this)+n);
		}
		
		protected:
		typename std::vector<IteratorPeakType*>* vector_;
		unsigned int position_; 
	};
	

	//Iterator for the DPeakArray
	template <class IteratorPeakType>
	class DPeakArrayIterator : public DPeakArrayConstIterator<IteratorPeakType>
	{
		
		friend class DPeakArray;
		
		public:
		
		typedef typename DPeakArrayConstIterator<IteratorPeakType>::value_type& reference;
		typedef typename DPeakArrayConstIterator<IteratorPeakType>::value_type* pointer;
		
		DPeakArrayIterator()
		{
			
		}

		DPeakArrayIterator(typename std::vector<IteratorPeakType*>* vec, unsigned int position): DPeakArrayConstIterator<IteratorPeakType>(vec,position)
		{

		}

		DPeakArrayIterator(const DPeakArrayIterator<IteratorPeakType>& it): DPeakArrayConstIterator<IteratorPeakType>(it)
		{

		}

		~DPeakArrayIterator()
		{
			
		}
		
		reference operator * ()
		{
			return *((*vector_)[DPeakArrayConstIterator<IteratorPeakType>::position_]);
		}

		pointer operator -> ()
		{
			return (*vector_)[DPeakArrayConstIterator<IteratorPeakType>::position_];
		}

		const pointer operator -> () const
		{
			return (*vector_)[DPeakArrayConstIterator<IteratorPeakType>::position_];
		}

		typename DPeakArrayIterator::reference operator [] (typename DPeakArrayIterator::difference_type n)
		{
			return *((*this)+n);
		}

		DPeakArrayIterator& operator ++ ()
		{
			DPeakArrayConstIterator<IteratorPeakType>::position_+=1;
			return *this;
		} 

		DPeakArrayIterator operator ++ (int)
		{
			DPeakArrayIterator tmp(*this);
			++(*this);
			return tmp;
		} 

		DPeakArrayIterator& operator -- ()
		{
			DPeakArrayConstIterator<IteratorPeakType>::position_-=1;
			return *this;
		} 

		DPeakArrayIterator operator -- (int)
		{
			DPeakArrayIterator tmp(*this);
			--(*this);
			return tmp;
		} 

		DPeakArrayIterator operator - (typename DPeakArrayIterator::difference_type n) const
		{
			DPeakArrayIterator tmp(*this);
			tmp.position_ -= n;
			return tmp;
		} 

		DPeakArrayIterator operator + (typename DPeakArrayIterator::difference_type n) const
		{
			DPeakArrayIterator tmp(*this);
			tmp.position_ += n;
			return tmp;
		} 

		friend DPeakArrayIterator operator + (typename DPeakArrayIterator::difference_type n, const DPeakArrayIterator& i )
		{
			DPeakArrayIterator tmp(i);
			tmp.position_ += n;
			return tmp;
		} 

		DPeakArrayIterator& operator += (typename DPeakArrayIterator::difference_type n)
		{
			DPeakArrayConstIterator<IteratorPeakType>::position_ += n;
			return *this;
		} 

		DPeakArrayIterator& operator -= (typename DPeakArrayIterator::difference_type n)
		{
			DPeakArrayConstIterator<IteratorPeakType>::position_ -= n;
			return *this;
		}
		
		friend void swap(DPeakArrayIterator& i1, DPeakArrayIterator& i2)
		{
			unsigned int tmp = i1.position_;
			i1.position_ = i2.position_;
			i2.position_ = tmp;
		}
		
		protected:
		
	};	
	*/
	
	typedef DPeak1 PeakType;
	typedef DPeak1& reference;
	typedef const DPeak1& const_reference;
	/*
	typedef DPeakArrayIterator<PeakType> Iterator;
	typedef DPeakArrayConstIterator<PeakType> ConstIterator;
	typedef std::reverse_iterator<Iterator> ReverseIterator;
	typedef std::reverse_iterator<ConstIterator> ConstReverseIterator;
	typedef typename std::vector<PeakType>::value_type value_type;
	typedef typename std::vector<PeakType>::size_type size_type;
	typedef typename std::vector<PeakType>::difference_type difference_type;
	///
	typedef typename std::vector<PeakType>::reference reference;
	typedef typename std::vector<PeakType>::const_reference const_reference;
	///
	typedef typename std::vector<PeakType>::pointer pointer;
	///
	typedef Iterator iterator;
	typedef ConstIterator const_iterator;
	typedef ReverseIterator reverse_iterator;
	typedef ConstReverseIterator const_reverse_iterator;
	*/

	void push_back(const PeakType& /x/);
	void pop_back();
	Size size() const;
	Size capacity() const;
	void reserve(Size /n/);
	Size max_size() const;

	/*	
	Iterator begin()
	{
		return Iterator((std::vector<PeakType*>*)&vector_,(unsigned int)0);
	}

	Iterator end()
	{
		return Iterator((std::vector<PeakType*>*)&vector_,(unsigned int)(vector_.size()));
	}

	ConstIterator begin() const
	{
		return ConstIterator((const std::vector<PeakType*>*)&vector_,(unsigned int)0);
	}

	ConstIterator end() const
	{
		return ConstIterator((const std::vector<PeakType*>*)&vector_,(unsigned int)(vector_.size()));
	}

	ReverseIterator rbegin()
	{
		return ReverseIterator(end());
	}

	ReverseIterator rend()
	{
		return ReverseIterator(begin());
	}

	ConstReverseIterator rbegin() const
	{
		return ConstReverseIterator(end());
	}

	ConstReverseIterator rend() const
	{
		return ConstReverseIterator(begin());
	}
	*/

	void resize(Size /new_size/);
	reference front();
	const_reference front() const;
	reference back();
	const_reference back() const;
	void clear();
	bool empty() const;

	// ???? reference operator [] (Size /n/);
	// ???? const_reference operator [](Size /n/) const;
	bool operator == (const DPeakArray1& /array/) const;
	bool operator !=(const DPeakArray1& /array/) const;
	bool operator < (const DPeakArray1& /array/) const;
	bool operator > (const DPeakArray1& /array/) const;
	bool operator <= (const DPeakArray1& /array/) const;
	bool operator >= (const DPeakArray1& /array/) const;
	void swap(DPeakArray1& /array/);

	//	friend void swap(DPeakArray& a1, DPeakArray& a2);
	//	Iterator insert(Iterator pos, const PeakType& peak);
	// 	void insert(Iterator pos, size_type n, const PeakType& peak);

	//	template <class InputIterator>
	//	void insert(Iterator pos, InputIterator f, InputIterator l)
	//	Iterator erase(Iterator pos)
	//	Iterator erase(Iterator first,Iterator last)

	DPeakArray1();
	DPeakArray1(Size /n/);
	DPeakArray1(const DPeakArray1& /p/);
	~DPeakArray1(); 

	//	template <class InputIterator>
	// DPeakArray(InputIterator f, InputIterator l) : PersistenceManagerInterface(), capacity_(0) 

	//	DPeakArray& operator = (const DPeakArray& rhs) 

	// template <class InputIterator>
	// void assign(InputIterator f , InputIterator l)

	// ???? void assign(Size /n/);
	void sortWrtIntensity();
	void sortWrtPosition(Index /i/ = -1);
	void sortWrtWidth(Index /i/ = -1);

	//	template < typename Cmp >
	// void sortWrtComparator ()

	//	void serialize(PersistenceManager& f)
	//	void deserialize_(const std::string& attribute_name, const std::string& value)
};
