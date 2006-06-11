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
// $Id: DPeakList.h,v 1.13 2006/06/01 15:26:18 groepl Exp $
// $Author: groepl $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_DPEAKLIST_H
#define OPENMS_KERNEL_DPEAKLIST_H

#include <OpenMS/config.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>

#include <list>
#include <algorithm>
#include <iterator>

namespace OpenMS
{

	/**	
		@brief Peak container implemented as a list.
		
		This class represents an list of D-dimensional peaks.
		It is fulfills the requirements of a STL list class, but provides more a more 
		convenient interface to manipulate these vectors, sort with
		respect to specific dimensions or intensity and 
		a convenient interface to the other OpenMS classes.

		<i><b>ATTENTION:</b> DPeakList can be used with DPeak and instances of
		classes derived from DPeak at a time, but do not use mutating STL
		algorithms on an inhomogenous DPeakList!  These algorithms may affect
		only the DPeak part of the contained Peaks. </i>
	
		@ingroup Kernel, Serialization
	*/
	
	template <Size D, typename PeakT = DPeak<D> >
	class DPeakList
		: public PersistentObject
	{
		public:

		template <class IteratorPeakType> class DPeakListConstIterator;

		/**
			@brief Iterator of the DPeakList.
			
			The template argument is the PeakType of the DPeakList.
		*/
		template <class IteratorPeakType>
		class DPeakListIterator
		{
			
			friend class DPeakList;
			friend class DPeakListConstIterator<IteratorPeakType>;
			
			public:
			typedef IteratorPeakType value_type;
			typedef typename std::list<IteratorPeakType*>::difference_type difference_type;
			typedef std::bidirectional_iterator_tag iterator_category;
			typedef typename DPeakListIterator<IteratorPeakType>::value_type& reference;
			typedef typename DPeakListIterator<IteratorPeakType>::value_type* pointer;
			
			DPeakListIterator()
			{
			}

			DPeakListIterator(typename std::list<IteratorPeakType*>::iterator it)
			{
				iterator_ = it;
			}

			DPeakListIterator(const DPeakListIterator& it)
			{
				iterator_ = it.iterator_;
			}

			~DPeakListIterator()
			{
			}

			DPeakListIterator& operator = (const DPeakListIterator& it)
			{
				if (this==&it) return *this;
				
				iterator_ = it.iterator_;
				
				return *this;
			}

			bool operator == (const DPeakListIterator& it) const
			{
				return (iterator_ == it.iterator_);
			}

			bool operator != (const DPeakListIterator& it) const
			{
				return (iterator_ != it.iterator_);
			}
			
			reference operator * ()
			{
				return **iterator_;
			}

			const reference operator * () const
			{
				return **iterator_;
			}

			pointer operator -> ()
			{
				return &**iterator_;
			}

			const pointer operator -> () const
			{
				return &**iterator_;
			}

			DPeakListIterator& operator ++ ()
			{
				++iterator_;
				return *this;
			} 

			DPeakListIterator operator ++ (int)
			{
	      DPeakListIterator tmp(*this);
	      ++(*this);
	      return tmp;
			} 

			DPeakListIterator& operator -- ()
			{
				--iterator_;
				return *this;
			} 

			DPeakListIterator operator -- (int)
			{
	      DPeakListIterator tmp(*this);
	      --(*this);
	      return tmp;
			} 

			friend void swap(DPeakListIterator& i1, DPeakListIterator& i2)
			{
				typename std::list<IteratorPeakType*>::iterator tmp =  i1.iterator_;
				i1.iterator_ = i2.iterator_;
				i2.iterator_ = tmp;
			}
			
			protected:
			typename std::list<IteratorPeakType*>::iterator iterator_;
		};	
		
		/**
			@brief Const iterator of the DPeakList.
			
			The template argument is the PeakType of the DPeakList.
		*/
		template <class IteratorPeakType>
		class DPeakListConstIterator
		{
			
			friend class DPeakList;
			
			public:
			typedef IteratorPeakType value_type;
			typedef typename std::list<IteratorPeakType*>::difference_type difference_type;
			typedef const value_type& reference;
			typedef const value_type* pointer;
			typedef std::bidirectional_iterator_tag iterator_category;
			
			DPeakListConstIterator()
			{
			}
			
			DPeakListConstIterator(typename std::list<IteratorPeakType*>::const_iterator it)
			{
				iterator_ = it;
			}

			DPeakListConstIterator(typename std::list<IteratorPeakType*>::iterator it)
			{
				iterator_ = it;
			}

			DPeakListConstIterator(const DPeakListConstIterator& it)
			{
				iterator_ = it.iterator_;
			}

			DPeakListConstIterator(const DPeakListIterator<IteratorPeakType>& it)
			{
				iterator_ = it.iterator_;
			}

			~DPeakListConstIterator()
			{
			}

			DPeakListConstIterator& operator = (const DPeakListConstIterator& it)
			{
				if (this==&it) return *this;
				
				iterator_ = it.iterator_;
				
				return *this;
			}

			DPeakListConstIterator& operator = (const DPeakListIterator<IteratorPeakType>& it)
			{
				iterator_ = it.iterator_;
				return *this;
			}

			bool operator == (const DPeakListIterator<IteratorPeakType>& it) const
			{
				return iterator_ == it.iterator_;
			}

			bool operator != (const DPeakListIterator<IteratorPeakType>& it) const
			{
				return iterator_ != it.iterator_;
			}
			
			bool operator == (const DPeakListConstIterator& it) const
			{
				return iterator_ == it.iterator_;
			}

			bool operator != (const DPeakListConstIterator& it) const
			{
				return iterator_ != it.iterator_;
			}
			
			DPeakListConstIterator& operator ++ ()
			{
				++iterator_;
				return *this;
			} 

			DPeakListConstIterator operator ++ (int)
			{
	      DPeakListConstIterator tmp(*this);
	      ++(*this);
	      return tmp;
			} 

			DPeakListConstIterator& operator -- ()
			{
				--iterator_;
				return *this;
			} 

			DPeakListConstIterator operator -- (int)
			{
	      DPeakListConstIterator tmp(*this);
	      --(*this);
	      return tmp;
			} 

			reference operator * () const
			{
				return **iterator_;
			}

			pointer operator -> () const
			{
				return &**iterator_;
			}

			friend void swap(DPeakListConstIterator& i1, DPeakListConstIterator& i2)
			{
				typename std::list<IteratorPeakType*>::iterator tmp =  i1.iterator_;
				i1.iterator_ = i2.iterator_;
				i2.iterator_ = tmp;
			}
			
			protected:
			typename std::list<IteratorPeakType*>::const_iterator iterator_;
		};

		

		/** @name Type definitions
		*/
		//@{
		typedef PeakT PeakType;
		typedef DPeakListIterator<PeakType> Iterator;
		typedef DPeakListConstIterator<PeakType> ConstIterator;
		typedef std::reverse_iterator<Iterator> ReverseIterator;
		typedef std::reverse_iterator<ConstIterator> ConstReverseIterator;
		
		// STL compliance definitions 
		typedef typename std::list<PeakType>::value_type value_type;
		typedef typename std::list<PeakType>::size_type size_type;
		typedef typename std::list<PeakType>::difference_type difference_type;
		typedef typename std::list<PeakType>::reference reference;
		typedef typename std::list<PeakType>::const_reference const_reference;
		typedef typename std::list<PeakType>::pointer pointer;
		typedef typename std::list<PeakType*>::allocator_type allocator_type;
		typedef Iterator iterator;
		typedef ConstIterator const_iterator;
		typedef ReverseIterator reverse_iterator;
		typedef ConstReverseIterator const_reverse_iterator;
		//@}

		/**	@name Constructors and Destructor
		*/
		//@{
		/// See std::list documentation.
		DPeakList(): PersistentObject() 
		{
		}
		/// See std::list documentation.
		DPeakList(size_type n) : PersistentObject()
		{
			list_=std::list<PeakType*>(n);
			for (typename std::list<PeakType*>::iterator it=list_.begin();it!=list_.end();++it)
			{
				PeakType* peak = new PeakType();
				*it = peak;
			}
		} 
		/// See std::list documentation.
		DPeakList(size_type n, const PeakType& peak) : PersistentObject()
		{
			list_=std::list<PeakType*>(n);
			for (typename std::list<PeakType*>::iterator it=list_.begin();it!=list_.end();++it)
			{
				PeakType* peak_p = reinterpret_cast<PeakType*>(peak.clone());
				*it = peak_p;
			}
		} 
		/// See std::list documentation.
		DPeakList(const DPeakList& p) : PersistentObject(p)
		{
			PeakType* peak;
			for (ConstIterator it=p.begin(); it!=p.end();++it)
			{
				peak = reinterpret_cast<PeakType*>(it->clone());
				list_.push_back(peak);
			}
		}
		/// See std::list documentation.
		template <class InputIterator>
		DPeakList(InputIterator f, InputIterator l) : PersistentObject()
		{
			PeakType* pointer;
			for (InputIterator it=f;it!=l;++it)
			{
				pointer = reinterpret_cast<PeakType*>(it->clone());
				list_.push_back(pointer);
			}
		}
		/// See std::list documentation.
		~DPeakList() 
		{
			for (typename std::list<PeakType*>::iterator it=list_.begin();it!=list_.end();++it)
			{
				delete(*it);
			}			
		}
		//@}
	
		/// See std::list documentation.
		DPeakList& operator = (const DPeakList& rhs)
		{
			if (this==&rhs) return *this;
			
			clear();
			PeakType* peak;
			for (ConstIterator it=rhs.begin(); it!=rhs.end();++it)
			{
				peak = reinterpret_cast<PeakType*>(it->clone());
				list_.push_back(peak);
			}
			
			return *this;
		}
		
		/// See std::list documentation.
		size_type size() const
		{
			return list_.size();
		}
		
		/// See std::list documentation.
		void resize(size_type new_size, const PeakType& t=PeakType())
		{
			size_type old_size = list_.size();
			if (new_size<old_size)
			{
				typename std::list<PeakType*>::iterator it = list_.end();
				for (size_type i=old_size ; i>new_size;i--)
				{
					--it;
					delete(*it);
				}					
				list_.resize(new_size);		
			}
			else if (new_size>old_size)
			{
				typename std::list<PeakType*>::iterator it = list_.end();
				--it;
				list_.resize(new_size);
				++it;
				for (typename std::list<PeakType*>::iterator i=it ; i!=list_.end();i++)
				{
					PeakType* peak = reinterpret_cast<PeakType*>(t.clone());
					*i = peak;
				}
			}
		}
		
		/// See std::list documentation.
		size_type max_size() const
		{
			return list_.max_size();
		}
		
		/// See std::list documentation.
		bool empty() const
		{
			return (list_.size()==0);
		}
		
		/// See std::list documentation.
		void clear()
		{
			for (typename std::list<PeakType*>::iterator it=list_.begin();it!=list_.end();++it)
			{
				delete(*it);
			}
			list_.clear();
		}
		
		/// See std::list documentation.
		Iterator begin()
		{
			return Iterator(list_.begin());
		}
		
		/// See std::list documentation.
		Iterator end()
		{
			return Iterator(list_.end());
		}
		
		/// See std::list documentation.
		ConstIterator begin() const
		{
			return ConstIterator(list_.begin());
		}
		
		/// See std::list documentation.
		ConstIterator end() const
		{
			return ConstIterator(list_.end());
		}
		
		/// See std::list documentation.
		ReverseIterator rbegin()
		{
			return ReverseIterator(end());
		}
		
		/// See std::list documentation.
		ReverseIterator rend()
		{
			return ReverseIterator(begin());
		}
		
		/// See std::list documentation.
		ConstReverseIterator rbegin() const
		{
			return ConstReverseIterator(end());
		}
		
		/// See std::list documentation.
		ConstReverseIterator rend() const
		{
			return ConstReverseIterator(begin());
		}
		
		/// See std::list documentation.
		reference front()
		{
			return *(begin());
		}
		
		/// See std::list documentation.
		const_reference front() const
		{
			return *(begin());
		}
		
		/// See std::list documentation.
		reference back()
		{
			Iterator it=end();
			--it;
			return *(it);
		}
		
		/// See std::list documentation.
		const_reference back() const
		{
			Iterator it=end();
			--it;
			return *(it);
		}
		
		/// See std::list documentation.
		void push_back(const PeakType& x)
		{
			PeakType* peak = reinterpret_cast<PeakType*>(x.clone());
			list_.push_back(peak);
		}
		
		/// See std::list documentation.
		void push_front(const PeakType& x)
		{
			PeakType* peak = reinterpret_cast<PeakType*>(x.clone());
			list_.push_front(peak);
		}
		
		/// See std::list documentation.
		void pop_back()
		{
			typename std::list<PeakType*>::iterator it = list_.end();
			--it;
			delete *(it);
			list_.pop_back();
		}
		
		/// See std::list documentation.
		void pop_front()
		{
			delete *(list_.begin());
			list_.pop_front();
		}
		
		/// See std::list documentation.
		void swap(DPeakList& list)
		{
			list_.swap(list.list_);
		}
		
		/// See std::list documentation.
		friend void swap(DPeakList& l1, DPeakList& l2)
		{
			l1.list_.swap(l2.list_);
		}
		
		/// See std::list documentation.
		Iterator insert(Iterator pos, const PeakType& peak)
		{
			PeakType* pointer = reinterpret_cast<PeakType*>(peak.clone());
			list_.insert(pos.iterator_,pointer);
			return pos;
		}
		
		/// See std::list documentation.
		void insert(Iterator pos, size_type n, const PeakType& peak)
		{
			//std::cout << "INSERT(pos,n,value)"<<std::endl;
			std::vector<PeakType*> tmp(n);
			for (size_type i=0;i<n;i++)
			{
				tmp[i] = reinterpret_cast<PeakType*>(peak.clone());
			}
			//std::cout << "tmp.size(): "<< tmp.size()<<std::endl;
			list_.insert(pos.iterator_,tmp.begin(),tmp.end());
		}
		
		/// See std::list documentation.
		template <class InputIterator>
		void insert(Iterator pos, InputIterator f, InputIterator l)
		{
			std::vector<PeakType*> tmp;
			for (InputIterator it=f;it!=l;++it)
			{
				tmp.push_back( reinterpret_cast<PeakType*>(it->clone()));
			}
			list_.insert(pos.iterator_,tmp.begin(),tmp.end());
		}
		
		/// See std::list documentation.
		Iterator erase(Iterator pos)
		{
			Iterator tmp=pos++;
			delete(&(*tmp));
			list_.erase(tmp.iterator_);
			return pos;
		}
		
		/// See std::list documentation.
		Iterator erase(Iterator first,Iterator last)
		{
			Iterator tmp=last++;
			for (Iterator it=first;it!=tmp;++it)
			{
				delete(&(*it));
			}
			list_.erase(first.iterator_,tmp.iterator_);
			return last;
		}
		
		/// See std::list documentation.
		bool operator == (const DPeakList& rhs) const
		{
			if (size()!=rhs.size())
			{
				return false;
			}
			typename std::list<PeakType*>::const_iterator it=list_.begin();
			typename std::list<PeakType*>::const_iterator rhs_it=rhs.list_.begin();
			for (unsigned int i=0;i<size();i++)
			{
				if (typeid(**it)!=typeid(**rhs_it))
				{
					return false;
				}
				if (**it!=**rhs_it)
				{
					return false;
				}
			}
			return true;
		}

		/// See std::list documentation.
		bool operator !=(const DPeakList& list) const
		{
			return !(operator==(list));
		}
		
		/// Comparison of container sizes
		bool operator < (const DPeakList& list) const
		{
			return size() < list.size();
		}
		
		/// Comparison of container sizes
		bool operator > (const DPeakList& list) const
		{
			return size() > list.size();
		}
		
		/// Comparison of container sizes
		bool operator <= (const DPeakList& list) const
		{
			return operator<(list) || operator==(list);
		}
		
		/// Comparison of container sizes
		bool operator >= (const DPeakList& list) const
		{
			return operator>(list) || operator==(list);
		}

		/// See std::list documentation.		
		void splice (iterator position, DPeakList& x)
		{
			list_.splice(position.iterator_,x.list_);
		}
		
		/// See std::list documentation.
		void splice (iterator position, DPeakList& x,iterator i)
		{
			list_.splice(position.iterator_,x.list_,i.iterator_);
		}

		/// See std::list documentation.
		void splice (iterator position, DPeakList& x,iterator f,iterator l)
		{
			list_.splice(position.iterator_,x.list_,f.iterator_,l.iterator_);
		}

		/// See std::list documentation.
		void remove(const PeakType& p)
		{
			typename std::list<PeakType*>::iterator it=list_.begin();
			while (it!=list_.end())
			{
				if (**it==p && typeid(p)==typeid(**it))
				{
					delete (*it);
					list_.erase (it++);
				}
				else
				{
					++it;
				}
			}
		}

		/// See std::list documentation.
		template <class Predicate>
		void remove_if(Predicate p)
		{
			typename std::list<PeakType*>::iterator it=list_.begin();
			while (it!=list_.end())
			{
				if (p(**it))
				{
					delete (*it);
					list_.erase (it++);
				}
				else
				{
					++it;
				}
			}
		}

		/// See std::list documentation.
		void unique()
		{
			typename std::list<PeakType*>::iterator it=list_.begin();
			typename std::list<PeakType*>::iterator it2=list_.begin();
			it2++;
			while (it2!=list_.end())
			{
				if (**it==**it2)
				{
					delete(*it2);
					list_.erase(it2);
					it2=it;
					++it2;
				}
				else
				{
				++it;
				++it2;
				}
			}
		}

		/// See std::list documentation.
		template <class BinaryPredicate>
		void unique(BinaryPredicate p)
		{
			typename std::list<PeakType*>::iterator it=list_.begin();
			typename std::list<PeakType*>::iterator it2=list_.begin();
			it2++;
			while (it2!=list_.end())
			{
				if (p(**it,**it2))
				{
					delete(*it2);
					list_.erase(it2);
					it2=it;
					++it2;
				}
				else
				{
				++it;
				++it2;
				}
			}
		}

		/// See std::list documentation.
		void merge(DPeakList& list)
		{
			list_.merge(list.list_,PointerComparator<typename PeakType::PositionLess>());
		}

		/// See std::list documentation.
		template <class BinaryPredicate>
		void merge(DPeakList& list , BinaryPredicate p)
		{
			list_.merge(list.list_,DPeakListBinaryPredicate<BinaryPredicate,PeakType>(p));
		}

		/// See std::list documentation.
		void sort()
		{
			list_.sort(PointerComparator<typename PeakType::PositionLess>());
		}

		/// See std::list documentation.
		template <class BinaryPredicate>
		void sort(BinaryPredicate p)
		{
			list_.sort(DPeakListBinaryPredicate<BinaryPredicate,PeakType>(p));
		}		

		/// See std::list documentation.
		void reverse()
		{
			list_.reverse();
		}

		/// See std::list documentation.
		template <class InputIterator>
		void assign(InputIterator f , InputIterator l)
		{
			clear();
			insert(end(),f,l);
		}

		/// See std::list documentation.
		void assign(size_type n , const PeakType& x)
		{
			clear();
			insert(end(),n,x);
		}
		
		/// See std::list documentation.
		allocator_type get_allocator() const
		{
			return list_.get_allocator();
		}

		/**	@name Sorting.
				These simplified sorting methods are supported in addition to	
				the standard sorting methods of std::list.
		*/
		//@{
		/// Sort the List according to pean intensities
		void sortByIntensity() 
		{ 
			sort(typename PeakType::IntensityLess()); 
		}

		/// Lexicographically sorts the peaks by their position.
		void sortByPosition() 
		{ 
			sort(typename PeakType::PositionLess() );
		}

		/**
			@brief Sorts the peaks by one dimension of their position.
			
			It is only sorted according to dimentsion @p i . 
		*/
		void sortByNthPosition(UnsignedInt i) throw (Exception::NotImplemented);
		//@}

		/** 
			@name Generic sorting function templates.
			Any peak comparator can be given as template argument.  You can also give the comparator as an
			argument to the function template (this is useful if the comparator is
			not default constructed, but keep in mind that STL copies comparators
			a lot).
			
			<p> Thus your can e.g. write <code>peaks.sortByComparator <
			DPeak<1>::IntensityLess > ()</code>, if peaks has type
			<code>DPeakArray < 1, DPeak <1> ></code>.
		*/
		//@{
		template < typename ComparatorType >
		void sortByComparator ( ComparatorType const & comparator )
		{ 
			sort(ComparatorType( comparator ) ); 
		}
		template < typename ComparatorType >
		void sortByComparator ()
		{ 
			sort(ComparatorType()); 
		}
		//@}
	
		///PersistentObject interface
		virtual void persistentWrite(PersistenceManager& pm, const char* name=0) const throw (Exception::Base)
		{
			pm.writeObjectHeader(this,name);
			//TODO Persistence
			pm.writeObjectTrailer(name);
		}
		
		///PersistentObject interface
		virtual void persistentRead(PersistenceManager& pm) throw (Exception::Base)
		{
			//TODO Persistence
			int dummy;
			pm.readPrimitive(dummy,"dummy_");
		}

		protected:
			///PersistentObject interface
	    virtual void clearChildIds_()
	    {
	    	//TODO	
	    };	


	  /// The actual list of peak pointers
		std::list<PeakType*> list_;

		///A wrapper class that adapts a binary predicate to the internal Structure of DPeakList
		template <class BinaryPredicate, class BinaryPredicateArgumentType>
		class DPeakListBinaryPredicate
		{
			public:
			
			DPeakListBinaryPredicate(BinaryPredicate& bp)
			{
					pred_ = &bp;
			}
			
			bool operator () (BinaryPredicateArgumentType* a,BinaryPredicateArgumentType* b)
			{
				return pred_->operator()(*a,*b);
			}
			
			protected:
			BinaryPredicate* pred_;
		};
	
		///@name Serialization
		//@{
	 private:
		/// Serialization interface
		template<class Archive>
		void serialize(Archive & ar, const unsigned int /* version */ )
		{
			ar & boost::serialization::make_nvp("list",list_);
		}
		//@}

		/// Serialization
		friend class boost::serialization::access;
		
	};

	///Print the contents to a stream.
	template <Size D, typename Peak>
	std::ostream& operator << (std::ostream& os, const DPeakList<D, Peak>& array)
	{
		os << "-- DPEAKLIST BEGIN --"<<std::endl;
		for (typename DPeakList<D, Peak>::const_iterator it = array.begin(); it!=array.end(); ++it)
		{
			os << *it << std::endl;
		}
		os << "-- DPEAKLIST END --"<<std::endl;
		return os;
	}

//---------------------------------------------------------------
//  Implementation of the inline / template functions
//---------------------------------------------------------------

	template <Size D, typename PeakT > 
	void DPeakList<D,PeakT>::sortByNthPosition(UnsignedInt i) throw (Exception::NotImplemented)
	{ 
		OPENMS_PRECONDITION(i < Index(D), "illegal dimension")
		if (i==0)
		{
			sort(typename PeakType::template NthPositionLess<0>() );
		}
		else if (i==1)
		{
			sort(typename PeakType::template NthPositionLess<1>() );
		}
		else if (i==2)
		{
			sort(typename PeakType::template NthPositionLess<2>() );
		}
		else
		{
			throw Exception::NotImplemented(__FILE__,__LINE__,__FUNCTION__);
		}
	}

} // namespace OpenMS

#endif // OPENMS_KERNEL_DPEAKLIST_H
