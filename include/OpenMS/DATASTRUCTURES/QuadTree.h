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
// $Maintainer: Thomas Kadauke $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_QUADTREE_H
#define OPENMS_DATASTRUCTURES_QUADTREE_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>

#include <vector>
#include <limits>
#include <iostream>

namespace OpenMS
{
	/**
		@brief Namespace to hide implementation details
		
		This namespace is to be used for classes that are only to be used by some
		parent class, but not by the user directly.
	*/
	namespace Internal
	{
		///Node of a QuadTree (internal use only)	
		template<typename Traits, typename Data>
		class QuadNode
		{
		private:
			inline bool isInf();
			
		public:
			typedef DPosition<2,Traits> PointType;
			typedef DRange<2,Traits> AreaType;
			typedef std::pair<PointType, Data*> ValueType;

			inline QuadNode();
			inline QuadNode(PointType pos, Data* data);

			inline bool isLeaf();
			inline bool isInner();
			inline bool isNil();

			inline bool isEmpty();

			inline QuadNode* getChildren();
			inline Data* getData();

			inline void resetPosition();

			// set point to (inf, inf) if this node is an inner node
			ValueType data;
		};
		
		///Chunk of a QuadTree (internal use only)			
		template<class T>
		class Chunk
		{
		public:
			inline Chunk(UnsignedInt size);
			inline ~Chunk();

			inline T* get(UnsignedInt num);

			inline void setNext(Chunk<T>* next);

			inline UnsignedInt size();

		private:
			Chunk<T>* next_;
			T* data_;
			T* current_;
			UnsignedInt size_;
		};

		///Allocator for the QuadTree (internal use only)			
		template<class T>
		class Alloc
		{
		public:
			inline Alloc();
			inline ~Alloc();

			inline T* getNodes();

		private:
			Chunk<T>* first_inner_;
			Chunk<T>* last_inner_;
		};

		///Pointer for the QuadTree (internal use only)
		template<class T>
		struct RemovePtr
		{
			typedef T Type;
		};


		///Pointer for the QuadTree (internal use only)
		template<class T>
		struct RemovePtr<T*>
		{
			typedef T Type;
		};
	}


	/**
		@brief Iterator on a QuadTree
		
		This class provides an interface to the contents of a
		QuadTree, which lie in a specified area. This interface
		aims to be compatible with the STL iterator interfaces.
		You should not create instances of this class directly;
		rather, use QuadTree::begin() and QuadTree::end() to
		construct an AreaIterator. You can use iterators like
		follows:
		
		@code
		QuadTree<double, double> tr(...);
		...
		for (QuadTree<double, double>::Iterator i
		     = tr.begin(AreaType(1, 1, 20, 20));
		     i != tr.end(); ++i) {
		  cout << *i << endl;
		}
		@endcode
		
		@param Type The items' value type.
		@param Ref The value type's reference type.
		@param Ptr The value type's pointer type.
	*/
	template<class Type, class Ref, class Ptr>
	class AreaIterator
	{
	private:
		typedef typename Type::first_type::CoordinateType Coord;
		typedef typename Internal::RemovePtr<typename Type::second_type>::Type Data;

	public:
	
		/// the default traits
		typedef typename Type::first_type::TraitsType Traits;
		
		/// the default area type
		typedef DRange<2,Traits> AreaType;
		/**
			Refers to an AreaIterator using non-const semantics.
			Used for STL-Compliance.
			@sa const_iterator, Iterator, ConstIterator
		*/
		typedef AreaIterator<Type, Type&, Type*> iterator;

		/**
			Refers to an AreaIterator using non-const semantics.
			Used for STL-Compliance.
			@sa iterator, Iterator, ConstIterator
		*/
		typedef AreaIterator<Type, const Type&, const Type*> const_iterator;

		/**
			This type is actually not used in this class, but
			required by some STL algorithms.
			@sa difference_type, iterator_category
		*/
		typedef std::size_t size_type;

		/**
			This type is actually not used in this class, but
			required by some STL algorithms.
			@sa size_type, iterator_category
		*/
		typedef std::ptrdiff_t difference_type;

		/**
			This type is actually not used in this class, but
			required by some STL algorithms. It hints algorithms
			what they can expect from this iterator class.
			@sa size_type, difference_type
		*/
		typedef std::forward_iterator_tag iterator_category;

		/**
			This is the standard interface to the Type template
			parameter. It is used by most STL algorithms.
			@sa reference, pointer, ValueType
		*/
		typedef Type value_type;

		/**
			This is the standard interface to the Ref template
			parameter, denoting the iterator's reference type.
			It is used by most STL algorithms.
			@sa value_type, pointer, Reference
		*/
		typedef Ref reference;

		/**
			This is the standard interface to the Type template
			parameter, denoting the iterator's pointer type.
			It is used by most STL algorithms.
			@sa reference, pointer, Pointer
		*/
		typedef Ptr pointer;

		/**
			This type is provided, because the OpenMS coding
			convention states that types start with a capital
			letter, in contrast to the STL types.
			@sa iterator, ConstIterator
		*/
		typedef iterator Iterator;

		/**
			This type is provided, because the OpenMS coding
			convention states that types start with a capital
			letter, in contrast to the STL types.
			@sa const_iterator, Iterator
		*/
		typedef const_iterator ConstIterator;

		/**
			This type is provided, because the OpenMS coding
			convention states that types start with a capital
			letter, in contrast to the STL types.
			@sa value_type
		*/
		typedef value_type ValueType;

		/**
			This type is provided, because the OpenMS coding
			convention states that types start with a capital
			letter, in contrast to the STL types.
			@sa reference
		*/
		typedef reference Reference;

		/**
			This type is provided, because the OpenMS coding
			convention states that types start with a capital
			letter, in contrast to the STL types.
			@sa pointer
		*/
		typedef pointer Pointer;

	private:
		typedef Internal::QuadNode<Traits, Data>* NodePointer;

	public:
		/**
			Default constructor. Constructs an iterator pointing
			to nowhere. This is used by the QuadTree::end() function.
			You should not use this constructor. Use the QuadTree
			iterator interface instead.
		*/
		inline AreaIterator();

		/**
			Constructor. Creates an iterator which can be used to
			extract all points from a QuadTree in a specified area.
			You should not use this constructor. Use the QuadTree
			iterator interface instead.
			@param tree_area The overall area of the quadtree.
			@param area The area where the desired points lie in.
			@param root_node The root of the QuadTree.
		*/
		inline AreaIterator(const AreaType& tree_area, const AreaType& area, NodePointer root_node);
		
		inline AreaIterator(const AreaIterator& it);
		
		inline AreaIterator& operator=(const AreaIterator& it);

		/**
			Returns a reference to the data pointed to by the iterator.
			Note that for ConstIterators, this is a const reference.
			This function's usage is the same as for STL iterators.
			@return A reference to the data pointed to by the iterator.
		*/
		inline Reference operator*() const;

		/**
			Returns a pointer to the data pointed to by the iterator.
			Note that for ConstIterators, this is a const pointer.
			This function's usage is the same as for STL iterators.
			This is used if the Type template parameter is an
			aggregate type. Example:
			*
			@code
			struct some_struct {
			  int data;
			  void print() { cout << data; }
			};
			...
			QuadTree<double, some_struct> tr(...);
			...
			// calls the some_struct::print() function
			// for every element in the tree, bounded by
			// the specified area.
			for (QuadTree<double, some_struct>::Iterator i
			     = tr.begin(AreaType(1, 1, 20, 20));
			     i != tr.end(); ++i) {
			  i->print();
			}
			@endcode
			*
			@return A pointer to the data pointed to by the iterator.
		*/
		inline Pointer operator->() const;

		/**
			Moves the iterator to next item stored in the tree, the
			position of which is in the specified area. Note that
			there are no guarantees about the order in which the
			items are discovered.
		*/
		inline AreaIterator& operator++();

#if 0
		/**
			This is the post-increment operator. It behaves like the
			pre-increment operator, except that it copies the iterator,
			which can be costly. Unless you really need to copy the
			iterator, use the pre-increment version instead.
		*/
		inline AreaIterator operator++(int);
#endif
		/**
			Checks this iterator and the @p it iterator for equality.
			Two iterators are considered equal, if they point to the
			same element in the tree, independent of the area which
			they iterate over. This means that it1 == it2 does not imply
			++it1 == ++it2.
		*/
		inline bool operator==(const AreaIterator& it);

		/**
			Checks this iterator and the @p it iterator for inequality.
			Two iterators are considered inequal, if they do not point
			to the same element in the tree, independent of the area which
			they iterate over. This means that it1 != it2 does not imply
			++it1 != ++it2.
		*/
		inline bool operator!=(const AreaIterator& it);

	private:
		void findNodes(const AreaType& area, NodePointer current, const AreaType& nodeArea);
		
		std::vector<NodePointer> nodes_;
		typename std::vector<NodePointer>::iterator current_;
	
#if 0
		AreaType area_;
		NodePointer current_;

		struct StackEntry_ {
			StackEntry_() : pointer(0), child(0) {}
			StackEntry_(NodePointer p, int ch, const AreaType& ar)
				: pointer(p), child(ch), area(ar) {}

			NodePointer pointer;
			int child;
			AreaType area;
		};

		std::stack<StackEntry_> stack_;
#endif
	};

	/**
		@brief QuadTree implementation for fast access to points in a plane.
		
		This class provides a QuadTree, which is essentially a data
		structure which stores points with 2D coordinates in the plane
		and provides fast access to the points.
		
		@param Data The type of the values associated with each point.
		@param Coord The points coordinate type.
		
		@ingroup Datastructures
	*/
	template<typename Traits, typename Data>
	class QuadTree
	{
	private:
		struct Cmp;

	public:
		///the default point type	
		typedef DPosition<2,Traits> PointType;

		///the default area type	
		typedef DRange<2,Traits> AreaType;
		
		/**
			This is the standard interface to the Type template
			parameter. It is used by most STL algorithms.
			@sa reference, pointer, ValueType
		*/
		typedef typename Internal::QuadNode<Traits, Data>::ValueType value_type;

		/**
			This is the standard interface to the Ref template
			parameter, denoting the iterator's reference type.
			It is used by most STL algorithms.
			@sa value_type, pointer, Reference
		*/
		typedef value_type& reference;

		/**
			This is the standard interface to the Type template
			parameter, denoting the iterator's pointer type.
			It is used by most STL algorithms.
			@sa reference, pointer, Pointer
		*/
		typedef value_type* pointer;

		/**
			Refers to an AreaIterator using non-const semantics.
			Used for STL-Compliance.
			@sa const_iterator, Iterator, ConstIterator
		*/
		typedef AreaIterator<value_type, value_type&, value_type*> iterator;

		/**
			Refers to an AreaIterator using const semantics.
			Used for STL-Compliance.
			@sa iterator, Iterator, ConstIterator
		*/
		typedef AreaIterator<value_type, const value_type&, const value_type*> const_iterator;

		/**
			This type is provided, because the OpenMS coding
			convention states that types start with a capital
			letter, in contrast to the STL types.
			@sa iterator, ConstIterator
		*/
		typedef iterator Iterator;

		/**
			This type is provided, because the OpenMS coding
			convention states that types start with a capital
			letter, in contrast to the STL types.
			@sa const_iterator, Iterator
		*/
		typedef const_iterator ConstIterator;

		/**
			This type is provided, because the OpenMS coding
			convention states that types start with a capital
			letter, in contrast to the STL types.
			@sa value_type
		*/
		typedef value_type ValueType;

		/**
			This type is provided, because the OpenMS coding
			convention states that types start with a capital
			letter, in contrast to the STL types.
			@sa reference
		*/
		typedef reference Reference;

		/**
			This type is provided, because the OpenMS coding
			convention states that types start with a capital
			letter, in contrast to the STL types.
			@sa pointer
		*/
		typedef pointer Pointer;

	private:
		///Comparator for QuadTree (internal use only)	
		struct Cmp
		{
			typedef typename QuadTree<Traits, Data>::value_type Param;

			inline bool operator()(const Param& a, const Param& b)
			{
				return a.second < b.second;
			}
		};

		typedef Internal::QuadNode<Traits, Data>* NodePointer;

	public:
		/**
			Constructs an empty QuadTree which is bordered by
			@p area. All points inserted to the tree MUST be
			inside the borders. Disobedience will be penalized
			with endless loops.
			@param area The QuadTree's borders.
		*/
		inline QuadTree(const AreaType& area);

		/**
			Destructor. Destructs all contents of the tree,
			so care must be taken when the construtor is
			called.
		*/
		inline ~QuadTree();

		/**
			@brief Inserts a new point into the tree.
		
			If the new point is already in the QuadTree, an IllegalTreeOperation exception is thrown.
		
			@param position The points position.
			@param data The data associated with the point.
		*/
		inline void insert(const PointType& position, /*const */Data* data) throw (Exception::IllegalTreeOperation);

		/**
			Constructs a mutable AreaIterator pointing to the first
			point in the tree which lies in the specified @p area.
			@param area The area to iterate over.
			@return The iterator pointing to the first point.
			@sa AreaIterator
		*/
		inline Iterator begin(const AreaType& area);

		/**
			Constructs a mutable AreaIterator pointing to the first
			point in the tree.
			@return The iterator pointing to the first point.
			@sa AreaIterator
		*/
		inline Iterator begin();
    
		/**
			Constructs a mutable AreaIterator pointing to nowhere
			(i.e. the element past the last element, in
			compliance to STL iterators). Every iterator which
			gets incremented over and over once points to
			nowhere and hence is equal to the iterator returned
			by this function.
			@return An iterator pointing to nowhere.
		*/
		inline Iterator end();

		/**
			This is the const version of begin(). Data associated
			with points that can be accessed through a ConstIterator
			is not mutable.
			@sa begin()
		*/
		inline ConstIterator begin(const AreaType& area) const;

		/**
			This is the const version of end().
			@sa end()
		*/
		inline ConstIterator end() const;

		/**
			Returns the QuadTree's area.
			@return The tree's area.
		*/
		inline const AreaType& getArea() const;

	private:
		void insert_(NodePointer node, const AreaType& area, const PointType& position, /*const */Data* data) throw (Exception::IllegalTreeOperation);

		Internal::Alloc<Internal::QuadNode<Traits, Data> > alloc_;
		NodePointer root_;
		AreaType area_;
	};

}

//---------------------------------------------------------------
//  Implementation of the inline / template functions
//---------------------------------------------------------------

// OpenMS::Internal::QuadNode

template<typename Traits, typename Data>
inline bool OpenMS::Internal::QuadNode<Traits, Data>::isInf()
{
	return data.first.X() == std::numeric_limits<typename Traits::CoordinateType>::infinity()
	    && data.first.Y() == std::numeric_limits<typename Traits::CoordinateType>::infinity();
}

template<typename Traits, typename Data>
inline OpenMS::Internal::QuadNode<Traits, Data>::QuadNode()
	: data(ValueType(PointType(std::numeric_limits<typename Traits::CoordinateType>::infinity(), std::numeric_limits<typename Traits::CoordinateType>::infinity()), 0))
{

}

template<typename Traits, typename Data>
inline OpenMS::Internal::QuadNode<Traits, Data>::QuadNode(PointType pos, Data* data)
	: data(ValueType(pos, data))
{

}

template<typename Traits, typename Data>
inline bool OpenMS::Internal::QuadNode<Traits, Data>::isLeaf()
{
	return !isInf() && data.second;
}

template<typename Traits, typename Data>
inline bool OpenMS::Internal::QuadNode<Traits, Data>::isInner()
{
	return isInf() && data.second;
}

template<typename Traits, typename Data>
inline bool OpenMS::Internal::QuadNode<Traits, Data>::isNil()
{
	return isInf() && !data.second;
}

template<typename Traits, typename Data>
inline bool OpenMS::Internal::QuadNode<Traits, Data>::isEmpty()
{
	for (int i = 0; i != 4; i++)
	{
		if (!getChildren()[i].isNil()) return false;
	}
	return true;
}

template<typename Traits, typename Data>
inline OpenMS::Internal::QuadNode<Traits, Data>* OpenMS::Internal::QuadNode<Traits, Data>::getChildren()
{
	return reinterpret_cast<QuadNode*>(data.second);
}

template<typename Traits, typename Data>
inline Data* OpenMS::Internal::QuadNode<Traits, Data>::getData()
{
	return data.second;
}

template<typename Traits , typename Data>
inline void OpenMS::Internal::QuadNode<Traits, Data>::resetPosition()
{
	data.first = PointType(std::numeric_limits<typename Traits::CoordinateType>::infinity(), std::numeric_limits<typename Traits::CoordinateType>::infinity());
}

// OpenMS::Internal::Chunk

template<class T>
inline OpenMS::Internal::Chunk<T>::Chunk(UnsignedInt size)
	: next_(0), size_(size)
{
	data_ = current_ = reinterpret_cast<T*>(std::malloc(size_ * sizeof(T)));
}

template<class T>
inline OpenMS::Internal::Chunk<T>::~Chunk()
{
	if (next_) delete next_;
	for (T* i = data_; i != current_; i++) i->~T();
	std::free(data_);
}

template<class T>
inline T* OpenMS::Internal::Chunk<T>::get(UnsignedInt num)
{
	if ((current_ - data_) + num <= size_)
	{
		for (UnsignedInt i = 0; i != num; i++)
		{
			new(current_++) T();
		}
		return current_ - num;
	}

	return 0;
}

template<class T>
inline void OpenMS::Internal::Chunk<T>::setNext(Chunk<T>* next)
{
	next_ = next;
}

template<class T>
inline OpenMS::UnsignedInt OpenMS::Internal::Chunk<T>::size()
{
	return size_;
}

// OpenMS::Internal::Alloc

template<class T>
inline OpenMS::Internal::Alloc<T>::Alloc()
{
	first_inner_ = last_inner_ = new Chunk<T>(1024);
}

template<class T>
inline OpenMS::Internal::Alloc<T>::~Alloc()
{
	delete first_inner_;
}

template<class T>
inline T* OpenMS::Internal::Alloc<T>::getNodes()
{
	T* nodes = last_inner_->get(4);
	if (!nodes)
	{
		Chunk<T>* c = new Chunk<T>(last_inner_->size() << 1);
		last_inner_->setNext(c);
		last_inner_ = c;
		return getNodes();
	}

	return nodes;
}

// OpenMS::AreaIterator

template<class Type, class Ref, class Ptr>
inline OpenMS::AreaIterator<Type, Ref, Ptr>::AreaIterator()
	: current_(nodes_.end())
{

}

template<class Type, class Ref, class Ptr>
inline OpenMS::AreaIterator<Type, Ref, Ptr>::AreaIterator(const DRange<2,AreaIterator<Type, Ref, Ptr>::Traits>& tree_area, const DRange<2,AreaIterator<Type, Ref, Ptr>::Traits>& area, NodePointer root_node)
{
#if 0
	// bad hack ;-)
	stack_.push(StackEntry_(root_node, -1, tree_area));
	operator++();
#endif

	findNodes(area, root_node, tree_area);
	current_ = nodes_.begin();
}

template<class Type, class Ref, class Ptr>
inline OpenMS::AreaIterator<Type, Ref, Ptr>::AreaIterator(const OpenMS::AreaIterator<Type, Ref, Ptr>& it)
	: nodes_(it.nodes_), current_(nodes_.begin())
{

}

template<class Type, class Ref, class Ptr>
inline OpenMS::AreaIterator<Type, Ref, Ptr>&
	OpenMS::AreaIterator<Type, Ref, Ptr>::operator=(const OpenMS::AreaIterator<Type, Ref, Ptr>& it)
{
	nodes_ = it.nodes_;
	current_ = nodes_.begin();
}

template<class Type, class Ref, class Ptr>
inline typename OpenMS::AreaIterator<Type, Ref, Ptr>::Reference
	OpenMS::AreaIterator<Type, Ref, Ptr>::operator*() const
{
	return (*current_)->data;
}

template<class Type, class Ref, class Ptr>
inline typename OpenMS::AreaIterator<Type, Ref, Ptr>::Pointer
	OpenMS::AreaIterator<Type, Ref, Ptr>::operator->() const
{
	return &(operator*());
}

#if 0
template<class Type, class Ref, class Ptr>
inline OpenMS::AreaIterator<Type, Ref, Ptr>&
	OpenMS::AreaIterator<Type, Ref, Ptr>::operator++()
{
	// OK, this function is a little tricky. The idea is
	// instead of using a recursive function to traverse
	// the tree and copying all discovered points into
	// a list, recursion is simulated by a stack, which
	// is stored in the iterator object, so it is persistent
	// between successive calls to this function.

	// The stack stores the path to the current leaf, i.e.
	// a node and the number of the node's child for the
	// next step of the path.
	StackEntry_ top;

	// This outermost loop is needed since not every leaf
	// discovered by the inner loops is in the area (see
	// below)
	do
	{
		// traverse upwards

		// each iteration in this loop corresponds to a return
		// statement in the recursive approach (hence, values
		// are popped from the stack)
		do
		{
			// if the stack is empty, all leafs in the area have
			// already been visited.
			if (stack_.empty())
			{
				current_ = 0;
				return *this;
			}

			top = stack_.top();

			if (top.pointer->isLeaf())
			{
				stack_.pop();
				top = stack_.top();
			}

			// follow the stack upwards until a path to an
			// unvisited leaf is discovered
			stack_.pop();

			// check all of the current node's children (if
			// it is not a leaf)
			++top.child;
			while (top.child != 4)
			{
				// break if all children are discovered or if an
				// undiscovered child whose area intersects with
				// the iterator's area is found
				if (!top.pointer->getChildren()[top.child].isNil() && area_.isIntersected(top.area)) break;
				++top.child;
			}
		} while (top.child == 4);

		// traverse downwards

		// there must exist another leaf in the tree
		// (otherwise the stack would have been empty
		// and this function would have been left)
		// So we are following the leftmost non-visited
		// path.

		typename Traits::CoordinateType mid_x = (top.area.minX() + top.area.maxX()) / 2.0;
		typename Traits::CoordinateType mid_y = (top.area.minY() + top.area.maxY()) / 2.0;

		// Each iteration in the following loop corresponds
		// to a recursive call (and hence, values are pushed
		// onto the stack)
		while (true)
		{
			stack_.push(top);

			// follow child # top.child and calculate area
			switch (top.child)
			{
			case 0:
				top = StackEntry_(top.pointer->getChildren(), 0, AreaType(top.area.minX(), top.area.minY(), mid_x, mid_y));
				break;
			case 1:
				top = StackEntry_(top.pointer->getChildren() + 1, 0, AreaType(mid_x, top.area.minY(), top.area.maxX(), mid_y));
				break;
			case 2:
				top = StackEntry_(top.pointer->getChildren() + 2, 0, AreaType(mid_x, mid_y, top.area.maxX(), top.area.maxY()));
				break;
			case 3:
				top = StackEntry_(top.pointer->getChildren() + 3, 0, AreaType(top.area.minX(), mid_y, mid_x, top.area.maxY()));
				break;

			}

			// is this a leaf? if so, we've found what we're
			// looking for.
			if (top.pointer->isLeaf()) break;

			// choose the current node's left-most non-discovered
			// child for the next iteration
			while (top.pointer->getChildren()[top.child].isNil())
			{
				++top.child;
				// the following can not happen since this would
				// mean that there is no leaf left in the tree
			}
		}

		// update the pointer to the current leaf
		current_ = top.pointer;
		stack_.push(top);

		// if the areas have been intersecting, but this
		// point lies outside of the area (which can actually
		// happen), we have to find the next node.
	} while (!area_.encloses(current_->data.first));

	return *this;
}
#endif

template<class Type, class Ref, class Ptr>
inline void OpenMS::AreaIterator<Type, Ref, Ptr>::findNodes(const AreaType& area, NodePointer node, const AreaType& nodeArea)
{
	// nodeArea can't be const ref because intersects() is not const and can't be as of now
	if (node->isNil()) return;
	if (!nodeArea.intersects(area)) return;
	
	if (node->isInner())
	{
		typename Traits::CoordinateType mid_x = (nodeArea.minX() + nodeArea.maxX()) / 2.0;
		typename Traits::CoordinateType mid_y = (nodeArea.minY() + nodeArea.maxY()) / 2.0;
		
		findNodes(area, node->getChildren() + 0, AreaType(nodeArea.minX(), nodeArea.minY(), mid_x, mid_y));
		findNodes(area, node->getChildren() + 1, AreaType(mid_x, nodeArea.minY(), nodeArea.maxX(), mid_y));
		findNodes(area, node->getChildren() + 2, AreaType(mid_x, mid_y, nodeArea.maxX(), nodeArea.maxY()));
		findNodes(area, node->getChildren() + 3, AreaType(nodeArea.minX(), mid_y, mid_x, nodeArea.maxY()));
	}
	else if (node->isLeaf())
	{
		if (area.encloses(node->data.first))
		{
			nodes_.push_back(node);
		}
	}
}

template<class Type, class Ref, class Ptr>
inline OpenMS::AreaIterator<Type, Ref, Ptr>&
	OpenMS::AreaIterator<Type, Ref, Ptr>::operator++()
{
	++current_;
	return *this;
}

#if 0
template<class Type, class Ref, class Ptr>
inline OpenMS::AreaIterator<Type, Ref, Ptr>
	OpenMS::AreaIterator<Type, Ref, Ptr>::operator++(int)
{
	// copy the iterator, move the original and return the copy
	AreaIterator tmp = *this;
	operator++();
	return tmp;
}
#endif

template<class Type, class Ref, class Ptr>
inline bool OpenMS::AreaIterator<Type, Ref, Ptr>::operator==(const AreaIterator& it)
{
	return current_ == nodes_.end() && it.current_ == it.nodes_.end();
}

template<class Type, class Ref, class Ptr>
inline bool OpenMS::AreaIterator<Type, Ref, Ptr>::operator!=(const AreaIterator& it)
{
	return !(operator==(it));
}

// OpenMS::QuadTree

template<typename Traits, typename Data>
inline OpenMS::QuadTree<Traits, Data>::QuadTree(const AreaType& area)
	: area_(area)
{
	root_ = new Internal::QuadNode<Traits, Data>();
	root_->data.second = reinterpret_cast<Data*>(alloc_.getNodes());
}

template<typename Traits, typename Data>
inline OpenMS::QuadTree<Traits, Data>::~QuadTree()
{
	delete root_;
}

template<typename Traits, typename Data>
inline void OpenMS::QuadTree<Traits, Data>::insert(const OpenMS::QuadTree<Traits, Data>::PointType& position, /*const */Data* data) throw (Exception::IllegalTreeOperation)
{
	//std::cout << std::endl << std::endl << std::endl;
	insert_(root_, area_, position, data);
}

template<typename Traits, typename Data>
inline typename OpenMS::QuadTree<Traits, Data>::Iterator OpenMS::QuadTree<Traits, Data>::begin(const AreaType& area)
{
	if (root_->isEmpty()) return Iterator();
	return Iterator(area_, area, root_);
}

template<typename Traits, typename Data>
inline typename OpenMS::QuadTree<Traits, Data>::Iterator OpenMS::QuadTree<Traits, Data>::begin()
{
	if (root_->isEmpty()) return Iterator();
	typename Traits::CoordinateType min = std::numeric_limits<typename Traits::CoordinateType>::min();
	typename Traits::CoordinateType max = std::numeric_limits<typename Traits::CoordinateType>::max(); 
	AreaType area(min, min, max, max);
	return Iterator(area_, area, root_);
}

template<typename Traits, typename Data>
inline typename OpenMS::QuadTree<Traits, Data>::Iterator OpenMS::QuadTree<Traits, Data>::end()
{
	return Iterator();
}

template<typename Traits, typename Data>
inline typename OpenMS::QuadTree<Traits, Data>::ConstIterator OpenMS::QuadTree<Traits, Data>::begin(const AreaType& area) const
{
	return ConstIterator(area_, area, root_);
}

template<typename Traits, typename Data>
inline typename OpenMS::QuadTree<Traits, Data>::ConstIterator OpenMS::QuadTree<Traits, Data>::end() const
{
	return ConstIterator();
}

template<typename Traits, typename Data>
inline const typename OpenMS::QuadTree<Traits, Data>::AreaType& OpenMS::QuadTree<Traits, Data>::getArea() const
{
	return area_;
}

template<typename Traits, typename Data>
void OpenMS::QuadTree<Traits, Data>::insert_(NodePointer node, const AreaType& area, const PointType& position, /*const */Data* data) throw (Exception::IllegalTreeOperation)
{
	// calculate the current node's mass point
	typename Traits::CoordinateType mid_x = (area.minX() + area.maxX()) / 2.0;
	typename Traits::CoordinateType mid_y = (area.minY() + area.maxY()) / 2.0;
  
  //std::cout << "Area: " << mid_x << " " << mid_y << std::endl;
  
	if (node->isLeaf())              
	{
		//std::cout << "Leaf" << std::endl;
		// if this is a leaf, we need to split the current node
		// store the leaf's position and data for later re-insertion
		PointType old_pos = node->data.first;
		Data* old_data = node->getData();

		// mark node as inner node
		node->resetPosition();
		// store pointer to the first child
		node->data.second = reinterpret_cast<Data*>(alloc_.getNodes());

		if (position==old_pos)
		{
			throw Exception::IllegalTreeOperation(__FILE__,__LINE__,__PRETTY_FUNCTION__);
		}
		// insert the node's original data in a child of node
		insert_(node, area, old_pos, old_data);
		// insert the new node's data in a child of node
		insert_(node, area, position, data);
	} 
	else 
	{
		//std::cout << "Not Leaf" << std::endl;
		// check which child node to follow (or create)
		// (I'm sure this can be optimized by size)
		if (position.Y() < mid_y)
		{
			//std::cout << "position.Y() < mid_y" << std::endl;
			if (position.X() < mid_x)
			{
				//std::cout << "position.X() < mid_x" << std::endl;
				// left-top child (node->getChildren()[0])
				if (!node->getChildren()[0].isNil())
				{
					// if the child exists, redirect the insert_() call to
					// the child.
					insert_(node->getChildren(), AreaType(area.minX(), area.minY(), mid_x, mid_y), position, data);
				}
				else
				{
					// if not, change node to a leaf
					node->getChildren()[0].data.first = position;
					node->getChildren()[0].data.second = data;
				}
			}
			else
			{
				//std::cout << "position.X() >= mid_x" << std::endl;
				// right-top child (node->getChildren()[1])
				if (!node->getChildren()[1].isNil())
				{
					insert_(node->getChildren() + 1, AreaType(mid_x, area.minY(), area.maxX(), mid_y), position, data);
				}
				else
				{
					node->getChildren()[1].data.first = position;
					node->getChildren()[1].data.second = data;
				}
			}
		}
		else
		{
			//std::cout << "position.Y() >= mid_y" << std::endl;
			if (position.X() > mid_x)
			{
				//std::cout << "position.X() < mid_x" << std::endl;
				// right-bottom child (node->getChildren()[2])
				if (!node->getChildren()[2].isNil())
				{
					insert_(node->getChildren() + 2, AreaType(mid_x, mid_y, area.maxX(), area.maxY()), position, data);
				}
				else
				{
					node->getChildren()[2].data.first = position;
					node->getChildren()[2].data.second = data;
				}
			}
			else
			{
				//std::cout << "position.X() >= mid_x" << std::endl;
				// left-bottom child (node->getChildren()[3])
				if (!node->getChildren()[3].isNil())
				{
					//std::cout << "NIL" << std::endl;
					insert_(node->getChildren() + 3, AreaType(area.minX(), mid_y, mid_x, area.maxY()), position, data);
				}
				else
				{
					//std::cout << "NOT NIL" << std::endl;
					node->getChildren()[3].data.first = position;
					node->getChildren()[3].data.second = data;
				}
			}
		}
	}
	//std::cout << std::endl;
}

#endif //OPENMS_DATASTRUCTURES_QUADTREE_H
