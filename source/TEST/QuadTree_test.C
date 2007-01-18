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
// $Maintainer: Thomas Kadauke $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/QuadTree.h>

///////////////////////////

#include <stdlib.h>
#include <vector>
#include <iostream>

template <class ForwardIterator>
bool is_sorted(ForwardIterator first, ForwardIterator last)
{
  if (first == last)
    return true;

  ForwardIterator next = first;
  for (++next; next != last; first = next, ++next) {
    if (*next < *first)
      return false;
  }

  return true;
}

// template<class ForwardIterator>
// bool is_sorted(ForwardIterator first, ForwardIterator last)
// {
// 	ForwardIterator it = first;
// 	for (++it; it != last; ++first, ++it) {
// 		if (*first > *it)
// 			return false;
// 	}
// 	return true;
// }

START_TEST(QuadTree<>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

typedef QuadTree<FloatKernelTraits, int> Tree;
typedef Tree::AreaType Area;
typedef Tree::PointType Point;
Tree* quadtree_ptr;

CHECK(QuadTree(const AreaType& area))
	quadtree_ptr = new Tree(Area(0, 0, 100, 100));
	TEST_NOT_EQUAL(quadtree_ptr, 0)
	TEST_EQUAL(quadtree_ptr->begin() == quadtree_ptr->end(), true)
	TEST_EQUAL(quadtree_ptr->getArea() == Area(0, 0, 100, 100), true)
RESULT

CHECK(~QuadTree())
	delete quadtree_ptr;
RESULT

CHECK(void insert(const PointType& position, Data* data) throw (Exception::IllegalTreeOperation))
	quadtree_ptr = new Tree(Area(0, 0, 100, 100));
	quadtree_ptr->insert(Point(10, 10), new int(10));
	TEST_EQUAL(quadtree_ptr->begin() == quadtree_ptr->end(), false)
	int size = 0;
	for (Tree::iterator i = quadtree_ptr->begin(); i != quadtree_ptr->end(); ++i)
		size++;
	TEST_EQUAL(size, 1)

	delete quadtree_ptr;
RESULT

CHECK([EXTRA] multiple insert())
	quadtree_ptr = new Tree(Area(0, 0, 100, 100));
	for (int x = 5; x < 100; x += 10) {
		for (int y = 5; y < 100; y += 10) {
			quadtree_ptr->insert(Point(x, y), new int(x * y));
		}
	}
	TEST_EQUAL(quadtree_ptr->begin() == quadtree_ptr->end(), false)
	
	int size = 0;
	for (Tree::iterator i = quadtree_ptr->begin(); i != quadtree_ptr->end(); ++i)
		size++;
	TEST_EQUAL(size, 100)
	
	size = 0;
	for (Tree::iterator i = quadtree_ptr->begin(Area(0, 0, 50, 50)); i != quadtree_ptr->end(); ++i)
		size++;
	TEST_EQUAL(size, 25)
	
	delete quadtree_ptr;
RESULT

CHECK([EXTRA] random insert())
	quadtree_ptr = new Tree(Area(0, 0, 100, 100));
	Area area(10, 10, 20, 20);
	int inarea = 0;
	for (int i = 0; i != 10000; i++) {
		double x = double(rand() % 1000000) / 10000.0;
		double y = double(rand() % 1000000) / 10000.0;	
		Point pos(x, y);
		if (area.encloses(pos)) inarea++;
		quadtree_ptr->insert(pos, new int(rand()));
	}
	TEST_EQUAL(quadtree_ptr->begin() == quadtree_ptr->end(), false)
	
	int size = 0;
	for (Tree::iterator i = quadtree_ptr->begin(); i != quadtree_ptr->end(); ++i)
		size++;
	TEST_EQUAL(size>9990, true)
	
	size = 0;
	for (Tree::iterator i = quadtree_ptr->begin(area); i != quadtree_ptr->end(); ++i)
		size++;
	TEST_EQUAL(size, inarea)
	delete quadtree_ptr;
RESULT

CHECK(Iterator begin(const AreaType& area))
	Area area(0, 0, 100, 100);
	quadtree_ptr = new Tree(area);
	quadtree_ptr->insert(Point(10, 10), new int(10));
	TEST_EQUAL(quadtree_ptr->begin(area) == quadtree_ptr->begin(), false); // equality only when both iterators are end()
	TEST_EQUAL(quadtree_ptr->begin(area) != quadtree_ptr->end(), true);
	delete quadtree_ptr;
RESULT

CHECK(Iterator begin())
	Area area(0, 0, 100, 100);
	quadtree_ptr = new Tree(area);
	TEST_EQUAL(quadtree_ptr->begin() == quadtree_ptr->end(), true);
	quadtree_ptr->insert(Point(10, 10), new int(10));
	TEST_EQUAL(quadtree_ptr->begin() == quadtree_ptr->begin(area), false); // equality only when both iterators are end()
	TEST_EQUAL(quadtree_ptr->begin() != quadtree_ptr->end(), true);
	delete quadtree_ptr;
RESULT

CHECK(Iterator end())
	quadtree_ptr = new Tree(Area(0, 0, 100, 100));
	Tree::Iterator end = quadtree_ptr->end();
	quadtree_ptr->insert(Point(10, 10), new int(10));
	TEST_EQUAL(quadtree_ptr->end() == end, true);
	delete quadtree_ptr;
RESULT

CHECK(ConstIterator begin(const AreaType& area) const)
	Area area(0, 0, 100, 100);
	quadtree_ptr = new Tree(area);
	
	const Tree* const const_tree = quadtree_ptr;
	TEST_EQUAL(const_tree->begin(area) == const_tree->end(), true);
	quadtree_ptr->insert(Point(10, 10), new int(10));
	
	TEST_EQUAL(const_tree->begin(area) == const_tree->begin(), false); // equality only when both iterators are end()
	TEST_EQUAL(const_tree->begin(area) != const_tree->end(), true);
	
	delete quadtree_ptr;
RESULT

CHECK(Iterator begin())
	Area area(0, 0, 100, 100);
	quadtree_ptr = new Tree(area);
	
	const Tree* const const_tree = quadtree_ptr;
	TEST_EQUAL(const_tree->begin() == const_tree->end(), true);
	
	quadtree_ptr->insert(Point(10, 10), new int(10));
	TEST_EQUAL(const_tree->begin() == const_tree->begin(area), false); // equality only when both iterators are end()
	TEST_EQUAL(const_tree->begin() != const_tree->end(), true);
	delete quadtree_ptr;
RESULT

CHECK(ConstIterator end() const)
	quadtree_ptr = new Tree(Area(0, 0, 100, 100));
	
	const Tree* const const_tree = quadtree_ptr;
	Tree::ConstIterator end = const_tree->end();
	
	quadtree_ptr->insert(Point(10, 10), new int(10));
	TEST_EQUAL(const_tree->end() == end, true);
	
	delete quadtree_ptr;
RESULT

CHECK(const AreaType& getArea() const)
	Area area(0, 0, 100, 100);
	quadtree_ptr = new Tree(area);
	TEST_EQUAL(quadtree_ptr->getArea() == area, true);
	delete quadtree_ptr;
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
