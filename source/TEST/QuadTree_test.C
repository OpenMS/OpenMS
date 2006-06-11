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
// $Id: QuadTree_test.C,v 1.7 2006/06/08 14:29:18 marc_sturm Exp $
// $Author: marc_sturm $
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

START_TEST(QuadTree<>, "$Id: QuadTree_test.C,v 1.7 2006/06/08 14:29:18 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

typedef QuadTree<OpenMS::FloatKernelTraits, int> Tree;
typedef Tree::AreaType Area;
typedef Tree::PointType Point;
Tree* quadtree_ptr;

CHECK(constructor QuadTree())
	quadtree_ptr = new Tree(Area(0, 0, 100, 100));
	TEST_NOT_EQUAL(quadtree_ptr, 0)
	TEST_EQUAL(quadtree_ptr->begin() == quadtree_ptr->end(), true)
	TEST_EQUAL(quadtree_ptr->sortedBegin(quadtree_ptr->getArea()) == quadtree_ptr->sortedEnd(), true)
	TEST_EQUAL(quadtree_ptr->getArea() == Area(0, 0, 100, 100), true)
RESULT

CHECK(~QuadTree())
	delete quadtree_ptr;
RESULT

CHECK(insert())
	quadtree_ptr = new Tree(Area(0, 0, 100, 100));
	quadtree_ptr->insert(Point(10, 10), new int(10));
	TEST_EQUAL(quadtree_ptr->begin() == quadtree_ptr->end(), false)
	TEST_EQUAL(quadtree_ptr->sortedBegin(quadtree_ptr->getArea()) == quadtree_ptr->sortedEnd(), false)
	int size = 0;
	for (Tree::iterator i = quadtree_ptr->begin(); i != quadtree_ptr->end(); ++i)
		size++;
	TEST_EQUAL(size, 1)
	
	size = 0;
	for (Tree::SortedIterator i = quadtree_ptr->sortedBegin(quadtree_ptr->getArea()); i != quadtree_ptr->sortedEnd(); ++i)
		size++;
	TEST_EQUAL(size, 1)
	delete quadtree_ptr;
RESULT

CHECK(multiple insert())
	quadtree_ptr = new Tree(Area(0, 0, 100, 100));
	for (int x = 5; x < 100; x += 10) {
		for (int y = 5; y < 100; y += 10) {
			quadtree_ptr->insert(Point(x, y), new int(x * y));
		}
	}
	TEST_EQUAL(quadtree_ptr->begin() == quadtree_ptr->end(), false)
	TEST_EQUAL(quadtree_ptr->sortedBegin(quadtree_ptr->getArea()) == quadtree_ptr->sortedEnd(), false)
	
	int size = 0;
	for (Tree::iterator i = quadtree_ptr->begin(); i != quadtree_ptr->end(); ++i)
		size++;
	TEST_EQUAL(size, 100)
	
	size = 0;
	for (Tree::iterator i = quadtree_ptr->begin(Area(0, 0, 50, 50)); i != quadtree_ptr->end(); ++i)
		size++;
	TEST_EQUAL(size, 25)
	
	size = 0;
	for (Tree::SortedIterator i = quadtree_ptr->sortedBegin(quadtree_ptr->getArea()); i != quadtree_ptr->sortedEnd(); ++i)
		size++;
	TEST_EQUAL(size, 100)
	
	size = 0;
	for (Tree::SortedIterator i = quadtree_ptr->sortedBegin(Area(0, 0, 50, 50)); i != quadtree_ptr->sortedEnd(); ++i)
		size++;
	TEST_EQUAL(size, 25)
	
	delete quadtree_ptr;
RESULT

CHECK(random insert())
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
	TEST_EQUAL(quadtree_ptr->sortedBegin(quadtree_ptr->getArea()) == quadtree_ptr->sortedEnd(), false)
	
	int size = 0;
	for (Tree::iterator i = quadtree_ptr->begin(); i != quadtree_ptr->end(); ++i)
		size++;
	TEST_EQUAL(size>9990, true)
	
	size = 0;
	for (Tree::iterator i = quadtree_ptr->begin(area); i != quadtree_ptr->end(); ++i)
		size++;
	TEST_EQUAL(size, inarea)

	size = 0;
	for (Tree::SortedIterator i = quadtree_ptr->sortedBegin(quadtree_ptr->getArea()); i != quadtree_ptr->sortedEnd(); ++i)
		size++;
	TEST_EQUAL(size, 10000)
	
	size = 0;
	for (Tree::SortedIterator i = quadtree_ptr->sortedBegin(area); i != quadtree_ptr->sortedEnd(); ++i)
		size++;
	TEST_EQUAL(size, inarea)
RESULT

// CHECK(SortedIterator)
// 	std::vector<int> vi;
// 	int num = 0;
// 	for (Tree::SortedIterator i = quadtree_ptr->sortedBegin(quadtree_ptr->getArea()); i != quadtree_ptr->sortedEnd(); ++i, ++num)
// 		vi.push_back(*(i->second));
// 	
// 	TEST_EQUAL(num, 10000)
// 	TEST_EQUAL(is_sorted(vi.begin(), vi.end()), true);
// 	
// 	for (int i = 0; i != 100; i++)
// 	  std::cout << vi[i] << std::endl;
// RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
