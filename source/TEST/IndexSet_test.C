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
// $Maintainer: Jens Joachim $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/IndexSet.h>
#include <sstream>


///////////////////////////

START_TEST(DataValue, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

// default ctor
IndexSet* ptr = 0;
CHECK(IndexSet())
	ptr = new IndexSet;
	TEST_NOT_EQUAL(ptr, 0)
	
	IndexSet is1;	
	TEST_EQUAL(is1.isEmpty(), 1)
RESULT

// destructor
CHECK(~IndexSet())
	delete ptr;
RESULT


// assignment operator
CHECK(IndexSet& operator = (const IndexSet& source))
	IndexSet is1(10,20);	
  is1.add(15).add(21);
  is1.sort();
  IndexSet is2(10,21);
  IndexSet is3;
  is3 = is1;
  is1.clear();
	TEST_EQUAL(is3, is2)
RESULT

// construct set with index
CHECK(IndexSet(Size index))
	IndexSet is2(4);	
	TEST_EQUAL(is2.isEmpty(), 0)
RESULT

// construct set [<i>from_index</i>,...,<i>to_index</to>]
CHECK(IndexSet(Size index_from, Size index_to))
	IndexSet s(1,10);
  IndexSet t;
  t.add(10,1);
  TEST_EQUAL(s,t)
RESULT

CHECK(bool operator == (const IndexSet& rhs) const)
	IndexSet s(1,10);
  IndexSet t;
  t.add(2,10).add(1);
  t.sort();
  TEST_EQUAL(s, t)
RESULT

CHECK(bool operator != (const IndexSet& rhs) const)
	IndexSet s(1,10);
  IndexSet t(5);
  TEST_NOT_EQUAL(s, t)
RESULT

CHECK(IndexSet& add(Size index))
	IndexSet s(1,6);
  IndexSet t(4);
  t.add(5).add(4).add(3).add(1).add(2).add(6);
  t.sort();
  TEST_EQUAL(s, t)
RESULT
	
	
// append indices [<i>from_index</i>,...,<i>to_index</to>] to set 
CHECK(IndexSet& add(Size index_from, Size index_to))
	IndexSet s(1,10);
  s.add(12,15).add(25,20).add(29);

  IndexSet t(10,20);
  t.add(15,25).add(26,26).add(11,13).add(12);
  t.sort();
  TEST_NOT_EQUAL(s, t)

	s = IndexSet(10,26);
  TEST_EQUAL(s, t) 
RESULT


// remove index from set
CHECK(IndexSet& remove(Size index))
	IndexSet s(9,19);
  s.add(21,25);

  IndexSet t(8,26);
  t.remove(20).remove(26).remove(8).remove(100);
  TEST_EQUAL(s, t) 
RESULT
 
// remove indices [<i>from_index</i>,...,<i>to_index</to>] from set 
CHECK(IndexSet& remove(Size index_from, Size index_to))
	IndexSet s(9,17);
  s.add(21,25);

  IndexSet t(1,30);
  t.remove(30,26).remove(18,20).remove(0,8).remove(100,50);
  TEST_EQUAL(s, t)
RESULT

 
CHECK(void clear())
  IndexSet t(54);
	TEST_NOT_EQUAL(t.isEmpty(), 1)
	t.clear();
  bool res = t.isEmpty();
	TEST_EQUAL(res, 1)  
  t.add(51,100);
  res = t.isEmpty();
	TEST_NOT_EQUAL(res, 1)
	t.clear();
  res = t.isEmpty();
	TEST_EQUAL(res, 1)  
RESULT

CHECK(ostream& operator << (ostream& os, const IndexSet& set);)
  IndexSet t(12,15);
  t.add(105,103).add(100);
  t.sort();
  
  stringstream stream;
  stream << t;
	TEST_EQUAL(stream.str(), "12..15,100,103..105")
RESULT

CHECK( ConstIterator begin() const)
	IndexSet t(12,13);
  t.add(105,103).add(100);
  t.sort();
  IndexSet::const_iterator it = t.begin();
  Size i = *it;
  TEST_EQUAL(i, 12)
	i = *++it;
  TEST_EQUAL(i, 13)
	i = *++it;
  TEST_EQUAL(i, 100)
	i = *++it;
  TEST_EQUAL(i, 103)		
	i = *++it;
  TEST_EQUAL(i, 104)
	i = *++it;
  TEST_EQUAL(i, 105)	
RESULT

CHECK( ConstIterator begin(Size index) const)
	IndexSet t(12,15);
  t.add(105,103).add(100);
  t.sort();
  IndexSet::const_iterator it = t.begin(102);
  Size i = *it;
  TEST_EQUAL(i, 103)		
	i = *++it;
  TEST_EQUAL(i, 104)
	i = *++it;
  TEST_EQUAL(i, 105)	
RESULT


CHECK( ConstIterator end() const)
	IndexSet t(12,12);
  IndexSet::const_iterator it = ++t.begin();
  bool res = it==t.end();
	TEST_EQUAL(res, 1)
RESULT
		


//==================================================================================

typedef IndexSet::IndexSetConstIterator Iterator;
Iterator* it_ptr = 0;

CHECK(IndexSetConstIterator())
	it_ptr = new Iterator;
	TEST_NOT_EQUAL(it_ptr, 0)
RESULT

CHECK( ~IndexSetConstIterator() )
	delete it_ptr;
RESULT

CHECK( IndexSetConstIterator(const Size index, const VecIndex pos, const IndexSet* ref) )
	IndexSet s(5,10);
  Iterator it(5,0,&s);
  Size i = *it;
  TEST_EQUAL( i, 5)

  it = Iterator(8,3,&s);
  i = *it;
  TEST_EQUAL( i, 8)
RESULT


CHECK(IndexSetConstIterator(const IndexSetConstIterator& it) )
	IndexSet s(5,10);
  Iterator it1(5,0,&s);
  Iterator it2(it1);

  Size i = *it2;
  TEST_EQUAL( i, 5)
	i = *++it2;
  TEST_EQUAL( i, 6)
RESULT

     
CHECK(IndexSetConstIterator& operator = (const IndexSetConstIterator& rhs) )
	IndexSet s(5,10);
  Iterator it1(5,0,&s);
  Iterator it2;
  it2 = it1;
  it1 = Iterator();
  Size i = *it2;
  TEST_EQUAL( i, 5)
	i = *++it2;
  TEST_EQUAL( i, 6)
RESULT
        
CHECK(bool operator < (const IndexSetConstIterator& it) const )
	IndexSet s(5,10);
  Iterator it1(5,0,&s);

  Iterator it2(8,3,&s);
  bool res = it1<it2;
  TEST_EQUAL(res , 1)
	res = it2<it1;
  TEST_EQUAL(res , 0)
	res = it1<it1;
  TEST_EQUAL(res , 0)
RESULT
      
CHECK(bool operator > (const IndexSetConstIterator& it) const )
	IndexSet s(5,10);
  Iterator it1(5,0,&s);

  Iterator it2(8,3,&s);
  bool res = it1>it2;
  TEST_EQUAL(res , 0)
	res = it2>it1;
  TEST_EQUAL(res , 1)
	res = it1>it1;
  TEST_EQUAL(res , 0)
RESULT

CHECK(bool operator <= (const IndexSetConstIterator& it) const )
	IndexSet s(5,10);
  Iterator it1(5,0,&s);

  Iterator it2(8,3,&s);
  bool res = it1<=it2;
  TEST_EQUAL(res , 1)
	res = it2<=it1;
  TEST_EQUAL(res , 0)
	res = it1<=it1;
  TEST_EQUAL(res , 1)
RESULT
      
CHECK(bool operator >= (const IndexSetConstIterator& it) const )
	IndexSet s(5,10);
  Iterator it1(5,0,&s);

  Iterator it2(8,3,&s);
  bool res = it1>=it2;
  TEST_EQUAL(res , 0)
	res = it2>=it1;
  TEST_EQUAL(res , 1)
	res = it1>=it1;
  TEST_EQUAL(res , 1)
RESULT

CHECK(bool operator == (const IndexSetConstIterator& it) const )
	IndexSet s(5,10);
  Iterator it1(5,0,&s);

  Iterator it2(6,1,&s);
  bool res = it1==it2;
  TEST_EQUAL( res, 0)
	res = it1==--it2;
  TEST_EQUAL( res, 1)
RESULT


CHECK(bool operator != (const IndexSetConstIterator& it) const )
	IndexSet s(5,10);
  Iterator it1(5,0,&s);

  Iterator it2(6,1,&s);
  bool res = it1!=it2;
  TEST_EQUAL( res, 1)
	res = it1!=--it2;
  TEST_EQUAL( res, 0)
RESULT

CHECK( const Size& operator * ())
	IndexSet s(5,10);
  Iterator it(5,0,&s);

  Size i = *it;
  TEST_EQUAL( i, 5)
RESULT

CHECK(IndexSetConstIterator& operator ++ ())
	IndexSet s(5,10);
  Iterator it(8,1,&s);

  Size i = *++it;
  TEST_EQUAL( i, 9)
	i = *++it;
  TEST_EQUAL( i, 10)

	bool res = ++it==s.end();
  TEST_EQUAL( res, 1)
RESULT

CHECK(IndexSetConstIterator operator ++ (int))
	IndexSet s(5,10);
  Iterator it(8,1,&s);

  Size i = *it++;
  TEST_EQUAL( i , 8)
	i = *it++;
  TEST_EQUAL( i , 9)
	it++;

  TEST_EQUAL( it==s.end(), true)
RESULT

 
CHECK(IndexSetConstIterator& operator --())
	IndexSet s(5,10);
  Iterator it = s.begin();
  it++;
  Size i = *--it;
  TEST_EQUAL(i, 5)

	it = s.end();
	i = *--it;
  TEST_EQUAL(i, 10)
RESULT


CHECK(IndexSetConstIterator operator -- (int))
	IndexSet s(5,10);
  Iterator it = s.end();

  bool res = ( (it--) == s.end());
  TEST_EQUAL( res, 1)
	Size i = *it--;
  TEST_EQUAL( i, 10)
RESULT
	 
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
