// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/Map.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Map, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Map<int, int>* map_ptr = nullptr;
Map<int, int>* map_nullPointer = nullptr;
START_SECTION((Map()))
	map_ptr = new Map<int, int>;
  TEST_NOT_EQUAL(map_ptr, map_nullPointer)
END_SECTION

START_SECTION((~Map()))
	delete map_ptr;
END_SECTION

START_SECTION((T& operator [] (const Key& key)))
	Map<int, int> hm;
	hm[0] = 0;
	hm[0] = 1;
	hm[1] = 2;
	hm[2] = 4;
	hm[3] = 8;
	hm[4] = 16;
	hm[5] = 32;
	TEST_EQUAL(hm.size(), 6)
	TEST_EQUAL(hm[0], 1)
	TEST_EQUAL(hm[1], 2)
	TEST_EQUAL(hm[2], 4)
	TEST_EQUAL(hm[3], 8)
	TEST_EQUAL(hm[4], 16)
	TEST_EQUAL(hm[5], 32)
END_SECTION

START_SECTION((const T & operator[](const Key &key) const ))
	Map<int, int> hm;
	hm[0] = 0;
	hm[0] = 1;
	hm[1] = 2;
	hm[2] = 4;
	hm[3] = 8;
	hm[4] = 16;
	hm[5] = 32;
	const Map<int, int>& const_map = const_cast<const Map<int, int>&>(hm);
	TEST_EQUAL(const_map.size(), 6)
	TEST_EQUAL(const_map[0], 1)
	TEST_EQUAL(const_map[1], 2)
	TEST_EQUAL(const_map[2], 4)
	TEST_EQUAL(const_map[3], 8)
	TEST_EQUAL(const_map[4], 16)
	TEST_EQUAL(const_map[5], 32)
	typedef Map<int,int> MyMap; // otherwise next line wont work
	TEST_EXCEPTION(MyMap::IllegalKey, const_map[6])
END_SECTION

START_SECTION((bool has(const Key& key) const))
	Map<int, int> hm;
	hm.insert(Map<int, int>::ValueType(0, 0));
	hm.insert(Map<int, int>::ValueType(1, 1));
	TEST_EQUAL(hm.has(0), true)
	TEST_EQUAL(hm.has(1), true)
	TEST_EQUAL(hm.has(2), false)
END_SECTION

START_SECTION(([Map::IllegalKey] IllegalKey(const char *file, int line, const char *function)))
	// already tested in const T & operator[](const Key &key) const
	NOT_TESTABLE
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
