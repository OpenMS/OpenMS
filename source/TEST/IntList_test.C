// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: David Wojnar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/IntList.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(IntList, "$$")

/////////////////////////////////////////////////////////////

IntList* ptr = 0;
IntList* nullPointer = 0;
START_SECTION(IntList())
	ptr = new IntList;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~IntList())
	delete ptr;
END_SECTION

START_SECTION(static IntList create(const String& list))
	IntList list = IntList::create("1,5");
	TEST_EQUAL(list.size(),2);
	TEST_EQUAL(list[0],1);
	TEST_EQUAL(list[1],5);

	IntList list2 = IntList::create("2");
	TEST_EQUAL(list2.size(),1);
	TEST_EQUAL(list2[0],2);

	IntList list3 = IntList::create("");
	TEST_EQUAL(list3.size(),0);
END_SECTION

START_SECTION(static IntList create(const StringList& list))
	IntList list = IntList::create(StringList::create("1,5"));
	TEST_EQUAL(list.size(),2);
	TEST_EQUAL(list[0],1);
	TEST_EQUAL(list[1],5);
	
	IntList list2 = IntList::create(StringList::create("2"));
	TEST_EQUAL(list2.size(),1);
	TEST_EQUAL(list2[0],2);

	IntList list3 = IntList::create(StringList::create(""));
	TEST_EQUAL(list3.size(),0);
	
	TEST_EXCEPTION(Exception::ConversionError,IntList::create(StringList::create("ein,exception")));
END_SECTION


START_SECTION(IntList(const IntList& rhs))
	IntList list = IntList::create("1,3");
	IntList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_EQUAL(list2[0],1);
	TEST_EQUAL(list2[1],3);
END_SECTION

START_SECTION(IntList(const std::vector<Int>& rhs))
	std::vector<Int> list;
	list.push_back(1);
	list.push_back(3);
	IntList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_EQUAL(list2[0],1);
	TEST_EQUAL(list2[1],3);
END_SECTION

START_SECTION(IntList(const std::vector<UInt>& rhs))
	std::vector<UInt> list;
	list.push_back(1);
	list.push_back(2);
	IntList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_EQUAL(list2[0],1);
	TEST_EQUAL(list2[1],2);

END_SECTION

START_SECTION(IntList& operator=(const IntList& rhs))
	IntList list = IntList::create("1,3");
	IntList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_EQUAL(list2[0],1);
	TEST_EQUAL(list2[1],3);

END_SECTION

START_SECTION(IntList& operator=(const std::vector<Int>& rhs))
	std::vector<Int> list;
	list.push_back(1);
	list.push_back(3);
	IntList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_EQUAL(list2[0],1);
	TEST_EQUAL(list2[1],3);

END_SECTION

START_SECTION(IntList& operator=(const std::vector<UInt>& rhs))
	std::vector<UInt> list;
	list.push_back(1);
	list.push_back(3);
	IntList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_EQUAL(list2[0],1);
	TEST_EQUAL(list2[1],3);
END_SECTION

START_SECTION((template<typename IntType> IntList& operator<<(IntType value)))
	IntList list;
	list << 1 << 2 << 3 << 1;
	TEST_EQUAL(list.size(),4);
	TEST_EQUAL(list[0],1);
	TEST_EQUAL(list[1],2);
	TEST_EQUAL(list[2],3);
	TEST_EQUAL(list[3],1);
END_SECTION

START_SECTION(bool contains(Int s) const)
	IntList list = IntList::create("1,3");
	TEST_EQUAL(list.contains(1),true)
	TEST_EQUAL(list.contains(3),true)
	TEST_EQUAL(list.contains(4),false)	
	TEST_EQUAL(list.contains(2),false)
	TEST_EQUAL(list.contains(0),false)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

