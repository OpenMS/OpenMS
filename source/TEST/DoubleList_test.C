// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: David Wojnar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DoubleList.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(DoubleList, "$Id$")

/////////////////////////////////////////////////////////////

DoubleList* ptr = 0;
DoubleList* nullPointer = 0;
START_SECTION(DoubleList())
	ptr = new DoubleList;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~DoubleList())
	delete ptr;
END_SECTION

START_SECTION(static DoubleList create(const String& list))
	DoubleList list = DoubleList::create("1.222,5.33789");
	TEST_EQUAL(list.size(),2);
	TEST_REAL_SIMILAR(list[0],1.222);
	TEST_REAL_SIMILAR(list[1],5.33789);

	DoubleList list2 = DoubleList::create("2.33334");
	TEST_EQUAL(list2.size(),1);
	TEST_REAL_SIMILAR(list2[0],2.33334);

	DoubleList list3 = DoubleList::create("");
	TEST_EQUAL(list3.size(),0);
END_SECTION

START_SECTION(static DoubleList create(const StringList& list))
	DoubleList list = DoubleList::create(StringList::create("1.222,5.33789"));
	TEST_EQUAL(list.size(),2);
	TEST_REAL_SIMILAR(list[0],1.222);
	TEST_REAL_SIMILAR(list[1],5.33789);

	DoubleList list2 = DoubleList::create(StringList::create("2.33334"));
	TEST_EQUAL(list2.size(),1);
	TEST_REAL_SIMILAR(list2[0],2.33334);

	DoubleList list3 = DoubleList::create(StringList::create(""));
	TEST_EQUAL(list3.size(),0);
	TEST_EXCEPTION(Exception::ConversionError,DoubleList::create(StringList::create("ein,exception")));
END_SECTION

START_SECTION(DoubleList(const DoubleList& rhs))
	DoubleList list = DoubleList::create("1.2,3.4");
	DoubleList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_REAL_SIMILAR(list2[0],1.2);
	TEST_REAL_SIMILAR(list2[1],3.4);
END_SECTION

START_SECTION(DoubleList(const std::vector<DoubleReal>& rhs))
	std::vector<DoubleReal> list;
	list.push_back(1.2345);
	list.push_back(3.45678);
	DoubleList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_REAL_SIMILAR(list2[0],1.2345);
	TEST_REAL_SIMILAR(list2[1],3.45678);
END_SECTION

START_SECTION(DoubleList(const std::vector<Real>& rhs))
	std::vector<Real> list;
	list.push_back(1.234f);
	list.push_back(2.345f);
	DoubleList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_REAL_SIMILAR(list2[0],1.234);
	TEST_REAL_SIMILAR(list2[1],2.345);

END_SECTION

START_SECTION(DoubleList& operator=(const DoubleList& rhs))
	DoubleList list = DoubleList::create("1.22,3.33");
	DoubleList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_REAL_SIMILAR(list2[0],1.22);
	TEST_REAL_SIMILAR(list2[1],3.33);

END_SECTION

START_SECTION(DoubleList& operator=(const std::vector<DoubleReal>& rhs))
	std::vector<DoubleReal> list;
	list.push_back(1.22);
	list.push_back(3.67);
	DoubleList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_REAL_SIMILAR(list2[0],1.22);
	TEST_REAL_SIMILAR(list2[1],3.67);
END_SECTION

START_SECTION(DoubleList& operator=(const std::vector<Real>& rhs))
	std::vector<Real> list;
	list.push_back(1.22f);
	list.push_back(3.67f);
	DoubleList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_REAL_SIMILAR(list2[0],1.22);
	TEST_REAL_SIMILAR(list2[1],3.67);
END_SECTION

START_SECTION((template<typename DoubleType> DoubleList& operator<<(DoubleType value)))
	DoubleList list;
	list << 1.2 << 2.3456 << 3.5678999 << 1.2;
	TEST_EQUAL(list.size(),4);
	TEST_EQUAL(list[0],1.2);
	TEST_EQUAL(list[1],2.3456);
	TEST_EQUAL(list[2],3.5678999);
	TEST_EQUAL(list[3],1.2);
END_SECTION

START_SECTION(bool contains(DoubleReal s, DoubleReal tolerance=0.00001) const)
	DoubleList list = DoubleList::create("1.2,3.4");
	TEST_EQUAL(list.contains(1.2),true)
	TEST_EQUAL(list.contains(1.21),false)
	TEST_EQUAL(list.contains(1.19),false)
	TEST_EQUAL(list.contains(1.21,0.02),true)
	TEST_EQUAL(list.contains(1.19,0.02),true)
	TEST_EQUAL(list.contains(3.4),true)
	TEST_EQUAL(list.contains(4.2),false)
	TEST_EQUAL(list.contains(2),false)
	TEST_EQUAL(list.contains(0),false)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

