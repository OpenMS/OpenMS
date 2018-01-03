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
// $Maintainer: Timo Sachsenberg,Andreas Bertsch$
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <iostream>
///////////////////////////
#include <OpenMS/DATASTRUCTURES/BigString.h>
//////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(BigString, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
typedef std::pair <String, String> FASTAEntry;

BigString* ptr = nullptr;
BigString* nullPointer = nullptr;

START_SECTION(BigString())
	ptr = new BigString();
	TEST_EQUAL (ptr->getBigString(),"$");
	TEST_EQUAL (ptr->size(),1);
	TEST_EQUAL (ptr->length(),1);
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~BigString())
	delete ptr;
END_SECTION

START_SECTION(void add(FASTAEntry const &new_entry))
	ptr = new BigString();
	const FASTAEntry fe ("ENTRY 1","AAAAA");
	ptr->add(fe);
	TEST_EQUAL (ptr->getBigString(),"$AAAAA$");
	const FASTAEntry fe2 ("ENTRY 2","BBBBB");
	ptr->add(fe2);
	TEST_EQUAL (ptr->getBigString(),"$AAAAA$BBBBB$");
	TEST_EQUAL (ptr->size(),3);
	TEST_EQUAL (ptr->length(),13);
END_SECTION

START_SECTION(void setSeparator(const char sep))
	ptr = new BigString();
	ptr->setSeparator('&');
	TEST_EQUAL(ptr->getSeparator(),'&');
END_SECTION

START_SECTION(char getSeparator())
	ptr = new BigString();
	TEST_EQUAL(ptr->getSeparator(),'$');
	ptr->setSeparator('&');
	TEST_EQUAL(ptr->getSeparator(),'&');
END_SECTION

START_SECTION(BigString(const BigString &bs))
	ptr = new BigString();
	const FASTAEntry fe ("ENTRY 1","AAAAA");
	ptr->add(fe);
	TEST_EQUAL (ptr->getBigString(),"$AAAAA$");
	const FASTAEntry fe2 ("ENTRY 2","BBBBB");
	ptr->add(fe2);
	ptr->setSeparator('&');
	BigString * new_big_string = new BigString (*ptr);
	TEST_EQUAL (ptr->getSeparator(),new_big_string->getSeparator());
	TEST_EQUAL (ptr->getBigString(),new_big_string->getBigString());
	TEST_EQUAL (ptr->size(),new_big_string->size());
	TEST_EQUAL (ptr->length(),new_big_string->length());	
	pair<String, String> result, ptr_result;
	new_big_string->getPeptide(result, 2, 2);
	ptr->getPeptide(ptr_result, 2, 2);
	TEST_EQUAL (ptr_result.first, result.first);
	TEST_EQUAL (ptr_result.second, result.second);
END_SECTION

START_SECTION(const String& getBigString() const )
	ptr = new BigString();
	const FASTAEntry fe ("ENTRY 1","AAAAA");
	ptr->add(fe);
	TEST_EQUAL (ptr->getBigString(),"$AAAAA$");
	const FASTAEntry fe2 ("ENTRY 2","BBBBB");
	ptr->add(fe2);
	TEST_EQUAL (ptr->getBigString(),"$AAAAA$BBBBB$");
END_SECTION

START_SECTION(Size size())
	ptr = new BigString();
	TEST_EQUAL (ptr->size(),1);
	const FASTAEntry fe ("ENTRY 1","AAAAA");
	for (Size i= 1; i < 10; i++)
	{
		ptr->add(fe);
		TEST_EQUAL (ptr->size(),i+1);
	}
END_SECTION

START_SECTION(Size length())
	ptr = new BigString();
	TEST_EQUAL (ptr->length(),1);
	const FASTAEntry fe ("ENTRY 1","AAAAA");
	for (Size i= 1; i < 10; i++)
	{
		ptr->add(fe);
		TEST_EQUAL (ptr->length(), i * 6 + 1);
	}
END_SECTION

START_SECTION(void getPeptide(FASTAEntry& entry, Size start, Size length))
	ptr = new BigString();
	const FASTAEntry fe ("ENTRY 1","AAAAA");
	ptr->add(fe);
	const FASTAEntry fe2 ("ENTRY 2","BBBBB");
	ptr->add(fe2);
	const FASTAEntry fe3 ("ENTRY 3","CCCCC");
	ptr->add(fe3);
	const FASTAEntry fe4 ("ENTRY 4","DDDDD");
	ptr->add(fe4);
	const FASTAEntry fe5 ("ENTRY 5","EEEEE");
	ptr->add(fe5);
	FASTAEntry res;
	ptr->getPeptide(res, 1, 3);
	TEST_EQUAL(res.first,"ENTRY 1");
	TEST_EQUAL(res.second,"AAA");
	ptr->getPeptide(res, 1, 5);
	TEST_EQUAL(res.first,"ENTRY 1");
	TEST_EQUAL(res.second,"AAAAA");
	ptr->getPeptide(res, 3, 2);
	TEST_EQUAL(res.first,"ENTRY 1");
	TEST_EQUAL(res.second,"AA");
	ptr->getPeptide(res, 7, 2);
	TEST_EQUAL(res.first,"ENTRY 2");
	TEST_EQUAL(res.second,"BB");
	ptr->getPeptide(res, 19, 2);
	TEST_EQUAL(res.first,"ENTRY 4");
	TEST_EQUAL(res.second,"DD");
	TEST_EXCEPTION(Exception::InvalidValue, ptr->getPeptide(res, 1,10));
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
