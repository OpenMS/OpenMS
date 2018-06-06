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

#include <OpenMS/DATASTRUCTURES/Date.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Date, "$Id$")

/////////////////////////////////////////////////////////////

Date* s_ptr = nullptr;
Date* s_nullPointer= nullptr;

START_SECTION((Date()))
	s_ptr = new Date();
  TEST_NOT_EQUAL(s_ptr, s_nullPointer)
END_SECTION

START_SECTION(([EXTRA]~Date()))
	delete s_ptr;
END_SECTION

START_SECTION(Date(const QDate &date))
	QDate qd(1999,12,24);
	Date d(qd);
	TEST_EQUAL(d.year(),1999)
	TEST_EQUAL(d.month(),12)
	TEST_EQUAL(d.day(),24)
END_SECTION

START_SECTION((void get(UInt& month, UInt& day, UInt& year) const))
  Date date;
  UInt d,m,y;
  date.set("2007-12-03");
  date.get(m,d,y);
  TEST_EQUAL(m,12);
  TEST_EQUAL(d,3);
  TEST_EQUAL(y,2007);
END_SECTION

START_SECTION((void set(UInt month, UInt day, UInt year) ))
  Date date;
  UInt d,m,y;
  date.set(12,1,1977);
  date.get(m,d,y);
  TEST_EQUAL(m,12);
  TEST_EQUAL(d,1);
  TEST_EQUAL(y,1977);
  
  //exceptions
  TEST_EXCEPTION(Exception::ParseError,date.set(0,12,1977));
  TEST_EXCEPTION(Exception::ParseError,date.set(12,0,1977));
  TEST_EXCEPTION(Exception::ParseError,date.set(1,32,1977));
  TEST_EXCEPTION(Exception::ParseError,date.set(13,1,1977));
	TEST_EXCEPTION(Exception::ParseError,date.set(02,29,2100));
END_SECTION

START_SECTION((Date& operator= (const Date& source)))
  Date date, date2;
  date.set(12,1,1977);
  TEST_EQUAL(date==date2,false);
	date2 = date;
	TEST_EQUAL(date==date2,true);
END_SECTION

START_SECTION((Date(const Date& date)))
  Date date;
  date.set(12,1,1977);
	Date date2(date);
	TEST_EQUAL(date==date2,true);
END_SECTION

START_SECTION((void set(const String& date) ))
  Date date;
  //german
  date.set("01.12.1977");
  UInt d,m,y;
  date.get(m,d,y);
  TEST_EQUAL(m,12);
  TEST_EQUAL(d,1);
  TEST_EQUAL(y,1977);  

  //english
  date.set("12/01/1977");
  date.get(m,d,y);
  TEST_EQUAL(m,12);
  TEST_EQUAL(d,1);
  TEST_EQUAL(y,1977);

  //iso/ansi
  date.set("1967-12-23");
  date.get(m,d,y);
  TEST_EQUAL(d,23);
  TEST_EQUAL(m,12);
  TEST_EQUAL(y,1967);
    
   //german short
  date.set("06.01.1688");
  date.get(m,d,y);
  TEST_EQUAL(m,1);
  TEST_EQUAL(d,6);
  TEST_EQUAL(y,1688);

	//exceptions
  TEST_EXCEPTION(Exception::ParseError,date.set("bla"));
  TEST_EXCEPTION(Exception::ParseError,date.set("01.01.01.2005"));
  TEST_EXCEPTION(Exception::ParseError,date.set("f1.01.1977"));
  TEST_EXCEPTION(Exception::ParseError,date.set("01.1x.1977"));
  TEST_EXCEPTION(Exception::ParseError,date.set("01.12.i135"));
  TEST_EXCEPTION(Exception::ParseError,date.set("1135-64-3"));
END_SECTION

START_SECTION((String get() const))
  Date d;
  TEST_EQUAL(d.get(),"0000-00-00");
  d.set("11.12.1977");
  TEST_EQUAL(d.get(),"1977-12-11");
	d.set("02.01.1999");
  TEST_EQUAL(d.get(),"1999-01-02");
END_SECTION

START_SECTION((void clear()))
  Date d;
  d.set("11.12.1977");
  TEST_EQUAL(d.get(),"1977-12-11");
	d.clear();
  TEST_EQUAL(d.get(),"0000-00-00");
END_SECTION

START_SECTION((static Date today()))
  TEST_EQUAL(Date::today().isValid(), true)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
