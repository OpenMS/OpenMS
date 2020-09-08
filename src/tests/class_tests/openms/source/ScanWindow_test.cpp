// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/ScanWindow.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

// static_assert(OpenMS::Test::fulfills_rule_of_5<ScanWindow>(), "Must fulfill rule of 5");
// static_assert(OpenMS::Test::fulfills_rule_of_6<ScanWindow>(), "Must fulfill rule of 6");
// static_assert(OpenMS::Test::fulfills_fast_vector<ScanWindow>(), "Must have fast vector semantics");
// static_assert(std::is_nothrow_move_constructible<ScanWindow>::value, "Must have nothrow move constructible");

START_TEST(ScanWindow, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ScanWindow* ptr = nullptr;
ScanWindow* nullPointer = nullptr;
START_SECTION((ScanWindow()))
	ptr = new ScanWindow();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~ScanWindow()))
	delete ptr;
END_SECTION

START_SECTION((ScanWindow(const ScanWindow& source)))
  ScanWindow tmp;
  tmp.begin = 1.0;
  tmp.end = 2.0;
  tmp.setMetaValue("label",String("label"));
  
  ScanWindow tmp2(tmp);
  TEST_REAL_SIMILAR(tmp2.begin, 1.0)
  TEST_REAL_SIMILAR(tmp2.end, 2.0)
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");  
END_SECTION

START_SECTION((ScanWindow& operator= (const ScanWindow& source)))
  ScanWindow tmp;
  tmp.begin = 1.0;
  tmp.end = 2.0;
  tmp.setMetaValue("label",String("label"));
  
  ScanWindow tmp2;
  tmp2 = tmp;
  TEST_REAL_SIMILAR(tmp2.begin, 1.0)
  TEST_REAL_SIMILAR(tmp2.end, 2.0)
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");  
END_SECTION

START_SECTION((bool operator==(const ScanWindow &source) const ))
  ScanWindow edit, empty;
  
  TEST_EQUAL(edit==empty,true);
  
  edit.begin = 1.0;
  TEST_EQUAL(edit==empty,false);
  
  edit = empty; 
  edit.end = 1.0;
  TEST_EQUAL(edit==empty,false);
  
	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit==empty,false);
END_SECTION

START_SECTION((bool operator!=(const ScanWindow &source) const ))
  ScanWindow edit, empty;
  
  TEST_EQUAL(edit!=empty,false);
  
  edit.begin = 1.0;
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty; 
  edit.end = 1.0;
  TEST_EQUAL(edit!=empty,true);
  
	edit = empty;
	edit.setMetaValue("label",String("label"));
	TEST_EQUAL(edit!=empty,true);
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



