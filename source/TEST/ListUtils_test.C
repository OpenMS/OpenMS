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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;
using namespace std;

START_TEST(ListUtils, "$Id$")

START_SECTION(static std::vector<T> create(const String& s))
{
  // Int
  std::vector<Int> iv;
  iv.push_back(1);
  iv.push_back(2);
  iv.push_back(3);
  iv.push_back(4);

	TEST_EQUAL(ListUtils::contains(iv, 1),true)
	TEST_EQUAL(ListUtils::contains(iv, 2),true)
	TEST_EQUAL(ListUtils::contains(iv, 3),true)
	TEST_EQUAL(ListUtils::contains(iv, 4),true)
	TEST_EQUAL(ListUtils::contains(iv, 5),false)
	TEST_EQUAL(ListUtils::contains(iv, 1011),false)

  //
  std::vector<DoubleReal> dv;
  dv.push_back(1.2);
  dv.push_back(3.4);
  TEST_EQUAL(ListUtils::contains(dv, 1.2),true)
  TEST_EQUAL(ListUtils::contains(dv, 1.21),false)
  TEST_EQUAL(ListUtils::contains(dv, 1.19),false)
  TEST_EQUAL(ListUtils::contains(dv, 1.21,0.02),true)
  TEST_EQUAL(ListUtils::contains(dv, 1.19,0.02),true)
  TEST_EQUAL(ListUtils::contains(dv, 3.4),true)
  TEST_EQUAL(ListUtils::contains(dv, 4.2),false)
  TEST_EQUAL(ListUtils::contains(dv, 2.0),false)
  TEST_EQUAL(ListUtils::contains(dv, 0.0),false)

  // String
  std::vector<String> sv;
  sv.push_back("yes");
  sv.push_back("no");
	TEST_EQUAL(ListUtils::contains(sv, "yes"),true)
	TEST_EQUAL(ListUtils::contains(sv, "no"),true)
	TEST_EQUAL(ListUtils::contains(sv, "jup"),false)
	TEST_EQUAL(ListUtils::contains(sv, ""),false)
	TEST_EQUAL(ListUtils::contains(sv, "noe"),false)
}
END_SECTION

START_SECTION((static std::vector<T> create(const String& s)))
{
  std::vector<String> sv = ListUtils::create<String>("yes,no");
  TEST_EQUAL(sv.size(), 2)
  ABORT_IF(sv.size() != 2)
	TEST_EQUAL(sv[0], "yes")
	TEST_EQUAL(sv[1], "no")

  std::vector<DoubleReal> dv = ListUtils::create<DoubleReal>("1.2,3.5");
  TEST_EQUAL(dv.size(), 2)
  ABORT_IF(dv.size() != 2)
	TEST_EQUAL(dv[0], 1.2)
	TEST_EQUAL(dv[1], 3.5)

	std::vector<Int> iv = ListUtils::create<Int>("1,5");
	TEST_EQUAL(iv.size(),2);
  ABORT_IF(iv.size() != 2)
	TEST_EQUAL(iv[0], 1)
	TEST_EQUAL(iv[1], 5)

	IntList iv2 = ListUtils::create<Int>("2");
	TEST_EQUAL(iv2.size(),1);
	TEST_EQUAL(iv2[0],2);

	IntList iv3 = ListUtils::create<Int>("");
	TEST_EQUAL(iv3.size(),0);

}
END_SECTION

START_SECTION((static String concatenate(const std::vector<ContainerType>& container, const String & glue = "")))
{
	std::vector<String> list;
  list.push_back("1");
  list.push_back("2");
  list.push_back("3");
  list.push_back("4");
  list.push_back("5");
	TEST_STRING_EQUAL(ListUtils::concatenate(list, "g"),"1g2g3g4g5");
	TEST_STRING_EQUAL(ListUtils::concatenate(list, ""),"12345");

	list.clear();
	TEST_STRING_EQUAL(ListUtils::concatenate(list, "g"),"");
	TEST_STRING_EQUAL(ListUtils::concatenate(list, ""),"");

	//test2 (from StringList)
	std::vector<String> tmp;
	TEST_EQUAL(ListUtils::concatenate(tmp),"")
	tmp.push_back("1\n");
	tmp.push_back("2\n");
	tmp.push_back("3\n");
	TEST_EQUAL(ListUtils::concatenate(tmp),"1\n2\n3\n")
}
END_SECTION

END_TEST
