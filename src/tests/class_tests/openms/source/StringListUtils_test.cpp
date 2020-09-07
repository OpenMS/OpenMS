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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <QStringList>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(StringList, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((static StringList fromQStringList(const QStringList &rhs)))
{
  QStringList q_str_list;
  q_str_list << "First Element" << "Second Element" << "Third Element";

  StringList str_list = StringListUtils::fromQStringList(q_str_list);
  TEST_EQUAL((int)str_list.size(), q_str_list.size())
  ABORT_IF((int)str_list.size() != q_str_list.size())
  for(Size i = 0 ; i < str_list.size() ; ++i)
  {
    TEST_EQUAL(str_list[i], String(q_str_list[(int)i]))
  }
}
END_SECTION


START_SECTION((static void toUpper(StringList &sl)))
{
	StringList list = ListUtils::create<String>("yes,no");
	StringListUtils::toUpper(list);
	TEST_EQUAL(list[0],"YES")
	TEST_EQUAL(list[1],"NO")
}
END_SECTION

START_SECTION((static void toLower(StringList &sl)))
{
	StringList list = ListUtils::create<String>("yES,nO");
	StringListUtils::toLower(list);
	TEST_EQUAL(list[0],"yes")
	TEST_EQUAL(list[1],"no")
}
END_SECTION

StringList tmp_list;
tmp_list.push_back("first_line");
tmp_list.push_back("");
tmp_list.push_back("");
tmp_list.push_back("middle_line");
tmp_list.push_back("");
tmp_list.push_back("  space_line");
tmp_list.push_back("	tab_line");
tmp_list.push_back("back_space_line   ");
tmp_list.push_back("back_tab_line			");
tmp_list.push_back("");
tmp_list.push_back("last_line");

StringList tmp_list2;
tmp_list2.push_back("first_line");
tmp_list2.push_back("");
tmp_list2.push_back("");
tmp_list2.push_back("middle_line");
tmp_list2.push_back("");
tmp_list2.push_back("space_line");
tmp_list2.push_back("tab_line");
tmp_list2.push_back("back_space_line");
tmp_list2.push_back("back_tab_line");
tmp_list2.push_back("");
tmp_list2.push_back("last_line");

START_SECTION((static Iterator searchPrefix(const Iterator &start, const Iterator &end, const String &text, bool trim=false)))
{
	StringList list(tmp_list);

	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "first_line") == list.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "middle_line") == (list.begin()+3), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "space_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "tab_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "last_line") == (list.end()-1), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "invented_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin()+1, list.end(), "first_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), " ") == (list.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "\t") == (list.begin()+6), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin()+9, list.end() ,"\t") == (list.end()), true)

	//trim
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "first_line",true) == list.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "space_line",true) == (list.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "tab_line",true) == (list.begin()+6), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "invented_line",true) == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin()+1, list.end(), "first_line",true) == list.end(), true)

	//Try it on the same file (but trimmed)
	list = tmp_list2;

	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "first_line") == list.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "middle_line") == (list.begin()+3), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "space_line",true) == (list.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "tab_line",true) == (list.begin()+6), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "last_line") == (list.end()-1), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "invented_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin()+1, list.end(), "first_line") == list.end(), true)

	//trim
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "first_line",true) == list.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "space_line",true) == (list.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "tab_line",true) == (list.begin()+6), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "invented_line",true) == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin()+1, list.end(), "first_line",true) == list.end(), true)	}

END_SECTION

START_SECTION((static Iterator searchPrefix(StringList &container, const String &text, bool trim=false)))
{
	StringList list(tmp_list);

	TEST_EQUAL(StringListUtils::searchPrefix(list, "first_line") == list.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "middle_line") == (list.begin()+3), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "space_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "tab_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "last_line") == (list.end()-1), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "invented_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, " ") == (list.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "\t") == (list.begin()+6), true)

	//trim
	TEST_EQUAL(StringListUtils::searchPrefix(list, "first_line",true) == list.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "space_line",true) == (list.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "tab_line",true) == (list.begin()+6), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "invented_line",true) == list.end(), true)

	//Try it on the same file (but trimmed)
	list = tmp_list2;

	TEST_EQUAL(StringListUtils::searchPrefix(list, "first_line") == list.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "middle_line") == (list.begin()+3), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "space_line",true) == (list.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "tab_line",true) == (list.begin()+6), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "last_line") == (list.end()-1), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "invented_line") == list.end(), true)

	//trim
	TEST_EQUAL(StringListUtils::searchPrefix(list, "first_line",true) == list.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "space_line",true) == (list.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "tab_line",true) == (list.begin()+6), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "invented_line",true) == list.end(), true)
}
END_SECTION

START_SECTION((static Iterator searchSuffix(const Iterator &start, const Iterator &end, const String &text, bool trim=false)))
{
	StringList list(tmp_list);

	TEST_EQUAL(StringListUtils::searchSuffix(list.begin(), list.end(), "invented_line",true) == list.end(), true)
	TEST_EQUAL(StringListUtils::searchSuffix(list.begin(), list.end(), "back_space_line",true) == list.begin()+7, true)
	TEST_EQUAL(StringListUtils::searchSuffix(list.begin(), list.end(), "back_tab_line",true) == list.begin()+8, true)
	TEST_EQUAL(StringListUtils::searchSuffix(list.begin()+8, list.end(), "back_space_line",true) == list.end(), true)

	TEST_EQUAL(StringListUtils::searchSuffix(list.begin(), list.end(), "invented_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchSuffix(list.begin(), list.end(), "back_space_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchSuffix(list.begin(), list.end(), "back_tab_line") == list.end(), true)
}
END_SECTION

START_SECTION((static Iterator searchSuffix(StringList &container, const String &text, bool trim=false)))
{
	StringList list(tmp_list);

	TEST_EQUAL(StringListUtils::searchSuffix(list, "invented_line",true) == list.end(), true)
	TEST_EQUAL(StringListUtils::searchSuffix(list, "back_space_line",true) == list.begin()+7, true)
	TEST_EQUAL(StringListUtils::searchSuffix(list, "back_tab_line",true) == list.begin()+8, true)

	TEST_EQUAL(StringListUtils::searchSuffix(list, "invented_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchSuffix(list, "back_space_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchSuffix(list, "back_tab_line") == list.end(), true)
}
END_SECTION

START_SECTION((static ConstIterator searchPrefix(const ConstIterator &start, const ConstIterator &end, const String &text, bool trim=false)))
{
	const StringList list(tmp_list);

	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "first_line") == list.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "middle_line") == (list.begin()+3), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "space_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "tab_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "last_line") == (list.end()-1), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "invented_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin() + 1, list.end(), "first_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), " ") == (list.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "\t") == (list.begin()+6), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin() + 9, list.end(), "\t") == (list.end()), true)

	//trim
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "first_line",true) == list.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "space_line",true) == (list.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "tab_line",true) == (list.begin()+6), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin(), list.end(), "invented_line",true) == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list.begin() + 1, list.end(), "first_line",true) == list.end(), true)

	//Try it on the same file (but trimmed)
	const StringList list2 = tmp_list2;

	TEST_EQUAL(StringListUtils::searchPrefix(list2.begin(), list2.end(), "first_line") == list2.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2.begin(), list2.end(), "middle_line") == (list2.begin()+3), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2.begin(), list2.end(), "space_line",true) == (list2.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2.begin(), list2.end(), "tab_line",true) == (list2.begin()+6), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2.begin(), list2.end(), "last_line") == (list2.end()-1), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2.begin(), list2.end(), "invented_line") == list2.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2.begin() + 1, list2.end(), "first_line") == list2.end(), true)

	//trim
	TEST_EQUAL(StringListUtils::searchPrefix(list2.begin(), list2.end(), "first_line",true) == list2.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2.begin(), list2.end(), "space_line",true) == (list2.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2.begin(), list2.end(), "tab_line",true) == (list2.begin()+6), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2.begin(), list2.end(), "invented_line",true) == list2.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2.begin()  + 1, list2.end(), "first_line",true) == list2.end(), true)
}
END_SECTION

START_SECTION((static ConstIterator searchPrefix(const StringList &container, const String &text, bool trim=false)))
{
	const StringList list(tmp_list);

	TEST_EQUAL(StringListUtils::searchPrefix(list, "first_line") == list.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "middle_line") == (list.begin()+3), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "space_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "tab_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "last_line") == (list.end()-1), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "invented_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, " ") == (list.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "\t") == (list.begin()+6), true)

	//trim
	TEST_EQUAL(StringListUtils::searchPrefix(list, "first_line",true) == list.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "space_line",true) == (list.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "tab_line",true) == (list.begin()+6), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list, "invented_line",true) == list.end(), true)

	//Try it on the same file (but trimmed)
	const StringList list2 = tmp_list2;

	TEST_EQUAL(StringListUtils::searchPrefix(list2, "first_line") == list2.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2, "middle_line") == (list2.begin()+3), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2, "space_line",true) == (list2.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2, "tab_line",true) == (list2.begin()+6), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2, "last_line") == (list2.end()-1), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2, "invented_line") == list2.end(), true)

	//trim
	TEST_EQUAL(StringListUtils::searchPrefix(list2, "first_line",true) == list2.begin(), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2, "space_line",true) == (list2.begin()+5), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2, "tab_line",true) == (list2.begin()+6), true)
	TEST_EQUAL(StringListUtils::searchPrefix(list2, "invented_line",true) == list2.end(), true)
}
END_SECTION

START_SECTION((static ConstIterator searchSuffix(const ConstIterator &start, const ConstIterator &end, const String &text, bool trim=false)))
{
	const StringList list(tmp_list);

	TEST_EQUAL(StringListUtils::searchSuffix(list.begin(), list.end(), "invented_line",true) == list.end(), true)
	TEST_EQUAL(StringListUtils::searchSuffix(list.begin(), list.end(), "back_space_line",true) == list.begin()+7, true)
	TEST_EQUAL(StringListUtils::searchSuffix(list.begin(), list.end(), "back_tab_line",true) == list.begin()+8, true)
	TEST_EQUAL(StringListUtils::searchSuffix(list.begin() + 8, list.end(), "back_space_line",true) == list.end(), true)

	TEST_EQUAL(StringListUtils::searchSuffix(list.begin(), list.end(), "invented_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchSuffix(list.begin(), list.end(), "back_space_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchSuffix(list.begin(), list.end(), "back_tab_line") == list.end(), true)
}
END_SECTION

START_SECTION((static ConstIterator searchSuffix(const StringList &container, const String &text, bool trim=false)))
{
	const StringList list(tmp_list);

	TEST_EQUAL(StringListUtils::searchSuffix(list, "invented_line",true) == list.end(), true)
	TEST_EQUAL(StringListUtils::searchSuffix(list, "back_space_line",true) == list.begin()+7, true)
	TEST_EQUAL(StringListUtils::searchSuffix(list, "back_tab_line",true) == list.begin()+8, true)

	TEST_EQUAL(StringListUtils::searchSuffix(list, "invented_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchSuffix(list, "back_space_line") == list.end(), true)
	TEST_EQUAL(StringListUtils::searchSuffix(list, "back_tab_line") == list.end(), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
