// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/SVOutStream.h>
#include <sstream>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SVOutStream, "$Id$")

/////////////////////////////////////////////////////////////

SVOutStream* sv_ptr = nullptr;
SVOutStream* sv_nullPointer = nullptr;

START_SECTION((SVOutStream(std::ostream& out, const String& sep="\t", const String& replacement="_", String::QuotingMethod quoting=String::DOUBLE)))
{
  stringstream strstr;
  sv_ptr = new SVOutStream(strstr);
  TEST_NOT_EQUAL(sv_ptr, sv_nullPointer);
}
END_SECTION

START_SECTION(([EXTRA] ~SVOutStream()))
{
  delete sv_ptr;
}
END_SECTION

START_SECTION((template <typename T> SVOutStream& operator<<(const T& value)))
{
  stringstream strstr;
  SVOutStream out(strstr, ",");
  out << 123 << 3.14 << -1.23e45 << nl;
  out << 456 << endl;
  // different cases for Unix/Windows:
  TEST_EQUAL((strstr.str() == "123,3.14,-1.23e+45\n456\n") ||
             (strstr.str() == "123,3.14,-1.23e+045\n456\n") ||
             (strstr.str() == "123,3.14,-1.23e45\n456\n"), true)
  std::cout << strstr.str() << std::endl;
}
{
  stringstream strstr;
  SVOutStream out(strstr, "_/_");
  out << 123 << 3.14 << -1.23e45 << endl;
  out << 456 << nl;
  // different cases for Unix/Windows:
  TEST_EQUAL((strstr.str() == "123_/_3.14_/_-1.23e+45\n456\n") ||
                 (strstr.str() == "123_/_3.14_/_-1.23e45\n456\n") ||
             (strstr.str() == "123_/_3.14_/_-1.23e+045\n456\n"), true);
}
END_SECTION

START_SECTION((SVOutStream& operator<<(String str)))
{
  stringstream strstr;
  SVOutStream out(strstr, ",", "_", String::NONE);
  out << String("a") << string("bc") << "d,f" << nl;
  out << String("g\"i\"k") << 'l' << endl;
  TEST_EQUAL(strstr.str(), "a,bc,d_f\ng\"i\"k,l\n");
}
{
  stringstream strstr;
  SVOutStream out(strstr, ",", "_", String::ESCAPE);
  out << string("a") << "bc" << String("d,f") << nl;
  out << "g\"i\"k" << 'l' << endl;
  TEST_EQUAL(strstr.str(), "\"a\",\"bc\",\"d,f\"\n\"g\\\"i\\\"k\",\"l\"\n");
}
{
  stringstream strstr;
  SVOutStream out(strstr, ",", "_", String::DOUBLE);
  out << "a" << String("bc") << string("d,f") << nl;
  out << string("g\"i\"k") << 'l' << endl;
  TEST_EQUAL(strstr.str(), "\"a\",\"bc\",\"d,f\"\n\"g\"\"i\"\"k\",\"l\"\n");
}
{
  stringstream strstr;
  SVOutStream out(strstr, "; ", ",_", String::NONE);
  out << String("a") << "bc" << string("d; f") << nl;
  out << "g\"i\"k" << 'l' << endl;
  TEST_EQUAL(strstr.str(), "a; bc; d,_f\ng\"i\"k; l\n");
}
END_SECTION

START_SECTION((SVOutStream& operator<<(const std::string& str)))
{
  NOT_TESTABLE // tested with "operator<<(String)"
}
END_SECTION

START_SECTION((SVOutStream& operator<<(const char* c_str)))
{
  NOT_TESTABLE // tested with "operator<<(String)"
}
END_SECTION

START_SECTION((SVOutStream& operator<<(const char c)))
{
  NOT_TESTABLE // tested with "operator<<(String)"
}
END_SECTION

START_SECTION((SVOutStream& operator<<(std::ostream& (*fp)(std::ostream&))))
{
  stringstream strstr;
  SVOutStream out(strstr, ",", "_", String::ESCAPE);
  out << endl << 123 << endl << "bla";
  TEST_EQUAL(strstr.str(), "\n123\n\"bla\"");
}
END_SECTION

START_SECTION((SVOutStream& operator<<(enum Newline)))
{
  stringstream strstr;
  SVOutStream out(strstr, ",", "_", String::ESCAPE);
  out << nl << 123 << nl << "bla";
  TEST_EQUAL(strstr.str(), "\n123\n\"bla\"");
}
END_SECTION

START_SECTION((SVOutStream& write(const String& str)))
{
  stringstream strstr;
  SVOutStream out(strstr, ",", "_", String::ESCAPE);
  out << "bla" << 123 << nl;
  out.write("#This, is, a, comment\n");
  out << 4.56 << "test" << endl;
  TEST_EQUAL(strstr.str(),
             "\"bla\",123\n#This, is, a, comment\n4.56,\"test\"\n");
}
END_SECTION

START_SECTION((bool modifyStrings(bool modify)))
{
  stringstream strstr;
  SVOutStream out(strstr, ",");
  out << "test";
  bool result = out.modifyStrings(false); // "true" by default
  TEST_EQUAL(result, true);
  out << "bla";
  result = out.modifyStrings(true);
  TEST_EQUAL(result, false);
  out << "laber" << endl;
  TEST_EQUAL(strstr.str(), "\"test\",bla,\"laber\"\n");
}
END_SECTION

START_SECTION((template <typename NumericT> SVOutStream& writeValueOrNan(NumericT thing)))
{
  stringstream strstr;
  SVOutStream out(strstr, ",");
  out.writeValueOrNan(123);
  out.writeValueOrNan(3.14);
  out << nl;
  out.writeValueOrNan(456);
  out.writeValueOrNan(std::numeric_limits<double>::quiet_NaN());
  out << nl;
  out.writeValueOrNan(std::numeric_limits<double>::infinity());
  out.writeValueOrNan(-std::numeric_limits<double>::infinity());
  out << endl;
  TEST_EQUAL(strstr.str(), "123,3.14\n456,nan\ninf,-inf\n");
}
END_SECTION

END_TEST
