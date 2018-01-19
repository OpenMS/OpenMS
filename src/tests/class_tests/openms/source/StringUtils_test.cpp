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
// $Maintainer: Timo Sachsenberg, Chris Bielow $
// $Authors: Marc Sturm, Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(StringUtils, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

StringUtils* ptr = nullptr;
StringUtils* null_ptr = nullptr;
START_SECTION(StringUtils())
{
	ptr = new StringUtils();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~StringUtils())
{
	delete ptr;
}
END_SECTION

START_SECTION((static String numberLength(double d, UInt n)))
{
  // TODO
}
END_SECTION

START_SECTION((static String number(double d, UInt n)))
{
  // TODO
}
END_SECTION

START_SECTION((static String& fillLeft(String &this_s, char c, UInt size)))
{
  // TODO
}
END_SECTION

START_SECTION((static String& fillRight(String &this_s, char c, UInt size)))
{
  // TODO
}
END_SECTION

START_SECTION((static bool hasPrefix(const String &this_s, const String &string)))
{
  // TODO
}
END_SECTION

START_SECTION((static bool hasSuffix(const String &this_s, const String &string)))
{
  // TODO
}
END_SECTION

START_SECTION((static bool hasSubstring(const String &this_s, const String &string)))
{
  // TODO
}
END_SECTION

START_SECTION((static bool has(const String &this_s, Byte byte)))
{
  // TODO
}
END_SECTION

START_SECTION((static String prefix(const String &this_s, size_t length)))
{
  // TODO
}
END_SECTION

START_SECTION((static String suffix(const String &this_s, size_t length)))
{
  // TODO
}
END_SECTION

START_SECTION((static String prefix(const String &this_s, Int length)))
{
  // TODO
}
END_SECTION

START_SECTION((static String suffix(const String &this_s, Int length)))
{
  // TODO
}
END_SECTION

START_SECTION((static String prefix(const String &this_s, char delim)))
{
  // TODO
}
END_SECTION

START_SECTION((static String suffix(const String &this_s, char delim)))
{
  // TODO
}
END_SECTION

START_SECTION((static String substr(const String &this_s, size_t pos, size_t n)))
{
  // TODO
}
END_SECTION

START_SECTION((static String chop(const String &this_s, Size n)))
{
  // TODO
}
END_SECTION

START_SECTION((static String& trim(String &this_s)))
{
  // TODO
}
END_SECTION

START_SECTION((static String& quote(String &this_s, char q, String::QuotingMethod method)))
{
  // TODO
}
END_SECTION

START_SECTION((static String& unquote(String &this_s, char q, String::QuotingMethod method)))
{
  // TODO
}
END_SECTION

START_SECTION((static String& simplify(String &this_s)))
{
  // TODO
}
END_SECTION

START_SECTION((static String random(UInt length)))
{
  // TODO
}
END_SECTION

START_SECTION((static String& reverse(String &this_s)))
{
  // TODO
}
END_SECTION

START_SECTION((static bool split(const String &this_s, const char splitter, std::vector< String > &substrings, bool quote_protect)))
{
  // TODO
}
END_SECTION

START_SECTION((static bool split(const String &this_s, const String &splitter, std::vector< String > &substrings)))
{
  // TODO
}
END_SECTION

START_SECTION((static bool split_quoted(const String &this_s, const String &splitter, std::vector< String > &substrings, char q, String::QuotingMethod method)))
{
  // TODO
}
END_SECTION

START_SECTION((static QString toQString(const String &this_s)))
{
  // TODO
}
END_SECTION

START_SECTION((static Int toInt(const String &this_s)))
{
  // easy case
  TEST_EQUAL(StringUtils::toInt("1234"), 1234)
  // with spaces (allowed)
  TEST_EQUAL(StringUtils::toInt("  1234"), 1234)
  TEST_EQUAL(StringUtils::toInt("1234 "), 1234)
  TEST_EQUAL(StringUtils::toInt("   1234  "), 1234)
  // with trailing chars (unexplained) --> error (because it means the input was not split correctly beforehand)!!!
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt("1234  moreText"))  // 'moreText' is not explained...
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt(" 1234 911.0"))     // '911.0' is not explained...
  // incorrect type
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt(" abc "))
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toInt(" 123.45 "))
}
END_SECTION

START_SECTION((static float toFloat(const String &this_s)))
{
  // easy case
  TEST_REAL_SIMILAR(StringUtils::toFloat("1234.45"), 1234.45)
  // with spaces (allowed)
  TEST_REAL_SIMILAR(StringUtils::toFloat("  1234.45"), 1234.45)
  TEST_REAL_SIMILAR(StringUtils::toFloat("1234.45 "), 1234.45)
  TEST_REAL_SIMILAR(StringUtils::toFloat("   1234.45  "), 1234.45)
  // with trailing chars (unexplained) --> error (because it means the input was not split correctly beforehand)!!!
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toFloat("1234.45  moreText"))  // 'moreText' is not explained...
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toFloat(" 1234.45 911.0"))     // '911.0' is not explained...
  // incorrect type
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toFloat(" abc "))
}
END_SECTION

START_SECTION((static double toDouble(const String &this_s)))
{
  // easy case
  TEST_REAL_SIMILAR(StringUtils::toDouble("1234.45"), 1234.45)
  // with spaces (allowed)
  TEST_REAL_SIMILAR(StringUtils::toDouble("  1234.45"), 1234.45)
  TEST_REAL_SIMILAR(StringUtils::toDouble("1234.45 "), 1234.45)
  TEST_REAL_SIMILAR(StringUtils::toDouble("   1234.45  "), 1234.45)
  // with trailing chars (unexplained) --> error (because it means the input was not split correctly beforehand)!!!
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toDouble("1234.45  moreText"))  // 'moreText' is not explained...
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toDouble(" 1234.45 911.0"))     // '911.0' is not explained...
  // incorrect type
  TEST_EXCEPTION(Exception::ConversionError, StringUtils::toDouble(" abc "))
}
END_SECTION

START_SECTION((static String& toUpper(String &this_s)))
{
  // TODO
}
END_SECTION

START_SECTION((static String& firstToUpper(String &this_s)))
{
  // TODO
}
END_SECTION

START_SECTION((static String& toLower(String &this_s)))
{
  // TODO
}
END_SECTION

START_SECTION((static String& substitute(String &this_s, char from, char to)))
{
  // TODO
}
END_SECTION

START_SECTION((static String& substitute(String &this_s, const String &from, const String &to)))
{
  // TODO
}
END_SECTION

START_SECTION((static String& remove(String &this_s, char what)))
{
  // TODO
}
END_SECTION

START_SECTION((static String& ensureLastChar(String &this_s, char end)))
{
  // TODO
}
END_SECTION

START_SECTION((static String& removeWhitespaces(String &this_s)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



