// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/PrecisionWrapper.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/DATASTRUCTURES/StringConversions.h>

#include <boost/functional/hash.hpp>

using namespace std;

namespace OpenMS
{
  const String String::EMPTY;

  String::String() :
    string()
  {
  }

  String::String(const string& s) :
    string(s)
  {
  }

  String::String(const std::string_view& sv) :
    string(sv)
  {
  }

  String::String(const char* s) :
    string(s)
  {
  }

  String::String(const QString& s) :
    string(s.toStdString())
  {
  }

  String::String(const char* s, SizeType length) :
    string(s, length)
  {
  }

  String::String(const char c) :
    string(1, c)
  {
  }

  String::String(size_t len, char c) :
    string(len, c)
  {
  }

  String::String(int i) :
    string()
  {
    StringConversions::append(i, *this);
  }

  String::String(unsigned int i) :
    string()
  {
    StringConversions::append(i, *this);
  }

  String::String(short int i) :
    string()
  {
    StringConversions::append(i, *this);
  }

  String::String(short unsigned int i) :
    string()
  {
    StringConversions::append(i, *this);
  }

  String::String(long int i) :
    string()
  {
    StringConversions::append(i, *this);
  }

  String::String(long unsigned int i) :
    string()
  {
    StringConversions::append(i, *this);
  }

  String::String(long long unsigned int i) :
    string()
  {
    StringConversions::append(i, *this);
  }

  String::String(long long signed int i) :
    string()
  {
    StringConversions::append(i, *this);
  }

  String::String(float f, bool full_precision) :
    string()
  {
    full_precision ? StringConversions::append(f, *this)
                   : StringConversions::appendLowP(f, *this);
  }

  String::String(double d, bool full_precision) :
    string()
  {
    full_precision ? StringConversions::append(d, *this)
                   : StringConversions::appendLowP(d, *this);
  }

  String::String(long double ld, bool full_precision) :
    string()
  {
    full_precision ? StringConversions::append(ld, *this)
                   : StringConversions::appendLowP(ld, *this);
  }

  String::String(const DataValue& d, bool full_precision) :
    string()
  {
    StringConversions::append(d, full_precision, *this);
  }

  String String::numberLength(double d, UInt n)
  {
    return StringUtils::numberLength(d, n);
  }

  String String::number(double d, UInt n)
  {
    return StringUtils::number(d, n);
  }

  String& String::fillLeft(char c, UInt size)
  {
    return StringUtils::fillLeft(*this, c, size);
  }

  String& String::fillRight(char c, UInt size)
  {
    return StringUtils::fillRight(*this, c, size);
  }

  bool String::hasPrefix(const String& string) const
  {
    return StringUtils::hasPrefix(*this, string);
  }

  bool String::hasSuffix(const String& string) const
  {
    return StringUtils::hasSuffix(*this, string);
  }

  bool String::hasSubstring(const String& string) const
  {
    return StringUtils::hasSubstring(*this, string);
  }

  bool String::has(Byte byte) const
  {
    return StringUtils::has(*this, byte);
  }

  String String::prefix(SizeType length) const
  {
    return StringUtils::prefix(*this, length);
  }

  String String::suffix(SizeType length) const
  {
    return StringUtils::suffix(*this, length);
  }

  String String::prefix(Int length) const
  {
    return StringUtils::prefix(*this, length);
  }

  String String::suffix(Int length) const
  {
    return StringUtils::suffix(*this, length);
  }

  String String::prefix(char delim) const
  {
    return StringUtils::prefix(*this, delim);
  }

  String String::suffix(char delim) const
  {
    return StringUtils::suffix(*this, delim);
  }

  String String::substr(size_t pos, size_t n) const
  {
    return StringUtils::substr(*this, pos, n);
  }

  String String::chop(Size n) const
  {
    return StringUtils::chop(*this, n);
  }

  String& String::trim()
  {
    return StringUtils::trim(*this);
  }

  bool String::isQuoted(char q)
  {
    return StringUtils::isQuoted(*this, q);
  }

  String& String::quote(char q, QuotingMethod method)
  {
    return StringUtils::quote(*this, q, method);
  }

  String& String::unquote(char q, QuotingMethod method)
  {
    return StringUtils::unquote(*this, q, method);
  }

  String& String::simplify()
  {
    return StringUtils::simplify(*this);
  }

  String String::random(UInt length)
  {
    return StringUtils::random(length);
  }

  String& String::reverse()
  {
    return StringUtils::reverse(*this);
  }

  bool String::split(const char splitter, vector<String>& substrings,
                     bool quote_protect) const
  {
    return StringUtils::split(*this, splitter, substrings, quote_protect);
  }

  bool String::split(const String& splitter, std::vector<String>& substrings)
  const
  {
    return StringUtils::split(*this, splitter, substrings);
  }

  bool String::split_quoted(const String& splitter, vector<String>& substrings,
                            char q, QuotingMethod method) const
  {
    return StringUtils::split_quoted(*this, splitter, substrings, q, method);
  }

  QString String::toQString() const
  {
    return StringUtils::toQString(*this);
  }

  Int String::toInt() const
  {
    if constexpr (is_same<Int, Int32>::value)
    {
      return StringUtils::toInt32(*this);
    }
    else
    {
      return StringUtils::toInt64(*this);
    }
  }

  Int32 String::toInt32() const
  {
    return StringUtils::toInt32(*this);
  }

  Int64 String::toInt64() const
  {
    return StringUtils::toInt64(*this);
  }

  float String::toFloat() const
  {
    return StringUtils::toFloat(*this);
  }

  double String::toDouble() const
  {
    return StringUtils::toDouble(*this);
  }

  String& String::toUpper()
  {
    return StringUtils::toUpper(*this);
  }

  String& String::firstToUpper()
  {
    return StringUtils::firstToUpper(*this);
  }

  String& String::toLower()
  {
    return StringUtils::toLower(*this);
  }

  String& String::substitute(char from, char to)
  {
    return StringUtils::substitute(*this, from, to);
  }

  String& String::substitute(const String& from, const String& to)
  {
    return StringUtils::substitute(*this, from, to);
  }

  String& String::remove(char what)
  {
    return StringUtils::remove(*this, what);
  }

  String& String::ensureLastChar(char end)
  {
    return StringUtils::ensureLastChar(*this, end);
  }

  String& String::removeWhitespaces()
  {
    return StringUtils::removeWhitespaces(*this);
  }

  ///
  ///// Operators
  ///

  String String::operator+(int i) const
  {
    String s(*this);
    StringConversions::append(i, s);
    return s;
  }

  String String::operator+(unsigned int i) const
  {
    String s(*this);
    StringConversions::append(i, s);
    return s;
  }

  String String::operator+(short int i) const
  {
    String s(*this);
    StringConversions::append(i, s);
    return s;
  }

  String String::operator+(short unsigned int i) const
  {
    String s(*this);
    StringConversions::append(i, s);
    return s;
  }

  String String::operator+(long int i) const
  {
    String s(*this);
    StringConversions::append(i, s);
    return s;
  }

  String String::operator+(long unsigned int i) const
  {
    String s(*this);
    StringConversions::append(i, s);
    return s;
  }

  String String::operator+(long long unsigned int i) const
  {
    String s(*this);
    StringConversions::append(i, s);
    return s;
  }

  String String::operator+(float f) const
  {
    String s(*this);
    StringConversions::append(f, s);
    return s;
  }

  String String::operator+(double d) const
  {
    String s(*this);
    StringConversions::append(d, s);
    return s;
  }

  String String::operator+(long double ld) const
  {
    String s(*this);
    StringConversions::append(ld, s);
    return s;
  }

  String String::operator+(char c) const
  {
    String tmp(*this);
    tmp.push_back(c);
    return tmp;
  }

  String String::operator+(const char* s) const
  {
    String tmp(*this);
    tmp.append(s);
    return tmp;
  }

  String String::operator+(const String& s) const
  {
    String tmp(*this);
    tmp.insert(tmp.end(), s.begin(), s.end());
    return tmp;
  }

  String String::operator+(const std::string& s) const
  {
    String tmp(*this);
    tmp.insert(tmp.end(), s.begin(), s.end());
    return tmp;
  }

  String& String::operator+=(int i)
  {
    StringConversions::append(i, *this);
    return *this;
  }

  String& String::operator+=(unsigned int i)
  {
    StringConversions::append(i, *this);
    return *this;
  }

  String& String::operator+=(short int i)
  {
    StringConversions::append(i, *this);
    return *this;
  }

  String& String::operator+=(short unsigned int i)
  {
    StringConversions::append(i, *this);
    return *this;
  }

  String& String::operator+=(long int i)
  {
    StringConversions::append(i, *this);
    return *this;
  }

  String& String::operator+=(long unsigned int i)
  {
    StringConversions::append(i, *this);
    return *this;
  }

  String& String::operator+=(long long unsigned int i)
  {
    StringConversions::append(i, *this);
    return *this;
  }

  String& String::operator+=(float f)
  {
    StringConversions::append(f, *this);
    return *this;
  }

  String& String::operator+=(double d)
  {
    StringConversions::append(d, *this);
    return *this;
  }

  String& String::operator+=(long double d)
  {
    StringConversions::append(d, *this);
    return *this;
  }

  String& String::operator+=(char c)
  {
    this->append(String(c));
    return *this;
  }

  String& String::operator+=(const char* s)
  {
    this->append(s);
    return *this;
  }

  String& String::operator+=(const String& s)
  {
    this->append(s);
    return *this;
  }

  String& String::operator+=(const std::string& s)
  {
    this->append(s);
    return *this;
  }
  
  ::size_t hash_value(String const& s)
  {
    boost::hash<std::string> hasher;
    return hasher(static_cast<std::string>(s));
  }

} // namespace OpenMS
