// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>

#include <string>
#include <cstring>
#include <vector>

class QString;

namespace OpenMS
{
  class DataValue;
  template <typename FloatingPointType>
  struct PrecisionWrapper;
  /**
      @brief A more convenient string class.

      It is based on std::string but adds a lot of methods for convenience.

      @ingroup Datastructures
  */
  class String :
    public std::string
  {
public:

    /// Empty string for comparisons
    OPENMS_DLLAPI static const String EMPTY;

    /** @name Type definitions
    */
    //@{
    /// Iterator
    typedef iterator    Iterator;
    /// Const Iterator
    typedef const_iterator  ConstIterator;
    /// Reverse Iterator
    typedef reverse_iterator    ReverseIterator;
    /// Const reverse Iterator
    typedef const_reverse_iterator  ConstReverseIterator;
    /// UInt type
    typedef size_type   SizeType;

    /// How to handle embedded quotes when quoting strings
    enum QuotingMethod {NONE, ESCAPE, DOUBLE};

    //@}

    /**	@name Constructors
    */
    //@{
    /// Default constructor
    OPENMS_DLLAPI String();
    /// Copy constructor
    OPENMS_DLLAPI String(const String&) = default;
    /// Move constructor
    OPENMS_DLLAPI String(String&&) = default;
    /// Constructor from std::string
    OPENMS_DLLAPI String(const std::string& s);
    /// Constructor from std::string_view
    OPENMS_DLLAPI String(const std::string_view& sv);
    /// Constructor from Qt QString
    OPENMS_DLLAPI String(const QString& s);
    /// Constructor from char*
    OPENMS_DLLAPI String(const char* s);
    /// Constructor from a char
    OPENMS_DLLAPI String(const char c);
    /// Constructor from char* (only @p length characters)
    OPENMS_DLLAPI String(const char* s, SizeType length);
    /// Constructor from char (repeats the char @p len times)
    OPENMS_DLLAPI String(size_t len, char c);
    /// Constructor from a char range
    template <class InputIterator>
    String(InputIterator first, InputIterator last) :
      std::string(first, last)
    {
    }

    /// Constructor from an integer
    OPENMS_DLLAPI String(int i);
    /// Constructor from an unsigned integer
    OPENMS_DLLAPI String(unsigned int i);
    /// Constructor from an integer
    OPENMS_DLLAPI String(short int i);
    /// Constructor from an unsigned integer
    OPENMS_DLLAPI String(short unsigned int i);
    /// Constructor from an integer
    OPENMS_DLLAPI String(long int i);
    /// Constructor from an unsigned integer
    OPENMS_DLLAPI String(long unsigned int i);
    /// Constructor from an unsigned integer
    OPENMS_DLLAPI String(long long unsigned int i);
    /// Constructor from an unsigned integer
    OPENMS_DLLAPI String(long long signed int i);
    /// Constructor from float (@p full_precision controls number of fractional digits, 3 digits when false, and 6 when true)
    OPENMS_DLLAPI String(float f, bool full_precision = true);
    /// Constructor from double (@p full_precision controls number of fractional digits, 3 digits when false, and 15 when true)
    OPENMS_DLLAPI String(double d, bool full_precision = true);
    /// Constructor from long double (@p full_precision controls number of fractional digits, 3 digits when false, and 15 when true)
    OPENMS_DLLAPI String(long double ld, bool full_precision = true);
    /// Constructor from DataValue (@p full_precision controls number of fractional digits for all double types or lists of double, 3 digits when false, and 15 when true)
    OPENMS_DLLAPI String(const DataValue& d, bool full_precision = true);

    //@}

    /** @name Predicates
    */
    //@{
    /// true if String begins with @p string, false otherwise
    OPENMS_DLLAPI bool hasPrefix(const String& string) const;

    /// true if String ends with @p string, false otherwise
    OPENMS_DLLAPI bool hasSuffix(const String& string) const;

    /// true if String contains the @p string, false otherwise
    OPENMS_DLLAPI bool hasSubstring(const String& string) const;

    /// true if String contains the @p byte, false otherwise
    OPENMS_DLLAPI bool has(Byte byte) const;
    //@}

    /// Assignment operator
    OPENMS_DLLAPI String& operator=(const String&) = default;
    /// Move assignment operator
    OPENMS_DLLAPI String& operator=(String&&) & = default;

    /** @name Accessors
    */
    //@{
    /**
      @brief returns the prefix of length @p length

      @exception Exception::IndexOverflow is thrown if @p length is bigger than the size
    */
    OPENMS_DLLAPI String prefix(SizeType length) const;

    /**
      @brief returns the suffix of length @p length

      @exception Exception::IndexOverflow is thrown if @p length is bigger than the size
    */
    OPENMS_DLLAPI String suffix(SizeType length) const;

    /**
      @brief returns the prefix of length @p length

      @exception Exception::IndexUnderflow is thrown if @p length is smaller than zero
      @exception Exception::IndexOverflow is thrown if @p length is bigger than the size
    */
    OPENMS_DLLAPI String prefix(Int length) const;

    /**
      @brief returns the suffix of length @p length

      @exception Exception::IndexUnderflow is thrown if @p length is smaller than zero
      @exception Exception::IndexOverflow is thrown if @p length is bigger than the size
    */
    OPENMS_DLLAPI String suffix(Int length) const;

    /**
      @brief returns the prefix up to the first occurrence of char @p delim (excluding it)

      @exception Exception::ElementNotFound is thrown if @p delim is not found
    */
    OPENMS_DLLAPI String prefix(char delim) const;

    /**
      @brief returns the suffix up to the last occurrence of char @p delim (excluding it)

      @exception Exception::ElementNotFound is thrown if @p delim is not found
    */
    OPENMS_DLLAPI String suffix(char delim) const;

    /**
     @brief Wrapper for the STL substr() method. Returns a String object with its contents initialized to a substring of the current object.

     @param pos Position of a character in the current string object to be used as starting character for the substring.
     If the @p pos is past the end of the string, it is set to the end of the string.

     @param n Length of the substring.
     If this value would make the substring to span past the end of the current string content, only those characters until the end of the string are used.
     npos is a static member constant value with the greatest possible value for an element of type size_t, therefore, when this value is used, all the
     characters between pos and the end of the string are used as the initialization substring.

     */
    OPENMS_DLLAPI String substr(size_t pos = 0, size_t n = npos) const;

    /**
      @brief Returns a substring where @p n characters were removed from the end of the string.

      If @p n is greater than size(), the result is an empty string.

      @param n Number of characters that will be removed from the end of the string.
     */
    OPENMS_DLLAPI String chop(Size n) const;

    //@}


    /**
        @name Mutators

        All these methods return a reference to the string in order to make them chainable
    */
    //@{
    /// inverts the direction of the string
    OPENMS_DLLAPI String& reverse();

    /// removes whitespaces (space, tab, line feed, carriage return) at the beginning and the end of the string
    OPENMS_DLLAPI String& trim();

    /**
        @brief Checks if the string is wrapped in quotation marks

        The quotation mark can be specified by parameter @p q (typically single or double quote).

        @see unquote()
    */
    OPENMS_DLLAPI bool isQuoted(char q = '"');

    /**
         @brief Wraps the string in quotation marks

         The quotation mark can be specified by parameter @p q (typically single or double quote); embedded quotation marks are handled according to @p method by backslash-escaping, doubling, or not at all.

         @see unquote()
    */
    OPENMS_DLLAPI String& quote(char q = '"', QuotingMethod method = ESCAPE);

    /**
         @brief Reverses changes made by the @p quote method

         Removes surrounding quotation marks (given by parameter @p q); handles embedded quotes according to @p method.

         @exception Exception::ConversionError is thrown if the string does not have the format produced by @p quote

         @see quote()
    */
    OPENMS_DLLAPI String& unquote(char q = '"', QuotingMethod method = ESCAPE);

    /// merges subsequent whitespaces to one blank character
    OPENMS_DLLAPI String& simplify();

    ///Adds @p c on the left side until the size of the string is @p size
    OPENMS_DLLAPI String& fillLeft(char c, UInt size);

    ///Adds @p c on the right side until the size of the string is @p size
    OPENMS_DLLAPI String& fillRight(char c, UInt size);

    ///Converts the string to uppercase
    OPENMS_DLLAPI String& toUpper();

    ///Converts the string to lowercase
    OPENMS_DLLAPI String& toLower();

    ///Converts the first letter of the string to uppercase
    OPENMS_DLLAPI String& firstToUpper();

    ///Replaces all occurrences of the character @p from by the character @p to.
    OPENMS_DLLAPI String& substitute(char from, char to);

    ///Replaces all occurrences of the string @p from by the string @p to.
    OPENMS_DLLAPI String& substitute(const String& from, const String& to);

    ///Remove all occurrences of the character @p what.
    OPENMS_DLLAPI String& remove(char what);

    ///Makes sure the string ends with the character @p end
    OPENMS_DLLAPI String& ensureLastChar(char end);

    ///removes whitespaces (space, tab, line feed, carriage return)
    OPENMS_DLLAPI String& removeWhitespaces();
    //@}

    /** @name Converters
    */
    //@{

    
    /**
        @brief Conversion to Int

        This method extracts only the integral part of the string.
        If you want the result rounded, use toFloat() and round the result.

        @exception Exception::ConversionError is thrown if the string could not be converted to Int
    */
    OPENMS_DLLAPI Int toInt() const;

    /**
        @brief Conversion to Int32

        This method extracts only the integral part of the string.
        If you want the result rounded, use toFloat() and round the result.

        @exception Exception::ConversionError is thrown if the string could not be converted to Int32
    */
    OPENMS_DLLAPI Int32 toInt32() const;

    /**
    @brief Conversion to Int64

    This method extracts only the integral part of the string.
    If you want the result rounded, use toFloat() and round the result.

    @exception Exception::ConversionError is thrown if the string could not be converted to Int64
    */
    OPENMS_DLLAPI Int64 toInt64() const;

    /**
      @brief Conversion to float

      @exception Exception::ConversionError is thrown if the string could not be converted to float
    */
    OPENMS_DLLAPI float toFloat() const;

    /**
      @brief Conversion to double

      @exception Exception::ConversionError is thrown if the string could not be converted to double
    */
    OPENMS_DLLAPI double toDouble() const;

    /// Conversion to Qt QString
    OPENMS_DLLAPI QString toQString() const;

    //@}

    /** @name Sum operator overloads
    */
    //@{
    /// Sum operator for an integer
    OPENMS_DLLAPI String operator+(int i) const;
    /// Sum operator for an unsigned integer
    OPENMS_DLLAPI String operator+(unsigned int i) const;
    /// Sum operator for an integer
    OPENMS_DLLAPI String operator+(short int i) const;
    /// Sum operator for an unsigned integer
    OPENMS_DLLAPI String operator+(short unsigned int i) const;
    /// Sum operator for an integer
    OPENMS_DLLAPI String operator+(long int i) const;
    /// Sum operator for an unsigned integer
    OPENMS_DLLAPI String operator+(long unsigned int i) const;
    /// Sum operator for an unsigned integer
    OPENMS_DLLAPI String operator+(long long unsigned int i) const;
    /// Sum operator for float
    OPENMS_DLLAPI String operator+(float f) const;
    /// Sum operator for double
    OPENMS_DLLAPI String operator+(double d) const;
    /// Sum operator for long double
    OPENMS_DLLAPI String operator+(long double ld) const;
    /// Sum operator for char
    OPENMS_DLLAPI String operator+(char c) const;
    /// Sum operator for char*
    OPENMS_DLLAPI String operator+(const char* s) const;
    /// Sum operator for String
    OPENMS_DLLAPI String operator+(const String& s) const;
    /// Sum operator for std::string
    OPENMS_DLLAPI String operator+(const std::string& s) const;
    //@}

    /** @name Append operator overloads
    */
    //@{
    /// Sum operator for an integer
    OPENMS_DLLAPI String& operator+=(int i);
    /// Sum operator for an unsigned integer
    OPENMS_DLLAPI String& operator+=(unsigned int i);
    /// Sum operator for an integer
    OPENMS_DLLAPI String& operator+=(short int i);
    /// Sum operator for an unsigned integer
    OPENMS_DLLAPI String& operator+=(short unsigned int i);
    /// Sum operator for an integer
    OPENMS_DLLAPI String& operator+=(long int i);
    /// Sum operator for an unsigned integer
    OPENMS_DLLAPI String& operator+=(long unsigned int i);
    /// Sum operator for an unsigned integer
    OPENMS_DLLAPI String& operator+=(long long unsigned int i);
    /// Sum operator for float
    OPENMS_DLLAPI String& operator+=(float f);
    /// Sum operator for double
    OPENMS_DLLAPI String& operator+=(double d);
    /// Sum operator for long double
    OPENMS_DLLAPI String& operator+=(long double d);
    /// Sum operator for char
    OPENMS_DLLAPI String& operator+=(char c);
    /// Sum operator for char*
    OPENMS_DLLAPI String& operator+=(const char* s);
    /// Sum operator for String
    OPENMS_DLLAPI String& operator+=(const String& s);
    /// Sum operator for std::string
    OPENMS_DLLAPI String& operator+=(const std::string& s);
    //@}

    ///returns a random string of the given length. It consists of [0-9a-zA-Z]
    OPENMS_DLLAPI static String random(UInt length);

    ///returns a string for @p d with exactly @p n decimal places
    OPENMS_DLLAPI static String number(double d, UInt n);

    /**
        @brief Returns a string with at maximum @p n characters for @p d

        If @p d is larger, scientific notation is used.
    */
    OPENMS_DLLAPI static String numberLength(double d, UInt n);

    /**
        @brief Splits a string into @p substrings using @p splitter as delimiter

        If @p splitter is not found, the whole string is put into @p substrings.
        If @p splitter is empty, the string is split into individual characters.
        If the invoking string is empty, @p substrings will also be empty.

        @p quote_protect (default: false) can be used to split only between quoted
        blocks e.g. ' "a string" , "another string with , in it" '
        results in only two substrings (with double quotation marks @em removed).
        Every returned substring is trimmed and then (if present) has surrounding quotation marks removed.

        @return @e true if one or more splits occurred, @e false otherwise

        @see concatenate().
    */
    OPENMS_DLLAPI bool split(const char splitter, std::vector<String>& substrings, bool quote_protect = false) const;

    /**
        @brief Splits a string into @p substrings using @p splitter (the whole string) as delimiter

        If @p splitter is not found,  the whole string is put into @p substrings.
        If @p splitter is empty, the string is split into individual characters.
        If the invoking string is empty, @p substrings will also be empty.

        @return @e true if one or more splits occurred, @e false otherwise

        @see concatenate().
    */
    OPENMS_DLLAPI bool split(const String& splitter, std::vector<String>& substrings) const;

    /**
        @brief Splits a string into @p substrings using @p splitter (the whole string) as delimiter, but does not split within quoted substrings

        A "quoted substring" has the format as produced by @p quote(q, method), where @p q is the quoting character and @p method defines the handling of embedded quotes. Substrings will not be "unquoted" or otherwise processed.

        If @p splitter is not found,  the whole string is put into @p substrings.
        If @p splitter or the invoking string is empty, @p substrings will also be empty.

        @return @e true if one or more splits occurred, @e false otherwise

        @exception Exception::ConversionError is thrown if quotation marks are not balanced

        @see concatenate(), quote().
    */
    OPENMS_DLLAPI bool split_quoted(const String& splitter, std::vector<String>& substrings,
                                    char q = '"', QuotingMethod method = ESCAPE) const;

    /**
        @brief Concatenates all elements from @p first to @p last-1 and inserts @p glue between the elements

        @see split().
    */
    template <class StringIterator>
    void concatenate(StringIterator first, StringIterator last, const String& glue = "")
    {
      //empty container
      if (first == last)
      {
        std::string::clear();
        return;
      }

      std::string::operator=(*first);
      for (StringIterator it = ++first; it != last; ++it)
      {
        std::string::operator+=(glue + (*it));
      }
    }
  };
  OPENMS_DLLAPI ::size_t hash_value(OpenMS::String const& s);
} // namespace OpenMS

namespace std
{
  template <> struct hash<OpenMS::String> //hash for String
  {
    std::size_t operator()( OpenMS::String const& s) const
    {
      return std::hash<string>()(static_cast<string>(s));
    }
  };
} // namespace std
