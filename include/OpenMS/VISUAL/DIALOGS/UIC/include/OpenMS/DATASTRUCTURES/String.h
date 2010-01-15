// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_STRING_H
#define OPENMS_DATASTRUCTURES_STRING_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <string>
#include <vector>

class QString;

namespace OpenMS
{
	class DataValue;
	/**	
		@brief A more convenient string class.
		
		It is based on std::string but adds a lot of methods for convenience.
		
		@ingroup Datastructures
	*/
	class OPENMS_DLLAPI String:
		public std::string
	{
		public:

		/// Empty string for comparisons
		static const String EMPTY;
		
		/** @name Type definitions
		*/
		//@{	
		/// Iterator
		typedef iterator	Iterator;
		/// Const Iterator
		typedef const_iterator	ConstIterator;
		/// Reverse Iterator
		typedef reverse_iterator	ReverseIterator;
		/// Const reverse Iterator
		typedef const_reverse_iterator	ConstReverseIterator;
		/// UInt type
		typedef size_type	SizeType;

		/// How to handle embedded quotes when quoting strings
		enum QuotingMethod {NONE, ESCAPE, DOUBLE};
		
		//@}		
		
		/**	@name Constructors
		*/
		//@{
		/// Default constructor
		String();
		/// Constructor from std::string
		String(const std::string& s);
		/// Constructor from Qt QString
		String(const QString& s);
		/// Constructor from char*
		String(const char* s);
		/// Constructor from a char
		String(const char c);
		/// Constructor from char* (only @p length chracters)
		String(const char* s, SizeType length);
		/// Constructor from char (repeats the char @p len times)
		String(size_t len, char c);
		/// Constructor from a char range
		template <class InputIterator>
		String(InputIterator first, InputIterator last)
		: std::string(first,last)
		{
			
		}
		/// Constructor from an integer
		String(int i);
		/// Constructor from an unsigned integer
		String(unsigned int i);
		/// Constructor from an integer
		String(short int i);
		/// Constructor from an unsigned integer
		String(short unsigned int i);
		/// Constructor from an integer
		String(long int i);
		/// Constructor from an unsigned integer
		String(long unsigned int i);
		/// Constructor from an unsigned integer
		String(long long unsigned int i);
		/// Constructor from an unsigned integer
		String(long long signed int i);
		/// Constructor from float
		String(float f);
		/// Constructor from double
		String(double d);
		/// Constructor from long double
		String(long double ld);
		/// Constructor from DataValue (casted to String)
		String(const DataValue& d);

		//@}

		/** @name Predicates
		*/
		//@{
		/// true if String begins with @p string, false otherwise
		bool hasPrefix(const String& string) const;

		/// true if String ends with @p string, false otherwise
		bool hasSuffix(const String& string) const;

		/// true if String contains the @p string, false otherwise
		bool hasSubstring(const String& string) const;

		/// true if String contains the @p byte, false otherwise
		bool has(Byte byte) const;
		//@}


		/** @name Accessors
		*/
		//@{
		/**
		  @brief returns the prefix of length @p length

		  @exception Exception::IndexOverflow is thrown if @p length is bigger than the size
		*/
		String prefix(SizeType length) const;
		
		/**
		  @brief returns the suffix of length @p length

		  @exception Exception::IndexOverflow is thrown if @p length is bigger than the size
		*/
		String suffix(SizeType length) const;

		/**
		  @brief returns the prefix of length @p length

		  @exception Exception::IndexUnderflow is thrown if @p length is smaller than zero
		  @exception Exception::IndexOverflow is thrown if @p length is bigger than the size
		*/
		String prefix(Int length) const;
		
		/**
		  @brief returns the suffix of length @p length

		  @exception Exception::IndexUnderflow is thrown if @p length is smaller than zero
		  @exception Exception::IndexOverflowis thrown if @p length is bigger than the size
		*/
		String suffix(Int length) const;
		
		/**
		  @brief returns the prefix up to the first occurence of char @p delim (excluding it)

		  @exception Exception::ElementNotFound is thrown if @p delim is not found
		*/
		String prefix(char delim) const;
		
		/**
		  @brief returns the suffix up to the last occurence of char @p delim (excluding it)

		  @exception Exception::ElementNotFound is thrown if @p delim is not found
		*/
		String suffix(char delim) const;
		/**
			@brief Returns a substring
			
			If start or end are out of bounds, they are automatically corrected - set to begin or end, respectively.
			
			If the end is positioned before the start an empty string is returned.
			
			@param start  start position of the substring.<br> 
			              If start is negative, the returned string will start at the start'th character from the end of string.
			@param n length of the substring.<br> 
			         If a negative length is given, then that many characters will be omitted from the end of string.
		*/
		String substr(SignedSize start, SignedSize n) const;
		/**
			@brief Returns the suffix of the string from position @p start

			If start is out of bounds, it is automatically corrected - set to begin or end, respectively.
			
			@param start  start position of the substring.<br> 
			              If start is negative, the returned string will start at the start'th character from the end of string.
		*/
		String substr(SignedSize start) const;
		
		//@}
		
		
		/** 
			@name Mutators
			
			All these methods return a reference to the string in order to make them chainable
		*/
		//@{
		/// inverts the direction of the string
		String& reverse();
		
		/// removes whitespaces (space, tab, line feed, carriage return) at the beginning and the end of the string
		String& trim();

		/**
			 @brief Wraps the string in quotation marks

			 The quotation mark can be specified by parameter @p q (typically single or double quote); embedded quotation marks are handled according to @p method by backslash-escaping, doubling, or not at all.

			 @see unquote()
		*/
		String& quote(char q = '"', QuotingMethod method = ESCAPE);

		/**
			 @brief Reverses changes made by the @p quote method

			 Removes surrounding quotation marks (given by parameter @p q); handles embedded quotes according to @p method.

			 @exception Exception::ConversionError is thrown if the string does not have the format produced by @p quote

			 @see quote()
		*/
		String& unquote(char q = '"', QuotingMethod method = ESCAPE);
		
		/// merges subsequent whitespaces to one blank character
		String& simplify();
		
		///Adds @p c on the left side until the size of the string is @p size
		String& fillLeft(char c, UInt size);
		
		///Adds @p c on the right side until the size of the string is @p size
		String& fillRight(char c, UInt size);

		///Converts the string to uppercase
		String& toUpper();
		
		///Converts the string to lowercase
		String& toLower();

		///Converts the first letter of the string to uppercase
		String& firstToUpper();

		///Replaces all occurences of the character @p from by the character @p to.
		String& substitute(char from, char to);

		///Replaces all occurences of the string @p from by the string @p to.
		String& substitute(const String& from, const String& to);

		///Remove all occurences of the character @p what.
		String& remove(char what);

		///Makes sure the string ends with the character @p end
		String& ensureLastChar(char end);

		///removes whitespaces (space, tab, line feed, carriage return)
		String& removeWhitespaces();
		//@}

		/** @name Converters
		*/
		//@{

		/**
			@brief Conversion to int
		
			This method extracts only the integral part of the string.
			If you want the result rounded, use toFloat() and round the result.
			
			@exception Exception::ConversionError is thrown if the string could not be converted to int
		*/
		Int toInt() const;

		/**
		  @brief Conversion to float

		  @exception Exception::ConversionError is thrown if the string could not be converted to float
		*/
		Real toFloat() const;

		/**
		  @brief Conversion to double

		  @exception Exception::ConversionError is thrown if the string could not be converted to double
		*/
		DoubleReal toDouble() const;

		/// Conversion to Qt QString
		QString toQString() const;

		//@}

		/** @name Sum operator overloads
		*/
		//@{
		/// Sum operator for an integer
		String operator+ (int i) const;
		/// Sum operator for an unsigned integer
		String operator+ (unsigned int i) const;
		/// Sum operator for an integer
		String operator+ (short int i) const;
		/// Sum operator for an unsigned integer
		String operator+ (short unsigned int i) const;
		/// Sum operator for an integer
		String operator+ (long int i) const;
		/// Sum operator for an unsigned integer
		String operator+ (long unsigned int i) const;
		/// Sum operator for an unsigned integer
		String operator+ (long long unsigned int i) const;
		/// Sum operator for float
		String operator+ (float f) const;
		/// Sum operator for double
		String operator+ (double d) const;
		/// Sum operator for long double
		String operator+ (long double ld) const;
		/// Sum operator for char
		String operator+ (char c) const;
		/// Sum operator for char*
		String operator+ (const char* s) const;
		/// Sum operatr for String
		String operator+ (const String& s) const;
		/// Sum operator for std::string
		String operator+ (const std::string& s) const;
		//@}

		/** @name Append operator overloads
		*/
		//@{
		/// Sum operator for an integer
		String& operator+= (int i);
		/// Sum operator for an unsigned integer
		String& operator+= (unsigned int i);
		/// Sum operator for an integer
		String& operator+= (short int i);
		/// Sum operator for an unsigned integer
		String& operator+= (short unsigned int i);
		/// Sum operator for an integer
		String& operator+= (long int i);
		/// Sum operator for an unsigned integer
		String& operator+= (long unsigned int i);
		/// Sum operator for an unsigned integer
		String& operator+= (long long unsigned int i);
		/// Sum operator for float
		String& operator+= (float f);
		/// Sum operator for double
		String& operator+= (double d);
		/// Sum operator for long double
		String& operator+= (long double d);
		/// Sum operator for char
		String& operator+= (char c);
		/// Sum operator for char*
		String& operator+= (const char* s);
		/// Sum operatr for String
		String& operator+= (const String& s);
		/// Sum operator for std::string
		String& operator+= (const std::string& s);
		//@}

		///returns a random string of the given length. It consists of [0-9a-zA-Z]
		static String random(UInt length);

		///returns a string for @p d with exactly @p n decimal places
		static String number(DoubleReal d, UInt n);
		/**
			@brief Returns a string with at maximum @p n characters for @p d
		
			If @p d is larger, scientific notation is used.
		*/
		static String numberLength(DoubleReal d, UInt n);

		
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
		bool split(const char splitter, std::vector<String>& substrings, bool quote_protect=false) const;

		/**
			@brief Splits a string into @p substrings using @p splitter (the whole string) as delimiter
			
			If @p splitter is not found,  the whole string is put into @p substrings.
			If @p splitter is empty, the string is split into individual characters.
			If the invoking string is empty, @p substrings will also be empty.
			
			@return @e true if one or more splits occurred, @e false otherwise

			@see concatenate().
		*/
		bool split(const String& splitter, std::vector<String>& substrings) const;

		/**
			@brief Splits a string into @p substrings using @p splitter (the whole string) as delimiter, but does not split within quoted substrings

			A "quoted substring" has the format as produced by @p quote(q, method), where @p q is the quoting character and @p method defines the handling of embedded quotes. Substrings will not be "unquoted" or otherwise processed.
			
			If @p splitter is not found,  the whole string is put into @p substrings.
			If @p splitter or the invoking string is empty, @p substrings will also be empty.		
			
			@return @e true if one or more splits occurred, @e false otherwise
			
			@exception Exception::ConversionError is thrown if quotation marks are not balanced
			
			@see concatenate(), quote().
		*/
		bool split_quoted(const String& splitter,	std::vector<String>& substrings,
											char q = '"', QuotingMethod method = ESCAPE) const;

		/**
			@brief Concatenates all elements from @p first to @p last-1 and inserts @p glue between the elements

			@see split().
		*/
		template<class StringIterator>
		void concatenate(StringIterator first, StringIterator last, const String& glue = "")
		{
			//empty container
			if (first==last)
			{
				std::string::clear();
				return;
			}
		
			std::string::operator=(*first);
			for (StringIterator it = ++first; it != last; ++it)
			{
				std::string::operator+=( glue + (*it));
			}
		}
		
	};
	
} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_HASHMAP_H
