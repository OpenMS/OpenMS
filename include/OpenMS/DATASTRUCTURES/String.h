// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
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
		
		It based on std::string but adds a lot of methods for convenience.
		
		@ingroup Datastructures
	*/
	class String:
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
		/// Constructor from float (precision is 7)
		String(float f);
		/// Constructor from double (precision is 10)
		String(double d);
		/// Constructor from long double (precision is 16)
		String(long double d);
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
			
			@param start  start position of the substring.<br> 
			              If start is negative, the returned string will start at the start'th character from the end of string.
			@param n length of the substring.<br> 
			         If a negative length is given, then that many characters will be omitted from the end of string.
		*/
		String substr(Int start, Int n) const;
		/**
			@brief Returns the suffix of the string from position @p start
			
			@param start  start position of the substring.<br> 
			              If start is negative, the returned string will start at the start'th character from the end of string.
		*/
		String substr(Int start) const;
		
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
		/// Sum operator for float (precision is 7)
		String operator+ (float f) const;
		/// Sum operator for double (precision is 10)
		String operator+ (double d) const;
		/// Sum operator for long double (precision is 16)
		String operator+ (long double d) const;
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
		/// Sum operator for float (precision is 7)
		String& operator+= (float f);
		/// Sum operator for double (precision is 10)
		String& operator+= (double d);
		/// Sum operator for long double (precision is 16)
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
			@brief returns a string with at maximum @p n characters for @p d
		
		 If @p d is larger, scientific notation is used.
		*/
		static String numberLength(DoubleReal d, UInt n);

		
		/**
			@brief splits a string into @p substrings using @p splitter as delimiter
			
			If the @p splitter is not found, @p substrings is empty.
			@return if the splitter was found (once or multiple times) in the string

			@see implode().
		*/
		bool split(char splitter, std::vector<String>& substrings) const;
		
		
		/**
			@brief Concatenates all elements from @p first to @p last-1 and inserts @p glue between the elements

			@see split().
		*/
		template<class StringIterator>
		void implode(StringIterator first, StringIterator last, const String& glue = "")
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
