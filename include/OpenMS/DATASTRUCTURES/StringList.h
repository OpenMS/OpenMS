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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_STRINGLIST_H
#define OPENMS_DATASTRUCTURES_STRINGLIST_H

#include <OpenMS/DATASTRUCTURES/String.h>

#ifdef OPENMS_COMPILER_MSVC
	#pragma warning( push )
	#pragma warning( disable : 4251 ) // disable MSVC dll-interface warning
#endif

namespace OpenMS
{
	/**
		@brief String list

		This class is based on std::vector<String> but adds some methods for convenience.

		@ingroup Datastructures
	*/
	class OPENMS_DLLAPI StringList
		: public std::vector<String>
	{
		public:

	 		/** @name Type definitions
			*/
			//@{
			/// Mutable iterator
			typedef iterator	Iterator;
			/// Non-mutable iterator
			typedef const_iterator	ConstIterator;
			/// Mutable reverse iterator
			typedef reverse_iterator	ReverseIterator;
			/// Non-mutable reverse iterator
			typedef const_reverse_iterator	ConstReverseIterator;
			//@}

			///@name Constructors and assignment operators
			//@{
			/// Default constructor
			StringList();
			/// Copy constructor
			StringList(const StringList& rhs);
			/// Constructor from vector<String>
			StringList(const std::vector<String>& rhs);
			/// Constructor from vector<string>
			StringList(const std::vector<std::string>& rhs);
			///  Assignment operator
			StringList& operator=(const StringList& rhs);
			///  Assignment operator from vector<String>
			StringList& operator=(const std::vector<String>& rhs);
			///  Assignment operator vector<string>
			StringList& operator=(const std::vector<std::string>& rhs);
			//@}

			///@name Search methods
			//@{
						/**
    		@brief Searches for the first line that starts with @p text beginning at line @p start
    		
    		@param start the line to start the search in
    		@param text the text to find
    		@param trim whether the line is trimmed before
    		@return returns an iterator to the matching line. If no line matches, end() is returned
    	*/
			Iterator search(const Iterator& start, const String& text, bool trim=false);

			/**
				@brief Searches for the first line that starts with @p text
				
				This is an overloaded member function, provided for convenience.<br>
				It behaves essentially like the above function but the search is start at the beginning of the file
    	*/
			Iterator search(const String& text, bool trim=false);

			/**
    		@brief Searches for the first line that ends with @p text beginning at line @p start
    		
    		@param start the line to start the search in
    		@param text the text to find
    		@param trim whether the line is trimmed before
    		@return returns an iterator to the matching line. If no line matches, end() is returned
    	*/
			Iterator searchSuffix(const Iterator& start, const String& text, bool trim=false);

			/**
				@brief Searches for the first line that ends with @p text
				
				This is an overloaded member function, provided for convenience.
				
				It behaves essentially like searchSuffix(const Iterator&, const String&, bool) but the search starts at the beginning of the file
    	*/
			Iterator searchSuffix(const String& text, bool trim=false);

      /**
        @brief Searches for the first line that starts with @p text beginning at line @p start

        @param start the line to start the search in
        @param text the text to find
        @param trim whether the line is trimmed before
        @return returns an iterator to the matching line. If no line matches, end() is returned
      */
      ConstIterator search(const ConstIterator& start, const String& text, bool trim=false) const;

      /**
        @brief Searches for the first line that starts with @p text

        This is an overloaded member function, provided for convenience.<br>
        It behaves essentially like the above function but the search is start at the beginning of the file
      */
      ConstIterator search(const String& text, bool trim=false) const;

      /**
        @brief Searches for the first line that ends with @p text beginning at line @p start

        @param start the line to start the search in
        @param text the text to find
        @param trim whether the line is trimmed before
        @return returns an iterator to the matching line. If no line matches, end() is returned
      */
      ConstIterator searchSuffix(const ConstIterator& start, const String& text, bool trim=false) const;

      /**
        @brief Searches for the first line that ends with @p text

        This is an overloaded member function, provided for convenience.

        It behaves essentially like searchSuffix(const Iterator&, const String&, bool) but the search starts at the beginning of the file
      */
      ConstIterator searchSuffix(const String& text, bool trim=false) const;
			//@}

			///Operator for appending entries with less code
			template<typename StringType>
			StringList& operator<<(const StringType& string)
			{
				this->push_back(string);
				return *this;
			}

			/// Returns a list that is created by splitting the given (comma-separated) string (String are not trimmed!)
			static StringList create(const String& list, const char splitter=',');
			/// Returns a list that is created from an array of char*
			static StringList create(const char * const * list, UInt size );
			/// Returns if a string is contained in the list
			bool contains(const String& s) const;
			/// Transforms all contained strings to upper case
			void toUpper();
			/// Transforms all contained strings to lower case
			void toLower();

			/// Concatenate the string elements and putting the @p glue string between elements
			String concatenate(const String& glue="") const;

			/// output stream operator
			friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const StringList& p);

	};

} // namespace OPENMS

#ifdef OPENMS_COMPILER_MSVC
	#pragma warning( pop )
#endif

#endif // OPENMS_DATASTRUCTURES_STRINGLIST_H
