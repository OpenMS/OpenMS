// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_FORMAT_TEXTFILE_H
#define OPENMS_FORMAT_TEXTFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <vector>

namespace OpenMS
{
	/**
		@brief This class provides some basic file handling methods and facilitates reading, writing and handling text files.
  
  	@ingroup FileIO
	*/
  class TextFile
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

    	///Default constructor
			TextFile();

			/// destructor
			virtual ~TextFile();
			
    	/**
    		@brief Constructor with filename
    	
    		@param filename the filename
    		@param trim_lines wether or not the lines are trimmed when reading them from file
    	*/
			TextFile(const String& filename, bool trim_lines=false) throw (Exception::FileNotFound);

    	/**
    		@brief Loads data from file
    	
    		@param filename the filename
    		@param trim_lines wether or not the lines are trimmed when reading them from file
    	*/
			void load(const String& filename, bool trim_lines=false) throw (Exception::FileNotFound);

    	/**
    		@brief Writes the data to a file
    		
    		Note that this function uses unix-style linebreaks
    		@param filename the filename
    	*/
			void save(const String& filename) throw (Exception::UnableToCreateFile);

			/**
    		@brief Searches for the first line that starts with @p text beginning at line @p start
    		
    		@param start the line to start the search in
    		@param text the text to find
    		@param trim wether the line is trimmed before
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
    		@param trim wether the line is trimmed before
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
        @param trim wether the line is trimmed before
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
        @param trim wether the line is trimmed before
        @return returns an iterator to the matching line. If no line matches, end() is returned
      */
      ConstIterator searchSuffix(const ConstIterator& start, const String& text, bool trim=false) const;

      /**
        @brief Searches for the first line that ends with @p text

        This is an overloaded member function, provided for convenience.

        It behaves essentially like searchSuffix(const Iterator&, const String&, bool) but the search starts at the beginning of the file
      */
      ConstIterator searchSuffix(const String& text, bool trim=false) const;
			
			/// Return the content as a single String
			String asString() const;
			
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_TEXTFILE_H
