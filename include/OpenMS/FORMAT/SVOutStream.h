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
// $Maintainer: $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef SVOUTSTREAM_H
#define SVOUTSTREAM_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <ostream>
#include <boost/math/special_functions/fpclassify.hpp> // for "isnan"

namespace OpenMS
{
	/**	
		@brief Stream class for writing to comma/tab/...-separated values files.
		
		Automatically inserts separators between items and handles quoting of strings. Requires @p std::endl as the line delimiter - @p "\n" won't be accepted.
		
		@ingroup Format
	*/
	class OPENMS_DLLAPI SVOutStream: public std::ostream
	{
	public:
		/**
			 @brief Constructor

			 @param out Output stream to write to (open file or @p cout)
			 @param sep Separator string (typically comma, semicolon, or tab)
			 @param replacement If @p quoting is @p NONE, used to replace occurrences of @p sep within strings before writing them
			 @param quoting Quoting method for strings (see @p String::quote)
		*/
		SVOutStream(std::ostream& out, const String& sep = "\t",
								const String& replacement = "_",
								String::QuotingMethod quoting = String::DOUBLE);

		
		/**
			 @brief Stream output operator for @p String

			 The argument is quoted before writing; it must not contain the newline character
		*/
		SVOutStream& operator<<(String str); // use call-by-value here

		
		/**
			 @brief Stream output operator for @p std::string

			 The argument is quoted before writing; it must not contain the newline character
		*/
		SVOutStream& operator<<(const std::string& str);

		
		/**
			 @brief Stream output operator for @p char*

			 The argument is quoted before writing; it must not contain the newline character
		*/
		SVOutStream& operator<<(const char* c_str);

		
		/**
			 @brief Stream output operator for @p char
			 
			 The argument is quoted before writing; it must not contain the newline character
		*/
		SVOutStream& operator<<(const char c);

		
		/// Stream output operator for manipulators (used to catch @p std::endl)
		SVOutStream& operator<<(std::ostream& (*fp)(std::ostream&));

		
		/// Generic stream output operator (for non-character-based types)
		template <typename T>
		SVOutStream& operator<<(const T& value)
			{
				if (!newline_) (std::ostream&)*this << sep_;
				else newline_ = false;
				(std::ostream&)*this << value;
				return *this;
			}

		
		/// Unformatted output (no quoting: useful for comments, but use only on a line of its own!)
		SVOutStream& write(const String& str); // write unmodified string

		
		/**
			 @brief Switch modification of strings (quoting/replacing of separators) on/off

			 @return previous modification state
		*/
		bool modifyStrings(bool modify);

		
		/// Write a numeric value or "nan", if applicable (would not be needed for Linux)
		template <typename NumericT>
		SVOutStream& writeValueOrNan(NumericT thing)
			{
				if (!boost::math::isnan(thing)) return operator<<(thing);
				bool old = modifyStrings(false);
				operator<<(nan_);
				modifyStrings(old);
				return *this;
			}

	protected:
		/// Separator string
		String sep_;

		/// Replacement for separator
		String replacement_;

		/// String to use for NaN values
		String nan_;
		
		/// String quoting method
		String::QuotingMethod quoting_;

		/// On/off switch for modification of strings
		bool modify_strings_;
		
		/// Are we at the beginning of a line? (Otherwise, insert separator before next item.)
		bool newline_;
	};

}

#endif
