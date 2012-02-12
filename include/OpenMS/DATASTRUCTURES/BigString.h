// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_DATASTRUCTURES_BIGSTRING_H
#define OPENMS_DATASTRUCTURES_BIGSTRING_H

#include <vector>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS {

/**
@brief concatenates Proteins given as FASTAEntry to one big string separated by a unique character (by default $)

Concatenates the strings given as FASTAEntry separating them with a unique character and storing the headers of FASTAEntry as well as the position of separator characters. So a substring can be accessed easily and the corresponding header can be found fast by using bineary search.
*/
class OPENMS_DLLAPI BigString
{

	public:

	typedef std::pair<String,String> FASTAEntry;

	/**
	@brief constructor
	*/
	BigString ();

	/**
	@brief copy constructor
	*/
	BigString (const BigString & bs);

	/**
	@brief desctructor
	*/
	virtual ~BigString ();

	/**
	@brief add new string to bigString
	@param new_entry FASTAEntry to be added to big_string
	*/
	void add (FASTAEntry const & new_entry);

	/**
	@brief setter for separator character by default $
	@param sep separator character
	*/
	void setSeparator (const char sep);

	/**
	@brief getter for separator character
	@return separator character
	*/
	char getSeparator ();

	/**
	@brief returns the number of strings
	@return int with number of strings
	*/
	Size size ();

	/**
	@brief length of bigString
	@return int with length of the created bigString
	*/
	Size length ();

	/**
	@brief getPeptide from start position with given length this includes FASTAHeader
	@param entry contains the entry of the given range after calling
	@param start start index
	@param length length of desired substring
	@return FASTAEntry describing the protein
	@throw InvalidValue if a peptide is part of two different fasta entrys
	*/
	void getPeptide (FASTAEntry& entry, Size start, Size length);

	/**
	@brief returns bigString
	@return const reference to bigString
	*/
	const String & getBigString () const;

	protected:

	/**
	@brief private function to implement binary search
	@param index
	@param start start index
	@param end end inxed
	@return int with index
	*/
	Size getIndex_ (Size index, Size start, Size end);

	/**
	@brief retrieves index of inserted protein by bigStringPosition
	@param index
	@return int with index
	*/
	Size getIndex_ (Size index);

	String big_string_; ///< concatenated String

	char separator_; ///< separator sign

	Size count_; ///< number of Strings added to big_string

	Size len_; ///< length of the big_string

	std::vector<Size> sep_indices_; ///< indices of separators

	std::vector<String> FASTA_header_; ///< vector with headers of FASTAEntry

};
}
#endif // OPENMS_DATASTRUCTURES_BIGSTRING_H
