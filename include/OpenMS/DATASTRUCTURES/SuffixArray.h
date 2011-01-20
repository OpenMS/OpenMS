// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------


#ifndef OPENMS_DATASTRUCTURES_SUFFIXARRAY_H
#define OPENMS_DATASTRUCTURES_SUFFIXARRAY_H

#include <vector>
#include <OpenMS/CONCEPT/Exception.h>


namespace OpenMS 
{
	class String;

/**
	@brief abstract class for suffix array
*/

class OPENMS_DLLAPI SuffixArray 
{
	
public:

	/**
	@brief constructor taking the string and the filename for writing or reading
	@param st the string as const reference with which the suffix array will be build
	@param filename the filename for writing or reading the suffix array
	@throw Exception::InvalidValue if string does not start with empty string ($)
	@throw Exception::FileNotFound is thrown if the given filename is not found
	*/
	SuffixArray(const String & st, const String & filename);

	/**
	@brief copy constructor
	*/
	SuffixArray(const SuffixArray & sa);

	/**
	@brief destructor
	*/
	virtual ~SuffixArray() = 0;

	/**
	@brief transforms suffix array to a printable String
	*/
	virtual String toString() = 0;

	/**
	@brief the function that will find all peptide candidates for a given spectrum
	@param spec const reference of DoubleReal vector describing the spectrum
	@param candidates the candidates which are returned for the masses given in spec
	@return a vector of SignedSize pairs.
	@throw InvalidValue if the spectrum is not sorted ascendingly
	
	*/
	virtual void findSpec(std::vector<std::vector<std::pair<std::pair<SignedSize, SignedSize>,DoubleReal > > >& candidates, const std::vector<DoubleReal> & spec) = 0;

	/**
	@brief saves the suffix array to disc
	@param filename const reference string describing the filename
	@return bool if operation was succesful
	@throw UnableToCreateFile if file could not be created (e.x. if you have no rigths)
	*/
	virtual bool save(const String & filename) = 0;
	/**
	@brief opens the suffix array
	@param filename const reference string describing the filename
	@return bool if operation was succesful
	@throw FileNotFound
	*/
	virtual bool open(const String & filename) = 0;

	/**
	@brief setter for tolerance
	@param t DoubleReal with tolerance
	@throw InvalidValue if tolerance is negative
	*/
	virtual void setTolerance (DoubleReal t) = 0;

	/**
	@brief getter for tolerance
	@return DoubleReal with tolerance
	*/
	virtual DoubleReal getTolerance () const  = 0;

	/**
	@brief returns if an enzyme will cut after first character
	@param aa1 const char as first aminoacid
	@param aa2 const char as second aminoacid
	@return bool describing if it is a digesting site
	*/
	virtual bool isDigestingEnd(const char aa1, const char aa2) const = 0;

	/**
	@brief setter for tags
	@param tags const vector of strings with tags with length 3 each
	@throw Exception::InvalidValue if at least one tag does not have size of 3
	*/
	virtual void setTags(const std::vector<String>& tags) = 0;

	/**
	@brief getter for tags
	@return const vector of string with tags
	*/
	virtual const std::vector<String>& getTags () = 0;

	/**
	@brief setter for use_tags
	@param use_tags indicating whether tags should be used or not
	*/
	virtual void setUseTags(bool use_tags) = 0;

	/**
	@brief getter for use_tags
	@return bool indicating whether tags are used or not
	*/
	virtual bool getUseTags() = 0;

	/**
	@brief setter for number of modifications
	@param number_of_mods
	*/
	virtual void setNumberOfModifications(Size number_of_mods) = 0;

	/**
	@brief getter for number of modifications
	@return Size describing number of modifications
	*/
	virtual Size getNumberOfModifications () = 0;

	/**
	@brief output for statistic
	*/
	virtual void printStatistic() = 0;

	/**
	@brief constructor
	*/
	SuffixArray ();	

	
};
}

#endif //OPENMS_DATASTRUCTURES_SUFARRAY_H
