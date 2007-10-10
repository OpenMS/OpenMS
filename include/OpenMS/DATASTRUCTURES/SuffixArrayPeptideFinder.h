// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bauer$
// --------------------------------------------------------------------------


#ifndef OPENMS_DATASTRUCTURES_SUFFIXARRAYPEPTIDEFINDER_H
#define OPENMS_DATASTRUCTURES_SUFFIXARRAYPEPTIDEFINDER_H

#include <vector>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/BigString.h>

namespace OpenMS 
{
	class String;
	class SuffixArray;
	/**
		@brief wrapper for easy use of sufArray

		@todo Test currently does not run (core dump), neither on 32 nor on 64 bit, crashes somewere in SeqAn...; TODO contact SeqAn developers or Chris Bauer (Clemens, Andreas)

	*/
class SuffixArrayPeptideFinder {

public:

	/**
	@brief
	*/
	typedef std::pair <String, String> FASTAEntry;

	/**
	@brief constructor
	@param fFile const string for location of FASTA File
	*/
	SuffixArrayPeptideFinder(const String & fFile, const String & method) throw (Exception::FileNotFound,Exception::ParseError,Exception::InvalidValue);

	/**
	@brief copy constructor
	*/
	SuffixArrayPeptideFinder(const SuffixArrayPeptideFinder & source);

	/**
	@brief destructor
	*/
	virtual ~SuffixArrayPeptideFinder();

	/**
	@brief finds all candidate for given spectrum in the sufix array
	@param spec const reference to float vector describing the MS spectrum
	@return	for every mass a entry with all Candidates as vector of FASTAEntrys
	@see sufArray.h
	*/
	std::vector<std::vector<std::pair<FASTAEntry, String > > > getCandidates (const std::vector<double> & spec);

	/**
	@brief finds all candidate for given DTA file
	@param DTAFile DTA file location
	@return	for every mass a entry with all Candidates as vector of FASTAEntrys
	@throw Exception::FileNotFound if DTA file does not exists
	@see sufArray.h
	*/
	std::vector<std::vector<std::pair<FASTAEntry, String > > > getCandidates (const String & DTA_file) throw (Exception::FileNotFound,Exception::ParseError);

	/**
	@brief setter for tolerance
	@param t const float tolerance
	*/
	void setTolerance(const float t);

	/**
	@brief getter for tolerance
	@return float with tolerance
	*/
	float getTolerance() const;

	/**
	@brief setter for number of modifications
	@param number_of_mods
	*/
	void setNumberOfModifications(UInt number_of_mods) const;

	/**
	@brief getter for number of modifications
	@return number of modifications
	*/
	UInt getNumberOfModifications() const;

	/**
	@brief setter for tags
	@param tags reference to vector of strings with tags
	@note sets use_tags = true
	*/
	void setTags(const std::vector<OpenMS::String> & tags) throw (OpenMS::Exception::InvalidValue);

	/**
	@brief getter for tags
	@return const reference to vector of strings
	*/
	const std::vector<OpenMS::String> & getTags();

	/**
	@brief setter for use_tags
	@param use_tags indicating whether tags should be used or not
	*/
	void setUseTags(bool use_tags);

	/**
	@brief getter for use_tags
	@return bool indicating whether tags are used or not
	*/
	bool getUseTags();

	/**
	@brief setter for modification output method
	@param s describing how modifications sould be given back
	*/
	void setModificationOutputMethod(const String & s) throw (OpenMS::Exception::InvalidValue);

	/**
	@brief getter for modification output method
	@return String
	*/
	String getModificationOutputMethod();

protected:

	String vToString_(std::vector<String> v);

	BigString big_string_; 	///< bigString object holding all peptides of fasta file

	SuffixArray* sa_; 	///< pointer to sufixarray

	String modification_output_method_; ///< output method for modifications

};

}

#endif //OPENMS_EXAMPLES_SuffixArrayPeptideFinder_H
