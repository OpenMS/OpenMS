// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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


#ifndef OPENMS_DATASTRUCTURES_SUFFIXARRAYPEPTIDEFINDER_H
#define OPENMS_DATASTRUCTURES_SUFFIXARRAYPEPTIDEFINDER_H

#include <vector>
#include <OpenMS/DATASTRUCTURES/BigString.h>
#include <OpenMS/CHEMISTRY/WeightWrapper.h>

namespace OpenMS 
{
	class String;
	class SuffixArray;
	/**
		@brief wrapper for easy use of sufArray
	*/
class OPENMS_DLLAPI SuffixArrayPeptideFinder
	: public WeightWrapper 
{

public:

	/**
	@brief
	*/
	typedef std::pair <String, String> FASTAEntry;

	/**
	@brief constructor
	@param filename FASTA File name
	@param method Name of the method used (trypticCompressed, seqan, trypticSeqan)
	@param weight_mode if not monoistopic weight should be used, this parameters can be set to AVERAGE
	@throw FileNotFound is thrown if the filename is not found
	@throw ParseError is thrown if a error in parsing of the fasta file occurs
	@throw InvalidValue is thrown if an unknown method is supplied 
	*/
	SuffixArrayPeptideFinder(const String& filename, const String& method, const WeightWrapper::WEIGHTMODE weight_mode = WeightWrapper::MONO);

	/**
	@brief copy constructor
	*/
	SuffixArrayPeptideFinder(const SuffixArrayPeptideFinder& source);

	/**
	@brief destructor
	*/
	virtual ~SuffixArrayPeptideFinder();

	/**
	@brief finds all candidates for given spectrum in the suffix array
	@param spec vector holding the mass values to query
	@param candidates Output holding the candidates for input masses (one vector per mass)
				 FASTAEntry contains the FASTA header and the peptide sequence
				 The String contains the modification (if any) in the format specified by getModificationOutputMethod()
	@see sufArray.h
	*/
	void getCandidates(std::vector<std::vector<std::pair<FASTAEntry, String > > >& candidates, const std::vector<DoubleReal> & spec);

	/**
	@brief finds all candidate for given DTA file
	@param DTA_file DTA file location
	@param candidates Output parameters which holds the candidates suitable for the mass given in the dta file
				 FASTAEntry contains the FASTA header and the peptide sequence
				 The String contains the modification (if any) in the format specified by getModificationOutputMethod()
	@throw FileNotFound if DTA file does not exists
	@throw ParseError is thrown if the dta file could not be parsed
	@see sufArray.h
	*/
	void getCandidates(std::vector<std::vector<std::pair<FASTAEntry, String > > >& candidates, const String & DTA_file);

	/**
	@brief allowed tolerance for mass match
	@param t Tolerance in u
	*/
	void setTolerance(const DoubleReal t);

	/**
	@brief allowed tolerance for mass match
	@return Tolerance in u
	*/
	DoubleReal getTolerance() const;

	/**
	@brief setter for number of modifications
	@param number_of_mods
	*/
	void setNumberOfModifications(Size number_of_mods) const;

	/**
	@brief getter for number of modifications
	@return number of modifications
	*/
	Size getNumberOfModifications() const;

	/**
	@brief setter for tags
	@param tags reference to vector of strings with tags
	@note sets use_tags = true
	*/
	void setTags(const std::vector<String>& tags);

	/**
	@brief getter for tags
	@return const reference to vector of strings
	*/
	const std::vector<String>& getTags();

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
	@brief set modification output method (valid are: "mass", "stringUnchecked", "stringChecked")
	@param s describing how modifications sould be given back
	@throw InvalidValue is thrown if method s is not known
	*/
	void setModificationOutputMethod(const String& s);

	/**
	@brief getter for modification output method
	@return String
	*/
	String getModificationOutputMethod();

protected:

	String vToString_(std::vector<String> v);

	BigString big_string_; 	///< bigString object holding all peptides of fasta file

	SuffixArray* sa_; 	///< pointer to suffixarray

	String modification_output_method_; ///< output method for modifications

};

}

#endif //OPENMS_EXAMPLES_SuffixArrayPeptideFinder_H
