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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------



#ifndef OPENMS_DATASTRUCTURES_SUFFIXARRAYTRYPTICSEQAN_H
#define OPENMS_DATASTRUCTURES_SUFFIXARRAYTRYPTICSEQAN_H

#include <OpenMS/DATASTRUCTURES/SuffixArraySeqan.h>

namespace OpenMS {

/**
	@brief Class that uses SEQAN library for a suffix array. It can be used to find peptide Candidates for a MS spectrum

	This class uses SEQAN suffix array. It can just be used for finding peptide Candidates for a given MS Spectrum within a certain mass tolerance. The suffix array can be saved to disc for reused so it has to be build just once.

*/

class OPENMS_DLLAPI SuffixArrayTrypticSeqan
	: public SuffixArraySeqan
{
	
public:
	
	/** @brief constructor for tryptic seqan array with a specially optimized implementation
	
			@param st the suffix array string, which is used to build the suffix array
			@param filename filename of fasta file

			@throw InvalidValue is thrown if string st if invalid
			@throw FileNotFound is thrown if given file is not found
	*/
	SuffixArrayTrypticSeqan(const String& st, const String& filename, const WeightWrapper::WEIGHTMODE weight_mode=WeightWrapper::MONO);
	
	/**
	@brief returns if an enzyme will cut after first character
	@param aa1 const char as first aminoacid
	@param aa2 const char as second aminoacid
	@return bool descibing if it is a digesting site
	*/
	bool isDigestingEnd(const char aa1, const char aa2) const;

	
};
}

#endif //OPENMS_DATASTRUCTURES_SUFFIXARRAYTRYPTICSEQAN_H
