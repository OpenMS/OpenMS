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
// $Authors: Ole Schulz-Trieglaff, Alexander Haupt$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_LCMSSAMPLE_H
#define OPENMS_SIMULATION_LCMSSAMPLE_H

#include <iostream>
#include <map>
#include <cmath>
#include <sstream>

#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>



#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>



namespace OpenMS
{

	/**
		@brief Representation of sample data, i.e. proteins and peptides.

		Reads the list of proteins from a FASTA file, digests and removes peptides
		with low detectability.

		Digestion parameters include <i>number of missed cleavages</i>,
		<i>minimum peptide length</i> and <i>maximum number of missed cleavages</i>.
	*/

	class OPENMS_DLLAPI LCMSSample
		: public DefaultParamHandler
	{

	public:

		/** @name Typedefs
		*/

		typedef std::map<OpenMS::String, unsigned int> PeptideSequences;

		typedef std::vector< std::pair<OpenMS::String,int> > SampleProteins;

		//@{
		/// Iterators
		typedef PeptideSequences::iterator Iterator;
		typedef PeptideSequences::iterator iterator;

		typedef PeptideSequences::const_iterator ConstIterator;
		typedef PeptideSequences::const_iterator const_iterator;

		typedef SampleProteins::const_iterator SampleProteinsConstIt;
		typedef SampleProteins::iterator SampleProteinsIt;
		//@}


		/** @name Constructors and Destructors
			*/
		//@{
		/// default constructor
		LCMSSample();

		/// Copy constructor
		LCMSSample(const LCMSSample& source);

		/// destructor
		virtual	~LCMSSample();

		LCMSSample& operator = (const LCMSSample& source);

		/// Load Proteins from FASTA file
		void loadFASTA(const String filename);
		//@}


		/** @name Accessors
			*/
		//@{
		/// Mutable accessor for peptides
		const PeptideSequences& getPeptideSequences() const { return peptides_; }
		//@}

		/// Digest proteins
		void digest();

		/// Print all sample proteins (for debugging)
		void printProteins() const;

		/// Print all sample peptides (for debugging)
		void printPeptides() const;

		/// Clear all sample proteins (for debugging)
		void clearProteins() { proteins_.clear(); }

		/// Set file name of SVM model 	for detectability prediction
		void setPdModelFile(OpenMS::String file) { DtModelFile_ = file; }

		/** @name Iterators
		*/
		//@{
		/// Iterators for accessing digested peptides
		inline Iterator begin() { return peptides_.begin(); }
		inline Iterator end()   { return peptides_.end();   }

		inline ConstIterator begin() const { return peptides_.begin(); }
		inline ConstIterator end()   const { return peptides_.end();   }

		inline size_t size() const { return peptides_.size(); }
		//@}

	private:

		/// Synchronize members with param class
		void updateMembers_();

		/// Filters peptides for detectability
		void filterForDetectability_(const std::vector<String>& all_peptides, std::vector<String>& filtered_peptides, UInt k_mer_length);

		/// Peptides
		PeptideSequences peptides_;
		/// Proteins
		SampleProteins proteins_;

		/// Minimum allowed detectability likelihood of a peptide
		double min_detect_;

		/// The SVM model for peptide detectability prediction
		String DtModelFile_;

		/// The support vector machine
		SVMWrapper svm_;

	}; // class definition

} // namespace OpenMS

#endif // OPENMS_SIMULATION_LCMSSAMPLE_H
