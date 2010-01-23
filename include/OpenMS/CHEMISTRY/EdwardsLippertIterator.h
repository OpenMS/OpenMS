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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_EDWARDSLIPPERTITERATOR_H
#define OPENMS_CHEMISTRY_EDWARDSLIPPERTITERATOR_H

#include <OpenMS/CHEMISTRY/PepIterator.h>

namespace OpenMS
{
	
	/**
	@brief finds all Peptide Candidates with given masses and given fasta file
	
	The used algorithm was described by Edwards and Lippert. The algorithm uses two pointers to iterator over a sequence. One pointer represents the start the other the end.
	*/
	class OPENMS_DLLAPI EdwardsLippertIterator : public PepIterator
	{
		
		public:
		
			typedef std::pair<String,String> FASTAEntry;	
			/**
			@brief Constructor
			*/
			EdwardsLippertIterator();
			/**
			@brief Copy Constructor
			*/
			EdwardsLippertIterator (const EdwardsLippertIterator &);
			/**
			@brief Destructor
			*/
			virtual ~EdwardsLippertIterator ();
			
			/**
			@brief * operator for getting the value of the iterator
			@return FASTAEntry with specific candidate
			@throw InvalidIterator if iterator has not been initialized
			*/
			virtual FASTAEntry operator*();
			
			/**
			@brief opperator ++ for postincrement
			@return Reference to PepIterator
			@throw InvalidIterator if iterator has not been initialized
			*/
			virtual PepIterator & operator++();
			
			/**
			@brief opperator ++ for preincrement
			@return pointer to PepIterator
			@throw InvalidIterator if iterator has not been initialized
			*/
			virtual PepIterator * operator++(int i);
		
			/**
			@brief setter for fasta file
			@param f String with fasta file location
			@throw FileNotFound if file could not be found
			*/
			virtual void setFastaFile (const String & f);
			
			/**
			@brief getter for FASTA file
			@return String with file location
			*/
			virtual String getFastaFile();
		
			/**
			@brief setter for tolerance
			@param t tolerance
			@throw InvalidValue if tolerance is negative
			*/
			virtual void setTolerance (DoubleReal t);
			
			/**
			@brief getter for tolerance
			@return tolerance
			*/
			virtual DoubleReal getTolerance ();
		
			/**
			@brief setter for spectrum
			@param s spectrum as a vector of DoubleReals
			@throw InvalidValue if spectrum is not sorted
			*/
			virtual void setSpectrum(const std::vector<DoubleReal> & s);
		
			/**
			@brief getter for spectrum
			@return the used spectrum
			*/
			virtual const std::vector<DoubleReal> & getSpectrum ();
		
			/**
			@brief initializing iterator
			@return true if everything was ok
			@throw InvalidIterator is thrown if the begin iterator is invalid
			*/
			virtual bool begin ();
		
			/**
			@brief indicates whether iterator is at end
			@return true if iterator is at end
			@see hasNext
			*/
			virtual bool isAtEnd ();
			
			/**
			@brief indicated if a digesting enzyme will cut at this position
			@return true if digenting enzym cuts the sequence
			*/
			virtual bool isDigestingEnd (char, char);
			
			/**
			@brief needed by Factory
			@return const string name of class
			*/
			static const String getProductName()
			{
				return "EdwardsLippertIterator";
			}
		
			/**
			@brief needed by Factory
			@return poiter to new object
			*/
			static PepIterator* create()
			{
				return new EdwardsLippertIterator;
			}
		
		
		
		protected:
			/**
			@brief getting the next candidate
			@return string with next sequence
			*/
			virtual std::string next_ ();
		
			/**
			@brief indicates if there will be a next element
			@return true if iterator has more elements
			*/
			bool hasNext_();
		
			/**
			@brief finds the next starting position where a digesting enzyme will cut the sequence
			*/
			void goToNextAA_ ();
	
			/**
			@brief indicates if a mass is in spectrum
			@return true if a given mass is in spectrum
			*/
			virtual bool isInSpectrum_ (DoubleReal & mass);
	
			String f_file_; ///< fasta file location
	
			std::string actual_pep_; ///< actual peptide
	
			std::vector<DoubleReal> spec_; ///< given spectrum
	
			DoubleReal tol_; ///< tolerance
	
			DoubleReal masse_[255]; ///< mass table
	
			bool is_at_end_; ///< indicates if iterator is at end
	
			PepIterator * f_iterator_; ///< FastaIterator
	
			FASTAEntry f_entry_; ///< actual fasta entry
	
			unsigned int b_,e_; ///< to ints representing a position within the actual string (b = begin, e = end)
	
			DoubleReal m_, massMax_; ///< mass and maximum masse
		
		
	};

}
#endif //OPENMS_CHEMISTRY_EDWARDSLIPPERTITERATOR_H
