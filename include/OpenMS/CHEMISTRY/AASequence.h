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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_AASEQUENCE_H
#define OPENMS_CHEMISTRY_AASEQUENCE_H

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/HashMap.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/Residue.h>

#include <vector>
#include <iostream>

namespace OpenMS
{
	class ResidueDB;

	/** 
		@brief Representation of a peptide/protein sequence
		
		@ingroup Chemistry
	*/
	class AASequence
	{
		public:
			
			/** @name Typedefs
			*/
			//@{
			typedef std::vector<const Residue*>::iterator iterator;
			typedef std::vector<const Residue*>::const_iterator const_iterator;
			typedef std::vector<const Residue*>::iterator Iterator;
			typedef std::vector<const Residue*>::const_iterator ConstIterator;
			//@}

			/** @name Constructors and Destructors
			*/
			//@{
			/// default constructor
			AASequence();
	
			/// copy constructor
			AASequence(const AASequence& rhs);

			/// copy constructor from a String
			AASequence(const String& rhs) throw(Exception::ParseError);

			/// constructor with given residue db pointer
			AASequence(ResidueDB* res_db);
			
			/// constructor with given a residue range
			AASequence(ConstIterator begin, ConstIterator end);
			
			/// destructor
			virtual ~AASequence();
			//@}

			/// assignment operator
			AASequence& operator = (const AASequence& rhs);
			
			/** @name Accessors
			*/
			//@{
			/// returns a pointer to the residue, which is at position index
			const Residue* getResidue(SignedInt index) const throw(Exception::IndexUnderflow, Exception::IndexOverflow);
			
			/// returns a pointer to the residue, which is at position index
			const Residue* getResidue(UnsignedInt index) const throw(Exception::IndexOverflow);
			
			/// returns the formula of the peptide
			EmpiricalFormula getFormula(Residue::ResidueType type = Residue::Full, SignedInt charge = 0) const;

			/// returns the average weight of the peptide
			Real getAverageWeight(Residue::ResidueType type = Residue::Full, SignedInt charge = 0) const;

			/// returns the mono isotopic weight of the peptide
			Real getMonoWeight(Residue::ResidueType type = Residue::Full, SignedInt charge = 0) const;

			/// fills in the map the neutral loss formulas associated with their occuring frequency
			HashMap<const EmpiricalFormula*, Size> getNeutralLosses() const;

			/// returns a pointer to the residue at given position
			const Residue* operator [] (SignedInt index) const throw(Exception::IndexUnderflow, Exception::IndexOverflow);
			
			/// returns a pointer to the residue at given position
			const Residue* operator [] (UnsignedInt index) const throw(Exception::IndexOverflow);
			
			/// adds the residues of the peptide
			AASequence operator + (const AASequence& peptide) const;

			/// adds the residues of the peptide, which is given as a string
			AASequence operator + (const String& peptide) const throw(Exception::ParseError);

			/// adds the residues of a peptide
			AASequence& operator += (const AASequence&);

			/// adds the residues of a peptide, which is given as a string
			AASequence& operator += (const String&) throw(Exception::ParseError);

			/** sets the residue db from an residue db; ATTENTION this affects all instances!
			 * 	calling with no argument resets to the default residues db usage
			 */
			void setResidueDB(ResidueDB* res_db = 0);

			/// returns the number of residues
			Size size() const;

			/// returns a peptide sequence of the first index residues
			AASequence getPrefix(Size index) const throw(Exception::IndexOverflow);

			/// returns a peptide sequence of the last index residues
			AASequence getSuffix(Size index) const throw(Exception::IndexOverflow);

			/// returns a peptide sequence of number residues, beginning at position index
			AASequence getSubsequence(Size index, Size number) const throw(Exception::IndexOverflow);
			//@}

			/** @name Predicates
			*/
			//@{
			/// returns true if the peptude contains the given residue
			bool has(const Residue* residue) const;

			/// returns true if the peptide contains the given residue
			bool has(const String& name) const;
			
			/// returns true if the peptide contains the given peptide
			bool hasSubsequence(const AASequence& peptide) const;

			// returns true if the peptide contains the given peptide
			bool hasSubsequence(const String& peptide) const throw(Exception::ParseError);
			
			/// returns true if the peptide has the given prefix
			bool hasPrefix(const AASequence& peptide) const;
			
			/// returns true if the peptide has the given prefix
			bool hasPrefix(const String& peptide) const throw(Exception::ParseError);

			/// returns true if the peptide has the given suffix
			bool hasSuffix(const AASequence& peptide) const;
			
			/// returns true if the peptide has the given suffix
			bool hasSuffix(const String& peptide) const throw(Exception::ParseError);

			/// equality operator
			bool operator == (const AASequence&) const;

			/// equality operator given the peptide as a string
			bool operator == (const String&) const throw(Exception::ParseError);

			/// inequality operator 
			bool operator != (const AASequence&) const;

			/// inequality operator given the peptide as a string
			bool operator != (const String&) const throw(Exception::ParseError);
			//@}

			/** @name Iterators
			*/
			//@{
			inline Iterator begin() { return peptide_.begin(); }

			inline ConstIterator begin() const { return peptide_.begin(); }

			inline Iterator end() { return peptide_.end(); }

			inline ConstIterator end() const { return peptide_.end(); }
			//@}
			
			/// writes a peptide to an output stream
			friend std::ostream& operator << (std::ostream& os, const AASequence& peptide);
			
			/// reads a peptide from an input stream
			friend std::istream& operator >> (std::istream& is, const AASequence& peptide);
			
		protected:
	
			ResidueDB* getResidueDB_() const;
			
			static ResidueDB* custom_res_db_;
			
			std::vector<const Residue*> peptide_;

			void parseString_(std::vector<const Residue*>& sequence, const String& peptide) const throw(Exception::ParseError);

	};			

	std::ostream& operator << (std::ostream& os, const AASequence& peptide);

	std::istream& operator >> (std::istream& os, const AASequence& peptide);
	
} // namespace OpenMS

#endif
