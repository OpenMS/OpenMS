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
#include <OpenMS/DATASTRUCTURES/Map.h>
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
			class Iterator;
						
			// iterator
			class ConstIterator
			{
				public: 
				
				// TODO Iterator constructor for ConstIterator
								
				typedef const Residue& 	const_reference;
				typedef Residue& 				reference;
				typedef const Residue* 	const_pointer;
				typedef std::vector<const Residue*>::difference_type difference_type;

				/** @name Constructors and destructors
				*/
				//@{
				/// default constructor
				ConstIterator()
				{
				}

				/// detailed constructor with pointer to the vector and offset position
				ConstIterator(const std::vector<const Residue*>* vec_ptr, int position)
				{
					vector_ = vec_ptr;
					position_ = position;
				}

				/// copy constructor
				ConstIterator(const ConstIterator& rhs)
					:	vector_(rhs.vector_),
						position_(rhs.position_)
				{
				}

				/// copy constructor from Iterator
				ConstIterator(const AASequence::Iterator& rhs)
					: vector_(rhs.vector_),
						position_(rhs.position_)
				{
				}

				/// destructor
				virtual ~ConstIterator()
				{
				}
				//@}

				/// assignment operator 
				ConstIterator& operator = (const ConstIterator& rhs)
				{
					if (this != &rhs)
					{
						position_ = rhs.position_;
						vector_ = rhs.vector_;
					}
					return *this;
				}

				/** @name Operators
				*/
				//@{
				/// dereference operator
				const_reference operator * () const
				{
					return *(*vector_)[position_];
				}

				/// dereference operator
				const_pointer operator -> () const
				{
					return (*vector_)[position_];
				}

				/// forward jump operator
				const ConstIterator operator + (difference_type diff) const
				{
					return ConstIterator(vector_, position_ + diff);
				}			

				difference_type operator - (ConstIterator rhs) const
				{
					return position_ - rhs.position_;
				}

				/// backward jump operator
				const ConstIterator operator - (difference_type diff) const
				{
					return ConstIterator(vector_, position_ - diff);
				}

				/// equality comparator
				bool operator == (const ConstIterator& rhs) const
				{
					return vector_ == rhs.vector_ && position_ == rhs.position_;
				}

				/// inequality operator
				bool operator != (const ConstIterator& rhs) const
				{
					return vector_ != rhs.vector_ || position_ != rhs.position_;
				}

				/// increment operator
				ConstIterator& operator ++ ()
				{
					++position_;
					return *this;
				}

				/// decrement operator
				ConstIterator& operator -- ()
				{
					--position_;
					return *this;
				}
				//@}

				protected:

				// pointer to the AASequence vector
				const std::vector<const Residue*>* vector_;

				// position in the AASequence vector
				difference_type position_;
			};

			class Iterator
			{
				public: 

				friend class AASequence::ConstIterator;
								
				typedef const Residue& const_reference;
				typedef Residue& reference;
				typedef const Residue* const_pointer;
				typedef const Residue* pointer;
				typedef std::vector<const Residue*>::difference_type difference_type;

				/** @name Constructors and destructors
				*/
				//@{
				/// default constructor
				Iterator()
				{
				}

				/// detailed constructor with pointer to the vector and offset position
				Iterator(std::vector<const Residue*>* vec_ptr, int position)
				{
					vector_ = vec_ptr;
					position_ = position;
				}

				/// copy constructor
				Iterator(const Iterator& rhs)
					:	vector_(rhs.vector_),
						position_(rhs.position_)
				{
				}

				/// destructor
				virtual ~Iterator()
				{
				}
				//@}

				/// assignment operator 
				Iterator& operator = (const Iterator& rhs)
				{
					if (this != &rhs)
					{
						position_ = rhs.position_;
						vector_ = rhs.vector_;
					}
					return *this;
				}

				/** @name Operators
				*/
				//@{
				/// dereference operator
				const_reference operator * () const
				{
					return *(*vector_)[position_];
				}

				/// dereference operator
				const_pointer operator -> () const
				{
					return (*vector_)[position_];
				}

				/// mutable dereference operator
				pointer operator -> ()
				{
					return (*vector_)[position_];
				}

				/// forward jump operator
				const Iterator operator + (difference_type diff) const
				{
					return Iterator(vector_, position_ + diff);
				}			

				difference_type operator - (Iterator rhs) const
				{
					return position_ - rhs.position_;
				}

				/// backward jump operator
				const Iterator operator - (difference_type diff) const
				{
					return Iterator(vector_, position_ - diff);
				}

				/// equality comparator
				bool operator == (const Iterator& rhs) const
				{
					return vector_ == rhs.vector_ && position_ == rhs.position_;
				}

				/// inequality operator
				bool operator != (const Iterator& rhs) const
				{
					return vector_ != rhs.vector_ || position_ != rhs.position_;
				}

				/// increment operator
				Iterator& operator ++ ()
				{
					++position_;
					return *this;
				}

				/// decrement operator
				Iterator& operator -- ()
				{
					--position_;
					return *this;
				}
				//@}

				protected:

				// pointer to the AASequence vector
				std::vector<const Residue*>* vector_;

				// position in the AASequence vector
				difference_type position_;
			};


			
			/** @name Constructors and Destructors
			*/
			//@{
			/// default constructor
			AASequence();
	
			/// copy constructor
			AASequence(const AASequence& rhs);

			/// copy constructor from a String
			AASequence(const String& rhs);

			/// copy consturctor from char* string
			AASequence(const char* rhs);

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
			/// returns the given string; if no string was given, the sequence is converted into a string
			String toString() const;

			String toUnmodifiedString() const;
			
			/// 
			void setModification(UInt index, const String& modification);

			void setNTerminalModification(const String& modification);

			const String& getNTerminalModification() const;
			
			void setCTerminalModification(const String& modification);

			const String& getCTerminalModification() const;
			
			/// sets the string of the sequence; returns true if the conversion to real AASequence was successful, false otherwise
			bool setStringSequence(const String& sequence);
			
			/// returns a pointer to the residue, which is at position index
			const Residue& getResidue(Int index) const;
			
			/// returns a pointer to the residue, which is at position index
			const Residue& getResidue(UInt index) const;
			
			/// returns the formula of the peptide
			EmpiricalFormula getFormula(Residue::ResidueType type = Residue::Full, Int charge = 0) const;

			/// returns the average weight of the peptide
			DoubleReal getAverageWeight(Residue::ResidueType type = Residue::Full, Int charge = 0) const;

			/// returns the mono isotopic weight of the peptide
			DoubleReal getMonoWeight(Residue::ResidueType type = Residue::Full, Int charge = 0) const;

			/// fills in the map the neutral loss formulas associated with their occuring frequency
			Map<const EmpiricalFormula*, UInt> getNeutralLosses() const;

			/// returns a pointer to the residue at given position
			const Residue& operator [] (Int index) const;
			
			/// returns a pointer to the residue at given position
			const Residue& operator [] (UInt index) const;
			
			/// adds the residues of the peptide
			AASequence operator + (const AASequence& peptide) const;

			/// adds the residues of the peptide, which is given as a string
			AASequence operator + (const String& peptide) const;

			AASequence operator + (const char* rhs) const;

			/// adds the residues of a peptide
			AASequence& operator += (const AASequence&);

			/// adds the residues of a peptide, which is given as a string
			AASequence& operator += (const String&);

			AASequence& operator += (const char* rhs);
			
			/// returns the number of residues
			UInt size() const;

			/// returns a peptide sequence of the first index residues
			AASequence getPrefix(UInt index) const;

			/// returns a peptide sequence of the last index residues
			AASequence getSuffix(UInt index) const;

			/// returns a peptide sequence of number residues, beginning at position index
			AASequence getSubsequence(UInt index, UInt number) const;
			//@}

			/** @name Predicates
			*/
			//@{
			/** @brief return true if the instance is valid
			 
			 		Valid means that a possible given sequence as string was successful
					converted into a real amino acid sequence which meaningful amino acids
					and modifications associated with it.
			*/
			bool isValid() const;
			
			/// returns true if the peptude contains the given residue
			bool has(const Residue& residue) const;

			/// returns true if the peptide contains the given residue
			bool has(const String& name) const;
			
			/// returns true if the peptide contains the given peptide
			bool hasSubsequence(const AASequence& peptide) const;

			/// returns true if the peptide contains the given peptide
			bool hasSubsequence(const String& peptide) const;
			
			/// returns true if the peptide has the given prefix
			bool hasPrefix(const AASequence& peptide) const;
			
			/// returns true if the peptide has the given prefix
			bool hasPrefix(const String& peptide) const;

			/// returns true if the peptide has the given suffix
			bool hasSuffix(const AASequence& peptide) const;
			
			/// returns true if the peptide has the given suffix
			bool hasSuffix(const String& peptide) const;

			bool hasNTerminalModification() const;

			bool hasCTerminalModification() const;
			
			// returns true if any of the residues is modified
			bool isModified() const;
			
			/// returns true if the residue at the position is modified
			bool isModified(UInt index) const;
			
			/// equality operator
			bool operator == (const AASequence&) const;

			/// equality operator given the peptide as a string
			bool operator == (const String&) const;

			/// 
			bool operator == (const char* rhs) const;

			bool operator < (const AASequence& rhs) const;

			bool operator <= (const AASequence& rhs) const;
			
			bool operator > (const AASequence& rhs) const;

			bool operator >= (const AASequence& rhs) const;
			
			/// inequality operator 
			bool operator != (const AASequence&) const;

			/// inequality operator given the peptide as a string
			bool operator != (const String&) const;

			bool operator != (const char* rhs) const;
			//@}

			/** @name Iterators
			*/
			//@{
			inline Iterator begin() { return Iterator(&peptide_, 0); }

			inline ConstIterator begin() const { return ConstIterator(&peptide_, 0); }

			inline Iterator end() { return Iterator(&peptide_, peptide_.size()); }

			inline ConstIterator end() const { return ConstIterator(&peptide_, peptide_.size()); }
			//@}
			
			/// writes a peptide to an output stream
			friend std::ostream& operator << (std::ostream& os, const AASequence& peptide);
			
			/// reads a peptide from an input stream
			friend std::istream& operator >> (std::istream& is, const AASequence& peptide);
			
		protected:

			std::vector<const Residue*> peptide_;

			String sequence_string_;
			
			void parseString_(std::vector<const Residue*>& sequence, const String& peptide);

			ResidueDB* getResidueDB_() const;

			bool valid_;

			String n_term_mod_;

			String c_term_mod_;
	};			

	std::ostream& operator << (std::ostream& os, const AASequence& peptide);

	std::istream& operator >> (std::istream& os, const AASequence& peptide);
	
} // namespace OpenMS

#endif
