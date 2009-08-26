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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_CHEMISTRY_EMPIRICALFORMULA_H
#define OPENMS_CHEMISTRY_EMPIRICALFORMULA_H

#include <iostream>
#include <vector>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
	class Element;
	class ElementDB;
	/**
			@ingroup Chemistry

			@brief Representation of an empirical formula

		 	The formula can be used as follows: elements are represented through its symbol or full name. 
			The symbol or name is followed by a number. If not, the frequency is set to one. Examples 
			are CH3OH or CarbonHydrogen3OH. The names must start with an capital letter (symbols always have
		 	an upper-case letter at the beginning). Additionally charges can be used with '+' followed by 
			a number, if no number follows the charge of +1 is set.

			This class also supports the usage of specific isotopes. By default "C" describes not one isotope
			but a natural distribution or different isotopes. This distribution can be accessed via the 
			member getIsotopeDistribution(). 

			If one wants only use a specific isotope, it can be specified using "(",")" brackets. For example,
			to specify 14C a heavy isotope of carbon it is expressed as "(14)C". The isotope distribution 
			of that instance contains only one isotope, 14C itself with a frequency of 100%. 

			Instances EmpiricalFormula support a (limited) set of mathmatical operations. Additions and subtractions
			are supported in different flavors. However, one must be careful, because this can lead to negative 
			frequencies. In most cases this might be misleading, however, the class therefore supports difference
			formulae. E.g. formula differences of reactions from post-translational modifications.
	*/

	class OPENMS_DLLAPI EmpiricalFormula
	{

		public:

			/** @name Typedefs
			*/
			//@{
			/// Iterators
			typedef Map<const Element*, SignedSize>::ConstIterator ConstIterator;
			typedef Map<const Element*, SignedSize>::ConstIterator const_iterator;

			/** @name Constructors and Destructors
			*/
			//@{
			/// default constructor
			EmpiricalFormula();

			/// copy constructor
			EmpiricalFormula(const EmpiricalFormula& rhs);

			/** constructor from an OpenMS String

					@throw throws ParseError if the formula cannot be parsed
			*/
			EmpiricalFormula(const String& rhs);

			/// constructor with element pointer and number
			EmpiricalFormula(SignedSize number, const Element* element, Size charge = 0);

			/// destructor
			virtual ~EmpiricalFormula();
			//@}

			/** @name Accessors
			*/
			//@{
			/// returns the mono isotopic weight of the formula
			DoubleReal getMonoWeight() const;

			/// returns the average weight of the formula
			DoubleReal getAverageWeight() const;

			/** @brief returns the isotope distribution of the formula
				*	The details of the calculation of the isotope distribution
				* are described in the doc to the IsotopeDistribution class.
				*	@param max_depth: this parameter gives the max isotope which is considered, if 0 all are reported
				*/
			IsotopeDistribution getIsotopeDistribution(UInt max_depth) const;

			/// returns a pointer to the element with name or symbol or 0 if no such element is fount
			const Element* getElement(const String& name) const;

			/// returns a pointer to the element with given atomic number or 0 if none if found
			const Element* getElement(UInt atomic_number) const;

			/// returns a pointer to the element db which is used with this class
			const ElementDB* getElementDB() const;

			/// returns the number of atoms with the given atomic_number
			Size getNumberOf(UInt atomic_number) const;

			/// returns the number of atoms with the given name
			Size getNumberOf(const String& name) const;

			/// returns the number of atoms
			Size getNumberOf(const Element* element) const;

			/// returns the atoms total
			Size getNumberOfAtoms() const;

			/// returns the charge
			Size getCharge() const;

			/// sets the charge
			void setCharge(Size charge);

			/// returns the formula as a string
			String getString() const;
			//@}

			/** Assignment
			*/
			//@{
			/// assignment operator
			EmpiricalFormula& operator = (const EmpiricalFormula& rhs);

			/** assignment operator which assigns an string to the formula

					@throw throws ParseError if the formula cannot be parsed
			*/
			EmpiricalFormula& operator = (const String& rhs);

			/// adds the elements of the given formula
			EmpiricalFormula& operator += (const EmpiricalFormula& rhs);

			/** adds the elements from the given formula, which is given as a OpenMS String
			
					@throw throws ParseError if the formula cannot be parsed
			*/
			EmpiricalFormula& operator += (const String& rhs);

			/// multiplies the elements and charge with a factor
			EmpiricalFormula operator * (const SignedSize& times) const;
			
			/// adds the elements of the given formula and returns a new formula
			EmpiricalFormula operator + (const EmpiricalFormula& rhs) const;

			/** adds the elements of the given formula (given as a String) and returns a new formula

					@throw throws ParseError if the formula cannot be parsed
			*/
			EmpiricalFormula operator + (const String& rhs) const;

			/** subtracts the elements of a formula

					@throw throws SizeUnderflow if one number of elements of right hand side is larger than left hand side
			*/
			EmpiricalFormula& operator -= (const EmpiricalFormula& rhs);

			/** subtracts the elements of a formula given as string
			
					@throw throws ParseError if the formula cannot be parsed
					@throw throws SizeUnderflow if one number of elements of right hand side is larger than left hand side
			*/
			EmpiricalFormula& operator -= (const String& rhs);

			/** subtracts the elements of a formula an returns a new formula

					@throw throws SizeUnderflow if one number of elements of right hand side is larger than left hand side
			*/
			EmpiricalFormula operator - (const EmpiricalFormula& rhs) const;

			/** subtracts the elements of a formula given as a String and returns a new formula

					@throw throws ParseError if the formula cannot be parsed
					@throw throws SizeUnderflow if one number of elements of right hand side is larger than left hand side
			*/
			EmpiricalFormula operator - (const String& rhs) const;
			//@}

			/**@name Predicates
			*/
			//@{
			/// returns true if the formula does not contain a element
			bool isEmpty() const;

			/// returns true if charge != 0
			bool isCharged() const;

			/// returns true if the formula contains the element
			bool hasElement(const Element* element) const;

			/// returns true if the formula contains the element, given with its name or symbol
			bool hasElement(const String& name) const;

			/// returns true if the formula contains the element with the given atomic number
			bool hasElement(UInt atomic_number) const;

			/// returns true if the formulas contain equal elements in equal quantities
			bool operator == (const EmpiricalFormula& rhs) const;

			/** returns true if the formulas contain equal elements in equal quantities
					
					@throw throws ParseError if the formula cannot be parsed
			*/
			bool operator == (const String& rhs) const;

			/// returns true if the formulas differ in elements composition
			bool operator != (const EmpiricalFormula& rhs) const;

			/** returns true if the formulas differ in elements composition

					@throw throws ParseError if the formula cannot be parsed
			*/
			bool operator != (const String& rhs) const;
			//@}

			/// writes the formula to a stream
			friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const EmpiricalFormula& formula);

			/** @name Iterators
			*/
			//@{
			inline ConstIterator begin() const { return formula_.begin(); }

			inline ConstIterator end() const { return formula_.end(); }
			//@}

		protected:

			Map<const Element*, SignedSize> formula_;

			Size charge_;

			void readElementsFromFile_(const String& file_name);

			Size parseFormula_(Map<const Element*, SignedSize>& ef,const String& formula) const;

			const ElementDB* element_db_;
	};

	OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const EmpiricalFormula& formula);

} // namespace OpenMS
#endif
