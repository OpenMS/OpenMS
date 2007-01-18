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
//
#ifndef OPENMS_CHEMISTRY_EMPIRICALFORMULA_H
#define OPENMS_CHEMISTRY_EMPIRICALFORMULA_H

#include <iostream>
#include <vector>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/HashMap.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
	class Element;
	class ElementDB;
	/** 
			@ingroup Chemistry
	
			@brief Representation of an empirical formula
	 	
		 	The formula can be used as follows: elements are 
		 	represented through its symbol or full name. The
		 	symbol or name is followed by a number. If not, the
		 	frequency is set to one. Examples are CH3OH or CarbonHydrogen3OH.
		 	The names must start with an capital letter (symbols always have
		 	an upper-case letter at the beginning). Additionally charges can
			be used with '+' or '-' followed by a number, if no number follows
		 	the charge of +/- 1 is set.
	*/
	
	class EmpiricalFormula
	{
	
		public:

			/** @name Typedefs
			*/
			//@{
			/// Iterators
			typedef HashMap<const Element*, Size>::ConstIterator ConstIterator;
			typedef HashMap<const Element*, Size>::ConstIterator const_iterator;

			/** @name Constructors and Destructors
			*/
			//@{
			/// default constructor
			EmpiricalFormula();
	
			/// copy constructor
			EmpiricalFormula(const EmpiricalFormula& rhs);
	
			/// constructor from an OpenMS String
			EmpiricalFormula(const String& rhs) throw(Exception::ParseError);

			/// constructor with element pointer and number
			EmpiricalFormula(Size number, const Element* element, SignedInt charge = 0);
			
			/// destructor
			virtual ~EmpiricalFormula();
			//@}

			/** @name Accessors
			*/
			//@{
			/// returns the mono isotopic weight of the formula
			Real getMonoWeight() const;
	
			/// returns the average weight of the formula
			Real getAverageWeight() const;

			/** @brief returns the isotope distribution of the formula
				*	The details of the calculation of the isotope distribution
				* are described in the doc to the IsotopeDistribution class.
				*	@param max_depth: this parameter gives the max isotope which is considered, if 0 all are reported
				*/
			IsotopeDistribution getIsotopeDistribution(Size max_depth = 20) const;
		
			/// sets the element db, the elements are read from the given file
			void setElementDB(const String& file_name) throw(Exception::FileNotFound, Exception::ParseError);

			/// returns a pointer to the element with name or symbol or 0 if no such element is fount
			const Element* getElement(const String& name) const;

			/// returns a pointer to the element with given atomic number or 0 if none if found
			const Element* getElement(Size atomic_number) const;
			
			/// returns a pointer to the element db which is used with this class
			const ElementDB* getElementDB() const;

			/// returns the number of atoms with the given atomic_number
			Size getNumberOf(Size atomic_number) const;

			/// returns the number of atoms with the given name
			Size getNumberOf(const String& name) const;

			/// returns the number of atoms 
			Size getNumberOf(const Element* element) const;

			/// returns the atoms total
			Size getNumberOfAtoms() const;

			/// returns the charge
			SignedInt getCharge() const;

			/// sets the charge
			void setCharge(SignedInt charge);

			/// returns the formula as a string
			String getString() const;
			//@}
			
			/** Assignment
			*/
			//@{
			/// assignment operator 
			EmpiricalFormula& operator = (const EmpiricalFormula& rhs);
	
			/// assignment operator which assigns an string to the formula
			EmpiricalFormula& operator = (const String& rhs) throw(Exception::ParseError);
	
			/// adds the elements of the given formula 
			EmpiricalFormula& operator += (const EmpiricalFormula& rhs);

			/// adds the elements from the given formula, which is given as a OpenMS String
			EmpiricalFormula& operator += (const String& rhs) throw(Exception::ParseError);
	
			/// adds the elements of the given formula and returns a new formula
			EmpiricalFormula operator + (const EmpiricalFormula& rhs) const;

			/// adds the elements of the given formula (given as a String) and returns a new formula
			EmpiricalFormula operator + (const String& rhs) const throw(Exception::ParseError);

			/// subtracts the elements of a formula
			EmpiricalFormula& operator -= (const EmpiricalFormula& rhs) throw(Exception::SizeUnderflow);

			/// subtracts the elements of a formula given as string
			EmpiricalFormula& operator -= (const String& rhs) throw(Exception::ParseError, Exception::SizeUnderflow);

			/// subtracts the elements of a formula an returns a new formula
			EmpiricalFormula operator - (const EmpiricalFormula& rhs) const throw(Exception::SizeUnderflow);

			/// subtracts the elements of a formula given as a String and returns a new formula
			EmpiricalFormula operator - (const String& rhs) const throw(Exception::ParseError, Exception::SizeUnderflow);
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
			bool hasElement(Size atomic_number) const;
			
			/// returns true if the formulas contain equal elements in equal quantities
			bool operator == (const EmpiricalFormula& rhs) const;

			/// returns true if the formulas contain equal elements in equal quantities
			bool operator == (const String& rhs) const throw(Exception::ParseError);

			/// returns true if the formulas differ in elements composition
			bool operator != (const EmpiricalFormula& rhs) const;
			
			/// returns true if the formulas differ in elements composition
			bool operator != (const String& rhs) const throw(Exception::ParseError);
			//@}

			/// writes the formula to a stream
			friend std::ostream& operator << (std::ostream&, const EmpiricalFormula&);
			
			/** @name Iterators
			*/
			//@{
			inline ConstIterator begin() const { return formula_.begin(); }

			inline ConstIterator end() const { return formula_.end(); }
			//@}
		
		protected:
		
			HashMap<const Element*, Size> formula_;

			SignedInt charge_;

			void readElementsFromFile_(const String& file_name) throw(Exception::FileNotFound, Exception::ParseError);

			SignedInt parseFormula_(HashMap<const Element*, Size>& ef,const String& formula) const throw(Exception::ParseError);
			
			const ElementDB* element_db_;
	};

	std::ostream& operator << (std::ostream&, const EmpiricalFormula::EmpiricalFormula&);

  /** 
      @ingroup Chemistry
  
      @brief some often used empirical formulas (which ist faster than creating instances of EmpiricalFormula from strings)

      To use i.e. an Hydrogen just use s.th. like Formulas::H
  */
	
	namespace Formulas
  {
    static const EmpiricalFormula& getH()
    {
      static const EmpiricalFormula H_("H");
      return H_;
    }

    static const EmpiricalFormula& getH2O()
    {
      static const EmpiricalFormula H2O_("H2O");
      return H2O_;
    }

    static const EmpiricalFormula& getNH()
    {
      static const EmpiricalFormula NH_("NH");
      return NH_;
    }

    static const EmpiricalFormula& getOH()
    {
      static const EmpiricalFormula OH_("OH");
      return OH_;
    }

    static const EmpiricalFormula& getNH3()
    {
      static const EmpiricalFormula NH3_("NH3");
      return NH3_;
    }

		const EmpiricalFormula H = getH();
		
		const EmpiricalFormula H2O = getH2O();
		
		const EmpiricalFormula NH = getNH();
		
		const EmpiricalFormula OH = getOH();
		
		const EmpiricalFormula NH3 = getNH3();
		
		const EmpiricalFormula Water = getH2O();

		const EmpiricalFormula Ammonia = getNH3();

		const EmpiricalFormula Hydrogen = getH();
  } // namespace Formulas

} // namespace OpenMS
#endif
