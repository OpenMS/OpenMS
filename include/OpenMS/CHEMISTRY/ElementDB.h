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

#ifndef OPENMS_CHEMISTRY_ELEMENTDB_H
#define OPENMS_CHEMISTRY_ELEMENTDB_H

	#include <OpenMS/DATASTRUCTURES/String.h>

	#include <OpenMS/DATASTRUCTURES/HashMap.h>

	#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>

#define OPENMS_CHEMISTRY_ELEMENTDB_DEFAULT_FILE "CHEMISTRY/Elements.xml"

namespace OpenMS
{
	class Element;

	/** @ingroup Chemistry
	
			@brief Stores elements

	  	The elements weights (in the default file) are taken from
	  	"Isotopic Compositions of the Elements 1997", Pure Appl. Chem., 70(1), 217-235, 1998.
	  	(http://haven.isb-sib.ch/tools/isotopident/htdocs/iso.pdf)
	 
	  	The elements isotope distributions (in the default file) are taken from
	  	"Atomic Weights of the Elements 1999", Pure Appl. Chem. 73(4), 667-683, 2001.
	 		(http://haven.isb-sib.ch/tools/isotopident/htdocs/atomic-weights.pdf)
	 
	 		ATOMIC WEIGHTS OF THE ELEMENTS: REVIEW 2000 (IUPAC Technical Report)
	 		Pure Appl. Chem., Vol. 75, No. 6, pp. 683 800, 2003.
	*/
	
	class ElementDB
	{
		public:
	
			/** @name Accessors
			*/
			//@{
			/// returns a pointer to the singleton instance of the element db
			inline static const ElementDB* getInstance()
			{
				static ElementDB* db_ = 0;
				if (db_ == 0)
				{
					db_ = new ElementDB;
				}
				return db_;
			}
			
			/// returns a hashmap that contains names mapped to pointers to the elements
			const HashMap<String, const Element*>& getNames() const;
			
			/// returns a hashmap that contains symbols mapped to pointers to the elements
			const HashMap<String, const Element*>& getSymbols() const;

			/// returns a hashmap that contains atomic numbers mapped to pointers of the elements
			const HashMap<Size, const Element*>& getAtomicNumbers() const;

			/** returns a pointer to the element with name or symbol given in parameter name;
				*	if no element exists with that name or symbol 0 is returned
				*	@param name: name or symbol of the element
			*/
			const Element* getElement(const String& name) const;

			/// returns a pointer to the element of atomic number; if no element is found 0 is returned
			const Element* getElement(Size atomic_number) const;
			//@}
			
			/** @name Predicates
			*/
			//@{
			/// returns true if the db contains an element with the given name
			bool hasElement(const String& name) const;

			/// returns true if the db contains an element with the given atomic_number
			bool hasElement(Size atomic_number) const;
			//@}

		protected:
	
			/*_ parses a Histogram given as a OpenMS String and return the distribution
			 */
			IsotopeDistribution parseIsotopeDistribution_(const HashMap<UnsignedInt, double>& distribution) 
				throw(Exception::ParseError);

			/*_ read elements from a XML file, formated as a Param file.
			 */
			void readFromFile_(const String& file_name) throw(Exception::FileNotFound, Exception::ParseError);

			/*_ resets all containers 
			 */
			void clear_();

			HashMap<String, const Element*> names_;

			HashMap<String, const Element*> symbols_;

			HashMap<Size, const Element*> atomic_numbers_;

		private:

			ElementDB();

			ElementDB(const ElementDB& db);

			ElementDB& operator = (const ElementDB& db);

			virtual ~ElementDB();
	};

} // namespace OpenMS
#endif

