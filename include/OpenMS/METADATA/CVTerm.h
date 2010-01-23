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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_CVTERM_H
#define OPENMS_METADATA_CVTERM_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>

namespace OpenMS
{
	/**
		@brief Representation of controlled vocabulary term
		
		This class simply stores a CV term, its value and unit if necessary. 
		
		@ingroup Metadata
	*/
			///Represenation of a CV term used by CVMappings
			class OPENMS_DLLAPI CVTerm
			{
				public:
			

				struct Unit
				{
					Unit()
					{
					}

					Unit(const String& p_accession, const String& p_name, const String& p_cv_ref)
						: accession(p_accession),
							name(p_name),
							cv_ref(p_cv_ref)							
					{
					}

					Unit(const Unit& rhs)
						: accession(rhs.accession),
							name(rhs.name),
							cv_ref(rhs.cv_ref)
					{
					}

					virtual ~Unit()
					{
					}

					Unit& operator = (const Unit& rhs)
					{
						if (this != &rhs)
						{
							accession = rhs.accession;
							name = rhs.name;
							cv_ref = rhs.cv_ref;
						}
						return *this;
					}

					bool operator == (const Unit& rhs) const
					{
						return  accession == rhs.accession &&
										name == rhs.name &&
										cv_ref == rhs.cv_ref;
					}

					bool operator != (const Unit& rhs) const
					{
						return !(*this == rhs);
					}

					String accession;
					String name;
					String cv_ref;
				};


				/// Default constructor
				CVTerm();

				/// Detailed constructor
				CVTerm(const String& accession, const String& name, const String& cv_identifier_ref, const String& value, const Unit& unit);

				/// Copy constructor
				CVTerm(const CVTerm& rhs);

				/// Destructor
				virtual ~CVTerm();

				/// Assignment operator
				CVTerm& operator = (const CVTerm& rhs);

				/** @name Accessors
				*/
				//@{
				/// sets the accession string of the term
				void setAccession(const String& accession);

				/// returns the accession string of the term
				const String& getAccession() const;

				/// sets the name of the term
				void setName(const String& name);

				/// returns the name of the term
				const String& getName() const;

				/// sets the cv identifier reference string, e.g. UO for unit obo 
				void setCVIdentifierRef(const String& cv_identifier_ref);

				/// returns the cv identifier reference string
				const String& getCVIdentifierRef() const;

				/// set the value of the term
				void setValue(const DataValue& value);

				/// returns the value of the term
				const DataValue& getValue() const;

				/// sets the unit of the term
				void setUnit(const Unit& unit);

				/// returns the unit
				const Unit& getUnit() const;
				//@}
			
				/** @name Predicates
				*/
				//@{
				/// equality operator 
				bool operator == (const CVTerm& rhs) const;
				
				/// inequality operator 
				bool operator != (const CVTerm& rhs) const;

				/// checks whether the term has a value
				bool hasValue() const;

				/// checks whether the term has a unit
				bool hasUnit() const;
				//}

				protected:
				
				String accession_;

				String name_;

				String cv_identifier_ref_;

				Unit unit_;

				DataValue value_;
			};
						
} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CVTERM_H
