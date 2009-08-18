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
// $Maintainer: Marc Sturm, Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_CONTROLLEDVOCABULARY_H
#define OPENMS_FORMAT_CONTROLLEDVOCABULARY_H

#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <set>

namespace OpenMS
{
	/**
		@brief Representation of a controlled vocabulary.
		
		This represenation only contains the information used for parsing and validation.
		All other lines are stored in the @em unparsed member of the the CVTerm struct.
  	
  	@ingroup Format
	*/
	class OPENMS_DLLAPI ControlledVocabulary
	{
		friend OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const ControlledVocabulary& cv);
			
		public:
			/// Representation of a CV term
			struct CVTerm
			{
				/// define xsd types allowed in cv term to specify their value-type
				enum XRefType
				{
					XSD_STRING=0, 						// xsd:string A string 
					XSD_INTEGER,  						// xsd:integer Any integer 
					XSD_DECIMAL,  						// xsd:decimal Any real number 
					XSD_NEGATIVE_INTEGER, 		// xsd:negativeInteger Any negative integer 
					XSD_POSITIVE_INTEGER, 		// xsd:positiveInteger Any integer > 0 
					XSD_NON_NEGATIVE_INTEGER, // xsd:nonNegativeInteger Any integer >= 0 
					XSD_NON_POSITIVE_INTEGER, // xsd:nonPositiveInteger Any integer < 0 
					XSD_BOOLEAN, 						 	// xsd:boolean True or false 
					XSD_DATE,								 	// xsd:date An XML-Schema date
					XSD_ANYURI,								// xsd:anyURI uniform resource identifier
					NONE
				};

				static String getXRefTypeName(XRefType type)
				{
					switch (type)
					{
						case XSD_STRING: return "xsd:string";
						case XSD_INTEGER: return "xsd:integer";
						case XSD_DECIMAL: return "xsd:decimal";
						case XSD_NEGATIVE_INTEGER: return "xsd:negativeInteger";
						case XSD_POSITIVE_INTEGER: return "xsd:positiveInteger";
						case XSD_NON_NEGATIVE_INTEGER: return "xsd:nonNegativeInteger";
						case XSD_NON_POSITIVE_INTEGER: return "xsd:nonPositiveInteger";
						case XSD_BOOLEAN: return "xsd:boolean";
						case XSD_DATE: return "xsd:date";
						case XSD_ANYURI: return "xsd:anyURI";
						default: return "none";
					}
					return "";
				}
							
				String name;									///< Text name
				String id;										///< Identifier
				std::set<String> parents;	    ///< The parent IDs
				std::set<String> children;    ///< The child IDs
				bool obsolete; 								///< Flag that indicates of the term is obsolete
				String description;       		///< Term description
				StringList synonyms;					///< List of synonyms
				StringList unparsed;					///< Unparsed lines from the definition file
				XRefType xref_type;						///< xref value-type for the CV-term
				StringList xref_binary;				///< xref binary-data-type for the CV-term (list of all allowed data value types for the current binary data array)
				std::set<String> units;       ///< unit accession ids, defined by relationship has units
				
				///Default constructor
				CVTerm()
					: name(),
						id(),
						parents(),
						children(),
						obsolete(false),
						description(),
						synonyms(),
						unparsed(),
						xref_type(NONE),
						xref_binary()
				{
				}

				CVTerm(const CVTerm& rhs)
					: name(rhs.name),
						id(rhs.id),
						parents(rhs.parents),
						children(rhs.children),
						obsolete(rhs.obsolete),
						description(rhs.description),
						synonyms(rhs.synonyms),
						unparsed(rhs.unparsed),
						xref_type(rhs.xref_type),
						xref_binary(rhs.xref_binary),
						units(rhs.units)
				{
				}

				CVTerm& operator = (const CVTerm& rhs)
				{
					if (this != &rhs)
					{
						name = rhs.name;
						id = rhs.id;
						parents = rhs.parents;
						children = rhs.children;
						obsolete = rhs.obsolete;
						description = rhs.description;
						synonyms = rhs.synonyms;
						unparsed = rhs.unparsed;
						xref_type = rhs.xref_type;
						xref_binary = rhs.xref_binary;
						units = rhs.units;
					}
					return *this;
				}
			};
			
			/// Constructor
			ControlledVocabulary();
			
			///Destructor
			virtual ~ControlledVocabulary();
			
			/// Returns the CV name (set in the load method)
			const String& name() const;
			
			/**
				@brief Loads the CV from an OBO file

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
			void loadFromOBO(const String& name, const String& filename);
			
			/// Returns true if the term is in the CV. Returns false otherwise.
			bool exists(const String& id) const;
			
			/**
				@brief Returns a term specified by ID
				
				@exception Exception::InvalidValue is thrown if the term is not present
			*/
			const CVTerm& getTerm(const String& id) const;


			/// returns all the terms stored in the CV
			const Map<String, CVTerm>& getTerms() const;

			/**
				@brief Writes all child terms recursively into terms

				If parent has child this method writes them recursively into the term object
				
				@exception Exception::InvalidValue is thrown if the term is not present
			*/
			void getAllChildTerms(std::set<String>& terms, const String& parent) const;
			
			/**
				@brief Returns if @p child is a child of @p parent
				
				@exception Exception::InvalidValue is thrown if one of the terms is not present
			*/
			bool isChildOf(const String& child, const String& parent) const;
			
		protected:
			/**
				@brief checks if a name corresponds to an id
				
				If the term is not known, 'true' is returned!				
			*/
			bool checkName_(const String& id, const String& name, bool ignore_case=true);
			
			///Map from ID to CVTerm
			Map<String, CVTerm> terms_;
			///Name set in the load method
			String name_;			
	};
	
	///Print the contents to a stream.
	OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const ControlledVocabulary& cv);
	
	
} // namespace OpenMS

#endif // OPENMS_FORMAT_CONTROLLEDVOCABULARY_H
