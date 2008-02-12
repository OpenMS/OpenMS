// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_SCHEMA2CV_H
#define OPENMS_FORMAT_SCHEMA2CV_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <vector>


namespace OpenMS
{
	namespace Internal
	{
		class Schema2CVHandler;
	}
	
	/**
		@brief A mapping of ControlledVocabulary terms to locations of an XML schema.
		
		@todo add check after reading the file -> ParseError (Marc)
		@todo add 'repeatable' (Marc)
		
		@experimental Perhaps used for mzML (Marc)
	*/
	class Schema2CV
	{
		friend class Internal::Schema2CVHandler;
		friend std::ostream& operator << (std::ostream& os, const Schema2CV& cv);
			
		public:
			
			///Helper struct for the CV definition
			struct CVDesc
			{
				String name;		///< Name of the CV
				String version;	///< Version of the CV
				String uri;			///< URL of the CV
				String id;			///< Identifier of the CV referenced by TermDesc
				String format;	///< Format of the CV (obo or owl)
			};
			
			///Helper struct for CV terms
			struct TermDesc
			{
				String accession;			///< Term accession
				String cv;						///< CV identifier
				bool allowSelf;				///< Indicates if the term itself is allowed
				bool allowChildren;	///< Indicates if child terms of the term are allowed
				bool repeatable;			///< Indicates if the term and its child terms can occur several times
			};
			
			///Helper struct for the location
			struct LocDesc
			{
				String location;									///< XPath location
				bool strict; 										///< Indicates if only the terms listed can be used
				std::vector<TermDesc> terms;	///< Allowed terms
			};
			
			/// Constructor
			Schema2CV();
			
			///Destructor
			~Schema2CV();
			
			///load the CV from a OBO file
			void loadFromFile(const String& filename) throw (Exception::FileNotFound, Exception::ParseError);
			
			///Returns the registered CVs
			const std::vector<CVDesc>& getCVs() const;

			///Returns the registered Paths
			const std::vector<LocDesc>& getLocations() const;
			
		protected:
			std::vector<CVDesc> cvs_;
			std::vector<LocDesc> locs_;
	};
	
	///Print the contents to a stream.
	std::ostream& operator << (std::ostream& os, const Schema2CV& mapping);
	
	namespace Internal
	{
		class Schema2CVHandler
			: public XMLHandler
		{
			public:
				/// Default constructor
				Schema2CVHandler(const String& filename, Schema2CV& mapping);
				
				/// Parsing method for opening tags
				virtual void startElement(const XMLCh* const uri, const XMLCh* const localname, const XMLCh* const qname, const xercesc::Attributes& attrs);
				
			protected:
				/// Reference to Schema2CV to fill
				Schema2CV& mapping_;

			private:
				///Not implemented
				Schema2CVHandler();
		};
	}
	
} // namespace OpenMS

#endif // OPENMS_FORMAT_SCHEMA2CV_H
