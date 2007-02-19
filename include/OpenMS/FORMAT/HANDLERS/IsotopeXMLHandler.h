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
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_IsotopeXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_IsotopeXMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/METADATA/Identification.h>

#include <vector>
#include <map>
#include <fstream>

namespace OpenMS
{
	namespace Internal
	{
		/**
			@brief Handler that is used for parsing IsotopeXML data
		*/
		class IsotopeXMLHandler:
			public XMLHandler
		{
			public:
				/// Constructor for loading
				IsotopeXMLHandler(std::map< String, std::vector< std::pair< DoubleReal, DoubleReal > > >& isotope_informations, const String& filename);
		
				/// Constructor for storing
				IsotopeXMLHandler(const std::map< String, std::vector< std::pair< DoubleReal, DoubleReal > > >& isotope_informations, const String& filename);
				
				/// Destructor
				~IsotopeXMLHandler();
				
				/// Writes the xml file to the ostream 'os'
				void writeTo(std::ostream& os); //throw (Exception::IndexOverflow);
				
				// Docu in base class
				virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
				
				// Docu in base class
				virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
				
				// Docu in base class
				virtual void characters(const XMLCh* const chars, const unsigned int /*length*/);
				
			protected:
				std::map< String, std::vector< std::pair< DoubleReal, DoubleReal > > >& isotope_informations_;
				const std::map< String, std::vector< std::pair< DoubleReal, DoubleReal > > >& const_isotope_informations_;
				String symbol_, tag_;
				bool open_tag_;
				DoubleReal mass_;
		};
		
	} // namespace Internal
	
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_IsotopeXMLHANDLER_H
