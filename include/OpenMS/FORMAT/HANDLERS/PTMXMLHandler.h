// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_PTMXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_PTMXMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <vector>
#include <map>
#include <fstream>

namespace OpenMS
{
	namespace Internal
	{
		/**
			@brief Handler that is used for parsing PTMXML data
		*/
		class OPENMS_DLLAPI PTMXMLHandler
			: public XMLHandler
		{
			public:
				/// Constructor for loading
				PTMXMLHandler(std::map< String, std::pair< String, String > >& ptm_informations, const String& filename);
				
				/// Destructor
				virtual ~PTMXMLHandler();
				
				/// Writes the xml file to the ostream 'os'
				void writeTo(std::ostream& os);
				
				// Docu in base class
				virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
				
				// Docu in base class
				virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
				
				// Docu in base class
				virtual void characters(const XMLCh* const chars, const XMLSize_t /*length*/);
				
			protected:
				std::map< String, std::pair< String, String > >& ptm_informations_;
				String name_, tag_, composition_;
				bool open_tag_;
		};
		
	} // namespace Internal
	
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_PTMXMLHANDLER_H
