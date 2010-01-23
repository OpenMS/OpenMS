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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_XTANDEMINFILEXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_XTANDEMINFILEXMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

namespace OpenMS
{
	class XTandemInfile;
	namespace Internal
	{

		/** @brief Note representation of bioml structure used by XTandem


		*/
	  struct XTandemInfileNote
    {
      String note_type;
      String note_label;
      String note_value;
    };

  /**
    @brief Handler that is used for parsing XTandemXML data
    
  */
  class OPENMS_DLLAPI XTandemInfileXMLHandler:
    public XMLHandler
  {
    public:

      /// Default constructor
      XTandemInfileXMLHandler(const String& filename, std::vector<XTandemInfileNote>& notes, XTandemInfile* infile);

      /// Destructor
      virtual ~XTandemInfileXMLHandler();
      
			// Docu in base class
      void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
      void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
		
			// Docu in base class
   		void characters(const XMLCh* const chars, const XMLSize_t /*length*/);
		  
    protected:
	
			std::vector<XTandemInfileNote>& notes_;
	
			XTandemInfile* infile_;

			XTandemInfileNote actual_note_;

			String tag_;
  };

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_XTANDEMXMLHANDLER_H
