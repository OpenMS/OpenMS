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

#ifndef OPENMS_FORMAT_HANDLERS_XTANDEMXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_XTANDEMXMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>
#include <map>
#include <set>

namespace OpenMS
{
	namespace Internal
	{
  /**
    @brief Handler that is used for parsing XTandemXML data
    
  */
  class XTandemXMLHandler:
    public XMLHandler
  {
    public:
      /// Default constructor
      XTandemXMLHandler(ProteinIdentification& protein_identification, std::map<UInt, std::vector<PeptideHit> >& peptide_hits, const String& filename);

      /// Destructor
      virtual ~XTandemXMLHandler();
      
			// Docu in base class
      void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
      void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
		
			// Docu in base class
   		void characters(const XMLCh* const chars, const unsigned int /*length*/);
		  
    private:
    	
      ProteinIdentification& protein_id_;
			std::map<UInt, std::vector<PeptideHit> >& peptide_hits_;
			String actual_protein_id_;
			std::set<String> accessions_;
			Int actual_charge_;
			String tag_;
  };

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_XTANDEMXMLHANDLER_H
