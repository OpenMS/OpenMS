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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_OMSSAXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_OMSSAXMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

namespace OpenMS
{
	namespace Internal
	{
  /**
    @brief Handler that is used for parsing OMSSAXML data
    
  */
  class OMSSAXMLHandler:
    public XMLHandler
  {
    public:

      /// Default constructor
      OMSSAXMLHandler(ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& peptide_identifications, const String& filename, bool load_proteins);

      /// Destructor
      virtual ~OMSSAXMLHandler();
      
			// Docu in base class
      void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
      void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
		
			// Docu in base class
   		void characters(const XMLCh* const chars, unsigned int /*length*/);
		  
    private:
    	
			/// the protein identifications
    	ProteinIdentification& protein_identification_;

			/// the identifications (storing the peptide hits)
      std::vector<PeptideIdentification>& peptide_identifications_;

      ProteinHit actual_protein_hit_;

      PeptideHit actual_peptide_hit_;

			PeptideIdentification actual_peptide_id_;

			ProteinIdentification actual_protein_id_;

			String tag_;

			/// site of the actual modification (simple position in the peptide)
			UInt actual_mod_site_;

			/// type of the modification
			String actual_mod_type_;

			/// modifications of the peptide defined by site and type
			std::vector<std::pair<UInt, String> > modifications_;

			/// should protein hits be read from the file?
			bool load_proteins_;
  };

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_OMSSAXMLHANDLER_H

