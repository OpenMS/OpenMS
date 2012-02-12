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
// $Maintainer: Nico Pfeifer $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MASCOTXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MASCOTXMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

namespace OpenMS
{
	namespace Internal
	{
  /**
    @brief Handler that is used for parsing MascotXML data
    
  */
  class OPENMS_DLLAPI MascotXMLHandler:
    public XMLHandler
  {
    public:
      /// Constructor
      MascotXMLHandler(ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& identifications, const String& filename, std::map<String, std::vector<AASequence> >& peptides);

      /// Destructor
      virtual ~MascotXMLHandler();
      
			// Docu in base class
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
      virtual void characters(const XMLCh* const chars, const XMLSize_t /*length*/);
		  
    private:
    	
    	ProteinIdentification& protein_identification_;	///< the protein identifications
      std::vector<PeptideIdentification>& id_data_;		///< the identifications (storing the peptide hits)
      ProteinHit actual_protein_hit_;												
      PeptideHit actual_peptide_hit_;
			UInt peptide_identification_index_;
			String tag_;
			DateTime date_;
			String date_time_string_;
			UInt actual_query_;
			ProteinIdentification::SearchParameters search_parameters_;
			String identifier_;
			String actual_title_;
			std::map<String, std::vector<AASequence> >& modified_peptides_;
      String warning_msg_;

			StringList tags_open_; ///< tracking the current XML tree
			String major_version_;
			String minor_version_;
  };

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_MASCOTXMLHANDLER_H
