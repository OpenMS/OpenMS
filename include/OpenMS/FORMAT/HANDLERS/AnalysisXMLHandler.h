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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_ANALYSISXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_ANALYSISXMLHANDLER_H

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
    @brief Handler that is used for parsing AnalysisXML data
    
    @todo do not work with that many pointers internally. Replace by references (Nico)
  */
  class AnalysisXMLHandler:
    public XMLHandler
  {
    public:
      /// Constructor for loading
      AnalysisXMLHandler(std::vector<ProteinIdentification>& protein_identifications, std::vector<IdentificationData>& id_data, const String& filename);
      /// Constructor for loading
      AnalysisXMLHandler(std::vector<ProteinIdentification>& protein_identifications, std::vector<IdentificationData>& id_data, std::map<String, DoubleReal>& predicted_retention_times, const String& filename);
      /// Constructor for storing
      AnalysisXMLHandler(const std::vector<ProteinIdentification>& protein_identifications, const std::vector<IdentificationData>& id_data, const String& filename);
      /// Constructor for storing
      AnalysisXMLHandler(const std::vector<ProteinIdentification>& protein_identifications, const std::vector<IdentificationData>& id_data, const std::map<String, DoubleReal>& predicted_retention_times, const String& filename);
      
      /// Destructor
      ~AnalysisXMLHandler();
      
      /// Writes the xml file to the ostream 'os'
      void writeTo(std::ostream& os);
      	
			// Docu in base class
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
      virtual void characters(const XMLCh* const chars, unsigned int /*length*/);


    protected:
      std::vector<ProteinIdentification>* protein_identifications_;
      std::vector<IdentificationData>* id_data_;
      ProteinHit actual_protein_hit_;
      std::vector<ProteinHit> actual_protein_hits_;
      PeptideHit actual_peptide_hit_;
      std::vector<PeptideHit> actual_peptide_hits_;
			UInt peptide_identification_index_;
			UInt protein_identification_index_;
			bool inside_peptide_;   	
      const std::vector<ProteinIdentification> const_protein_identifications_;
      const std::vector<IdentificationData> const_id_data_;
			const std::map<String, DoubleReal> const_predicted_retention_times_;
			String tag_;      	
			UInt charge_identification_index_;
			bool inside_protein_;
			bool inside_global_protein_;
			std::vector<UInt> actual_peptide_indices_;
			std::map<String, DoubleReal>* predicted_retention_times_;
			std::vector< String > date_times_temp_;
			UInt date_times_counter_;
			String actual_date_time_;
			
		private:
				/// determines the date group index
		    UInt getDateGroupIndex(DateTime 												date_time,
			  															std::map< String , UInt> date_times);
			  																
				/// writes a peptide to the ostream 'os'			  																
				void writePeptideHit(std::ostream& os, 
		 													String shift,
		 													PeptideHit hit,
		 													Real significance_threshold,
		 													UInt identification_index,
		 													Int charge, 
		 													Real precursor_retention_time,
		 													Real precursor_mz,
		 													DateTime date_time,
					  									std::map< String , UInt> date_times,
														  DoubleReal predicted_rt_p_value = -1);


  };

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_ANALYSISXMLHANDLER_H
