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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_IDXMLFILE_H
#define OPENMS_FORMAT_IDXMLFILE_H

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief Used to load and store IdXML files
    
    This class is used to load and store documents that implement 
    the schema of IdXML files.
		
		A documented schema for this format can be found at http://open-ms.sourceforge.net/schemas/. 
		
		One file can contain several ProteinIdentification runs. Each run consists of peptide hits stored in 
		PeptideIdentification and (optional) protein hits stored in Identification. Peptide and protein
		hits are connected via a string identifier. We use the search engine and the date as identifier.
		
  	@note This format will eventually be replaced by the HUPO-PSI (mzIdentML and mzQuantML)) AnalysisXML formats!
  	
  	@ingroup FileIO
  */
  class OPENMS_DLLAPI IdXMLFile
  	: protected Internal::XMLHandler,
  		public Internal::XMLFile
  {
		public:
		
			/// Constructor
			IdXMLFile();
		


			/**
				@brief Loads the identifications of an IdXML file without identifier
				
				The information is read in and the information is stored in the
				corresponding variables

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
			void load(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids);

			/**
				@brief Loads the identifications of an IdXML file
				
				The information is read in and the information is stored in the
				corresponding variables

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
			void load(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids, String& document_id);
			 			 
			/**
				@brief Stores the data in an IdXML file
				
				The data is read in and stored in the file 'filename'.

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
			*/
			void store(String filename, const std::vector<ProteinIdentification>& protein_ids, const std::vector<PeptideIdentification>& peptide_ids ,const String& document_id=""); 
  	
  	protected:
			// Docu in base class
			virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
			virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			/// @name members for loading data
			//@{
			/// Pointer to fill in protein identifications
			std::vector<ProteinIdentification>* prot_ids_;
			/// Pointer to fill in peptide identifications
			std::vector<PeptideIdentification>* pep_ids_;
			/// Pointer to last read object with MetaInfoInterface
			MetaInfoInterface* last_meta_;
			/// Search parameters map (key is the "id")
			std::map<String,ProteinIdentification::SearchParameters> parameters_;
			/// Temporary search parameters variable
			ProteinIdentification::SearchParameters param_;
			/// Temporary id
			String id_;
			/// Temporary protein ProteinIdentification
			ProteinIdentification prot_id_;
			/// Temporary peptide ProteinIdentification
			PeptideIdentification pep_id_;
			/// Temporary protein hit
			ProteinHit prot_hit_;
			/// Temporary peptide hit
			PeptideHit pep_hit_;
			/// Map from protein id to accession
			std::map<String,String> proteinid_to_accession_;
			/// Document identitifier
			String* document_id_;
			/// true if a prot id is contained in the current run
			bool prot_id_in_run_;
			//@}
  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_IDXMLFILE_H
