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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PROTXMLFILE_H
#define OPENMS_FORMAT_PROTXMLFILE_H

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief Used to load (storing not supported, yet) ProtXML files
    
    This class is used to load (storing not supported, yet) documents that implement 
    the schema of ProtXML files.
		
		A documented schema for this format comes with the TPP.
		
		OpenMS can only read parts of the protein_summary subtree to extract protein-peptide associations. All other parts are silently ignored.

  	@note This format will eventually be replaced by the HUPO-PSI (mzIdentML and mzQuantML) AnalysisXML formats!
  	
  	@ingroup FileIO
  */
  class OPENMS_DLLAPI ProtXMLFile
  	: protected Internal::XMLHandler,
  		public Internal::XMLFile
  {
		public:
		
			/// Constructor
			ProtXMLFile();
		
			/**
				@brief Loads the identifications of an ProtXML file without identifier
				
				The information is read in and the information is stored in the
				corresponding variables

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
			void load(const String& filename, std::vector<ProteinIdentification>& protein_ids, std::vector<PeptideIdentification>& peptide_ids);

			/**
				@brief [not implemented yet!] Stores the data in an ProtXML file
				
				[not implemented yet!]
				The data is stored in the file 'filename'.

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
			*/
			void store(const String& filename, const std::vector<ProteinIdentification>& protein_ids, const std::vector<PeptideIdentification>& peptide_ids ,const String& document_id=""); 
  	
  	protected:

			/// reset members after reading/writing
			void resetMembers_();

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

			/// probability of a protein group
			DoubleReal p_protein_group_;
			/// probability of a protein
			DoubleReal protein_p_;
			/// protein identifier
			String protein_name_;
			/// peptide sequence
			String peptide_seq_;
			/// peptide nsp adjusted probability
			DoubleReal peptide_nspp_;

			//@}
  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_PROTXMLFILE_H
