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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_OMSSAXMLFILE_H
#define OPENMS_FORMAT_OMSSAXMLFILE_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>

#include <vector>

namespace OpenMS 
{
	class String;
	class ModificationsDB;
	
  /**
    @brief Used to load OMSSAXML files
    
    This class is used to load documents that implement 
    the schema of OMSSAXML files.

  	@ingroup FileIO
  */
  class OPENMS_DLLAPI OMSSAXMLFile 
		: protected Internal::XMLHandler, 
			public Internal::XMLFile
  {
    public:
						
      /// Default constructor
      OMSSAXMLFile();
			
			/// Destructor
			virtual ~OMSSAXMLFile();
		  /**
		    @brief loads data from a OMSSAXML file
		    
		    @param filename the file to be loaded
		    @param protein_identification protein identifications belonging to the whole experiment
		    @param id_data the identifications with m/z and RT
				@param load_proteins if this flag is set to false, the protein identifications are not loaded

		    This class serves to read in a OMSSAXML file. The information can be 
		    retrieved via the load function.      
		  	
				@exception FileNotFound
				@exception ParseError
				
		  	@ingroup FileIO
		  */
	    void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data, bool load_proteins = true);

			/// sets the valid modifications
			void setModificationDefinitionsSet(const ModificationDefinitionsSet& rhs);

		protected:
			// Docu in base class
      void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
      void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
		
			// Docu in base class
   		void characters(const XMLCh* const chars, const XMLSize_t /*length*/);

		private:

			OMSSAXMLFile(const OMSSAXMLFile& rhs);

			OMSSAXMLFile& operator = (const OMSSAXMLFile& rhs);
 	
			/// reads the mapping file needed for modifications
			void readMappingFile_();
			
			/// the identifications (storing the peptide hits)
      std::vector<PeptideIdentification>* peptide_identifications_;

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

			/// modifications mapping file from OMSSA mod num to UniMod accession
			Map<UInt, std::vector<ResidueModification> > mods_map_;

			/// modification mapping reverse, from the modification to the mod_num
			Map<String, UInt> mods_to_num_;

			/// modification definitions set of the search, needed to annotate fixed modifications
			ModificationDefinitionsSet mod_def_set_;
  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_OMSSAXMLFILE_H
