// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_FORMAT_XTANDEMXMLFILE_H
#define OPENMS_FORMAT_XTANDEMXMLFILE_H

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>

#include <vector>

namespace OpenMS 
{
	class String;
  class ProteinIdentification;

  /**
    @brief Used to load XTandemXML files
    
    This class is used to load documents that implement 
    the schema of XTandemXML files.
  
  	@ingroup FileIO
  */
  class OPENMS_DLLAPI XTandemXMLFile
		: protected Internal::XMLHandler,
			public Internal::XMLFile
  {
    public:
						
      /// Default constructor
      XTandemXMLFile();
			
			/// Destructor
			virtual ~XTandemXMLFile();
		  /**
		    @brief loads data from a XTandemXML file
		    
		    @param filename the file to be loaded
		    @param protein_identification protein identifications belonging to the whole experiment
		    @param id_data the identifications with m/z and RT

		    This class serves to read in a XTandemXML file. The information can be 
		    retrieved via the load function.      
		  	
		  	@ingroup FileIO
		  */
	    void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data);


			/// sets the valid modifications
			void setModificationDefinitionsSet(const ModificationDefinitionsSet& rhs);
			
		protected:

			// Docu in base class
			void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
			void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);

      // Docu in base class
      void characters(const XMLCh* const chars, const XMLSize_t /*length*/);

			XTandemXMLFile(const XTandemXMLFile& rhs);

			XTandemXMLFile& operator = (const XTandemXMLFile& rhs);

		private:

			ProteinIdentification* protein_identification_;

			// used to indicate that an protein tag is open
			bool protein_open_;
     
			// true if actual 
			bool is_description_;
			
			// peptide hits of one spectrum
			Map<UInt, std::vector<PeptideHit> > peptide_hits_;

			// protein hits, sorted by the id
			Map<String, ProteinHit> protein_hits_;
			
			// id of the actual protein
			String actual_protein_id_;
			
			// charge of actual peptide
      Int actual_charge_;
			
			// id of actual peptide
      Int actual_id_;
			
			// tag
      String tag_;
			
			// actual start position of peptide in protein sequence
			UInt actual_start_;
			
			// actual stop position of peptide in protein sequence
			UInt actual_stop_;
	
			// modification definitions
			ModificationDefinitionsSet mod_def_set_;
      					 
  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_XTANDEMXMLFILE_H
