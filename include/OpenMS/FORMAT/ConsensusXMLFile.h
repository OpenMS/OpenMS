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
// $Maintainer: Clemens Groepl $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_CONSENSUSXMLFILE_H
#define OPENMS_FORMAT_CONSENSUSXMLFILE_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{
  /**
	@brief This class provides Input functionality for ConsensusMaps and Output functionality for
	alignments and quantitation.

	This class can be used to load the content of a consensusXML file into a ConsensusMap
	or to save the content of an ConsensusMap object into an XML file.

	A documented schema for this format can be found at http://open-ms.sourceforge.net/schemas/.

  @todo Take care that unique ids are assigned properly by TOPP tools before calling ConsensusXMLFile::store().  There will be a message on LOG_INFO but we will make no attempt to fix the problem in this class.  (all developers)

	@ingroup FileIO
  */
  class OPENMS_DLLAPI ConsensusXMLFile
  	: public Internal::XMLHandler,
  		public Internal::XMLFile,
			public ProgressLogger
  {
		public:
			///Default constructor
			ConsensusXMLFile();
			///Destructor
			~ConsensusXMLFile();


			/**
			@brief Loads a consensus map from file

			@exception Exception::FileNotFound is thrown if the file could not be opened
			@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
			void load(const String& filename, ConsensusMap& map);

			/**
			@brief Stores a consensus map to file

			@exception Exception::UnableToCreateFile is thrown if the file name is not writable
			@exception Exception::IllegalArgument is thrown if the consensus map is not valid
			*/
			void store(const String& filename, const ConsensusMap& consensus_map);

			/// Mutable access to the options for loading/storing
			PeakFileOptions& getOptions();

			/// Non-mutable access to the options for loading/storing
			const PeakFileOptions& getOptions() const;

		protected:

			// Docu in base class
			virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);

			// Docu in base class
			virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);

			// Docu in base class
			virtual void characters(const XMLCh* const chars, const XMLSize_t length);


			/// Writes a peptide identification to a stream (for assigned/unassigned peptide identifications)
			void writePeptideIdentification_(const String& filename, std::ostream& os, const PeptideIdentification& id, const String& tag_name, UInt indentation_level);


			/// Options that can be set
			PeakFileOptions options_;

			///@name Temporary variables for parsing
			//@{
			ConsensusMap* consensus_map_;
			ConsensusFeature act_cons_element_;
			DPosition<2> pos_;
			DoubleReal it_;
			UInt last_map_;
			//@}

			/// Pointer to last read object as a MetaInfoInterface, or null.
			MetaInfoInterface* last_meta_;
			/// Temporary protein ProteinIdentification
			ProteinIdentification prot_id_;
			/// Temporary peptide ProteinIdentification
			PeptideIdentification pep_id_;
			/// Temporary protein hit
			ProteinHit prot_hit_;
			/// Temporary peptide hit
			PeptideHit pep_hit_;
			/// Map from protein id to accession
			Map<String,String> proteinid_to_accession_;
			/// Map from search identifier concatenated with protein accession to id
			Map<String,Size> accession_to_id_;
			/// Map from identification run identifier to file xs:id (for linking peptide identifications to the corresponding run)
			Map<String,String> identifier_id_;
			/// Map from file xs:id to identification run identifier (for linking peptide identifications to the corresponding run)
			Map<String,String> id_identifier_;
			/// Temporary search parameters file
			ProteinIdentification::SearchParameters search_param_;

			UInt progress_;

  };
} // namespace OpenMS

#endif // OPENMS_FOMAT_CONSENSUSXMLFILE_H
