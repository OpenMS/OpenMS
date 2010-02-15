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
// $Maintainer: Andreas Bertsch$
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MZIDENTMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZIDENTMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/METADATA/IdentificationHit.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

namespace OpenMS
{
	class ProgressLogger;

	namespace Internal
	{

		/**
			@brief XML handler for MzIdentMLFile
			
			@note Do not use this class. It is only needed in MzIdentMLFile.
		*/
		class OPENMS_DLLAPI MzIdentMLHandler
			: public XMLHandler
		{
		 public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler
      MzIdentMLHandler(const Identification& id, const String& filename, const String& version, const ProgressLogger& logger);

      /// Constructor for a read-only handler
      MzIdentMLHandler(Identification& id, const String& filename, const String& version, const ProgressLogger& logger);

      /// Destructor
      virtual ~MzIdentMLHandler();
      //@}


			// Docu in base class
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
      virtual void characters(const XMLCh* const chars, const XMLSize_t length);
			
			//Docu in base class
			virtual void writeTo(std::ostream& os);
			
		 protected:
      
			/// Progress logger
			const ProgressLogger& logger_;
			
			///Controlled vocabulary (psi-pi from OpenMS/share/OpenMS/CV/psi-pi.obo)
			ControlledVocabulary cv_;
		
			String tag_;

			Identification* id_;

			const Identification* cid_;

			SpectrumIdentification current_spectrum_id_;

			IdentificationHit current_id_hit_;

			/// Handles CV terms
			void handleCVParam_(const String& parent_parent_tag, const String& parent_tag, const String& accession, const String& name, const String& value, const xercesc::Attributes& attributes, const String& cv_ref,  const String& unit_accession="");

			/// Handles user terms
			void handleUserParam_(const String& parent_parent_tag, const String& parent_tag, const String& name, const String& type, const String& value);
			
			/// Writes user terms
			void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, UInt indent) const;
			
			/// Looks up a child CV term of @p parent_accession with the name @p name. If no such term is found, an empty term is returned.
			ControlledVocabulary::CVTerm getChildWithName_(const String& parent_accession, const String& name) const;
			
			/// Helper method that writes a source file
			//void writeSourceFile_(std::ostream& os, const String& id, const SourceFile& software);


			private:

				MzIdentMLHandler();
				MzIdentMLHandler(const MzIdentMLHandler& rhs);
				MzIdentMLHandler& operator = (const MzIdentMLHandler& rhs);
				Map<String, AASequence> pep_sequences_;
				AASequence actual_peptide_;
				Int current_mod_location_;
				ProteinHit actual_protein_;
				
		};
	} // namespace Internal
} // namespace OpenMS

#endif
