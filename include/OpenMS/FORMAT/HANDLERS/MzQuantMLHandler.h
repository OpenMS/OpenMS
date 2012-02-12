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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MZQUANTMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZQUANTMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{
	class ProgressLogger;

	namespace Internal
	{

		/**
			@brief XML handler for MzQuantMLFile

			@note Do not use this class. It is only needed in MzQuantMLFile.
		*/
		class OPENMS_DLLAPI MzQuantMLHandler
			: public XMLHandler
		{
		 public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler
      MzQuantMLHandler(const ConsensusMap& consensus_map, /* const FeatureMap& feature_map, */ const String& filename, const String& version, const ProgressLogger& logger);

      /// Constructor for a read-only handler
      MzQuantMLHandler(ConsensusMap& consensus_map, /* FeatureMap& feature_map, */ const String& filename, const String& version, const ProgressLogger& logger);

      /// Destructor
      virtual ~MzQuantMLHandler();
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

			/// Controlled vocabulary (hopefully the psi-pi from OpenMS/share/OpenMS/CV/psi-pi.obo)
			ControlledVocabulary cv_;

			String tag_;

			//~ FeatureMap* fm_;
			ConsensusMap* cm_;

			//~ const FeatureMap* cfm_;
			const ConsensusMap* ccm_;

			/// Handles CV terms
			void handleCVParam_(const String& parent_parent_tag, const String& parent_tag, const String& accession, /* const String& name, */ /* const String& value, */ const xercesc::Attributes& attributes, const String& cv_ref/* ,  const String& unit_accession="" */);

			/// Handles user terms
			void handleUserParam_(const String& parent_parent_tag, const String& parent_tag, const String& name, const String& type, const String& value);

			/// Writes user terms
			void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, UInt indent);

			void writeUserParam_(String& s, const MetaInfoInterface& meta, UInt indent);

			/// Looks up a child CV term of @p parent_accession with the name @p name. If no such term is found, an empty term is returned.
			ControlledVocabulary::CVTerm getChildWithName_(const String& parent_accession, const String& name) const;


			/// Helper method that writes a feature
			void writeFeature_(std::ostream& os, const String& identifier_prefix, UInt64 identifier, UInt indentation_level);

			/// Helper method that writes a source file
			//void writeSourceFile_(std::ostream& os, const String& id, const SourceFile& software);


			private:
				MzQuantMLHandler();
				MzQuantMLHandler(const MzQuantMLHandler& rhs);
				MzQuantMLHandler& operator = (const MzQuantMLHandler& rhs);
				enum QUANT_TYPES {MS1LABEL=0, MS2LABEL, LABELFREE, SIZE_OF_QUANT_TYPES};

				//~ Double ratio_;

		};
	} // namespace Internal
} // namespace OpenMS

#endif
