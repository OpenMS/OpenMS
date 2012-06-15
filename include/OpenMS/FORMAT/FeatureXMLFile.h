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
// $Maintainer: Chris Bielow, Clemens Groepl $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_FEATUREXMLFILE_H
#define OPENMS_FORMAT_FEATUREXMLFILE_H

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <iostream>

namespace OpenMS
{
	/**
  	@brief This class provides Input/Output functionality for feature maps

		A documented schema for this format can be found at http://open-ms.sourceforge.net/schemas/.

    @todo Take care that unique ids are assigned properly by TOPP tools before calling FeatureXMLFile::store().  There will be a message on LOG_INFO but we will make no attempt to fix the problem in this class.  (all developers)

  	@note This format will eventually be replaced by the HUPO-PSI AnalysisXML (mzIdentML and mzQuantML) formats!

  	@ingroup FileIO
  */
  class OPENMS_DLLAPI FeatureXMLFile
  	: protected Internal::XMLHandler,
  		public Internal::XMLFile, 
  		public ProgressLogger
  {

		public:

			/** @name Constructors and Destructor */
			//@{
			///Default constructor
			FeatureXMLFile();
			///Destructor
			~FeatureXMLFile();
			//@}

			/**
				@brief loads the file with name @p filename into @p map.

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
			void load(const String& filename, FeatureMap<>& feature_map);

      Size loadSize(const String& filename);

			/**
				@brief stores the map @p feature_map in file with name @p filename.

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
			*/
			void store(const String& filename, const FeatureMap<>& feature_map);

      /// Mutable access to the options for loading/storing
      FeatureFileOptions& getOptions();

      /// Non-mutable access to the options for loading/storing
      const FeatureFileOptions& getOptions() const;

		protected:

      // restore default state for next load/store operation
      void resetMembers_();

			// Docu in base class
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);

			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);

			// Docu in base class
      virtual void characters(const XMLCh* const chars, const XMLSize_t length);

			/// Writes a feature to a stream
			void writeFeature_(const String& filename, std::ostream& os, const Feature& feat, const String& identifier_prefix, UInt64 identifier, UInt indentation_level);

			/// Writes a peptide identification to a stream (for assigned/unassigned peptide identifications)
			void writePeptideIdentification_(const String& filename, std::ostream& os, const PeptideIdentification& id, const String& tag_name, UInt indentation_level);


			/**
				@brief update the pointer to the current feature

				@param create If true, a new (empty) Feature is added at the appropriate subordinate_feature_level_
			*/
			void updateCurrentFeature_(bool create);

      /// allows for early return in parsing functions when certain sections should be ignored
      /// <=0 - parsing ON
      ///  >0 - this number of tags have been entered that forbid parsing and need to be exited before parsing continues
      Int disable_parsing_;

			/// points to the last open &lt;feature&gt; tag (possibly a subordinate feature)
			Feature* current_feature_;
			/// Feature map pointer for reading
			FeatureMap<Feature>* map_;
			/// Options that can be set
			FeatureFileOptions options_;
      /// only parse until "count" tag is reached (used in loadSize())
      bool size_only_;
      /// holds the putative size given in count
      Size expected_size_;

			/**@name temporary data structures to hold parsed data */
	    //@{
			ModelDescription<2> model_desc_;
			Param param_;
			ConvexHull2D::PointArrayType current_chull_;
			DPosition<2> hull_position_;
	    //@}

			/// current dimension of the feature position, quality, or convex hull point
	 		UInt dim_;

			/// for downward compatibility, all tags in the old description must be ignored
			bool in_description_;

			/// level in Feature stack during parsing
			Int subordinate_feature_level_;

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

	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_FeatureXMLFile_H
