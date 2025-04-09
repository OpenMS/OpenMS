// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Andreas Bertsch, Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <iosfwd>
#include <map>

namespace OpenMS
{
  namespace Internal
  {

    /**
      @brief XML handler for TraMLFile

      @note Do not use this class. It is only needed in TraMLFile.
    */
    class OPENMS_DLLAPI TraMLHandler :
      public XMLHandler
    {
public:

      TraMLHandler() = delete;
      TraMLHandler(const TraMLHandler & rhs) = delete;
      TraMLHandler & operator=(const TraMLHandler & rhs) = delete;

      typedef std::vector<ReactionMonitoringTransition::Product> ProductListType;
      typedef std::vector<ReactionMonitoringTransition::Configuration> ConfigurationListType;

      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler
      TraMLHandler(const TargetedExperiment & exp, const String & filename, const String & version, const ProgressLogger & logger);

      /// Constructor for a read-only handler
      TraMLHandler(TargetedExperiment & exp, const String & filename, const String & version, const ProgressLogger & logger);

      /// Destructor
      ~TraMLHandler() override;
      //@}


      // Docu in base class
      void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname) override;

      // Docu in base class
      void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes) override;

      // Docu in base class
      void characters(const XMLCh * const chars, const XMLSize_t length) override;

      //Docu in base class
      void writeTo(std::ostream & os) override;

protected:

      /// Progress logger
      const ProgressLogger& logger_;

      ///Controlled vocabulary (psi-ms from OpenMS/share/OpenMS/CV/psi-ms.obo)
      ControlledVocabulary cv_;

      String tag_;

      TargetedExperiment* exp_;

      const TargetedExperiment* cexp_;

      TargetedExperiment::Publication actual_publication_;

      TargetedExperiment::Contact actual_contact_;

      TargetedExperiment::Instrument actual_instrument_;

      TargetedExperiment::Prediction actual_prediction_;

      Software actual_software_;

      TargetedExperiment::Protein actual_protein_;

      TargetedExperiment::RetentionTime actual_rt_;

      TargetedExperiment::Peptide actual_peptide_;

      TargetedExperiment::Compound actual_compound_;

      ReactionMonitoringTransition actual_transition_;

      IncludeExcludeTarget actual_target_;

      CVTermList actual_validation_;

      TargetedExperiment::Interpretation actual_interpretation_;

      std::vector<ReactionMonitoringTransition::Product> actual_intermediate_products_;

      ReactionMonitoringTransition::Product actual_product_;

      ReactionMonitoringTransition::Configuration actual_configuration_;

      SourceFile actual_sourcefile_;

      /// Handles CV terms
      void handleCVParam_(const String & parent_parent_tag, const String & parent_tag, const CVTerm & cv_term);

      /// Handles user terms
      void handleUserParam_(const String & parent_parent_tag, const String & parent_tag, const String & name, const String & type, const String & value);

      /// Writes user terms
      void writeUserParam_(std::ostream & os, const MetaInfoInterface & meta, UInt indent) const;

      void writeUserParams_(std::ostream & os, const std::vector<MetaInfoInterface> & meta, UInt indent) const;

      void writeCVParams_(std::ostream & os, const CVTermList & cv_terms, UInt indent) const;
      void writeCVParams_(std::ostream & os, const CVTermListInterface & cv_terms, UInt indent) const;

      void writeCVList_(std::ostream & os, const std::map<String, std::vector<CVTerm>> & cv_terms, UInt indent) const;

      // subfunctions of write
      void writeTarget_(std::ostream & os, const std::vector<IncludeExcludeTarget>::const_iterator & it) const;

      void writeRetentionTime_(std::ostream& os, const TargetedExperimentHelper::RetentionTime& rt) const;

      void writeProduct_(std::ostream & os, const std::vector<ReactionMonitoringTransition::Product>::const_iterator & prod_it) const;

      void writeConfiguration_(std::ostream & os, const std::vector<ReactionMonitoringTransition::Configuration>::const_iterator & cit) const;

      /// Looks up a child CV term of @p parent_accession with the name @p name. If no such term is found, an empty term is returned.
      ControlledVocabulary::CVTerm getChildWithName_(const String & parent_accession, const String & name) const;

      /// Helper method that writes a source file
      //void writeSourceFile_(std::ostream& os, const String& id, const SourceFile& software);

    };
  } // namespace Internal
} // namespace OpenMS


