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
// $Maintainer: Hannes Roest $
// $Authors: Andreas Bertsch, Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_TRAMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_TRAMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{
  namespace Internal
  {

    /**
      @brief XML handler for TraMLFile
      
      @note Do not use this class. It is only needed in TraMLFile.
    */
    class OPENMS_DLLAPI TraMLHandler
      : public XMLHandler
    {
     public:

      typedef std::vector<ReactionMonitoringTransition::Product> ProductListType;
      typedef std::vector<ReactionMonitoringTransition::Configuration> ConfigurationListType;

      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler
      TraMLHandler(const TargetedExperiment& exp, const String& filename, const String& version, const ProgressLogger& logger);

      /// Constructor for a read-only handler
      TraMLHandler(TargetedExperiment& exp, const String& filename, const String& version, const ProgressLogger& logger);

      /// Destructor
      virtual ~TraMLHandler();
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
  
      CVTermList actual_interpretation_;

      std::vector<ReactionMonitoringTransition::Product> actual_intermediate_products_;

      ReactionMonitoringTransition::Product actual_product_;
      
      ReactionMonitoringTransition::Configuration actual_configuration_;

      SourceFile actual_sourcefile_;

      /// Handles CV terms
      void handleCVParam_(const String& parent_parent_tag, const String& parent_tag, const CVTerm& cv_term);

      /// Handles user terms
      void handleUserParam_(const String& parent_parent_tag, const String& parent_tag, const String& name, const String& type, const String& value);
      
      /// Writes user terms
      void writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, UInt indent) const;

      void writeUserParams_(std::ostream& os, const std::vector<MetaInfoInterface>& meta, UInt indent) const;

      void writeCVParams_(std::ostream& os, const CVTermList& cv_terms, UInt indent) const;

      // subfunctions of write
      void write_target_ (std::ostream& os, const std::vector<IncludeExcludeTarget>::const_iterator& it) const;

      void write_product_ (std::ostream& os, const std::vector<ReactionMonitoringTransition::Product>::const_iterator& prod_it) const;

      void write_configuration_ (std::ostream& os, const std::vector<ReactionMonitoringTransition::Configuration>::const_iterator& cit) const;
      
      /// Looks up a child CV term of @p parent_accession with the name @p name. If no such term is found, an empty term is returned.
      ControlledVocabulary::CVTerm getChildWithName_(const String& parent_accession, const String& name) const;
      
      /// Helper method that writes a source file
      //void writeSourceFile_(std::ostream& os, const String& id, const SourceFile& software);

      private:

        TraMLHandler();
        TraMLHandler(const TraMLHandler& rhs);
        TraMLHandler& operator = (const TraMLHandler& rhs);
    };
  } // namespace Internal
} // namespace OpenMS

#endif
