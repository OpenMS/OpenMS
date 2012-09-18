// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
