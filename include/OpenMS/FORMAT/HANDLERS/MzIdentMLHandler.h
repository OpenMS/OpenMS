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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer, Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MZIDENTMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZIDENTMLHANDLER_H

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/METADATA/IdentificationHit.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

namespace OpenMS
{
  class ProgressLogger;

  namespace Internal
  {

    /**
        @brief XML handler for MzIdentMLFile

        @note Do not use this class. It is only needed in MzIdentMLFile.
    */
    class OPENMS_DLLAPI MzIdentMLHandler :
      public XMLHandler
    {
public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler
      MzIdentMLHandler(const Identification & id, const String & filename, const String & version, const ProgressLogger & logger);
      /// Constructor for a write-only handler for internal identification structures
      MzIdentMLHandler(const std::vector<ProteinIdentification> & pro_id, const std::vector<PeptideIdentification> & pep_id, const String & filename, const String & version, const ProgressLogger & logger);

      /// Constructor for a read-only handler
      MzIdentMLHandler(Identification & id, const String & filename, const String & version, const ProgressLogger & logger);
      /// Constructor for a read-only handler for internal identification structures
      MzIdentMLHandler(std::vector<ProteinIdentification> & pro_id, std::vector<PeptideIdentification> & pep_id, const String & filename, const String & version, const ProgressLogger & logger);

      /// Destructor
      virtual ~MzIdentMLHandler();
      //@}


      // Docu in base class
      virtual void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname);

      // Docu in base class
      virtual void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes);

      // Docu in base class
      virtual void characters(const XMLCh * const chars, const XMLSize_t length);

      //Docu in base class
      virtual void writeTo(std::ostream & os);

protected:

      /// Progress logger
      const ProgressLogger & logger_;

      ///Controlled vocabulary (psi-ms from OpenMS/share/OpenMS/CV/psi-ms.obo)
      ControlledVocabulary cv_;
      ///Controlled vocabulary for modifications (unimod from OpenMS/share/OpenMS/CV/unimod.obo)
      ControlledVocabulary unimod_;

      //~ MSExperiment<>* ms_exp_;

      ///XML tag parse element
      String tag_;

      ///Identification Item
      Identification * id_;
      ///internal Identification Item for proteins
      std::vector<ProteinIdentification> * pro_id_;
      ///Identification Item for peptides
      std::vector<PeptideIdentification> * pep_id_;

      const Identification * cid_;
      const std::vector<ProteinIdentification> * cpro_id_;
      const std::vector<PeptideIdentification> * cpep_id_;

      ///SpectrumIdentification Item
      SpectrumIdentification current_spectrum_id_;

      ///IdentificationHit Item
      IdentificationHit current_id_hit_;

      /// Handles CV terms
      void handleCVParam_(const String & parent_parent_tag, const String & parent_tag, const String & accession, /* const String& name, */ /* const String& value, */ const xercesc::Attributes & attributes, const String & cv_ref /* ,  const String& unit_accession="" */);

      /// Handles user terms
      void handleUserParam_(const String & parent_parent_tag, const String & parent_tag, const String & name, const String & type, const String & value);

      /// Writes user terms
      void writeUserParam_(std::ostream & os, const MetaInfoInterface & meta, UInt indent) const;

      /// Looks up a child CV term of @p parent_accession with the name @p name. If no such term is found, an empty term is returned.
      ControlledVocabulary::CVTerm getChildWithName_(const String & parent_accession, const String & name) const;

      /// Helper method that writes a source file
      //void writeSourceFile_(std::ostream& os, const String& id, const SourceFile& software);


private:
      MzIdentMLHandler();
      MzIdentMLHandler(const MzIdentMLHandler & rhs);
      MzIdentMLHandler & operator=(const MzIdentMLHandler & rhs);
      Map<String, AASequence> pep_sequences_;
      AASequence actual_peptide_;
      Int current_mod_location_;
      ProteinHit actual_protein_;

    };
  }   // namespace Internal
} // namespace OpenMS

#endif
