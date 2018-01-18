// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MASCOTXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MASCOTXMLHANDLER_H

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>

#include <vector>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Handler that is used for parsing MascotXML data
    */
    class OPENMS_DLLAPI MascotXMLHandler :
      public XMLHandler
    {
public:
      /// Constructor
      MascotXMLHandler(ProteinIdentification& protein_identification,
                       std::vector<PeptideIdentification>& identifications,
                       const String& filename,
                       std::map<String, std::vector<AASequence> >& peptides,
                       const SpectrumMetaDataLookup& lookup);

      /// Destructor
      ~MascotXMLHandler() override;

      // Docu in base class
      void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

      // Docu in base class
      void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

      // Docu in base class
      void characters(const XMLCh* const chars, const XMLSize_t /*length*/) override;
      
      /// Split modification search parameter if for more than one amino acid specified e.g. Phospho (ST)
      static std::vector<String> splitModificationBySpecifiedAA(String mod);

private:

      ProteinIdentification& protein_identification_; ///< the protein identifications
      std::vector<PeptideIdentification>& id_data_; ///< the identifications (storing the peptide hits)
      ProteinHit actual_protein_hit_;
      PeptideHit actual_peptide_hit_;
      PeptideEvidence actual_peptide_evidence_;
      UInt peptide_identification_index_;
      String tag_;
      DateTime date_;
      String date_time_string_;
      UInt actual_query_;
      ProteinIdentification::SearchParameters search_parameters_;
      String identifier_;
      String actual_title_;
      std::map<String, std::vector<AASequence> >& modified_peptides_;

      StringList tags_open_; ///< tracking the current XML tree
      String character_buffer_; ///< filled by MascotXMLHandler::characters
      String major_version_;
      String minor_version_;
      
      // list of modifications, which cannot be set as fixed and needs
      // to be removed, because added from mascot as variable modification
      std::vector<String> remove_fixed_mods_;

      /// Helper object for looking up RT information
      const SpectrumMetaDataLookup& lookup_;

      /// Error for missing RT information already reported?
      bool no_rt_error_;
    };

  } // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_MASCOTXMLHANDLER_H
