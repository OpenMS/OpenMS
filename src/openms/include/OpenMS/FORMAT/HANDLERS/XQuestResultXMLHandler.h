// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann $
// --------------------------------------------------------------------------
#ifndef OPENMS_FORMAT_HANDLERS_XQUESTRESULTXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_XQUESTRESULTXMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/XQuestResultMeta.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>

namespace OpenMS
{
  namespace Internal
  {

    /** @brief XMLHandler for the result files of XQuest
     */
    class OPENMS_DLLAPI XQuestResultXMLHandler :
      public XMLHandler
    {
    public:

      XQuestResultXMLHandler(const String & /* filename */,
                             std::vector< XQuestResultMeta > & /* metas */,
                             std::vector< std::vector< PeptideIdentification > > & /* csms */,
                             int & n_hits,
                             std::vector< int > * cum_hits,
                             size_t min_n_ions_per_spectrum,
                             bool load_to_peptideHit_);
      virtual ~XQuestResultXMLHandler();

      // Docu in base class
      void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname);

      // Docu in base class
      void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes);

      // Docu in base class
      void characters(const XMLCh * const chars, const XMLSize_t /*length*/);

    private:

      std::vector< XQuestResultMeta > & metas_;
      std::vector< std::vector< PeptideIdentification > > & csms_;
      // The current spectrum search
      std::vector< PeptideIdentification > current_spectrum_search;
      // The current meta value
      XQuestResultMeta current_meta_;
      int & n_hits_; // Total no. of hits found in the result XML file
      std::vector< int > * cum_hits_;
      size_t min_n_ions_per_spectrum_;
      bool load_to_peptideHit_;

      // Stores the attributes of a record (peptide identification)
      std::map<String, DataValue> peptide_id_meta_values;

      // Assign all attributes in the peptide_id_attributes map to the MetaInfoInterface object
      void add_meta_values(MetaInfoInterface & meta_info_interface);

      // Retrieves the link location for cross-links and loop links
      void getLinkPosition_(const xercesc::Attributes &, std::pair<SignedSize, SignedSize> &);

      /*
       * Sets the meta data for one or both peptide hits
       */
      void setMetaValue(const String & /* key */, const DataValue & /* datavalue */, PeptideIdentification & /* pep_id */, PeptideHit & /* alpha */);
      void setMetaValue(const String & /* key */, const DataValue & /* datavalue */, PeptideIdentification & /* pep_id */, PeptideHit & /* alpha */, PeptideHit & /* beta */);
    };
  } // namespace Internal
} // namespace OpenMS
#endif // OPENMS_FORMAT_HANDLERS_XQUESTRESULTXMLHANDLER_H
