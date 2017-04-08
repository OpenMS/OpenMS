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
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/CHEMISTRY/EnzymesDB.h>

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

      // Maps enzyme_num in xQuest result file to the enzyme name used by OpenMS
      static std::map< Size, String > enzymes;

      // Maps String encoding month to the numeric value
      static std::map<String, UInt> months;

      // Decoy string used by xQuest
      static const String decoy_string;

      /**
       * @brief Removes a substring from a larger string
       * @param large String from which `small` should be removed
       * @param small Substring to be removed from `large`
       */
      static void removeSubstring(String & large, const String & small);


      XQuestResultXMLHandler(const String & /* filename */,
                             std::vector< XQuestResultMeta > & /* metas */,
                             std::vector< std::vector< PeptideIdentification > > & /* csms */,
                             std::vector< ProteinIdentification > &,
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
      //void characters(const XMLCh * const chars, const XMLSize_t /*length*/);

      //Docu in base class
      //virtual void writeTo(std::ostream & os);


    private:
   
      // Main data structures that are populated during loading the file
      std::vector< XQuestResultMeta > & metas_;
      std::vector< std::vector< PeptideIdentification > > & csms_;
      std::vector< ProteinIdentification > & prot_ids_; 
      
      // Whether or not current xquest result tag comes from OpenProXL (xQuest otherwise)
      bool is_openproxl_;

      // Set of all protein accessions that are within the ProteinHits.
      std::set< String > accessions_;
      
      // The enzyme database for enzyme lookup
      EnzymesDB * enzymes_db_;

      // Keeps track of the charges of the hits
      std::set< UInt > charges_;
      UInt min_precursor_charge_;
      UInt max_precursor_charge_;

      // Current Retention time of light spectrum
      double rt_light_;

      // The masses of the Monolinks
      std::set< double > monolinks_masses_;
      
      // The current spectrum search
      std::vector< PeptideIdentification > current_spectrum_search_;
      XQuestResultMeta current_meta_; // The current meta value
      int & n_hits_; // Total no. of hits found in the result XML file
      std::vector< int > * cum_hits_;
      size_t min_n_ions_per_spectrum_;
      bool load_to_peptideHit_;  // Whether Meta data of peptide identification should also be loaded to peptide hit
      
      // Stores the attributes of a record (peptide identification)
      std::map<String, DataValue> peptide_id_meta_values_;

      /**
       * @brief Extracts the DateTime from datetime string from xQuest
       * @param xquest_datetime_string The DateTime String to be processed
       * @param date_time DateTime that reflects the value given in the `xquest_datetime_string`
       */
      inline void extractDateTime_(const String & xquest_datetime_string, DateTime & date_time);

      /**
       * @brief Assignes all meta values stored in the peptide_id_attributes member to an meta info interface
       * @param meta_info_interface Where the meta values from the peptide_id_attributes member should be assigned to
       */
      void addMetaValues_(MetaInfoInterface & meta_info_interface);

      /**
       * @brief Gets the link location of a xQuest xlinkPositionString
       */
      void getLinkPosition_(const xercesc::Attributes &, std::pair<SignedSize, SignedSize> &);
      
      /**
       * @brief Sets the peptide Evidence for Alpha and Beta
       */
      void setPeptideEvidence_(const String &, PeptideHit &);

      /*
       * Sets the meta data for one or both peptide hits
       */
      void setMetaValue_(const String & /* key */, const DataValue & /* datavalue */, PeptideIdentification & /* pep_id */, PeptideHit & /* alpha */);
      void setMetaValue_(const String & /* key */, const DataValue & /* datavalue */, PeptideIdentification & /* pep_id */, PeptideHit & /* alpha */, PeptideHit & /* beta */);
    };
  } // namespace Internal
} // namespace OpenMS
#endif // OPENMS_FORMAT_HANDLERS_XQUESTRESULTXMLHANDLER_H
