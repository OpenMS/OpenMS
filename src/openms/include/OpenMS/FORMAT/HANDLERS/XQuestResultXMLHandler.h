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
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann $
// --------------------------------------------------------------------------
#ifndef OPENMS_FORMAT_HANDLERS_XQUESTRESULTXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_XQUESTRESULTXMLHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

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

      XQuestResultXMLHandler(const String & filename,
                             std::vector< std::vector< PeptideIdentification > > & csms,
                             std::vector< ProteinIdentification > & prot_ids,
                             Size min_n_ions_per_spectrum,
                             bool load_to_peptideHit_);
      virtual ~XQuestResultXMLHandler();

      // Docu in base class
      void endElement(const XMLCh * const uri, const XMLCh * const local_name, const XMLCh * const qname);

      // Docu in base class
      void startElement(const XMLCh * const uri, const XMLCh * const local_name, const XMLCh * const qname, const xercesc::Attributes & attributes);

      /**
       * @brief Returns the minimum score encountered in the file.
       * @return Minimum score encountered in the file.
       */
      double getMinScore() const;

      /**
       * @brief Returns the maximum score encountered in the file.
       * @return Maximum score encountered in the file.
       */
      double getMaxScore() const;

      /**
       * @brief Returns the total number of hits in the file.
       * @return Total number of hits in the file.
       */
      UInt getNumberOfHits() const;

    private:

      // Main data structures that are populated during loading the file
      std::vector< std::vector< PeptideIdentification > > & csms_;
      std::vector< ProteinIdentification > & prot_ids_;

      UInt n_hits_; // Total no. of hits found in the result XML file

      // Keeps track of the minscore and maxscore encountered
      double min_score_;
      double max_score_;

      Size min_n_ions_per_spectrum_;
      bool load_to_peptideHit_;  // Whether Meta data of peptide identification should also be loaded to peptide hit

      // Whether or not current xquest result tag comes from OpenProXL (xQuest otherwise)
      bool is_openproxl_;

      // Set of all protein accessions that are within the ProteinHits.
      std::set< String > accessions_;

      // The enzyme database for enzyme lookup
      ProteaseDB* enzymes_db_;

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
       * @brief Gets the link location of a xQuest xlinkPositionString.
       * @param attributes XML attributes of Xerces.
       * @param pair Pair to be populated with the xlinkposition in xQuest.
       */
      void getLinkPosition_(const xercesc::Attributes & attributes, std::pair<SignedSize, SignedSize> & pair);

      /**
       * @brief Sets the peptide evidence for Alpha and Beta.
       * @param prot_string Protein string of the xquest file the peptide evidence should be populated from.
       * @param pep_hit For which peptide hit the peptide evidence should be set.
       */
      void setPeptideEvidence_(const String & prot_string, PeptideHit & pep_hit);

      /**
       * @brief Sets the meta value of the peptide identification for alpha hit.
       * @param key Which meta value to set
       * @param datavalue Value to be set
       * @param pep_id For which peptide identification the meta value should be set.
       * @param alpha Alpha peptide hit for which the meta value should be set.
       */
      void setMetaValue_(const String & key, const DataValue & datavalue, PeptideIdentification & pep_id, PeptideHit & alpha);

      /**
       * @brief Sets the meta value of the peptide identification for alpha hit.
       * @param key Which meta value to set
       * @param datavalue Value to be set
       * @param pep_id For which peptide identification the meta value should be set.
       * @param alpha Alpha peptide hit for which the meta value should be set.
       * @param beta Beta peptide hit for which the meta value should be set.
       */
      void setMetaValue_(const String & key, const DataValue & datavalue, PeptideIdentification & pep_id, PeptideHit & alpha, PeptideHit & beta);
    };
  } // namespace Internal
} // namespace OpenMS
#endif // OPENMS_FORMAT_HANDLERS_XQUESTRESULTXMLHANDLER_H
