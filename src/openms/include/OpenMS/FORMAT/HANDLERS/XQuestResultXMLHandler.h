// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Lukas Zimmermann $
// --------------------------------------------------------------------------
#pragma once

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

      /// Maps enzyme_num in xQuest result file to the enzyme name used by OpenMS
      static std::map< Size, String > enzymes;

      /// Maps String encoding month to the numeric value
      static std::map<String, UInt> months;

      /// Constructor for a read-only handler for internal identification structures
      XQuestResultXMLHandler(const String & filename,
                             std::vector< PeptideIdentification > & pep_ids,
                             std::vector< ProteinIdentification > & prot_ids
                             );

      /// Constructor for a write-only handler for internal identification structures
      XQuestResultXMLHandler(const std::vector<ProteinIdentification>& pro_id,
                             const std::vector<PeptideIdentification>& pep_id,
                             const String& filename,
                             const String& version
                           );

      ~XQuestResultXMLHandler() override;

      // Docu in base class
      void endElement(const XMLCh * const uri, const XMLCh * const local_name, const XMLCh * const qname) override;

      // Docu in base class
      void startElement(const XMLCh * const uri, const XMLCh * const local_name, const XMLCh * const qname, const xercesc::Attributes & attributes) override;

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

      //Docu in base class
      void writeTo(std::ostream& os) override;

      // TODO move these to StringUtils?
      /**
        * @brief splits the @p input string at the nth occurrence of the @p separator

          If the separator does not occur in the input string n times, then the first output string will be the entire input string
          and the second one will be empty.

        * @return StringList with two elements, the two parts of the input without the nth separator
      */
      static StringList splitByNth(const String& input, const char separator, const Size n);

      /**
        * @brief counts occurrences of the @p separator and splits the @p string into two at the middle

          If the separator occurs 5 times in the input string, the string will be split at the 3rd occurrence.
          If 7 times, then at the 4th.
          The separator has to occur in the string an uneven number of times.
          If the separator occurs once, the string will be split at this one instance.
          If this one occurrence is at the beginning or end, one of the result strings will be empty.

        * @exception Exception::IllegalArgument is thrown if the @p separator does not occur in the @p input string an uneven number of times and at least once
        * @return StringList with two elements, the two halves of the input without the middle separator
      */
      static StringList splitByMiddle(const String& input, const char separator);

    private:


      // Decoy string used by xQuest, initialize to a default value
      String decoy_string_ = "decoy_";
      int spectrum_index_light_;
      int spectrum_index_heavy_;
      String cross_linker_name_;

      // Main data structures that are populated during loading the file
      std::vector< PeptideIdentification >* pep_ids_;
      std::vector< ProteinIdentification >* prot_ids_;

      // internal ID items for writing files
      const std::vector<ProteinIdentification>* cpro_id_;
      const std::vector<PeptideIdentification>* cpep_id_;

      UInt n_hits_; ///< Total no. of hits found in the result XML file

      // Keeps track of the minscore and maxscore encountered
      double min_score_;
      double max_score_;

      /// Whether or not current xquest result tag comes from OpenPepXL (xQuest otherwise)
      bool is_openpepxl_;

      /// Set of all protein accessions that are within the ProteinHits.
      std::set< String > accessions_;

      /// The enzyme database for enzyme lookup
      ProteaseDB* enzymes_db_;

      /// Keeps track of the charges of the hits
      std::set< UInt > charges_;
      UInt min_precursor_charge_;
      UInt max_precursor_charge_;

      // Current Retention time of spectrum pair
      double rt_light_;
      double rt_heavy_;

      // Current experimental m/z of spectrum pair
      double mz_light_;
      double mz_heavy_;

      // primary MS run path
      StringList ms_run_path_;
      String spectrum_input_file_;

      /// The current spectrum search
      std::vector< PeptideIdentification > current_spectrum_search_;

      /// Stores the attributes of a record (peptide identification)
      std::map<String, DataValue> peptide_id_meta_values_;

      /**
       * @brief Extracts the DateTime from datetime string from xQuest
       * @param xquest_datetime_string The DateTime String to be processed
       * @param date_time DateTime that reflects the value given in the `xquest_datetime_string`
       */
      inline void extractDateTime_(const String & xquest_datetime_string, DateTime & date_time) const;

      /**
       * @brief Assigns all meta values stored in the peptide_id_attributes
       * member to an meta info interface
       *
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

    };
  } // namespace Internal
} // namespace OpenMS
