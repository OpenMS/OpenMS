// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <vector>

namespace OpenMS
{

  /**
     @brief Class for storing Percolator tab-delimited input files.

  */
  class OPENMS_DLLAPI PercolatorInfile
  {
    public:
      static void store(const String& pin_file, 
        const std::vector<PeptideIdentification>& peptide_ids, 
        const StringList& feature_set, 
        const std::string& enz, 
        int min_charge, 
        int max_charge);

      /** @brief load pin file and convert to a vector of PeptideIdentification using the given score column @p score_name and orientation @p higher_score_better.
          If a decoy prefix is provided, the decoy status is set from the protein accessions.
          Otherwise, it assumes that the pin file already contains the correctly annotated decoy status.
          If @p extra_scores is not empty, the scores are added to the PeptideHit as MetaValues.
          If a filename column is encountered the set of @p filenames is filled in the order of appearance and PeptideIdentifications annotated with the id_merge_index meta value to link them to the filename (similar to a merged idXML file). 
          TODO: implement something similar to PepXMLFile().setPreferredFixedModifications(getModifications_(fixed_modifications_names));
          **/
      static std::vector<PeptideIdentification> load(const String& pin_file, 
        bool higher_score_better, 
        const String& score_name, 
        const StringList& extra_scores,
        StringList& filenames, 
        String decoy_prefix = "");

      // uses spectrum_reference, if empty uses spectrum_id, if also empty fall back to using index
      static String getScanIdentifier(const PeptideIdentification& pid, size_t index);
      
    protected:

      //id <tab> label <tab> scannr <tab> calcmass <tab> expmass <tab> feature1 <tab> ... <tab> featureN <tab> peptide <tab> proteinId1 <tab> .. <tab> proteinIdM
      static TextFile preparePin_(
        const std::vector<PeptideIdentification>& peptide_ids, 
        const StringList& feature_set, 
        const std::string& enz, 
        int min_charge, 
        int max_charge);

      static bool isEnz_(const char& n, const char& c, const std::string& enz);

      static Size countEnzymatic_(const String& peptide, const std::string& enz);

  };
} // namespace OpenMS
