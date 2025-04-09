// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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


      /**
      * @brief Loads peptide identifications from a Percolator input file.
      * 
      * This function reads a Percolator input file (`pin_file`) and returns a vector of `PeptideIdentification` objects. 
      * It extracts relevantinformation such as peptide sequences, scores, charges, annotations, and protein accessions, applying
      * specified thresholds and handling decoy targets as needed.
      * Note: If a filename column is encountered the set of @p filenames is filled in the order of appearance and PeptideIdentifications annotated with the id_merge_index meta value to link them to the filename (similar to a merged idXML file). 
      * 
      * @param pin_file he path to the Percolator input file with a `.pin` extension.
      * 
      * @param higher_score_better A boolean flag indicating whether higher scores are considered better (`true`) or lower scores are better (`false`).
      * 
      * @param score_name The name of the primary score to be used for ranking peptide hits.
      * 
      * @param extra_scores A list of additional score names that should be extracted and stored in each `PeptideHit`.
      * 
      * @param filenames Will be populated with the unique raw file names extracted from the input data.
      * 
      * @param decoy_prefix The prefix used to identify decoy protein accessions. Proteins with accessions starting with this prefix are marked as decoys. Otherwise, it assumes that the pin file already contains the correctly annotated decoy status.
      * @param threshold A double value representing the threshold for the `spectrum_q` value. Only spectra with `spectrum_q` below this threshold are processed.
                         Implemented to allow prefiltering of Sage results.
      * @param SageAnnotation A boolean value used to determine if the pin file is coming from Sage or not 
      * @return A `std::vector` of `PeptideIdentification` objects containing the peptide identifications.
      
      * @throws `Exception::ParseError` if any line in the input file does not have the expected number of columns.
      * TODO: implement something similar to PepXMLFile().setPreferredFixedModifications(getModifications_(fixed_modifications_names));      
      */
      static std::vector<PeptideIdentification> load(const String& pin_file, 
        bool higher_score_better, 
        const String& score_name, 
        const StringList& extra_scores,
        StringList& filenames, 
        String decoy_prefix = "",
        double threshold = 0.01, 
        bool SageAnnotation = false);

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
