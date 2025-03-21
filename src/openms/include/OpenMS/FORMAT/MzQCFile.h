// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>
#include <vector>

namespace OpenMS
{
  class FeatureMap;
  
  /**
      @brief File adapter for mzQC files used to load and store mzQC files

      This class collects the data for the mzQC File

      @ingroup FileIO
  */
  class OPENMS_DLLAPI MzQCFile
  {
  public:
    // Default constructor
    MzQCFile() = default;

    /**
      @brief Stores QC data in mzQC file with JSON format
      @param input_file mzML input file name
      @param output_file mzQC output file name
      @param exp MSExperiment to extract QC data from, prior sortSpectra() and updateRanges() required
      @param contact_name name of the person creating the mzQC file
      @param contact_address contact address (mail/e-mail or phone) of the person creating the mzQC file
      @param description description and comments about the mzQC file contents
      @param label unique and informative label for the run
      @param feature_map FeatureMap from feature file (featureXML)
      @param prot_ids protein identifications from ID file (idXML)
      @param pep_ids protein identifications from ID file (idXML)
    */
    void store(const String& input_file,
               const String& output_file,
               const MSExperiment& exp,
               const String& contact_name,
               const String& contact_address,
               const String& description,
               const String& label,
               const FeatureMap& feature_map,
               std::vector<ProteinIdentification>& prot_ids,
               std::vector<PeptideIdentification>& pep_ids) const;
  };
}