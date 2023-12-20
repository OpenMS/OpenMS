// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Svetlana Kutuzova, Douglas McCloskey $
// $Authors: Svetlana Kutuzova, Douglas McCloskey $
// --------------------------------------------------------------------------
 
#pragma once

#include <OpenMS/FORMAT/CsvFile.h>
#include <map>

namespace OpenMS
{
class MzMLFile;
/*
    @brief Scheduler for FIA-MS data batches. Works with FIAMSDataProcessor.

    The class is initialised with the path to the csv file that must contain the following columns:
    - *filename* - the mzML filename for the sample. Must not contains the extension and path. The results filenames follow the pattern @filename_@n_seconds.
    - *dir_input* - the location of the directory with the mzML file
    - *dir_output* - the location of the directory where the results will be stored
    - *resolution* - the resolution of the instrument
    - *polarity* - the charge of the instrument, accepts "positive" or "negative"
    - *db:mapping* - database input file containing three tab-separated columns of mass, formula, identifier for the accurate mass search
    - *db:struct* - database input file containing four tab-separated columns of identifier, name, SMILES, INCHI for the accurate mass search
    - *positive_adducts* - file containing the list of potential positive adducts for the accurate mass search
    - *negative_adducts* - file containing the list of potential negative adducts for the accurate mass search
    - *time* - ";"-separated numbers of seconds to process, f.e. "30;60;90;180"
*/
class OPENMS_DLLAPI FIAMSScheduler 
  {
public:
    FIAMSScheduler() = default;

    /// Default constructor
    FIAMSScheduler(
      String filename,
      String base_dir = "/",
      bool load_cached_ = true
    );

    /// Default destructor
    ~FIAMSScheduler() = default;

    /// Copy constructor
    FIAMSScheduler(const FIAMSScheduler& cp) = default;

    /// Assignment
    FIAMSScheduler& operator=(const FIAMSScheduler& fdp) = default;

    /**
      @brief Run the FIA-MS data analysis for the batch defined in the @p filename_
    */
    void run();

    /**
      @brief Get the batch
    */
    const std::vector<std::map<String, String>>& getSamples();

    /**
      @brief Get the base directory for the relevant paths from the csv file
    */
    const String& getBaseDir();

private:
    /**
      @brief Load the batch from the csv file and store as the vector of maps
    */
    void loadSamples_();

    String filename_;
    String base_dir_;
    bool load_cached_;
    std::vector<std::map<String, String>> samples_;
  };
} // namespace OpenMS
