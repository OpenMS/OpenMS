// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>

namespace OpenMS
{
  class OnDiscMSExperiment;
  typedef OpenMS::OnDiscMSExperiment OnDiscPeakMap;

  /**
    @brief A class to load an indexedmzML file.

    Providing the same interface as the other classes such as MzMLFile,
    MzXMLFile etc. to load and store a file. Reading a file from disk will load
    the file into a OnDiscMSExperiment while the class can write to disk both,
    a MSExperiment and a OnDiscMSExperiment.

  */
  class OPENMS_DLLAPI IndexedMzMLFileLoader
  {
    public:

    /// Constructor
    IndexedMzMLFileLoader();

    /// Destructor
    ~IndexedMzMLFileLoader();

    /// Mutable access to the options for loading/storing
    PeakFileOptions& getOptions();

    /// Non-mutable access to the options for loading/storing
    const PeakFileOptions& getOptions() const;

    /// set options for loading/storing
    void setOptions(const PeakFileOptions &);

    /**
      @brief Load a file 

      Tries to parse the file, success needs to be checked with the return value.

      @param filename Filename determines where the file is located
      @param exp Object which will contain the data after the call

      @return Indicates whether parsing was successful (if it is false, the file most likely was not an mzML or not indexed).
    */
    bool load(const String& filename, OnDiscPeakMap& exp);

    /**
      @brief Store a file from an on-disc data-structure

      @param filename Filename determines where the file will be stored 
      @param exp MS data to be stored
    */
    void store(const String& filename, OnDiscPeakMap& exp);

    /**
      @brief Store a file from an in-memory data-structure

      @param filename Filename determines where the file will be stored 
      @param exp MS data to be stored
    */
    void store(const String& filename, PeakMap& exp);

private:

    /// Options for storing
    PeakFileOptions options_;

  };
}


