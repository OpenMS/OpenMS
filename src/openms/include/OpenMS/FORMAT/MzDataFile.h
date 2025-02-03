// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/FORMAT/HANDLERS/MzDataHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>

namespace OpenMS
{
  /**
      @brief File adapter for MzData files

      @ingroup FileIO
  */
  class OPENMS_DLLAPI MzDataFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
	typedef PeakMap MapType;

public:

    ///Default constructor
    MzDataFile();
    ///Destructor
    ~MzDataFile() override;

    /// Mutable access to the options for loading/storing
    PeakFileOptions & getOptions();

    /// Non-mutable access to the options for loading/storing
    const PeakFileOptions & getOptions() const;

    /// set options for loading/storing
    void setOptions(const PeakFileOptions &);

    /**
        @brief Loads a map from a MzData file.

        @p map has to be a MSExperiment or have the same interface.

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String & filename, MapType & map);

    /**
        @brief Stores a map in a MzData file.

        @p map has to be a MSExperiment or have the same interface.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String & filename, const MapType & map) const;

    /**
        @brief Checks if a file is valid with respect to the mapping file and the controlled vocabulary.

        @param filename File name of the file to be checked.
        @param errors Errors during the validation are returned in this output parameter.
        @param warnings Warnings during the validation are returned in this output parameter.

        @exception Exception::FileNotFound is thrown if the file could not be opened
    */
    bool isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings);

private:

    /// Options for loading / storing
    PeakFileOptions options_;
  };

} // namespace OpenMS


