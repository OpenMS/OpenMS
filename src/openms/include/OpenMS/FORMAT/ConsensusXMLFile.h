// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{
  class ConsensusMap;
  /**
    @brief This class provides Input functionality for ConsensusMaps and Output functionality for
    alignments and quantitation.

    This class can be used to load the content of a consensusXML file into a ConsensusMap
    or to save the content of an ConsensusMap object into an XML file.

    A documented schema for this format can be found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS

  @todo Take care that unique ids are assigned properly by TOPP tools before calling ConsensusXMLFile::store().  There will be a message on OPENMS_LOG_INFO but we will make no attempt to fix the problem in this class.  (all developers)

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ConsensusXMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    ///Default constructor
    ConsensusXMLFile();
    ///Destructor
    ~ConsensusXMLFile() override;

    /**
    @brief Loads a consensus map from file and calls updateRanges

    @exception Exception::FileNotFound is thrown if the file could not be opened
    @exception Exception::ParseError is thrown if an error occurs during parsing
    @exception Exception::MissingInformation is thrown if source files are missing/duplicated or map-IDs are referencing non-existing maps
    */
    void load(const String& filename, ConsensusMap& map);

    /**
    @brief Stores a consensus map to file

    @exception Exception::UnableToCreateFile is thrown if the file name is not writable
    @exception Exception::IllegalArgument is thrown if the consensus map is not valid
    @exception Exception::MissingInformation is thrown if source files are missing/duplicated or map-IDs are referencing non-existing maps
    */
    void store(const String& filename, const ConsensusMap& consensus_map);

    /// Mutable access to the options for loading/storing
    PeakFileOptions& getOptions();

    /// Non-mutable access to the options for loading/storing
    const PeakFileOptions& getOptions() const;

protected:

    /// Options that can be set
    PeakFileOptions options_;
  };
} // namespace OpenMS

