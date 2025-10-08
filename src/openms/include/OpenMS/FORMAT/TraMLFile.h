// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

namespace OpenMS
{
  class Identification;
  class TargetedExperiment;

  /**
      @brief File adapter for HUPO PSI TraML files

      TraML files contain information about transitions used for targeted
      proteomics and metabolomics experiments:

      Deutsch et al. "TraML--a standard format for exchange of selected reaction monitoring transition lists."
      Mol Cell Proteomics. 2012 Apr;11(4):R111.015040. doi: 10.1074/mcp.R111.015040. 

      In OpenMS, TraML files can be generated from TSV or CSV files using the
      @ref OpenMS::TransitionTSVFile "TransitionTSVFile class" or the @ref
      TOPP_TargetedFileConverter "TargetedFileConverter TOPP Tool". For more information on the TSV format required by the TOPP tool, see
      see also the documentation of @ref
      OpenMS::TransitionTSVFile "TransitionTSVFile".

      @ingroup FileIO
  */
  class OPENMS_DLLAPI TraMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    ///Default constructor
    TraMLFile();
    ///Destructor
    ~TraMLFile() override;

    /**
        @brief Loads a map from a TraML file.

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String & filename, TargetedExperiment & id);

    /**
        @brief Stores a map in a TraML file.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String & filename, const TargetedExperiment & id) const;

    /**
        @brief Checks if a file is valid with respect to the mapping file and the controlled vocabulary.

        @param filename File name of the file to be checked.
        @param errors Errors during the validation are returned in this output parameter.
        @param warnings Warnings during the validation are returned in this output parameter.

        @exception Exception::FileNotFound is thrown if the file could not be opened
    */
    bool isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings);

private:

  };

} // namespace OpenMS


