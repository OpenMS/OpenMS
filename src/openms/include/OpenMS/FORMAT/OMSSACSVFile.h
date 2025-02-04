// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>


namespace OpenMS
{
  class String;
  /**
      @brief File adapter for OMSSACSV files.

  The files contain the results of the OMSSA algorithm in a comma separated manner. This file adapter is able to
      load the data from such a file into the structures of OpenMS

  @ingroup FileIO
*/
  class OPENMS_DLLAPI OMSSACSVFile
  {
public:

    /// Default constructor
    OMSSACSVFile();
    /// Destructor
    virtual ~OMSSACSVFile();

    /**
              @brief Loads a OMSSA file

              @param filename the name of the file to read from
              @param protein_identification the protein ProteinIdentification data
              @param id_data the peptide ids of the file

              @throw FileNotFound is thrown if the given file could not be found
              @throw ParseError is thrown if the given file could not be parsed
    */
    void load(const String & filename, ProteinIdentification & protein_identification, std::vector<PeptideIdentification> & id_data) const;
  };
} // namespace OpenMS

