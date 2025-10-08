// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/XMLFile.h>

#include <map>
#include <vector>

namespace OpenMS
{
  /**
      @brief Used to load and store PTMXML files

      This class is used to load and store documents that implement
      the schema of PTMXML files.

      @ingroup FileIO
  */
  class OPENMS_DLLAPI PTMXMLFile :
    public Internal::XMLFile
  {
public:
    /// Constructor
    PTMXMLFile();

    /**
        @brief Loads the information of a PTMXML file

        @param filename The name of the file
        @param ptm_informations the PTM information from the file are stored herein
        @throw FileNotFound is thrown if the given file could not be found
        @throw ParseError is thrown if the given file could not be parsed
        The information is read in and stored in the corresponding variables
    */
    void load(const String & filename, std::map<String, std::pair<String, String> > & ptm_informations);

    /**
        @brief Stores the data in an PTMXML file

        @throw UnableToCreateFile is thrown if the given filename could not be created

        The data is read in and stored in the file 'filename'.
    */
    void store(const String& filename, std::map<String, std::pair<String, String> > & ptm_informations) const;
  };

} // namespace OpenMS

