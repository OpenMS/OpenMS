// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/DATASTRUCTURES/ToolDescription.h>

namespace OpenMS
{
  /**
      @brief File adapter for ToolDescriptor files

      If a critical error occurs due to the missing functionality, Exception::NotImplemented is thrown.

      @ingroup FileIO
  */
  class OPENMS_DLLAPI ToolDescriptionFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    ///Default constructor
    ToolDescriptionFile();
    ///Destructor
    ~ToolDescriptionFile() override;

    /**
        @brief Loads a map from a ToolDescriptor file.

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String & filename, std::vector<Internal::ToolDescription> & tds);

    /**
        @brief Stores a map in a ToolDescriptor file.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String & filename, const std::vector<Internal::ToolDescription> & tds) const;

private:

  };

} // namespace OpenMS

