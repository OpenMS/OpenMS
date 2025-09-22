// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/ToolDescription.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>
#include <vector>

namespace OpenMS
{
  /**
    @brief Load TOPP tool configurations from JSON file.

    This class provides functionality to load tool definitions from a JSON configuration file
    instead of hardcoding them in the source code. The JSON file should have the following structure:

    @code
    {
      "categories": {
        "cat_id": "Display Name",
        ...
      },
      "tools": {
        "ToolName": {
          "category": "cat_id",
          "requires_gui": true,
          "requires_openswath": true,
          "requires_parquet": true
        },
        ...
      }
    }
    @endcode

    The optional flags (requires_gui, requires_openswath, requires_parquet) indicate
    build-time dependencies that affect tool availability.
  */
  class OPENMS_DLLAPI ToolJSONFile
  {
  public:
    /**
      @brief Load tool configurations from JSON file.

      @param filename The JSON file to load tool configurations from
      @param tools Output vector of tool descriptions
      @param categories Output map of category IDs to display names
      @return true if file was successfully loaded, false otherwise

      @exception Exception::FileNotFound is thrown if the file could not be found
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    static bool load(const String& filename, 
                     std::vector<Internal::ToolDescription>& tools,
                     std::map<String, String>& categories);

    /**
      @brief Get the default path for the tools JSON configuration file.
      
      @return Path to the default tools.json file in the OpenMS data directory
    */
    static String getDefaultConfigPath();

  private:
    /// Not implemented (static class)
    ToolJSONFile();
  };

} // namespace OpenMS