// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Ruben Grünberg $
// $Authors: Ruben Grünberg $
// --------------------------------------------------------------------------

#pragma once

#include <string>
#include <OpenMS/DATASTRUCTURES/Param.h>

namespace OpenMS
{

  /**
   @brief A struct to pass information about the tool as one parameter
   */
  struct ToolInfo
  {
    std::string version_;
    std::string name_;
    std::string docurl_;
    std::string category_;
    std::string description_;
    std::vector<std::string> citations_;
  };

  /**
  @brief Serializes a Param class in paramCTD file format.
         Note: only storing is currently possible

*/
  class OPENMS_DLLAPI ParamCTDFile
  {
  public:
    ParamCTDFile() = default; ///Constructor
    
    ~ParamCTDFile() = default; ///Destructor

    /**
       @brief Write CTD file

       @param filename The name of the file the param data structure should be stored in.
       @param param The param data structure that should be stored.
       @param tool_info Additional information about the Tool for which the param data should be stored.

       @exception std::ios::failure is thrown if the file could not be created
     */
    void store(const std::string& filename, const Param& param, const ToolInfo& tool_info) const;

    /**
       @brief Write CTD to output stream.

       @param os_ptr The stream to which the param data should be written.
       @param param The param data structure that should be writte to stream.
       @param tool_info Additional information about the Tool for which the param data should be written.
     */
    void writeCTDToStream(std::ostream* os_ptr, const Param& param, const ToolInfo& tool_info) const;

  private:
    /**
      @brief Escapes certain characters in a string that are not allowed in XML
             Escaped characters are: & < > " '

      @param to_escape The string in which the characters should be escaped

      @returns The escaped string
     */
    static std::string escapeXML(const std::string& to_escape);

    /**
      @brief Replace all occurrences of a character in a string with a string

      @param replace_in The string in which the characters should be replaced.
      @param to_replace The character that should be replaced.
      @param replace_with The string the character should be replaced with.
     */
    static void replace(std::string& replace_in, char to_replace, const std::string& replace_with);

    const std::string schema_location_ = "/SCHEMAS/Param_1_8_0.xsd";
    const std::string schema_version_ = "1.8.0";
  };
}
