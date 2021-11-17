// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Authors: Simon Gene Gottlieb $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/ParamCTDFile.h>

namespace OpenMS
{

  /**
  @brief Load from JSON (in a Common Workflow Language (CWL) compatible way) into the Param class.
         Exports .cwl files.

  @param FlatHierarchy If set to true, all parameters will be listed on store without nesting.
                       The names will be expanded to include the nesting hierarchy.


      The JSON file must contain one top level mapping of param value names to actual values.
      These values can be one of the following types:
          - null
          - boolean
          - int
          - long
          - float
          - double
          - string
          - a CWL style file path ({ "class": "File", "path": "./myFolder/myFile.txt" })
          - an array of these

      param value names match the command line option without the leading '-'. Optionally the ':'
      can be replaced with a double underscore "__".
  @code
  {
      "in": {
          "class": "File",
          "path": "./myFolder/myFile.txt"
      },
      "out_prefix": "test_cwl_",
      "algorithm:threshold": 5,
      "algorithm:score_type": "ID"
  }
  @endcode

  Same file with "__" instead of ':' as the section separator.
  @code
  {
      "in": {
          "class": "File",
          "path": "./myFolder/myFile.txt"
      },
      "out_prefix": "test_cwl_",
      "algorithm__threshold": 5,
      "algorithm__score_type": "ID"
  }
  @endcode
  */
  template <bool FlatHierarchy>
  class OPENMS_DLLAPI ParamCWLFile
  {
  public:
    /**
       @brief Write CWL file

       @param filename The name of the file the param data structure should be stored in.
       @param param The param data structure that should be stored.
       @param ToolInfo Additional information about the Tool for which the param data should be stored.

       @exception std::ios::failure is thrown if the file could not be created
     */
    void store(const std::string& filename, const Param& param, const ToolInfo& tool_info) const;

    /**
       @brief Write CWL to output stream.

       @param os_ptr The stream to which the param data should be written.
       @param param The param data structure that should be writte to stream.
       @param tool_info Additional information about the Tool for which the param data should be written.
     */
    void writeCWLToStream(std::ostream* os_ptr, const Param& param, const ToolInfo& tool_info) const;

    /**
      @brief Read JSON file that is formatted in CWL conforming style.

      @param filename The file from where to read the Param object.
      @param param A param object with pre-filled defaults, which are updated by the values in the JSON file
      @return returns true if file was successfully loaded; false if an unknown (non-default) parameter name was encountered in the JSON file

      @exception Exception::FileNotFound is thrown if the file could not be found
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    static bool load(const std::string& filename, Param& param);
  };

  // Delaying instanciation of the ParamCWLFile classes, so we can write the function bodies into the .cpp file
  extern template class ParamCWLFile<true>;
  extern template class ParamCWLFile<false>;

} // namespace OpenMS
