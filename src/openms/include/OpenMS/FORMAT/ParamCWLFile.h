// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
  @brief Exports .cwl files.

        If Names include ':' it will be replaced with "__";
  */
  class OPENMS_DLLAPI ParamCWLFile
  {
  public:
    /**
       \brief If set to true, all parameters will be listed without nesting when writing the CWL File.
              The names will be expanded to include the nesting hierarchy.
     */
    bool flatHierarchy{};

    /**
       @brief Write CWL file

       @param filename The name of the file the param data structure should be stored in.
       @param param The param data structure that should be stored.
       @param tool_info Additional information about the Tool for which the param data should be stored.

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
  };
} // namespace OpenMS
