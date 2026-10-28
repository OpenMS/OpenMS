// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Simon Gene Gottlieb $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ToolInfo.h>

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
       @brief Whether this build can write CWL files.

       Writing CWL needs the TDL library, which is only linked when OpenMS is configured with
       ENABLE_TDL=ON (the default is OFF). Without it, store() and writeCWLToStream() throw.
     */
    static bool isSupported();

    /**
       @brief Write CWL file

       @param[out] filename The name of the file the param data structure should be stored in.
       @param[in] param The param data structure that should be stored.
       @param[out] tool_info Additional information about the Tool for which the param data should be stored.

       @exception Exception::NotImplemented is thrown if this build has no CWL support (see isSupported()); no file is created then
       @exception Exception::UnableToCreateFile is thrown if the file could not be created
     */
    void store(const std::string& filename, const Param& param, const ToolInfo& tool_info) const;

    /**
       @brief Write CWL to output stream.

       @param[out] os_ptr The stream to which the param data should be written.
       @param[out] param The param data structure that should be writte to stream.
       @param[out] tool_info Additional information about the Tool for which the param data should be written.

       @exception Exception::NotImplemented is thrown if this build has no CWL support (see isSupported())
     */
    void writeCWLToStream(std::ostream* os_ptr, const Param& param, const ToolInfo& tool_info) const;
  };
} // namespace OpenMS
