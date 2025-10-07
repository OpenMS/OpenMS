// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

namespace OpenMS
{
  /**
    @brief The file pendant of the Param class used to load and store the param
           datastructure as paramXML (i.e. INI files).

    A documented schema for this format can be found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS

  */
  class OPENMS_DLLAPI ParamXMLFile :
    public Internal::XMLFile
  {
public:
    /// Constructor.
    ParamXMLFile();

    /**
      @brief Write XML file.

      @param filename The filename where the param data structure should be stored.
      @param param The Param class that should be stored in the file.

      @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String& filename, const Param& param) const;

    /**
      @brief Write XML to output stream.

      @param os_ptr The stream where the param class should be written to.
      @param param The Param class that should be written to the stream.
    */
    void writeXMLToStream(std::ostream* os_ptr, const Param& param) const;

    /**
      @brief Read XML file.

      @param filename The file from where to read the Param object.
      @param param The param object where the read data should be stored.

      @exception Exception::FileNotFound is thrown if the file could not be found
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, Param& param);
  };

} // namespace OpenMS


