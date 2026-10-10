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

      The content is serialised into an exclusively-created temporary next to @p filename and then
      renamed over the target, so a failure part way through (full disk, I/O error) leaves the
      previous, valid file untouched instead of truncating it. Directories and write-protected
      targets are rejected up front; symlinks are resolved so the link target is replaced rather
      than the link. The special filename "-" writes to stdout.

      On POSIX the final rename is atomic. On Windows it is atomic where the platform replaces an
      existing target in a single step; in the rare case it does not, a best-effort fallback removes
      the target and retries, and if that retry fails the freshly written content is preserved under
      the temporary name (reported in the exception message) rather than both copies being lost.
      Note that replacing by rename gives the target a new inode, so ownership, ACLs, extended
      attributes and hard links are not preserved (Qt's QSaveFile behaves the same way).

      @param[in] filename The filename where the param data structure should be stored.
      @param[in] param The Param class that should be stored in the file.

      @exception Exception::UnableToCreateFile is thrown if the target cannot be written: it is a
                 directory or write-protected, the temporary file cannot be created, the write or
                 flush fails, or the temporary cannot be renamed over the target.
    */
    void store(const std::string& filename, const Param& param) const;

    /**
      @brief Write XML to output stream.

      @param[out] os_ptr The stream where the param class should be written to.
      @param[in] param The Param class that should be written to the stream.
    */
    void writeXMLToStream(std::ostream* os_ptr, const Param& param) const;

    /**
      @brief Read XML file.

      @param[in] filename The file from where to read the Param object.
      @param[out] param The param object where the read data should be stored.

      @exception Exception::FileNotFound is thrown if the file could not be found
      @exception Exception::ParseError is thrown if an error occurs during parsing, including a
                 well-formed document that contains no @em PARAMETERS element (i.e. is not a
                 parameter file at all); @p param is cleared in that case.
    */
    void load(const std::string& filename, Param& param);
  };

} // namespace OpenMS


