// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/FileTypes.h>

#include <cstddef>
#include <string>

namespace OpenMS
{
  /**
    @brief Filename-only format recognition and extension manipulation.

    These operations preserve FileHandler's filename rules and never inspect file
    contents or require the path to exist. Extensions are resolved through the
    FileTypes registry, so registered aliases (e.g. '.fa' for FASTA) are recognized
    alongside the preferred extension. Both slash styles are recognized on
    every platform. FileHandler retains its existing forwarding entry points.
    @ingroup FileIO
  */
  class OPENMS_DLLAPI FileNameUtils
  {
  public:
    /** @brief Determine a file type from its name, including compression suffixes.
        @param[in] filename Filename to inspect.
        @return Recognized type, or FileTypes::UNKNOWN.
    */
    static FileTypes::Type getTypeByFileName(const std::string& filename);

    /** @brief Accept the expected extension or an unknown extension.
        @param[in] filename Filename to inspect.
        @param[in] type Expected type.
    */
    static bool hasValidExtension(const std::string& filename, FileTypes::Type type);

    /** @brief Remove the extension using the established FileHandler rules.
        @param[in] filename Filename whose extension should be removed.
        @return Filename with its recognized extension removed, or the last
        basename suffix removed for an unknown type.
    */
    static std::string stripExtension(const std::string& filename);

    /** @brief Strip the existing extension and append the requested extension.
        @param[in] filename Original filename.
        @param[in] new_type Type of the new extension.
        @return Filename with the replacement extension.
    */
    static std::string swapExtension(const std::string& filename, FileTypes::Type new_type);

  private:
    /** @brief Match the longest known extension at the end of @p filename.

        Every dot-delimited suffix of the basename is looked up in the FileTypes registry (preferred
        extensions and aliases alike) and the longest match wins, so '.pep.xml' is preferred over
        '.xml'. A trailing compression suffix is seen through and counted as part of the extension.

        @param[in] filename Filename to inspect.
        @param[out] ext_start Offset of the '.' that starts the extension, or npos if the basename has
        no extension at all. Set even when the type is UNKNOWN, in which case it points at the last dot.
        @return Recognized type, or FileTypes::UNKNOWN.
    */
    static FileTypes::Type matchExtension_(const std::string& filename, size_t& ext_start);
  };
}
