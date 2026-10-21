// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/InputStream.h>
#include <string>

namespace OpenMS
{
  /** @brief Internal source selecting bzip2, ZIP or gzip decompression by file magic. */
  class OPENMS_DLLAPI CompressedInputSource
  {
  public:
    /// Relative UTF-8 paths are resolved against the current working directory.
    CompressedInputSource(const std::string& file_path, const std::string& header);
    CompressedInputSource(const char16_t* file_path, const std::string& header);
    /// Create a stream owned by the caller. BZ selects bzip2, PK ZIP, otherwise gzip.
    InputStream* makeStream() const;
    const std::string& getSystemId() const { return file_path_; }

    CompressedInputSource() = delete;
    CompressedInputSource(const CompressedInputSource&) = delete;
    CompressedInputSource& operator=(const CompressedInputSource&) = delete;
  private:
    std::string file_path_;
    std::string head_;
  };
}
