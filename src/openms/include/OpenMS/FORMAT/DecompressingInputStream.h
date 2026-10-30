// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Julianus Pfeuffer, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/InputStream.h>

#include <memory>
#include <optional>
#include <string>

namespace OpenMS
{
  /**
    @brief Decompresses gzip or bzip2 data read from another InputStream, e.g. a remote object.

    Unlike GzipInputStream and Bzip2InputStream, which open local files, this decompresses an arbitrary stream on the
    fly. Concatenated members (as written by pigz or pbzip2) are read in sequence. Bytes after the last member that do
    not start another member are ignored, as gzread() does.

    Errors of the source stream propagate. Corrupt or truncated compressed data throws Exception::ParseError.
  */
  class OPENMS_DLLAPI DecompressingInputStream :
    public InputStream
  {
public:
    enum class Format { GZIP, BZIP2 };

    /// The format whose magic bytes start @p head (1F 8B: gzip, "BZh": bzip2), if any
    static std::optional<Format> detect(const unsigned char* head, Size size);

    /**
      @brief Decompresses @p source, which holds data in @p format.

      @param[in] source The compressed stream; it is read from its current position
      @param[in] format The compression format of @p source
      @param[in] name The name of the data, for error messages (e.g. its URI)
    */
    DecompressingInputStream(std::unique_ptr<InputStream> source, Format format, const std::string& name);

    ~DecompressingInputStream() override;

    /// The number of decompressed bytes read so far
    UInt64 curPos() const override;

    /// Reads up to @p max_to_read decompressed bytes into @p to_fill; 0 at the end of the data
    Size readBytes(unsigned char* to_fill, Size max_to_read) override;

    /// Cancels a blocking read of the source stream
    void cancel() override;

    DecompressingInputStream(const DecompressingInputStream&) = delete;
    DecompressingInputStream& operator=(const DecompressingInputStream&) = delete;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_; ///< zlib/bzip2 state (their headers stay out of this one)
  };

} // namespace OpenMS
