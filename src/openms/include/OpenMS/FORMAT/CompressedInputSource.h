// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <xercesc/sax/InputSource.hpp>

namespace OpenMS
{
  /**
    @brief Xerces @c InputSource that opens a bzip2-, zip- or gzip-compressed file on demand.

    Resolves the file path the same way as @c xercesc::LocalFileInputSource
    (relative paths are completed against the current working directory),
    but @ref makeStream picks an OpenMS decompression stream
    (@ref Bzip2InputStream, @ref ZipInputStream or @ref GzipInputStream)
    based on the first two bytes of the file as passed to the
    constructor.

    Copy construction, copy assignment and default construction are
    deliberately suppressed.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI CompressedInputSource :
    public xercesc::InputSource
  {
public:
    /**
      @brief Construct from a UTF-8 @c String path.

      @param[in] file_path Path to the compressed file; relative paths
                           are completed against Xerces' current
                           working directory.
      @param[in] header    First two bytes of the file; @ref makeStream
                           uses them to pick the decompressor.
                           @c "BZ" selects bzip2, @c "PK" selects zip,
                           anything else (including a header shorter
                           than two characters, which is silently
                           replaced with two NUL bytes) selects gzip.
      @param[in] manager   Memory manager used for Xerces allocations.
                           Defaults to the platform-wide manager.
    */
    CompressedInputSource(const   String & file_path, const String & header, xercesc::MemoryManager * const manager = xercesc::XMLPlatformUtils::fgMemoryManager);

    /**
      @brief Construct from a wide-character (@c XMLCh) path.

      @param[in] file_path Path to the compressed file; relative paths
                           are completed against Xerces' current
                           working directory.
      @param[in] header    Same semantics as the @c header argument of the @c String overload of this constructor.
      @param[in] manager   Memory manager used for Xerces allocations.
                           Defaults to the platform-wide manager.
    */
    CompressedInputSource(const   XMLCh * const file_path, const String & header, xercesc::MemoryManager * const manager = xercesc::XMLPlatformUtils::fgMemoryManager);

    /// Destructor.
    ~CompressedInputSource() override;

    /**
      @brief Create the decompression input stream for this source.

      Inspects the @c header passed to the constructor and returns a
      newly-allocated @c BinInputStream of the matching type:
        - @c "BZ" -> @ref Bzip2InputStream
        - @c "PK" -> @ref ZipInputStream
        - anything else (including the synthetic NUL pair injected when
          the constructor received a header shorter than two characters)
          -> @ref GzipInputStream

      Ownership of the returned stream is transferred to the caller
      (typically the Xerces parser).

      @return A newly allocated decompression stream, or @c nullptr
              when the chosen stream fails to open the file.
    */
    xercesc::BinInputStream * makeStream() const override;

private:
    String head_; ///< Header bytes captured at construction; @ref makeStream uses the first two characters to pick a decompressor.

    /// Default construction is deliberately suppressed (declared but not defined).
    CompressedInputSource();
    /// Copy construction is deliberately suppressed (declared but not defined).
    CompressedInputSource(const CompressedInputSource & source);
    /// Assignment is deliberately suppressed (declared but not defined).
    CompressedInputSource & operator=(const CompressedInputSource & source);
  };

} // namespace OpenMS

