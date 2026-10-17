// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Julianus Pfeuffer, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_S3

#include <xercesc/sax/InputSource.hpp>
#include <xercesc/util/BinInputStream.hpp>

#include <memory>
#include <string>

namespace arrow::io { class InputStream; }
namespace arrow::fs { class FileSystem; }

namespace OpenMS
{
  /**
    @brief Bounded, sequential Arrow-to-Xerces stream with gzip/bzip2 magic detection.

    Owns and closes the supplied stream. Positions refer to decompressed bytes.
    Arrow failures, including truncated compressed input, become Exception::IOException.
    @ingroup FileIO
  */
  class OPENMS_DLLAPI ArrowInputStream : public xercesc::BinInputStream
  {
  public:
    /**
      @brief Take ownership of an Arrow input stream.
      @param[in] stream Input stream, positioned at its beginning.
      @param[in] filesystem Optional filesystem kept alive until the stream is closed.
    */
    explicit ArrowInputStream(std::shared_ptr<arrow::io::InputStream> stream,
                              std::shared_ptr<arrow::fs::FileSystem> filesystem = {});
    ~ArrowInputStream() override;
    ArrowInputStream(const ArrowInputStream&) = delete;
    ArrowInputStream& operator=(const ArrowInputStream&) = delete;

    /// Current position in decompressed bytes.
    XMLFilePos curPos() const override;
    /**
      @brief Read at most 64 KiB, returning zero only at EOF or for a zero-size request.
      @param[out] toFill Destination buffer.
      @param[in] maxToRead Destination capacity.
      @return Number of bytes read (partial reads are permitted).
    */
    XMLSize_t readBytes(XMLByte* const toFill, XMLSize_t maxToRead) override;
    /// No externally supplied content type.
    const XMLCh* getContentType() const override;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
  };

  /**
    @brief Experimental streaming XML input from an Arrow S3 URI.

    Arrow's URI factory preserves its default credential chain and supported URI
    options (for example region and endpoint). Temporary credentials can be supplied
    through Arrow's standard environment/credential providers. Anonymous access and
    session_token query parameters are not supported by this URI-based interface.
    gzip and bzip2 are detected by content, not by filename.

    The URI factory initializes S3 safely on concurrent first use. This reader
    borrows that process-wide lifecycle and never finalizes S3. The host must
    coordinate shutdown: quiesce all Arrow S3 users, release their resources, then
    call arrow::fs::EnsureS3Finalized() before process teardown. Finalization is
    irreversible and must not be deferred to a reader destructor or atexit handler.
    @warning Prefer Arrow's standard credential providers to URI credentials.
    Userinfo, query and fragment are removed from parser identifiers and document
    provenance; the original URI is retained only for opening the resource.
    @ingroup FileIO
  */
  class OPENMS_DLLAPI S3InputSource : public xercesc::InputSource
  {
  public:
    /** @brief Construct a source without making a network request.
        @param[in] uri Arrow S3 URI, including any Arrow-supported query options.
        @throws Exception::ParseError if the scheme is not s3://.
    */
    explicit S3InputSource(const std::string& uri);
    /// Open a fresh stream. The caller (normally Xerces) owns the result.
    xercesc::BinInputStream* makeStream() const override;

  private:
    std::string uri_;
  };
}

#endif
