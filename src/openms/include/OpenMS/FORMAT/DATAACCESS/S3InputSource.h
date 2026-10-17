// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Julianus Pfeuffer, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/config.h>

#include <xercesc/sax/InputSource.hpp>

#include <memory>
#include <string>

namespace OpenMS
{
  /**
    @brief Xerces InputSource for reading XML files from Amazon S3 (s3://bucket/key), with Arrow's S3 filesystem.

    The object is streamed to the parser; gzip and bzip2 compressed objects (recognized by their magic bytes) are
    decompressed on the fly.

    Credentials, region and endpoint are those of Arrow's S3 filesystem: the AWS default credential chain
    (e.g. the environment variables AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY and AWS_DEFAULT_REGION, ~/.aws) and the
    query parameters of the URI, e.g. s3://bucket/key?region=eu-central-1, or
    s3://bucket/key?endpoint_override=localhost:9000&scheme=http for an S3-compatible server.

    Needs Arrow with S3 support (as in OpenMS' vcpkg build, arrow[s3]); see isAvailable().

    @ingroup Format
  */
  class OPENMS_DLLAPI S3InputSource :
    public xercesc::InputSource
  {
public:
    /**
      @brief Opens the S3 object @p uri.

      @param[in] uri An S3 URI (s3://bucket/key, optionally with query parameters, see above)

      @throws Exception::IllegalArgument if @p uri is no S3 URI or OpenMS has no S3 support (see isAvailable())
      @throws Exception::FileNotFound if the object cannot be opened (e.g. it does not exist, or access is denied)
    */
    explicit S3InputSource(const std::string& uri);

    ~S3InputSource() override;

    /**
      @brief A stream of the (decompressed) content of the object, for the parser.

      @return A new BinInputStream (the caller takes ownership), or nullptr if the object cannot be read
    */
    xercesc::BinInputStream* makeStream() const override;

    /// Whether @p uri is an S3 URI (starts with "s3://")
    static bool isS3URI(const std::string& uri);

    /// Whether OpenMS can read S3 URIs: it was built with an Arrow library that has S3 support
    static bool isAvailable();

    /**
      @brief The size of the S3 object @p uri, in bytes.

      @throws Exception::IllegalArgument, Exception::FileNotFound as the constructor
    */
    static Int64 objectSize(const std::string& uri);

    /**
      @brief The bytes [@p offset, @p offset + @p length) of the S3 object @p uri (fewer at its end).

      @throws Exception::IllegalArgument, Exception::FileNotFound as the constructor
      @throws Exception::IOException if the bytes cannot be read
    */
    static std::string readRange(const std::string& uri, Int64 offset, Int64 length);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_; ///< The opened object (Arrow types stay out of this header)
  };

} // namespace OpenMS
