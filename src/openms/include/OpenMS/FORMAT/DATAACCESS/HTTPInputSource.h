// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/InputStream.h>
#include <OpenMS/config.h>

#include <memory>
#include <string>

namespace OpenMS
{
  /**
    @brief Source for reading XML files from HTTP, HTTPS and FTP URLs, with libcurl.

    The response body is streamed to the parser; gzip and bzip2 compressed files (recognized by their magic bytes)
    are decompressed on the fly. Redirects are followed. Proxies, CA certificates etc. are libcurl's defaults, e.g.
    the environment variables https_proxy and CURL_CA_BUNDLE. Credentials can be part of the URL
    (ftp://user:password\@host/file).

    The counterpart of S3InputSource for s3:// URIs.

    @ingroup Format
  */
  class OPENMS_DLLAPI HTTPInputSource
  {
public:
    /**
      @brief The source of the URL @p url; no request is made yet.

      @throws Exception::IllegalArgument if @p url is no HTTP, HTTPS or FTP URL
    */
    explicit HTTPInputSource(const std::string& url);

    ~HTTPInputSource();

    /**
      @brief Starts the download and returns a stream of the (decompressed) content, for the parser.

      Read errors of the stream throw Exception::IOException; corrupt compressed data throws Exception::ParseError.

      @throws Exception::FileNotFound if the server reports that the file does not exist (HTTP 404 or 410, FTP 550)
      @throws Exception::IOException if the download fails otherwise
    */
    std::unique_ptr<InputStream> makeStream() const;

    /// Whether @p url is an HTTP, HTTPS or FTP URL (starts with "http://", "https://" or "ftp://", in any letter case)
    static bool isHTTPURI(const std::string& url);

    /**
      @brief The size of the file at @p url, in bytes, as reported by the server.

      @throws Exception::IllegalArgument, Exception::FileNotFound, Exception::IOException as makeStream(); also
      Exception::IOException if the server does not report the size
    */
    static Int64 objectSize(const std::string& url);

    /**
      @brief The bytes [@p offset, @p offset + @p length) of the file at @p url (fewer at its end).

      Uses a range request; from a server that ignores it, the bytes before @p offset are downloaded and discarded.

      @throws Exception::IllegalArgument, Exception::FileNotFound, Exception::IOException as makeStream()
    */
    static std::string readRange(const std::string& url, Int64 offset, Int64 length);

    HTTPInputSource(const HTTPInputSource&) = delete;
    HTTPInputSource& operator=(const HTTPInputSource&) = delete;

private:
    std::string url_;
  };

} // namespace OpenMS
