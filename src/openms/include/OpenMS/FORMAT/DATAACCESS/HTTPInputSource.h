// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Julianus Pfeuffer, Timo Sachsenberg $
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

    The resource is streamed to the parser; gzip and bzip2 compressed resources (recognized by their magic bytes) are
    decompressed on the fly. HTTP transfer compression (Content-Encoding) is decoded by libcurl.

    Redirects to HTTP, HTTPS and FTP URLs are followed. Proxies come from the usual environment variables
    (http_proxy, https_proxy, no_proxy). file:// URLs are accepted too, but XMLFile reads local files directly.

    @ingroup Format
  */
  class OPENMS_DLLAPI HTTPInputSource
  {
public:
    /**
      @brief A source for @p url.

      @param[in] url An HTTP, HTTPS, FTP or file URL

      @throws Exception::IllegalArgument if @p url is no such URL
    */
    explicit HTTPInputSource(const std::string& url);

    /**
      @brief Starts a transfer of the resource and returns a stream of its (decompressed) content, for the parser.

      Transfer errors of the stream throw Exception::IOException; corrupt compressed data throws Exception::ParseError.

      @throws Exception::FileNotFound if the server reports that the resource does not exist (e.g. HTTP 404)
      @throws Exception::IOException if the transfer fails (e.g. the host cannot be reached)
    */
    std::unique_ptr<InputStream> makeStream() const;

    /**
      @brief The size of the resource at @p url, in bytes, as reported by the server.

      @throws Exception::IllegalArgument, Exception::FileNotFound, Exception::IOException as the constructor and
      makeStream(); also Exception::IOException if the server does not report the size
    */
    static Int64 objectSize(const std::string& url);

    /**
      @brief The bytes [@p offset, @p offset + @p length) of the resource at @p url (fewer at its end), with a range request.

      From an HTTP server that ignores the range, the bytes before @p offset are received and discarded.

      @throws Exception::IllegalArgument, Exception::FileNotFound, Exception::IOException as the constructor and makeStream()
    */
    static std::string readRange(const std::string& url, Int64 offset, Int64 length);

    /// Whether @p url is an HTTP, HTTPS or FTP URL (starts with "http://", "https://" or "ftp://")
    static bool isURL(const std::string& url);

    /**
      @brief @p url without user info ("user:password@"), query and fragment, for log and error messages.

      These parts can hold credentials, e.g. ftp://user:password@host/run.mzML or the signature of a presigned URL.
      Applies to any "scheme://" URL (e.g. s3:// too). Local file names are returned unchanged, also when they
      contain '?', '#' or '@'.
    */
    static std::string redactURL(const std::string& url);

private:
    std::string url_;
  };

} // namespace OpenMS
