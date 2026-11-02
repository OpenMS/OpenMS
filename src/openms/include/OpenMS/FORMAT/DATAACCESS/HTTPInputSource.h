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
    decompressed on the fly. HTTP transfer compression is requested as gzip (Content-Encoding) and decompressed the
    same way, as the stream is read, so a highly compressed response is never held in memory at once. A compressed
    resource (e.g. .mzML.gz) sent with transfer compression is decompressed in both layers.

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

      For HTTP and HTTPS, from the response to a request of the first byte, not a HEAD request: a presigned URL (e.g. of
      S3) is signed for one request method, usually GET.

      @param[in] url An HTTP, HTTPS, FTP or file URL

      @throws Exception::IllegalArgument, Exception::FileNotFound, Exception::IOException as the constructor and
      makeStream(); also Exception::IOException if the server does not report the size
    */
    static Int64 objectSize(const std::string& url);

    /**
      @brief The bytes [@p offset, @p offset + @p length) of the resource at @p url (fewer at its end), with a range request.

      From an HTTP server that ignores the range, the bytes before @p offset are received and discarded.

      @param[in] url An HTTP, HTTPS, FTP or file URL
      @param[in] offset The first byte to read
      @param[in] length The number of bytes to read

      @throws Exception::IllegalArgument, Exception::FileNotFound, Exception::IOException as the constructor and makeStream()
    */
    static std::string readRange(const std::string& url, Int64 offset, Int64 length);

    /**
      @brief Whether @p url is an HTTP, HTTPS or FTP URL (starts with "http://", "https://" or "ftp://")

      @param[in] url A URL or local file name
    */
    static bool isURL(const std::string& url);

    /**
      @brief @p url without user info ("user:password@"), query and fragment, for log and error messages.

      These parts can hold credentials, e.g. ftp://user:password@host/run.mzML or the signature of a presigned URL.
      Applies to any "scheme://" URL (e.g. s3:// too). Local file names are returned unchanged, also when they
      contain '?', '#' or '@'.

      @param[in] url A URL or local file name
    */
    static std::string redactURL(const std::string& url);

    /**
      @brief @p text with every "scheme://" URL in it redacted as by redactURL(), for diagnostics of other libraries.

      Arrow's or libcurl's messages can repeat the URL they were given, credentials included. A URL ends at
      whitespace, a quote or an angle bracket.

      @param[in] text A message that may contain URLs
    */
    static std::string redactURLs(const std::string& text);

    /**
      @brief @p text redacted as by redactURLs(text), and also without the parts of @p url that redactURL() removes.

      A library can also repeat a part of the URL it was given on its own, e.g. Arrow's "Invalid S3 connection scheme
      '...'" for the value of the query parameter "scheme". The user name, the password, each query value and the
      fragment of @p url (as given and percent-decoded) are replaced by "***" wherever they are a whole word of @p text,
      i.e. not part of a longer word.

      @param[in] text A message about @p url
      @param[in] url The URL the message is about (local file names have no such parts)
    */
    static std::string redactURLs(const std::string& text, const std::string& url);

private:
    std::string url_;
  };

} // namespace OpenMS
