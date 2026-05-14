// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief Synchronous HTTP GET request backed by libcurl.

    Typical lifecycle: set the URL via @ref setUrl, optionally set a
    timeout via @ref setTimeout, call @ref run to perform the request,
    then read the result via @ref getResponse / @ref getResponseBinary
    or check @ref hasError / @ref getErrorString. The instance can be
    reused for further requests; each call to @ref run replaces the
    previously observed state. Redirects are followed automatically.

    Failures (libcurl initialisation, transport errors, HTTP status
    @c >= @c 400) are reported through @ref hasError and
    @ref getErrorString; @ref run does not throw on those conditions.

    Instances are not copyable.

    @ingroup System
  */
  class OPENMS_DLLAPI NetworkGetRequest
  {
  public:
    NetworkGetRequest();
    ~NetworkGetRequest();

    /// Set the URL to request. The URL is consumed by the next @ref run call; no validation is performed here.
    void setUrl(const std::string& url);

    /// Set the request timeout in seconds. The default of @c 0 leaves the timeout unset (request blocks until the server responds or libcurl gives up on its own).
    void setTimeout(int seconds);

    /**
      @brief Execute the GET request synchronously.

      Performs the request configured by the last @ref setUrl /
      @ref setTimeout call, blocking the caller until the response is
      complete or libcurl reports a timeout / error. Redirects are
      followed. The previous response and error state are cleared at
      the start of the call.

      The call sets @ref hasError to @c true when:
        - libcurl could not be initialised (error string
          @c "Failed to initialize libcurl"),
        - a transport-level libcurl error occurred (error string from
          @c curl_easy_strerror), or
        - the server responded with HTTP status @c >= @c 400 (error
          string @c "HTTP error <code>").

      No exception is thrown for any of these.
    */
    void run();

    /// Response body as a string (copy of @ref getResponseBinary).
    std::string getResponse() const;

    /// Raw response body. Reference is valid until the next @ref run call or destruction.
    const std::vector<char>& getResponseBinary() const;

    /// @c true if the last @ref run call did not produce a usable response (see @ref run for the conditions).
    bool hasError() const;

    /// Human-readable error message; empty when @ref hasError is @c false.
    std::string getErrorString() const;

  private:
    NetworkGetRequest(const NetworkGetRequest&) = delete;
    NetworkGetRequest& operator=(const NetworkGetRequest&) = delete;

    std::vector<char> response_bytes_;
    std::string url_;
    int timeout_ = 0;
    bool has_error_ = false;
    std::string error_string_;
  };
}
