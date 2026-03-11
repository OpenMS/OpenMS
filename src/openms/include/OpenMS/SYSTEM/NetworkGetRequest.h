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
  /// Synchronous HTTP GET request using libcurl.
  class OPENMS_DLLAPI NetworkGetRequest
  {
  public:
    NetworkGetRequest();
    ~NetworkGetRequest();

    /// set the URL to request
    void setUrl(const std::string& url);

    /// set the timeout in seconds (0 = no timeout)
    void setTimeout(int seconds);

    /// execute the GET request (blocks until complete or timeout)
    void run();

    /// returns the response as a string
    std::string getResponse() const;

    /// returns the raw response bytes
    const std::vector<char>& getResponseBinary() const;

    /// returns true if an error occurred during the query
    bool hasError() const;

    /// returns the error message
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
