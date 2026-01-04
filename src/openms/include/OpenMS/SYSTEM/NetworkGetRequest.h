// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <string>
#include <vector>

namespace OpenMS
{

  /**
    @brief Simple HTTP GET request handler (replaces Qt networking)

    This class provides basic HTTP GET functionality without Qt dependencies.
    Uses libcurl if available, otherwise falls back to platform-specific implementation.

    @ingroup System
  */
  class OPENMS_DLLAPI NetworkGetRequest
  {
  public:

    /** @name Constructors and destructors
    */
    //@{
    /// default constructor
    NetworkGetRequest();

    /// destructor
    ~NetworkGetRequest();
    //@}

    // set request parameters
    void setUrl(const String& url);

    /// returns the response as string
    String getResponse() const;

    /// returns the response as binary data
    const std::vector<char>& getResponseBinary() const;

    /// returns true if an error occurred during the query
    bool hasError() const;

    /// returns the error message, if hasError can be used to check whether an error has occurred
    String getErrorString() const;

    /// execute the request (blocking)
    void run();

    /// set timeout in milliseconds (default: 10000ms = 10s)
    void setTimeout(int timeout_ms);

  private:
    /// assignment operator
    NetworkGetRequest& operator=(const NetworkGetRequest& rhs);
    /// copy constructor
    NetworkGetRequest(const NetworkGetRequest& rhs);

    std::vector<char> response_bytes_;
    String url_;
    String error_string_;
    bool has_error_;
    int timeout_ms_;
  };
}

