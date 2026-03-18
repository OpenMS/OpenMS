// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/NetworkGetRequest.h>
#include <OpenMS/SYSTEM/CurlInit.h>

#include <curl/curl.h>

using namespace std;

namespace OpenMS
{
  namespace
  {
    size_t writeCallback(char* ptr, size_t size, size_t nmemb, void* userdata)
    {
      auto* buf = static_cast<vector<char>*>(userdata);
      size_t total = size * nmemb;
      buf->insert(buf->end(), ptr, ptr + total);
      return total;
    }
  }

  NetworkGetRequest::NetworkGetRequest() = default;
  NetworkGetRequest::~NetworkGetRequest() = default;

  void NetworkGetRequest::setUrl(const string& url)
  {
    url_ = url;
  }

  void NetworkGetRequest::setTimeout(int seconds)
  {
    timeout_ = seconds;
  }

  void NetworkGetRequest::run()
  {
    CurlInit::ensure();

    has_error_ = false;
    error_string_.clear();
    response_bytes_.clear();

    CURL* curl = curl_easy_init();
    if (!curl)
    {
      has_error_ = true;
      error_string_ = "Failed to initialize libcurl";
      return;
    }

    curl_easy_setopt(curl, CURLOPT_URL, url_.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writeCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response_bytes_);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_NOSIGNAL, 1L); // thread-safe DNS timeout
    if (timeout_ > 0)
    {
      curl_easy_setopt(curl, CURLOPT_TIMEOUT, static_cast<long>(timeout_));
    }

    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK)
    {
      has_error_ = true;
      error_string_ = curl_easy_strerror(res);
    }
    else
    {
      long http_code = 0;
      curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
      if (http_code >= 400)
      {
        has_error_ = true;
        error_string_ = "HTTP error " + std::to_string(http_code);
      }
    }

    curl_easy_cleanup(curl);
  }

  string NetworkGetRequest::getResponse() const
  {
    return string(response_bytes_.begin(), response_bytes_.end());
  }

  const vector<char>& NetworkGetRequest::getResponseBinary() const
  {
    return response_bytes_;
  }

  bool NetworkGetRequest::hasError() const
  {
    return has_error_;
  }

  string NetworkGetRequest::getErrorString() const
  {
    return error_string_;
  }
}
