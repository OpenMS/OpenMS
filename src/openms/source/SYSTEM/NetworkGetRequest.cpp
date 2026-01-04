// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/NetworkGetRequest.h>

#include <OpenMS/CONCEPT/LogStream.h>

#ifdef OPENMS_HAS_CURL
#include <curl/curl.h>
#else
// Fallback implementation without libcurl
#include <cstdio>
#include <cstring>
#include <array>
#ifndef OPENMS_WINDOWSPLATFORM
#include <sys/wait.h>
#endif
#endif

using namespace std;

namespace OpenMS
{

#ifdef OPENMS_HAS_CURL
  // Callback function for libcurl to write received data
  static size_t WriteCallback(void* contents, size_t size, size_t nmemb, void* userp)
  {
    size_t totalSize = size * nmemb;
    std::vector<char>* buffer = static_cast<std::vector<char>*>(userp);
    buffer->insert(buffer->end(), static_cast<char*>(contents), static_cast<char*>(contents) + totalSize);
    return totalSize;
  }
#endif

  NetworkGetRequest::NetworkGetRequest() :
    response_bytes_(),
    url_(),
    error_string_(),
    has_error_(false),
    timeout_ms_(10000)
  {
  }

  NetworkGetRequest::~NetworkGetRequest() = default;

  void NetworkGetRequest::setUrl(const String& url)
  {
    url_ = url;
  }

  String NetworkGetRequest::getResponse() const
  {
    if (response_bytes_.empty())
    {
      return "";
    }
    return String(response_bytes_.data(), response_bytes_.size());
  }

  const std::vector<char>& NetworkGetRequest::getResponseBinary() const
  {
    return response_bytes_;
  }

  bool NetworkGetRequest::hasError() const
  {
    return has_error_;
  }

  String NetworkGetRequest::getErrorString() const
  {
    return error_string_;
  }

  void NetworkGetRequest::setTimeout(int timeout_ms)
  {
    timeout_ms_ = timeout_ms;
  }

  void NetworkGetRequest::run()
  {
    response_bytes_.clear();
    error_string_.clear();
    has_error_ = false;

#ifdef OPENMS_HAS_CURL
    CURL* curl = curl_easy_init();
    if (!curl)
    {
      has_error_ = true;
      error_string_ = "Failed to initialize libcurl";
      return;
    }

    // Set URL
    curl_easy_setopt(curl, CURLOPT_URL, url_.c_str());

    // Set write callback
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response_bytes_);

    // Follow redirects
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);

    // Set timeout (in seconds)
    curl_easy_setopt(curl, CURLOPT_TIMEOUT_MS, (long)timeout_ms_);

    // Set user agent
    curl_easy_setopt(curl, CURLOPT_USERAGENT, "OpenMS");

    // Perform the request
    CURLcode res = curl_easy_perform(curl);

    if (res != CURLE_OK)
    {
      has_error_ = true;
      error_string_ = String("libcurl error: ") + curl_easy_strerror(res);
    }
    else
    {
      // Check HTTP response code
      long response_code = 0;
      curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &response_code);
      if (response_code >= 400)
      {
        has_error_ = true;
        error_string_ = String("HTTP error: ") + String(response_code);
      }
    }

    curl_easy_cleanup(curl);

#else
    // Fallback: No libcurl available - use platform-specific approach
    // For Unix/Mac: use curl command line if available
    #ifndef OPENMS_WINDOWSPLATFORM

    // First check if curl command exists
    FILE* check_curl = popen("command -v curl 2>/dev/null", "r");
    if (!check_curl)
    {
      has_error_ = true;
      error_string_ = "curl command not available. Please install curl or build OpenMS with libcurl support.";
      return;
    }
    char curl_path[256];
    bool curl_exists = (fgets(curl_path, sizeof(curl_path), check_curl) != nullptr);
    int close_status = pclose(check_curl);
    if (close_status == -1)
    {
      // Log warning but continue - this is just a check
      OPENMS_LOG_WARN << "pclose() failed while checking for curl, but continuing" << std::endl;
    }

    if (!curl_exists)
    {
      has_error_ = true;
      error_string_ = "curl command not found in PATH. Please install curl or build OpenMS with libcurl support.";
      return;
    }

    // SECURITY: Validate URL structure - must start with http:// or https://
    if (!url_.hasPrefix("http://") && !url_.hasPrefix("https://"))
    {
      has_error_ = true;
      error_string_ = "URL must start with http:// or https://";
      return;
    }

    // SECURITY: Validate URL doesn't contain dangerous shell metacharacters
    // Allow common URL characters including #, +, ~, @, ;, , etc. but reject shell metacharacters
    for (size_t i = 0; i < url_.size(); ++i)
    {
      char c = url_[i];
      // Reject dangerous shell metacharacters: backtick, pipe, dollar, redirect, etc.
      // Note: We use single quotes around URL in shell command, so single quote is also dangerous
      if (c == '`' || c == '|' || c == '$' || c == '<' || c == '>' ||
          c == '\'' || c == '"' || c == '\\' || c == '\n' || c == '\r')
      {
        has_error_ = true;
        error_string_ = String("URL contains dangerous shell metacharacter: ") + c;
        return;
      }
      // Optionally, ensure it's printable ASCII (0x20-0x7E) to catch other weird characters
      if (c < 0x20 || c > 0x7E)
      {
        has_error_ = true;
        error_string_ = String("URL contains non-printable character");
        return;
      }
    }

    // Build curl command with safe URL
    // Ensure timeout is at least 1 second to avoid 0 for sub-second values
    int timeout_seconds = std::max(1, timeout_ms_ / 1000);
    String command = "curl -s -m " + String(timeout_seconds) + " '" + url_ + "' 2>&1";
    std::array<char, 128> buffer;
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe)
    {
      has_error_ = true;
      error_string_ = "Failed to execute curl command";
      return;
    }

    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr)
    {
      size_t len = strlen(buffer.data());
      response_bytes_.insert(response_bytes_.end(), buffer.data(), buffer.data() + len);
    }

    // Decode pclose return value properly (Unix only - already inside #ifndef OPENMS_WINDOWSPLATFORM)
    int status = pclose(pipe);
    if (status == -1)
    {
      has_error_ = true;
      error_string_ = "pclose() failed";
    }
    else if (WIFEXITED(status))
    {
      int exit_code = WEXITSTATUS(status);
      if (exit_code != 0)
      {
        has_error_ = true;
        error_string_ = String("curl command failed with exit code: ") + String(exit_code);
      }
    }
    else if (WIFSIGNALED(status))
    {
      has_error_ = true;
      error_string_ = String("curl command terminated by signal: ") + String(WTERMSIG(status));
    }
    #else
    // Windows without curl: error
    has_error_ = true;
    error_string_ = "Network requests require libcurl on Windows. Please build OpenMS with libcurl support.";
    #endif
#endif
  }

} // namespace OpenMS
