// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/NetworkGetRequest.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <cstdio>
#include <cstring>
#include <array>
#ifndef OPENMS_WINDOWSPLATFORM
#include <sys/wait.h>
#endif

using namespace std;

namespace OpenMS
{

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

#ifndef OPENMS_WINDOWSPLATFORM
    // Unix/Mac: use curl command line

    // First check if curl command exists
    FILE* check_curl = popen("command -v curl 2>/dev/null", "r");
    if (!check_curl)
    {
      has_error_ = true;
      error_string_ = "curl command not available. Please install curl.";
      return;
    }
    char curl_path[256];
    bool curl_exists = (fgets(curl_path, sizeof(curl_path), check_curl) != nullptr);
    int close_status = pclose(check_curl);
    if (close_status == -1)
    {
      // Log warning but continue - this is just a check
      OPENMS_LOG_WARN << "pclose() failed while checking for curl, but continuing\n";
    }

    if (!curl_exists)
    {
      has_error_ = true;
      error_string_ = "curl command not found in PATH. Please install curl.";
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
    for (Size i = 0; i < url_.size(); ++i)
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
      Size len = strlen(buffer.data());
      response_bytes_.insert(response_bytes_.end(), buffer.data(), buffer.data() + len);
    }

    // Decode pclose return value properly
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
    // Windows: error - no curl command line fallback
    has_error_ = true;
    error_string_ = "Network requests are not supported on Windows in this build.";
#endif
  }

} // namespace OpenMS
