// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <string>

namespace OpenMS
{
  /**
    @brief Convenience wrapper around @c NetworkGetRequest for one-shot HTTP/HTTPS downloads.

    Single-method utility used by tooling that needs to fetch a remote URL to disk without
    composing a @c NetworkGetRequest pipeline by hand. Synchronous (the call blocks until
    the transfer finishes or fails).

    @ingroup System
  */
  class OPENMS_DLLAPI Network
  {
  public:
    /**
      @brief Download @p url to @p download_folder and return when the file is fully written.

      Internally creates a @ref NetworkGetRequest with a 10-minute (@c 600 second) timeout
      and runs it synchronously. The destination filename is derived from the URL:
        - any @c "?query" and @c "#fragment" suffixes are stripped first;
        - the basename of the remaining path component is used (e.g. @c "https://host/a/b/c.tsv" → @c "c.tsv");
        - if the resulting basename is empty, @c "download" is used as a fallback.
      If a file with the same name already exists in @p download_folder, the suffixes
      @c ".0", @c ".1", @c ".2", ... are tried in order until an unused name is found —
      existing files are never overwritten.

      An empty @p download_folder is treated as @c "./" (the current working directory).
      No trailing-slash normalisation is performed on the directory string.

      Success path: writes the response bytes to the destination in binary mode and
      emits two @c OPENMS_LOG_INFO lines (download successful + final on-disk path).
      Failure paths throw:
        - download error (HTTP error or timeout in @c NetworkGetRequest),
        - opening the destination file fails (e.g. permission denied, missing folder),
        - writing the response bytes fails partway.
      A partially-written file is not cleaned up on write failure.

      @param[in] url             Source URL.
      @param[in] download_folder Destination directory; @c "" is treated as @c "./".

      @throws OpenMS::Exception::IOException on any of the failure paths described above.
    */
    static void downloadFile(const std::string& url, const std::string& download_folder);
  };
}
