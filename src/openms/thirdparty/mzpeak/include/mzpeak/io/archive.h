/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <filesystem>
#include <vector>

#include "mzpeak/io/file.h"

namespace MzPeak::IO {

namespace fs = std::filesystem;

/**
 * An interface for read-only archives.
 */
class Archive {
public:
  /// Destructor.
  virtual ~Archive() {}

  /**
   * Return a list of files names.
   */
  virtual std::vector<fs::path> list() = 0;

  /**
   * Open a file and gain random read-only access.
   */
  virtual std::unique_ptr<File> read_file(const fs::path&) = 0;

protected:
  Archive() = default;
};

} // namespace MzPeak::IO
