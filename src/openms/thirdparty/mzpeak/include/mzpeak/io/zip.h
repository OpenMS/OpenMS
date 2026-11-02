/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include "mzpeak/io/archive.h"

namespace MzPeak::IO {

/**
 * Access files stored in a ZIP archive.
 */
class Zip final : public Archive {
public:
  /**
   * Open a zip archive at the given path.
   */
  Zip(const fs::path&);

  /**
   * Close the zip archive.
   */
  ~Zip();

  /**
   * Retrieve a list of files in the zip archive.
   */
  std::vector<fs::path> list() override;

  /**
   * Open a file from within the zip archive for reading.
   *
   * NOTE: The path given must be one returned from the `list`
   * method.
   */
  std::unique_ptr<File> read_file(const fs::path&) override;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace MzPeak::IO
