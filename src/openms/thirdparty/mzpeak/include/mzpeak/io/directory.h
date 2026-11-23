/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include "mzpeak/io/archive.h"

namespace MzPeak::IO {

/**
 * Access files from a directory.
 */
class Directory final : public Archive {
public:
  /// Constructor.
  Directory(const fs::path& path)
      : path_(path.lexically_normal())
  {
  }

  /**
   * Return a list of files in the directory.
   *
   * NOTE: This is not a recursive listing.
   *
   * NOTE: Only files will be returned.
   */
  std::vector<fs::path> list() override;

  /**
   * Open a file for reading.
   *
   * The path *must* be a file name relative to the directory as
   * returned by the `list` method.
   */
  std::unique_ptr<File> read_file(const fs::path& path) override;

private:
  fs::path path_;
};

} // namespace MzPeak::IO
