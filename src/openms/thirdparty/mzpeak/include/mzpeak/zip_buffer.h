/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <cstddef>
#include <vector>

#include "mzpeak/archive.h"

namespace MzPeak {

/**
 * RDR-22 — access files in a ZIP archive whose bytes live in memory.
 *
 * Mirrors the Rust reader's `from_buf`: the archive OWNS the byte buffer
 * (moved in at construction) and that buffer outlives every `File` it
 * hands out.  Like `Zip` (RDR-26), each `read_file`/`list` opens a fresh
 * libzip handle over a `zip_source_buffer` view of the owned bytes, so the
 * reader can keep several members open at once with independent read
 * positions.
 */
class ZipBuffer final : public Archive {
public:
  /**
   * Construct from an in-memory `.mzpeak` archive.  The bytes are moved in
   * and owned for the lifetime of the archive.
   */
  ZipBuffer(std::vector<std::byte>);

  /**
   * Close the archive.
   */
  ~ZipBuffer();

  /**
   * Retrieve a list of files in the zip archive.
   */
  std::vector<fs::path> list();

  /**
   * Open a file from within the zip archive for reading.
   *
   * NOTE: The path given must be one returned from the `list` method.
   */
  std::unique_ptr<File> read_file(const fs::path&);

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace MzPeak
