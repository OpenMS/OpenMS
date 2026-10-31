/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

#include <filesystem>

#include "mzpeak/index.h"

namespace MzPeak {

namespace fs = std::filesystem;

/**
 * What to check while opening, beyond what reading the archive requires.
 */
enum class Validate {
  /// Open without reading anything that decoding does not need.  The default,
  /// because verification costs a full pass over every indexed member.
  Nothing,

  /// Recompute every SHA-512 the index records and refuse the archive if one
  /// disagrees.  An archive that records NO digests -- anything written before
  /// the specification required them -- opens without complaint; use
  /// `Index::verify_checksums()` directly to tell that case apart from a fully
  /// verified one.
  Checksums,
};

/**
 * Open a MzPeak file for reading.
 *
 * @throws ChecksumError when @p validate is `Validate::Checksums` and a
 *         member's bytes do not match the digest the index records for it.
 */
MzPeak::Index open(const fs::path&, Validate validate = Validate::Nothing);

} // namespace MzPeak
