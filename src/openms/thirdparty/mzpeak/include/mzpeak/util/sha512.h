/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <cstdint>
#include <string>
#include <string_view>

namespace MzPeak::Util {

/**
 * SHA-512 (FIPS 180-4), as the mzPeak specification requires for the
 * `checksum` of every file named in `mzpeak_index.json`.
 *
 * Implemented here rather than pulled from a crypto library because nothing
 * else in this project's dependency set (Arrow, Parquet, Boost, libzip)
 * provides it, and making OpenSSL a hard requirement of a mass-spectrometry
 * I/O library for one hash function is a poor trade for its consumers.  This
 * is an integrity check against bit rot and truncated transfers, which is
 * exactly what the specification scopes it to -- tamper resistance is deferred
 * there to a separate provenance mechanism, so nothing here is load bearing
 * for security.
 *
 * Correctness is pinned by the published FIPS vectors AND by the digests the
 * reference implementation wrote into its own fixture, so agreement is
 * cross-checked against the other implementation rather than self-asserted.
 */
class Sha512 final {
public:
  /// Constructor.
  Sha512();

  /// Absorb more of the message.  May be called any number of times; the
  /// result depends only on the concatenation, not on how it was split.
  void update(const void* data, std::size_t size);

  /// Finish and return the digest as 128 lowercase hex characters with no
  /// separators, which is the spelling the index requires.
  std::string hex_digest();

private:
  void compress(const unsigned char* block);

  std::uint64_t state_[8];
  std::uint64_t length_; ///< Message length in BYTES (< 2^61, so no 128-bit
                         ///< counter is needed for any real file).
  unsigned char buffer_[128];
  std::size_t buffered_;
};

/// Convenience: the hex digest of a contiguous byte range.
std::string sha512_hex(std::string_view bytes);

} // namespace MzPeak::Util
