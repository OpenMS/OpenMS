/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <variant>

#include "mzpeak/schema/cv.h"

namespace MzPeak::Schema::PSI {

/**
 * Encoding methods used in the `chunk_encoding` column.
 */
class ChunkEncoding final {
public:
  enum Type {
    /// MS:1000576
    ///
    /// Values are not encoded nor compressed.
    NoCompression,

    /// MS:1003089
    ///
    /// Data array compression using mantissa bit truncation, delta
    /// prediction and zlib compression.
    ///
    /// NOTE: Parquet takes care of everything except the delta
    /// encoding.
    Delta,

    /// MS:1002312
    ///
    /// Compression using MS-Numpress linear prediction compression.
    NumpressLinear,

    /// MS:1003826
    ///
    /// Coordinate grid encoding: the chunk stores integer indices into a
    /// parametric model carried in a sibling `<array>_grid` struct column,
    /// and no coordinate values at all.  Decoded by Transformer::Grid.
    Grid,
  };

  // Internal storage type.
  using value_type = std::variant<Type, CV>;

  /**
   * Construct an encoding type from a CV.
   */
  ChunkEncoding(const CV&);

  /**
   * Convert an encoding type to a CV.
   */
  CV to_cv() const;

  /**
   * Return the encoding type if it is known.
   */
  std::optional<Type> type() const;

private:
  value_type val_;
};

} // namespace MzPeak::Schema::PSI
