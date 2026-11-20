/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <string>

namespace MzPeak::Schema {

/**
 * How arrays are encoded in a Parquet file.
 */
enum class BufferFormat {
  /// This array is stored using the `point` layout. The `point`
  /// layout is all-or-nothing, every array **MUST** be in that
  /// format.
  Point,

  /// This array is part of the `chunked` layout.  It contains the
  /// starting value of the "main" axis for the chunk, inclusive.
  ChunkStart,

  /// This array is part of the `chunked` layout.  It contains the
  /// ending value of the "main" axis for the chunk, inclusive.  It
  /// should be less than the next chunk's `chunk_start` value.
  ChunkEnd,

  /// This array is part of the `chunked` layout.  It contains the
  /// list of values of the "main" axis bounded between the chunk's
  /// start and end point.
  ChunkValues,

  /// This array is part of the `chunked` layout.  It contains a CURIE
  /// indicating how `chunk_values` were encoded.
  ChunkEncoding,

  /// This array is part of the `chunked` layout.  It contains the
  /// list of values of an array other than the main axis for the
  /// chunk.
  ChunkSecondary,

  /// This array is part of the `chunked` layout.  It contains the
  /// list of raw byte contents of an array in the chunk that was
  /// opaquely transformed, e.g. using MS-Numpress.  It may be present
  /// in addition to a referenced `chunk_values` or `chunk_secondary`
  /// column.
  ChunkTransform,

};

/**
 * Convert a BufferFormat to a string.
 */
std::string buffer_format_to_string(BufferFormat);

/**
 * Parse an BufferFormat from a string.
 */
BufferFormat buffer_format_from_string(std::string_view);

} // namespace MzPeak::Schema
