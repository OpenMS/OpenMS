/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <stdexcept>
#include <string>

namespace MzPeak {

/**
 * Exception base class for the MzPeak library.
 */
class Exception : public std::runtime_error {
public:
  /// Constructor.
  Exception(const std::string& msg)
      : std::runtime_error(msg)
  {
  }

  /// Destructor.
  ~Exception() = default;
};

/**
 * Error thrown when a JSON file could not be parsed.
 */
class JsonError final : public Exception {
public:
  /// Constructor.
  JsonError(const std::string& msg)
      : Exception(msg)
  {
  }

  /// Destructor.
  ~JsonError() = default;
};

/**
 * Thrown when an Arrow or Parquet error is encountered.
 */
class ParquetError final : public Exception {
public:
  /// Constructor.
  ParquetError(const std::string& msg)
      : Exception(msg)
  {
  }

  /// Destructor.
  ~ParquetError() = default;
};

/**
 * The mzPeak file is malformed and the error could not be recovered
 * from.
 */
class InvalidFormatError final : public Exception {
public:
  /// Constructor.
  InvalidFormatError(const std::string& msg)
      : Exception(msg)
  {
  }

  /// Destructor.
  ~InvalidFormatError() = default;
};

/**
 * A member's bytes do not hash to the SHA-512 the index records for it, so the
 * archive has been altered or damaged since it was written.  Only raised when
 * the caller asked for verification: it is a MAY in the specification, and
 * most archives predate the requirement entirely.
 */
class ChecksumError final : public Exception {
public:
  /// Constructor.
  ChecksumError(const std::string& msg)
      : Exception(msg)
  {
  }

  /// Destructor.
  ~ChecksumError() = default;
};

/**
 * Automatic decoding of signal data is only supported for standard
 * file layouts such as point and chunked.
 */
class UnknownLayoutError final : public Exception {
public:
  /// Constructor.
  UnknownLayoutError(const std::string& msg)
      : Exception(msg)
  {
  }

  /// Destructor.
  ~UnknownLayoutError() = default;
};

/**
 * Attempt to access an invalid iterator.
 */
class InvalidIteratorError final : public Exception {
public:
  /// Constructor.
  InvalidIteratorError(const std::string& msg)
      : Exception(msg)
  {
  }

  /// Destructor.
  ~InvalidIteratorError() = default;
};

/**
 * Thrown when an unexpected type is encountered at run time.
 */
class TypeError final : public Exception {
public:
  /// Constructor.
  TypeError(const std::string& msg)
      : Exception(msg)
  {
  }

  /// Destructor.
  ~TypeError() = default;
};

/**
 * Failed to allocate memory.
 */
class AllocationError final : public Exception {
public:
  /// Constructor.
  AllocationError(const std::string& msg)
      : Exception(msg)
  {
  }

  /// Destructor.
  ~AllocationError() = default;
};

} // namespace MzPeak
