/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <variant>

#include "mzpeak/schema/cv.h"

namespace MzPeak::Util {
enum class Type : int;
}

namespace MzPeak::Schema::PSI {

/**
 * Children of the MS:1000518 type:
 *
 * Encoding type of binary data specifying the binary representation
 * and precision, e.g. 64-bit float.
 */
class DataType {
public:
  enum Type {
    /// MS:1000519
    ///
    /// Signed 32-bit little-endian integer.
    Int32,

    /// MS:1000522
    ///
    /// Signed 64-bit little-endian integer.
    Int64,

    /// MS:1000521
    ///
    /// 32-bit precision little-endian floating point conforming to
    /// IEEE-754.
    Float32,

    /// MS:1000523
    ///
    /// 64-bit precision little-endian floating point conforming to
    /// IEEE-754.
    Float64,

    /// MS:1001479
    ///
    /// Sequence of zero or more non-zero ASCII characters terminated by
    /// a single null (0) byte.
    ASCII,
  };

  // Internal storage type.
  using value_type = std::variant<Type, CV>;

  /**
   * Create a DataType from a CV value.
   */
  explicit DataType(const CV&);

  /**
   * Create a DataType from the given enum Type.
   */
  explicit DataType(Type);

  /**
   * Convert a DataType to a CV.
   */
  CV to_cv() const;

  /**
   * Return the internal Type if the CV mapped to a standard type.
   */
  std::optional<Type> data_type() const;

  /**
   * Convert a DataType to a Util::Type.
   */
  std::optional<Util::Type> as_type() const;

  /**
   * Ensure this DataType is treated like the given Util::Type.
   *
   * This is mandatory if the internal CV term is non-standard and
   * therefore can't be mapped to a Util::Type.
   */
  void as_type(Util::Type);

  /// Equality.
  bool operator==(const DataType& other) const;

  /// Less than (for sorting).
  bool operator<(const DataType& other) const;

private:
  value_type val_;
  std::optional<Util::Type> as_type_;
};

} // namespace MzPeak::Schema::PSI
