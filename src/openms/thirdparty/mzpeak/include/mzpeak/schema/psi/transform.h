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
 * Possible transformation types.
 */
class Transform final {
public:
  /**
   * Known transform types.
   */
  enum Type {
    /// MS:1003901
    ///
    /// Apply an algorithm to remove excess zero intensity value data
    /// points from a spectrum. Data may be retained for
    /// interperatbility such as retaining only zeros that flank
    /// non-zero intensity value data points from a profile spectrum.
    ZeroIntensityTrim,

    /// MS:1003902
    ///
    /// A zero intensity point trimming algorithm that interpolates
    /// the m/z coordinate values from the local data or an estimated
    /// model.
    ZeroIntensityInterpolation,

    /// MS:1002314
    ///
    /// Compression using MS-Numpress short logged float compression.
    NumpressSLOF,

    /// MS:1002313
    ///
    /// Compression using MS-Numpress positive integer compression.
    NumpressPIC,
  };

  /// The transform value can be a type as described by the Type enum,
  /// or a raw curie when the value was outside of the enum.
  using value_type = std::variant<Type, CV>;

  /// Constructor from a CV value..
  Transform(const CV&);

  /// Constructor from a `Transform::Type`.
  Transform(Type);

  /// Destructor.
  ~Transform() = default;

  /**
   * Return the Type for this Transform if it was decoded
   * successfully.  Otherwise return the raw CV term.
   */
  value_type value() const;

  /**
   * Return the Transform::Type if it is set.
   */
  std::optional<Type> type() const;

  /**
   * Convert this Transform to a CV term.
   */
  CV to_cv() const;

  /**
   * Return `true` if this transformation needs a delta model.
   */
  bool needs_delta_model() const noexcept;

  /// Equality.
  bool operator==(const Transform& other) const { return val_ == other.val_; }

private:
  value_type val_;
};

} // namespace MzPeak::Schema::PSI
