/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <cmath>
#include <concepts>
#include <ranges>
#include <vector>

namespace MzPeak::Util {

/// Two types that can be converted back and forth.
template <class A, class B>
concept convertible_between = std::convertible_to<A, B> && std::convertible_to<B, A>;

/**
 * A type-safe wrapper for the NULL decoding regression model.
 */
template <typename T> class DeltaEstimator final {
public:
  /// Constructor.
  explicit DeltaEstimator(const std::vector<T> model)
      : model_(model)
  {
  }

  /**
   * Predict a delta given a real value.
   *
   * Type `T` is the data type for the internal regression model
   * parameters and type `U` is the value type stored in an Arrow
   * array.  An example of where they may differ is when the m/z
   * values are 32-bit floats since the m/z data model is always
   * stored as a vector of doubles.
   */
  template <typename U>
    requires convertible_between<T, U>
  U predict(U val) const;

private:
  template <typename U> friend class DeltaEstimator;
  std::vector<T> model_;
};

/******************************************************************************/
template <typename T>
template <typename U>
  requires convertible_between<T, U>
U DeltaEstimator<T>::predict(U base) const
{
  T acc{};
  T base_t = static_cast<T>(base);

  for (std::size_t i : std::views::iota(0ul, model_.size())) {
    if (i == 0ul) {
      acc += model_[i];
    } else {
      acc += model_[i] * std::pow(base_t, static_cast<T>(i));
    }
  }

  return static_cast<U>(acc);
}

} // namespace MzPeak::Util
