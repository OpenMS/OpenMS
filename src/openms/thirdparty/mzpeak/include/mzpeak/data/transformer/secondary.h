/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include "mzpeak/data/array_index.h"
#include "mzpeak/util/types.h"

namespace MzPeak::Data::Transformer::Secondary {

/**
 * Decoding class for secondary dimensions.
 */
template <Util::supported_type T> class Decoder final {
public:
  /// The type of result this function object returns.
  using result_type = std::variant<std::shared_ptr<arrow::Array>,
                                   std::pair<T, std::shared_ptr<arrow::Array>>,
                                   std::shared_ptr<std::vector<T>>>;
  /// Constructor.
  Decoder(const ArrayIndex::Dimension&);

  /**
   * Transform the given array.
   *
   * This array could come from the `chunk_values` column, or from the
   * `chunk_transform` column.
   */
  result_type operator()(int64_t, const std::shared_ptr<arrow::Array>&) const;

private:
  std::string dim_name_;
  std::optional<Schema::PSI::Transform> transform_;
};

/******************************************************************************/
template <Util::supported_type T>
Decoder<T>::Decoder(const ArrayIndex::Dimension& dim)
    : dim_name_(dim.name)
    , transform_(dim.transform)
{
}

/******************************************************************************/
template <Util::supported_type T>
Decoder<T>::result_type
Decoder<T>::operator()(int64_t, const std::shared_ptr<arrow::Array>& src) const
{
  if (!transform_.has_value()) return src;

  std::optional<Schema::PSI::Transform::Type> type = transform_.value().type();

  if (!type.has_value()) {
    std::string msg("while decoding " + dim_name_ + ": ");
    msg += "unknown transform method: " + transform_.value().to_cv().to_string();
    throw InvalidFormatError(msg);
  }

  switch (type.value()) {
  case Schema::PSI::Transform::ZeroIntensityTrim:
  case Schema::PSI::Transform::ZeroIntensityInterpolation:
    // Handled by the null decoding code.
    return src;

  case Schema::PSI::Transform::NumpressSLOF:
    return Util::Numpress::decode_slof_convert<T>(src);

  case Schema::PSI::Transform::NumpressPIC:
    return Util::Numpress::decode_pic_convert<T>(src);
  }

  std::unreachable();
}

} // namespace MzPeak::Data::Transformer::Secondary
