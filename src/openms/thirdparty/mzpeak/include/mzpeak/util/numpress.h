/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <arrow/array.h>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <vector>

namespace MzPeak::Util::Numpress {

/**
 * The kind of Numpress compression.
 */
enum Type {
  Linear,
  SLOF,
  PIC,
};

/**
 * Attempt to infer the numpress method type from the name of a
 * column.
 */
std::optional<Type> type_from_column_name(const std::string&);

/**
 * Return the number of elements that should be reserved in order to
 * decode `n` bytes encoding with `Type` `t`.
 */
std::size_t decoding_space_needed(std::size_t n, Type t);

/**
 * Possibly create a new vector and copy all of the elements from
 * `doubles` with a static cast to `T`.
 *
 * If `T` is `double` then return the input vector unchanged.
 */
template <typename T>
std::shared_ptr<std::vector<T>>
cast(const std::shared_ptr<std::vector<double>>& doubles)
{
  if constexpr (std::is_same_v<T, double>) {
    return doubles;
  } else {
    std::shared_ptr<std::vector<T>> result = std::make_shared<std::vector<T>>();
    result->reserve(doubles->size());

    for (const auto& d : *doubles) {
      result->push_back(static_cast<T>(d));
    }

    return result;
  }
}

/**
 * Decode a vector of bytes into a vector of doubles.
 *
 * The bytes need to be encoded using the MS-Numpress Linear encoding.
 */
void decode_linear(const std::vector<uint8_t>&, std::vector<double>&);

/**
 * Decode an arrow array of `uint8_t` values.
 */
std::shared_ptr<std::vector<double>>
decode_linear(const std::shared_ptr<arrow::Array>&);

/**
 * Decode and perform type conversion if necessary.
 */
template <typename T>
std::shared_ptr<std::vector<T>>
decode_linear_convert(const std::shared_ptr<arrow::Array>& src)
{
  std::shared_ptr<std::vector<double>> doubles = decode_linear(src);
  return cast<T>(doubles);
}

/**
 * Decode a vector of bytes into a vector of doubles.
 *
 * The bytes need to be encoded using the MS-Numpress short logged
 * float compression encoding.
 */
void decode_slof(const std::vector<uint8_t>&, std::vector<double>&);

/**
 * Decode an arrow array of `uint8_t` values.
 */
std::shared_ptr<std::vector<double>>
decode_slof(const std::shared_ptr<arrow::Array>&);

/**
 * Decode and perform type conversion if necessary.
 */
template <typename T>
std::shared_ptr<std::vector<T>>
decode_slof_convert(const std::shared_ptr<arrow::Array>& src)
{
  std::shared_ptr<std::vector<double>> doubles = decode_slof(src);
  return cast<T>(doubles);
}

/**
 * Decode a vector of bytes into a vector of doubles.
 *
 * The bytes need to be encoding using MS-Numpress positive integer
 * compression encoding.
 */
void decode_pic(const std::vector<uint8_t>&, std::vector<double>&);

/**
 * Decode an arrow array of `uint8_t` values.
 */
std::shared_ptr<std::vector<double>>
decode_pic(const std::shared_ptr<arrow::Array>&);

/**
 * Decode and perform type conversion if necessary.
 */
template <typename T>
std::shared_ptr<std::vector<T>>
decode_pic_convert(const std::shared_ptr<arrow::Array>& src)
{
  std::shared_ptr<std::vector<double>> doubles = decode_pic(src);
  return cast<T>(doubles);
}

} // namespace MzPeak::Util::Numpress
