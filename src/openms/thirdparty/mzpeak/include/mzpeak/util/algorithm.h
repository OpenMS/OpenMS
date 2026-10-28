/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

#include <arrow/array.h>
#include <arrow/builder.h>
#include <boost/math/statistics/univariate_statistics.hpp>
#include <functional>
#include <iterator>
#include <ranges>
#include <utility>
#include <vector>

#include "mzpeak/util/decoders.h"
#include "mzpeak/util/types.h"

namespace MzPeak::Util::Algorithm {

/**
 * Compute all runs of true values in the given vector.  For each run,
 * invoke the given function with index to the first and last true
 * value in the current span.
 */
template <typename Fn, typename... Args>
void spans(const std::vector<bool>& bins, Fn&& func, Args&&... args)
{
  auto i = std::ranges::find(bins.begin(), bins.end(), true);
  auto j = bins.end();

  for (; i != bins.end(); i = std::ranges::find(j, bins.end(), true)) {
    j = std::ranges::find(i + 1, bins.end(), false);
    auto first = std::distance(bins.begin(), i);
    auto last = std::distance(bins.begin(), j - 1);
    std::invoke(func, static_cast<std::size_t>(first), static_cast<size_t>(last),
                args...);
  }
}

/**
 * Return a vector of deltas computed from the given input vector.
 */
template <typename T> std::vector<T> deltas(const std::vector<T>& values)
{
  std::vector<T> result;

  if (values.size() > 1) {
    result.reserve(values.size() - 1);
  }

  for (std::size_t i : std::views::iota(1ul, values.size())) {
    result.push_back(values[i] - values[i - 1]);
  }

  return result;
}

/**
 * Compute the "median-below-median".
 *
 * That is, compute deltas for the input vector, then return the
 * median value from the deltas that are themselves below their
 * median.
 */
template <typename T> T median_delta(const std::vector<T>& values, T or_else)
{
  std::vector<T> ds(deltas(values));
  if (ds.empty()) return or_else;

  // NOTE: Explicit sorting is not needed to compute the median.
  T median = boost::math::statistics::median(ds.begin(), ds.end());

  auto [first, last] =
      std::ranges::remove_if(ds, [&median](auto& v) { return v > median; });

  ds.erase(first, last);
  if (ds.empty()) return or_else;

  return boost::math::statistics::median(ds.begin(), ds.end());
}

/**
 * Intersect an absolute row range with one record batch.
 *
 * @param range_offset,range_length  the range, in rows relative to the start of the
 * row group.
 * @param batch_start                the batch's first row within the row group.
 * @param batch_rows                 the batch's row count.
 * @return `{offset, length}` RELATIVE to the batch; length is 0 when the two
 *         do not overlap.
 *
 * Extracted so it can be tested directly: the previous inline form computed
 * `range_length - batch_start`, subtracting an absolute position from a length.
 * That is correct only for the first batch (`batch_start == 0`) and negative
 * afterwards, so no single-batch fixture could catch it.
 */
constexpr std::pair<int64_t, int64_t> intersect_range(int64_t range_offset,
                                                      int64_t range_length,
                                                      int64_t batch_start,
                                                      int64_t batch_rows)
{
  const int64_t begin = std::max(range_offset, batch_start);
  const int64_t end =
      std::min(range_offset + range_length, batch_start + batch_rows);
  if (end <= begin) return {0, 0};
  return {begin - batch_start, end - begin};
}

/**
 * Decode an Arrow array that was encoded with delta encoding.
 *
 * Parameters:
 *
 *   - start: The excluded starting value.
 *
 *   - src: The array to decode.
 *
 * Returns a newly allocated array containing the decoded values.
 * Null values transferred from the source array to the returned array
 * unchanged.
 */
template <Type T>
std::shared_ptr<arrow::Array>
null_delta_decode(typename type_traits<T>::value_type start,
                  const std::shared_ptr<arrow::Array>& src)
{
  using ValueType = typename type_traits<T>::value_type;
  ValueType zero = {};

  std::optional<ValueType> last(start);
  int64_t length = src->length();

  using Builder = type_traits<T>::builder_type;
  Builder builder;

  // N.B.: "The start point is *excluded* from the chunk-values array."
  if (!builder.Reserve(length + 1).ok()) {
    throw AllocationError("unable to decode delta encoded array");
  }

  using ArrayType = typename type_traits<T>::array_type;
  std::shared_ptr<ArrayType> casted = std::static_pointer_cast<ArrayType>(src);

  auto append = [&](const std::optional<ValueType>& v) -> void {
    arrow::Status status;

    if (v.has_value()) {
      status = builder.Append(v.value());
    } else {
      status = builder.AppendNull();
    }

    if (!status.ok()) {
      throw AllocationError("unable to decode delta encoded element");
    }
  };

  // N.B.: "The start point is *excluded* from the chunk-values array."
  //
  // If the first value in the list is NULL then the initial starting
  // value is ignored and the delta starts at 0.  However, if the
  // first *two* values are NULL then we still treat the starting
  // point as 0 but add the starting value as the first element in the
  // results.
  //
  // This is due to delta encoding coming after null marking.
  if (length > 0 && casted->IsNull(0)) {
    if (length > 1 && casted->IsNull(1)) {
      append({start});
    }

    last = {};
  } else {
    append({start});
  }

  for (int64_t index : std::views::iota(0, length)) {
    if (casted->IsValid(index)) {
      ValueType delta = Decoders::unsafe_array_value<T>(casted, index);
      last = last.value_or(zero) + delta;
      append(last);
    } else {
      last = {};
      append(last);
    }
  }

  std::shared_ptr<arrow::Array> result;

  if (!builder.Finish(&result).ok()) {
    throw AllocationError("failed to decode delta encoded array");
  }

  return result;
}

} // namespace MzPeak::Util::Algorithm
