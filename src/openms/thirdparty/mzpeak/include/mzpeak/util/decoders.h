/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

#include <arrow/array.h>
#include <functional>
#include <memory>
#include <ranges>

#include "mzpeak/schema/group.h"
#include "mzpeak/util/compat.h" // IWYU pragma: keep
#include "mzpeak/util/types.h"

namespace MzPeak::Util::Decoders {

/******************************************************************************/
/**
 * `C` is either a container holding values of `V`, or `C` and `V` are
 * both scalar values of the same type.
 */
template <typename C, typename V>
concept scalar_or_container_of =
    std::same_as<std::remove_cvref_t<C>, std::remove_cvref_t<V>> ||
    (std::ranges::range<C> && std::convertible_to<std::ranges::range_value_t<C>, V>);

/******************************************************************************/
/**
 * `T` is a type that has a `decode` function that can decode values
 * from an `arrow::Array` and place the result in `R`.  The `R` type
 * can be a container or scalar value.
 */
template <typename T, typename R>
concept from_arrow_array =
    requires { typename T::value_type; } &&
    scalar_or_container_of<R, typename T::value_type> &&
    requires(T t, const std::shared_ptr<arrow::Array>& a, R& r) {
      { t.decode(a, r) } -> std::same_as<void>;
    };

/******************************************************************************/
/**
 * Read a value from an array without checking bounds or if it is NULL.
 *
 * This is needed because Parquet/Arrow uses std::string_view for byte
 * arrays, but that means that the original arrow array needs to
 * remain resident in memory.  Therefore we need to copy the memory
 * referenced by a std::string_view into a std::string.
 */
template <Type T>
type_traits<T>::value_type
unsafe_array_value(const std::shared_ptr<typename type_traits<T>::array_type>& ary,
                   int64_t index)
{
  using A = type_traits<T>::array_type;

  if constexpr (std::is_same_v<A, arrow::StringArray>) {
    return ary->GetString(index);
  } else {
    return ary->Value(index);
  }
}

/******************************************************************************/
/**
 * If the given array is a "list of lists" then visit each element of
 * the outer list.  The given function is called on non-null elements
 * and given the index to the list element.
 *
 * Returns the length of the outer list.
 */
template <typename F> int64_t visit(const std::shared_ptr<arrow::Array>& ary, F f)
{
  auto go = [&f]<typename L>(const std::shared_ptr<L>& list) -> int64_t {
    for (int64_t index : std::views::iota(0, list->length())) {
      if (list->IsValid(index)) {
        std::invoke(f, index, list->value_slice(index));
      }
    }

    return list->length();
  };

  auto type = ary->type_id();

  if (type == arrow::Type::LIST) {
    return go(std::static_pointer_cast<arrow::ListArray>(ary));
  } else if (type == arrow::Type::FIXED_SIZE_LIST) {
    return go(std::static_pointer_cast<arrow::FixedSizeListArray>(ary));
  } else if (type == arrow::Type::LARGE_LIST) {
    return go(std::static_pointer_cast<arrow::LargeListArray>(ary));
  } else if (type == arrow::Type::LIST_VIEW) {
    return go(std::static_pointer_cast<arrow::ListViewArray>(ary));
  } else if (type == arrow::Type::LARGE_LIST_VIEW) {
    return go(std::static_pointer_cast<arrow::LargeListViewArray>(ary));
  } else {
    std::string msg("expected an arrow list array but found: ");
    msg += ary->type()->name();
    throw TypeError(msg);
  }
}

/******************************************************************************/
/**
 * Try to figure out how many bytes should be reserved to decode the
 * given array.
 */
std::size_t guess_array_length(const Schema::Column&,
                               const std::shared_ptr<arrow::Array>&);

/******************************************************************************/
/**
 * A NULL decoder that always skips NULL values.
 */
template <typename T> struct NullSkip final {
  /// Skip this null;
  std::optional<T> operator()(int64_t) { return std::nullopt; }
};

/******************************************************************************/
/**
 * A NULL decoder that replaces NULL values with zero.
 */
template <typename T> struct NullToZero final {
  /// Replace the given NULL with zero.
  std::optional<T> operator()(int64_t)
  {
    T zero{};
    return zero;
  }
};

/******************************************************************************/
template <typename Child> struct Helper {

  // Helper function to push a value into a container, or when `V` and
  // `C` are both the same type do assignment.
  template <typename V, typename C>
    requires scalar_or_container_of<C, V>
  inline void push(C& dst, V&& v)
  {
    if constexpr (requires { dst.push_back(v); }) {
      dst.push_back(v);
    } else if constexpr (std::is_same_v<C, std::remove_cvref_t<V>>) {
      dst = std::move(v);
    } else {
      static_assert(false_type<C, V>, "bad destination");
    }
  }
};

/******************************************************************************/
/**
 * A decoder that produces a vector of scalar values.
 */
template <typename V, typename C = std::vector<V>, typename N = NullSkip<V>>
  requires scalar_or_container_of<C, V>
class Scalar final : Helper<Scalar<V, C, N>> {
public:
  /// The types of values this decoder can decode.
  using value_type = V;

  /// The range or scalar type.
  using range_type = C;

  /// The null decoder type.
  using null_decoder_type = N;

  /// Default constructor.
  Scalar()
      : null_decoder_({})
  {
  }

  /// Constructor.
  Scalar(const null_decoder_type& decoder)
      : null_decoder_(decoder)
  {
  }

  /// Decoding function.
  void decode(const std::shared_ptr<arrow::Array>& src, C& dst)
  {
    using array_type = type_traits<enum_type_v<V>>::array_type;
    using array_ptr_type = std::shared_ptr<array_type>;

    array_ptr_type casted = std::static_pointer_cast<array_type>(src);

    if constexpr (requires { null_decoder_.chunk(casted); }) {
      null_decoder_.chunk(casted);
    }

    for (int64_t i : std::views::iota(0, casted->length())) {
      if (casted->IsNull(i)) {
        std::optional<V> value = std::invoke(null_decoder_, i);
        if (value.has_value()) this->push(dst, std::move(*value));
      } else {
        V value = unsafe_array_value<enum_type_v<V>>(casted, i);
        this->push(dst, std::move(value));
      }
    }
  }

private:
  null_decoder_type null_decoder_;
};

/******************************************************************************/
/**
 * A decoder where array elements are lists.
 *
 * NULL `ListArray` elements are skipped.
 *
 * NULL elements inside the `ListArray` elements are decoded using the
 * (optionally) provided null decoder.
 */
template <typename V, typename C = std::vector<V>, typename N = NullSkip<V>>
  requires Decoders::scalar_or_container_of<C, V>
class List final : Helper<List<V, C>> {
public:
  /// Decodes vectors of type T.
  using value_type = std::vector<V>;

  /// The range or scalar type.
  using range_type = C;

  /// The null decoder type.
  using null_decoder_type = N;

  /// Constructor.
  List() {}

  // Constructor where you can pass a null decoder to the scalar decoder.
  List(const null_decoder_type& null_decoder)
      : scalar_decoder_(null_decoder)
  {
  }

  /// Decoding function.
  void decode(const std::shared_ptr<arrow::Array>& src, C& dst)
  {
    visit(src, [&](int64_t, const std::shared_ptr<arrow::Array>& values) {
      value_type res;
      res.reserve(values->length());
      scalar_decoder_.decode(values, res);
      this->push(dst, res);
    });
  }

private:
  Scalar<V, C, N> scalar_decoder_;
};

/******************************************************************************/
/**
 * An array transformer that returns its argument unchanged.
 */
struct IdentityTransform {
  const std::shared_ptr<arrow::Array>&
  operator()(int64_t, const std::shared_ptr<arrow::Array>& a) const
  {
    return a;
  }
};

/******************************************************************************/
/**
 * A decoder that handles array elements that are lists, and the
 * destination object is a list of scalar values.
 *
 * This class can decode null values in the lists, and also transform
 * the lists using a helper object.  Once use of the transformer
 * object is to decode delta encoding prior to null decoding.
 *
 * Transformers are called with two arguments:
 *
 *   - The index of the array element currently being decoded.
 *
 *   - The array element itself, as an Arrow Array.
 *
 * The transformer can return one of the following types:
 *
 *   - `std::shared_ptr<arrow::Array>` which contains the transformed
 *      values that can then be decoded.
 *
 *   - A pair where the first element is a starting value that should
 *     be inserted into the destination and the second element is an
 *     Arrow array to decode.
 *
 *   - A container of decoded values when can be inserted into the
 *     destination vector and further decoding can be skipped.
 */
template <typename Value,
          typename Container = std::vector<Value>,
          typename NullDecoder = NullSkip<Value>,
          typename Transformer = IdentityTransform>
  requires Decoders::scalar_or_container_of<Container, Value>
class Flattened final : Helper<Flattened<Value, Container>> {
public:
  /// The types of values this decoder can decode.
  using value_type = Value;

  /// The type of return value allowed from transformers.
  using transform_result_type =
      std::variant<std::shared_ptr<arrow::Array>,
                   std::pair<Value, std::shared_ptr<arrow::Array>>,
                   std::shared_ptr<Container>>;

  // Constructor where you can pass a null decoder to the scalar decoder.
  Flattened(const NullDecoder& null_decoder, Transformer transformer = {})
      : scalar_decoder_(null_decoder)
      , transformer_(transformer)
  {
  }

  /// Decoding function.
  void decode(const std::shared_ptr<arrow::Array>& src, Container& dst)
  {
    index_ += visit(src, [&](int64_t i, const std::shared_ptr<arrow::Array>& elm) {
      transform_result_type values(transformer_(index_ + i, elm));

      std::visit(
          [&](auto&& v) -> void {
            using U = std::decay_t<decltype(v)>;
            using P = std::pair<Value, std::shared_ptr<arrow::Array>>;

            if constexpr (std::is_same_v<U, std::shared_ptr<arrow::Array>>) {
              scalar_decoder_.decode(v, dst);
            } else if constexpr (std::is_same_v<U, P>) {
              dst.push_back(v.first);
              scalar_decoder_.decode(v.second, dst);
            } else if constexpr (std::is_same_v<U, std::shared_ptr<Container>>) {
              dst.insert(dst.end(), v->begin(), v->end());
            } else {
              static_assert(false_type<U>, "invalid transform result");
            }
          },
          values);
    });
  }

private:
  Scalar<Value, Container, NullDecoder> scalar_decoder_;
  Transformer transformer_;
  int64_t index_ = 0;
};

} // namespace MzPeak::Util::Decoders
