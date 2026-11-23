/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

#include <arrow/array.h>

#include "mzpeak/schema/group.h"
#include "mzpeak/util/compat.h" // IWYU pragma: keep
#include "mzpeak/util/decoders.h"

namespace MzPeak::Util {

// Forward declarations.
class Executor;

/******************************************************************************/
using Column = MzPeak::Schema::Column;

/******************************************************************************/
/**
 * Represents a subset of a Parquet file resulting from executing a
 * query.
 */
class Slice final {
public:
  /// Raw, chunked arrays from Parquet.
  using Raw = std::vector<std::shared_ptr<arrow::Array>>;

  /// Destructor.
  ~Slice();

  /**
   * Return a list of fields that can be extracted from this slice.
   */
  const std::vector<Column>& fields() const;

  /**
   * Return true if the given column is in the slice.
   */
  bool has_column(const Column&) const;

  /**
   * Return the raw array for the given field.
   *
   * The returned raw array is removed from the internal storage
   * therefore calling this method again with the same column will
   * fail.
   *
   * NOTE: If you request a field that does not exist in the slice
   * this function will return a nullptr.
   */
  std::shared_ptr<Raw> raw(const Column&);

  /**
   * Decode the first non-null value from the given column.  The
   * column is then removed from internal storage.
   *
   * Template Parameters:
   *
   *   - T: The Decoder class to use (see MzPeak::Util::Decoders)
   *
   *   - R: The destination object to update with the decoded value
   */
  template <typename T, typename R = std::optional<typename T::value_type>>
  void singleton(const Column&, R&, T&& = {});

  /**
   * A helper function that calls the `singleton` method with the
   * scalar decoder.
   */
  template <typename T> void scalar(const Schema::Column&, T&);

  /// Overload that resets the destination type if column is nullopt.
  template <typename T> void scalar(const std::optional<Schema::Column>&, T&);

  /**
   * A helper function that calls the `singleton` method with the list
   * decoder.
   */
  template <typename T> void list(const Schema::Column&, T&);

  /// Overload that resets the destination type if column is nullopt.
  template <typename T> void list(const std::optional<Schema::Column>&, T&);

  /**
   * Extract and decode an array.  The array is then removed from the
   * internal storage.
   *
   * Use one of the decoders defined in `decoders.h`, or write your own.
   *
   * Template Parameters:
   *
   *   - T: The Decoder class to use (see MzPeak::Util::Decoders)
   *
   *   - V: The destination object to fill with decoded values
   */
  template <typename T, typename V = std::vector<typename T::value_type>>
    requires Decoders::from_arrow_array<T, V>
  void array(const Column&, V&, T&& = {});

  /// Override that takes an optional column.
  template <typename T, typename V = std::vector<typename T::value_type>>
    requires Decoders::from_arrow_array<T, V>
  void array(const std::optional<Column>&, V&, T&& = {});

  /// Override that takes a reference to a decoder.
  template <typename T, typename V = std::vector<typename T::value_type>>
    requires Decoders::from_arrow_array<T, V>
  void array(const Column&, V&, T&);

  /**
   * A helper function that calls `array` with a scalar decoder that
   * skips null values.
   */
  template <typename T> void non_null(const Column&, std::vector<T>&);

  /// Override that takes an optional column.
  template <typename T> void non_null(const std::optional<Column>&, std::vector<T>&);

private:
  friend class MzPeak::Util::Executor;

  /// Constructor.
  Slice(const std::vector<Column>&);

  /// Add an array chunk.
  void append(const Column&, std::shared_ptr<arrow::Array>);

  struct Impl;
  std::unique_ptr<Impl> impl_;
};

/******************************************************************************/
template <typename T, typename R>
void Slice::singleton(const Column& field, R& dst, T&& decoder)
{
  std::shared_ptr<Raw> chunks = raw(field);

  if (chunks != nullptr) {
    for (const auto& chunk : *chunks) {
      for (int64_t i : std::views::iota(0, chunk->length())) {
        if (!chunk->IsNull(i)) {
          std::shared_ptr<arrow::Array> sliced = chunk->Slice(i, 1);

          if constexpr (std::is_same_v<R, typename T::value_type>) {
            decoder.decode(sliced, dst);
            return;
          } else {
            typename T::value_type val{};
            decoder.decode(sliced, val);
            dst = {std::move(val)};
            return;
          }
        }
      }
    }
  }

  dst = {};
}

/******************************************************************************/
template <typename T> void Slice::scalar(const Schema::Column& field, T& dst)
{
  if constexpr (requires { typename T::value_type; }) {
    using V = typename T::value_type;
    singleton<Decoders::Scalar<V, V>, T>(field, dst);
  } else {
    singleton<Decoders::Scalar<T, T>, T>(field, dst);
  }
}

/******************************************************************************/
template <typename T>
void Slice::scalar(const std::optional<Schema::Column>& field, T& dst)
{
  if (field.has_value()) {
    scalar(field.value(), dst);
  } else {
    dst = {};
  }
}

/******************************************************************************/
template <typename T> void Slice::list(const Schema::Column& field, T& dst)
{
  using V = T::value_type;
  using D = Decoders::List<V>;
  singleton<D, T>(field, dst);
}

/******************************************************************************/
template <typename T>
void Slice::list(const std::optional<Schema::Column>& field, T& dst)
{
  if (field.has_value()) {
    list(field.value(), dst);
  } else {
    dst = {};
  }
}

/******************************************************************************/
template <typename T, typename V>
  requires Decoders::from_arrow_array<T, V>
void Slice::array(const Column& field, V& v, T&& t)
{
  array(field, v, t);
}

/******************************************************************************/
template <typename T, typename V>
  requires Decoders::from_arrow_array<T, V>
void Slice::array(const std::optional<Column>& field, V& v, T&& t)
{
  if (field.has_value()) {
    array(field.value(), v, t);
  }
}

/******************************************************************************/
template <typename T, typename V>
  requires Decoders::from_arrow_array<T, V>
void Slice::array(const Column& field, V& v, T& t)
{
  std::shared_ptr<Raw> chunks = raw(field);
  if (chunks == nullptr) return;

  if constexpr (requires { v.reserve(std::size_t{}); }) {
    std::size_t size{};

    for (const auto& chunk : *chunks) {
      size += Decoders::guess_array_length(field, chunk);
    }

    v.reserve(v.size() + size);
  }

  for (const auto& chunk : *chunks) {
    t.decode(chunk, v);
  }
}

/******************************************************************************/
template <typename T> void Slice::non_null(const Column& field, std::vector<T>& dst)
{
  using D = Decoders::Scalar<T, std::vector<T>, Decoders::NullSkip<T>>;
  array<D>(field, dst);
}

/******************************************************************************/
template <typename T>
void Slice::non_null(const std::optional<Column>& field, std::vector<T>& dst)
{
  if (field.has_value()) {
    non_null(field.value(), dst);
  }
}

} // namespace MzPeak::Util
