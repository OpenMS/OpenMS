/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <arrow/builder.h>
#include <arrow/type_fwd.h>
#include <boost/type_index.hpp>
#include <cstdint>
#include <optional>
#include <parquet/types.h>
#include <type_traits>
#include <variant>

#include "mzpeak/exception.h"
#include "mzpeak/util/compat.h" // IWYU pragma: keep

namespace parquet::schema {
class PrimitiveNode;
}

namespace MzPeak::Util {

/**
 * List of supported data types.
 */
template <typename T>
concept supported_type = std::same_as<std::remove_cvref_t<T>, int8_t> ||
                         std::same_as<std::remove_cvref_t<T>, uint8_t> ||
                         std::same_as<std::remove_cvref_t<T>, int32_t> ||
                         std::same_as<std::remove_cvref_t<T>, uint32_t> ||
                         std::same_as<std::remove_cvref_t<T>, int64_t> ||
                         std::same_as<std::remove_cvref_t<T>, uint64_t> ||
                         std::same_as<std::remove_cvref_t<T>, float> ||
                         std::same_as<std::remove_cvref_t<T>, double> ||
                         std::same_as<std::remove_cvref_t<T>, std::string>;

/// A variant that can hold any supported type.
using any_value_type = std::variant<int8_t,
                                    uint8_t,
                                    int32_t,
                                    uint32_t,
                                    int64_t,
                                    uint64_t,
                                    float,
                                    double,
                                    std::string>;

/// A variant that can hold pairs of any supported type.
using any_pair_type = std::variant<std::pair<int8_t, int8_t>,
                                   std::pair<uint8_t, uint8_t>,
                                   std::pair<int32_t, int32_t>,
                                   std::pair<uint32_t, uint32_t>,
                                   std::pair<int64_t, int64_t>,
                                   std::pair<uint64_t, uint64_t>,
                                   std::pair<float, float>,
                                   std::pair<double, double>,
                                   std::pair<std::string, std::string>>;

/// Enum of supported data types for tracking at run time.
enum class Type {
  Int8,
  UInt8,
  Int32,
  UInt32,
  Int64,
  UInt64,
  Float32,
  Float64,
  ByteArray,
};

/**
 * Compile-time type traits.
 *
 * This is a compile-time map from the Type enum to the following
 * traits:
 *
 *   - `name`:
 *
 *     The name of the type for error messages.
 *
 *   - `value_type`:
 *
 *     The associated C++ type.
 *
 *   - `parquet_type`:
 *
 *     The Parquet storage type.  NOTE: This type can be larger than
 *     the C++ type and does not have unsigned versions.  Therefore
 *     casting to the C++ type is required.
 *
 *   - `array_type`:
 *
 *     The Arrow array type that holds this Type.
 */
template <Type T> struct type_traits;

template <> struct type_traits<Type::Int8> {
  static constexpr const char* name = "int8";
  using value_type = int8_t;
  using parquet_type = parquet::Int32Type;
  using array_type = arrow::Int8Array;
  using builder_type = arrow::Int8Builder;
};

template <> struct type_traits<Type::UInt8> {
  static constexpr const char* name = "uint8";
  using value_type = uint8_t;
  using parquet_type = parquet::Int32Type;
  using array_type = arrow::UInt8Array;
  using builder_type = arrow::UInt8Builder;
};

template <> struct type_traits<Type::Int32> {
  static constexpr const char* name = "int32";
  using value_type = int32_t;
  using parquet_type = parquet::Int32Type;
  using array_type = arrow::Int32Array;
  using builder_type = arrow::Int32Builder;
};

template <> struct type_traits<Type::UInt32> {
  static constexpr const char* name = "uint32";
  using value_type = uint32_t;
  using parquet_type = parquet::Int32Type;
  using array_type = arrow::UInt32Array;
  using builder_type = arrow::UInt32Builder;
};

template <> struct type_traits<Type::Int64> {
  static constexpr const char* name = "int64";
  using value_type = int64_t;
  using parquet_type = parquet::Int64Type;
  using array_type = arrow::Int64Array;
  using builder_type = arrow::Int64Builder;
};

template <> struct type_traits<Type::UInt64> {
  static constexpr const char* name = "uint64";
  using value_type = uint64_t;
  using parquet_type = parquet::Int64Type;
  using array_type = arrow::UInt64Array;
  using builder_type = arrow::UInt64Builder;
};

template <> struct type_traits<Type::Float32> {
  static constexpr const char* name = "float";
  using value_type = float;
  using parquet_type = parquet::FloatType;
  using array_type = arrow::FloatArray;
  using builder_type = arrow::FloatBuilder;
};

template <> struct type_traits<Type::Float64> {
  static constexpr const char* name = "double";
  using value_type = double;
  using parquet_type = parquet::DoubleType;
  using array_type = arrow::DoubleArray;
  using builder_type = arrow::DoubleBuilder;
};

template <> struct type_traits<Type::ByteArray> {
  static constexpr const char* name = "bytes";
  using value_type = std::string;
  using parquet_type = parquet::ByteArrayType;
  using array_type = arrow::StringArray;
  using builder_type = arrow::StringBuilder;
};

/******************************************************************************/
/**
 * Compile-type map from C++ types back to the Type enum.
 */
template <supported_type T> struct type_from_value_type;

template <> struct type_from_value_type<int8_t> {
  static constexpr Type enum_type = Type::Int8;
};

template <> struct type_from_value_type<uint8_t> {
  static constexpr Type enum_type = Type::UInt8;
};

template <> struct type_from_value_type<int32_t> {
  static constexpr Type enum_type = Type::Int32;
};

template <> struct type_from_value_type<uint32_t> {
  static constexpr Type enum_type = Type::UInt32;
};

template <> struct type_from_value_type<int64_t> {
  static constexpr Type enum_type = Type::Int64;
};

template <> struct type_from_value_type<uint64_t> {
  static constexpr Type enum_type = Type::UInt64;
};

template <> struct type_from_value_type<float> {
  static constexpr Type enum_type = Type::Float32;
};

template <> struct type_from_value_type<double> {
  static constexpr Type enum_type = Type::Float64;
};

template <> struct type_from_value_type<std::string> {
  static constexpr Type enum_type = Type::ByteArray;
};

/// Helper to use the above map.
template <supported_type T>
inline constexpr Type enum_type_v = type_from_value_type<T>::enum_type;

/**
 * Cast or convert a value from one type (usually Type::parquet_type)
 * to another (usually Type::value_type).
 *
 * This is needed because Parquet/Arrow uses std::string_view for byte
 * arrays, but that means that the original arrow array needs to
 * remain resident in memory.  Therefore we need to copy the memory
 * referenced by a std::string_view into a std::string.
 */
template <typename To, typename From> To safe_cast_or_copy(From value)
{
  using from_t = std::remove_cvref_t<From>;
  using to_t = std::remove_cvref_t<To>;

  if constexpr (std::is_convertible_v<From, To>) {
    return value;
  } else if constexpr (std::is_same_v<to_t, std::string> &&
                       std::is_same_v<from_t, parquet::ByteArray>) {
    return parquet::ByteArrayToString(value);
  } else {
    static_assert(false_type<From, To>, "no conversion available");
  }
}

/**
 * Return a type for the given parquet node.
 */
std::optional<Type> type_from_parquet(const parquet::schema::PrimitiveNode&);

/**
 * Lift a run-time type to a template parameter for a template
 * function.  Very useful when combined with generic/template lambdas.
 */
template <typename Fn, typename... Args>
decltype(auto) lift_type(Type t, Fn&& func, Args&&... args)
{
  using enum Type;

  switch (t) {
  case Int8:
    return std::forward<Fn>(func).template operator()<Int8>(
        std::forward<Args>(args)...);
  case UInt8:
    return std::forward<Fn>(func).template operator()<UInt8>(
        std::forward<Args>(args)...);
  case Int32:
    return std::forward<Fn>(func).template operator()<Int32>(
        std::forward<Args>(args)...);
  case UInt32:
    return std::forward<Fn>(func).template operator()<UInt32>(
        std::forward<Args>(args)...);
  case Int64:
    return std::forward<Fn>(func).template operator()<Int64>(
        std::forward<Args>(args)...);
  case UInt64:
    return std::forward<Fn>(func).template operator()<UInt64>(
        std::forward<Args>(args)...);
  case Float32:
    return std::forward<Fn>(func).template operator()<Float32>(
        std::forward<Args>(args)...);
  case Float64:
    return std::forward<Fn>(func).template operator()<Float64>(
        std::forward<Args>(args)...);
  case ByteArray:
    return std::forward<Fn>(func).template operator()<ByteArray>(
        std::forward<Args>(args)...);
  }

  throw TypeError("unknown type (type index " + std::to_string(static_cast<int>(t)) +
                  ")");
}

/**
 * Return `true` if the template type is the same as the given types
 * `value_type` from it's type traits.
 */
template <typename T> bool is_same_type(Type type)
{
  using TT = std::decay_t<T>;

  return lift_type(type, []<Type U>() -> bool {
    return std::is_same_v<TT, typename type_traits<U>::value_type>;
  });
}

/**
 * Throw an exception if `T` and `type_traits<Type>::value_type` are
 * not the same type.
 */
template <typename T> void ensure_same_type(Type type)
{
  if (!is_same_type<T>(type)) {
    lift_type(type, []<Type U>() -> void {
      std::string msg("run time type mismatch: ");
      msg += type_traits<U>::name;
      msg += " != ";
      msg += boost::typeindex::type_id<T>().pretty_name();
      throw TypeError(msg);
    });
  }
}

} // namespace MzPeak::Util
