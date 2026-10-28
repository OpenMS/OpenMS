/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include <parquet/schema.h>
#include <parquet/types.h>

#include "mzpeak/util/types.h"

namespace MzPeak::Util {

/******************************************************************************/
std::optional<Type> type_from_parquet(const parquet::schema::PrimitiveNode& node)
{
  using enum Type;

  // Logical types are the most recent addition.
  const std::shared_ptr<const parquet::LogicalType>& logical = node.logical_type();

  if (logical->is_int()) {
    auto logical_int =
        std::static_pointer_cast<const parquet::IntLogicalType>(logical);

    switch (logical_int->bit_width()) {
    case 8:
      return (logical_int->is_signed() ? Int8 : UInt8);
    case 32:
      return (logical_int->is_signed() ? Int32 : UInt32);
    case 64:
      return (logical_int->is_signed() ? Int64 : UInt64);
    default:
      break;
    }
  } else if (logical->is_string()) {
    return ByteArray;
  }

  // Try to use the converted type which is more accurate.
  switch (node.converted_type()) {
  case parquet::ConvertedType::UINT_8:
    return UInt8;
  case parquet::ConvertedType::UINT_32:
    return UInt32;
  case parquet::ConvertedType::UINT_64:
    return UInt64;
  case parquet::ConvertedType::INT_8:
    return Int8;
  case parquet::ConvertedType::INT_32:
    return Int32;
  case parquet::ConvertedType::INT_64:
    return Int64;

  // What we don't handle:
  case parquet::ConvertedType::NONE:
  case parquet::ConvertedType::UTF8:
  case parquet::ConvertedType::MAP:
  case parquet::ConvertedType::MAP_KEY_VALUE:
  case parquet::ConvertedType::LIST:
  case parquet::ConvertedType::ENUM:
  case parquet::ConvertedType::DECIMAL:
  case parquet::ConvertedType::DATE:
  case parquet::ConvertedType::TIME_MILLIS:
  case parquet::ConvertedType::TIME_MICROS:
  case parquet::ConvertedType::TIMESTAMP_MILLIS:
  case parquet::ConvertedType::TIMESTAMP_MICROS:
  case parquet::ConvertedType::INT_16:
  case parquet::ConvertedType::UINT_16:
  case parquet::ConvertedType::JSON:
  case parquet::ConvertedType::BSON:
  case parquet::ConvertedType::INTERVAL:
  case parquet::ConvertedType::NA:
  case parquet::ConvertedType::UNDEFINED:
    break;
  }

  // Fall back to the physical type.
  switch (node.physical_type()) {
  case parquet::Type::INT32:
    return Int32;
  case parquet::Type::INT64:
    return Int64;
  case parquet::Type::FLOAT:
    return Float32;
  case parquet::Type::DOUBLE:
    return Float64;
  case parquet::Type::BYTE_ARRAY:
    return ByteArray;

  // What we don't handle:
  case parquet::Type::FIXED_LEN_BYTE_ARRAY:
  case parquet::Type::INT96:
  case parquet::Type::BOOLEAN:
  case parquet::Type::UNDEFINED:
    break;
  }

  // Fail.
  return {};
}

} // namespace MzPeak::Util
