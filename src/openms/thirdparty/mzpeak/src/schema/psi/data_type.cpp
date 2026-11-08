/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include <type_traits>
#include <utility>

#include "mzpeak/schema/psi/data_type.h"
#include "mzpeak/util/types.h"

namespace MzPeak::Schema::PSI {

/******************************************************************************/
DataType::DataType(const CV& cv)
    : val_(cv)
{
  if (cv.code() == "MS") {
    const auto& accession = cv.accession();

    if (accession == "1000519") {
      val_ = Int32;
    } else if (accession == "1000521") {
      val_ = Float32;
    } else if (accession == "1000522") {
      val_ = Int64;
    } else if (accession == "1000523") {
      val_ = Float64;
    } else if (accession == "1001479") {
      val_ = ASCII;
    }
  }
}

/******************************************************************************/
DataType::DataType(Type t)
    : val_(t)
{
}

/******************************************************************************/
CV DataType::to_cv() const
{
  return std::visit(
      [](auto&& v) -> CV {
        using T = std::decay_t<decltype(v)>;
        if constexpr (std::is_same_v<T, CV>) {
          return v;
        } else if constexpr (std::is_same_v<T, Type>) {
          switch (v) {
          case Int32:
            return CV("MS", "1000519");
          case Int64:
            return CV("MS", "1000522");
          case Float32:
            return CV("MS", "1000521");
          case Float64:
            return CV("MS", "1000523");
          case ASCII:
            return CV("MS", "1001479");
          }

          std::unreachable();
        }
      },
      val_);
}

/******************************************************************************/
std::optional<DataType::Type> DataType::data_type() const
{
  if (std::holds_alternative<Type>(val_)) {
    return std::get<Type>(val_);
  } else {
    return {};
  }
}

/******************************************************************************/
std::optional<Util::Type> DataType::as_type() const
{
  if (as_type_.has_value()) {
    return as_type_;
  } else if (std::holds_alternative<Type>(val_)) {
    switch (std::get<Type>(val_)) {
    case Int32:
      return Util::Type::Int32;
    case Int64:
      return Util::Type::Int64;
    case Float32:
      return Util::Type::Float32;
    case Float64:
      return Util::Type::Float64;
    case ASCII:
      return Util::Type::ByteArray;
    }
  }

  return {};
}

/******************************************************************************/
void DataType::as_type(Util::Type t) { as_type_ = t; }

/******************************************************************************/
bool DataType::operator==(const DataType& other) const
{
  return val_ == other.val_ && as_type_ == other.as_type_;
}

/******************************************************************************/
bool DataType::operator<(const DataType& other) const
{
  return as_type() < other.as_type();
}

} // namespace MzPeak::Schema::PSI
