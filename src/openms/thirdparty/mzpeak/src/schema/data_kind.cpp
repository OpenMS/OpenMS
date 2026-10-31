/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include <utility>

#include "mzpeak/schema/data_kind.h"

namespace MzPeak::Schema {

/******************************************************************************/
DataKind::value_type data_kind_from_string(std::string_view s)
{
  // The two string variants of `DataArray` are for backwards
  // compatibility. See https://github.com/HUPO-PSI/mzPeak/issues/26

  using enum DataKind::Type;

  if (s == "data_arrays") {
    return DataArray;
  } else if (s == "data arrays") {
    return DataArray;
  } else if (s == "peaks") {
    return Peaks;
  } else if (s == "metadata") {
    return Metadata;
  } else if (s == "scans") {
    return Scans;
  } else if (s == "precursors") {
    return Precursors;
  } else if (s == "selected_ions") {
    return SelectedIons;
  } else if (s == "products") {
    return Products;
  } else if (s == "proprietary") {
    return Proprietary;
  } else {
    return std::string(s);
  }
}

/******************************************************************************/
DataKind::DataKind(std::string_view s)
    : val_(data_kind_from_string(s))
{
}

/******************************************************************************/
DataKind::DataKind(Type t)
    : val_(t)
{
}

/******************************************************************************/
std::string DataKind::to_string() const
{
  return std::visit(
      [](auto&& v) -> std::string {
        using T = std::decay_t<decltype(v)>;

        if constexpr (std::is_same_v<T, std::string>) {
          return v;
        } else {
          switch (v) {
          case DataArray:
            return "data_arrays";
          case Peaks:
            return "peaks";
          case Metadata:
            return "metadata";
          case Scans:
            return "scans";
          case Precursors:
            return "precursors";
          case SelectedIons:
            return "selected_ions";
          case Products:
            return "products";
          case Proprietary:
            return "proprietary";
          }
        }

        std::unreachable();
      },
      val_);
}

/******************************************************************************/
std::optional<DataKind::Type> DataKind::type() const
{
  if (std::holds_alternative<Type>(val_)) {
    return std::get<Type>(val_);
  } else {
    return {};
  }
}

/******************************************************************************/
bool DataKind::is_metadata() const
{
  auto enum_type = type();

  if (enum_type.has_value()) {
    switch (enum_type.value()) {
    case DataArray:
      return false;
    case Peaks:
      return false;
    case Metadata:
      return true;
    case Scans:
      return true;
    case Precursors:
      return true;
    case SelectedIons:
      return true;
    case Products:
      return true;
    case Proprietary:
      return false;
    }

    std::unreachable();
  } else {
    return false;
  }
}

} // namespace MzPeak::Schema
