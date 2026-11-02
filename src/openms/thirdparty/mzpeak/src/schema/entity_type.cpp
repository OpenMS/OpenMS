/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/schema/entity_type.h"

#include <utility>

namespace MzPeak::Schema {

/******************************************************************************/
std::string entity_type_to_string(EntityType::Type t)
{
  using enum EntityType::Type;

  switch (t) {
  case Spectrum:
    return "spectrum";
  case Chromatogram:
    return "chromatogram";
  case WavelengthSpectrum:
    // UNDERSCORE, not "wavelength spectrum": index_column_name(),
    // array_index_name() and metadata_count_key() are all built by appending
    // to this string, and the space form yields "wavelength spectrum_index",
    // which matches no column in any archive.  Upstream converged on the same
    // spelling in d70ce3a.  The space form is still ACCEPTED when parsing.
    return "wavelength_spectrum";
  }

  std::unreachable();
}

/******************************************************************************/
EntityType::value_type entity_type_from_string(std::string_view s)
{
  using enum EntityType::Type;

  // The two spellings of `WavelengthSpectrum` are for backwards compatibility:
  // upstream renamed the space form to the underscore form, mirroring the
  // `data arrays` -> `data_arrays` rename handled in data_kind.cpp.  Accepting
  // only one silently degrades every entry to `Other`, which makes the whole
  // wavelength collection come back empty (has_uv.mzpeak: 520 rows).

  if (s == "spectrum") {
    return Spectrum;
  } else if (s == "chromatogram") {
    return Chromatogram;
  } else if (s == "wavelength_spectrum") {
    return WavelengthSpectrum;
  } else if (s == "wavelength spectrum") {
    // Kept deliberately after upstream dropped it in d70ce3a: archives written
    // before the rename still carry the space form, and refusing it degrades
    // the entry to `Other` and empties the collection without an error.
    return WavelengthSpectrum;
  } else {
    return std::string(s);
  }
}

/******************************************************************************/
EntityType::EntityType(std::string_view s)
    : val_(entity_type_from_string(s))
{
}

/******************************************************************************/
EntityType::EntityType(Type t)
    : val_(t)
{
}

/******************************************************************************/
std::string EntityType::to_string() const
{
  return std::visit(
      [](auto&& v) {
        using T = std::decay_t<decltype(v)>;

        if constexpr (std::is_same_v<T, Type>) {
          return entity_type_to_string(v);
        } else {
          return v;
        }
      },
      val_);
}

/******************************************************************************/
std::optional<EntityType::Type> EntityType::type() const
{
  return std::visit(
      [](auto&& v) -> std::optional<EntityType::Type> {
        using T = std::decay_t<decltype(v)>;

        if constexpr (std::is_same_v<T, Type>) {
          return v;
        } else {
          return std::nullopt;
        }
      },
      val_);
}

/******************************************************************************/
std::string EntityType::index_column_name() const { return to_string() + "_index"; }

/******************************************************************************/
std::string EntityType::array_index_name() const
{
  return to_string() + "_array_index";
}

/******************************************************************************/
std::string EntityType::metadata_count_key() const { return to_string() + "_count"; }

} // namespace MzPeak::Schema
