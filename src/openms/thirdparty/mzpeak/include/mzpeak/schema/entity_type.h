/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <optional>
#include <string>
#include <variant>

namespace MzPeak::Schema {

/**
 * The type of data entity stored in the file.
 */
class EntityType {
public:
  enum Type {
    /// Mass spectra.
    Spectrum,

    /// Chromatograms or other measurements over time.
    Chromatogram,

    /// Similar to Spectrum except the unit of measure is a wavelength
    /// measurement.
    WavelengthSpectrum,
  };

  // Internal storage type.
  using value_type = std::variant<Type, std::string>;

  /// Constructor.
  EntityType(std::string_view);

  /// Constructor.
  EntityType(Type);

  /**
   * Return the string representation of an entity type.
   */
  std::string to_string() const;

  /**
   * Return the enumerated type if it is known.
   */
  std::optional<Type> type() const;

  /**
   * The name of the index column matching this entity.
   */
  std::string index_column_name() const;

  /**
   * The name of the array index metadata key.
   */
  std::string array_index_name() const;

  /**
   * The metadata key for the entity count.
   */
  std::string metadata_count_key() const;

private:
  value_type val_;
};

} // namespace MzPeak::Schema
