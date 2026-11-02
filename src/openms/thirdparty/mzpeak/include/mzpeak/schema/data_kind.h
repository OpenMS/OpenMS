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
 * Indicates how data is encoded in a parquet file.
 */
class DataKind final {
public:
  enum Type {
    /// Files that contain signal data.
    DataArray,

    /// Processed data (e.g. centroided).  It is expected that the
    /// unprocessed version of the data is present in the same MzPeak
    /// file as a DataArray.
    Peaks,

    /// Metadata relating to one of the other files.
    Metadata,

    /// A scan or acquisition from the original raw file used to
    /// create a spectrum.
    Scans,

    /// The method of precursor-ion selection and activation.
    Precursors,

    /// An ion isolated for dissociation.
    SelectedIons,

    /// When describing single reaction monitoring (SRM) or multiple
    /// reaction monitoring (MRM) experiments, each product ion is
    /// isolated separately with a different isolation window. This
    /// table is usually empty or absent
    Products,

    /// Non-standard file that can't be decoded by this library.
    /// However, users of this library can access the raw bytes of
    /// this file.
    Proprietary,
  };

  // Internal storage type.
  using value_type = std::variant<Type, std::string>;

  /// Constructor.
  DataKind(std::string_view);

  /// Constructor.
  DataKind(Type);

  /**
   * Return the string representation of the data kind.
   */
  std::string to_string() const;

  /**
   * Return the enumerated type of the data kind if it is known.
   */
  std::optional<Type> type() const;

  /**
   * Return `true` if this data kind is a known holder of metadata.
   * For example, the `Scans` and `Precursors` type are both
   * considered to be metadata.
   */
  bool is_metadata() const;

private:
  value_type val_;
};

} // namespace MzPeak::Schema
