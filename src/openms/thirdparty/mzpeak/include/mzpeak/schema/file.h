/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <boost/json.hpp>
#include <optional>
#include <string>
#include <vector>

#include "mzpeak/schema/data_kind.h"
#include "mzpeak/schema/entity_type.h"

namespace MzPeak::Schema {
namespace json = boost::json;

/**
 * A description of a file in the mzPeak archive.
 */
class File final {
public:
  /**
   * One entry of a file's `column_mapping`: the binding between a plain Parquet
   * column path and the CV term it carries.
   *
   * Newer writers name metadata columns plainly (`scan_start_time`) and declare
   * the CV term here, instead of encoding it in the column name
   * (`MS_1000016_scan_start_time_unit_UO_0000031`).  Resolving fields through
   * this mapping is therefore the stable way to find them, and it is the ONLY
   * place the unit is stated -- which matters for retention time, where
   * mistaking minutes for seconds is a 60x error.
   */
  struct Column {
    /// Human-readable term name, e.g. "scan start time".
    std::string name;

    /// Parquet column path within the file, e.g. "scan_start_time".
    std::string path;

    /// CV accession of the term, e.g. "MS:1000016".  Absent if not declared.
    std::optional<std::string> accession;

    /// CV accession of the unit, e.g. "UO:0000031" (minutes).  Absent when the
    /// term is unitless or the writer did not declare one.
    std::optional<std::string> unit;

    /// This column marks the presence of a term rather than carrying a value.
    ///
    /// Two forms, per the specification: a BOOLEAN column where `true` means
    /// this mapping's own accession applies to the row, and a STRING column
    /// whose value is the CURIE of a CHILD of this mapping's accession (so
    /// `MS:1000525` "spectrum representation" yields `MS:1000127` centroid or
    /// `MS:1000128` profile).  A standardised column may leave the flag unset,
    /// so its absence does not mean the column carries a value.
    bool term_marker = false;

    /// Equality operator (File's operator== needs it).
    bool operator==(const Column&) const = default;
  };

  /// Constructor from a file name.
  explicit File(const std::string& name);

  /// Conversion from JSON.
  explicit File(const json::object&);

  /// Equality operator.
  bool operator==(const File&) const;

  /// The name of this file.
  const std::string& file_name() const { return file_name_; }

  /// This file's data kind.
  DataKind data_kind() const { return data_kind_; }

  /// The SHA-512 the index records for this file, lowercase hex with no
  /// separators, or absent when the archive declares none.  The specification
  /// requires one per indexed file, but archives written before that rule are
  /// common and are not rejected over it.
  const std::optional<std::string>& checksum() const { return checksum_; }

  /// This file's entity type.
  EntityType entity_type() const { return entity_type_; }

  /// Declared column/CV bindings for this file.  Empty for signal-data files
  /// (which describe their arrays through the array index instead) and for
  /// older writers that encoded the CV term in the column name.
  const std::vector<Column>& columns() const { return columns_; }

  /// Return the column path bound to @p accession, or nullopt when this file
  /// declares no such term.  Prefer this over hardcoding a column name.
  std::optional<std::string> path_for(const std::string_view& accession) const;

  /// Return the unit accession declared for @p accession's column, or nullopt
  /// when the term is absent or carries no declared unit.
  std::optional<std::string> unit_for(const std::string_view& accession) const;

private:
  std::string file_name_;
  DataKind data_kind_ = DataKind("other");
  EntityType entity_type_ = EntityType("other");
  std::vector<Column> columns_;

  // SHA-512 hex digest from the index, if the archive declares one.
  std::optional<std::string> checksum_;
};

} // namespace MzPeak::Schema
