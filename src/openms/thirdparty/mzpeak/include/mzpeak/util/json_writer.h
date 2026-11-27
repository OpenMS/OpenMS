/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <boost/json.hpp>
#include <string>
#include <vector>

namespace MzPeak::Util {

/**
 * Emit the file-level array index JSON for a point-layout spectra_data
 * table whose columns are point.spectrum_index (uint64), point.mz
 * (double) and point.intensity (float).
 *
 * The returned JSON mirrors the `spectrum_array_index` blob stored in
 * the Parquet key/value metadata: an object with a `prefix` and an
 * `entries` array. The index column itself is *not* serialized; it is
 * synthesized by the reader (see ArrayIndex parsing).
 */
std::string point_spectra_array_index_json();

/**
 * Emit the file-level array index for a point-layout chromatograms_data table
 * whose columns are point.chromatogram_index (uint64), point.time (double,
 * MINUTES) and point.intensity (float).
 *
 * @param intensity_unit CURIE for the intensity column.  A chromatogram's
 *        intensity is not always detector counts -- an absorbance trace from a
 *        diode-array detector is in UO:0000269 -- and the unit is the only
 *        thing distinguishing them once decoded.
 */
std::string point_chromatograms_array_index_json(
    const std::string& intensity_unit = "MS:1000131");

/**
 * Emit the file-level array index for a point-layout wavelength_spectra_data
 * table whose columns are point.wavelength_spectrum_index (uint64),
 * point.wavelength (float, NANOMETRES) and point.intensity (float).
 */
std::string
point_wavelength_array_index_json(const std::string& intensity_unit = "MS:1000131");

/**
 * A single entry in the mzpeak_index.json `files` array.
 */
/// One `column_mapping` entry: binds a plain column path to its CV term.  The
/// reference reader resolves metadata columns THROUGH this mapping — a column
/// that is absent from it is logged as "unspecified" and effectively ignored,
/// so omitting the mapping silently loses the field.
struct IndexColumnMapping {
  std::string name;
  std::string path;
  std::string accession;
  std::string unit; // empty => null
};

struct IndexFileEntry {
  std::string name;
  std::string entity_type;
  std::string data_kind;
  std::vector<IndexColumnMapping> column_mapping;

  /// Lowercase, separator-free SHA-512 of this member's bytes.  The
  /// specification requires one per indexed file; empty is serialised as JSON
  /// null, which is what an archive written before that rule looks like.
  std::string checksum;
};

/**
 * Emit the contents of mzpeak_index.json: a `files` array and a minimal
 * but valid `metadata` object carrying at least `version`.
 */
std::string mzpeak_index_json(const std::vector<IndexFileEntry>& files,
                              const std::string& version = "0.9.0");

/**
 * WRT-2 — emit mzpeak_index.json with run-level metadata blocks.
 *
 * Identical to the @ref mzpeak_index_json overload above, but merges the
 * members of `run_metadata` (already serialized run-level blocks, e.g. from
 * `RunMetadata::to_json()`) into the emitted `metadata{}` object alongside
 * `version`.  `version` always wins over any `version` key in `run_metadata`.
 */
std::string mzpeak_index_json(const std::vector<IndexFileEntry>& files,
                              const std::string& version,
                              const boost::json::object& run_metadata);

} // namespace MzPeak::Util
