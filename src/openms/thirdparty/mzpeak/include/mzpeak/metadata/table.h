/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <cstdint>
#include <map>
#include <optional>
#include <memory>

#include "mzpeak/schema/data_kind.h"
#include "mzpeak/schema/group.h"
#include "mzpeak/spectrum_metadata.h"
#include "mzpeak/util/index_map.h"
#include "mzpeak/util/parquet.h"
#include "mzpeak/util/projection.h"

namespace MzPeak::Metadata {

/**
 * Low-level access to metadata files in a mzpeak file.
 */
class Table final {
public:
  /// Constructor.
  Table(std::unique_ptr<Util::Parquet>);

  /// Attach a facet file (scans / precursors / selected ions) that belongs to
  /// this entity's metadata.  Newer writers split the facets into their own
  /// files; older ones nest them as struct columns and never call this.
  void add_facet(Schema::DataKind::Type, std::unique_ptr<Util::Parquet>);

  /// Destructor.
  ~Table();

  /**
   * Return a group with the given name.
   */
  std::shared_ptr<Schema::Group> group(std::string_view) const;

  /**
   * Read all rows from the given group where the index column
   * matches the given value.
   *
   * Returns a slice with the given column projection.
   */
  std::unique_ptr<Util::Slice> indexed(uint64_t,
                                       const std::shared_ptr<Schema::Group>&,
                                       const Util::Projection&) const;

  /**
   * Read the per-spectrum descriptive metadata map (RT, precursors, selected
   * ions, scan windows, ion mobility) keyed by spectrum.index.
   *
   * Reads the whole metadata Parquet table once — the caller is expected to
   * cache the result.  It is NOT small: measured on a 7,534-spectrum Thermo
   * run it costs 26.7 MB, about 3.5 KB per spectrum, and it is built before
   * the first peak is read.  Pass `MetadataDetail::Lean` when the CV-parameter
   * lists, scan windows and auxiliary arrays are not wanted.
   */
  Util::IndexMap<SpectrumMetadata>
  read_spectrum_metadata(MetadataDetail = MetadataDetail::Full) const;

  /// The map at MetadataDetail::Minimal, or nullopt when this archive does
  /// not fit it; see Util::read_minimal_spectra_metadata().
  std::optional<Util::IndexMap<MinimalSpectrumMetadata>>
  read_minimal_spectrum_metadata() const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace MzPeak::Metadata
