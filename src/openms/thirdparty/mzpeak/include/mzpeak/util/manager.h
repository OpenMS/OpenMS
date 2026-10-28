/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "mzpeak/ims_calibration.h"
#include "mzpeak/io/archive.h"
#include "mzpeak/run_metadata.h"
#include "mzpeak/schema/file.h"
#include "mzpeak/spectrum_metadata.h"
#include "mzpeak/util/index_map.h"
#include "mzpeak/util/parquet.h"

namespace MzPeak::Util {

/**
 * The archive itself: the parsed `mzpeak_index.json` plus the means to open any
 * file it names.
 *
 * Held by `shared_ptr` and handed to every reader built from one archive, so
 * that the index is parsed once and the archive outlives the readers that read
 * through it.
 */
class Manager final {
public:
  /// Constructor.
  Manager(std::unique_ptr<MzPeak::IO::Archive>);

  /**
   * Return a vector of files that are located in the mzPeak archive.
   */
  const std::vector<Schema::File>& files() const;

  /**
   * Find a file given its name.
   */
  std::vector<Schema::File>::const_iterator find_file(std::string_view) const;

  /**
   * Find a file given its `EntityType` and `DataKind`.
   */
  std::vector<Schema::File>::const_iterator find_file(Schema::EntityType::Type,
                                                      Schema::DataKind::Type) const;

  /**
   * Open a Parquet file from the mzPeak archive.
   */
  /// A reader over @p file.  Every Parquet from one Manager shares the
  /// row-group cache below, so a group is decoded once for all of them, and
  /// every Parquet over one member shares its parsed footer.  Thread-safe.
  std::unique_ptr<Util::Parquet> parquet(const Schema::File&) const;

  /**
   * Open a member for RAW byte access, without interpreting it as Parquet.
   *
   * Exists for integrity checking, which must hash exactly the bytes the
   * archive stores; decoding and re-encoding would digest a different
   * sequence.
   */
  std::unique_ptr<MzPeak::IO::File> read_member(const std::string&) const;

  /// The archive-wide decoded-row-group cache: size its budget to the number
  /// of readers you run, read its stats to see what they cost.  Expect about
  /// one decoded group per live reader (two with a separate peaks file) on
  /// top of the budget: each Parquet from parquet() holds the last group it
  /// read.  See RowGroupCache.
  RowGroupCache& row_group_cache() const { return *row_group_cache_; }

  /**
   * The mzPeak format version from the index `metadata.version`, or an empty
   * string if the index does not declare one.
   */
  const std::string& version() const { return version_; }

  /**
   * TOF -> m/z calibration declared by this archive, for the Bruker TDF
   * "ims-compact" layout.  `valid` is false when the archive declares none.
   */
  const ImsCalibration& ims_calibration() const { return ims_; }

  /**
   * Run-level metadata from the index `metadata{}` block (run, file
   * description, software, instrument configurations, ...).  Empty when the
   * index carries none.
   */
  const RunMetadata& metadata() const { return metadata_; }

  /// The per-spectrum descriptive metadata map, keyed by spectrum index.
  using SpectrumMetadataMap = IndexMap<SpectrumMetadata>;

  /**
   * The cached descriptive metadata for this archive at the given detail,
   * building it with @p build on first request.
   *
   * WHY IT LIVES HERE.  The map is the single largest fixed cost of opening a
   * run (26.7 MB on a 7,534-spectrum Thermo file) and it is immutable once
   * built, so every reader over one archive should share one copy.  Before
   * this, each `Index::spectra()` built its own -- and since the only safe way
   * to read one archive from N threads is one `Spectra` per thread, that was N
   * copies of the same table.  A shared `Index` now pays for it once.
   *
   * Thread-safe: concurrent callers serialise on the build and the later ones
   * find it cached.  The returned map is const and shared, so readers on
   * different threads may hold it simultaneously.
   *
   * @param build  called at most once per detail level, under the lock.
   */
  std::shared_ptr<const SpectrumMetadataMap>
  spectrum_metadata(MetadataDetail detail,
                    const std::function<SpectrumMetadataMap()>& build) const;

  /// The same at MetadataDetail::Minimal.  Null when the archive does not
  /// fit it (@p build returned nullopt), which is cached too.
  using MinimalSpectrumMetadataMap = IndexMap<MinimalSpectrumMetadata>;
  std::shared_ptr<const MinimalSpectrumMetadataMap> minimal_spectrum_metadata(
      const std::function<std::optional<MinimalSpectrumMetadataMap>()>& build) const;

private:
  std::shared_ptr<MzPeak::IO::Archive> archive_;
  std::shared_ptr<RowGroupCache> row_group_cache_ = std::make_shared<RowGroupCache>();
  std::vector<Schema::File> files_;

  // The mzPeak format version from metadata.version (empty if absent).
  std::string version_;

  // TOF -> m/z calibration for the ims-compact layout (valid=false if absent).
  ImsCalibration ims_;

  // Typed run-level metadata blocks from metadata{} (empty if absent).
  RunMetadata metadata_;

  // Built on demand, at most once per detail level, and shared from there on.
  // Mutable because caching is not an observable state change: every accessor
  // that reaches it is logically const.
  mutable std::map<MetadataDetail, std::shared_ptr<const SpectrumMetadataMap>>
      spectrum_metadata_;
  mutable std::optional<std::shared_ptr<const MinimalSpectrumMetadataMap>>
      minimal_spectrum_metadata_;
  mutable std::mutex spectrum_metadata_mutex_;

  // Each member's parsed Parquet footer, by file name, from the first Parquet
  // opened over it; see parquet().  Mutable for the same reason as above.
  mutable std::map<std::string, Parquet::file_metadata_t, std::less<>> footers_;
  mutable std::mutex footers_mutex_;
};

} // namespace MzPeak::Util
