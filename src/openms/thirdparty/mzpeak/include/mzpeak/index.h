/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "mzpeak/chromatograms.h"
#include "mzpeak/ims_calibration.h"
#include "mzpeak/run_metadata.h"
#include "mzpeak/schema/file.h"
#include "mzpeak/spectra.h"
#include "mzpeak/wavelength_spectra.h"

namespace MzPeak {

namespace IO {
class Archive;
}

namespace Schema {
class File;
}

namespace Util {
class Manager;
}

class Spectra;

/**
 * Read-only access to the index inside a MzPeak archive.
 */
/// One member whose bytes do not hash to the digest the index records.
struct ChecksumMismatch {
  std::string file_name;
  std::string expected; ///< what `mzpeak_index.json` claims
  std::string actual;   ///< what the stored bytes actually hash to
};

/**
 * The outcome of verifying an archive's recorded SHA-512 digests.
 *
 * `unchecked` is reported separately and deliberately: an archive written
 * before the digests became mandatory records none at all, and such an archive
 * would otherwise "pass" verification having proved nothing.  A caller that
 * cares about integrity must look at `verified`, not only at `ok()`.
 */
struct ChecksumReport {
  std::size_t verified = 0;  ///< members whose digest matched their bytes
  std::size_t unchecked = 0; ///< members the index records no digest for
  std::vector<ChecksumMismatch> mismatches;

  /// True when nothing disagreed.  Vacuously true for an archive that records
  /// no digests -- see `verified`.
  bool ok() const { return mismatches.empty(); }
};

class Index {
public:
  /// Constructor.
  Index(std::unique_ptr<MzPeak::IO::Archive>);

  /**
   * Recompute the SHA-512 of every member the index records one for.
   *
   * Reads every byte of every indexed member, so the cost is the size of the
   * archive; that is why it is a separate call rather than part of opening.
   * Members are streamed, not held.
   *
   * Reports rather than throws, so a caller can distinguish a damaged member
   * from an archive that simply predates the requirement, and can decide for
   * itself whether to continue.  `MzPeak::open(path, Validate::Checksums)`
   * wraps this and throws instead.
   */
  ChecksumReport verify_checksums() const;

  /**
   * Return a list of files found in the index.
   */
  const std::vector<Schema::File>& files() const;

  /**
   * Find a file in the mzPeak archive with the given name.
   */
  std::vector<Schema::File>::const_iterator find(std::string_view) const;

  /**
   * Find a file by its `EntityType` and `DataKind` (upstream's API).
   */
  std::vector<Schema::File>::const_iterator find_file(Schema::EntityType::Type,
                                                      Schema::DataKind::Type) const;

  /**
   * The mzPeak format version from the index `metadata.version`, or an
   * empty string if the index does not declare one.
   */
  const std::string& version() const;

  /**
   * Indicates which source to fetch spectra data from (upstream's API).
   */
  enum class SpectraSource {
    /// Use `spectra_data.parquet`, which may hold profile or centroid data.
    Data,

    /// Use the optional `spectra_peaks.parquet`.
    Peaks,
  };

  /**
   * Returns `true` if the given source file exists.
   */
  bool has_spectra(SpectraSource) const;

  /**
   * Access the spectra in the file.
   *
   * @param detail  how much per-spectrum metadata to materialise.  The map is
   *   built once per archive and per detail level and SHARED by every `Spectra`
   *   this `Index` hands out, so calling this once per thread over one shared
   *   `Index` costs the metadata once, not once per thread.  See
   *   @ref MzPeak::MetadataDetail for what `Lean` leaves out.
   * @param source  which table is primary.  With the default `Data` this reader
   *   also attaches the peaks file when present and picks the right table PER
   *   SPECTRUM from its declared representation, so callers rarely need
   *   `Peaks` explicitly.  Throws if `Peaks` is requested and no peaks file
   *   exists -- check `has_spectra` first.
   */
  Spectra spectra(MetadataDetail detail = MetadataDetail::Full,
                  SpectraSource source = SpectraSource::Data) const;

  /**
   * TOF -> m/z calibration declared by this archive, for the Bruker TDF
   * "ims-compact" layout.  `valid` is false when the archive declares none.
   */
  const ImsCalibration& ims_calibration() const;

  /**
   * Access the chromatograms in the file.
   */
  Chromatograms chromatograms() const;

  /**
   * Access the wavelength spectra in the file.
   */
  WavelengthSpectra wavelength_spectra() const;

  /**
   * Run-level metadata from the index `metadata{}` block.  Empty when the
   * archive declares none.
   */
  const RunMetadata& metadata() const;

  /**
   * Access the low-level MzPeak Manager object.
   */
  std::shared_ptr<Util::Manager> manager() const;

protected:
  std::shared_ptr<Util::Manager> manager_;
};

} // namespace MzPeak
