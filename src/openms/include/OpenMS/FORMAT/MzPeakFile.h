// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <memory>

namespace OpenMS
{

/**
  @brief File adapter for the mzPeak format (Apache Parquet payloads packaged in a ZIP container).

  The mzPeak format stores mass spectrometric data as one or more Apache Parquet tables bundled
  inside a ZIP archive.  Three usage modes are provided:

  - **Bulk load/store**: load() / store() read or write a complete MSExperiment.
  - **Streaming**: transform() pushes individual spectra through an IMSDataConsumer without
    building a full MSExperiment in memory.
  - **On-disc random access**: openFile() + getSpectrum() give indexed access to individual
    spectra without loading the whole file.  Only the (small) per-spectrum metadata table is
    read into memory; the large data/peaks Parquet tables stay on disc and are accessed via
    Parquet row-group statistics seeks.

  @note The three modes are independent: load() and getSpectrum() can both be called on the
  same instance without interfering with each other.

  @note **Thread safety**: the on-disc interface (openFile / getSpectrum) is @em not thread-safe.
  Concurrent calls to getSpectrum() race on the shared file readers and the lazy row-group index.
  Give each thread its own MzPeakFile instance.

  @note **Copy semantics**: MzPeakFile is move-only (not copyable) because on-disc state holds
  unique file handles.  Use a separate instance per thread.

  @note **Known limitations**:
  - The on-disc interface assumes the @c spectrum_index values in the Parquet point tables are
    contiguous starting at 0 (exactly as written by store()).  Non-contiguous third-party files
    may cause getSpectrum() to throw InvalidParameter or return a different spectrum.
  - PeakFileOptions (MS-level / RT / m/z filters) are honoured by load() and transform() but
    @em not by getSpectrum().  getSpectrum() always returns the full spectrum regardless of
    options.
  - The point-column prefix is assumed to be @c "point" (the mzPeak 0.x default); files using
    a non-standard prefix will lose peak data silently (same limitation applies to load() /
    transform()).
*/
class OPENMS_DLLAPI MzPeakFile
{
public:
  typedef MSExperiment MapType;

  /** @name Constructors and Destructor
   */
  //@{
  /// Default constructor
  MzPeakFile();

  /// Default destructor
  ~MzPeakFile();
  //@}

  /** @name PeakFileOptions accessors
   *
   * Options apply to load() and transform(). getSpectrum() ignores them.
   */
  //@{
  /// Returns a mutable reference to the current file options.
  PeakFileOptions& getOptions()
  { return options_; }

  /// Returns a const reference to the current file options.
  const PeakFileOptions& getOptions() const
  { return options_; }

  /// Replaces the current file options.
  void setOptions(const PeakFileOptions& opts)
  { options_ = opts; }
  //@}

  /** @name Read / Write a complete mass spectrometric experiment
   */
  //@{

  /**
    @brief Load an mzPeak file into an MSExperiment.
  */
  void load(const String& filename, MapType& map) const;

  /**
    @brief Store an MSExperiment in mzPeak format.

    Writes the experiment as point-layout Parquet tables (profile spectra to
    @c spectra_data.parquet, centroid spectra to @c spectra_peaks.parquet) plus
    a per-spectrum @c spectra_metadata.parquet and an @c mzpeak_index.json,
    bundled (STORED) into a ZIP archive. The output round-trips through load()
    to an equivalent experiment (spectrum count, ms_level, type, RT, peaks).

    @note Run-level metadata and precursor facets are not yet emitted.
  */
  void store(const String& filename, const MapType& map) const;

  /**
    @brief Stream the contents of an mzPeak file into a data consumer.

    Performs an optional first pass to count spectra that pass the active
    PeakFileOptions filters and notify the consumer via setExpectedSize() and
    setExperimentalSettings(), then a second pass that calls
    consumeSpectrum() for each passing spectrum.

    @param filename_in    Path to the input .mzpeak archive.
    @param consumer       Receiver for streamed spectra; must not be null.
    @param skip_full_count Ignored (kept for API compatibility).
    @param skip_first_pass When true, the first pass (setExpectedSize /
                           setExperimentalSettings) is skipped.
  */
  void transform(const String& filename_in, Interfaces::IMSDataConsumer* consumer, bool skip_full_count = false, bool skip_first_pass = false) const;
  //@}

  /** @name On-disc random-access interface
   *
   * Provides indexed access to individual spectra without loading the whole file.
   * See the class-level documentation for threading and limitation notes.
   */
  //@{

  /**
    @brief Open an mzPeak file for random-access streaming.

    Reads only the (small) metadata Parquet table into memory; the data and
    peaks Parquet tables are kept open on disc for on-demand access via
    getSpectrum().  Calling openFile() again on the same object closes the
    previously open file first (old state is released atomically only after
    the new file is fully open; a failed re-open leaves the old file intact).

    @param filename  Path to the .mzpeak archive.
    @throw Exception::FileNotFound if the file does not exist.
    @throw Exception::ParseError   if the archive is malformed.
  */
  void openFile(const String& filename);

  /// Returns the number of spectra in the currently open file, or 0 if no file is open.
  Size getNrSpectra() const;

  /**
    @brief Fetch a single spectrum by index without loading the full file.

    Uses Parquet row-group statistics to read only the row group(s) that
    contain spectrum @p index, then filters within those row groups.  Profile
    spectra have null-marked m/z values reconstructed identically to load().

    @note PeakFileOptions are @em not applied — the full spectrum is always
    returned regardless of MS-level or m/z range filters.

    @note This method is @em not @c const because it lazily builds a
    row-group index on first access.

    @param index  Zero-based spectrum index.
    @throw Exception::IllegalArgument  if no file is open (call openFile() first).
    @throw Exception::InvalidParameter if @p index is not in [0, getNrSpectra()).
    @throw Exception::ParseError       if Parquet reading fails.
  */
  MSSpectrum getSpectrum(Size index);

  //@}

private:
  struct OnDiscState; ///< PIMPL: on-disc reader state (defined in .cpp to avoid exposing Arrow types).

  PeakFileOptions options_;              ///< Active read options (MS-level / RT / m/z / metadata-only filters).
  std::unique_ptr<OnDiscState> on_disc_; ///< Non-null while a file is open via openFile().
};
} // namespace OpenMS
