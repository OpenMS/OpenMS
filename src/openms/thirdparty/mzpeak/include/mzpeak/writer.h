/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "mzpeak/run_metadata.h"

namespace MzPeak {

/**
 * One selected ion of a precursor.  Every field is optional because every one
 * is nullable in the format; a DIA window commonly has a target and no ion.
 */
struct SelectedIonData {
  std::optional<double> mz;
  std::optional<int> charge;
  std::optional<float> intensity;
};

/**
 * One precursor of an MS2+ spectrum: its isolation window and the ions
 * selected from it.  A DIA frame carries several of these per spectrum.
 *
 * The window bounds are float because that is what the format stores
 * (MS:1000827-9 are float32 columns in every archive this library has seen);
 * a double here would promise precision the file cannot keep.
 *
 * Activation parameters and a precursor id are not carried yet -- nothing
 * this writer serves needs them, and an empty list is honest where an
 * invented one is not.
 */
struct PrecursorData {
  std::optional<float> isolation_target_mz;
  std::optional<float> isolation_lower_offset;
  std::optional<float> isolation_upper_offset;
  std::vector<SelectedIonData> selected_ions;
};

/**
 * The m/z and intensity arrays of a single spectrum to be written, with the
 * per-spectrum metadata the reader hands back.
 */
struct SpectrumData {
  std::vector<double> mz;
  std::vector<float> intensity;

  /// `true` for a centroid (peak) spectrum (written to the peaks table),
  /// `false` for a profile spectrum (written to the data-arrays table).
  bool centroid = false;

  /// MS level (1 = MS1, 2 = MS2, etc). Defaults to 1.
  uint8_t ms_level = 1;

  /// Optional retention time in seconds.
  std::optional<double> retention_time;

  /// Optional ion mobility of the whole spectrum (scan level; e.g. a timsTOF frame's 1/K0 or a FAIMS
  /// compensation voltage), with the CV accession of its kind ("MS:1002815" inverse reduced ion
  /// mobility, "MS:1002476" ion mobility drift time, "MS:1001581" FAIMS compensation voltage).
  /// Written to the scans facet as ion_mobility_value / ion_mobility_type, where the reader finds it
  /// (Spectrum::ion_mobility(), ion_mobility_type()). (OpenMS addition.)
  std::optional<double> ion_mobility;
  std::optional<std::string> ion_mobility_type;

  /// Optional scan polarity: 1 = positive, -1 = negative.
  std::optional<int> polarity;

  /// Optional spectrum ID string; auto-generated ("index=N") when absent.
  std::optional<std::string> id;

  /// Precursors, in order; empty for MS1.  Written to the precursor and
  /// selected-ion facets, keyed by this spectrum's index and the precursor's
  /// position here, which is how the reader joins them back.
  std::vector<PrecursorData> precursors;
};

/**
 * One chromatogram to be written.
 *
 * @note There is no `products` member.  An SRM/MRM transition's Q3 isolation
 * window belongs in a product facet that no reference implementation writes and
 * whose reader path is unimplemented, so writing a chromatogram whose type
 * implies a product would produce a file that silently lacks the transition.
 * The writer REJECTS those types rather than emitting one; see
 * @ref write_run_directory.
 */
struct ChromatogramData {
  /// Time values in SECONDS.  Stored as minutes, which the format expects; the
  /// conversion happens in the writer so callers use one time base throughout.
  std::vector<double> time;

  std::vector<float> intensity;

  /// MS:1000626 chromatogram type CURIE, e.g. "MS:1000235" (total ion current)
  /// or "MS:1000812" (absorption chromatogram).  Empty is permitted but leaves
  /// the chromatogram unidentifiable.
  std::string chromatogram_type;

  /// Unit CURIE for @ref intensity.  Defaults to detector counts; an optical
  /// absorbance trace is UO:0000269 and a reader cannot tell them apart
  /// afterwards without this.
  std::string intensity_unit = "MS:1000131";

  /// Scan polarity: 1 positive, -1 negative.  Absent means unknown -- unlike
  /// the reference writer, which records 0 for unknown.
  std::optional<int> polarity;

  /// Native identifier; auto-generated ("chromatogram=N") when absent.
  std::optional<std::string> id;
};

/**
 * One UV/Vis (wavelength) spectrum to be written.
 */
struct WavelengthSpectrumData {
  /// Wavelength values in NANOMETRES.
  std::vector<float> wavelength;

  std::vector<float> intensity;

  /// Acquisition time in SECONDS; stored as minutes.
  std::optional<double> time;

  /// Unit CURIE for @ref intensity.  UV detectors commonly report absorbance
  /// (UO:0000269) rather than detector counts.
  std::string intensity_unit = "MS:1000131";

  /// Native identifier; auto-generated ("wavelength_spectrum=N") when absent.
  std::optional<std::string> id;
};

/**
 * Everything one mzPeak archive holds.  A real run mixes entity types -- mass
 * spectra alongside a TIC and a diode-array trace -- so they are written
 * together rather than through one entry point per type.
 */
struct RunContents {
  std::vector<SpectrumData> spectra;
  std::vector<ChromatogramData> chromatograms;
  std::vector<WavelengthSpectrumData> wavelength_spectra;
};

/**
 * Write a whole run as an unpacked DIRECTORY.
 *
 * Emits only the members the contents require, plus `mzpeak_index.json`.
 * Metadata uses the split (flat column + `column_mapping`) layout, which is
 * what the current reference reader resolves through; the nested layout is
 * readable but the reference cannot open it as metadata.
 *
 * Summary fields for wavelength spectra (lowest/highest observed wavelength,
 * lambda max, base peak intensity, total ion current) are computed from the
 * data actually written.  The reference writer derives its range from the
 * unsorted input after sorting a copy for output, and seeds its maximum at
 * zero so an all-negative absorbance spectrum records a maximum of 0; neither
 * is reproduced here.
 *
 * @throws ParquetError on I/O or encoding errors, if any entity's parallel
 *         arrays differ in length, or if a chromatogram's type implies an
 *         SRM/MRM product selection this format cannot carry.
 */
/// The streaming counterpart of write_run_archive() for spectra: append one
/// at a time and they land on disk row group by row group, so memory is one
/// row group of points plus a few hundred bytes of metadata per spectrum,
/// never the run.  finish() writes the metadata tables and the index, seals
/// the archive and removes the working files; a writer destroyed unfinished
/// removes them too and leaves no archive.  Spectra only: chromatograms and
/// wavelength spectra still go through write_run_archive().
class RunArchiveWriter final {
public:
  /// @param points_per_row_group  points per Parquet row group, the unit a
  ///   reader decodes; 1,048,576 matches the reference converter.
  explicit RunArchiveWriter(const std::filesystem::path& zip_path,
                            std::size_t points_per_row_group = std::size_t(1) << 20);
  ~RunArchiveWriter();
  RunArchiveWriter(const RunArchiveWriter&) = delete;
  RunArchiveWriter& operator=(const RunArchiveWriter&) = delete;

  /// Append the next spectrum; its index is the number appended before it.
  /// Validated like write_run_archive(): equal mz/intensity lengths, no NaN
  /// m/z; points are stored in ascending m/z.
  void add(SpectrumData spectrum);
  /// Run-level metadata for the index; any time before finish().
  void set_metadata(const RunMetadata& metadata);
  /// Spectra appended so far.
  std::size_t size() const;
  /// Seal the archive.  Idempotent.
  void finish();

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

void write_run_directory(const std::filesystem::path& dir,
                         const RunContents& contents,
                         const RunMetadata* run_metadata = nullptr);

/**
 * Write a whole run as a single ZIP ARCHIVE (`.mzpeak`), with STORED members
 * as the specification mandates.  See @ref write_run_directory.
 */
void write_run_archive(const std::filesystem::path& zip_path,
                       const RunContents& contents,
                       const RunMetadata* run_metadata = nullptr);

/**
 * Write a point-layout mzPeak file as an unpacked DIRECTORY.
 *
 * Produces, inside `dir`: `spectra_data.parquet` (point-layout profile
 * spectra), `spectra_peaks.parquet` (centroid spectra, only when any are
 * present), `spectra_metadata.parquet` (per-spectrum scalar metadata) and
 * `mzpeak_index.json`.  The result is readable by `MzPeak::open(dir)`.
 *
 * Still point-layout only: the chunked layout, numpress and other transforms
 * are not yet emitted.  Run-level index metadata blocks (`run`,
 * `software_list`, …) ARE emitted when the @ref RunMetadata overload is used
 * (WRT-2); this overload emits a `metadata{}` carrying only `version`.
 *
 * @throws on I/O or encoding errors, or if any spectrum's mz and
 *         intensity arrays differ in length.
 */
void write_spectra_directory(const std::filesystem::path& dir,
                             const std::vector<SpectrumData>& spectra);

/**
 * WRT-2 — as @ref write_spectra_directory, additionally emitting the run-level
 * metadata blocks of `metadata` into the index's `metadata{}` object
 * (alongside `version`).  The emitted blocks round-trip back through
 * `Index::metadata()`.
 */
void write_spectra_directory(const std::filesystem::path& dir,
                             const std::vector<SpectrumData>& spectra,
                             const RunMetadata& metadata);

/**
 * Write a point-layout mzPeak file as a single ZIP ARCHIVE (`.mzpeak`).
 *
 * Produces an archive at `zip_path` containing `spectra_data.parquet`
 * (point layout) and `mzpeak_index.json` as STORED (uncompressed)
 * members, as the mzPeak spec mandates.  The result is readable by
 * `MzPeak::open(zip_path)`.
 *
 * The spectra are flattened, validated and per-spectrum m/z sorted
 * identically to @ref write_spectra_directory; only the container differs.
 * Any existing file at `zip_path` is truncated.
 *
 * @throws ParquetError on I/O, encoding or libzip errors, or if any
 *         spectrum's mz and intensity arrays differ in length.
 */
void write_spectra_archive(const std::filesystem::path& zip_path,
                           const std::vector<SpectrumData>& spectra);

/**
 * WRT-2 — as @ref write_spectra_archive, additionally emitting the run-level
 * metadata blocks of `metadata` into the index's `metadata{}` object
 * (alongside `version`).  The emitted blocks round-trip back through
 * `Index::metadata()`.
 */
void write_spectra_archive(const std::filesystem::path& zip_path,
                           const std::vector<SpectrumData>& spectra,
                           const RunMetadata& metadata);

} // namespace MzPeak
