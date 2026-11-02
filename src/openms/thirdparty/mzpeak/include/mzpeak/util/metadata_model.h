/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "mzpeak/chromatogram_metadata.h"
#include "mzpeak/spectrum_metadata.h"
#include "mzpeak/util/index_map.h"
#include "mzpeak/util/parquet.h"
#include "mzpeak/wavelength_spectrum_metadata.h"

namespace MzPeak::Util {

/**
 * Read the per-spectrum `mz_delta_model` (the m/z null-marking delta model)
 * from a spectra_metadata Parquet table, keyed by `spectrum.index`.
 *
 * Spectra without a model are omitted.  Used to reconstruct null-marked
 * profile m/z values (see null_fill.h).
 */
std::map<uint64_t, std::vector<double>> read_mz_delta_models(Parquet& metadata);

/**
 * Read per-spectrum descriptive metadata (top-level `spectrum` struct fields,
 * plus precursor / selected_ion / scan facets) from a spectra_metadata Parquet
 * table, keyed by `spectrum.index`.
 *
 * RDR-10b: exposes spectrum_type (MS:1000559), lowest/highest observed m/z,
 * data_processing_ref, and the per-spectrum parameters CvParam list.
 *
 * RDR-10c: reads the Parquet table ONCE and joins precursor, selected_ion, and
 * scan columns into each SpectrumMetadata entry by `source_index` VALUE (H1
 * fix — never by chunk-local/positional row index; chunk boundaries differ per
 * column).  Selected ions are attached to their owning PrecursorInfo by
 * explicit (source_index, precursor_index) matching (H2 fix — DIA-safe for
 * multiple precursors per spectrum and multiple selected ions per precursor).
 * Unmatched source_index values are logged to stderr (no silent loss).
 *
 * Ion-mobility fields (ion_mobility_value/type) are read from both the
 * selected_ion struct (SelectedIonInfo) and the scan struct (SpectrumMetadata::
 * ion_mobility/ion_mobility_type).  They are NULL in every bundled fixture, so
 * the decode is exercised for null-safety only; value-level correctness is
 * fixture-gated on a real diaPASEF/timsTOF run.
 *
 * Retention time (SpectrumMetadata::retention_time) is in SECONDS: sourced from
 * scan.MS_1000016_scan_start_time (UO_0000031 = minutes) ×60, with spectrum.time
 * as a fallback (also minutes).
 *
 * Nullable source values map to absent std::optional fields or empty
 * strings/vectors as documented on SpectrumMetadata.
 */
IndexMap<SpectrumMetadata> read_spectra_metadata(Parquet& metadata);

/**
 * The Parquet files that together carry one entity's metadata.
 *
 * Older writers nest every facet as a struct column of @ref primary, so the
 * facet members stay null.  Newer writers give each facet its own file with
 * flat columns; those are joined to @ref primary by `source_index` VALUE.
 *
 * Non-owning: the caller keeps the Parquet objects alive for the call.
 */
struct SpectraMetadataFiles {
  Parquet* primary = nullptr;
  Parquet* scans = nullptr;
  Parquet* precursors = nullptr;
  Parquet* selected_ions = nullptr;
};

/**
 * Read per-spectrum metadata from a primary file plus any separate facet files.
 * Handles both the nested single-table layout and the split-file layout; see
 * @ref SpectraMetadataFiles.
 */
IndexMap<SpectrumMetadata>
read_spectra_metadata(const SpectraMetadataFiles&);

/**
 * As above, but materialising only what @p detail asks for.  See
 * @ref MzPeak::MetadataDetail: `Lean` omits the CV-parameter lists, scan
 * windows and auxiliary arrays and keeps everything a decode or a selection
 * query reads.
 */
IndexMap<SpectrumMetadata>
read_spectra_metadata(const SpectraMetadataFiles&, MetadataDetail detail);

/**
 * The same read at @ref MzPeak::MetadataDetail::Minimal, into compact
 * records.  nullopt when a spectrum has more than one precursor row, more
 * than one selected ion, or an ion whose precursor_index differs from that of
 * its precursor row: read the archive as `Lean` then.  An ion on a spectrum
 * with no precursor row is held, under a precursor of its own.
 */
std::optional<IndexMap<MinimalSpectrumMetadata>>
read_minimal_spectra_metadata(const SpectraMetadataFiles&);

/**
 * From this many spectra on, the spectrum metadata readers above attach
 * precursors and selected ions in four contiguous index ranges, three of them
 * on threads of their own; below it, in one range on the calling thread.  The
 * result is the same either way.
 */
inline constexpr std::size_t kParallelAttachMinSpectra = 1024;

/**
 * The Parquet files carrying one chromatogram set's metadata.  As with
 * @ref SpectraMetadataFiles, the facets are either struct columns of @ref
 * primary or files of their own.  Non-owning.
 *
 * There is no `products` member: see @ref ChromatogramMetadata.
 */
struct ChromatogramMetadataFiles {
  Parquet* primary = nullptr;
  Parquet* precursors = nullptr;
  Parquet* selected_ions = nullptr;
};

/**
 * Read per-chromatogram metadata, keyed by `chromatogram.index`.  Handles the
 * nested single-table and split-file layouts.
 */
std::map<uint64_t, ChromatogramMetadata>
read_chromatogram_metadata(const ChromatogramMetadataFiles&);

/**
 * The Parquet files carrying one wavelength-spectrum set's metadata.
 * Non-owning.  There is no precursor or selected-ion facet for this entity.
 */
struct WavelengthMetadataFiles {
  Parquet* primary = nullptr;
  Parquet* scans = nullptr;
};

/**
 * Read per-wavelength-spectrum metadata, keyed by `spectrum.index`.  Handles
 * the nested single-table and split-file layouts.
 */
std::map<uint64_t, WavelengthSpectrumMetadata>
read_wavelength_spectrum_metadata(const WavelengthMetadataFiles&);

/**
 * Read the native-id → entity-index map from a metadata Parquet table whose
 * top-level column is named @p col_name (e.g. "chromatogram" or
 * "wavelength_spectrum").  Only rows with a non-null, non-empty `id` field are
 * inserted.  Used by Chromatograms::by_id() and WavelengthSpectra::by_id().
 */
std::unordered_map<std::string, std::size_t>
read_entity_id_map(Parquet& metadata, const std::string& col_name);

} // namespace MzPeak::Util
