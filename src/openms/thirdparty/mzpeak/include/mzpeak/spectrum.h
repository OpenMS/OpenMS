/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

#include "mzpeak/data/array_index.h"
#include "mzpeak/data/encoding.h"
#include "mzpeak/data/signals.h"
#include "mzpeak/ims_calibration.h"
#include "mzpeak/metadata/spectrum.h"
#include "mzpeak/metadata/table.h"
#include "mzpeak/spectrum_metadata.h"
#include "mzpeak/util/index_map.h"
#include "mzpeak/util/slice.h"

namespace MzPeak {

// Forward declaration.
class Spectra;

/**
 * Access to a single spectrum in an MzPeak file.
 *
 * @note THREAD SAFETY.  A `Spectrum` is safe to *read* concurrently: the lazy
 * peak decode behind mz()/intensity() is serialised with `std::call_once`, and
 * copies of a `Spectrum` share one decode, so the arrays are materialised
 * exactly once no matter how many copies or threads ask for them.
 *
 * What is NOT safe is driving one `Spectra` (or the `Data::Signals` behind it)
 * from several threads at once: decoding two DIFFERENT spectra concurrently
 * still funnels through the same Parquet reader, whose thread-safety is not
 * established.  **Use one `Spectra` per thread.**
 */
class Spectrum final {
public:
  /// The type of decoder used.
  using decoder_type = Data::Encoding::Decoder<double>;

  /// Default constructor: an empty spectrum with no backing file.  mz() and
  /// intensity() return empty arrays and metadata() the default record; used
  /// for out-of-range slots in Spectra::get_spectra_batch().
  Spectrum();

  // Other special members are implicit (rule of zero): copyable AND movable.
  // The decoded peak arrays live in a shared_ptr, so copying is cheap and
  // copies share the decode rather than each repeating the file read.

  /**
   * Mass-to-charge values.  Decodes the peak arrays on first access (lazy):
   * calling only the metadata accessors below never touches the peak data.
   *
   * @warning The returned reference is owned by this `Spectrum`.  Because
   * `Spectra::operator[]` yields a Spectrum BY VALUE, writing
   * `const auto& mz = spectra[0].mz();` dangles — bind the spectrum to a named
   * local first: `auto s = spectra[0]; const auto& mz = s.mz();`
   */
  const std::vector<double>& mz() const;

  /**
   * Intensity values.  Lazily decoded (see mz()).
   */
  const std::vector<float>& intensity() const;

  /**
   * Per-peak ion mobility, parallel to mz() and intensity().  Lazily decoded.
   *
   * EMPTY unless the file stores a mobility array.  diaPASEF files keep a frame
   * as ONE spectrum carrying several isolation windows over disjoint mobility
   * ranges plus this array, so a per-spectrum scalar cannot represent it:
   * collapsing to one value gives every peak in the frame the same drift time
   * and lets the windows bleed together.  Use ion_mobility() for the older
   * one-scalar-per-spectrum layout.
   */
  const std::vector<double>& ion_mobility_array() const;

  /**
   * Stage number achieved in a multi stage mass spectrometry
   * acquisition.  Read from metadata; does NOT decode peaks.
   */
  uint8_t ms_level() const;

  /**
   * Full per-spectrum descriptive metadata (RT, precursors, isolation windows,
   * selected ions, scan windows, ion mobility).  Resolved from the cached
   * spectra_metadata table WITHOUT decoding the peak arrays.
   *
   * @warning Same lifetime rule as mz(): the reference is kept alive by this
   * `Spectrum` (and the `Spectra` that produced it), so
   * `const auto& md = run.spectra()[3].metadata();` dangles — both temporaries
   * die at the end of the statement.  Bind a named local first.
   */
  const SpectrumMetadata& metadata() const;

  /// Retention time in SECONDS (handoff P0).  See SpectrumMetadata::retention_time.
  std::optional<double> retention_time() const { return metadata().retention_time; }

  /// Precursor list (empty for MS1); isolation window is target_mz ± offsets.
  const std::vector<PrecursorInfo>& precursors() const
  {
    return metadata().precursors;
  }

  /// Scan-level ion mobility value (1/K0 for diaPASEF/timsTOF) — handoff P1.
  std::optional<double> ion_mobility() const { return metadata().ion_mobility; }

  /// Scan-level ion mobility type — handoff P1.
  std::optional<std::string> ion_mobility_type() const
  {
    return metadata().ion_mobility_type;
  }

protected:
  friend class Spectra;

  /// Internal constructor.  The peak arrays (and the signal-file read that backs
  /// them) are deferred to first mz()/intensity(); @p md_map is the shared,
  /// cached spectra_metadata map (may be null when no metadata file).
  Spectrum(uint64_t index,
           std::shared_ptr<Data::Signals>,
           std::vector<Data::ArrayIndex::Dimension>,
           std::shared_ptr<Metadata::Table>,
           std::shared_ptr<const Util::IndexMap<SpectrumMetadata>> md_map,
           ImsCalibration ims = {},
           std::shared_ptr<const SpectrumMetadata> md_own = nullptr);

private:
  /// The lazily-decoded peak arrays, plus the flag that serialises the decode.
  /// Held behind a shared_ptr so that copying a Spectrum is cheap, copies share
  /// one decode, and `std::once_flag` (neither copyable nor movable) does not
  /// leak into Spectrum's own special members.
  struct Peaks {
    std::once_flag once;
    std::vector<double> mz;
    std::vector<float> intensity;
    std::vector<double> mobility;
  };

  /// Read + decode the peak arrays into peaks_ on first access.  Runs at most
  /// once per Spectrum (and its copies) even under concurrent callers.  This is
  /// the ONLY place the signal (peak) file is touched, so metadata-only
  /// iteration never reads peak data.
  void decode_() const;

  uint64_t index_;
  std::shared_ptr<Metadata::Table> md_table_;
  std::shared_ptr<const Util::IndexMap<SpectrumMetadata>> md_map_;
  /// This spectrum's own metadata, when its Spectra holds a Minimal map.
  std::shared_ptr<const SpectrumMetadata> md_own_;

  // Ingredients retained for the lazy peak read + decode.
  std::shared_ptr<Data::Signals> signals_;
  std::vector<Data::ArrayIndex::Dimension> dims_;

  // TOF -> m/z calibration; only used by the ims-compact layout, which stores
  // no m/z array at all.
  ImsCalibration ims_;

  std::shared_ptr<Peaks> peaks_;
};

} // namespace MzPeak
