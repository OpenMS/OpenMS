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

#include "mzpeak/ims_calibration.h"
#include "mzpeak/spectrum.h"
#include "mzpeak/spectrum_metadata.h"
#include "mzpeak/util/enumerable_proxy.h"
#include "mzpeak/util/index_map.h"

// Forward declarations:
namespace MzPeak::Data {
class Signals;
} // namespace MzPeak::Data

namespace MzPeak::Metadata {
class Table;
}

namespace MzPeak {

/**
 * One sample of an extracted-ion chromatogram: the summed intensity within an
 * m/z window for a single spectrum, tagged with that spectrum's retention time
 * and index.  An EIC is the ascending-time sequence of these points.
 */
struct EicPoint final {
  /// The spectrum's retention time, in SECONDS.
  double time = 0.0;

  /// Summed intensity of the points falling inside the m/z window.
  double intensity = 0.0;

  /// Index of the spectrum this sample came from.
  std::size_t spectrum_index = 0;
};

/**
 * Access all spectra in a MzPeak file.
 *
 * Iteration is streaming: spectra are fetched on demand and their peak arrays
 * decode lazily, so a run is traversed in bounded memory.  The selection
 * helpers below resolve from the cached descriptive metadata and therefore do
 * NOT decode peaks — only extract_ion_chromatogram() does, and only for the
 * spectra it actually selects.
 */
class Spectra final : public Util::EnumerableProxy<Spectrum> {
public:
  /// Default constructor: empty (zero spectra).
  Spectra() = default;

  /// Neither copyable nor movable.
  ///
  /// The base class stores a fetch callback that binds `this`.  With implicit
  /// copy/move that callback keeps pointing at the ORIGINAL object, so after
  /// `auto b = a;` every `b[i]` dispatches through `a` — and dangles once `a`
  /// dies, or silently reads a different run after `a` is reassigned.  Nothing
  /// about the copy looks wrong at the call site.
  ///
  /// Deleting them costs nothing in practice: `Spectra s = index.spectra();`
  /// initialises directly from a prvalue and is elided, which is how every
  /// caller obtains one.  It also makes the "use one Spectra per thread" rule
  /// enforceable rather than advisory, since a Spectra can no longer be copied
  /// into a thread.
  Spectra(const Spectra&) = delete;
  Spectra& operator=(const Spectra&) = delete;
  Spectra(Spectra&&) = delete;
  Spectra& operator=(Spectra&&) = delete;

  /// The cached per-spectrum descriptive metadata, keyed by spectrum index.
  using MetadataMap = Util::IndexMap<SpectrumMetadata>;
  /// MetadataDetail::Minimal's compact map; a Spectra holds this or the above.
  using MinimalMetadataMap = Util::IndexMap<MinimalSpectrumMetadata>;

  /// Constructor for profile-only or centroid-only data.
  ///
  /// @param md  a metadata map already built for this archive, to be SHARED
  ///   rather than re-read.  It is the largest fixed cost of opening a run
  ///   (26.7 MB on a 7,534-spectrum Thermo file), it is immutable once built,
  ///   and the only safe way to read one archive from N threads is one
  ///   `Spectra` per thread -- so without sharing, N threads mean N copies of
  ///   one table.  Null means "read it from @p meta", which is what a caller
  ///   with a single reader wants.
  /// @param minimal  the MetadataDetail::Minimal map, held in place of @p md.
  explicit Spectra(std::unique_ptr<Data::Signals>,
                   std::unique_ptr<Metadata::Table>,
                   ImsCalibration ims = {},
                   std::shared_ptr<const MetadataMap> md = nullptr,
                   std::shared_ptr<const MinimalMetadataMap> minimal = nullptr);

  /// Constructor for mixed profile+centroid data (data = profile file, peaks =
  /// centroid file).  @p md as above.
  explicit Spectra(std::unique_ptr<Data::Signals> data,
                   std::unique_ptr<Data::Signals> peaks,
                   std::unique_ptr<Metadata::Table>,
                   ImsCalibration ims = {},
                   std::shared_ptr<const MetadataMap> md = nullptr,
                   std::shared_ptr<const MinimalMetadataMap> minimal = nullptr);

  /**
   * Resolve a native spectrum id (e.g. "controllerType=0 controllerNumber=1
   * scan=1") to its index, or nullopt when no spectrum carries that id.
   */
  std::optional<std::size_t> index_for_id(const std::string& id) const;

  /**
   * Fetch the spectrum with the given native id.
   * @throws ParquetError when the id is unknown.
   */
  Spectrum by_id(const std::string& id) const;

  /**
   * Indices of every spectrum whose retention time lies in the INCLUSIVE range
   * [rt_low, rt_high], in ascending time order.  Spectra with no retention time
   * are skipped.
   *
   * @note Times are SECONDS, matching SpectrumMetadata::retention_time — the
   * file stores minutes, and the reader converts.  Passing minutes here is a
   * 60x error that silently selects the wrong scans.
   *
   * @note An inverted range (`rt_low > rt_high`) is NORMALISED rather than
   * rejected: the bounds are swapped, so the call returns the same window as
   * the correctly ordered one.
   */
  std::vector<std::size_t> indices_in_time_range(double rt_low,
                                                 double rt_high) const;

  /**
   * Extracted-ion chromatogram over an m/z window and an RT window (both
   * inclusive, RT in SECONDS), optionally restricted to one MS level.
   *
   * Emits one point per selected spectrum in ascending time order.  A spectrum
   * with no signal in the m/z window still emits a point with intensity 0 — an
   * EIC is a dense trace over the selected scans, so gaps must stay visible
   * rather than being dropped.
   *
   * This is the only selection helper that decodes peaks, and it decodes only
   * the spectra that survive the RT and ms-level filters.
   *
   * @note An inverted range is NORMALISED by swapping rather than rejected —
   * for m/z exactly as indices_in_time_range() does for time.
   */
  std::vector<EicPoint>
  extract_ion_chromatogram(double mz_low,
                           double mz_high,
                           double rt_low,
                           double rt_high,
                           std::optional<int> ms_level = std::nullopt) const;

  /**
   * Read several spectra by index.  Reads proceed in ascending index order so
   * the underlying file access is sequential, but results are returned in the
   * SAME order as @p indices.  An out-of-range index yields a default Spectrum.
   */
  std::vector<Spectrum>
  get_spectra_batch(const std::vector<std::size_t>& indices) const;

private:
  // Internal data access.
  std::shared_ptr<Data::Signals> data_;
  std::shared_ptr<Data::Signals> peaks_; // centroid-only file (optional)
  std::shared_ptr<Metadata::Table> meta_;

  // TOF -> m/z calibration, passed to each Spectrum (ims-compact layout only).
  ImsCalibration ims_;

  // Cached per-spectrum descriptive metadata (RT, precursors, ion mobility, …),
  // read once at construction and shared into every Spectrum.  Stays null when
  // there is no metadata table.
  std::shared_ptr<const MetadataMap> md_map_;
  std::shared_ptr<const MinimalMetadataMap> md_minimal_;

  // Native id -> index, for by_id().  Ids are not guaranteed unique; the FIRST
  // spectrum carrying an id wins.
  //
  // Built on first use, not with md_map_: it duplicates every native id string
  // (~48 bytes each plus a tree node), and iteration, RT selection and EIC
  // extraction -- what a reader that streams a run actually does -- never touch
  // it.  `call_once` keeps that lazy publication safe, which is why the map
  // itself is still eager (a lazily-published md_map_ WAS a data race).
  mutable std::map<std::string, std::size_t> id_to_index_;
  mutable std::once_flag id_index_once_;

  // Populate md_map_ from meta_, or adopt the shared map the caller passed.
  // Called from the constructors so that fetch() never publishes it lazily.
  void load_metadata_(std::shared_ptr<const MetadataMap> md);

  /// f(index, entry) for every entry of whichever map this holds.  The two
  /// entry types share the fields these loops read, by name and type.
  template <typename F> void for_each_metadata_(F&& f) const
  {
    if (md_map_)
      for (const auto& [index, md] : *md_map_) f(index, md);
    else if (md_minimal_)
      for (const auto& [index, md] : *md_minimal_) f(index, md);
  }

  /// f(entry) for @p index's entry, if either map has one; whether it did.
  template <typename F> bool with_metadata_(uint64_t index, F&& f) const
  {
    if (md_map_) {
      auto it = md_map_->find(index);
      if (it == md_map_->end()) return false;
      f(it->second);
      return true;
    }
    if (md_minimal_) {
      auto it = md_minimal_->find(index);
      if (it == md_minimal_->end()) return false;
      f(it->second);
      return true;
    }
    return false;
  }

  // Grow the run to cover the metadata's highest spectrum index.  Called from
  // the constructors AFTER load_metadata_ and on BOTH of its paths -- a caller
  // supplying a cached map takes an early return inside it, which is the
  // ordinary path through Index::spectra().  See the definition: a signal
  // table's per-file spectrum_count is a floor, not the run's size.
  void resize_from_metadata_();

  // Build id_to_index_ once, on first by_id()/index_for_id().
  void build_id_index_() const;

  // Function to fetch a specific spectrum.  Const because it mutates no Spectra
  // state — the metadata cache is built in the constructor — which lets the
  // selection helpers above be const without casting.
  Spectrum fetch_(uint64_t) const;

  /// EnumerableProxy's fetch is a non-const pure virtual; ours is const so the
  /// const public API (by_id, range queries) still works.  Forward.
  Spectrum fetch(uint64_t n) override { return fetch_(n); }
};

} // namespace MzPeak
