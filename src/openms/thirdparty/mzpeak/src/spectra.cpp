/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#include "mzpeak/spectra.h"

#include <algorithm>
#include <iterator>
#include <limits>
#include <memory>

#include "mzpeak/data/signals.h"
#include "mzpeak/exception.h"
#include "mzpeak/metadata/table.h"
#include "mzpeak/schema/group.h"
#include "mzpeak/spectrum.h"
#include "mzpeak/util/decoders.h"
#include "mzpeak/util/enumerable_proxy.h"
#include "mzpeak/util/projection.h"

namespace MzPeak {

/******************************************************************************/
Spectra::Spectra(std::unique_ptr<Data::Signals> data,
                 std::unique_ptr<Metadata::Table> meta,
                 ImsCalibration ims,
                 std::shared_ptr<const MetadataMap> md,
                 std::shared_ptr<const MinimalMetadataMap> minimal)
    : EnumerableProxy(0)
    , data_(std::move(data))
    , peaks_()
    , meta_(std::move(meta))
    , ims_(ims)
    , md_minimal_(std::move(minimal))
{
  resize(data_->record_count());
  load_metadata_(std::move(md));
  resize_from_metadata_();
}

/******************************************************************************/
Spectra::Spectra(std::unique_ptr<Data::Signals> data,
                 std::unique_ptr<Data::Signals> peaks,
                 std::unique_ptr<Metadata::Table> meta,
                 ImsCalibration ims,
                 std::shared_ptr<const MetadataMap> md,
                 std::shared_ptr<const MinimalMetadataMap> minimal)
    : EnumerableProxy(0)
    , data_(std::move(data))
    , peaks_(std::move(peaks))
    , meta_(std::move(meta))
    , ims_(ims)
    , md_minimal_(std::move(minimal))
{
  // Size from BOTH tables, not just the profile one.  The two normally share a
  // spectrum_count KV, but when that is absent the count falls back to the
  // table's own maximum index — and a run whose last spectra are centroid-only
  // then reports a size that stops short of them.  Iteration and batch reads
  // would skip those spectra while by_id() still found them.
  std::size_t count = data_->record_count();
  if (peaks_) count = std::max(count, peaks_->record_count());
  resize(count);
  load_metadata_(std::move(md));
  resize_from_metadata_();
}

/******************************************************************************/
void Spectra::resize_from_metadata_()
{
  // The signal tables' own counts are a FLOOR, not the answer.
  //
  // `spectrum_count` is a per-FILE key, and a converter is free to stamp each
  // member with the number of spectra IN THAT MEMBER: mzpeak-convert does
  // exactly that from 0.11 on, so a run of 14 profile and 34 centroid spectra
  // carries "14" on spectra_data and "34" on spectra_peaks while holding 48
  // spectra.  Taking the larger of the two then reports 34 and everything
  // after it silently disappears from iteration and from get_spectra_batch --
  // measured on a converted small.mzML, and on a 371,187-spectrum ddaPASEF run
  // that reported 371,184.
  //
  // The metadata table is the one place that counts spectra rather than rows
  // of one array, so a metadata index beyond both tables extends the run.  It
  // can only ever grow the count: a metadata table that is SHORT (a facet
  // misread, a truncated file) must not shrink the reachable set.
  //
  // Indices are dense and 0-based in every archive this reader has seen, and
  // the enumeration contract here is numeric, so a sparse index set yields
  // default-constructed spectra in the gaps -- the same behaviour a file with
  // missing rows already gets.
  // The map is sorted by index and free of duplicates (IndexMap::sort(), which
  // every map is built through before any lookup), so its last entry holds the
  // largest index.  Walking all of them instead read one cache line per
  // spectrum of a map shared by every reader: 3 ms per Spectra on a
  // 717,924-spectrum run, paid again by each per-thread copy.
  uint64_t last = 0;
  if (md_map_ && md_map_->begin() != md_map_->end())
    last = std::prev(md_map_->end())->first;
  else if (md_minimal_ && md_minimal_->begin() != md_minimal_->end())
    last = std::prev(md_minimal_->end())->first;
  else
    return;
  if (last == std::numeric_limits<uint64_t>::max()) return; // last + 1 would wrap
  const std::size_t from_metadata = static_cast<std::size_t>(last) + 1;
  if (from_metadata > size()) resize(from_metadata);
}

/******************************************************************************/
void Spectra::load_metadata_(std::shared_ptr<const MetadataMap> md)
{
  // Adopting a map the caller already has is the whole point of the parameter:
  // it is the difference between one copy of the metadata per ARCHIVE and one
  // per Spectra, and a multi-threaded reader needs one Spectra per thread.
  if (md) {
    md_map_ = std::move(md);
    return;
  }
  if (md_minimal_) return;

  // Built here rather than lazily on first fetch(): almost every access path
  // needs it, and a lazily-published cache was a data race between concurrent
  // fetch() calls.  Peak decode stays lazy — that is the expensive part.
  if (!meta_) return;
  md_map_ =
      std::make_shared<const MetadataMap>(meta_->read_spectrum_metadata());
}

/******************************************************************************/
void Spectra::build_id_index_() const
{
  // First id wins: ids SHOULD be unique, but nothing enforces it, and silently
  // remapping an id to a later spectrum would be worse than ignoring the
  // duplicate.
  std::call_once(id_index_once_, [this] {
    for_each_metadata_([this](uint64_t index, const auto& md) {
      if (md.id.empty()) return;
      id_to_index_.emplace(md.id, static_cast<std::size_t>(index));
    });
  });
}

/******************************************************************************/
std::optional<std::size_t> Spectra::index_for_id(const std::string& id) const
{
  build_id_index_();
  auto it = id_to_index_.find(id);
  if (it == id_to_index_.end()) return std::nullopt;
  return it->second;
}

/******************************************************************************/
Spectrum Spectra::by_id(const std::string& id) const
{
  auto index = index_for_id(id);
  if (!index) throw ParquetError("no spectrum with id '" + id + "'");
  return fetch_(static_cast<uint64_t>(*index));
}

/******************************************************************************/
std::vector<std::size_t> Spectra::indices_in_time_range(double rt_low,
                                                        double rt_high) const
{
  std::vector<std::size_t> result;
  if (!md_map_ && !md_minimal_) return result;
  if (rt_low > rt_high) std::swap(rt_low, rt_high);

  // Collect (time, index) then sort by time: the map is keyed by index, and
  // acquisition order is not guaranteed to be time order.
  std::vector<std::pair<double, std::size_t>> selected;
  for_each_metadata_([&](uint64_t index, const auto& md) {
    if (!md.retention_time.has_value()) return;
    const double t = *md.retention_time;
    if (t < rt_low || t > rt_high) return; // inclusive on both ends
    selected.emplace_back(t, static_cast<std::size_t>(index));
  });

  std::ranges::sort(selected);
  result.reserve(selected.size());
  for (const auto& [t, index] : selected)
    result.push_back(index);
  return result;
}

/******************************************************************************/
std::vector<EicPoint>
Spectra::extract_ion_chromatogram(double mz_low,
                                  double mz_high,
                                  double rt_low,
                                  double rt_high,
                                  std::optional<int> ms_level) const
{
  std::vector<EicPoint> result;
  if (mz_low > mz_high) std::swap(mz_low, mz_high);

  const std::vector<std::size_t> indices = indices_in_time_range(rt_low, rt_high);
  result.reserve(indices.size());

  for (std::size_t index : indices) {
    double time = 0.0;
    if (md_map_ || md_minimal_) {
      bool wrong_level = false;
      const bool found = with_metadata_(static_cast<uint64_t>(index), [&](const auto& md) {
        // ms-level filter is decided from metadata, before any peak decode.
        wrong_level = ms_level.has_value() && md.ms_level != ms_level;
        time = md.retention_time.value_or(0.0);
      });
      if (wrong_level || (!found && ms_level.has_value())) continue;
    }

    Spectrum spectrum = fetch_(static_cast<uint64_t>(index));
    const auto& mz = spectrum.mz();
    const auto& intensity = spectrum.intensity();

    // Scan the whole array rather than breaking at the first out-of-window
    // point: m/z is usually ascending, but the layout does not guarantee it
    // (sorting_rank may be absent), and an early break would undercount.
    const std::size_t n = std::min(mz.size(), intensity.size());
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      if (mz[i] >= mz_low && mz[i] <= mz_high) {
        sum += static_cast<double>(intensity[i]);
      }
    }

    // Dense trace: a scan with no signal in the window still contributes a
    // zero sample, so gaps stay visible instead of silently closing up.
    result.push_back(EicPoint{time, sum, index});
  }

  return result;
}

/******************************************************************************/
std::vector<Spectrum>
Spectra::get_spectra_batch(const std::vector<std::size_t>& indices) const
{
  std::vector<Spectrum> result(indices.size());

  // Read in ascending index order so file access is sequential, but write each
  // result back to its ORIGINAL position so callers keep positional
  // correspondence with the indices they passed.
  std::vector<std::size_t> order(indices.size());
  for (std::size_t i = 0; i < order.size(); ++i)
    order[i] = i;
  std::ranges::sort(order, [&indices](std::size_t a, std::size_t b) {
    return indices[a] < indices[b];
  });

  for (std::size_t slot : order) {
    const std::size_t index = indices[slot];
    if (index >= size()) continue; // out of range -> default Spectrum
    result[slot] = fetch_(static_cast<uint64_t>(index));
  }

  return result;
}

/******************************************************************************/
Spectrum Spectra::fetch_(uint64_t index) const
{
  // Choose the table this spectrum lives in.  Read from the cached metadata map
  // rather than re-querying the metadata file per spectrum: the old projection
  // needed a `spectrum` struct group, which the flat layout does not have, so it
  // silently fell through to the profile table.
  //
  // The spectrum REPRESENTATION decides, not merely "does it have peaks".  A
  // spectrum may carry both a profile array and a centroid array — in
  // small.chunked.mzpeak spectrum 0 declares number_of_data_points 13589 AND
  // number_of_peaks 1612 while being profile (MS:1000128).  Dispatching on the
  // peak count alone handed back 1612 centroids for a profile spectrum: the
  // wrong array, silently, and with a plausible length.
  std::shared_ptr<Data::Signals> signals = data_;
  if (peaks_) {
    with_metadata_(index, [&](const auto& md) {
      const bool centroid = md.representation == "MS:1000127";
      const bool profile = md.representation == "MS:1000128";

      if (centroid) {
        signals = peaks_;
      } else if (!profile) {
        // No usable representation: fall back to "wherever the data is".
        if (md.number_of_data_points.value_or(0) == 0 &&
            md.number_of_peaks.value_or(0) > 0) {
          signals = peaks_;
        }
      }
    });
  }

  // A plain copy_if rather than views::filter | ranges::to: ranges::to is a
  // C++23 LIBRARY feature Apple's clang 15 libc++ lacks.
  std::vector<Data::ArrayIndex::Dimension> dims;
  std::ranges::copy_if(signals->array_index()->dimensions(), std::back_inserter(dims),
                       [](auto& d) {
        // Mobility must be selected here too, or the column is never
        // projected and ion_mobility_array() comes back empty however well the
        // decoder handles it.
        return d.array_type == Schema::PSI::ArrayType::Mz ||
               d.array_type == Schema::PSI::ArrayType::Intensity ||
               Schema::PSI::is_ion_mobility(d.array_type) ||
               d.name.find("mobility") != std::string::npos ||
               // ims-compact keeps m/z in a non-standard `tof` column; without
               // selecting it here the spectrum has no coordinate at all.
               d.name.find("tof") != std::string::npos;
                       });

  // A Minimal entry is expanded here, once per Spectrum: that is the
  // SpectrumMetadata its metadata() hands out, and what the decode reads.
  std::shared_ptr<const SpectrumMetadata> own;
  if (md_minimal_) {
    auto it = md_minimal_->find(index);
    if (it != md_minimal_->end())
      own = std::make_shared<const SpectrumMetadata>(it->second.expand());
  }
  return Spectrum(index, signals, std::move(dims), meta_, md_map_, ims_, std::move(own));
}

} // namespace MzPeak
