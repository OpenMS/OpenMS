/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/index.h"

#include <arrow/util/key_value_metadata.h>
#include <charconv>
#include <memory>
#include <string_view>

#include "mzpeak/data/signals.h"
#include "mzpeak/exception.h"
#include "mzpeak/io/archive.h"
#include "mzpeak/metadata/table.h"
#include "mzpeak/util/manager.h"
#include "mzpeak/util/sha512.h"
#include "mzpeak/util/metadata_model.h"

namespace MzPeak {

/******************************************************************************/
Index::Index(std::unique_ptr<MzPeak::IO::Archive> archive)
    : manager_(std::make_shared<Util::Manager>(std::move(archive)))
{
}

/******************************************************************************/
const std::vector<Schema::File>& Index::files() const { return manager_->files(); }

/******************************************************************************/
ChecksumReport Index::verify_checksums() const
{
  ChecksumReport report;

  for (const Schema::File& file : manager_->files()) {
    const std::optional<std::string>& expected = file.checksum();
    if (!expected.has_value() || expected->empty()) {
      ++report.unchecked;
      continue;
    }

    std::unique_ptr<IO::File> member(manager_->read_member(file.file_name()));
    if (member == nullptr) {
      // The index names a member the archive does not contain.  That is a
      // broken archive, and silently counting it as unchecked would hide it.
      report.mismatches.push_back({file.file_name(), *expected, std::string()});
      continue;
    }

    Util::Sha512 hash;
    std::vector<uint8_t> buffer(64 * 1024);
    for (;;) {
      std::optional<std::size_t> got = member->read(buffer.data(), buffer.size());
      if (!got.has_value() || *got == 0) break;
      hash.update(buffer.data(), *got);
    }

    std::string actual(hash.hex_digest());
    if (actual == *expected) {
      ++report.verified;
    } else {
      report.mismatches.push_back({file.file_name(), *expected, std::move(actual)});
    }
  }

  return report;
}

/******************************************************************************/
std::vector<Schema::File>::const_iterator Index::find(std::string_view name) const
{
  return manager_->find_file(name);
}

/******************************************************************************/
std::vector<Schema::File>::const_iterator
Index::find_file(Schema::EntityType::Type et, Schema::DataKind::Type dk) const
{
  return manager_->find_file(et, dk);
}

/******************************************************************************/
const std::string& Index::version() const { return manager_->version(); }

/******************************************************************************/
const ImsCalibration& Index::ims_calibration() const
{
  return manager_->ims_calibration();
}

/******************************************************************************/
const RunMetadata& Index::metadata() const { return manager_->metadata(); }

/******************************************************************************/
std::shared_ptr<Util::Manager> Index::manager() const { return manager_; }

/******************************************************************************/
bool Index::has_spectra(SpectraSource source) const
{
  const auto kind = source == SpectraSource::Peaks ? Schema::DataKind::Peaks
                                                   : Schema::DataKind::DataArray;
  return manager_->find_file(Schema::EntityType::Spectrum, kind) !=
         manager_->files().end();
}

/******************************************************************************/
Spectra Index::spectra(MetadataDetail detail, SpectraSource source) const
{
  using enum Schema::DataKind::Type;

  if (source == SpectraSource::Peaks && !has_spectra(SpectraSource::Peaks)) {
    throw ParquetError("no spectra_peaks file in this archive");
  }

  const Schema::File* data_file = nullptr;
  const Schema::File* meta_file = nullptr;

  for (const auto& file : manager_->files()) {
    if (file.entity_type().type() != Schema::EntityType::Spectrum) continue;
    if (file.data_kind().type() == DataArray || file.data_kind().type() == Peaks) {
      // Prefer DataArray (profile); fall back to Peaks (centroid-only).  When
      // the caller explicitly asked for Peaks, that file is the primary one.
      const bool wanted = source == SpectraSource::Peaks
                              ? file.data_kind().type() == Peaks
                              : file.data_kind().type() == DataArray;
      if (!data_file || wanted) data_file = &file;
    } else if (file.data_kind().type() == Metadata) {
      meta_file = &file;
    }
  }

  if (!data_file) return Spectra(); // no spectrum data

  const Schema::File* peaks_file = nullptr;
  for (const auto& file : manager_->files()) {
    if (file.entity_type().type() == Schema::EntityType::Spectrum &&
        file.data_kind().type() == Peaks &&
        data_file->data_kind().type() == DataArray) {
      peaks_file = &file;
    }
  }

  std::unique_ptr<Metadata::Table> meta = nullptr;
  if (meta_file) {
    meta = std::make_unique<Metadata::Table>(manager_->parquet(*meta_file));

    // Newer writers put each metadata facet in its own file; attach them so the
    // reader can join them by source_index.  Older writers have none of these.
    for (const auto& file : manager_->files()) {
      if (file.entity_type().type() != Schema::EntityType::Spectrum) continue;
      if (file.data_kind().type() == Scans ||
          file.data_kind().type() == Precursors ||
          file.data_kind().type() == SelectedIons) {
        meta->add_facet(file.data_kind().type().value(), manager_->parquet(file));
      }
    }
  }

  // Read the descriptive metadata through the archive's cache, so a second
  // Spectra over this Index -- the per-thread reader the library's own thread
  // safety rules require -- shares the map instead of re-reading it.  The
  // Table is still handed to the Spectra: Spectrum keeps it for the queries
  // that go past the cached map.
  //
  // Minimal falls back to Lean, shared under Lean, for an archive it cannot
  // hold (Manager caches that verdict too).
  std::shared_ptr<const Spectra::MetadataMap> md;
  std::shared_ptr<const Spectra::MinimalMetadataMap> minimal;
  if (meta) {
    Metadata::Table* table = meta.get();
    if (detail == MetadataDetail::Minimal) {
      minimal = manager_->minimal_spectrum_metadata(
          [table] { return table->read_minimal_spectrum_metadata(); });
    }
    if (!minimal) {
      const MetadataDetail d =
          detail == MetadataDetail::Minimal ? MetadataDetail::Lean : detail;
      md = manager_->spectrum_metadata(
          d, [table, d] { return table->read_spectrum_metadata(d); });
    }
  }

  std::unique_ptr<Data::Signals> data =
      std::make_unique<Data::Signals>(manager_->parquet(*data_file));

  if (peaks_file) {
    std::unique_ptr<Data::Signals> peaks =
        std::make_unique<Data::Signals>(manager_->parquet(*peaks_file));
    return Spectra(std::move(data), std::move(peaks), std::move(meta),
                   manager_->ims_calibration(), std::move(md), std::move(minimal));
  }

  return Spectra(std::move(data), std::move(meta), manager_->ims_calibration(),
                 std::move(md), std::move(minimal));
}

/******************************************************************************/
WavelengthSpectra Index::wavelength_spectra() const
{
  using enum Schema::DataKind::Type;

  const Schema::File* data = nullptr;
  const Schema::File* metadata = nullptr;
  for (const auto& file : manager_->files()) {
    if (file.entity_type().type() != Schema::EntityType::WavelengthSpectrum)
      continue;
    if (file.data_kind().type() == DataArray) data = &file;
    // FIRST match wins.  Several members can carry data_kind "metadata" for one
    // entity -- the specification's own chromatogram example labels the
    // precursor facet that way -- and the primary table is written first.
    // Taking the last would hand a facet to the primary reader, which finds no
    // index column and returns an empty map with no error.
    else if (file.data_kind().type() == Metadata && metadata == nullptr)
      metadata = &file;
  }

  if (!data) return WavelengthSpectra();

  // Read the declared count and the descriptive metadata from the metadata
  // table.  The Parquet objects must outlive the read, so they are held here.
  std::optional<std::size_t> count;
  std::map<uint64_t, WavelengthSpectrumMetadata> md;
  if (metadata) {
    auto md_parquet(manager_->parquet(*metadata));
    if (auto fmd = md_parquet->file_metadata()) {
      if (auto kv = fmd->key_value_metadata()) {
        auto result(kv->Get("wavelength_spectrum_count"));
        if (result.ok()) {
          std::size_t r{};
          const std::string& value(result.ValueOrDie());
          auto [ptr,
                ec]{std::from_chars(value.data(), value.data() + value.size(), r)};
          if (ec == std::errc()) count = r;
        }
      }
    }

    Util::WavelengthMetadataFiles files;
    files.primary = md_parquet.get();

    // A split-layout writer puts the scan facet in its own file.
    std::unique_ptr<Util::Parquet> scans;
    for (const auto& file : manager_->files()) {
      if (file.entity_type().type() != Schema::EntityType::WavelengthSpectrum)
        continue;
      if (file.data_kind().type() != Scans) continue;
      scans = manager_->parquet(file);
      files.scans = scans.get();
    }

    md = Util::read_wavelength_spectrum_metadata(files);
  }

  return WavelengthSpectra(std::make_unique<Data::Signals>(manager_->parquet(*data)),
                           count, std::move(md));
}

/******************************************************************************/
Chromatograms Index::chromatograms() const
{
  using enum Schema::DataKind::Type;

  const Schema::File* data = nullptr;
  const Schema::File* metadata = nullptr;
  for (const auto& file : manager_->files()) {
    if (file.entity_type().type() != Schema::EntityType::Chromatogram) continue;
    if (file.data_kind().type() == DataArray) data = &file;
    // FIRST match wins; see the wavelength case above.
    else if (file.data_kind().type() == Metadata && metadata == nullptr)
      metadata = &file;
  }

  if (!data) return Chromatograms();

  std::optional<std::size_t> count;
  std::map<uint64_t, ChromatogramMetadata> md;
  if (metadata) {
    auto md_parquet(manager_->parquet(*metadata));
    if (auto fmd = md_parquet->file_metadata()) {
      if (auto kv = fmd->key_value_metadata()) {
        auto result(kv->Get("chromatogram_count"));
        if (result.ok()) {
          std::size_t r{};
          const std::string& value(result.ValueOrDie());
          auto [ptr,
                ec]{std::from_chars(value.data(), value.data() + value.size(), r)};
          if (ec == std::errc()) count = r;
        }
      }
    }

    Util::ChromatogramMetadataFiles files;
    files.primary = md_parquet.get();

    std::unique_ptr<Util::Parquet> precursors;
    std::unique_ptr<Util::Parquet> selected_ions;
    for (const auto& file : manager_->files()) {
      if (file.entity_type().type() != Schema::EntityType::Chromatogram) continue;
      if (file.data_kind().type() == Precursors) {
        precursors = manager_->parquet(file);
        files.precursors = precursors.get();
      } else if (file.data_kind().type() == SelectedIons) {
        selected_ions = manager_->parquet(file);
        files.selected_ions = selected_ions.get();
      }
    }

    md = Util::read_chromatogram_metadata(files);
  }

  return Chromatograms(std::make_unique<Data::Signals>(manager_->parquet(*data)),
                       count, std::move(md));
}

} // namespace MzPeak
