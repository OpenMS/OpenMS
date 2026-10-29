/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/writer.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <numeric>
#include <string>
#include <zip.h>
#ifndef ZIP_LENGTH_TO_END // libzip < 1.10.1: a length of -1 (or 0) reads the file to its end as well
#define ZIP_LENGTH_TO_END -1
#endif

#include "mzpeak/exception.h"
#include "mzpeak/util/compat.h"
#include "mzpeak/schema/data_kind.h"
#include "mzpeak/util/json_writer.h"
#include "mzpeak/util/parquet_writer.h"
#include "mzpeak/util/sha512.h"

namespace MzPeak {

namespace fs = std::filesystem;

namespace {

/******************************************************************************/
// The flattened parallel point columns of one table (data or peaks).
struct PointColumns {
  std::vector<uint64_t> spectrum_index;
  std::vector<double> mz;
  std::vector<float> intensity;
};

/******************************************************************************/
// Validate all spectra up front so a bad input cannot leave a partial
// artifact behind.
void validate_one(const SpectrumData& s, std::size_t i, const char* context)
{
  if (s.mz.size() != s.intensity.size()) {
    throw ParquetError(std::string(context) + ": spectrum " + std::to_string(i) +
                       " has mismatched mz/intensity lengths");
  }
  if (std::ranges::any_of(s.mz, [](double m) { return std::isnan(m); })) {
    throw ParquetError(std::string(context) + ": spectrum " + std::to_string(i) +
                       " has a NaN m/z; points are written in ascending order "
                       "and NaN has no place in that order");
  }
}

// Points are written in ascending m/z under a sorting_rank:0 declaration.  A
// NaN makes the sort comparator ill-defined (undefined behaviour) and could
// leave the stored coordinate array non-ascending, contradicting the
// declaration a reader may binary search on -- hence the NaN check above.
void validate(const std::vector<SpectrumData>& spectra, const char* context)
{
  for (std::size_t i = 0; i < spectra.size(); ++i) validate_one(spectra[i], i, context);
}

/******************************************************************************/
// Flatten the spectra whose `centroid` flag matches `want_centroid` into the
// parallel point columns, carrying the GLOBAL spectrum index and emitting each
// spectrum's points in ascending m/z order (so sorting_rank:0 holds).
PointColumns flatten(const std::vector<SpectrumData>& spectra, bool want_centroid)
{
  PointColumns c;
  for (std::size_t i = 0; i < spectra.size(); ++i) {
    const SpectrumData& s = spectra[i];
    if (s.centroid != want_centroid) continue;

    std::vector<std::size_t> order(s.mz.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    // stable: points of equal m/z keep the order the caller gave (OpenMS change; an unstable sort permuted them,
    // which changed what a consumer that merges or deisotopes neighbouring peaks made of the spectrum)
    std::ranges::stable_sort(
        order, [&](std::size_t a, std::size_t b) { return s.mz[a] < s.mz[b]; });

    for (std::size_t k : order) {
      c.spectrum_index.push_back(static_cast<uint64_t>(i));
      c.mz.push_back(s.mz[k]);
      c.intensity.push_back(s.intensity[k]);
    }
  }
  return c;
}

/******************************************************************************/
// File-level key/value metadata for a point table.  `total_spectra` is the
// full spectrum count (drives the reader's record_count, the same on every
// table); `n_points` is the number of points in THIS table.
std::map<std::string, std::string> point_file_kv(std::size_t total_spectra,
                                                 std::size_t n_points)
{
  return {
      {"spectrum_array_index", Util::point_spectra_array_index_json()},
      {"spectrum_count", std::to_string(total_spectra)},
      {"spectrum_data_point_count", std::to_string(n_points)},
  };
}

/******************************************************************************/
bool any_centroid(const std::vector<SpectrumData>& spectra)
{
  return std::ranges::any_of(spectra,
                             [](const SpectrumData& s) { return s.centroid; });
}

bool any_profile(const std::vector<SpectrumData>& spectra)
{
  return std::ranges::any_of(spectra,
                             [](const SpectrumData& s) { return !s.centroid; });
}

/******************************************************************************/
// The three spectrum metadata facets and the column_mapping each declares.
// One definition: both index writers below emit exactly these, and the reader
// resolves MS:1000827 and friends through these mappings when the column name
// alone does not tell it.
std::vector<Util::IndexFileEntry> spectra_facet_entries()
{
  using Schema::DataKind;
  return {
      {"spectra_metadata_scans.parquet",
       "spectrum",
       DataKind(DataKind::Scans).to_string(),
       {{"scan start time", "scan_start_time", "MS:1000016", "UO:0000031"}}},
      {"spectra_metadata_precursors.parquet",
       "spectrum",
       DataKind(DataKind::Precursors).to_string(),
       {{"isolation window target m/z", "isolation_window.isolation_window_target",
         "MS:1000827", "MS:1000040"},
        {"isolation window lower offset",
         "isolation_window.isolation_window_lower_offset", "MS:1000828",
         "MS:1000040"},
        {"isolation window upper offset",
         "isolation_window.isolation_window_upper_offset", "MS:1000829",
         "MS:1000040"}}},
      {"spectra_metadata_selected_ions.parquet",
       "spectrum",
       DataKind(DataKind::SelectedIons).to_string(),
       {{"selected ion m/z", "selected_ion_mz", "MS:1000744", "MS:1000040"},
        {"charge state", "charge_state", "MS:1000041", ""},
        {"intensity", "peak_intensity", "MS:1000042", "MS:1000131"}}},
  };
}

/******************************************************************************/
// mzpeak_index.json describing the data + metadata tables (+ peaks if present).
// WRT-2: when `run_metadata` is non-null, its serialized run-level blocks are
// merged into the emitted `metadata{}` alongside `version`.
std::string
spectra_index_json(bool with_data,
                   bool with_peaks,
                   const RunMetadata* run_metadata = nullptr,
                   const std::map<std::string, std::string>& checksums = {})
{
  using Schema::DataKind;
  std::vector<Util::IndexFileEntry> files;
  if (with_data) {
    files.push_back({"spectra_data.parquet", "spectrum",
                     Schema::DataKind(DataKind::DataArray).to_string()});
  }
  files.push_back(
      {"spectra_metadata.parquet",
       "spectrum",
       Schema::DataKind(DataKind::Metadata).to_string(),
       {{"ms level", "ms_level", "MS:1000511", ""},
        {"scan polarity", "scan_polarity", "MS:1000465", ""},
        {"spectrum representation", "spectrum_representation", "MS:1000525", ""},
        {"number of data points", "number_of_data_points", "MS:1003060", ""},
        {"number of peaks", "number_of_peaks", "MS:1003059", ""}}});
  for (auto& facet : spectra_facet_entries()) files.push_back(std::move(facet));
  if (with_peaks) {
    files.push_back({"spectra_peaks.parquet", "spectrum",
                     Schema::DataKind(DataKind::Peaks).to_string()});
  }
  // The index cannot carry its own digest, so only the members written beside
  // it are stamped; a name absent from the map stays null rather than empty.
  for (auto& file : files) {
    if (auto it = checksums.find(file.name); it != checksums.end()) {
      file.checksum = it->second;
    }
  }

  if (run_metadata != nullptr) {
    return Util::mzpeak_index_json(files, "0.9.0", run_metadata->to_json());
  }
  return Util::mzpeak_index_json(files, "0.9.0");
}

/******************************************************************************/
// Per-spectrum metadata rows.  The data-point / peak counts gate which table
// a reader loads a spectrum from, so they must match where the points went.
Util::SpectrumMetaRow metadata_row(std::size_t i, const SpectrumData& s)
{
  return {/*index=*/static_cast<uint64_t>(i),
          /*id=*/s.id.value_or("index=" + std::to_string(i)),
          /*ms_level=*/s.ms_level,
          /*retention_time=*/s.retention_time,
          /*polarity=*/s.polarity,
          /*number_of_data_points=*/s.centroid ? uint64_t{0} : s.mz.size(),
          /*number_of_peaks=*/s.centroid ? s.mz.size() : uint64_t{0},
          /*representation=*/s.centroid ? "MS:1000127" : "MS:1000128",
          /*precursors=*/s.precursors,
          /*ion_mobility=*/s.ion_mobility,
          /*ion_mobility_type=*/s.ion_mobility_type};
}

std::vector<Util::SpectrumMetaRow>
build_metadata_rows(const std::vector<SpectrumData>& spectra)
{
  std::vector<Util::SpectrumMetaRow> rows;
  rows.reserve(spectra.size());
  for (std::size_t i = 0; i < spectra.size(); ++i) {
    rows.push_back(metadata_row(i, spectra[i]));
  }
  return rows;
}

/******************************************************************************/
// A file to be emitted, held as bytes.  Both the directory and archive writers
// consume the same list, so an archive and an unpacked directory of the same
// run differ only in how the bytes land.
struct Member {
  std::string name;
  std::string bytes;
};

/******************************************************************************/
/// SHA-512 of a file already written to disk, streamed rather than slurped: a
/// point table can be gigabytes, and hashing must not be the thing that decides
/// how much memory writing a run takes.
std::string checksum_of_file(const fs::path& path)
{
  std::ifstream in(path, std::ios::binary);
  if (!in) throw ParquetError("cannot reopen to checksum: " + path.string());

  Util::Sha512 hash;
  char buffer[64 * 1024];
  while (in.read(buffer, sizeof(buffer)) || in.gcount() > 0) {
    hash.update(buffer, static_cast<std::size_t>(in.gcount()));
    if (!in) break;
  }
  if (in.bad()) throw ParquetError("error reading to checksum: " + path.string());
  return hash.hex_digest();
}

/******************************************************************************/
// Chromatogram types that select ions this writer cannot record.
//
// SRM/MRM names a Q1 AND a Q3; SIM names only a Q1.  Neither survives here:
// ChromatogramData carries no precursor or product at all, the product facet
// has no writer anywhere and an unimplemented reference reader path, and a file
// claiming to be one of these would look valid while lacking the very selection
// that defines it.  Both are refused, each for its own reason.
//
// Both CURIE spellings are listed: mzdata emits MS:1000472/MS:1000473 for
// SIM/SRM while its own reader expects MS:1001472/MS:1001473.
bool selects_ions(const std::string& type, const char** what)
{
  if (type == "MS:1001473" || type == "MS:1000473") {
    *what = "a precursor (Q1) and a product (Q3)";
    return true;
  }
  if (type == "MS:1001472" || type == "MS:1000472") {
    *what = "a precursor (Q1)";
    return true;
  }
  return false;
}

/******************************************************************************/
void validate_run(const RunContents& contents)
{
  validate(contents.spectra, "write_run");

  for (std::size_t i = 0; i < contents.chromatograms.size(); ++i) {
    const ChromatogramData& c = contents.chromatograms[i];
    if (c.time.size() != c.intensity.size()) {
      throw ParquetError("write_run: chromatogram " + std::to_string(i) +
                         " has mismatched time/intensity lengths");
    }
    const char* selects = nullptr;
    if (selects_ions(c.chromatogram_type, &selects)) {
      throw ParquetError("write_run: chromatogram " + std::to_string(i) +
                         " has type '" + c.chromatogram_type + "', which names " +
                         selects +
                         " that this writer cannot record; writing it would "
                         "silently drop the selection that defines it");
    }

    // An empty unit becomes "unit":"" in the array index, which violates the
    // schema and the reference's CURIE parser rejects it.  It is also what one
    // of the decode paths hands back when a file declares no unit, so this is
    // reachable by round-tripping rather than only by a careless caller.
    if (c.intensity_unit.empty()) {
      throw ParquetError("write_run: chromatogram " + std::to_string(i) +
                         " has an empty intensity unit; a unit CURIE is "
                         "required in the array index");
    }

    if (std::ranges::any_of(c.time, [](double t) { return std::isnan(t); })) {
      throw ParquetError("write_run: chromatogram " + std::to_string(i) +
                         " has a NaN time; the points are written in ascending "
                         "order and NaN has no place in that order");
    }
  }

  for (std::size_t i = 0; i < contents.wavelength_spectra.size(); ++i) {
    const WavelengthSpectrumData& w = contents.wavelength_spectra[i];
    if (w.wavelength.size() != w.intensity.size()) {
      throw ParquetError("write_run: wavelength spectrum " + std::to_string(i) +
                         " has mismatched wavelength/intensity lengths");
    }
    if (w.intensity_unit.empty()) {
      throw ParquetError("write_run: wavelength spectrum " + std::to_string(i) +
                         " has an empty intensity unit; a unit CURIE is "
                         "required in the array index");
    }
    if (std::ranges::any_of(w.wavelength, [](float x) { return std::isnan(x); })) {
      throw ParquetError("write_run: wavelength spectrum " + std::to_string(i) +
                         " has a NaN wavelength; the points are written in "
                         "ascending order and NaN has no place in that order");
    }
  }
}

/******************************************************************************/
// One intensity unit per file: the coalesced multi-unit layout needs a column
// per unit, which this writer does not emit.  Mixing them silently would label
// absorbance as detector counts.
template <typename Entity>
std::string single_intensity_unit(const std::vector<Entity>& entities,
                                  const char* what)
{
  const std::string* unit = nullptr;
  for (std::size_t i = 0; i < entities.size(); ++i) {
    const std::string& u = entities[i].intensity_unit;
    if (unit == nullptr) {
      unit = &u;
    } else if (*unit != u) {
      throw ParquetError(std::string("write_run: ") + what + " " +
                         std::to_string(i) + " uses intensity unit '" + u +
                         "' but an earlier one uses '" + *unit +
                         "'; this writer emits a single intensity column per "
                         "file and cannot carry both");
    }
  }
  return unit ? *unit : std::string("MS:1000131");
}

/******************************************************************************/
// The order that sorts one entity's points along its primary axis.  Both entity
// types need it and neither cares how it is obtained.
template <typename Axis>
std::vector<std::size_t> ascending_order(const std::vector<Axis>& axis)
{
  std::vector<std::size_t> order(axis.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::ranges::stable_sort(order, // stable: equal coordinates keep their order (OpenMS change)
                           [&](std::size_t a, std::size_t b) { return axis[a] < axis[b]; });
  return order;
}

/******************************************************************************/
// Flatten chromatograms into parallel point columns, converting the API's
// seconds to the minutes the format stores, and emitting each chromatogram's
// points in ascending time order so sorting_rank 0 holds.
struct ChromatogramColumns {
  std::vector<uint64_t> index;
  std::vector<double> time;
  std::vector<float> intensity;
};

ChromatogramColumns flatten_chromatograms(const std::vector<ChromatogramData>& in)
{
  ChromatogramColumns c;
  for (std::size_t i = 0; i < in.size(); ++i) {
    for (std::size_t k : ascending_order(in[i].time)) {
      c.index.push_back(static_cast<uint64_t>(i));
      c.time.push_back(Util::minutes_from_seconds(in[i].time[k])); // seconds -> minutes
      c.intensity.push_back(in[i].intensity[k]);
    }
  }
  return c;
}

/******************************************************************************/
struct WavelengthColumns {
  std::vector<uint64_t> index;
  std::vector<float> wavelength;
  std::vector<float> intensity;
};

WavelengthColumns flatten_wavelength(const std::vector<WavelengthSpectrumData>& in)
{
  WavelengthColumns c;
  for (std::size_t i = 0; i < in.size(); ++i) {
    for (std::size_t k : ascending_order(in[i].wavelength)) {
      c.index.push_back(static_cast<uint64_t>(i));
      c.wavelength.push_back(in[i].wavelength[k]);
      c.intensity.push_back(in[i].intensity[k]);
    }
  }
  return c;
}

/******************************************************************************/
std::vector<Util::ChromatogramMetaRow>
build_chromatogram_rows(const std::vector<ChromatogramData>& in)
{
  std::vector<Util::ChromatogramMetaRow> rows;
  rows.reserve(in.size());
  for (std::size_t i = 0; i < in.size(); ++i) {
    rows.push_back({/*index=*/static_cast<uint64_t>(i),
                    /*id=*/in[i].id.value_or("chromatogram=" + std::to_string(i)),
                    /*chromatogram_type=*/in[i].chromatogram_type,
                    /*polarity=*/in[i].polarity,
                    /*number_of_data_points=*/in[i].time.size()});
  }
  return rows;
}

/******************************************************************************/
// Summary fields are computed from the values actually written.
//
// The maximum is a real maximum, seeded from the first element rather than from
// zero: a baseline-subtracted absorbance spectrum is entirely negative, and a
// zero-seeded maximum records lambda max 0 and base peak 0 for it -- plausible,
// and wrong.  The observed range comes from the sorted output, not the input
// order, for the same reason.
std::vector<Util::WavelengthMetaRow>
build_wavelength_rows(const std::vector<WavelengthSpectrumData>& in)
{
  std::vector<Util::WavelengthMetaRow> rows;
  rows.reserve(in.size());

  for (std::size_t i = 0; i < in.size(); ++i) {
    const WavelengthSpectrumData& w = in[i];

    Util::WavelengthMetaRow row;
    row.index = static_cast<uint64_t>(i);
    row.id = w.id.value_or("wavelength_spectrum=" + std::to_string(i));
    row.time =
        w.time.has_value() ? std::optional<double>(Util::minutes_from_seconds(*w.time)) : std::nullopt;
    row.number_of_data_points = w.wavelength.size();

    if (!w.wavelength.empty()) {
      auto [lo, hi] = std::ranges::minmax_element(w.wavelength);
      row.lowest_observed_wavelength = static_cast<double>(*lo);
      row.highest_observed_wavelength = static_cast<double>(*hi);
    }

    if (!w.intensity.empty()) {
      auto peak = std::ranges::max_element(w.intensity);
      row.base_peak_intensity = *peak;
      row.lambda_max = static_cast<double>(
          w.wavelength[static_cast<std::size_t>(peak - w.intensity.begin())]);
      // Accumulated in double: summing a thousand small samples into a float
      // alongside one large one loses the small ones entirely.
      row.total_ion_current = static_cast<float>(
          std::accumulate(w.intensity.begin(), w.intensity.end(), 0.0));
    }

    rows.push_back(std::move(row));
  }
  return rows;
}

/******************************************************************************/
// Encode every member of a run, and produce the matching mzpeak_index.json.
std::vector<Member> build_run_members(const RunContents& contents,
                                      const RunMetadata* run_metadata)
{
  using Schema::DataKind;
  std::vector<Member> members;
  std::vector<Util::IndexFileEntry> files;

  // ---- mass spectra ------------------------------------------------------
  const std::size_t total_spectra = contents.spectra.size();
  if (!contents.spectra.empty()) {
    if (any_profile(contents.spectra)) {
      PointColumns data(flatten(contents.spectra, /*want_centroid=*/false));
      members.push_back({"spectra_data.parquet",
                         Util::point_spectra_data_bytes(
                             data.spectrum_index, data.mz, data.intensity,
                             point_file_kv(total_spectra, data.mz.size()))});
      files.push_back({"spectra_data.parquet", "spectrum",
                       Schema::DataKind(DataKind::DataArray).to_string()});
    }
    if (any_centroid(contents.spectra)) {
      PointColumns peaks(flatten(contents.spectra, /*want_centroid=*/true));
      members.push_back({"spectra_peaks.parquet",
                         Util::point_spectra_data_bytes(
                             peaks.spectrum_index, peaks.mz, peaks.intensity,
                             point_file_kv(total_spectra, peaks.mz.size()))});
      files.push_back({"spectra_peaks.parquet", "spectrum",
                       Schema::DataKind(DataKind::Peaks).to_string()});
    }

    auto rows(build_metadata_rows(contents.spectra));
    members.push_back(
        {"spectra_metadata.parquet", Util::spectra_metadata_bytes(rows)});
    files.push_back(
        {"spectra_metadata.parquet",
         "spectrum",
         Schema::DataKind(DataKind::Metadata).to_string(),
         {{"ms level", "ms_level", "MS:1000511", ""},
          {"scan polarity", "scan_polarity", "MS:1000465", ""},
          {"spectrum representation", "spectrum_representation", "MS:1000525", ""},
          {"number of data points", "number_of_data_points", "MS:1003060", ""},
          {"number of peaks", "number_of_peaks", "MS:1003059", ""}}});

    std::array<std::string, 3> facets(Util::spectra_metadata_facet_bytes(rows));
    members.push_back({"spectra_metadata_scans.parquet", std::move(facets[0])});
    members.push_back({"spectra_metadata_precursors.parquet", std::move(facets[1])});
    members.push_back(
        {"spectra_metadata_selected_ions.parquet", std::move(facets[2])});
    for (auto& facet : spectra_facet_entries()) files.push_back(std::move(facet));
  }

  // ---- chromatograms -----------------------------------------------------
  if (!contents.chromatograms.empty()) {
    const std::string unit(
        single_intensity_unit(contents.chromatograms, "chromatogram"));
    ChromatogramColumns c(flatten_chromatograms(contents.chromatograms));

    members.push_back(
        {"chromatograms_data.parquet",
         Util::point_chromatograms_data_bytes(
             c.index, c.time, c.intensity,
             {{"chromatogram_array_index",
               Util::point_chromatograms_array_index_json(unit)},
              {"chromatogram_count", std::to_string(contents.chromatograms.size())},
              {"chromatogram_data_point_count", std::to_string(c.time.size())}})});
    files.push_back({"chromatograms_data.parquet", "chromatogram",
                     Schema::DataKind(DataKind::DataArray).to_string()});

    members.push_back({"chromatograms_metadata.parquet",
                       Util::chromatograms_metadata_bytes(
                           build_chromatogram_rows(contents.chromatograms),
                           {{"chromatogram_count",
                             std::to_string(contents.chromatograms.size())}})});
    files.push_back(
        {"chromatograms_metadata.parquet",
         "chromatogram",
         Schema::DataKind(DataKind::Metadata).to_string(),
         {{"chromatogram type", "chromatogram_type", "MS:1000626", ""},
          {"scan polarity", "scan_polarity", "MS:1000465", ""},
          {"number of data points", "number_of_data_points", "MS:1003060", ""}}});

    // The current reference reader opens both facets unconditionally when it
    // loads chromatogram metadata, so a plain TIC without them fails to load
    // there.  Emitted with their schema and zero rows.
    std::array<std::string, 2> facets(Util::chromatogram_facet_bytes(
        {{"chromatogram_count", std::to_string(contents.chromatograms.size())}}));
    members.push_back(
        {"chromatograms_metadata_precursors.parquet", std::move(facets[0])});
    members.push_back(
        {"chromatograms_metadata_selected_ions.parquet", std::move(facets[1])});
    files.push_back({"chromatograms_metadata_precursors.parquet", "chromatogram",
                     Schema::DataKind(DataKind::Precursors).to_string()});
    files.push_back({"chromatograms_metadata_selected_ions.parquet", "chromatogram",
                     Schema::DataKind(DataKind::SelectedIons).to_string()});
  }

  // ---- wavelength spectra ------------------------------------------------
  if (!contents.wavelength_spectra.empty()) {
    const std::string unit(
        single_intensity_unit(contents.wavelength_spectra, "wavelength spectrum"));
    WavelengthColumns w(flatten_wavelength(contents.wavelength_spectra));
    const std::size_t count = contents.wavelength_spectra.size();

    members.push_back({"wavelength_spectra_data.parquet",
                       Util::point_wavelength_data_bytes(
                           w.index, w.wavelength, w.intensity,
                           {{"wavelength_spectrum_array_index",
                             Util::point_wavelength_array_index_json(unit)},
                            {"wavelength_spectrum_count", std::to_string(count)},
                            {"wavelength_spectrum_data_point_count",
                             std::to_string(w.wavelength.size())}})});
    files.push_back({"wavelength_spectra_data.parquet", "wavelength_spectrum",
                     Schema::DataKind(DataKind::DataArray).to_string()});

    members.push_back({"wavelength_spectra_metadata.parquet",
                       Util::wavelength_metadata_bytes(
                           build_wavelength_rows(contents.wavelength_spectra),
                           {{"wavelength_spectrum_count", std::to_string(count)}})});
    files.push_back(
        {"wavelength_spectra_metadata.parquet",
         "wavelength_spectrum",
         Schema::DataKind(DataKind::Metadata).to_string(),
         {{"number of data points", "number_of_data_points", "MS:1003060", ""},
          {"lowest observed wavelength", "lowest_observed_wavelength", "MS:1000619",
           "UO:0000018"},
          {"highest observed wavelength", "highest_observed_wavelength",
           "MS:1000618", "UO:0000018"},
          {"lambda max", "lambda_max", "MS:1003812", "UO:0000018"},
          {"spectrum type", "spectrum_type", "MS:1000559", ""},
          {"spectrum representation", "spectrum_representation", "MS:1000525", ""},
          // These summarise the intensity array, so they carry ITS unit.  Left
          // null, an absorbance base peak is indistinguishable from a detector
          // count to anything reading the file.
          {"base peak intensity", "base_peak_intensity", "MS:1000505", unit},
          {"total ion current", "total_ion_current", "MS:1000285", unit}}});

    // The reference reader ignores the primary `time` column outright, so
    // without a scan facet the acquisition time is invisible to it.
    members.push_back({"wavelength_spectra_metadata_scans.parquet",
                       Util::wavelength_scans_bytes(
                           build_wavelength_rows(contents.wavelength_spectra),
                           {{"wavelength_spectrum_count", std::to_string(count)}})});
    files.push_back(
        {"wavelength_spectra_metadata_scans.parquet",
         "wavelength_spectrum",
         Schema::DataKind(DataKind::Scans).to_string(),
         {{"scan start time", "scan_start_time", "MS:1000016", "UO:0000031"}}});
  }

  // Run-level metadata goes in unconditionally.  The reference writer copies it
  // only in its chromatogram close path, so a wavelength-only archive can lose
  // the version, CV list and run description entirely.
  // Stamp each index entry with the digest of the member it names, which must
  // happen before the index is serialised and cannot cover the index itself.
  for (auto& file : files) {
    const auto member = std::ranges::find(members, file.name, &Member::name);
    if (member != members.end()) file.checksum = Util::sha512_hex(member->bytes);
  }

  std::string index_json(
      run_metadata != nullptr
          ? Util::mzpeak_index_json(files, "0.9.0", run_metadata->to_json())
          : Util::mzpeak_index_json(files, "0.9.0"));
  members.push_back({"mzpeak_index.json", std::move(index_json)});

  return members;
}

/******************************************************************************/
// Add a STORED (uncompressed) member to a zip archive from an in-memory
// buffer.  The buffer must outlive zip_close (libzip reads the source
// lazily), so callers keep the backing std::string in scope.
void add_stored_member(zip_t* archive, const char* name, const std::string& data)
{
  zip_source_t* source =
      zip_source_buffer(archive, data.data(), data.size(), /*freep=*/0);
  if (source == nullptr) {
    throw ParquetError(std::string("write_spectra_archive: zip_source_buffer "
                                   "failed for ") +
                       name + ": " + zip_strerror(archive));
  }

  zip_int64_t idx =
      zip_file_add(archive, name, source, ZIP_FL_ENC_UTF_8 | ZIP_FL_OVERWRITE);
  if (idx < 0) {
    // On failure libzip does not take ownership of the source.
    zip_source_free(source);
    throw ParquetError(std::string("write_spectra_archive: zip_file_add failed "
                                   "for ") +
                       name + ": " + zip_strerror(archive));
  }

  // The mzPeak spec mandates ZIP_CM_STORE; the reader rejects compressed
  // members.
  if (zip_set_file_compression(archive, static_cast<zip_uint64_t>(idx), ZIP_CM_STORE,
                               0) != 0) {
    throw ParquetError(std::string("write_spectra_archive: "
                                   "zip_set_file_compression failed for ") +
                       name + ": " + zip_strerror(archive));
  }
}

/******************************************************************************/
// Add a STORED member backed by a file on disk.  libzip reads it at
// zip_close, so the file must stay put until then.
void add_stored_file_member(zip_t* archive, const char* name, const fs::path& path)
{
  zip_source_t* source =
      zip_source_file(archive, Util::narrow(path).c_str(), 0, ZIP_LENGTH_TO_END);
  if (source == nullptr) {
    throw ParquetError(std::string("RunArchiveWriter: zip_source_file failed for ") +
                       name + ": " + zip_strerror(archive));
  }
  zip_int64_t idx =
      zip_file_add(archive, name, source, ZIP_FL_ENC_UTF_8 | ZIP_FL_OVERWRITE);
  if (idx < 0) {
    zip_source_free(source);
    throw ParquetError(std::string("RunArchiveWriter: zip_file_add failed for ") +
                       name + ": " + zip_strerror(archive));
  }
  if (zip_set_file_compression(archive, static_cast<zip_uint64_t>(idx), ZIP_CM_STORE,
                               0) != 0) {
    throw ParquetError(std::string("RunArchiveWriter: zip_set_file_compression "
                                   "failed for ") +
                       name + ": " + zip_strerror(archive));
  }
}

/******************************************************************************/
void write_index_file(const fs::path& path, const std::string& json)
{
  std::ofstream out(path, std::ios::binary);
  if (!out) {
    throw ParquetError("cannot open index for writing: " + path.string());
  }
  out << json;
  out.close();
  if (!out) {
    throw ParquetError("error writing index: " + path.string());
  }
}

/******************************************************************************/
// WRT-2: shared directory writer; `run_metadata` is null for the no-metadata
// overload and otherwise carries the run-level blocks to emit.
void write_spectra_directory_impl(const fs::path& dir,
                                  const std::vector<SpectrumData>& spectra,
                                  const RunMetadata* run_metadata)
{
  validate(spectra, "write_spectra_directory");

  fs::path tmp_dir = dir;
  tmp_dir += ".tmp";

  try {
    fs::create_directories(tmp_dir);
  } catch (const fs::filesystem_error& e) {
    throw ParquetError(std::string("write_spectra_directory: ") + e.what());
  }

  try {
    const std::size_t total = spectra.size();
    const bool with_data = any_profile(spectra);
    const bool with_peaks = any_centroid(spectra);

    if (with_data) {
      PointColumns data(flatten(spectra, /*want_centroid=*/false));
      Util::write_point_spectra_data((tmp_dir / "spectra_data.parquet").string(),
                                     data.spectrum_index, data.mz, data.intensity,
                                     point_file_kv(total, data.mz.size()));
    }

    if (with_peaks) {
      PointColumns peaks(flatten(spectra, /*want_centroid=*/true));
      Util::write_point_spectra_data((tmp_dir / "spectra_peaks.parquet").string(),
                                     peaks.spectrum_index, peaks.mz, peaks.intensity,
                                     point_file_kv(total, peaks.mz.size()));
    }

    Util::write_spectra_metadata((tmp_dir / "spectra_metadata.parquet").string(),
                                 build_metadata_rows(spectra));
    Util::write_spectra_metadata_facets(tmp_dir.string(),
                                        build_metadata_rows(spectra));

    // Hashed back off disk rather than kept in memory: these tables were
    // written straight through, and a run large enough to matter is a run too
    // large to hold twice.
    std::map<std::string, std::string> checksums;
    for (const auto& entry : fs::directory_iterator(tmp_dir)) {
      if (!entry.is_regular_file()) continue;
      checksums.emplace(entry.path().filename().string(),
                        checksum_of_file(entry.path()));
    }

    write_index_file(
        tmp_dir / "mzpeak_index.json",
        spectra_index_json(with_data, with_peaks, run_metadata, checksums));

    // Remove any existing target before the atomic rename (M10: POSIX rename
    // returns ENOTEMPTY for non-empty directories).
    if (fs::exists(dir)) fs::remove_all(dir);
    fs::rename(tmp_dir, dir);
  } catch (...) {
    if (fs::exists(tmp_dir)) {
      fs::remove_all(tmp_dir);
    }
    throw;
  }
}

/******************************************************************************/
// WRT-2: shared archive writer; `run_metadata` is null for the no-metadata
// overload and otherwise carries the run-level blocks to emit.
void write_spectra_archive_impl(const fs::path& zip_path,
                                const std::vector<SpectrumData>& spectra,
                                const RunMetadata* run_metadata)
{
  validate(spectra, "write_spectra_archive");

  const std::size_t total = spectra.size();
  const bool with_data = any_profile(spectra);
  const bool with_peaks = any_centroid(spectra);

  // Encode all members in memory first.  These strings back the libzip
  // sources and MUST stay alive until zip_close returns.
  std::string data_bytes;
  if (with_data) {
    PointColumns data(flatten(spectra, /*want_centroid=*/false));
    data_bytes =
        Util::point_spectra_data_bytes(data.spectrum_index, data.mz, data.intensity,
                                       point_file_kv(total, data.mz.size()));
  }

  std::string peaks_bytes;
  if (with_peaks) {
    PointColumns peaks(flatten(spectra, /*want_centroid=*/true));
    peaks_bytes = Util::point_spectra_data_bytes(
        peaks.spectrum_index, peaks.mz, peaks.intensity,
        point_file_kv(total, peaks.mz.size()));
  }

  std::string metadata_bytes(
      Util::spectra_metadata_bytes(build_metadata_rows(spectra)));
  std::array<std::string, 3> facet_bytes(
      Util::spectra_metadata_facet_bytes(build_metadata_rows(spectra)));
  std::map<std::string, std::string> checksums;
  if (with_data) {
    checksums.emplace("spectra_data.parquet", Util::sha512_hex(data_bytes));
  }
  if (with_peaks) {
    checksums.emplace("spectra_peaks.parquet", Util::sha512_hex(peaks_bytes));
  }
  checksums.emplace("spectra_metadata.parquet", Util::sha512_hex(metadata_bytes));
  checksums.emplace("spectra_metadata_scans.parquet",
                    Util::sha512_hex(facet_bytes[0]));
  checksums.emplace("spectra_metadata_precursors.parquet",
                    Util::sha512_hex(facet_bytes[1]));
  checksums.emplace("spectra_metadata_selected_ions.parquet",
                    Util::sha512_hex(facet_bytes[2]));

  std::string index_json(
      spectra_index_json(with_data, with_peaks, run_metadata, checksums));

  int errnum = 0;
  zip_t* archive = zip_open(Util::narrow(zip_path).c_str(), ZIP_CREATE | ZIP_TRUNCATE, &errnum);
  if (archive == nullptr) {
    zip_error_t error;
    zip_error_init_with_code(&error, errnum);
    std::string msg("write_spectra_archive: failed to create archive " +
                    zip_path.string() + ": " + zip_error_strerror(&error));
    zip_error_fini(&error);
    throw ParquetError(msg);
  }

  // Add the members as STORED.  On any error, discard the (unwritten) archive
  // so no partial file is left behind, then rethrow.
  try {
    if (with_data) {
      add_stored_member(archive, "spectra_data.parquet", data_bytes);
    }
    if (with_peaks) {
      add_stored_member(archive, "spectra_peaks.parquet", peaks_bytes);
    }
    add_stored_member(archive, "spectra_metadata.parquet", metadata_bytes);
    add_stored_member(archive, "spectra_metadata_scans.parquet", facet_bytes[0]);
    add_stored_member(archive, "spectra_metadata_precursors.parquet",
                      facet_bytes[1]);
    add_stored_member(archive, "spectra_metadata_selected_ions.parquet",
                      facet_bytes[2]);
    add_stored_member(archive, "mzpeak_index.json", index_json);
  } catch (...) {
    zip_discard(archive);
    throw;
  }

  // zip_close flushes the archive and reads the source buffers; the backing
  // strings above are still in scope here.
  if (zip_close(archive) != 0) {
    std::string msg("write_spectra_archive: failed to finalize archive " +
                    zip_path.string() + ": " + zip_strerror(archive));
    zip_discard(archive);
    throw ParquetError(msg);
  }
}

} // namespace

/******************************************************************************/
void write_spectra_directory(const fs::path& dir,
                             const std::vector<SpectrumData>& spectra)
{
  write_spectra_directory_impl(dir, spectra, /*run_metadata=*/nullptr);
}

/******************************************************************************/
void write_spectra_directory(const fs::path& dir,
                             const std::vector<SpectrumData>& spectra,
                             const RunMetadata& metadata)
{
  write_spectra_directory_impl(dir, spectra, &metadata);
}

/******************************************************************************/
void write_spectra_archive(const fs::path& zip_path,
                           const std::vector<SpectrumData>& spectra)
{
  write_spectra_archive_impl(zip_path, spectra, /*run_metadata=*/nullptr);
}

/******************************************************************************/
void write_spectra_archive(const fs::path& zip_path,
                           const std::vector<SpectrumData>& spectra,
                           const RunMetadata& metadata)
{
  write_spectra_archive_impl(zip_path, spectra, &metadata);
}

/******************************************************************************/
void write_run_directory(const fs::path& dir,
                         const RunContents& contents,
                         const RunMetadata* run_metadata)
{
  validate_run(contents);

  fs::path tmp_dir = dir;
  tmp_dir += ".tmp";

  try {
    fs::create_directories(tmp_dir);
  } catch (const fs::filesystem_error& e) {
    throw ParquetError(std::string("write_run_directory: ") + e.what());
  }

  try {
    for (const Member& member : build_run_members(contents, run_metadata)) {
      std::ofstream out(tmp_dir / member.name, std::ios::binary);
      if (!out) {
        throw ParquetError("cannot open for writing: " +
                           (tmp_dir / member.name).string());
      }
      out.write(member.bytes.data(),
                static_cast<std::streamsize>(member.bytes.size()));
      out.close();
      if (!out) {
        throw ParquetError("error writing: " + (tmp_dir / member.name).string());
      }
    }

    // Remove any existing target before the atomic rename (POSIX rename
    // returns ENOTEMPTY for non-empty directories).
    if (fs::exists(dir)) fs::remove_all(dir);
    fs::rename(tmp_dir, dir);
  } catch (...) {
    if (fs::exists(tmp_dir)) fs::remove_all(tmp_dir);
    throw;
  }
}

/******************************************************************************/
void write_run_archive(const fs::path& zip_path,
                       const RunContents& contents,
                       const RunMetadata* run_metadata)
{
  validate_run(contents);

  // The member bytes back the libzip sources and MUST stay alive until
  // zip_close returns.
  std::vector<Member> members(build_run_members(contents, run_metadata));

  int errnum = 0;
  zip_t* archive = zip_open(Util::narrow(zip_path).c_str(), ZIP_CREATE | ZIP_TRUNCATE, &errnum);
  if (archive == nullptr) {
    zip_error_t error;
    zip_error_init_with_code(&error, errnum);
    std::string msg("write_run_archive: failed to create archive " +
                    zip_path.string() + ": " + zip_error_strerror(&error));
    zip_error_fini(&error);
    throw ParquetError(msg);
  }

  try {
    for (const Member& member : members) {
      add_stored_member(archive, member.name.c_str(), member.bytes);
    }
  } catch (...) {
    zip_discard(archive);
    throw;
  }

  if (zip_close(archive) != 0) {
    std::string msg("write_run_archive: failed to finalize archive " +
                    zip_path.string() + ": " + zip_strerror(archive));
    zip_discard(archive);
    throw ParquetError(msg);
  }
}

/******************************************************************************/
struct RunArchiveWriter::Impl {
  /// One point table (data or peaks): the pending row group and, once the
  /// first group has been flushed, the stream writing to the working file.
  struct Stream {
    fs::path path;
    PointColumns pending;
    std::unique_ptr<Util::PointTableStream> table;
    std::size_t points = 0;
  };

  fs::path zip_path;
  fs::path tmp_dir;
  std::size_t points_per_group;
  std::optional<RunMetadata> metadata;
  Stream data;
  Stream peaks;
  std::vector<Util::SpectrumMetaRow> rows;
  bool finished = false;

  Impl(fs::path path, std::size_t per_group)
      : zip_path(std::move(path)), points_per_group(std::max<std::size_t>(1, per_group))
  {
    tmp_dir = zip_path;
    tmp_dir += ".tmp";
    std::error_code ec;
    fs::remove_all(tmp_dir, ec);
    fs::create_directories(tmp_dir);
    data.path = tmp_dir / "spectra_data.parquet";
    peaks.path = tmp_dir / "spectra_peaks.parquet";
  }

  ~Impl()
  {
    // Finished: the working files are already gone.  Abandoned: take them
    // with us, and the streams' handles before them.
    data.table.reset();
    peaks.table.reset();
    std::error_code ec;
    fs::remove_all(tmp_dir, ec);
  }

  void flush(Stream& s)
  {
    if (s.pending.mz.empty()) return;
    if (!s.table) s.table = std::make_unique<Util::PointTableStream>(s.path.string());
    s.table->write_row_group(s.pending.spectrum_index, s.pending.mz, s.pending.intensity);
    s.points += s.pending.mz.size();
    s.pending.spectrum_index.clear();
    s.pending.mz.clear();
    s.pending.intensity.clear();
  }
};

/******************************************************************************/
RunArchiveWriter::RunArchiveWriter(const fs::path& zip_path,
                                   std::size_t points_per_row_group)
    : impl_(std::make_unique<Impl>(zip_path, points_per_row_group))
{
}

RunArchiveWriter::~RunArchiveWriter() = default;

/******************************************************************************/
void RunArchiveWriter::add(SpectrumData s)
{
  if (impl_->finished) throw ParquetError("RunArchiveWriter: add after finish");
  const std::size_t index = impl_->rows.size();
  validate_one(s, index, "RunArchiveWriter");
  impl_->rows.push_back(metadata_row(index, s));

  // Points in ascending m/z, as flatten() does, so sorting_rank:0 holds.
  Impl::Stream& stream = s.centroid ? impl_->peaks : impl_->data;
  std::vector<std::size_t> order(s.mz.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  // stable: points of equal m/z keep the caller's order (OpenMS change, see flatten())
  std::ranges::stable_sort(order, [&](std::size_t a, std::size_t b) { return s.mz[a] < s.mz[b]; });
  PointColumns& p = stream.pending;
  for (std::size_t k : order) {
    p.spectrum_index.push_back(static_cast<uint64_t>(index));
    p.mz.push_back(s.mz[k]);
    p.intensity.push_back(s.intensity[k]);
  }
  // Whole spectra per group, so a group boundary never splits one: the
  // reader's straddling case is allowed but costs it a second group.
  if (p.mz.size() >= impl_->points_per_group) impl_->flush(stream);
}

/******************************************************************************/
void RunArchiveWriter::set_metadata(const RunMetadata& metadata)
{
  impl_->metadata = metadata;
}

std::size_t RunArchiveWriter::size() const { return impl_->rows.size(); }

/******************************************************************************/
void RunArchiveWriter::finish()
{
  if (impl_->finished) return;
  Impl& im = *impl_;
  im.flush(im.data);
  im.flush(im.peaks);

  const std::size_t total = im.rows.size();
  // A run with spectra but no points at all gets an empty peaks table: the reader locates the spectra through a
  // signal table and would report none without one (OpenMS change).
  if (!im.rows.empty() && !im.data.table && !im.peaks.table) {
    im.peaks.table = std::make_unique<Util::PointTableStream>(im.peaks.path.string());
  }
  const bool with_data = im.data.table != nullptr;
  const bool with_peaks = im.peaks.table != nullptr;
  if (with_data) im.data.table->close(point_file_kv(total, im.data.points));
  if (with_peaks) im.peaks.table->close(point_file_kv(total, im.peaks.points));

  // Same members as write_run_archive(): the point tables from their working
  // files, the metadata tables and the index from memory (small: a few
  // hundred bytes per spectrum).
  std::vector<Member> members;
  if (!im.rows.empty()) {
    members.push_back({"spectra_metadata.parquet", Util::spectra_metadata_bytes(im.rows)});
    std::array<std::string, 3> facets(Util::spectra_metadata_facet_bytes(im.rows));
    members.push_back({"spectra_metadata_scans.parquet", std::move(facets[0])});
    members.push_back({"spectra_metadata_precursors.parquet", std::move(facets[1])});
    members.push_back({"spectra_metadata_selected_ions.parquet", std::move(facets[2])});
  }
  std::map<std::string, std::string> checksums;
  for (const Member& member : members) {
    checksums.emplace(member.name, Util::sha512_hex(member.bytes));
  }
  // The point tables never passed through memory here; they are streamed from
  // the working files they were written to.
  if (with_data) {
    checksums.emplace("spectra_data.parquet", checksum_of_file(im.data.path));
  }
  if (with_peaks) {
    checksums.emplace("spectra_peaks.parquet", checksum_of_file(im.peaks.path));
  }

  std::string index_json;
  if (im.rows.empty()) {
    index_json = Util::mzpeak_index_json({}, "0.9.0");
  } else {
    index_json = spectra_index_json(
        with_data, with_peaks, im.metadata ? &*im.metadata : nullptr, checksums);
  }
  members.push_back({"mzpeak_index.json", std::move(index_json)});

  int errnum = 0;
  zip_t* archive = zip_open(Util::narrow(im.zip_path).c_str(), ZIP_CREATE | ZIP_TRUNCATE, &errnum);
  if (archive == nullptr) {
    zip_error_t error;
    zip_error_init_with_code(&error, errnum);
    std::string msg("RunArchiveWriter: failed to create archive " + im.zip_path.string() +
                    ": " + zip_error_strerror(&error));
    zip_error_fini(&error);
    throw ParquetError(msg);
  }
  try {
    if (with_data) add_stored_file_member(archive, "spectra_data.parquet", im.data.path);
    if (with_peaks) add_stored_file_member(archive, "spectra_peaks.parquet", im.peaks.path);
    for (const Member& member : members) {
      add_stored_member(archive, member.name.c_str(), member.bytes);
    }
  } catch (...) {
    zip_discard(archive);
    throw;
  }
  if (zip_close(archive) != 0) {
    std::string msg("RunArchiveWriter: failed to finalize archive " + im.zip_path.string() +
                    ": " + zip_strerror(archive));
    zip_discard(archive);
    throw ParquetError(msg);
  }
  im.finished = true;
  im.data.table.reset();
  im.peaks.table.reset();
  std::error_code ec;
  fs::remove_all(im.tmp_dir, ec);
}

} // namespace MzPeak
