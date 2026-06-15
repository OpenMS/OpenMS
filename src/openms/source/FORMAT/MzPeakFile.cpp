// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/MzPeakFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/FORMAT/ZipRandomAccessFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/SYSTEM/File.h>
#include <algorithm>
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <parquet/arrow/reader.h>
#include <string>
#include <vector>

namespace OpenMS
{

namespace
{
  // ==========================================================================
  // Null-marked m/z reconstruction (RDR-5).
  //
  // Ported near-verbatim from the mzpeak reference (src/util/null_fill.cpp,
  // PORT-01 verified C++20-clean). Profile m/z columns carry interior NULLs
  // where flanking zero-intensity points were removed by null-marking; each
  // maximal run of valid m/z values emits one reconstructed value into the
  // adjacent null slots, offset by a local delta (the run's median spacing,
  // or the per-spectrum mz_delta_model polynomial for a single-value run).
  // ==========================================================================

  /// Median of an ALREADY-SORTED vector, matching the reference @c median:
  /// <=2 elements -> the first; odd -> the middle; even -> the mean of the two
  /// upper-middle elements (mid and mid+1, where mid = n/2).
  double median_sorted_(const std::vector<double>& d)
  {
    std::size_t n = d.size();
    if (n == 0) return 0.0;
    if (n <= 2) return d[0];
    std::size_t mid = n / 2;
    if (n % 2 == 1) return d[mid];
    return (d[mid] + d[mid + 1]) / 2.0;
  }

  /// Evaluate the m/z delta model @p beta at @p mz: a polynomial
  /// beta[0] + beta[1]*mz + beta[2]*mz^2 + .... A single coefficient is a
  /// constant model. Mirrors the reference @c MZDeltaModel::predict.
  double predict_delta_(const std::vector<double>& beta, double mz)
  {
    double acc = 0.0;
    for (std::size_t i = 0; i < beta.size(); ++i)
    {
      acc += (i == 0) ? beta[i] : beta[i] * std::pow(mz, static_cast<int>(i));
    }
    return acc;
  }

  /// The "median-below-median" delta of the values [begin, end): the median of
  /// the consecutive differences that are <= the overall median difference.
  /// Mirrors the reference @c MedianDeltaEstimator::estimate_median_delta.
  double estimate_median_delta_(const std::vector<double>& values, std::size_t begin, std::size_t end)
  {
    std::vector<double> deltas;
    if (end > begin + 1) deltas.reserve(end - begin - 1);
    for (std::size_t k = begin + 1; k < end; ++k)
    {
      deltas.push_back(values[k] - values[k - 1]);
    }
    std::sort(deltas.begin(), deltas.end());
    if (deltas.empty()) return 0.0;

    double m = median_sorted_(deltas);
    // Keep only deltas <= the median; the vector stays sorted.
    std::erase_if(deltas, [m](double v) { return v > m; });
    if (deltas.empty()) return m;
    return median_sorted_(deltas);
  }

  /// Reconstruct the null-marked m/z values of a profile spectrum.
  ///
  /// @p values[i] is meaningful only where @p valid[i] is true; null (invalid)
  /// positions are flanking zero-intensity points removed by null-marking.
  /// Returns the full reconstructed array, or an EMPTY vector if the layout is
  /// not the expected interior-paired-null form (the caller then leaves nulls
  /// as 0). Mirrors the reference @c fill_nulls_for.
  std::vector<double> reconstruct_null_mz_(const std::vector<double>& values, const std::vector<bool>& valid, const std::vector<double>& beta)
  {
    const std::size_t n = values.size();
    // valid is the per-position null bitmap for values; mismatched sizes would
    // read out of bounds. Decline rather than risk UB (callers treat {} as
    // "no reconstruction").
    if (valid.size() != n) return {};

    std::vector<double> out;
    out.reserve(n);

    std::size_t i = 0;
    while (i < n)
    {
      if (! valid[i])
      {
        ++i; // null slots are filled by the adjacent runs
        continue;
      }

      // Maximal run of valid values [s, e).
      std::size_t s = i;
      std::size_t e = s;
      while (e < n && valid[e])
        ++e;

      std::size_t len = e - s;
      double delta = (len > 1) ? estimate_median_delta_(values, s, e) : predict_delta_(beta, values[s]);

      // Runs are separated by nulls, so a run starting after position 0 has a
      // null immediately before it, and a run ending before n has one after.
      if (s > 0) out.push_back(values[s] - delta);
      for (std::size_t k = s; k < e; ++k)
        out.push_back(values[k]);
      if (e < n) out.push_back(values[e - 1] + delta);

      i = e;
    }

    // Only the interior-paired-null layout reconstructs to the original length.
    // Anything else (e.g. unpaired leading/trailing nulls) is left to the caller.
    if (out.size() != n) return {};
    return out;
  }

  // ==========================================================================
  // Arrow / Parquet helpers.
  // ==========================================================================

  /// Per-spectrum metadata needed to build an MSSpectrum.
  struct SpectrumMeta
  {
    int ms_level = 0;
    double retention_time = 0.0; ///< seconds (mzpeak stores RT in seconds)
    std::vector<double> mz_delta_model;
  };

  /// Read a whole Parquet table from a named archive entry via Arrow.
  std::shared_ptr<arrow::Table> readTableFromArchive_(const String& archive, const String& entry, std::unique_ptr<File::TempDir>& temp_dir)
  {
    auto raf_result = ZipRandomAccessFile::Open(archive, entry, temp_dir);
    if (! raf_result.ok())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, archive,
                                  "Failed to open mzPeak archive entry '" + entry + "': " + raf_result.status().ToString());
    }
    std::shared_ptr<arrow::io::RandomAccessFile> raf = raf_result.ValueOrDie();

    auto reader_result = parquet::arrow::OpenFile(raf, arrow::default_memory_pool());
    if (! reader_result.ok())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, archive,
                                  "Failed to open Parquet reader for '" + entry + "': " + reader_result.status().ToString());
    }
    std::unique_ptr<parquet::arrow::FileReader> reader = std::move(reader_result.ValueOrDie());

    std::shared_ptr<arrow::Table> table;
    auto status = reader->ReadTable(&table);
    if (! status.ok())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, archive,
                                  "Failed to read Parquet table '" + entry + "': " + status.ToString());
    }
    return table;
  }

  /// Read raw bytes of a (small) archive entry into a string (for JSON members).
  std::string readEntryBytes_(const String& archive, const String& entry, std::unique_ptr<File::TempDir>& temp_dir)
  {
    auto raf_result = ZipRandomAccessFile::Open(archive, entry, temp_dir);
    if (! raf_result.ok())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, archive,
                                  "Failed to open mzPeak archive entry '" + entry + "': " + raf_result.status().ToString());
    }
    std::shared_ptr<arrow::io::RandomAccessFile> raf = raf_result.ValueOrDie();

    auto size_result = raf->GetSize();
    if (! size_result.ok())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, archive,
                                  "Failed to size archive entry '" + entry + "': " + size_result.status().ToString());
    }
    int64_t size = size_result.ValueOrDie();

    auto buf_result = raf->ReadAt(0, size);
    if (! buf_result.ok())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, archive,
                                  "Failed to read archive entry '" + entry + "': " + buf_result.status().ToString());
    }
    std::shared_ptr<arrow::Buffer> buf = buf_result.ValueOrDie();
    return std::string(reinterpret_cast<const char*>(buf->data()), static_cast<std::size_t>(buf->size()));
  }

  /// Resolve the top-level "point" struct column to a StructArray. The mzpeak
  /// point layout stores {spectrum_index uint64, mz double, intensity float}
  /// under a single struct column. Returns nullptr if absent or not a struct.
  std::shared_ptr<arrow::StructArray> pointStruct_(const std::shared_ptr<arrow::Table>& table)
  {
    auto col = table->GetColumnByName("point");
    if (! col) return nullptr;
    if (col->num_chunks() == 0) return nullptr;
    if (col->num_chunks() == 1)
    {
      if (col->chunk(0)->type_id() != arrow::Type::STRUCT) return nullptr;
      return std::static_pointer_cast<arrow::StructArray>(col->chunk(0));
    }
    // Multiple row groups: concatenate the chunks into one contiguous array so
    // the spectrum_index grouping scan sees the whole table.
    auto concat = arrow::Concatenate(col->chunks(), arrow::default_memory_pool());
    if (! concat.ok()) return nullptr;
    auto arr = concat.ValueOrDie();
    if (arr->type_id() != arrow::Type::STRUCT) return nullptr;
    return std::static_pointer_cast<arrow::StructArray>(arr);
  }

  /// Read the per-spectrum metadata map (index -> {ms_level, RT, betas}) from
  /// the spectra_metadata.parquet "spectrum" struct column.
  std::map<std::uint64_t, SpectrumMeta> readSpectraMeta_(const std::shared_ptr<arrow::Table>& table)
  {
    std::map<std::uint64_t, SpectrumMeta> out;

    auto col = table->GetColumnByName("spectrum");
    if (! col) return out;

    for (const auto& chunk : col->chunks())
    {
      if (chunk->type_id() != arrow::Type::STRUCT) continue;
      auto spectrum = std::static_pointer_cast<arrow::StructArray>(chunk);

      auto index_field = spectrum->GetFieldByName("index");
      if (! index_field || index_field->type_id() != arrow::Type::UINT64) continue;
      auto index = std::static_pointer_cast<arrow::UInt64Array>(index_field);

      auto level_field = spectrum->GetFieldByName("MS_1000511_ms_level");
      auto time_field = spectrum->GetFieldByName("time");
      auto model_field = spectrum->GetFieldByName("mz_delta_model");

      auto large_list = std::dynamic_pointer_cast<arrow::LargeListArray>(model_field);
      auto list = std::dynamic_pointer_cast<arrow::ListArray>(model_field);
      std::shared_ptr<arrow::DoubleArray> model_values;
      if (large_list) model_values = std::dynamic_pointer_cast<arrow::DoubleArray>(large_list->values());
      else if (list)
        model_values = std::dynamic_pointer_cast<arrow::DoubleArray>(list->values());

      for (int64_t r = 0; r < spectrum->length(); ++r)
      {
        if (index->IsNull(r)) continue;
        SpectrumMeta m;

        if (level_field && ! level_field->IsNull(r))
        {
          switch (level_field->type_id())
          {
            case arrow::Type::UINT8:
              m.ms_level = static_cast<int>(std::static_pointer_cast<arrow::UInt8Array>(level_field)->Value(r));
              break;
            case arrow::Type::INT8:
              m.ms_level = static_cast<int>(std::static_pointer_cast<arrow::Int8Array>(level_field)->Value(r));
              break;
            case arrow::Type::UINT16:
              m.ms_level = static_cast<int>(std::static_pointer_cast<arrow::UInt16Array>(level_field)->Value(r));
              break;
            case arrow::Type::INT32:
              m.ms_level = static_cast<int>(std::static_pointer_cast<arrow::Int32Array>(level_field)->Value(r));
              break;
            case arrow::Type::INT64:
              m.ms_level = static_cast<int>(std::static_pointer_cast<arrow::Int64Array>(level_field)->Value(r));
              break;
            default:
              break;
          }
        }

        if (time_field && ! time_field->IsNull(r) && time_field->type_id() == arrow::Type::DOUBLE)
        {
          m.retention_time = std::static_pointer_cast<arrow::DoubleArray>(time_field)->Value(r);
        }

        if (model_values)
        {
          bool null_model = large_list ? large_list->IsNull(r) : (list ? list->IsNull(r) : true);
          if (! null_model)
          {
            int64_t offset = large_list ? large_list->value_offset(r) : list->value_offset(r);
            int64_t length = large_list ? large_list->value_length(r) : list->value_length(r);
            m.mz_delta_model.reserve(static_cast<std::size_t>(length));
            for (int64_t k = 0; k < length; ++k)
            {
              m.mz_delta_model.push_back(model_values->Value(offset + k));
            }
          }
        }

        out[index->Value(r)] = std::move(m);
      }
    }

    return out;
  }

  /// Decode all point-layout spectra from a data/peaks table and add them to
  /// @p exp. Rows are grouped by spectrum_index. When @p reconstruct is true
  /// (profile data), interior null m/z values are reconstructed from the
  /// per-spectrum mz_delta_model; centroid peaks carry no nulls.
  void addSpectraFromTable_(const std::shared_ptr<arrow::Table>& table,
                            const std::map<std::uint64_t, SpectrumMeta>& meta,
                            bool reconstruct,
                            MSExperiment& exp)
  {
    auto pts = pointStruct_(table);
    if (! pts) return;

    auto si_field = pts->GetFieldByName("spectrum_index");
    auto mz_field = pts->GetFieldByName("mz");
    auto int_field = pts->GetFieldByName("intensity");
    if (! si_field || ! mz_field || ! int_field) return;
    if (si_field->type_id() != arrow::Type::UINT64) return;
    if (mz_field->type_id() != arrow::Type::DOUBLE) return;
    if (int_field->type_id() != arrow::Type::FLOAT) return;

    auto si = std::static_pointer_cast<arrow::UInt64Array>(si_field);
    auto mz = std::static_pointer_cast<arrow::DoubleArray>(mz_field);
    auto intensity = std::static_pointer_cast<arrow::FloatArray>(int_field);

    const int64_t n = pts->length();
    int64_t i = 0;
    while (i < n)
    {
      // Rows are grouped by spectrum_index; gather the contiguous run.
      std::uint64_t cur = si->Value(i);
      int64_t j = i;

      std::vector<double> mz_values;
      std::vector<bool> mz_valid;
      std::vector<float> int_values;
      while (j < n && si->Value(j) == cur)
      {
        bool ok = ! mz->IsNull(j);
        mz_valid.push_back(ok);
        mz_values.push_back(ok ? mz->Value(j) : 0.0);
        int_values.push_back(intensity->IsNull(j) ? 0.0f : intensity->Value(j));
        ++j;
      }
      i = j;

      // Reconstruct null-marked m/z for profile spectra.
      std::vector<double> final_mz = mz_values;
      if (reconstruct && std::find(mz_valid.begin(), mz_valid.end(), false) != mz_valid.end())
      {
        std::vector<double> betas;
        if (auto it = meta.find(cur); it != meta.end()) betas = it->second.mz_delta_model;
        auto rec = reconstruct_null_mz_(mz_values, mz_valid, betas);
        if (! rec.empty()) final_mz = std::move(rec);
      }

      MSSpectrum spec;
      spec.reserve(final_mz.size());
      for (std::size_t k = 0; k < final_mz.size(); ++k)
      {
        spec.push_back(Peak1D(final_mz[k], int_values[k]));
      }

      if (auto it = meta.find(cur); it != meta.end())
      {
        if (it->second.ms_level > 0) spec.setMSLevel(static_cast<UInt>(it->second.ms_level));
        spec.setRT(it->second.retention_time);
      }
      spec.setNativeID("index=" + String(cur));
      spec.sortByPosition();
      exp.addSpectrum(std::move(spec));
    }
  }

} // namespace

MzPeakFile::MzPeakFile() = default;

MzPeakFile::~MzPeakFile() = default;

void MzPeakFile::load(const String& filename, MapType& map) const
{
  if (! File::exists(filename)) { throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename); }

  map.clear(true);

  // Keep a single TempDir alive for all archive-entry reads (ZipRandomAccessFile
  // may extract entries to disk on some platforms).
  std::unique_ptr<File::TempDir> temp_dir;

  // ------------------------------------------------------------------
  // 1. Parse mzpeak_index.json to locate the data / peaks / metadata members.
  // ------------------------------------------------------------------
  std::string index_json = readEntryBytes_(filename, "mzpeak_index.json", temp_dir);

  String data_entry;     // profile data arrays
  String peaks_entry;    // centroid peaks
  String metadata_entry; // per-spectrum metadata
  try
  {
    nlohmann::json idx = nlohmann::json::parse(index_json);
    for (const auto& f : idx.at("files"))
    {
      std::string name = f.value("name", "");
      std::string entity = f.value("entity_type", "");
      std::string kind = f.value("data_kind", "");
      if (entity != "spectrum") continue; // chromatograms ignored in this plan
      if (kind == "data arrays") data_entry = name;
      else if (kind == "peaks")
        peaks_entry = name;
      else if (kind == "metadata")
        metadata_entry = name;
    }
  }
  catch (const std::exception& e)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, std::string("Failed to parse mzpeak_index.json: ") + e.what());
  }

  // ------------------------------------------------------------------
  // 2. Read per-spectrum metadata (ms_level, RT seconds, mz_delta_model).
  // ------------------------------------------------------------------
  std::map<std::uint64_t, SpectrumMeta> meta;
  if (! metadata_entry.empty())
  {
    auto meta_table = readTableFromArchive_(filename, metadata_entry, temp_dir);
    meta = readSpectraMeta_(meta_table);
  }

  // ------------------------------------------------------------------
  // 3. Decode profile data arrays (null-marked m/z reconstructed).
  // ------------------------------------------------------------------
  if (! data_entry.empty())
  {
    auto data_table = readTableFromArchive_(filename, data_entry, temp_dir);
    addSpectraFromTable_(data_table, meta, /*reconstruct=*/true, map);
  }

  // ------------------------------------------------------------------
  // 4. Decode centroid peaks (no null complications).
  // ------------------------------------------------------------------
  if (! peaks_entry.empty())
  {
    auto peaks_table = readTableFromArchive_(filename, peaks_entry, temp_dir);
    addSpectraFromTable_(peaks_table, meta, /*reconstruct=*/false, map);
  }

  // Spectra are appended data-table-first then peaks-table; order them by RT so
  // the experiment follows acquisition order.
  map.sortSpectra(false);
  map.updateRanges();
}

void MzPeakFile::store(const String& /* filename */, const MapType& /* map */) const
{ throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }

void MzPeakFile::transform(const String& /* filename_in */,
                           Interfaces::IMSDataConsumer* /* consumer */,
                           bool /* skip_full_count */,
                           bool /* skip_first_pass */) const
{ throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }

} // namespace OpenMS
