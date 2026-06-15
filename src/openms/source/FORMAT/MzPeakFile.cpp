// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/FORMAT/MzPeakFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/FORMAT/ZipRandomAccessFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/CVTerm.h>
#include <OpenMS/METADATA/Instrument.h>
#include <OpenMS/METADATA/IonDetector.h>
#include <OpenMS/METADATA/IonSource.h>
#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/METADATA/Software.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <algorithm>
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>
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

  /// A flat PSI-MS CV parameter tuple, mirroring mzPeak's {accession, name,
  /// value, unit} CvParam. The value is normalised to a string so the CV
  /// dispatch can apply typed setters (toDouble/toInt) uniformly, exactly like
  /// MzMLHandler consumes its XML cvParam @c value attribute.
  struct CvParam
  {
    String accession;
    String name;
    String value; ///< stringified union value ("" if absent)
    String unit;  ///< unit accession ("" if absent)
  };

  /// One precursor (RDR-10c): an isolation window (target m/z + offsets),
  /// activation params, and one selected ion. mzPeak emits these as separate
  /// row-aligned facet columns joined on @c source_index; we map every selected
  /// ion to its own OpenMS Precursor so MSn/DIA precursor info is never dropped.
  struct PrecursorData
  {
    bool has_isolation = false;
    double isolation_target_mz = 0.0;
    double isolation_lower_offset = 0.0;
    double isolation_upper_offset = 0.0;
    std::vector<CvParam> activation;

    bool has_selected_ion = false;
    double selected_ion_mz = 0.0;
    bool has_charge = false;
    int charge = 0;
    bool has_intensity = false;
    double intensity = 0.0;
  };

  /// Per-spectrum metadata needed to build an MSSpectrum.
  struct SpectrumMeta
  {
    int ms_level = 0;
    double retention_time = 0.0; ///< seconds (mzpeak stores RT in seconds)
    std::vector<double> mz_delta_model;
    String native_id;
    String representation;           ///< MS:1000127 centroid / MS:1000128 profile
    int polarity = 0;                ///< +1 positive, -1 negative, 0 unknown
    std::vector<CvParam> parameters; ///< per-spectrum flat CV params
    std::vector<PrecursorData> precursors;
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

  /// Stringify one element of a mzPeak value union
  /// struct<integer,float,string,boolean>. mzPeak stores a CvParam value in
  /// exactly one populated sub-field; we normalise it to a String so the CV
  /// dispatch can apply typed setters uniformly. Returns "" if all are null.
  String valueUnionToString_(const std::shared_ptr<arrow::StructArray>& value_struct, int64_t row)
  {
    if (! value_struct || value_struct->IsNull(row)) return "";

    if (auto f = value_struct->GetFieldByName("string"))
    {
      if (! f->IsNull(row) && f->type_id() == arrow::Type::LARGE_STRING)
      {
        return String(std::static_pointer_cast<arrow::LargeStringArray>(f)->GetString(row));
      }
    }
    if (auto f = value_struct->GetFieldByName("float"))
    {
      if (! f->IsNull(row) && f->type_id() == arrow::Type::DOUBLE) { return String(std::static_pointer_cast<arrow::DoubleArray>(f)->Value(row)); }
    }
    if (auto f = value_struct->GetFieldByName("integer"))
    {
      if (! f->IsNull(row) && f->type_id() == arrow::Type::INT64) { return String(std::static_pointer_cast<arrow::Int64Array>(f)->Value(row)); }
    }
    if (auto f = value_struct->GetFieldByName("boolean"))
    {
      if (! f->IsNull(row) && f->type_id() == arrow::Type::BOOL)
      {
        return std::static_pointer_cast<arrow::BooleanArray>(f)->Value(row) ? String("true") : String("false");
      }
    }
    return "";
  }

  /// Read a string sub-field of a struct array at @p row, accepting either
  /// arrow string or large_string. Returns "" if absent/null.
  String structString_(const std::shared_ptr<arrow::StructArray>& s, const String& field, int64_t row)
  {
    if (! s) return "";
    auto f = s->GetFieldByName(field);
    if (! f || f->IsNull(row)) return "";
    if (f->type_id() == arrow::Type::LARGE_STRING) { return String(std::static_pointer_cast<arrow::LargeStringArray>(f)->GetString(row)); }
    if (f->type_id() == arrow::Type::STRING) { return String(std::static_pointer_cast<arrow::StringArray>(f)->GetString(row)); }
    return "";
  }

  /// Decode a CvParam list column (a large_list<struct<value, accession, name,
  /// unit>>) at the given outer @p row into a flat vector<CvParam>. Mirrors the
  /// mzpeak reader's parameter-facet decode. Empty/null lists yield {}.
  std::vector<CvParam> readCvParamList_(const std::shared_ptr<arrow::Array>& list_field, int64_t row)
  {
    std::vector<CvParam> out;
    auto large_list = std::dynamic_pointer_cast<arrow::LargeListArray>(list_field);
    if (! large_list || large_list->IsNull(row)) return out;

    auto items = std::dynamic_pointer_cast<arrow::StructArray>(large_list->values());
    if (! items) return out;

    int64_t offset = large_list->value_offset(row);
    int64_t length = large_list->value_length(row);

    auto value_struct = std::dynamic_pointer_cast<arrow::StructArray>(items->GetFieldByName("value"));

    out.reserve(static_cast<std::size_t>(length));
    for (int64_t k = 0; k < length; ++k)
    {
      int64_t r = offset + k;
      CvParam p;
      p.accession = structString_(items, "accession", r);
      p.name = structString_(items, "name", r);
      p.unit = structString_(items, "unit", r);
      p.value = valueUnionToString_(value_struct, r);
      out.push_back(std::move(p));
    }
    return out;
  }

  /// Read the precursor / selected_ion facet columns and join them by
  /// @c source_index (the owning spectrum's index). mzPeak stores these as
  /// row-aligned top-level struct columns separate from the @c spectrum column;
  /// a single spectrum may own several precursor rows (MSn/DIA). Returns a map
  /// source_index -> vector<PrecursorData> (one entry per precursor row).
  std::map<std::uint64_t, std::vector<PrecursorData>> readPrecursors_(const std::shared_ptr<arrow::Table>& table)
  {
    std::map<std::uint64_t, std::vector<PrecursorData>> out;

    // Index precursor rows by (source_index, precursor_index) so the selected
    // ion facet can be matched back onto the right precursor.
    std::map<std::pair<std::uint64_t, std::uint64_t>, PrecursorData> joined;

    auto prec_col = table->GetColumnByName("precursor");
    if (prec_col)
    {
      for (const auto& chunk : prec_col->chunks())
      {
        auto prec = std::dynamic_pointer_cast<arrow::StructArray>(chunk);
        if (! prec) continue;

        auto si = std::dynamic_pointer_cast<arrow::UInt64Array>(prec->GetFieldByName("source_index"));
        auto pi = std::dynamic_pointer_cast<arrow::UInt64Array>(prec->GetFieldByName("precursor_index"));
        auto iso = std::dynamic_pointer_cast<arrow::StructArray>(prec->GetFieldByName("isolation_window"));
        auto act = std::dynamic_pointer_cast<arrow::StructArray>(prec->GetFieldByName("activation"));
        if (! si || ! pi) continue;

        std::shared_ptr<arrow::FloatArray> tgt, lo, hi;
        if (iso)
        {
          tgt = std::dynamic_pointer_cast<arrow::FloatArray>(iso->GetFieldByName("MS_1000827_isolation_window_target_mz"));
          lo = std::dynamic_pointer_cast<arrow::FloatArray>(iso->GetFieldByName("MS_1000828_isolation_window_lower_offset"));
          hi = std::dynamic_pointer_cast<arrow::FloatArray>(iso->GetFieldByName("MS_1000829_isolation_window_upper_offset"));
        }
        auto act_params = act ? act->GetFieldByName("parameters") : nullptr;

        for (int64_t r = 0; r < prec->length(); ++r)
        {
          if (si->IsNull(r) || pi->IsNull(r)) continue;
          PrecursorData pd;
          if (tgt && ! tgt->IsNull(r))
          {
            pd.has_isolation = true;
            pd.isolation_target_mz = tgt->Value(r);
            pd.isolation_lower_offset = (lo && ! lo->IsNull(r)) ? lo->Value(r) : 0.0;
            pd.isolation_upper_offset = (hi && ! hi->IsNull(r)) ? hi->Value(r) : 0.0;
          }
          if (act_params) pd.activation = readCvParamList_(act_params, r);
          joined[{si->Value(r), pi->Value(r)}] = std::move(pd);
        }
      }
    }

    auto si_col = table->GetColumnByName("selected_ion");
    if (si_col)
    {
      for (const auto& chunk : si_col->chunks())
      {
        auto sion = std::dynamic_pointer_cast<arrow::StructArray>(chunk);
        if (! sion) continue;

        auto si = std::dynamic_pointer_cast<arrow::UInt64Array>(sion->GetFieldByName("source_index"));
        auto pi = std::dynamic_pointer_cast<arrow::UInt64Array>(sion->GetFieldByName("precursor_index"));
        auto mz = std::dynamic_pointer_cast<arrow::DoubleArray>(sion->GetFieldByName("MS_1000744_selected_ion_mz_unit_MS_1000040"));
        auto charge = std::dynamic_pointer_cast<arrow::Int32Array>(sion->GetFieldByName("MS_1000041_charge_state"));
        auto inten = std::dynamic_pointer_cast<arrow::FloatArray>(sion->GetFieldByName("MS_1000042_intensity_unit_MS_1000131"));
        if (! si || ! pi) continue;

        for (int64_t r = 0; r < sion->length(); ++r)
        {
          if (si->IsNull(r) || pi->IsNull(r)) continue;
          // Attach to the matching precursor row (or create a bare one if the
          // precursor facet had no isolation window for this index pair).
          PrecursorData& pd = joined[{si->Value(r), pi->Value(r)}];
          if (mz && ! mz->IsNull(r))
          {
            pd.has_selected_ion = true;
            pd.selected_ion_mz = mz->Value(r);
          }
          if (charge && ! charge->IsNull(r))
          {
            pd.has_charge = true;
            pd.charge = charge->Value(r);
          }
          if (inten && ! inten->IsNull(r))
          {
            pd.has_intensity = true;
            pd.intensity = inten->Value(r);
          }
        }
      }
    }

    for (auto& [key, pd] : joined)
    {
      out[key.first].push_back(std::move(pd));
    }
    return out;
  }

  /// Read the per-spectrum metadata map (index -> {ms_level, RT, betas, ...})
  /// from the spectra_metadata.parquet "spectrum" struct column.
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
      auto polarity_field = spectrum->GetFieldByName("MS_1000465_scan_polarity");
      auto params_field = spectrum->GetFieldByName("parameters");

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

        m.native_id = structString_(spectrum, "id", r);
        m.representation = structString_(spectrum, "MS_1000525_spectrum_representation", r);

        if (polarity_field && ! polarity_field->IsNull(r) && polarity_field->type_id() == arrow::Type::INT8)
        {
          m.polarity = static_cast<int>(std::static_pointer_cast<arrow::Int8Array>(polarity_field)->Value(r));
        }

        if (params_field) m.parameters = readCvParamList_(params_field, r);

        out[index->Value(r)] = std::move(m);
      }
    }

    return out;
  }

  // ==========================================================================
  // CV-param dispatch (§4a M2): mirror the MzMLHandler::handleCVParam_ subset
  // mzPeak actually emits. Recognised accession -> typed setter; everything
  // else -> setMetaValue keyed by the bare accession (mirrors
  // handleUserParam_, so the accession is never lost).
  // ==========================================================================

  /// Spectrum-context CV dispatch (parent_tag "spectrum"/"scan" in mzML).
  void applyCVParamToSpectrum_(const CvParam& p, MSSpectrum& spec)
  {
    const String& a = p.accession;
    if (a.empty()) // no accession (e.g. Thermo trailer) -> keep name->value
    {
      if (! p.name.empty()) spec.setMetaValue(p.name, p.value);
      return;
    }
    else if (a == "MS:1000127") { spec.setType(SpectrumSettings::SpectrumType::CENTROID); }
    else if (a == "MS:1000128") { spec.setType(SpectrumSettings::SpectrumType::PROFILE); }
    else if (a == "MS:1000129") { spec.getInstrumentSettings().setPolarity(IonSource::Polarity::NEGATIVE); }
    else if (a == "MS:1000130") { spec.getInstrumentSettings().setPolarity(IonSource::Polarity::POSITIVE); }
    else
    {
      spec.setMetaValue(a, p.value);
    } // keep accession-keyed (handleUserParam_)
  }

  /// Activation-context CV dispatch: map dissociation methods + collision
  /// energy onto an OpenMS Precursor (parent_tag "activation" in mzML).
  void applyActivationCVParam_(const CvParam& p, Precursor& prec)
  {
    const String& a = p.accession;
    if (a == "MS:1000133") { prec.getActivationMethods().insert(Precursor::ActivationMethod::CID); }
    else if (a == "MS:1000134") { prec.getActivationMethods().insert(Precursor::ActivationMethod::PD); }
    else if (a == "MS:1000135") { prec.getActivationMethods().insert(Precursor::ActivationMethod::PSD); }
    else if (a == "MS:1000136") { prec.getActivationMethods().insert(Precursor::ActivationMethod::SID); }
    else if (a == "MS:1000242") { prec.getActivationMethods().insert(Precursor::ActivationMethod::BIRD); }
    else if (a == "MS:1000250") { prec.getActivationMethods().insert(Precursor::ActivationMethod::ECD); }
    else if (a == "MS:1000262") { prec.getActivationMethods().insert(Precursor::ActivationMethod::IMD); }
    else if (a == "MS:1000282") { prec.getActivationMethods().insert(Precursor::ActivationMethod::SORI); }
    else if (a == "MS:1000422") { prec.getActivationMethods().insert(Precursor::ActivationMethod::HCD); }
    else if (a == "MS:1000598") { prec.getActivationMethods().insert(Precursor::ActivationMethod::ETD); }
    else if (a == "MS:1000599") { prec.getActivationMethods().insert(Precursor::ActivationMethod::PQD); }
    else if (a == "MS:1002472") { prec.getActivationMethods().insert(Precursor::ActivationMethod::TRAP); }
    else if (a == "MS:1000045") // collision energy
    {
      if (! p.value.empty()) prec.setActivationEnergy(p.value.toDouble());
    }
    else if (a == "MS:1000509") // activation energy
    {
      if (! p.value.empty()) prec.setActivationEnergy(p.value.toDouble());
    }
    else if (! a.empty()) { prec.setMetaValue(a, p.value); } // keep otherwise
  }

  /// Build OpenMS Precursors from the joined precursor/selected_ion facets.
  /// Isolation window is stored as OFFSETS from the target m/z (the OpenMS
  /// convention): setMZ(target), setIsolationWindow{Lower,Upper}Offset(offset).
  /// Every selected ion becomes its own Precursor so MSn/DIA info is preserved.
  void buildPrecursors_(const std::vector<PrecursorData>& src, MSSpectrum& spec)
  {
    for (const PrecursorData& pd : src)
    {
      Precursor prec;
      // Isolation target m/z is the primary m/z; offsets are widths.
      if (pd.has_isolation)
      {
        prec.setMZ(pd.isolation_target_mz);
        prec.setIsolationWindowLowerOffset(pd.isolation_lower_offset);
        prec.setIsolationWindowUpperOffset(pd.isolation_upper_offset);
      }
      // A selected-ion m/z overrides the isolation target as the precursor m/z
      // (the mzML "precursor m/z from selected ion" behaviour); keep the
      // isolation target as a meta value so it is not lost.
      if (pd.has_selected_ion)
      {
        if (pd.has_isolation) prec.setMetaValue("isolation window target m/z", pd.isolation_target_mz);
        prec.setMZ(pd.selected_ion_mz);
      }
      if (pd.has_charge) prec.setCharge(pd.charge);
      if (pd.has_intensity) prec.setIntensity(pd.intensity);
      for (const CvParam& p : pd.activation)
        applyActivationCVParam_(p, prec);
      spec.getPrecursors().push_back(std::move(prec));
    }
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
        const SpectrumMeta& sm = it->second;
        if (sm.ms_level > 0) spec.setMSLevel(static_cast<UInt>(sm.ms_level));
        spec.setRT(sm.retention_time);

        // Spectrum type from the representation accession (MS:1000127/1000128).
        if (sm.representation == "MS:1000127") spec.setType(SpectrumSettings::SpectrumType::CENTROID);
        else if (sm.representation == "MS:1000128")
          spec.setType(SpectrumSettings::SpectrumType::PROFILE);

        // Polarity (mzPeak scan polarity: +1 positive, -1 negative).
        if (sm.polarity > 0) spec.getInstrumentSettings().setPolarity(IonSource::Polarity::POSITIVE);
        else if (sm.polarity < 0)
          spec.getInstrumentSettings().setPolarity(IonSource::Polarity::NEGATIVE);

        // Per-spectrum flat CV params via the dispatch helper.
        for (const CvParam& p : sm.parameters)
          applyCVParamToSpectrum_(p, spec);

        // Precursors (isolation offsets + activation + every selected ion).
        buildPrecursors_(sm.precursors, spec);

        // Keep the mzPeak native id (controllerType=...) without breaking the
        // stable "index=N" id used for lookups.
        if (! sm.native_id.empty()) spec.setMetaValue("mzpeak_native_id", sm.native_id);
      }
      spec.setNativeID("index=" + String(cur));
      spec.sortByPosition();
      exp.addSpectrum(std::move(spec));
    }
  }

  // ==========================================================================
  // Run-level metadata{} -> ExperimentalSettings (§4 mapping).
  // ==========================================================================

  /// Convenience accessors for nlohmann JSON CvParam tuples ({accession, name,
  /// value, unit}); value may be a string/number/null.
  String jsonAccession_(const nlohmann::json& p)
  { return p.value("accession", String()); }
  String jsonName_(const nlohmann::json& p)
  { return p.value("name", String()); }
  String jsonParamValue_(const nlohmann::json& p)
  {
    if (! p.contains("value") || p.at("value").is_null()) return "";
    const auto& v = p.at("value");
    if (v.is_string()) return String(v.get<std::string>());
    if (v.is_number_integer()) return String(v.get<long long>());
    if (v.is_number_float()) return String(v.get<double>());
    if (v.is_boolean()) return v.get<bool>() ? String("true") : String("false");
    return "";
  }

  /// Map a source-file CvParam (subset MzMLHandler recognises) onto a SourceFile.
  /// SourceFile IS-A CVTermList; unrecognised accessions are kept via addCVTerm.
  void applySourceFileParam_(const nlohmann::json& p, SourceFile& sf)
  {
    String a = jsonAccession_(p);
    String value = jsonParamValue_(p);
    if (a == "MS:1000569") { sf.setChecksum(value, SourceFile::ChecksumType::SHA1); }
    else if (a == "MS:1000568") { sf.setChecksum(value, SourceFile::ChecksumType::MD5); }
    else if (a == "MS:1000563") { sf.setFileType(jsonName_(p)); } // Thermo RAW format
    else if (a == "MS:1000584") { sf.setFileType(jsonName_(p)); } // mzML format
    else if (a == "MS:1000768")                                   // Thermo nativeID format
    {
      sf.setNativeIDType(jsonName_(p));
      sf.setNativeIDTypeAccession(a);
    }
    else if (! a.empty()) { sf.addCVTerm(CVTerm(a, jsonName_(p), "MS", value)); }
  }

  /// Map the run-level metadata{} object onto the MSExperiment's
  /// ExperimentalSettings: run date/ids, instrument (+ components), source
  /// files, sample, and software. Best-effort; missing blocks are skipped.
  void applyRunMetadata_(const nlohmann::json& md, MSExperiment& exp)
  {
    // ---- run: date + identifier ------------------------------------------
    if (md.contains("run") && md.at("run").is_object())
    {
      const auto& run = md.at("run");
      if (run.contains("start_time") && run.at("start_time").is_string())
      {
        // ISO 8601 (e.g. 2005-07-20T19:44:22Z); DateTime::set tolerates the 'Z'.
        DateTime dt;
        String iso = run.value("start_time", String());
        iso.substitute("Z", "");
        iso.substitute("T", " ");
        try
        {
          dt.set(iso);
          exp.setDateTime(dt);
        }
        catch (const std::exception&)
        { /* leave default date if unparseable */
        }
      }
      if (run.contains("id")) exp.setMetaValue("mzpeak_run_id", run.value("id", String()));
    }

    // ---- instrument_configuration_list -> Instrument ---------------------
    if (md.contains("instrument_configuration_list") && md.at("instrument_configuration_list").is_array()
        && ! md.at("instrument_configuration_list").empty())
    {
      const auto& ic = md.at("instrument_configuration_list").front();
      Instrument instrument;

      // Instrument-level params: a named model term -> setName; others kept.
      if (ic.contains("parameters"))
      {
        for (const auto& p : ic.at("parameters"))
        {
          String a = jsonAccession_(p);
          String value = jsonParamValue_(p);
          if (a == "MS:1000529") { instrument.setMetaValue("instrument serial number", value); }
          else if (! jsonName_(p).empty() && value.empty())
          {
            instrument.setName(jsonName_(p)); // model term (e.g. "LTQ FT")
          }
          else if (! a.empty()) { instrument.setMetaValue(a, value); }
        }
      }

      // Components -> IonSource / MassAnalyzer / IonDetector (order preserved;
      // CV accessions kept as meta values since these hosts are not CVTermList).
      if (ic.contains("components"))
      {
        std::vector<IonSource> sources;
        std::vector<MassAnalyzer> analyzers;
        std::vector<IonDetector> detectors;
        for (const auto& c : ic.at("components"))
        {
          String ctype = c.value("component_type", String());
          int order = c.value("order", 0);
          auto first_acc = [&](String& acc, String& nm) {
            if (c.contains("parameters") && ! c.at("parameters").empty())
            {
              acc = jsonAccession_(c.at("parameters").front());
              nm = jsonName_(c.at("parameters").front());
            }
          };
          String acc, nm;
          first_acc(acc, nm);
          if (ctype == "ionsource")
          {
            IonSource s;
            s.setOrder(order);
            if (! acc.empty()) s.setMetaValue(acc, nm);
            sources.push_back(s);
          }
          else if (ctype == "analyzer")
          {
            MassAnalyzer m;
            m.setOrder(order);
            if (! acc.empty()) m.setMetaValue(acc, nm);
            analyzers.push_back(m);
          }
          else if (ctype == "detector")
          {
            IonDetector d;
            d.setOrder(order);
            if (! acc.empty()) d.setMetaValue(acc, nm);
            detectors.push_back(d);
          }
        }
        if (! sources.empty()) instrument.setIonSources(sources);
        if (! analyzers.empty()) instrument.setMassAnalyzers(analyzers);
        if (! detectors.empty()) instrument.setIonDetectors(detectors);
      }

      exp.setInstrument(instrument);
    }

    // ---- file_description.source_files -> setSourceFiles -----------------
    if (md.contains("file_description") && md.at("file_description").is_object())
    {
      const auto& fd = md.at("file_description");
      if (fd.contains("source_files") && fd.at("source_files").is_array())
      {
        std::vector<SourceFile> source_files;
        for (const auto& sfj : fd.at("source_files"))
        {
          SourceFile sf;
          sf.setNameOfFile(sfj.value("name", String()));
          sf.setPathToFile(sfj.value("location", String()));
          if (sfj.contains("parameters"))
          {
            for (const auto& p : sfj.at("parameters"))
              applySourceFileParam_(p, sf);
          }
          source_files.push_back(sf);
        }
        if (! source_files.empty()) exp.setSourceFiles(source_files);
      }
    }

    // ---- sample_list -> setSample (first sample) -------------------------
    if (md.contains("sample_list") && md.at("sample_list").is_array() && ! md.at("sample_list").empty())
    {
      const auto& sj = md.at("sample_list").front();
      Sample sample;
      sample.setName(sj.value("name", String()));
      if (sj.contains("id")) sample.setNumber(sj.value("id", String()));
      exp.setSample(sample);
    }

    // ---- software_list -> Instrument::setSoftware (first entry) ----------
    // Software has no run-level home in OpenMS; attach to the instrument.
    if (md.contains("software_list") && md.at("software_list").is_array() && ! md.at("software_list").empty())
    {
      const auto& swj = md.at("software_list").front();
      Software sw;
      sw.setName(swj.value("id", String()));
      sw.setVersion(swj.value("version", String()));
      if (swj.contains("parameters"))
      {
        for (const auto& p : swj.at("parameters"))
        {
          String a = jsonAccession_(p);
          if (! a.empty()) sw.addCVTerm(CVTerm(a, jsonName_(p), "MS", jsonParamValue_(p)));
        }
      }
      exp.getInstrument().setSoftware(sw);
    }
  }

  // ==========================================================================
  // Store path (03-01): build the point-layout Arrow tables + metadata table
  // from an MSExperiment, write them to Parquet members in a temp dir, and zip
  // them into a .mzpeak archive that MzPeakFile::load re-reads into an
  // equivalent experiment. Schemas mirror the known-good mzpeak-lib writer and,
  // crucially, the exact field names the load path above expects.
  // ==========================================================================

  /// Throw a uniform ParseError on a failed Arrow status during write.
  void checkArrowStatus_(const arrow::Status& status, const String& what)
  {
    if (! status.ok()) { throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "mzPeak store: " + what, status.ToString()); }
  }

  /// The flattened point columns of one table (data or peaks): the spectrum
  /// index repeated per point plus the parallel m/z and intensity vectors.
  struct PointColumns
  {
    std::vector<std::uint64_t> spectrum_index;
    std::vector<double> mz;
    std::vector<float> intensity;
  };

  /// Write the point-layout table (top-level `point` struct {spectrum_index
  /// uint64, mz float64, intensity float32}) to a Parquet file at @p path.
  /// Matches @c pointStruct_ / @c addSpectraFromTable_ in the load path: a
  /// single struct column whose three children carry the literal values.
  void writePointTable_(const String& path, const PointColumns& cols)
  {
    arrow::UInt64Builder index_builder;
    arrow::DoubleBuilder mz_builder;
    arrow::FloatBuilder intensity_builder;
    checkArrowStatus_(index_builder.AppendValues(cols.spectrum_index), "append spectrum_index");
    checkArrowStatus_(mz_builder.AppendValues(cols.mz), "append mz");
    checkArrowStatus_(intensity_builder.AppendValues(cols.intensity), "append intensity");

    std::shared_ptr<arrow::Array> index_array, mz_array, intensity_array;
    checkArrowStatus_(index_builder.Finish(&index_array), "finish spectrum_index");
    checkArrowStatus_(mz_builder.Finish(&mz_array), "finish mz");
    checkArrowStatus_(intensity_builder.Finish(&intensity_array), "finish intensity");

    arrow::FieldVector point_fields {
      arrow::field("spectrum_index", arrow::uint64(), /*nullable=*/true),
      arrow::field("mz", arrow::float64(), /*nullable=*/true),
      arrow::field("intensity", arrow::float32(), /*nullable=*/true),
    };
    auto point_result = arrow::StructArray::Make({index_array, mz_array, intensity_array}, point_fields);
    checkArrowStatus_(point_result.status(), "build point struct");
    std::shared_ptr<arrow::Array> point_array = point_result.ValueOrDie();

    auto schema = arrow::schema({arrow::field("point", point_array->type(), /*nullable=*/true)});
    auto table = arrow::Table::Make(schema, {point_array});

    auto outfile_result = arrow::io::FileOutputStream::Open(std::string(path));
    checkArrowStatus_(outfile_result.status(), "open output file " + path);
    std::shared_ptr<arrow::io::FileOutputStream> outfile = outfile_result.ValueOrDie();

    // One bounded row group; a positive size is required even for empty tables.
    constexpr int64_t kMaxRowGroup = 1 << 20;
    int64_t row_group_size = table->num_rows() > 0 ? std::min<int64_t>(table->num_rows(), kMaxRowGroup) : kMaxRowGroup;
    checkArrowStatus_(parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, row_group_size), "write point table " + path);
  }

  /// One row of the spectra_metadata table.
  struct MetaRow
  {
    std::uint64_t index = 0;
    String id;
    std::uint8_t ms_level = 0;
    double time = 0.0; ///< RT in seconds
    std::int8_t polarity = 0;
    String representation; ///< MS:1000127 centroid / MS:1000128 profile
    std::uint64_t n_points = 0;
    std::uint64_t n_peaks = 0;
  };

  /// Write the spectra_metadata table: a single `spectrum` struct column whose
  /// FIRST child is the uint64 index (the load path accesses it positionally /
  /// reads `time`, `MS_1000511_ms_level`, `MS_1000465_scan_polarity`,
  /// `MS_1000525_spectrum_representation`).
  void writeMetadataTable_(const String& path, const std::vector<MetaRow>& rows)
  {
    arrow::UInt64Builder index_b;
    arrow::LargeStringBuilder id_b;
    arrow::UInt8Builder level_b;
    arrow::DoubleBuilder time_b;
    arrow::Int8Builder polarity_b;
    arrow::StringBuilder repr_b;
    arrow::UInt64Builder npts_b;
    arrow::UInt64Builder npks_b;

    for (const MetaRow& r : rows)
    {
      checkArrowStatus_(index_b.Append(r.index), "append index");
      checkArrowStatus_(id_b.Append(std::string(r.id)), "append id");
      checkArrowStatus_(level_b.Append(r.ms_level), "append ms_level");
      checkArrowStatus_(time_b.Append(r.time), "append time");
      checkArrowStatus_(polarity_b.Append(r.polarity), "append polarity");
      checkArrowStatus_(repr_b.Append(std::string(r.representation)), "append representation");
      checkArrowStatus_(npts_b.Append(r.n_points), "append number_of_data_points");
      checkArrowStatus_(npks_b.Append(r.n_peaks), "append number_of_peaks");
    }

    std::shared_ptr<arrow::Array> index_a, id_a, level_a, time_a, polarity_a, repr_a, npts_a, npks_a;
    checkArrowStatus_(index_b.Finish(&index_a), "finish index");
    checkArrowStatus_(id_b.Finish(&id_a), "finish id");
    checkArrowStatus_(level_b.Finish(&level_a), "finish ms_level");
    checkArrowStatus_(time_b.Finish(&time_a), "finish time");
    checkArrowStatus_(polarity_b.Finish(&polarity_a), "finish polarity");
    checkArrowStatus_(repr_b.Finish(&repr_a), "finish representation");
    checkArrowStatus_(npts_b.Finish(&npts_a), "finish number_of_data_points");
    checkArrowStatus_(npks_b.Finish(&npks_a), "finish number_of_peaks");

    // `index` MUST be the first child (load reads it positionally as child 0).
    arrow::FieldVector spectrum_fields {
      arrow::field("index", arrow::uint64(), /*nullable=*/true),
      arrow::field("id", arrow::large_utf8(), /*nullable=*/true),
      arrow::field("MS_1000511_ms_level", arrow::uint8(), /*nullable=*/true),
      arrow::field("time", arrow::float64(), /*nullable=*/true),
      arrow::field("MS_1000465_scan_polarity", arrow::int8(), /*nullable=*/true),
      arrow::field("MS_1000525_spectrum_representation", arrow::utf8(), /*nullable=*/true),
      arrow::field("MS_1003060_number_of_data_points", arrow::uint64(), /*nullable=*/true),
      arrow::field("MS_1003059_number_of_peaks", arrow::uint64(), /*nullable=*/true),
    };
    auto spectrum_result = arrow::StructArray::Make({index_a, id_a, level_a, time_a, polarity_a, repr_a, npts_a, npks_a}, spectrum_fields);
    checkArrowStatus_(spectrum_result.status(), "build spectrum struct");
    std::shared_ptr<arrow::Array> spectrum_array = spectrum_result.ValueOrDie();

    auto schema = arrow::schema({arrow::field("spectrum", spectrum_array->type(), /*nullable=*/true)});
    auto table = arrow::Table::Make(schema, {spectrum_array});

    auto outfile_result = arrow::io::FileOutputStream::Open(std::string(path));
    checkArrowStatus_(outfile_result.status(), "open output file " + path);
    std::shared_ptr<arrow::io::FileOutputStream> outfile = outfile_result.ValueOrDie();

    constexpr int64_t kMaxRowGroup = 1 << 20;
    int64_t row_group_size = table->num_rows() > 0 ? std::min<int64_t>(table->num_rows(), kMaxRowGroup) : kMaxRowGroup;
    checkArrowStatus_(parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, row_group_size), "write metadata table " + path);
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
  // 2. Read per-spectrum metadata (ms_level, RT seconds, mz_delta_model,
  //    type/polarity/native-id/CV-params) + precursor facets joined by
  //    source_index.
  // ------------------------------------------------------------------
  std::map<std::uint64_t, SpectrumMeta> meta;
  if (! metadata_entry.empty())
  {
    auto meta_table = readTableFromArchive_(filename, metadata_entry, temp_dir);
    meta = readSpectraMeta_(meta_table);

    auto precursors = readPrecursors_(meta_table);
    for (auto& [source_index, plist] : precursors)
    {
      meta[source_index].precursors = std::move(plist);
    }
  }

  // ------------------------------------------------------------------
  // 2a. Run-level metadata{} -> ExperimentalSettings (instrument, source files,
  //     sample, software, run date/ids). map IS-A ExperimentalSettings.
  // ------------------------------------------------------------------
  try
  {
    nlohmann::json idx = nlohmann::json::parse(index_json);
    if (idx.contains("metadata")) applyRunMetadata_(idx.at("metadata"), map);
  }
  catch (const std::exception&)
  {
    // Run-level metadata is best-effort; a malformed block must not abort the
    // peak load (the spectra are the primary payload).
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

void MzPeakFile::store(const String& filename, const MapType& map) const
{
  // ------------------------------------------------------------------
  // 1. Split spectra into profile (-> data arrays) and centroid (-> peaks),
  //    flatten each into the point columns, and build the metadata rows. The
  //    point order per spectrum is ascending m/z (sortByPosition) so the load
  //    path reads strictly ascending values. The spectrum_index is the global
  //    spectrum position i, matching the "index=i" native id load emits.
  // ------------------------------------------------------------------
  PointColumns data;  // profile points
  PointColumns peaks; // centroid points
  std::vector<MetaRow> meta_rows;
  meta_rows.reserve(map.size());

  bool any_centroid = false;
  for (Size i = 0; i < map.size(); ++i)
  {
    MSSpectrum spec = map[i]; // copy so we can sort without mutating the input
    spec.sortByPosition();

    const bool centroid = spec.getType() == SpectrumSettings::SpectrumType::CENTROID;
    PointColumns& target = centroid ? peaks : data;
    if (centroid) any_centroid = true;

    for (const Peak1D& p : spec)
    {
      target.spectrum_index.push_back(static_cast<std::uint64_t>(i));
      target.mz.push_back(p.getMZ());
      target.intensity.push_back(p.getIntensity());
    }

    MetaRow row;
    row.index = static_cast<std::uint64_t>(i);
    // Prefer the spectrum native id; fall back to the stable "index=i" form.
    row.id = spec.getNativeID().empty() ? ("index=" + String(i)) : String(spec.getNativeID());
    row.ms_level = static_cast<std::uint8_t>(spec.getMSLevel());
    row.time = spec.getRT();
    // Representation accession: profile MS:1000128, centroid MS:1000127.
    row.representation = centroid ? "MS:1000127" : "MS:1000128";
    // Scan polarity as the mzPeak signed int8 (+1 positive, -1 negative).
    switch (spec.getInstrumentSettings().getPolarity())
    {
      case IonSource::Polarity::POSITIVE:
        row.polarity = 1;
        break;
      case IonSource::Polarity::NEGATIVE:
        row.polarity = -1;
        break;
      default:
        row.polarity = 0;
        break;
    }
    row.n_points = centroid ? 0 : static_cast<std::uint64_t>(spec.size());
    row.n_peaks = centroid ? static_cast<std::uint64_t>(spec.size()) : 0;
    meta_rows.push_back(std::move(row));
  }

  // ------------------------------------------------------------------
  // 2. Write the Parquet members into a temp dir.
  // ------------------------------------------------------------------
  File::TempDir temp_dir;
  const String dir = temp_dir.getPath();
  const String data_path = dir + "/spectra_data.parquet";
  const String peaks_path = dir + "/spectra_peaks.parquet";
  const String meta_path = dir + "/spectra_metadata.parquet";
  const String index_path = dir + "/mzpeak_index.json";

  // The data table is always present (may be empty); peaks only when a centroid
  // spectrum exists, mirroring the mzpeak-lib writer.
  writePointTable_(data_path, data);
  if (any_centroid) writePointTable_(peaks_path, peaks);
  writeMetadataTable_(meta_path, meta_rows);

  // ------------------------------------------------------------------
  // 3. Write mzpeak_index.json (files[] + minimal metadata{version}). The load
  //    path matches members by data_kind, so the names + kinds must be exact.
  // ------------------------------------------------------------------
  nlohmann::json idx;
  nlohmann::json files = nlohmann::json::array();
  files.push_back({{"name", "spectra_data.parquet"}, {"entity_type", "spectrum"}, {"data_kind", "data arrays"}});
  if (any_centroid) { files.push_back({{"name", "spectra_peaks.parquet"}, {"entity_type", "spectrum"}, {"data_kind", "peaks"}}); }
  files.push_back({{"name", "spectra_metadata.parquet"}, {"entity_type", "spectrum"}, {"data_kind", "metadata"}});
  idx["files"] = std::move(files);
  idx["metadata"] = {{"version", "0.9.0"}};

  {
    std::string index_json = idx.dump();
    std::ofstream out(std::string(index_path), std::ios::binary);
    if (! out) { throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index_path); }
    out << index_json;
    out.close();
    if (! out) { throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index_path); }
  }

  // ------------------------------------------------------------------
  // 4. Zip the temp dir into the output .mzpeak (STORE members; parquet is
  //    already internally compressed and the format mandates ZIP_CM_STORE so
  //    the members stay randomly accessible).
  // ------------------------------------------------------------------
  ZipArchiveFile::zipDirectory(dir, filename);
}

void MzPeakFile::transform(const String& /* filename_in */,
                           Interfaces::IMSDataConsumer* /* consumer */,
                           bool /* skip_full_count */,
                           bool /* skip_first_pass */) const
{ throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }

} // namespace OpenMS
