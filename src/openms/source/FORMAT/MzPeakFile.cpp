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

  /// Write a Parquet table with Rust-compat writer properties (ZSTD, page index,
  /// SortingColumn on the first column, store_schema) and optional file-level
  /// key-value metadata. Column 0 is used as the sort key (the `point` struct
  /// for data/peaks tables and the `spectrum` struct for the metadata table).
  void writeTableWithProps_(const String& path, const std::shared_ptr<arrow::Table>& table, const std::map<std::string, std::string>& file_kv)
  {
    auto outfile_result = arrow::io::FileOutputStream::Open(std::string(path));
    checkArrowStatus_(outfile_result.status(), "open output file " + path);
    std::shared_ptr<arrow::io::FileOutputStream> outfile = outfile_result.ValueOrDie();

    parquet::WriterProperties::Builder props_b;
    props_b.compression(arrow::Compression::ZSTD);
    props_b.enable_statistics();
    props_b.enable_write_page_index();
    props_b.set_sorting_columns({parquet::SortingColumn {0, false, false}});
    auto writer_props = props_b.build();
    auto arrow_props = parquet::ArrowWriterProperties::Builder().store_schema()->build();

    auto writer_result = parquet::arrow::FileWriter::Open(*table->schema(), arrow::default_memory_pool(), outfile, writer_props, arrow_props);
    checkArrowStatus_(writer_result.status(), "open FileWriter for " + path);
    auto writer = std::move(writer_result).ValueOrDie();

    constexpr int64_t kMaxRowGroup = 1 << 20;
    int64_t chunk = table->num_rows() > 0 ? std::min<int64_t>(table->num_rows(), kMaxRowGroup) : kMaxRowGroup;
    checkArrowStatus_(writer->WriteTable(*table, chunk), "write table " + path);

    if (! file_kv.empty())
    {
      std::vector<std::string> keys, values;
      for (const auto& [k, v] : file_kv)
      {
        keys.push_back(k);
        values.push_back(v);
      }
      auto kvm = std::make_shared<arrow::KeyValueMetadata>(std::move(keys), std::move(values));
      checkArrowStatus_(writer->AddKeyValueMetadata(kvm), "add kv metadata " + path);
    }

    checkArrowStatus_(writer->Close(), "close FileWriter for " + path);
  }

  /// Write the point-layout table (top-level `point` struct {spectrum_index
  /// uint64, mz float64, intensity float32}) to a Parquet file at @p path.
  /// Matches @c pointStruct_ / @c addSpectraFromTable_ in the load path: a
  /// single struct column whose three children carry the literal values.
  void writePointTable_(const String& path, const PointColumns& cols, const std::map<std::string, std::string>& file_kv)
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

    writeTableWithProps_(path, table, file_kv);
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

  /// Per-spectrum precursor data for the store path.
  struct PrecursorOut
  {
    bool has_precursor = false; ///< false → null source_index (MS1)
    std::uint64_t source_index = 0;
    std::uint64_t precursor_index = 0;
    bool has_isolation = false;
    float isolation_target = 0.0f;
    float isolation_lower = 0.0f;
    float isolation_upper = 0.0f;
    /// Activation CvParams: {accession, name, has_float_value, float_val}
    struct ActParam
    {
      String accession;
      String name;
      bool has_float = false;
      double float_val = 0.0;
    };
    std::vector<ActParam> activation;
    bool has_selected_ion = false;
    double selected_ion_mz = 0.0;
    bool has_charge = false;
    int charge = 0;
    bool has_intensity = false;
    float intensity = 0.0f;
  };

  /// Map an OpenMS ActivationMethod to its PSI-MS accession and name.
  /// Returns {"", ""} for unknown methods.
  std::pair<String, String> activationMethodCV_(Precursor::ActivationMethod method)
  {
    switch (method)
    {
      case Precursor::ActivationMethod::CID:
        return {"MS:1000133", "collision-induced dissociation"};
      case Precursor::ActivationMethod::PD:
        return {"MS:1000134", "plasma desorption"};
      case Precursor::ActivationMethod::PSD:
        return {"MS:1000135", "post-source decay"};
      case Precursor::ActivationMethod::SID:
        return {"MS:1000136", "surface-induced dissociation"};
      case Precursor::ActivationMethod::BIRD:
        return {"MS:1000242", "blackbody infrared radiative dissociation"};
      case Precursor::ActivationMethod::ECD:
        return {"MS:1000250", "electron capture dissociation"};
      case Precursor::ActivationMethod::IMD:
        return {"MS:1000262", "infrared multiphoton dissociation"};
      case Precursor::ActivationMethod::SORI:
        return {"MS:1000282", "sustained off-resonance irradiation"};
      case Precursor::ActivationMethod::HCD:
        return {"MS:1000422", "beam-type collision-induced dissociation"};
      case Precursor::ActivationMethod::ETD:
        return {"MS:1000598", "electron transfer dissociation"};
      case Precursor::ActivationMethod::PQD:
        return {"MS:1000599", "pulsed q dissociation"};
      case Precursor::ActivationMethod::TRAP:
        return {"MS:1002472", "trap-type collision-induced dissociation"};
      default:
        return {"", ""};
    }
  }

  /// Build a PrecursorOut from one MSSpectrum.
  /// Returns {has_precursor=false} for MS1 spectra (no precursors).
  PrecursorOut buildPrecursorOut_(const MSSpectrum& spec, std::uint64_t spectrum_index, std::uint64_t precursor_index)
  {
    PrecursorOut po;
    if (spec.getPrecursors().empty()) return po; // has_precursor stays false
    po.has_precursor = true;
    po.source_index = spectrum_index;
    po.precursor_index = precursor_index;

    const Precursor& prec = spec.getPrecursors()[0];

    // Isolation window: target from meta value if selected-ion overrode getMZ()
    po.has_isolation = true;
    if (prec.metaValueExists("isolation window target m/z"))
    {
      po.isolation_target = static_cast<float>(static_cast<double>(prec.getMetaValue("isolation window target m/z")));
    }
    else
    {
      po.isolation_target = static_cast<float>(prec.getMZ());
    }
    po.isolation_lower = static_cast<float>(prec.getIsolationWindowLowerOffset());
    po.isolation_upper = static_cast<float>(prec.getIsolationWindowUpperOffset());

    // Activation methods
    for (const auto& method : prec.getActivationMethods())
    {
      auto [acc, name] = activationMethodCV_(method);
      if (! acc.empty())
      {
        PrecursorOut::ActParam ap;
        ap.accession = acc;
        ap.name = name;
        ap.has_float = false;
        po.activation.push_back(std::move(ap));
      }
    }
    // Collision energy as an additional ActParam
    if (prec.getActivationEnergy() > 0.0)
    {
      PrecursorOut::ActParam ap;
      ap.accession = "MS:1000045";
      ap.name = "collision energy";
      ap.has_float = true;
      ap.float_val = prec.getActivationEnergy();
      po.activation.push_back(std::move(ap));
    }

    // Selected ion
    po.has_selected_ion = true;
    po.selected_ion_mz = prec.getMZ(); // getMZ() is already the selected-ion mz
    if (prec.getCharge() != 0)
    {
      po.has_charge = true;
      po.charge = prec.getCharge();
    }
    if (prec.getIntensity() != 0.0f)
    {
      po.has_intensity = true;
      po.intensity = prec.getIntensity();
    }

    return po;
  }

  /// Build the run-level metadata JSON from ExperimentalSettings.
  /// Best-effort: if the settings are empty, returns {{"version","0.9.0"}}.
  nlohmann::json buildRunMetadataJson_(const MSExperiment& map)
  {
    nlohmann::json md;
    bool has_anything = false;

    // ---- run: start_time + id -------------------------------------------
    {
      nlohmann::json run_obj = nlohmann::json::object();
      bool has_run = false;
      if (! map.getDateTime().isNull())
      {
        // Convert DateTime to ISO 8601 string: "YYYY-MM-DD HH:MM:SS" -> replace space with T, append Z
        String dt_str = map.getDateTime().get();
        dt_str.substitute(" ", "T");
        run_obj["start_time"] = std::string(dt_str) + "Z";
        has_run = true;
      }
      if (map.metaValueExists("mzpeak_run_id"))
      {
        run_obj["id"] = std::string(static_cast<String>(map.getMetaValue("mzpeak_run_id")));
        has_run = true;
      }
      if (has_run)
      {
        md["run"] = std::move(run_obj);
        has_anything = true;
      }
    }

    // ---- instrument_configuration_list ----------------------------------
    {
      const Instrument& instr = map.getInstrument();
      bool has_instr
        = ! instr.getName().empty() || ! instr.getIonSources().empty() || ! instr.getMassAnalyzers().empty() || ! instr.getIonDetectors().empty();
      if (has_instr)
      {
        nlohmann::json ic = nlohmann::json::object();
        nlohmann::json params = nlohmann::json::array();
        if (! instr.getName().empty()) { params.push_back({{"accession", ""}, {"name", std::string(instr.getName())}, {"value", nullptr}}); }
        ic["parameters"] = std::move(params);

        nlohmann::json components = nlohmann::json::array();
        int order = 1;
        for (const IonSource& s : instr.getIonSources())
        {
          nlohmann::json comp = nlohmann::json::object();
          comp["component_type"] = "ionsource";
          comp["order"] = s.getOrder() > 0 ? s.getOrder() : order;
          comp["parameters"] = nlohmann::json::array();
          components.push_back(std::move(comp));
          ++order;
        }
        for (const MassAnalyzer& a : instr.getMassAnalyzers())
        {
          nlohmann::json comp = nlohmann::json::object();
          comp["component_type"] = "analyzer";
          comp["order"] = a.getOrder() > 0 ? a.getOrder() : order;
          comp["parameters"] = nlohmann::json::array();
          components.push_back(std::move(comp));
          ++order;
        }
        for (const IonDetector& d : instr.getIonDetectors())
        {
          nlohmann::json comp = nlohmann::json::object();
          comp["component_type"] = "detector";
          comp["order"] = d.getOrder() > 0 ? d.getOrder() : order;
          comp["parameters"] = nlohmann::json::array();
          components.push_back(std::move(comp));
          ++order;
        }
        ic["components"] = std::move(components);

        md["instrument_configuration_list"] = nlohmann::json::array({std::move(ic)});
        has_anything = true;
      }
    }

    // ---- file_description.source_files ----------------------------------
    {
      const std::vector<SourceFile>& sfs = map.getSourceFiles();
      if (! sfs.empty())
      {
        nlohmann::json sf_array = nlohmann::json::array();
        for (const SourceFile& sf : sfs)
        {
          nlohmann::json sfj = nlohmann::json::object();
          sfj["name"] = std::string(sf.getNameOfFile());
          sfj["location"] = std::string(sf.getPathToFile());
          nlohmann::json sf_params = nlohmann::json::array();
          // SHA-1 checksum
          if (! sf.getChecksum().empty() && sf.getChecksumType() == SourceFile::ChecksumType::SHA1)
          {
            sf_params.push_back({{"accession", "MS:1000569"}, {"name", "SHA-1"}, {"value", std::string(sf.getChecksum())}});
          }
          else if (! sf.getChecksum().empty() && sf.getChecksumType() == SourceFile::ChecksumType::MD5)
          {
            sf_params.push_back({{"accession", "MS:1000568"}, {"name", "MD5"}, {"value", std::string(sf.getChecksum())}});
          }
          // Native ID type
          if (! sf.getNativeIDTypeAccession().empty())
          {
            sf_params.push_back(
              {{"accession", std::string(sf.getNativeIDTypeAccession())}, {"name", std::string(sf.getNativeIDType())}, {"value", nullptr}});
          }
          // File type
          if (! sf.getFileType().empty()) { sf_params.push_back({{"accession", ""}, {"name", std::string(sf.getFileType())}, {"value", nullptr}}); }
          sfj["parameters"] = std::move(sf_params);
          sf_array.push_back(std::move(sfj));
        }
        md["file_description"] = {{"source_files", std::move(sf_array)}};
        has_anything = true;
      }
    }

    // ---- sample_list ----------------------------------------------------
    {
      const Sample& sample = map.getSample();
      if (! sample.getName().empty() || ! sample.getNumber().empty())
      {
        nlohmann::json sj = nlohmann::json::object();
        sj["name"] = std::string(sample.getName());
        sj["id"] = std::string(sample.getNumber());
        md["sample_list"] = nlohmann::json::array({std::move(sj)});
        has_anything = true;
      }
    }

    // ---- software_list --------------------------------------------------
    {
      const Software& sw = map.getInstrument().getSoftware();
      if (! sw.getName().empty())
      {
        nlohmann::json swj = nlohmann::json::object();
        swj["id"] = std::string(sw.getName());
        swj["version"] = std::string(sw.getVersion());
        swj["parameters"] = nlohmann::json::array();
        md["software_list"] = nlohmann::json::array({std::move(swj)});
        has_anything = true;
      }
    }

    if (! has_anything) return nlohmann::json {{"version", "0.9.0"}};
    return md;
  }

  /// Write the spectra_metadata table: a `spectrum` struct column (load path),
  /// plus `precursor` and `selected_ion` struct columns (round-trip path).
  /// The precursor/selected_ion schemas mirror the mzpeak-lib reference writer
  /// exactly so that readPrecursors_() can join them by source_index.
  void writeMetadataTable_(const String& path, const std::vector<MetaRow>& rows, const std::vector<PrecursorOut>& prec_rows)
  {
    // ---- spectrum column ------------------------------------------------
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

    // ---- precursor column -----------------------------------------------
    // Schema:
    //  precursor: struct<
    //    source_index: uint64,
    //    precursor_index: uint64,
    //    precursor_id: large_string,
    //    isolation_window: struct<tgt float, lo float, hi float, parameters large_list<CvParam>>,
    //    activation: struct<parameters large_list<CvParam>>
    //  >
    //
    // CvParam value union: struct<integer int64, float double, string large_string, boolean bool>
    // CvParam: struct<value ..., accession string, name large_string, unit string>
    //
    // We build all builders manually so we can control null/non-null per row.

    auto* mp = arrow::default_memory_pool();

    // ---- CvParam value-union sub-builder (shared for iso_params and act_params) ---
    // We need two independent sets of value-union builders (iso-window params and
    // activation params) because they are children of different LargeListBuilders.

    // Helper lambda to build one CvParam struct array given a list of rows
    // (rows encoded as a vector-of-vectors of ActParam). Each sub-call appends
    // exactly one row per element of prec_rows. The outer struct is always non-null
    // (Append(true)); source_index / precursor_index children are null for MS1 rows.

    // ---- Build value-union type once (reused for both param lists) --------
    auto val_union_type = arrow::struct_({
      arrow::field("integer", arrow::int64(), true),
      arrow::field("float", arrow::float64(), true),
      arrow::field("string", arrow::large_utf8(), true),
      arrow::field("boolean", arrow::boolean(), true),
    });

    auto cvparam_type = arrow::struct_({
      arrow::field("value", val_union_type, true),
      arrow::field("accession", arrow::utf8(), true),
      arrow::field("name", arrow::large_utf8(), true),
      arrow::field("unit", arrow::utf8(), true),
    });

    // isolation_window.parameters type
    auto iso_params_type = arrow::large_list(cvparam_type);

    auto iso_window_type = arrow::struct_({
      arrow::field("MS_1000827_isolation_window_target_mz", arrow::float32(), true),
      arrow::field("MS_1000828_isolation_window_lower_offset", arrow::float32(), true),
      arrow::field("MS_1000829_isolation_window_upper_offset", arrow::float32(), true),
      arrow::field("parameters", iso_params_type, true),
    });

    // activation.parameters type (same CvParam)
    auto act_params_type = arrow::large_list(cvparam_type);
    auto act_type = arrow::struct_({arrow::field("parameters", act_params_type, true)});

    auto prec_type = arrow::struct_({
      arrow::field("source_index", arrow::uint64(), true),
      arrow::field("precursor_index", arrow::uint64(), true),
      arrow::field("precursor_id", arrow::large_utf8(), true),
      arrow::field("isolation_window", iso_window_type, true),
      arrow::field("activation", act_type, true),
    });

    // ---- Precursor column builder ----------------------------------------
    // We build children manually then call StructArray::Make at the end.
    arrow::UInt64Builder prec_si_b(mp);
    arrow::UInt64Builder prec_pi_b(mp);
    arrow::LargeStringBuilder prec_id_b(mp);

    // isolation_window children
    arrow::FloatBuilder iso_tgt_b(mp);
    arrow::FloatBuilder iso_lo_b(mp);
    arrow::FloatBuilder iso_hi_b(mp);

    // isolation_window.parameters — empty list for every row
    // We build a LargeListBuilder whose values are CvParam structs.
    // Since we always emit empty lists (no CV params in iso window),
    // we just accumulate offsets.
    std::vector<int64_t> iso_param_list_offsets; // one per row + sentinel
    iso_param_list_offsets.reserve(prec_rows.size() + 1);
    iso_param_list_offsets.push_back(0); // initial offset

    // No CvParam values for isolation_window.parameters (always empty list).

    // activation.parameters LargeList builder
    // Each row: for MS2 with activation, append the ActParam entries; for MS1, empty list.
    std::vector<int64_t> act_param_list_offsets;
    act_param_list_offsets.reserve(prec_rows.size() + 1);
    act_param_list_offsets.push_back(0);

    // Activation CvParam children arrays (all actual values):
    // value union children
    arrow::Int64Builder act_val_int_b(mp);
    arrow::DoubleBuilder act_val_float_b(mp);
    arrow::LargeStringBuilder act_val_str_b(mp);
    arrow::BooleanBuilder act_val_bool_b(mp);

    arrow::StringBuilder act_acc_b(mp);
    arrow::LargeStringBuilder act_name_b(mp);
    arrow::StringBuilder act_unit_b(mp);

    // outer struct null-bitmaps for precursor and isolation_window and activation
    // We build using StructArray::Make so we need validity bitmaps as uint8 vectors.
    // Simpler: use explicit validity vectors.
    std::vector<bool> prec_outer_valid; // always true
    std::vector<bool> prec_si_valid;
    std::vector<bool> prec_pi_valid;
    std::vector<bool> prec_id_null; // always null
    std::vector<bool> iso_valid;    // always true (outer struct is always non-null)
    std::vector<bool> iso_tgt_valid;
    std::vector<bool> iso_lo_valid;
    std::vector<bool> iso_hi_valid;
    std::vector<bool> iso_param_null; // always emit list (non-null, empty)
    std::vector<bool> act_valid;      // always true

    for (const PrecursorOut& po : prec_rows)
    {
      prec_outer_valid.push_back(true); // outer struct always non-null

      if (po.has_precursor)
      {
        prec_si_valid.push_back(true);
        prec_pi_valid.push_back(true);
        checkArrowStatus_(prec_si_b.Append(po.source_index), "append prec source_index");
        checkArrowStatus_(prec_pi_b.Append(po.precursor_index), "append prec precursor_index");
      }
      else
      {
        prec_si_valid.push_back(false);
        prec_pi_valid.push_back(false);
        checkArrowStatus_(prec_si_b.AppendNull(), "append null prec source_index");
        checkArrowStatus_(prec_pi_b.AppendNull(), "append null prec precursor_index");
      }
      prec_id_null.push_back(false); // always null
      checkArrowStatus_(prec_id_b.AppendNull(), "append null prec id");

      // isolation_window outer struct — always non-null to match reference
      iso_valid.push_back(true);
      if (po.has_isolation)
      {
        iso_tgt_valid.push_back(true);
        iso_lo_valid.push_back(true);
        iso_hi_valid.push_back(true);
        checkArrowStatus_(iso_tgt_b.Append(po.isolation_target), "append iso tgt");
        checkArrowStatus_(iso_lo_b.Append(po.isolation_lower), "append iso lo");
        checkArrowStatus_(iso_hi_b.Append(po.isolation_upper), "append iso hi");
      }
      else
      {
        iso_tgt_valid.push_back(false);
        iso_lo_valid.push_back(false);
        iso_hi_valid.push_back(false);
        checkArrowStatus_(iso_tgt_b.AppendNull(), "append null iso tgt");
        checkArrowStatus_(iso_lo_b.AppendNull(), "append null iso lo");
        checkArrowStatus_(iso_hi_b.AppendNull(), "append null iso hi");
      }
      iso_param_null.push_back(true);                                  // iso params always an empty (non-null) list
      iso_param_list_offsets.push_back(iso_param_list_offsets.back()); // empty list: offset doesn't advance

      // activation outer struct — always non-null
      act_valid.push_back(true);

      // activation.parameters list
      std::size_t n_act = po.activation.size();
      act_param_list_offsets.push_back(act_param_list_offsets.back() + static_cast<int64_t>(n_act));

      for (const PrecursorOut::ActParam& ap : po.activation)
      {
        // value union: if has_float, fill float child, otherwise all null
        if (ap.has_float)
        {
          checkArrowStatus_(act_val_int_b.AppendNull(), "append null act val int");
          checkArrowStatus_(act_val_float_b.Append(ap.float_val), "append act val float");
          checkArrowStatus_(act_val_str_b.AppendNull(), "append null act val str");
          checkArrowStatus_(act_val_bool_b.AppendNull(), "append null act val bool");
        }
        else
        {
          checkArrowStatus_(act_val_int_b.AppendNull(), "append null act val int");
          checkArrowStatus_(act_val_float_b.AppendNull(), "append null act val float");
          checkArrowStatus_(act_val_str_b.AppendNull(), "append null act val str");
          checkArrowStatus_(act_val_bool_b.AppendNull(), "append null act val bool");
        }
        checkArrowStatus_(act_acc_b.Append(std::string(ap.accession)), "append act accession");
        checkArrowStatus_(act_name_b.Append(std::string(ap.name)), "append act name");
        checkArrowStatus_(act_unit_b.AppendNull(), "append null act unit");
      }
    }

    // ---- Finish arrays ---------------------------------------------------
    std::shared_ptr<arrow::Array> prec_si_a, prec_pi_a, prec_id_a;
    checkArrowStatus_(prec_si_b.Finish(&prec_si_a), "finish prec_si");
    checkArrowStatus_(prec_pi_b.Finish(&prec_pi_a), "finish prec_pi");
    checkArrowStatus_(prec_id_b.Finish(&prec_id_a), "finish prec_id");

    std::shared_ptr<arrow::Array> iso_tgt_a, iso_lo_a, iso_hi_a;
    checkArrowStatus_(iso_tgt_b.Finish(&iso_tgt_a), "finish iso_tgt");
    checkArrowStatus_(iso_lo_b.Finish(&iso_lo_a), "finish iso_lo");
    checkArrowStatus_(iso_hi_b.Finish(&iso_hi_a), "finish iso_hi");

    // Build isolation_window.parameters as a LargeListArray of always-empty lists.
    // iso_param_list_offsets has N+1 entries, all equal to 0 (empty lists only).
    // Use MakeArrayOfNull to get a zero-length typed CvParam values array:
    auto empty_cvparam_result = arrow::MakeArrayOfNull(cvparam_type, 0, mp);
    checkArrowStatus_(empty_cvparam_result.status(), "build empty cvparam");
    auto empty_cvparam_a = empty_cvparam_result.ValueOrDie();

    int64_t n_rows = static_cast<int64_t>(prec_rows.size());

    // Manually build the LargeListArray for iso params (always empty lists):
    auto iso_params_offsets_buf = std::make_shared<arrow::Buffer>(reinterpret_cast<const uint8_t*>(iso_param_list_offsets.data()),
                                                                  static_cast<int64_t>(iso_param_list_offsets.size() * sizeof(int64_t)));
    auto iso_params_arr = std::make_shared<arrow::LargeListArray>(iso_params_type, n_rows, iso_params_offsets_buf, empty_cvparam_a);

    // ---- Build activation.parameters arrays ------------------------------
    std::shared_ptr<arrow::Array> act_val_int_a, act_val_float_a, act_val_str_a, act_val_bool_a;
    checkArrowStatus_(act_val_int_b.Finish(&act_val_int_a), "finish act_val_int");
    checkArrowStatus_(act_val_float_b.Finish(&act_val_float_a), "finish act_val_float");
    checkArrowStatus_(act_val_str_b.Finish(&act_val_str_a), "finish act_val_str");
    checkArrowStatus_(act_val_bool_b.Finish(&act_val_bool_a), "finish act_val_bool");

    auto act_val_result
      = arrow::StructArray::Make({act_val_int_a, act_val_float_a, act_val_str_a, act_val_bool_a},
                                 {arrow::field("integer", arrow::int64(), true), arrow::field("float", arrow::float64(), true),
                                  arrow::field("string", arrow::large_utf8(), true), arrow::field("boolean", arrow::boolean(), true)});
    checkArrowStatus_(act_val_result.status(), "build act value union struct");
    auto act_val_a = act_val_result.ValueOrDie();

    std::shared_ptr<arrow::Array> act_acc_a, act_name_a, act_unit_a;
    checkArrowStatus_(act_acc_b.Finish(&act_acc_a), "finish act_acc");
    checkArrowStatus_(act_name_b.Finish(&act_name_a), "finish act_name");
    checkArrowStatus_(act_unit_b.Finish(&act_unit_a), "finish act_unit");

    auto act_cvparam_result = arrow::StructArray::Make({act_val_a, act_acc_a, act_name_a, act_unit_a},
                                                       {arrow::field("value", val_union_type, true), arrow::field("accession", arrow::utf8(), true),
                                                        arrow::field("name", arrow::large_utf8(), true), arrow::field("unit", arrow::utf8(), true)});
    checkArrowStatus_(act_cvparam_result.status(), "build act cvparam struct");
    auto act_cvparam_a = act_cvparam_result.ValueOrDie();

    // Build the activation params LargeListArray from offsets
    auto act_params_offsets_buf = std::make_shared<arrow::Buffer>(reinterpret_cast<const uint8_t*>(act_param_list_offsets.data()),
                                                                  static_cast<int64_t>(act_param_list_offsets.size() * sizeof(int64_t)));
    auto act_params_arr = std::make_shared<arrow::LargeListArray>(act_params_type, n_rows, act_params_offsets_buf, act_cvparam_a);

    // ---- Build isolation_window struct -----------------------------------
    auto iso_result = arrow::StructArray::Make({iso_tgt_a, iso_lo_a, iso_hi_a, iso_params_arr},
                                               {arrow::field("MS_1000827_isolation_window_target_mz", arrow::float32(), true),
                                                arrow::field("MS_1000828_isolation_window_lower_offset", arrow::float32(), true),
                                                arrow::field("MS_1000829_isolation_window_upper_offset", arrow::float32(), true),
                                                arrow::field("parameters", iso_params_type, true)});
    checkArrowStatus_(iso_result.status(), "build isolation_window struct");
    auto iso_a = iso_result.ValueOrDie();

    // ---- Build activation struct -----------------------------------------
    auto act_result = arrow::StructArray::Make({act_params_arr}, {arrow::field("parameters", act_params_type, true)});
    checkArrowStatus_(act_result.status(), "build activation struct");
    auto act_a = act_result.ValueOrDie();

    // ---- Build precursor outer struct ------------------------------------
    auto prec_result
      = arrow::StructArray::Make({prec_si_a, prec_pi_a, prec_id_a, iso_a, act_a},
                                 {arrow::field("source_index", arrow::uint64(), true), arrow::field("precursor_index", arrow::uint64(), true),
                                  arrow::field("precursor_id", arrow::large_utf8(), true), arrow::field("isolation_window", iso_window_type, true),
                                  arrow::field("activation", act_type, true)});
    checkArrowStatus_(prec_result.status(), "build precursor struct");
    auto prec_array = prec_result.ValueOrDie();

    // ---- selected_ion column --------------------------------------------
    // Schema: struct<source_index uint64, precursor_index uint64,
    //   MS_1000744_selected_ion_mz_unit_MS_1000040 double,
    //   MS_1000041_charge_state int32,
    //   MS_1000042_intensity_unit_MS_1000131 float>
    arrow::UInt64Builder sion_si_b(mp);
    arrow::UInt64Builder sion_pi_b(mp);
    arrow::DoubleBuilder sion_mz_b(mp);
    arrow::Int32Builder sion_charge_b(mp);
    arrow::FloatBuilder sion_int_b(mp);

    std::vector<bool> sion_si_valid;
    std::vector<bool> sion_pi_valid;
    std::vector<bool> sion_mz_valid;
    std::vector<bool> sion_charge_valid;
    std::vector<bool> sion_int_valid;

    for (const PrecursorOut& po : prec_rows)
    {
      if (po.has_precursor && po.has_selected_ion)
      {
        sion_si_valid.push_back(true);
        sion_pi_valid.push_back(true);
        sion_mz_valid.push_back(true);
        checkArrowStatus_(sion_si_b.Append(po.source_index), "append sion si");
        checkArrowStatus_(sion_pi_b.Append(po.precursor_index), "append sion pi");
        checkArrowStatus_(sion_mz_b.Append(po.selected_ion_mz), "append sion mz");
      }
      else
      {
        sion_si_valid.push_back(false);
        sion_pi_valid.push_back(false);
        sion_mz_valid.push_back(false);
        checkArrowStatus_(sion_si_b.AppendNull(), "append null sion si");
        checkArrowStatus_(sion_pi_b.AppendNull(), "append null sion pi");
        checkArrowStatus_(sion_mz_b.AppendNull(), "append null sion mz");
      }

      if (po.has_precursor && po.has_charge)
      {
        sion_charge_valid.push_back(true);
        checkArrowStatus_(sion_charge_b.Append(po.charge), "append sion charge");
      }
      else
      {
        sion_charge_valid.push_back(false);
        checkArrowStatus_(sion_charge_b.AppendNull(), "append null sion charge");
      }

      if (po.has_precursor && po.has_intensity)
      {
        sion_int_valid.push_back(true);
        checkArrowStatus_(sion_int_b.Append(po.intensity), "append sion intensity");
      }
      else
      {
        sion_int_valid.push_back(false);
        checkArrowStatus_(sion_int_b.AppendNull(), "append null sion intensity");
      }
    }

    std::shared_ptr<arrow::Array> sion_si_a, sion_pi_a, sion_mz_a, sion_charge_a, sion_int_a;
    checkArrowStatus_(sion_si_b.Finish(&sion_si_a), "finish sion_si");
    checkArrowStatus_(sion_pi_b.Finish(&sion_pi_a), "finish sion_pi");
    checkArrowStatus_(sion_mz_b.Finish(&sion_mz_a), "finish sion_mz");
    checkArrowStatus_(sion_charge_b.Finish(&sion_charge_a), "finish sion_charge");
    checkArrowStatus_(sion_int_b.Finish(&sion_int_a), "finish sion_int");

    auto sion_result = arrow::StructArray::Make(
      {sion_si_a, sion_pi_a, sion_mz_a, sion_charge_a, sion_int_a},
      {arrow::field("source_index", arrow::uint64(), true), arrow::field("precursor_index", arrow::uint64(), true),
       arrow::field("MS_1000744_selected_ion_mz_unit_MS_1000040", arrow::float64(), true),
       arrow::field("MS_1000041_charge_state", arrow::int32(), true), arrow::field("MS_1000042_intensity_unit_MS_1000131", arrow::float32(), true)});
    checkArrowStatus_(sion_result.status(), "build selected_ion struct");
    auto sion_array = sion_result.ValueOrDie();

    // ---- Assemble table with 3 columns ----------------------------------
    auto schema = arrow::schema({
      arrow::field("spectrum", spectrum_array->type(), /*nullable=*/true),
      arrow::field("precursor", prec_array->type(), /*nullable=*/true),
      arrow::field("selected_ion", sion_array->type(), /*nullable=*/true),
    });
    auto table = arrow::Table::Make(schema, {spectrum_array, prec_array, sion_array});

    writeTableWithProps_(path, table, {});
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
  std::vector<PrecursorOut> prec_rows;
  meta_rows.reserve(map.size());
  prec_rows.reserve(map.size());

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

    // Precursor data for this spectrum
    prec_rows.push_back(buildPrecursorOut_(spec, static_cast<std::uint64_t>(i), 0));
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

  // The spectrum_array_index JSON is required by the Rust reference reader to
  // locate the column layout (prefix "point", mz/intensity array type info).
  // It mirrors the exact structure emitted by the mzpeak-lib Rust writer.
  static const std::string kSpectrumArrayIndex
    = R"({"prefix":"point","entries":[{"context":"spectrum","path":"point.mz","data_type":"MS:1000523","array_type":"MS:1000514","array_name":"m/z array","unit":"MS:1000040","buffer_format":"point","transform":"MS:1003901","data_processing_id":null,"buffer_priority":"primary","sorting_rank":0},{"context":"spectrum","path":"point.intensity","data_type":"MS:1000521","array_type":"MS:1000515","array_name":"intensity array","unit":"MS:1000131","buffer_format":"point","transform":"MS:1003902","data_processing_id":null,"buffer_priority":"primary","sorting_rank":null}]})";

  // File-kv for point tables: spectrum_count and spectrum_data_point_count.
  const std::string n_spectra = std::to_string(map.size());
  const std::map<std::string, std::string> data_kv {
    {"spectrum_count", n_spectra}, {"spectrum_data_point_count", std::to_string(data.mz.size())}, {"spectrum_array_index", kSpectrumArrayIndex}};
  const std::map<std::string, std::string> peaks_kv {
    {"spectrum_count", n_spectra}, {"spectrum_data_point_count", std::to_string(peaks.mz.size())}, {"spectrum_array_index", kSpectrumArrayIndex}};

  // The data table is always present (may be empty); peaks only when a centroid
  // spectrum exists, mirroring the mzpeak-lib writer.
  writePointTable_(data_path, data, data_kv);
  if (any_centroid) writePointTable_(peaks_path, peaks, peaks_kv);
  writeMetadataTable_(meta_path, meta_rows, prec_rows);

  // ------------------------------------------------------------------
  // 3. Write mzpeak_index.json (files[] + run-level metadata). The load
  //    path matches members by data_kind, so the names + kinds must be exact.
  // ------------------------------------------------------------------
  nlohmann::json idx;
  nlohmann::json files = nlohmann::json::array();
  files.push_back({{"name", "spectra_data.parquet"}, {"entity_type", "spectrum"}, {"data_kind", "data arrays"}});
  if (any_centroid) { files.push_back({{"name", "spectra_peaks.parquet"}, {"entity_type", "spectrum"}, {"data_kind", "peaks"}}); }
  files.push_back({{"name", "spectra_metadata.parquet"}, {"entity_type", "spectrum"}, {"data_kind", "metadata"}});
  idx["files"] = std::move(files);

  // Build run-level metadata from ExperimentalSettings (best-effort)
  nlohmann::json run_md = buildRunMetadataJson_(map);
  run_md["version"] = "0.9.0";
  idx["metadata"] = std::move(run_md);

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
