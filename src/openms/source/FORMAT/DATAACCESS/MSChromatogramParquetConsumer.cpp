// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/ZlibCompression.h>

#include <exception>
#include <limits>
#include <memory>
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace OpenMS
{
  namespace
  {
    void appendOrThrow_(const arrow::Status& status, const char* column)
    {
      if (!status.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      String("Failed to append value for ") + column, status.ToString());
      }
    }

    void reserveOrThrow_(const arrow::Status& status, const char* column)
    {
      if (!status.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      String("Failed to reserve capacity for ") + column, status.ToString());
      }
    }

    template <typename Builder>
    std::shared_ptr<arrow::Array> finishArray_(Builder& builder, const char* name)
    {
      std::shared_ptr<arrow::Array> array;
      auto status = builder.Finish(&array);
      if (!status.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      String("Failed to finish array for ") + name, status.ToString());
      }
      return array;
    }

    void appendOptionalInt_(arrow::Int64Builder& builder, bool has_value, int64_t value, const char* column)
    {
      if (!has_value)
      {
        appendOrThrow_(builder.AppendNull(), column);
      }
      else
      {
        appendOrThrow_(builder.Append(value), column);
      }
    }

    void appendOptionalString_(arrow::StringBuilder& builder, const String& value, const char* column)
    {
      if (value.empty())
      {
        appendOrThrow_(builder.AppendNull(), column);
      }
      else
      {
        appendOrThrow_(builder.Append(std::string(value)), column);
      }
    }

    void appendBinary_(arrow::BinaryBuilder& builder, const String& value, const char* column)
    {
      if (value.size() > static_cast<size_t>(std::numeric_limits<int32_t>::max()))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      String("Binary value for column ") + column + " exceeds Arrow int32 size limit (2 GB).",
                                      String(value.size()));
      }
      appendOrThrow_(builder.Append(reinterpret_cast<const uint8_t*>(value.c_str()),
                                    static_cast<int32_t>(value.size())), column);
    }

    struct CompoundInfo
    {
      int64_t precursor_id = 0;
      String modified_sequence;
      int64_t precursor_charge = 0;
      int64_t precursor_decoy = 0;
    };

    struct TransitionInfo
    {
      int64_t transition_id = 0;
      int64_t precursor_id = 0;
      String modified_sequence;
      int64_t precursor_charge = 0;
      int64_t product_charge = 0;
      int64_t detecting_transition = 0;
      int64_t precursor_decoy = 0;
      int64_t product_decoy = 0;
      int64_t transition_ordinal = 0;
      String transition_type;
      String annotation;
    };

    int64_t parseOrAssignId_(const String& text, int64_t& next_id, std::unordered_set<int64_t>& used_ids)
    {
      try
      {
        int64_t value = text.toInt64();
        if (used_ids.insert(value).second)
        {
          if (value >= next_id)
          {
            next_id = value + 1;
          }
          return value;
        }
      }
      catch (Exception::ConversionError&)
      {
        // fall through to auto-assigned ID
      }
      while (used_ids.count(next_id))
      {
        ++next_id;
      }
      used_ids.insert(next_id);
      return next_id++;
    }

    String buildAnnotation_(const String& transition_type, int64_t ordinal, int64_t charge)
    {
      if (transition_type.empty() || ordinal < 0)
      {
        return "";
      }
      String annotation = transition_type + String(ordinal);
      if (charge > 0)
      {
        annotation += "^" + String(charge);
      }
      return annotation;
    }

    String buildPrecursorAnnotation_(const String& native_id)
    {
      const String tag = "_Precursor_i";
      const Size pos = native_id.rfind(tag);
      if (pos != String::npos)
      {
        return native_id.substr(pos + 1);
      }
      return "";
    }
  } // namespace

  class MSChromatogramParquetConsumerImpl
  {
  public:
    MSChromatogramParquetConsumerImpl(const String& filename,
                                      UInt64 run_id,
                                      const String& source_file,
                                      const OpenSwath::LightTargetedExperiment& transition_exp) :
      filename_(filename),
      run_id_(run_id),
      source_file_(source_file)
    {
      buildTransitionMaps_(transition_exp);
      // Build the invariant schema once; reuse on every flushBuffers_() call.
      schema_ = XICSchema::schema();
    }

    ~MSChromatogramParquetConsumerImpl()
    {
      try
      {
        finalize();
      }
      catch (const Exception::BaseException& e)
      {
        OPENMS_LOG_ERROR << "Failed to write chromatogram parquet file '" << filename_
                         << "': " << e.what() << "\n";
      }
      catch (const std::exception& e)
      {
        OPENMS_LOG_ERROR << "Failed to write chromatogram parquet file '" << filename_
                         << "': " << e.what() << "\n";
      }
      catch (...)
      {
        OPENMS_LOG_ERROR << "Failed to write chromatogram parquet file '" << filename_
                         << "': unknown exception.\n";
      }
    }

    void consumeSpectrum(MSChromatogramParquetConsumer::SpectrumType& s)
    {
      s.clear(false);
    }

    void consumeChromatogram(MSChromatogramParquetConsumer::ChromatogramType& c)
    {
      const String native_id = c.getNativeID();
      const bool is_precursor = native_id.hasSubstring("_Precursor_i");
      const int64_t ms_level = is_precursor ? 1 : 2;

      appendOrThrow_(run_id_builder_.Append(static_cast<int64_t>(run_id_)), "RUN_ID");
      appendOptionalString_(source_file_builder_, source_file_, "SOURCE_FILE");
      appendOrThrow_(ms_level_builder_.Append(ms_level), "MS_LEVEL");
      if (is_precursor)
      {
        const String precursor_annotation = buildPrecursorAnnotation_(native_id);
        const String group_id = OpenSwathHelper::computeTransitionGroupId(native_id);
        auto comp_it = compound_info_.find(group_id);
        if (comp_it != compound_info_.end())
        {
          const CompoundInfo& comp = comp_it->second;
          appendOptionalInt_(precursor_id_builder_, true, comp.precursor_id, "PRECURSOR_ID");
          appendOptionalString_(modified_sequence_builder_, comp.modified_sequence, "MODIFIED_SEQUENCE");
          appendOptionalInt_(precursor_charge_builder_, comp.precursor_charge > 0, comp.precursor_charge, "PRECURSOR_CHARGE");
          appendOptionalInt_(precursor_decoy_builder_, true, comp.precursor_decoy, "PRECURSOR_DECOY");
        }
        else
        {
          // If we have an extracted ion chromatogram to write out to disk, we likely would want to know what precursor the XIC belongs to. Otherwise we could end up with a parquet file with a log of XICs but without any meta data to identify them.
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        String("Chromatogram native ID '") + native_id + "' looks like a precursor but no matching compound metadata entry was found. Please ensure the native ID is correct and that the transition library contains metadata for this precursor.", String());
        }

        appendOrThrow_(transition_id_builder_.AppendNull(), "TRANSITION_ID");
        appendOrThrow_(product_charge_builder_.AppendNull(), "PRODUCT_CHARGE");
        appendOrThrow_(detecting_transition_builder_.AppendNull(), "DETECTING_TRANSITION");
        appendOrThrow_(product_decoy_builder_.AppendNull(), "PRODUCT_DECOY");
        appendOrThrow_(transition_ordinal_builder_.AppendNull(), "TRANSITION_ORDINAL");
        appendOrThrow_(transition_type_builder_.AppendNull(), "TRANSITION_TYPE");
        appendOptionalString_(annotation_builder_, precursor_annotation, "ANNOTATION");
      }
      else
      {
        auto tr_it = transition_info_.find(native_id);
        if (tr_it != transition_info_.end())
        {
          const TransitionInfo& tr = tr_it->second;
          appendOptionalInt_(precursor_id_builder_, true, tr.precursor_id, "PRECURSOR_ID");
          appendOptionalInt_(transition_id_builder_, true, tr.transition_id, "TRANSITION_ID");
          appendOptionalString_(modified_sequence_builder_, tr.modified_sequence, "MODIFIED_SEQUENCE");
          appendOptionalInt_(precursor_charge_builder_, tr.precursor_charge > 0, tr.precursor_charge, "PRECURSOR_CHARGE");
          appendOptionalInt_(product_charge_builder_, tr.product_charge > 0, tr.product_charge, "PRODUCT_CHARGE");
          appendOptionalInt_(detecting_transition_builder_, true, tr.detecting_transition, "DETECTING_TRANSITION");
          appendOptionalInt_(precursor_decoy_builder_, true, tr.precursor_decoy, "PRECURSOR_DECOY");
          appendOptionalInt_(product_decoy_builder_, true, tr.product_decoy, "PRODUCT_DECOY");
          appendOptionalInt_(transition_ordinal_builder_, tr.transition_ordinal >= 0, tr.transition_ordinal, "TRANSITION_ORDINAL");
          appendOptionalString_(transition_type_builder_, tr.transition_type, "TRANSITION_TYPE");
          appendOptionalString_(annotation_builder_, tr.annotation, "ANNOTATION");
        }
        else
        {
          // Same comment as for the precursor case applies here: if we have a chromatogram to write out to disk, we likely would want to know what transition it belongs to. Otherwise we could end up with a parquet file with a log of chromatograms but without any meta data to identify them.
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        String("Chromatogram native ID '") + native_id + "' does not have a matching transition metadata entry. Please ensure the native ID is correct and that the transition library contains metadata for this transition.", String());
        }
      }

      encodeChromatogram_(c);
      c.clear(false);
    }

    void setExpectedSize(Size, Size expectedChromatograms)
    {
      if (expectedChromatograms > 0)
      {
        reserveOrThrow_(run_id_builder_.Reserve(expectedChromatograms), "RUN_ID");
        reserveOrThrow_(source_file_builder_.Reserve(expectedChromatograms), "SOURCE_FILE");
        reserveOrThrow_(ms_level_builder_.Reserve(expectedChromatograms), "MS_LEVEL");
        reserveOrThrow_(precursor_id_builder_.Reserve(expectedChromatograms), "PRECURSOR_ID");
        reserveOrThrow_(transition_id_builder_.Reserve(expectedChromatograms), "TRANSITION_ID");
        reserveOrThrow_(modified_sequence_builder_.Reserve(expectedChromatograms), "MODIFIED_SEQUENCE");
        reserveOrThrow_(precursor_charge_builder_.Reserve(expectedChromatograms), "PRECURSOR_CHARGE");
        reserveOrThrow_(product_charge_builder_.Reserve(expectedChromatograms), "PRODUCT_CHARGE");
        reserveOrThrow_(detecting_transition_builder_.Reserve(expectedChromatograms), "DETECTING_TRANSITION");
        reserveOrThrow_(precursor_decoy_builder_.Reserve(expectedChromatograms), "PRECURSOR_DECOY");
        reserveOrThrow_(product_decoy_builder_.Reserve(expectedChromatograms), "PRODUCT_DECOY");
        reserveOrThrow_(transition_ordinal_builder_.Reserve(expectedChromatograms), "TRANSITION_ORDINAL");
        reserveOrThrow_(transition_type_builder_.Reserve(expectedChromatograms), "TRANSITION_TYPE");
        reserveOrThrow_(annotation_builder_.Reserve(expectedChromatograms), "ANNOTATION");
        reserveOrThrow_(rt_data_builder_.Reserve(expectedChromatograms), "RT_DATA");
        reserveOrThrow_(intensity_data_builder_.Reserve(expectedChromatograms), "INTENSITY_DATA");
        reserveOrThrow_(rt_compression_builder_.Reserve(expectedChromatograms), "RT_COMPRESSION");
        reserveOrThrow_(intensity_compression_builder_.Reserve(expectedChromatograms), "INTENSITY_COMPRESSION");
      }
    }

    void setExperimentalSettings(const ExperimentalSettings&)
    {
    }

    void finalize()
    {
      if (wrote_)
      {
        return;
      }
      write_();
    }

  private:
    String filename_;
    UInt64 run_id_{0};
    String source_file_;
    bool wrote_{false};
    void buildTransitionMaps_(const OpenSwath::LightTargetedExperiment& transition_exp)
    {
      // Build lookup tables for precursor- and transition-level metadata.

      std::unordered_map<String, int64_t> precursor_ids;
      precursor_ids.reserve(transition_exp.getCompounds().size());

      std::unordered_map<String, int64_t> precursor_decoy;
      precursor_decoy.reserve(transition_exp.getCompounds().size());

      for (const auto& compound : transition_exp.getCompounds())
      {
        const String compound_id = compound.id;
        const int64_t precursor_id = parseOrAssignId_(compound_id, next_precursor_id_, used_precursor_ids_);
        precursor_ids[compound_id] = precursor_id;
        precursor_decoy[compound_id] = 0;

        CompoundInfo info;
        info.precursor_id = precursor_id;
        info.modified_sequence = compound.sequence;
        info.precursor_charge = compound.charge;
        info.precursor_decoy = 0;
        compound_info_.emplace(compound_id, std::move(info));
      }

      for (const auto& transition : transition_exp.getTransitions())
      {
        const String transition_name = transition.transition_name;
        const String peptide_ref = transition.peptide_ref;

        int64_t transition_id = 0;
        auto it = transition_ids_.find(transition_name);
        if (it == transition_ids_.end())
        {
          transition_id = parseOrAssignId_(transition_name, next_transition_id_, used_transition_ids_);
          transition_ids_[transition_name] = transition_id;
        }
        else
        {
          transition_id = it->second;
        }

        auto precursor_it = precursor_ids.find(peptide_ref);
        if (precursor_it == precursor_ids.end())
        {
          const int64_t precursor_id = parseOrAssignId_(peptide_ref, next_precursor_id_, used_precursor_ids_);
          precursor_ids[peptide_ref] = precursor_id;
          precursor_decoy[peptide_ref] = 0;

          CompoundInfo info;
          info.precursor_id = precursor_id;
          info.modified_sequence = "";
          info.precursor_charge = 0;
          info.precursor_decoy = 0;
          compound_info_.emplace(peptide_ref, std::move(info));
          precursor_it = precursor_ids.find(peptide_ref);
        }

        const bool is_decoy = transition.getDecoy();
        if (is_decoy)
        {
          precursor_decoy[peptide_ref] = 1;
          auto comp_it = compound_info_.find(peptide_ref);
          if (comp_it != compound_info_.end())
          {
            comp_it->second.precursor_decoy = 1;
          }
        }

        TransitionInfo info;
        info.transition_id = transition_id;
        info.precursor_id = precursor_it->second;
        auto comp_info_it = compound_info_.find(peptide_ref);
        if (comp_info_it != compound_info_.end())
        {
          info.modified_sequence = comp_info_it->second.modified_sequence;
          info.precursor_charge = comp_info_it->second.precursor_charge;
          info.precursor_decoy = comp_info_it->second.precursor_decoy;
        }
        info.product_charge = transition.fragment_charge;
        info.detecting_transition = transition.isDetectingTransition() ? 1 : 0;
        info.product_decoy = is_decoy ? 1 : 0;
        info.transition_ordinal = transition.fragment_nr;
        info.transition_type = transition.getFragmentType();
        info.annotation = buildAnnotation_(info.transition_type, info.transition_ordinal, info.product_charge);

        transition_info_.emplace(transition_name, std::move(info));
      }
    }

    void encodeChromatogram_(const MSChromatogramParquetConsumer::ChromatogramType& c)
    {
      // Marshal chromatogram data into RT/intensity arrays for encoding.
      std::vector<double> rt_data;
      std::vector<double> intensity_data;
      rt_data.reserve(c.size());
      intensity_data.reserve(c.size());
      for (const auto& peak : c)
      {
        rt_data.push_back(peak.getRT());
        intensity_data.push_back(peak.getIntensity());
      }

      const bool use_lossy_compression = true;
      const int64_t rt_compression = use_lossy_compression ? 5 : 1;
      const int64_t intensity_compression = use_lossy_compression ? 6 : 1;

      String rt_encoded;
      String int_encoded;
      if (rt_data.empty())
      {
        appendBinary_(rt_data_builder_, rt_encoded, "RT_DATA");
        appendBinary_(intensity_data_builder_, int_encoded, "INTENSITY_DATA");
        appendOrThrow_(rt_compression_builder_.Append(rt_compression), "RT_COMPRESSION");
        appendOrThrow_(intensity_compression_builder_.Append(intensity_compression), "INTENSITY_COMPRESSION");
        return;
      }
      if (use_lossy_compression)
      {
        // MSNumpress encodes doubles to byte strings, then zlib compresses.
        MSNumpressCoder::NumpressConfig npconfig_mz;
        npconfig_mz.estimate_fixed_point = true;
        npconfig_mz.numpressErrorTolerance = -1.0;
        npconfig_mz.setCompression("linear");
        npconfig_mz.linear_fp_mass_acc = 0.05;

        MSNumpressCoder::NumpressConfig npconfig_int;
        npconfig_int.estimate_fixed_point = true;
        npconfig_int.numpressErrorTolerance = -1.0;
        npconfig_int.setCompression("slof");

        String rt_uncompressed;
        MSNumpressCoder().encodeNPRaw(rt_data, rt_uncompressed, npconfig_mz);
        ZlibCompression::compressString(rt_uncompressed, rt_encoded);

        String int_uncompressed;
        MSNumpressCoder().encodeNPRaw(intensity_data, int_uncompressed, npconfig_int);
        ZlibCompression::compressString(int_uncompressed, int_encoded);
      }
      else
      {
        // Raw doubles are zlib-compressed without MSNumpress.
        std::string rt_bytes(reinterpret_cast<const char*>(rt_data.data()), rt_data.size() * sizeof(double));
        ZlibCompression::compressString(rt_bytes, rt_encoded);
        std::string int_bytes(reinterpret_cast<const char*>(intensity_data.data()), intensity_data.size() * sizeof(double));
        ZlibCompression::compressString(int_bytes, int_encoded);
      }

      appendBinary_(rt_data_builder_, rt_encoded, "RT_DATA");
      appendBinary_(intensity_data_builder_, int_encoded, "INTENSITY_DATA");
      appendOrThrow_(rt_compression_builder_.Append(rt_compression), "RT_COMPRESSION");
      appendOrThrow_(intensity_compression_builder_.Append(intensity_compression), "INTENSITY_COMPRESSION");
      
      // Track accumulated binary size and row count. Flush BEFORE the next
      // append if the thresholds would be exceeded, so we never exceed Arrow's
      // int32 capacity limits within a single batch.
      accumulated_rows_++;
      accumulated_binary_bytes_ += rt_encoded.size();
      accumulated_binary_bytes_ += int_encoded.size();

      if (accumulated_rows_ >= flush_row_limit_ || accumulated_binary_bytes_ >= flush_byte_limit_)
      {
        flushBuffers_();
      }
    }

    void write_()
    {
      // Finalize builders into arrays and assemble the table and write any
      // remaining data.
      flushBuffers_();

      // Close writer and outfile so Parquet metadata is written.
      // Null out each pointer immediately after closing so a subsequent call
      // from the destructor (if the first call threw) cannot double-close.
      if (parquet_writer_)
      {
        auto s = parquet_writer_->Close();
        parquet_writer_.reset(); // prevent double-close even if outfile_ close throws
        if (!s.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to close parquet writer", s.ToString());
        }
      }
      // Mark as written before closing the file so the destructor's retry path
      // (triggered if outfile_->Close() throws) does not re-enter write_().
      wrote_ = true;
      if (outfile_)
      {
        auto s2 = outfile_->Close();
        outfile_.reset(); // prevent double-close
        if (!s2.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to close parquet outfile", s2.ToString());
        }
      }
      return;
    }

    void flushBuffers_()
    {
      // If there is nothing to flush, we usually return early.
      // However, if no outfile/writer has been created yet (i.e. no
      // flush has ever happened) we still want to create an empty
      // parquet file so consumers that expect the file to exist
      // (even if it contains zero rows) can open it. Only return
      // when there are truly no rows to flush and a writer/outfile
      // already exists.
      if (accumulated_rows_ == 0)
      {
        if (outfile_ || parquet_writer_)
        {
          return;
        }
        // otherwise, fall through and create an empty file below
        OPENMS_LOG_WARN << "No chromatograms written to '" << filename_
                         << "': creating empty Parquet file with zero rows.\n";
      }

      // Build Arrow schema for the chromatogram table.
      auto schema = schema_;

      // Finalize builders into arrays and assemble the table.
      // Finish arrays from the current builders, then immediately reset the
      // builders so they are fresh for the next batch regardless of whether
      // WriteTable succeeds or throws.
      auto arr_run_id         = finishArray_(run_id_builder_,                "RUN_ID");
      auto arr_source_file    = finishArray_(source_file_builder_,           "SOURCE_FILE");
      auto arr_ms_level       = finishArray_(ms_level_builder_,              "MS_LEVEL");
      auto arr_prec_id        = finishArray_(precursor_id_builder_,          "PRECURSOR_ID");
      auto arr_tr_id          = finishArray_(transition_id_builder_,         "TRANSITION_ID");
      auto arr_mod_seq        = finishArray_(modified_sequence_builder_,     "MODIFIED_SEQUENCE");
      auto arr_prec_charge    = finishArray_(precursor_charge_builder_,      "PRECURSOR_CHARGE");
      auto arr_prod_charge    = finishArray_(product_charge_builder_,        "PRODUCT_CHARGE");
      auto arr_detecting      = finishArray_(detecting_transition_builder_,  "DETECTING_TRANSITION");
      auto arr_prec_decoy     = finishArray_(precursor_decoy_builder_,       "PRECURSOR_DECOY");
      auto arr_prod_decoy     = finishArray_(product_decoy_builder_,         "PRODUCT_DECOY");
      auto arr_tr_ordinal     = finishArray_(transition_ordinal_builder_,    "TRANSITION_ORDINAL");
      auto arr_tr_type        = finishArray_(transition_type_builder_,       "TRANSITION_TYPE");
      auto arr_annotation     = finishArray_(annotation_builder_,            "ANNOTATION");
      auto arr_rt_data        = finishArray_(rt_data_builder_,               "RT_DATA");
      auto arr_int_data       = finishArray_(intensity_data_builder_,        "INTENSITY_DATA");
      auto arr_rt_comp        = finishArray_(rt_compression_builder_,        "RT_COMPRESSION");
      auto arr_int_comp       = finishArray_(intensity_compression_builder_, "INTENSITY_COMPRESSION");

      // Reset builders and counters immediately — before WriteTable — so the
      // destructor's retry path always finds builders in a usable state.
      run_id_builder_               = arrow::Int64Builder();
      source_file_builder_          = arrow::StringBuilder();
      ms_level_builder_             = arrow::Int64Builder();
      precursor_id_builder_         = arrow::Int64Builder();
      transition_id_builder_        = arrow::Int64Builder();
      modified_sequence_builder_    = arrow::StringBuilder();
      precursor_charge_builder_     = arrow::Int64Builder();
      product_charge_builder_       = arrow::Int64Builder();
      detecting_transition_builder_ = arrow::Int64Builder();
      precursor_decoy_builder_      = arrow::Int64Builder();
      product_decoy_builder_        = arrow::Int64Builder();
      transition_ordinal_builder_   = arrow::Int64Builder();
      transition_type_builder_      = arrow::StringBuilder();
      annotation_builder_           = arrow::StringBuilder();
      rt_data_builder_              = arrow::BinaryBuilder();
      intensity_data_builder_       = arrow::BinaryBuilder();
      rt_compression_builder_       = arrow::Int64Builder();
      intensity_compression_builder_= arrow::Int64Builder();
      accumulated_rows_             = 0;
      accumulated_binary_bytes_     = 0;

      auto table = arrow::Table::Make(schema, {
        arr_run_id, arr_source_file, arr_ms_level,
        arr_prec_id, arr_tr_id, arr_mod_seq,
        arr_prec_charge, arr_prod_charge, arr_detecting,
        arr_prec_decoy, arr_prod_decoy, arr_tr_ordinal,
        arr_tr_type, arr_annotation,
        arr_rt_data, arr_int_data, arr_rt_comp, arr_int_comp
      });

      // Validate table against registry schema
      auto validation = ArrowSchemaValidation::validate(table, XICSchema::schema());
      if (!validation.valid)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "XIC table schema validation failed: " + validation.toString(), "");
      }

      // Open output file on first flush and prepare writer properties.
      if (!outfile_)
      {
        auto outfile_result = arrow::io::FileOutputStream::Open(std::string(filename_));
        if (!outfile_result.ok())
        {
          throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename_);
        }
        outfile_ = std::move(outfile_result).ValueOrDie();
        parquet::WriterProperties::Builder builder;
        builder.compression(parquet::Compression::ZSTD);
        builder.compression_level(11);
        props_ = builder.build();
      }

      // Write the table as row-group(s). Use parquet::arrow::FileWriter for
      // incremental writes; create the writer on first flush and reuse it for
      // subsequent flushes. FileWriter writes metadata at Close(), so we keep
      // it open until finalize().
      if (!parquet_writer_)
      {
        // Create a FileWriter (returns arrow::Result)
        auto writer_result = parquet::arrow::FileWriter::Open(*schema,
                             arrow::default_memory_pool(),
                             outfile_,
                             props_,
                             parquet::default_arrow_writer_properties());
        if (!writer_result.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to create parquet FileWriter", writer_result.status().ToString());
        }
        parquet_writer_ = std::move(writer_result).ValueOrDie();
      }
      auto status = parquet_writer_->WriteTable(*table, 1024);
      if (!status.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to write chromatogram parquet table", status.ToString());
      }
    }

    std::unordered_map<String, CompoundInfo> compound_info_;
    std::unordered_map<String, TransitionInfo> transition_info_;
    std::unordered_map<String, int64_t> transition_ids_;
    // Auto-assignment state for when no transition experiment is provided.
    int64_t next_precursor_id_ = 1;
    int64_t next_transition_id_ = 1;
    std::unordered_set<int64_t> used_precursor_ids_;
    std::unordered_set<int64_t> used_transition_ids_;

    arrow::Int64Builder run_id_builder_;
    arrow::StringBuilder source_file_builder_;
    arrow::Int64Builder ms_level_builder_;
    arrow::Int64Builder precursor_id_builder_;
    arrow::Int64Builder transition_id_builder_;
    arrow::StringBuilder modified_sequence_builder_;
    arrow::Int64Builder precursor_charge_builder_;
    arrow::Int64Builder product_charge_builder_;
    arrow::Int64Builder detecting_transition_builder_;
    arrow::Int64Builder precursor_decoy_builder_;
    arrow::Int64Builder product_decoy_builder_;
    arrow::Int64Builder transition_ordinal_builder_;
    arrow::StringBuilder transition_type_builder_;
    arrow::StringBuilder annotation_builder_;
    arrow::BinaryBuilder rt_data_builder_;
    arrow::BinaryBuilder intensity_data_builder_;
    arrow::Int64Builder rt_compression_builder_;
    arrow::Int64Builder intensity_compression_builder_;
    // Streaming/flush state to avoid building giant arrays that exceed Arrow's
    // int32 capacity limits.
    size_t accumulated_rows_ = 0;
    size_t accumulated_binary_bytes_ = 0;
    // Flush when this many chromatograms are buffered
    const size_t flush_row_limit_ = 1000;
    // Flush when this many binary bytes are buffered (safety margin below 2^31)
    const size_t flush_byte_limit_ = 1900000000ULL; // ~1.9 GB: well below Arrow's int32 limit
    std::shared_ptr<arrow::io::FileOutputStream> outfile_;
    std::shared_ptr<parquet::WriterProperties> props_;
    std::unique_ptr<parquet::arrow::FileWriter> parquet_writer_;
    std::shared_ptr<arrow::Schema> schema_;
  };

  MSChromatogramParquetConsumer::MSChromatogramParquetConsumer(const String& filename,
                                                               UInt64 run_id,
                                                               const String& source_file,
                                                               const OpenSwath::LightTargetedExperiment& transition_exp)
  {
    impl_ = std::make_unique<MSChromatogramParquetConsumerImpl>(filename, run_id, source_file, transition_exp);
  }

  MSChromatogramParquetConsumer::~MSChromatogramParquetConsumer() = default;

  void MSChromatogramParquetConsumer::consumeSpectrum(SpectrumType& s)
  {
    impl_->consumeSpectrum(s);
  }

  void MSChromatogramParquetConsumer::consumeChromatogram(ChromatogramType& c)
  {
    impl_->consumeChromatogram(c);
  }

  void MSChromatogramParquetConsumer::setExpectedSize(Size expectedSpectra, Size expectedChromatograms)
  {
    (void)expectedSpectra;
    impl_->setExpectedSize(expectedSpectra, expectedChromatograms);
  }

  void MSChromatogramParquetConsumer::setExperimentalSettings(const ExperimentalSettings& exp)
  {
    (void)exp;
    impl_->setExperimentalSettings(exp);
  }

  void MSChromatogramParquetConsumer::finalize()
  {
    impl_->finalize();
  }

} // namespace OpenMS
