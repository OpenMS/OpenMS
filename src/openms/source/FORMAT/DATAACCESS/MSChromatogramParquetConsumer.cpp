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

#include <condition_variable>
#include <deque>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <thread>
#include <utility>
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
                                      std::string("Failed to append value for ") + column, status.ToString());
      }
    }

    void reserveOrThrow_(const arrow::Status& status, const char* column)
    {
      if (!status.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      std::string("Failed to reserve capacity for ") + column, status.ToString());
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
                                      std::string("Failed to finish array for ") + name, status.ToString());
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

    void appendOptionalString_(arrow::StringBuilder& builder, const std::string& value, const char* column)
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

    void appendBinary_(arrow::BinaryBuilder& builder, const std::string& value, const char* column)
    {
      if (value.size() > static_cast<size_t>(std::numeric_limits<int32_t>::max()))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      std::string("Binary value for column ") + column + " exceeds Arrow int32 size limit (2 GB).",
                                      StringUtils::toStr(value.size()));
      }
      appendOrThrow_(builder.Append(reinterpret_cast<const uint8_t*>(value.c_str()),
                                    static_cast<int32_t>(value.size())), column);
    }

    struct CompoundInfo
    {
      int64_t precursor_id = 0;
      std::string modified_sequence;
      int64_t precursor_charge = 0;
      int64_t precursor_decoy = 0;
    };

    struct TransitionInfo
    {
      int64_t transition_id = 0;
      int64_t precursor_id = 0;
      std::string modified_sequence;
      int64_t precursor_charge = 0;
      int64_t product_charge = 0;
      int64_t detecting_transition = 0;
      int64_t precursor_decoy = 0;
      int64_t product_decoy = 0;
      int64_t transition_ordinal = 0;
      std::string transition_type;
      std::string annotation;
    };

    /// One encoded chromatogram: what consumeChromatogram() computes before it touches the builders.
    /// The metadata pointers refer to the maps of the consumer that encoded the row.
    struct EncodedRow
    {
      const CompoundInfo* compound = nullptr;     ///< set for precursor (MS1) traces
      const TransitionInfo* transition = nullptr; ///< set for fragment (MS2) traces
      std::string precursor_annotation;           ///< annotation of a precursor trace
      std::string rt_data;                        ///< encoded RT array (empty for an empty chromatogram)
      std::string intensity_data;                 ///< encoded intensity array
    };

    int64_t parseOrAssignId_(const std::string& text, int64_t& next_id, std::unordered_set<int64_t>& used_ids)
    {
      try
      {
        int64_t value = StringUtils::toInt64(text);
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
      while (used_ids.contains(next_id))
      {
        ++next_id;
      }
      used_ids.insert(next_id);
      return next_id++;
    }

    std::string buildAnnotation_(const std::string& transition_type, int64_t ordinal, int64_t charge)
    {
      if (transition_type.empty() || ordinal < 0)
      {
        return "";
      }
      std::string annotation = transition_type + StringUtils::toStr(ordinal);
      if (charge > 0)
      {
        annotation += "^" + StringUtils::toStr(charge);
      }
      return annotation;
    }

    std::string buildPrecursorAnnotation_(const std::string& native_id)
    {
      const std::string tag = "_Precursor_i";
      const Size pos = native_id.rfind(tag);
      if (pos != std::string::npos)
      {
        return StringUtils::substr(native_id, pos + 1);
      }
      return "";
    }
  } // namespace

  class MSChromatogramParquetConsumerImpl
  {
  public:
    MSChromatogramParquetConsumerImpl(const std::string& filename,
                                      UInt64 run_id,
                                      const std::string& source_file,
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
      stopWriterThread_(); // no-op if finalize() got that far
    }

    void consumeSpectrum(MSChromatogramParquetConsumer::SpectrumType& s)
    {
      s.clear(false);
    }

    void consumeChromatogram(MSChromatogramParquetConsumer::ChromatogramType& c)
    {
      EncodedRow row = encodeRow(c);
      appendRow_(row);
      c.clear(false);
    }

    /// Looks up the metadata of @p c and compresses its arrays. Does not modify the consumer, so it
    /// can run concurrently with other calls; appendRow_() does the serialized part.
    EncodedRow encodeRow(const MSChromatogramParquetConsumer::ChromatogramType& c) const
    {
      EncodedRow row;
      const std::string native_id = c.getNativeID();
      const bool is_precursor = StringUtils::hasSubstring(native_id, "_Precursor_i");
      if (is_precursor)
      {
        row.precursor_annotation = buildPrecursorAnnotation_(native_id);
        const std::string group_id = OpenSwathHelper::computeTransitionGroupId(native_id);
        auto comp_it = compound_info_.find(group_id);
        if (comp_it == compound_info_.end())
        {
          // If we have an extracted ion chromatogram to write out to disk, we likely would want to know what precursor the XIC belongs to. Otherwise we could end up with a parquet file with a log of XICs but without any meta data to identify them.
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        std::string("Chromatogram native ID '") + native_id + "' looks like a precursor but no matching compound metadata entry was found. Please ensure the native ID is correct and that the transition library contains metadata for this precursor.", std::string{});
        }
        row.compound = &comp_it->second;
      }
      else
      {
        auto tr_it = transition_info_.find(native_id);
        if (tr_it == transition_info_.end())
        {
          // Same comment as for the precursor case applies here: if we have a chromatogram to write out to disk, we likely would want to know what transition it belongs to. Otherwise we could end up with a parquet file with a log of chromatograms but without any meta data to identify them.
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        std::string("Chromatogram native ID '") + native_id + "' does not have a matching transition metadata entry. Please ensure the native ID is correct and that the transition library contains metadata for this transition.", std::string{});
        }
        row.transition = &tr_it->second;
      }
      encodeArrays_(c, row.rt_data, row.intensity_data);
      return row;
    }

    /// Appends the rows to the builders (flushing as needed) and leaves @p rows empty.
    void appendRows(std::vector<EncodedRow>& rows)
    {
      for (EncodedRow& row : rows)
      {
        appendRow_(row);
      }
      rows.clear();
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
    std::string filename_;
    UInt64 run_id_{0};
    std::string source_file_;
    bool wrote_{false};
    void buildTransitionMaps_(const OpenSwath::LightTargetedExperiment& transition_exp)
    {
      // Build lookup tables for precursor- and transition-level metadata.

      std::unordered_map<std::string, bool> precursor_has_target_transition;
      std::unordered_map<std::string, bool> precursor_has_decoy_transition;
      precursor_has_target_transition.reserve(transition_exp.getTransitions().size());
      precursor_has_decoy_transition.reserve(transition_exp.getTransitions().size());
      for (const auto& transition : transition_exp.getTransitions())
      {
        if (transition.getDecoy())
        {
          precursor_has_decoy_transition[transition.peptide_ref] = true;
        }
        else
        {
          precursor_has_target_transition[transition.peptide_ref] = true;
        }
      }

      // IPF libraries can contain target detecting transitions and decoy
      // identifying transitions for the same precursor. Treat the precursor as
      // decoy only if it has no target transitions at all.
      const auto resolve_precursor_decoy = [&](const std::string& peptide_ref) -> int64_t
      {
        const bool has_target = precursor_has_target_transition.find(peptide_ref) != precursor_has_target_transition.end();
        const bool has_decoy = precursor_has_decoy_transition.find(peptide_ref) != precursor_has_decoy_transition.end();
        return (!has_target && has_decoy) ? 1 : 0;
      };

      std::unordered_map<std::string, int64_t> precursor_ids;
      precursor_ids.reserve(transition_exp.getCompounds().size());
      for (const auto& compound : transition_exp.getCompounds())
      {
        const std::string compound_id = compound.id;
        const int64_t precursor_id = parseOrAssignId_(compound_id, next_precursor_id_, used_precursor_ids_);
        precursor_ids[compound_id] = precursor_id;

        CompoundInfo info;
        info.precursor_id = precursor_id;
        info.modified_sequence = compound.sequence;
        info.precursor_charge = compound.charge;
        info.precursor_decoy = resolve_precursor_decoy(compound_id);
        compound_info_.emplace(compound_id, std::move(info));
      }

      for (const auto& transition : transition_exp.getTransitions())
      {
        const std::string transition_name = transition.transition_name;
        const std::string peptide_ref = transition.peptide_ref;

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

          CompoundInfo info;
          info.precursor_id = precursor_id;
          info.modified_sequence = "";
          info.precursor_charge = 0;
          info.precursor_decoy = resolve_precursor_decoy(peptide_ref);
          compound_info_.emplace(peptide_ref, std::move(info));
          precursor_it = precursor_ids.find(peptide_ref);
        }

        const bool is_decoy = transition.getDecoy();

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

    void appendRow_(const EncodedRow& row)
    {
      appendOrThrow_(run_id_builder_.Append(static_cast<int64_t>(run_id_)), "RUN_ID");
      appendOptionalString_(source_file_builder_, source_file_, "SOURCE_FILE");
      appendOrThrow_(ms_level_builder_.Append(static_cast<int64_t>(row.compound != nullptr ? 1 : 2)), "MS_LEVEL");
      if (row.compound != nullptr)
      {
        const CompoundInfo& comp = *row.compound;
        appendOptionalInt_(precursor_id_builder_, true, comp.precursor_id, "PRECURSOR_ID");
        appendOptionalString_(modified_sequence_builder_, comp.modified_sequence, "MODIFIED_SEQUENCE");
        appendOptionalInt_(precursor_charge_builder_, comp.precursor_charge > 0, comp.precursor_charge, "PRECURSOR_CHARGE");
        appendOptionalInt_(precursor_decoy_builder_, true, comp.precursor_decoy, "PRECURSOR_DECOY");

        appendOrThrow_(transition_id_builder_.AppendNull(), "TRANSITION_ID");
        appendOrThrow_(product_charge_builder_.AppendNull(), "PRODUCT_CHARGE");
        appendOrThrow_(detecting_transition_builder_.AppendNull(), "DETECTING_TRANSITION");
        appendOrThrow_(product_decoy_builder_.AppendNull(), "PRODUCT_DECOY");
        appendOrThrow_(transition_ordinal_builder_.AppendNull(), "TRANSITION_ORDINAL");
        appendOrThrow_(transition_type_builder_.AppendNull(), "TRANSITION_TYPE");
        appendOptionalString_(annotation_builder_, row.precursor_annotation, "ANNOTATION");
      }
      else
      {
        const TransitionInfo& tr = *row.transition;
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

      appendBinary_(rt_data_builder_, row.rt_data, "RT_DATA");
      appendBinary_(intensity_data_builder_, row.intensity_data, "INTENSITY_DATA");
      appendOrThrow_(rt_compression_builder_.Append(rt_compression_), "RT_COMPRESSION");
      appendOrThrow_(intensity_compression_builder_.Append(intensity_compression_), "INTENSITY_COMPRESSION");

      // Track accumulated binary size and row count for both populated and
      // empty chromatograms so trailing empty rows are flushed as well.
      accumulated_rows_++;
      accumulated_binary_bytes_ += row.rt_data.size();
      accumulated_binary_bytes_ += row.intensity_data.size();

      if (accumulated_rows_ >= flush_row_limit_ || accumulated_binary_bytes_ >= flush_byte_limit_)
      {
        flushBuffers_();
      }
    }

    /// Compresses the RT and intensity arrays of @p c (both stay empty for an empty chromatogram).
    static void encodeArrays_(const MSChromatogramParquetConsumer::ChromatogramType& c,
                              std::string& rt_encoded, std::string& int_encoded)
    {
      if (c.empty())
      {
        return;
      }
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

      if (use_lossy_compression_)
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

        std::string rt_uncompressed;
        MSNumpressCoder().encodeNPRaw(rt_data, rt_uncompressed, npconfig_mz);
        ZlibCompression::compressString(rt_uncompressed, rt_encoded);

        std::string int_uncompressed;
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
    }

    void write_()
    {
      // Finalize builders into arrays and assemble the table and write any
      // remaining data.
      flushBuffers_();

      // All batches are handed over: wait until the writer thread has written them.
      stopWriterThread_();
      if (writer_error_)
      {
        std::rethrow_exception(std::exchange(writer_error_, nullptr));
      }

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
        if (table_submitted_)
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
      const size_t batch_bytes = accumulated_binary_bytes_;
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

      submitTable_(std::move(table), batch_bytes);
    }

    /**
      @brief Hands a batch to the writer thread.

      Writing a batch means compressing it (ZSTD), and the caller holds the lock that serializes
      all chromatogram output. A single thread writes the batches in the order in which they are
      submitted, so the file is the same as when every batch is written in place. The queue is
      bounded by the size of the encoded arrays it holds, so a slow disk still holds the
      producers back; the bound is large enough for the burst at the end of a scheduler wave of
      OpenSwathWorkflow, when all threads deliver their chromatograms at once (a bound of 64
      batches made them wait for the compression again).
    */
    void submitTable_(std::shared_ptr<arrow::Table> table, const size_t batch_bytes)
    {
      std::unique_lock<std::mutex> lock(writer_mutex_);
      writer_cv_.wait(lock, [this]() { return writer_error_ != nullptr || pending_tables_.empty() || pending_bytes_ < max_pending_bytes_; });
      if (writer_error_)
      {
        std::rethrow_exception(writer_error_); // kept, so that finalize() reports it as well
      }
      pending_tables_.push_back({std::move(table), batch_bytes});
      pending_bytes_ += batch_bytes;
      table_submitted_ = true;
      if (!writer_thread_.joinable())
      {
        writer_stop_ = false;
        writer_thread_ = std::thread([this]() { writerLoop_(); });
      }
      lock.unlock();
      writer_cv_.notify_all();
    }

    void writerLoop_()
    {
      while (true)
      {
        PendingTable batch;
        {
          std::unique_lock<std::mutex> lock(writer_mutex_);
          writer_cv_.wait(lock, [this]() { return writer_stop_ || !pending_tables_.empty(); });
          if (pending_tables_.empty())
          {
            return; // stop requested and everything is written
          }
          batch = std::move(pending_tables_.front());
          pending_tables_.pop_front();
        }
        try
        {
          writeTable_(*batch.table);
        }
        catch (...)
        {
          std::lock_guard<std::mutex> lock(writer_mutex_);
          writer_error_ = std::current_exception();
          pending_tables_.clear();
          pending_bytes_ = 0;
          writer_cv_.notify_all();
          return;
        }
        batch.table.reset();
        {
          std::lock_guard<std::mutex> lock(writer_mutex_);
          pending_bytes_ -= batch.bytes;
        }
        writer_cv_.notify_all(); // room for further batches
      }
    }

    /// Lets the writer thread finish the queued batches and joins it. Does not throw.
    void stopWriterThread_() noexcept
    {
      {
        std::lock_guard<std::mutex> lock(writer_mutex_);
        writer_stop_ = true;
      }
      writer_cv_.notify_all();
      if (writer_thread_.joinable())
      {
        writer_thread_.join();
      }
    }

    /// Writes one batch as row group(s); runs on the writer thread only.
    void writeTable_(const arrow::Table& table)
    {
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
        auto writer_result = parquet::arrow::FileWriter::Open(*schema_,
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
      auto status = parquet_writer_->WriteTable(table, 1024);
      if (!status.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to write chromatogram parquet table", status.ToString());
      }
    }

    std::unordered_map<std::string, CompoundInfo> compound_info_;
    std::unordered_map<std::string, TransitionInfo> transition_info_;
    std::unordered_map<std::string, int64_t> transition_ids_;
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
    // Array encoding and the scheme ids stored with every row
    static constexpr bool use_lossy_compression_ = true;
    static constexpr int64_t rt_compression_ = use_lossy_compression_ ? 5 : 1;
    static constexpr int64_t intensity_compression_ = use_lossy_compression_ ? 6 : 1;
    // Flush when this many chromatograms are buffered
    const size_t flush_row_limit_ = 1000;
    // Flush when this many binary bytes are buffered (safety margin below 2^31)
    const size_t flush_byte_limit_ = 1900000000ULL; // ~1.9 GB: well below Arrow's int32 limit
    std::shared_ptr<arrow::io::FileOutputStream> outfile_;
    std::shared_ptr<parquet::WriterProperties> props_;
    std::unique_ptr<parquet::arrow::FileWriter> parquet_writer_;
    // Writer thread and its queue of finished batches (see submitTable_()). outfile_, props_ and
    // parquet_writer_ belong to the writer thread while it runs.
    std::thread writer_thread_;
    std::mutex writer_mutex_;
    std::condition_variable writer_cv_;
    struct PendingTable
    {
      std::shared_ptr<arrow::Table> table;
      size_t bytes = 0; ///< size of the encoded RT and intensity arrays of the batch
    };
    std::deque<PendingTable> pending_tables_;
    size_t pending_bytes_ = 0; ///< encoded arrays queued or being written
    bool writer_stop_ = false;
    bool table_submitted_ = false; ///< a batch (possibly empty) was handed over, so the file will exist
    std::exception_ptr writer_error_;
    static constexpr size_t max_pending_bytes_ = 512ULL * 1024 * 1024;
    std::shared_ptr<arrow::Schema> schema_;
  };

  MSChromatogramParquetConsumer::MSChromatogramParquetConsumer(const std::string& filename,
                                                               UInt64 run_id,
                                                               const std::string& source_file,
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

  struct MSChromatogramParquetConsumer::EncodedChromatograms::Rows
  {
    std::vector<EncodedRow> rows;
  };

  MSChromatogramParquetConsumer::EncodedChromatograms::EncodedChromatograms() :
    rows_(std::make_unique<Rows>())
  {
  }

  MSChromatogramParquetConsumer::EncodedChromatograms::~EncodedChromatograms() = default;
  MSChromatogramParquetConsumer::EncodedChromatograms::EncodedChromatograms(EncodedChromatograms&&) noexcept = default;
  MSChromatogramParquetConsumer::EncodedChromatograms& MSChromatogramParquetConsumer::EncodedChromatograms::operator=(EncodedChromatograms&&) noexcept = default;

  Size MSChromatogramParquetConsumer::EncodedChromatograms::size() const
  {
    return rows_ ? rows_->rows.size() : 0;
  }

  void MSChromatogramParquetConsumer::encodeChromatogram(ChromatogramType& c, EncodedChromatograms& rows) const
  {
    if (!rows.rows_)
    {
      rows.rows_ = std::make_unique<EncodedChromatograms::Rows>(); // moved-from batch
    }
    rows.rows_->rows.push_back(impl_->encodeRow(c));
    c.clear(false);
  }

  void MSChromatogramParquetConsumer::consumeEncodedChromatograms(EncodedChromatograms& rows)
  {
    if (rows.rows_)
    {
      impl_->appendRows(rows.rows_->rows);
    }
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
