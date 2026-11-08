// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/XIPMParquetConsumer.h>
#include "ParquetConsumerHelpers.h"
#include <OpenMS/FORMAT/ParquetFile.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/ZlibCompression.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace OpenMS
{
  namespace
  {
    using Internal::ParquetConsumerHelpers::appendOrThrow;
    using Internal::ParquetConsumerHelpers::reserveOrThrow;
    using Internal::ParquetConsumerHelpers::appendOptionalInt;
    using Internal::ParquetConsumerHelpers::appendOptionalString;
    using Internal::ParquetConsumerHelpers::appendBinary;
    using Internal::ParquetConsumerHelpers::parseOrAssignId;
    using Internal::ParquetConsumerHelpers::buildAnnotation;
    using Internal::ParquetConsumerHelpers::buildPrecursorAnnotation;
    using Internal::ParquetConsumerHelpers::CompoundInfo;
    using Internal::ParquetConsumerHelpers::TransitionInfo;

    void appendOptionalDouble_(arrow::DoubleBuilder& builder, bool has_value, double value, const char* column)
    {
      if (!has_value || !std::isfinite(value))
      {
        appendOrThrow(builder.AppendNull(), column);
      }
      else
      {
        appendOrThrow(builder.Append(value), column);
      }
    }

  } // namespace

  class XIPMParquetConsumerImpl
  {
  public:
    XIPMParquetConsumerImpl(const std::string& filename,
                           const OpenSwath::LightTargetedExperiment& transition_exp) :
      filename_(filename)
    {
      buildTransitionMaps_(transition_exp);
    }

    ~XIPMParquetConsumerImpl()
    {
      try
      {
        finalize();
      }
      catch (const Exception::BaseException& e)
      {
        OPENMS_LOG_ERROR << "Failed to write XIPM parquet file '" << filename_
                         << "': " << e.what() << "\n";
      }
      catch (const std::exception& e)
      {
        OPENMS_LOG_ERROR << "Failed to write XIPM parquet file '" << filename_
                         << "': " << e.what() << "\n";
      }
    }

    void consumePeakMap(const PeakMapExtractor::ExtractedPeakMap& peak_map,
                        UInt64 run_id,
                        const std::string& source_file,
                        Int64 ms_level)
    {
      if (peak_map.mz.size() != peak_map.rt.size() ||
          peak_map.mz.size() != peak_map.ion_mobility.size() ||
          peak_map.mz.size() != peak_map.intensity.size())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "All peak-map arrays need to have the same size.");
      }

      EncodedPeakMap encoded = encodePeakMap_(peak_map);
      const int64_t next_binary_bytes = pending_binary_bytes_ + static_cast<int64_t>(
        encoded.mz_encoded.size() + encoded.rt_encoded.size() +
        encoded.mobility_encoded.size() + encoded.intensity_encoded.size());

      if (pending_rows_ > 0 && next_binary_bytes > MAX_BUFFERED_BINARY_BYTES_)
      {
        flushChunk_();
      }

      if (wrote_)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "XIPMParquetConsumer cannot accept writes after finalize().");
      }

      const bool is_precursor = peak_map.native_id.find("_Precursor_i") != std::string::npos;
      appendOrThrow(run_id_builder_.Append(static_cast<int64_t>(run_id)), "RUN_ID");
      appendOptionalString(source_file_builder_, source_file, "SOURCE_FILE");
      appendOrThrow(ms_level_builder_.Append(ms_level), "MS_LEVEL");
      appendOrThrow(peakmap_type_builder_.Append(is_precursor ? "precursor" : "transition"), "PEAKMAP_TYPE");

      if (is_precursor)
      {
        const std::string group_id = OpenSwathHelper::computeTransitionGroupId(peak_map.native_id);
        auto comp_it = compound_info_.find(group_id);
        if (comp_it == compound_info_.end())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        std::string("Peak map native ID '") + peak_map.native_id +
                                        "' looks like a precursor but no matching compound metadata entry was found.", "");
        }

        const CompoundInfo& comp = comp_it->second;
        appendOptionalInt(precursor_id_builder_, true, comp.precursor_id, "PRECURSOR_ID");
        appendOrThrow(transition_id_builder_.AppendNull(), "TRANSITION_ID");
        appendOptionalString(modified_sequence_builder_, comp.modified_sequence, "MODIFIED_SEQUENCE");
        appendOptionalInt(precursor_charge_builder_, comp.precursor_charge > 0, comp.precursor_charge, "PRECURSOR_CHARGE");
        appendOrThrow(product_charge_builder_.AppendNull(), "PRODUCT_CHARGE");
        appendOrThrow(detecting_transition_builder_.AppendNull(), "DETECTING_TRANSITION");
        appendOptionalInt(precursor_decoy_builder_, true, comp.precursor_decoy, "PRECURSOR_DECOY");
        appendOrThrow(product_decoy_builder_.AppendNull(), "PRODUCT_DECOY");
        appendOrThrow(transition_ordinal_builder_.AppendNull(), "TRANSITION_ORDINAL");
        appendOrThrow(transition_type_builder_.AppendNull(), "TRANSITION_TYPE");
        appendOptionalString(annotation_builder_, buildPrecursorAnnotation(peak_map.native_id), "ANNOTATION");
      }
      else
      {
        auto tr_it = transition_info_.find(peak_map.native_id);
        if (tr_it == transition_info_.end())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        std::string("Peak map native ID '") + peak_map.native_id +
                                        "' does not have a matching transition metadata entry.", "");
        }

        const TransitionInfo& tr = tr_it->second;
        appendOptionalInt(precursor_id_builder_, true, tr.precursor_id, "PRECURSOR_ID");
        appendOptionalInt(transition_id_builder_, true, tr.transition_id, "TRANSITION_ID");
        appendOptionalString(modified_sequence_builder_, tr.modified_sequence, "MODIFIED_SEQUENCE");
        appendOptionalInt(precursor_charge_builder_, tr.precursor_charge > 0, tr.precursor_charge, "PRECURSOR_CHARGE");
        appendOptionalInt(product_charge_builder_, tr.product_charge > 0, tr.product_charge, "PRODUCT_CHARGE");
        appendOptionalInt(detecting_transition_builder_, true, tr.detecting_transition, "DETECTING_TRANSITION");
        appendOptionalInt(precursor_decoy_builder_, true, tr.precursor_decoy, "PRECURSOR_DECOY");
        appendOptionalInt(product_decoy_builder_, true, tr.product_decoy, "PRODUCT_DECOY");
        appendOptionalInt(transition_ordinal_builder_, tr.transition_ordinal >= 0, tr.transition_ordinal, "TRANSITION_ORDINAL");
        appendOptionalString(transition_type_builder_, tr.transition_type, "TRANSITION_TYPE");
        appendOptionalString(annotation_builder_, tr.annotation, "ANNOTATION");
      }

      appendOrThrow(target_mz_builder_.Append(peak_map.target_mz), "TARGET_MZ");
      const bool has_rt_window = (peak_map.rt_end - peak_map.rt_start) > 0.0;
      const bool has_target_rt = std::isfinite(peak_map.target_rt) || has_rt_window;
      const double target_rt = std::isfinite(peak_map.target_rt) ? peak_map.target_rt :
        (has_rt_window ? (peak_map.rt_start + peak_map.rt_end) / 2.0 : 0.0);
      appendOptionalDouble_(target_rt_builder_, has_target_rt, target_rt, "TARGET_RT");
      appendOptionalDouble_(target_im_builder_, peak_map.target_ion_mobility >= 0.0,
                            peak_map.target_ion_mobility, "TARGET_ION_MOBILITY");
      appendOptionalDouble_(rt_start_builder_, has_rt_window, peak_map.rt_start, "RT_START");
      appendOptionalDouble_(rt_end_builder_, has_rt_window, peak_map.rt_end, "RT_END");

      appendEncodedPeakMap_(encoded);
      if (pending_binary_bytes_ >= MAX_BUFFERED_BINARY_BYTES_)
      {
        flushChunk_();
      }
    }

    void setExpectedSize(Size expectedPeakMaps)
    {
      if (expectedPeakMaps == 0)
      {
        return;
      }

      reserveOrThrow(run_id_builder_.Reserve(expectedPeakMaps), "RUN_ID");
      reserveOrThrow(source_file_builder_.Reserve(expectedPeakMaps), "SOURCE_FILE");
      reserveOrThrow(ms_level_builder_.Reserve(expectedPeakMaps), "MS_LEVEL");
      reserveOrThrow(peakmap_type_builder_.Reserve(expectedPeakMaps), "PEAKMAP_TYPE");
      reserveOrThrow(precursor_id_builder_.Reserve(expectedPeakMaps), "PRECURSOR_ID");
      reserveOrThrow(transition_id_builder_.Reserve(expectedPeakMaps), "TRANSITION_ID");
      reserveOrThrow(modified_sequence_builder_.Reserve(expectedPeakMaps), "MODIFIED_SEQUENCE");
      reserveOrThrow(precursor_charge_builder_.Reserve(expectedPeakMaps), "PRECURSOR_CHARGE");
      reserveOrThrow(product_charge_builder_.Reserve(expectedPeakMaps), "PRODUCT_CHARGE");
      reserveOrThrow(detecting_transition_builder_.Reserve(expectedPeakMaps), "DETECTING_TRANSITION");
      reserveOrThrow(precursor_decoy_builder_.Reserve(expectedPeakMaps), "PRECURSOR_DECOY");
      reserveOrThrow(product_decoy_builder_.Reserve(expectedPeakMaps), "PRODUCT_DECOY");
      reserveOrThrow(transition_ordinal_builder_.Reserve(expectedPeakMaps), "TRANSITION_ORDINAL");
      reserveOrThrow(transition_type_builder_.Reserve(expectedPeakMaps), "TRANSITION_TYPE");
      reserveOrThrow(annotation_builder_.Reserve(expectedPeakMaps), "ANNOTATION");
      reserveOrThrow(target_mz_builder_.Reserve(expectedPeakMaps), "TARGET_MZ");
      reserveOrThrow(target_rt_builder_.Reserve(expectedPeakMaps), "TARGET_RT");
      reserveOrThrow(target_im_builder_.Reserve(expectedPeakMaps), "TARGET_ION_MOBILITY");
      reserveOrThrow(rt_start_builder_.Reserve(expectedPeakMaps), "RT_START");
      reserveOrThrow(rt_end_builder_.Reserve(expectedPeakMaps), "RT_END");
      reserveOrThrow(mz_data_builder_.Reserve(expectedPeakMaps), "MZ_DATA");
      reserveOrThrow(rt_data_builder_.Reserve(expectedPeakMaps), "RT_DATA");
      reserveOrThrow(mobility_data_builder_.Reserve(expectedPeakMaps), "MOBILITY_DATA");
      reserveOrThrow(intensity_data_builder_.Reserve(expectedPeakMaps), "INTENSITY_DATA");
      reserveOrThrow(mz_compression_builder_.Reserve(expectedPeakMaps), "MZ_COMPRESSION");
      reserveOrThrow(rt_compression_builder_.Reserve(expectedPeakMaps), "RT_COMPRESSION");
      reserveOrThrow(mobility_compression_builder_.Reserve(expectedPeakMaps), "MOBILITY_COMPRESSION");
      reserveOrThrow(intensity_compression_builder_.Reserve(expectedPeakMaps), "INTENSITY_COMPRESSION");
    }

    void finalize()
    {
      if (wrote_)
      {
        return;
      }

      if (pending_rows_ > 0)
      {
        flushChunk_();
      }
      else if (!writer_)
      {
        openWriterIfNeeded_();
        auto empty_table = buildTableFromBuilders_();
        auto status = writer_->WriteTable(*empty_table, PARQUET_ROW_GROUP_SIZE_);
        if (!status.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to write empty XIPM parquet table", status.ToString());
        }
      }

      if (writer_)
      {
        auto status = writer_->Close();
        writer_.reset();
        if (!status.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to close XIPM parquet writer", status.ToString());
        }
      }
      // Mark as written before closing the file so the destructor cannot
      // reopen and truncate a successfully written file if outfile_->Close()
      // throws afterwards.
      wrote_ = true;
      if (outfile_)
      {
        auto status = outfile_->Close();
        outfile_.reset();
        if (!status.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to close XIPM parquet output stream", status.ToString());
        }
      }
    }

  private:
    struct EncodedPeakMap
    {
      std::string mz_encoded;
      std::string rt_encoded;
      std::string mobility_encoded;
      std::string intensity_encoded;
      int64_t mz_compression = 0;
      int64_t rt_compression = 0;
      int64_t mobility_compression = 0;
      int64_t intensity_compression = 0;
    };

    static constexpr int64_t PARQUET_ROW_GROUP_SIZE_ = 1024;
    static constexpr int64_t MAX_BUFFERED_BINARY_BYTES_ = 256LL * 1024LL * 1024LL;
    static constexpr int64_t COMPRESSION_NP_LINEAR_ZLIB_ = 5;
    static constexpr int64_t COMPRESSION_NP_SLOF_ZLIB_ = 6;

    void buildTransitionMaps_(const OpenSwath::LightTargetedExperiment& transition_exp)
    {
      int64_t next_precursor_id = 1;
      int64_t next_transition_id = 1;
      std::unordered_set<int64_t> used_precursor_ids;
      std::unordered_set<int64_t> used_transition_ids;

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
        const int64_t precursor_id = parseOrAssignId(compound_id, next_precursor_id, used_precursor_ids);
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
        auto id_it = transition_ids_.find(transition_name);
        if (id_it == transition_ids_.end())
        {
          transition_id = parseOrAssignId(transition_name, next_transition_id, used_transition_ids);
          transition_ids_[transition_name] = transition_id;
        }
        else
        {
          transition_id = id_it->second;
        }

        auto precursor_it = precursor_ids.find(peptide_ref);
        if (precursor_it == precursor_ids.end())
        {
          const int64_t precursor_id = parseOrAssignId(peptide_ref, next_precursor_id, used_precursor_ids);
          precursor_ids[peptide_ref] = precursor_id;

          CompoundInfo info;
          info.precursor_id = precursor_id;
          info.precursor_decoy = resolve_precursor_decoy(peptide_ref);
          compound_info_.emplace(peptide_ref, std::move(info));
          precursor_it = precursor_ids.find(peptide_ref);
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
        info.product_decoy = transition.getDecoy() ? 1 : 0;
        info.transition_ordinal = transition.fragment_nr;
        info.transition_type = transition.getFragmentType();
        info.annotation = buildAnnotation(info.transition_type, info.transition_ordinal, info.product_charge);

        transition_info_.emplace(transition_name, std::move(info));
      }
    }

    std::shared_ptr<arrow::Schema> buildSchema_() const
    {
      return XIPMSchema::schema();
    }

    void openWriterIfNeeded_()
    {
      if (writer_)
      {
        return;
      }

      schema_ = buildSchema_();
      auto outfile_result = arrow::io::FileOutputStream::Open(std::string(filename_));
      if (!outfile_result.ok())
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename_);
      }
      outfile_ = outfile_result.ValueOrDie();

      parquet::WriterProperties::Builder builder;
      builder.compression(parquet::Compression::ZSTD);
      builder.compression_level(11);
      auto props = builder.build();

      auto writer_result = parquet::arrow::FileWriter::Open(*schema_,
                                                            arrow::default_memory_pool(),
                                                            outfile_,
                                                            props);
      if (!writer_result.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to open XIPM parquet writer", writer_result.status().ToString());
      }
      writer_ = std::move(writer_result.ValueOrDie());
    }

    std::shared_ptr<arrow::Table> buildTableFromBuilders_()
    {
      if (!schema_)
      {
        schema_ = buildSchema_();
      }

      const auto& fields = schema_->fields();
      int i = 0;
      std::vector<std::shared_ptr<arrow::Array>> arrays = {
        ParquetFile::finishArray(run_id_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(source_file_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(ms_level_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(peakmap_type_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(precursor_id_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(transition_id_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(modified_sequence_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(precursor_charge_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(product_charge_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(detecting_transition_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(precursor_decoy_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(product_decoy_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(transition_ordinal_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(transition_type_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(annotation_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(target_mz_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(target_rt_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(target_im_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(rt_start_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(rt_end_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(mz_data_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(rt_data_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(mobility_data_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(intensity_data_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(mz_compression_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(rt_compression_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(mobility_compression_builder_, fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(intensity_compression_builder_, fields.at(i++)->name().c_str()),
      };
      OPENMS_PRECONDITION(i == schema_->num_fields(),
                          "Column count mismatch: buildSchema_ and buildTableFromBuilders_ are out of sync");
      return arrow::Table::Make(schema_, arrays);
    }

    void flushChunk_()
    {
      if (pending_rows_ == 0)
      {
        return;
      }

      openWriterIfNeeded_();
      auto table = buildTableFromBuilders_();
      auto validation = ArrowSchemaValidation::validate(table, XIPMSchema::schema());
      if (!validation.valid)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "XIPM table schema validation failed: " + validation.toString(), "");
      }

      auto status = writer_->WriteTable(*table, PARQUET_ROW_GROUP_SIZE_);
      if (!status.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to write XIPM parquet chunk", status.ToString());
      }

      pending_rows_ = 0;
      pending_binary_bytes_ = 0;
    }

    static std::string encodeLinear_(const std::vector<double>& input)
    {
      if (input.empty())
      {
        return "";
      }

      MSNumpressCoder::NumpressConfig npcfg;
      npcfg.estimate_fixed_point = true;
      npcfg.numpressErrorTolerance = -1.0;
      npcfg.setCompression("linear");

      std::string uncompressed;
      MSNumpressCoder().encodeNPRaw(input, uncompressed, npcfg);
      std::string encoded;
      ZlibCompression::compressString(uncompressed, encoded);
      return encoded;
    }

    static std::string encodeSlof_(const std::vector<double>& input)
    {
      if (input.empty())
      {
        return "";
      }

      MSNumpressCoder::NumpressConfig npcfg;
      npcfg.estimate_fixed_point = true;
      npcfg.numpressErrorTolerance = -1.0;
      npcfg.setCompression("slof");

      std::string uncompressed;
      MSNumpressCoder().encodeNPRaw(input, uncompressed, npcfg);
      std::string encoded;
      ZlibCompression::compressString(uncompressed, encoded);
      return encoded;
    }

    EncodedPeakMap encodePeakMap_(const PeakMapExtractor::ExtractedPeakMap& peak_map) const
    {
      EncodedPeakMap encoded;
      encoded.mz_compression = COMPRESSION_NP_LINEAR_ZLIB_;
      encoded.rt_compression = COMPRESSION_NP_LINEAR_ZLIB_;
      encoded.mobility_compression = COMPRESSION_NP_LINEAR_ZLIB_;
      encoded.intensity_compression = COMPRESSION_NP_SLOF_ZLIB_;

      encoded.mz_encoded = encodeLinear_(peak_map.mz);
      encoded.rt_encoded = encodeLinear_(peak_map.rt);
      encoded.mobility_encoded = encodeLinear_(peak_map.ion_mobility);
      encoded.intensity_encoded = encodeSlof_(peak_map.intensity);
      return encoded;
    }

    void appendEncodedPeakMap_(const EncodedPeakMap& encoded)
    {
      appendBinary(mz_data_builder_, encoded.mz_encoded, "MZ_DATA");
      appendBinary(rt_data_builder_, encoded.rt_encoded, "RT_DATA");
      appendBinary(mobility_data_builder_, encoded.mobility_encoded, "MOBILITY_DATA");
      appendBinary(intensity_data_builder_, encoded.intensity_encoded, "INTENSITY_DATA");
      appendOrThrow(mz_compression_builder_.Append(encoded.mz_compression), "MZ_COMPRESSION");
      appendOrThrow(rt_compression_builder_.Append(encoded.rt_compression), "RT_COMPRESSION");
      appendOrThrow(mobility_compression_builder_.Append(encoded.mobility_compression), "MOBILITY_COMPRESSION");
      appendOrThrow(intensity_compression_builder_.Append(encoded.intensity_compression), "INTENSITY_COMPRESSION");

      ++pending_rows_;
      pending_binary_bytes_ += static_cast<int64_t>(
        encoded.mz_encoded.size() + encoded.rt_encoded.size() +
        encoded.mobility_encoded.size() + encoded.intensity_encoded.size());
    }

    std::string filename_;
    bool wrote_{false};

    arrow::Int64Builder run_id_builder_;
    arrow::StringBuilder source_file_builder_;
    arrow::Int64Builder ms_level_builder_;
    arrow::StringBuilder peakmap_type_builder_;
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
    arrow::DoubleBuilder target_mz_builder_;
    arrow::DoubleBuilder target_rt_builder_;
    arrow::DoubleBuilder target_im_builder_;
    arrow::DoubleBuilder rt_start_builder_;
    arrow::DoubleBuilder rt_end_builder_;
    arrow::BinaryBuilder mz_data_builder_;
    arrow::BinaryBuilder rt_data_builder_;
    arrow::BinaryBuilder mobility_data_builder_;
    arrow::BinaryBuilder intensity_data_builder_;
    arrow::Int64Builder mz_compression_builder_;
    arrow::Int64Builder rt_compression_builder_;
    arrow::Int64Builder mobility_compression_builder_;
    arrow::Int64Builder intensity_compression_builder_;

    std::shared_ptr<arrow::Schema> schema_;
    std::shared_ptr<arrow::io::FileOutputStream> outfile_;
    std::unique_ptr<parquet::arrow::FileWriter> writer_;
    int64_t pending_rows_{0};
    int64_t pending_binary_bytes_{0};

    std::unordered_map<std::string, CompoundInfo> compound_info_;
    std::unordered_map<std::string, TransitionInfo> transition_info_;
    std::unordered_map<std::string, int64_t> transition_ids_;
  };

  XIPMParquetConsumer::XIPMParquetConsumer(const std::string& filename,
                                           const OpenSwath::LightTargetedExperiment& transition_exp)
  {
    impl_ = std::make_unique<XIPMParquetConsumerImpl>(filename, transition_exp);
  }

  XIPMParquetConsumer::~XIPMParquetConsumer() = default;

  void XIPMParquetConsumer::consumePeakMap(const PeakMapExtractor::ExtractedPeakMap& peak_map,
                                           UInt64 run_id,
                                           const std::string& source_file,
                                           Int64 ms_level)
  {
    impl_->consumePeakMap(peak_map, run_id, source_file, ms_level);
  }

  void XIPMParquetConsumer::finalize()
  {
    impl_->finalize();
  }

  void XIPMParquetConsumer::setExpectedSize(Size expectedPeakMaps)
  {
    impl_->setExpectedSize(expectedPeakMaps);
  }
} // namespace OpenMS
