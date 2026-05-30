// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MobilogramParquetConsumer.h>

#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/ZlibCompression.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <memory>
#include <cmath>
#include <limits>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

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
      if (value.size() > static_cast<Size>(std::numeric_limits<int32_t>::max()))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      String("Failed to append value for ") + column,
                                      "Single binary payload exceeds Arrow binary cell limit");
      }
      appendOrThrow_(builder.Append(reinterpret_cast<const uint8_t*>(value.c_str()),
                                    static_cast<int32_t>(value.size())), column);
    }

    void appendOptionalDouble_(arrow::DoubleBuilder& builder, double value, const char* column)
    {
      if (!std::isfinite(value))
      {
        appendOrThrow_(builder.AppendNull(), column);
      }
      else
      {
        appendOrThrow_(builder.Append(value), column);
      }
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
        // Reject negative parsed IDs: negative values are used as "unset" sentinels elsewhere
        // and must not be preserved as real IDs in the output. Only accept non-negative parsed IDs.
        if (value >= 0 && used_ids.insert(value).second)
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
  } // namespace

  class MobilogramParquetConsumerImpl
  {
  public:
    MobilogramParquetConsumerImpl(const String& filename,
                                  UInt64 run_id,
                                  const String& source_file,
                                  const OpenSwath::LightTargetedExperiment& transition_exp) :
      filename_(filename), run_id_(run_id & ~(1ULL << 63)), source_file_(source_file)
    {
      buildTransitionMaps_(transition_exp);
    }

    ~MobilogramParquetConsumerImpl()
    {
      try
      {
        finalize();
      }
      catch (const Exception::BaseException& e)
      {
        OPENMS_LOG_ERROR << "Failed to write mobilogram parquet file '" << filename_
                         << "': " << e.what() << "\n";
      }
    }

    void consumeMobilogram(const Mobilogram& m,
                           const String& mobilogram_type,
                           Int64 ms_level,
                           Int64 transition_id,
                           const String& transition_native_id,
                           double feature_rt,
                           Int64 feature_id)
    {
      // Do expensive encoding outside the writer mutex to avoid serializing
      // the hot Numpress+zlib work when multiple threads share this consumer.
      EncodedMobilogram encoded = encodeMobilogram_(m);

      {
        std::lock_guard<std::mutex> lock(write_mutex_);

        // Fail fast if the consumer has been finalized. Writes after finalize() are not allowed
        // and would otherwise be silently dropped or reopen/truncate the file.
        if (wrote_)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                           "MobilogramParquetConsumer cannot accept writes after finalize()");
        }

        const int64_t next_binary_bytes = pending_binary_bytes_ +
          static_cast<int64_t>(encoded.mobility_encoded.size() + encoded.intensity_encoded.size());
        if (pending_rows_ > 0 && next_binary_bytes > MAX_BUFFERED_BINARY_BYTES_)
        {
          flushChunk_();
        }

        // Store contextual information for each mobilogram row.
        appendOrThrow_(run_id_builder_.Append(static_cast<int64_t>(run_id_)), "RUN_ID");
        appendOptionalString_(source_file_builder_, source_file_, "SOURCE_FILE");

        appendOptionalInt_(ms_level_builder_, ms_level >= 0, ms_level, "MS_LEVEL");
        appendOptionalString_(mobilogram_type_builder_, mobilogram_type, "MOBILOGRAM_TYPE");
        auto tr_it = transition_info_.find(transition_native_id);
        if (tr_it != transition_info_.end())
        {
          const TransitionInfo& tr = tr_it->second;
          const Int64 resolved_transition_id = transition_id >= 0 ? transition_id : tr.transition_id;
          appendOptionalInt_(precursor_id_builder_, true, tr.precursor_id, "PRECURSOR_ID");
          appendOptionalInt_(transition_id_builder_, true, resolved_transition_id, "TRANSITION_ID");
          appendOptionalInt_(feature_id_builder_, feature_id >= 0, feature_id, "FEATURE_ID");
          appendOptionalDouble_(feature_rt_builder_, feature_rt, "FEATURE_RT");
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
          appendOrThrow_(precursor_id_builder_.AppendNull(), "PRECURSOR_ID");
          appendOptionalInt_(transition_id_builder_, transition_id >= 0, transition_id, "TRANSITION_ID");
          appendOptionalInt_(feature_id_builder_, feature_id >= 0, feature_id, "FEATURE_ID");
          appendOptionalDouble_(feature_rt_builder_, feature_rt, "FEATURE_RT");
          appendOrThrow_(modified_sequence_builder_.AppendNull(), "MODIFIED_SEQUENCE");
          appendOrThrow_(precursor_charge_builder_.AppendNull(), "PRECURSOR_CHARGE");
          appendOrThrow_(product_charge_builder_.AppendNull(), "PRODUCT_CHARGE");
          appendOrThrow_(detecting_transition_builder_.AppendNull(), "DETECTING_TRANSITION");
          appendOrThrow_(precursor_decoy_builder_.AppendNull(), "PRECURSOR_DECOY");
          appendOrThrow_(product_decoy_builder_.AppendNull(), "PRODUCT_DECOY");
          appendOrThrow_(transition_ordinal_builder_.AppendNull(), "TRANSITION_ORDINAL");
          appendOrThrow_(transition_type_builder_.AppendNull(), "TRANSITION_TYPE");
          appendOrThrow_(annotation_builder_.AppendNull(), "ANNOTATION");
        }

        appendEncodedMobilogram_(encoded);
        if (pending_binary_bytes_ >= MAX_BUFFERED_BINARY_BYTES_)
        {
          flushChunk_();
        }
      }
    }

    void setExpectedSize(Size expectedMobilograms)
    {
      if (expectedMobilograms > 0)
      {
        reserveOrThrow_(run_id_builder_.Reserve(expectedMobilograms), "RUN_ID");
        reserveOrThrow_(source_file_builder_.Reserve(expectedMobilograms), "SOURCE_FILE");
        reserveOrThrow_(ms_level_builder_.Reserve(expectedMobilograms), "MS_LEVEL");
        reserveOrThrow_(mobilogram_type_builder_.Reserve(expectedMobilograms), "MOBILOGRAM_TYPE");
        reserveOrThrow_(precursor_id_builder_.Reserve(expectedMobilograms), "PRECURSOR_ID");
        reserveOrThrow_(transition_id_builder_.Reserve(expectedMobilograms), "TRANSITION_ID");
        reserveOrThrow_(feature_id_builder_.Reserve(expectedMobilograms), "FEATURE_ID");
        reserveOrThrow_(feature_rt_builder_.Reserve(expectedMobilograms), "FEATURE_RT");
        reserveOrThrow_(modified_sequence_builder_.Reserve(expectedMobilograms), "MODIFIED_SEQUENCE");
        reserveOrThrow_(precursor_charge_builder_.Reserve(expectedMobilograms), "PRECURSOR_CHARGE");
        reserveOrThrow_(product_charge_builder_.Reserve(expectedMobilograms), "PRODUCT_CHARGE");
        reserveOrThrow_(detecting_transition_builder_.Reserve(expectedMobilograms), "DETECTING_TRANSITION");
        reserveOrThrow_(precursor_decoy_builder_.Reserve(expectedMobilograms), "PRECURSOR_DECOY");
        reserveOrThrow_(product_decoy_builder_.Reserve(expectedMobilograms), "PRODUCT_DECOY");
        reserveOrThrow_(transition_ordinal_builder_.Reserve(expectedMobilograms), "TRANSITION_ORDINAL");
        reserveOrThrow_(transition_type_builder_.Reserve(expectedMobilograms), "TRANSITION_TYPE");
        reserveOrThrow_(annotation_builder_.Reserve(expectedMobilograms), "ANNOTATION");
        reserveOrThrow_(mobility_data_builder_.Reserve(expectedMobilograms), "MOBILITY_DATA");
        reserveOrThrow_(intensity_data_builder_.Reserve(expectedMobilograms), "INTENSITY_DATA");
        reserveOrThrow_(mobility_compression_builder_.Reserve(expectedMobilograms), "MOBILITY_COMPRESSION");
        reserveOrThrow_(intensity_compression_builder_.Reserve(expectedMobilograms), "INTENSITY_COMPRESSION");
      }
    }

    void finalize()
    {
      std::lock_guard<std::mutex> lock(write_mutex_);
      if (wrote_) return;

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
                                        "Failed to write empty mobilogram parquet table", status.ToString());
        }
      }

      if (writer_)
      {
        auto status = writer_->Close();
        if (!status.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to close mobilogram parquet writer", status.ToString());
        }
        writer_.reset();
      }
      if (outfile_)
      {
        auto status = outfile_->Close();
        if (!status.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to close mobilogram parquet output stream", status.ToString());
        }
        outfile_.reset();
      }
      wrote_ = true;
    }

  private:
    String filename_;
    UInt64 run_id_{0};
    String source_file_;
    bool wrote_{false};
    std::mutex write_mutex_;
    struct EncodedMobilogram
    {
      String mobility_encoded;
      String intensity_encoded;
      int64_t mobility_compression = 0;
      int64_t intensity_compression = 0;
    };

    static constexpr int64_t PARQUET_ROW_GROUP_SIZE_ = 1024;
    static constexpr int64_t MAX_BUFFERED_BINARY_BYTES_ = 256LL * 1024LL * 1024LL;
    // Compression codes matching the decoder in XIMParquetFile.cpp decodeBinary_():
    //   5 = np-linear + zlib (used for mobility data)
    //   6 = np-slof + zlib   (used for intensity data)
    static constexpr int64_t COMPRESSION_NP_LINEAR_ZLIB_ = 5;
    static constexpr int64_t COMPRESSION_NP_SLOF_ZLIB_   = 6;

    void buildTransitionMaps_(const OpenSwath::LightTargetedExperiment& transition_exp)
    {
      int64_t next_precursor_id = 1;
      int64_t next_transition_id = 1;
      std::unordered_set<int64_t> used_precursor_ids;
      std::unordered_set<int64_t> used_transition_ids;

      std::unordered_map<String, int64_t> precursor_ids;
      precursor_ids.reserve(transition_exp.getCompounds().size());

      for (const auto& compound : transition_exp.getCompounds())
      {
        const String compound_id = compound.id;
        const int64_t precursor_id = parseOrAssignId_(compound_id, next_precursor_id, used_precursor_ids);
        precursor_ids[compound_id] = precursor_id;

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
        auto id_it = transition_ids_.find(transition_name);
        if (id_it == transition_ids_.end())
        {
          transition_id = parseOrAssignId_(transition_name, next_transition_id, used_transition_ids);
          transition_ids_[transition_name] = transition_id;
        }
        else
        {
          transition_id = id_it->second;
        }

        auto precursor_it = precursor_ids.find(peptide_ref);
        if (precursor_it == precursor_ids.end())
        {
          const int64_t precursor_id = parseOrAssignId_(peptide_ref, next_precursor_id, used_precursor_ids);
          precursor_ids[peptide_ref] = precursor_id;

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

    std::shared_ptr<arrow::Schema> buildSchema_() const
    {
      return XIMSchema::schema();
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
                                      "Failed to open mobilogram parquet writer", writer_result.status().ToString());
      }
      writer_ = std::move(writer_result.ValueOrDie());
    }

    std::shared_ptr<arrow::Table> buildTableFromBuilders_()
    {
      if (!schema_) schema_ = buildSchema_();

      // Array order MUST match field order in buildSchema_().
      // Use schema field names as the error-message label for each finish call so that
      // any count mismatch is caught immediately by the precondition below.
      const auto& fields = schema_->fields();
      int i = 0;
      std::vector<std::shared_ptr<arrow::Array>> arrays = {
        finishArray_(run_id_builder_,                fields.at(i++)->name().c_str()),
        finishArray_(source_file_builder_,           fields.at(i++)->name().c_str()),
        finishArray_(ms_level_builder_,              fields.at(i++)->name().c_str()),
        finishArray_(mobilogram_type_builder_,       fields.at(i++)->name().c_str()),
        finishArray_(precursor_id_builder_,          fields.at(i++)->name().c_str()),
        finishArray_(transition_id_builder_,         fields.at(i++)->name().c_str()),
        finishArray_(feature_id_builder_,            fields.at(i++)->name().c_str()),
        finishArray_(feature_rt_builder_,            fields.at(i++)->name().c_str()),
        finishArray_(modified_sequence_builder_,     fields.at(i++)->name().c_str()),
        finishArray_(precursor_charge_builder_,      fields.at(i++)->name().c_str()),
        finishArray_(product_charge_builder_,        fields.at(i++)->name().c_str()),
        finishArray_(detecting_transition_builder_,  fields.at(i++)->name().c_str()),
        finishArray_(precursor_decoy_builder_,       fields.at(i++)->name().c_str()),
        finishArray_(product_decoy_builder_,         fields.at(i++)->name().c_str()),
        finishArray_(transition_ordinal_builder_,    fields.at(i++)->name().c_str()),
        finishArray_(transition_type_builder_,       fields.at(i++)->name().c_str()),
        finishArray_(annotation_builder_,            fields.at(i++)->name().c_str()),
        finishArray_(mobility_data_builder_,         fields.at(i++)->name().c_str()),
        finishArray_(intensity_data_builder_,        fields.at(i++)->name().c_str()),
        finishArray_(mobility_compression_builder_,  fields.at(i++)->name().c_str()),
        finishArray_(intensity_compression_builder_, fields.at(i++)->name().c_str()),
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

      // Validate table against registry schema
      auto validation = ArrowSchemaValidation::validate(table, XIMSchema::schema());
      if (!validation.valid)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "XIM table schema validation failed: " + validation.toString(), "");
      }

      auto status = writer_->WriteTable(*table, PARQUET_ROW_GROUP_SIZE_);
      if (!status.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to write mobilogram parquet chunk", status.ToString());
      }

      pending_rows_ = 0;
      pending_binary_bytes_ = 0;
    }

    EncodedMobilogram encodeMobilogram_(const Mobilogram& m)
    {
      std::vector<double> mobility;
      std::vector<double> intensity;
      mobility.reserve(m.size());
      intensity.reserve(m.size());
      for (const auto& p : m)
      {
        mobility.push_back(p.getMobility());
        intensity.push_back(p.getIntensity());
      }

      EncodedMobilogram encoded;
      encoded.mobility_compression = COMPRESSION_NP_LINEAR_ZLIB_;
      encoded.intensity_compression = COMPRESSION_NP_SLOF_ZLIB_;

      if (mobility.empty())
      {
        return encoded;
      }

      MSNumpressCoder::NumpressConfig npcfg_m;
      npcfg_m.estimate_fixed_point = true;
      npcfg_m.numpressErrorTolerance = -1.0;
      npcfg_m.setCompression("linear");

      MSNumpressCoder::NumpressConfig npcfg_i;
      npcfg_i.estimate_fixed_point = true;
      npcfg_i.numpressErrorTolerance = -1.0;
      npcfg_i.setCompression("slof");

      String mob_uncomp;
      MSNumpressCoder().encodeNPRaw(mobility, mob_uncomp, npcfg_m);
      ZlibCompression::compressString(mob_uncomp, encoded.mobility_encoded);

      String int_uncomp;
      MSNumpressCoder().encodeNPRaw(intensity, int_uncomp, npcfg_i);
      ZlibCompression::compressString(int_uncomp, encoded.intensity_encoded);

      return encoded;
    }

    void appendEncodedMobilogram_(const EncodedMobilogram& encoded)
    {
      appendBinary_(mobility_data_builder_, encoded.mobility_encoded, "MOBILITY_DATA");
      appendBinary_(intensity_data_builder_, encoded.intensity_encoded, "INTENSITY_DATA");
      appendOrThrow_(mobility_compression_builder_.Append(encoded.mobility_compression), "MOBILITY_COMPRESSION");
      appendOrThrow_(intensity_compression_builder_.Append(encoded.intensity_compression), "INTENSITY_COMPRESSION");
      ++pending_rows_;
      pending_binary_bytes_ += static_cast<int64_t>(encoded.mobility_encoded.size() + encoded.intensity_encoded.size());
    }

    // Builders
    arrow::Int64Builder run_id_builder_;
    arrow::StringBuilder source_file_builder_;
    arrow::Int64Builder ms_level_builder_;
    arrow::StringBuilder mobilogram_type_builder_;
    arrow::Int64Builder precursor_id_builder_;
    arrow::Int64Builder transition_id_builder_;
    arrow::Int64Builder feature_id_builder_;
    arrow::DoubleBuilder feature_rt_builder_;
    arrow::StringBuilder modified_sequence_builder_;
    arrow::Int64Builder precursor_charge_builder_;
    arrow::Int64Builder product_charge_builder_;
    arrow::Int64Builder detecting_transition_builder_;
    arrow::Int64Builder precursor_decoy_builder_;
    arrow::Int64Builder product_decoy_builder_;
    arrow::Int64Builder transition_ordinal_builder_;
    arrow::StringBuilder transition_type_builder_;
    arrow::StringBuilder annotation_builder_;
    arrow::BinaryBuilder mobility_data_builder_;
    arrow::BinaryBuilder intensity_data_builder_;
    arrow::Int64Builder mobility_compression_builder_;
    arrow::Int64Builder intensity_compression_builder_;
    std::shared_ptr<arrow::Schema> schema_;
    std::shared_ptr<arrow::io::FileOutputStream> outfile_;
    std::unique_ptr<parquet::arrow::FileWriter> writer_;
    int64_t pending_rows_{0};
    int64_t pending_binary_bytes_{0};

    std::unordered_map<String, CompoundInfo> compound_info_;
    std::unordered_map<String, TransitionInfo> transition_info_;
    std::unordered_map<String, int64_t> transition_ids_;
  };

  MobilogramParquetConsumer::MobilogramParquetConsumer(const String& filename,
                                                       UInt64 run_id,
                                                       const String& source_file,
                                                       const OpenSwath::LightTargetedExperiment& transition_exp)
  {
    impl_ = std::make_unique<MobilogramParquetConsumerImpl>(filename, run_id, source_file, transition_exp);
  }

  MobilogramParquetConsumer::~MobilogramParquetConsumer() = default;

  void MobilogramParquetConsumer::consumeMobilogram(const Mobilogram& m,
                                                    const String& mobilogram_type,
                                                    Int64 ms_level,
                                                    Int64 transition_id,
                                                    const String& transition_native_id,
                                                    double feature_rt,
                                                    Int64 feature_id)
  {
    impl_->consumeMobilogram(m, mobilogram_type, ms_level, transition_id, transition_native_id, feature_rt, feature_id);
  }

  void MobilogramParquetConsumer::finalize()
  {
    impl_->finalize();
  }

  void MobilogramParquetConsumer::setExpectedSize(Size expectedMobilograms)
  {
    impl_->setExpectedSize(expectedMobilograms);
  }

} // namespace OpenMS
