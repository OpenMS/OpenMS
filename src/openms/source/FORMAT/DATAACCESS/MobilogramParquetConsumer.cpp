// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MobilogramParquetConsumer.h>
#include "ParquetConsumerHelpers.h"
#include <OpenMS/FORMAT/ParquetFile.h>

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
    using Internal::ParquetConsumerHelpers::appendOrThrow;
    using Internal::ParquetConsumerHelpers::reserveOrThrow;
    using Internal::ParquetConsumerHelpers::appendOptionalInt;
    using Internal::ParquetConsumerHelpers::appendOptionalString;
    using Internal::ParquetConsumerHelpers::appendBinary;
    using Internal::ParquetConsumerHelpers::IdAssigner;
    using Internal::ParquetConsumerHelpers::buildAnnotation;
    using Internal::ParquetConsumerHelpers::CompoundInfo;
    using Internal::ParquetConsumerHelpers::TransitionInfo;

    void appendOptionalDouble_(arrow::DoubleBuilder& builder, double value, const char* column)
    {
      if (!std::isfinite(value))
      {
        appendOrThrow(builder.AppendNull(), column);
      }
      else
      {
        appendOrThrow(builder.Append(value), column);
      }
    }

  } // namespace

  class MobilogramParquetConsumerImpl
  {
  public:
    MobilogramParquetConsumerImpl(const std::string& filename,
                                  UInt64 run_id,
                                  const std::string& source_file,
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
                           const std::string& mobilogram_type,
                           Int64 ms_level,
                           Int64 transition_id,
                           const std::string& transition_native_id,
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
        appendOrThrow(run_id_builder_.Append(static_cast<int64_t>(run_id_)), "RUN_ID");
        appendOptionalString(source_file_builder_, source_file_, "SOURCE_FILE");

        appendOptionalInt(ms_level_builder_, ms_level >= 0, ms_level, "MS_LEVEL");
        appendOptionalString(mobilogram_type_builder_, mobilogram_type, "MOBILOGRAM_TYPE");
        auto tr_it = transition_info_.find(transition_native_id);
        if (tr_it != transition_info_.end())
        {
          const TransitionInfo& tr = tr_it->second;
          const Int64 resolved_transition_id = transition_id >= 0 ? transition_id : tr.transition_id;
          appendOptionalInt(precursor_id_builder_, true, tr.precursor_id, "PRECURSOR_ID");
          appendOptionalInt(transition_id_builder_, true, resolved_transition_id, "TRANSITION_ID");
          appendOptionalInt(feature_id_builder_, feature_id >= 0, feature_id, "FEATURE_ID");
          appendOptionalDouble_(feature_rt_builder_, feature_rt, "FEATURE_RT");
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
        else
        {
          appendOrThrow(precursor_id_builder_.AppendNull(), "PRECURSOR_ID");
          appendOptionalInt(transition_id_builder_, transition_id >= 0, transition_id, "TRANSITION_ID");
          appendOptionalInt(feature_id_builder_, feature_id >= 0, feature_id, "FEATURE_ID");
          appendOptionalDouble_(feature_rt_builder_, feature_rt, "FEATURE_RT");
          appendOrThrow(modified_sequence_builder_.AppendNull(), "MODIFIED_SEQUENCE");
          appendOrThrow(precursor_charge_builder_.AppendNull(), "PRECURSOR_CHARGE");
          appendOrThrow(product_charge_builder_.AppendNull(), "PRODUCT_CHARGE");
          appendOrThrow(detecting_transition_builder_.AppendNull(), "DETECTING_TRANSITION");
          appendOrThrow(precursor_decoy_builder_.AppendNull(), "PRECURSOR_DECOY");
          appendOrThrow(product_decoy_builder_.AppendNull(), "PRODUCT_DECOY");
          appendOrThrow(transition_ordinal_builder_.AppendNull(), "TRANSITION_ORDINAL");
          appendOrThrow(transition_type_builder_.AppendNull(), "TRANSITION_TYPE");
          appendOrThrow(annotation_builder_.AppendNull(), "ANNOTATION");
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
        reserveOrThrow(run_id_builder_.Reserve(expectedMobilograms), "RUN_ID");
        reserveOrThrow(source_file_builder_.Reserve(expectedMobilograms), "SOURCE_FILE");
        reserveOrThrow(ms_level_builder_.Reserve(expectedMobilograms), "MS_LEVEL");
        reserveOrThrow(mobilogram_type_builder_.Reserve(expectedMobilograms), "MOBILOGRAM_TYPE");
        reserveOrThrow(precursor_id_builder_.Reserve(expectedMobilograms), "PRECURSOR_ID");
        reserveOrThrow(transition_id_builder_.Reserve(expectedMobilograms), "TRANSITION_ID");
        reserveOrThrow(feature_id_builder_.Reserve(expectedMobilograms), "FEATURE_ID");
        reserveOrThrow(feature_rt_builder_.Reserve(expectedMobilograms), "FEATURE_RT");
        reserveOrThrow(modified_sequence_builder_.Reserve(expectedMobilograms), "MODIFIED_SEQUENCE");
        reserveOrThrow(precursor_charge_builder_.Reserve(expectedMobilograms), "PRECURSOR_CHARGE");
        reserveOrThrow(product_charge_builder_.Reserve(expectedMobilograms), "PRODUCT_CHARGE");
        reserveOrThrow(detecting_transition_builder_.Reserve(expectedMobilograms), "DETECTING_TRANSITION");
        reserveOrThrow(precursor_decoy_builder_.Reserve(expectedMobilograms), "PRECURSOR_DECOY");
        reserveOrThrow(product_decoy_builder_.Reserve(expectedMobilograms), "PRODUCT_DECOY");
        reserveOrThrow(transition_ordinal_builder_.Reserve(expectedMobilograms), "TRANSITION_ORDINAL");
        reserveOrThrow(transition_type_builder_.Reserve(expectedMobilograms), "TRANSITION_TYPE");
        reserveOrThrow(annotation_builder_.Reserve(expectedMobilograms), "ANNOTATION");
        reserveOrThrow(mobility_data_builder_.Reserve(expectedMobilograms), "MOBILITY_DATA");
        reserveOrThrow(intensity_data_builder_.Reserve(expectedMobilograms), "INTENSITY_DATA");
        reserveOrThrow(mobility_compression_builder_.Reserve(expectedMobilograms), "MOBILITY_COMPRESSION");
        reserveOrThrow(intensity_compression_builder_.Reserve(expectedMobilograms), "INTENSITY_COMPRESSION");
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
        writer_.reset();
        if (!status.ok())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Failed to close mobilogram parquet writer", status.ToString());
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
                                        "Failed to close mobilogram parquet output stream", status.ToString());
        }
      }
    }

  private:
    std::string filename_;
    UInt64 run_id_{0};
    std::string source_file_;
    bool wrote_{false};
    std::mutex write_mutex_;
    struct EncodedMobilogram
    {
      std::string mobility_encoded;
      std::string intensity_encoded;
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
      IdAssigner precursor_id_assigner;
      IdAssigner transition_id_assigner;

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

      // Claim every explicit ID first: a fresh ID handed out below must not take a value that an
      // entry later in the list carries itself.
      for (const auto& compound : transition_exp.getCompounds())
      {
        precursor_id_assigner.reserve(compound.id);
      }
      for (const auto& transition : transition_exp.getTransitions())
      {
        transition_id_assigner.reserve(transition.transition_name);
        precursor_id_assigner.reserve(transition.peptide_ref);
      }

      std::unordered_map<std::string, int64_t> precursor_ids;
      precursor_ids.reserve(transition_exp.getCompounds().size());

      for (const auto& compound : transition_exp.getCompounds())
      {
        const std::string compound_id = compound.id;
        const int64_t precursor_id = precursor_id_assigner.assign(compound_id);
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
          transition_id = transition_id_assigner.assign(transition_name);
          transition_ids_[transition_name] = transition_id;
        }
        else
        {
          transition_id = id_it->second;
        }

        auto precursor_it = precursor_ids.find(peptide_ref);
        if (precursor_it == precursor_ids.end())
        {
          const int64_t precursor_id = precursor_id_assigner.assign(peptide_ref);
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
        info.annotation = buildAnnotation(info.transition_type, info.transition_ordinal, info.product_charge);

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
        ParquetFile::finishArray(run_id_builder_,                fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(source_file_builder_,           fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(ms_level_builder_,              fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(mobilogram_type_builder_,       fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(precursor_id_builder_,          fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(transition_id_builder_,         fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(feature_id_builder_,            fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(feature_rt_builder_,            fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(modified_sequence_builder_,     fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(precursor_charge_builder_,      fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(product_charge_builder_,        fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(detecting_transition_builder_,  fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(precursor_decoy_builder_,       fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(product_decoy_builder_,         fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(transition_ordinal_builder_,    fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(transition_type_builder_,       fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(annotation_builder_,            fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(mobility_data_builder_,         fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(intensity_data_builder_,        fields.at(i++)->name().c_str()),
        ParquetFile::finishArray(mobility_compression_builder_,  fields.at(i++)->name().c_str()),
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

      std::string mob_uncomp;
      MSNumpressCoder().encodeNPRaw(mobility, mob_uncomp, npcfg_m);
      ZlibCompression::compressString(mob_uncomp, encoded.mobility_encoded);

      std::string int_uncomp;
      MSNumpressCoder().encodeNPRaw(intensity, int_uncomp, npcfg_i);
      ZlibCompression::compressString(int_uncomp, encoded.intensity_encoded);

      return encoded;
    }

    void appendEncodedMobilogram_(const EncodedMobilogram& encoded)
    {
      appendBinary(mobility_data_builder_, encoded.mobility_encoded, "MOBILITY_DATA");
      appendBinary(intensity_data_builder_, encoded.intensity_encoded, "INTENSITY_DATA");
      appendOrThrow(mobility_compression_builder_.Append(encoded.mobility_compression), "MOBILITY_COMPRESSION");
      appendOrThrow(intensity_compression_builder_.Append(encoded.intensity_compression), "INTENSITY_COMPRESSION");
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

    std::unordered_map<std::string, CompoundInfo> compound_info_;
    std::unordered_map<std::string, TransitionInfo> transition_info_;
    std::unordered_map<std::string, int64_t> transition_ids_;
  };

  MobilogramParquetConsumer::MobilogramParquetConsumer(const std::string& filename,
                                                       UInt64 run_id,
                                                       const std::string& source_file,
                                                       const OpenSwath::LightTargetedExperiment& transition_exp)
  {
    impl_ = std::make_unique<MobilogramParquetConsumerImpl>(filename, run_id, source_file, transition_exp);
  }

  MobilogramParquetConsumer::~MobilogramParquetConsumer() = default;

  void MobilogramParquetConsumer::consumeMobilogram(const Mobilogram& m,
                                                    const std::string& mobilogram_type,
                                                    Int64 ms_level,
                                                    Int64 transition_id,
                                                    const std::string& transition_native_id,
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
