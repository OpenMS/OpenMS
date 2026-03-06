// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MobilogramParquetConsumer.h>

#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/FORMAT/ZlibCompression.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/METADATA/ExperimentalSettings.h>

#include <memory>
#include <cmath>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef WITH_PARQUET
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#endif

namespace OpenMS
{
  namespace
  {
#ifdef WITH_PARQUET
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
#endif
  } // namespace

  class MobilogramParquetConsumerImpl
  {
  public:
    MobilogramParquetConsumerImpl(const String& filename,
                                  UInt64 run_id,
                                  const String& source_file,
                                  const OpenSwath::LightTargetedExperiment& transition_exp) :
      filename_(filename), run_id_(Internal::SqliteHelper::clearSignBit(run_id)), source_file_(source_file)
    {
#ifndef WITH_PARQUET
      (void)transition_exp;
      (void)source_file_;
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#else
      buildTransitionMaps_(transition_exp);
#endif
    }

    ~MobilogramParquetConsumerImpl()
    {
#ifdef WITH_PARQUET
      try
      {
        finalize();
      }
      catch (const Exception::BaseException& e)
      {
        OPENMS_LOG_ERROR << "Failed to write mobilogram parquet file '" << filename_
                         << "': " << e.what() << "\n";
      }
#endif
    }

    void consumeMobilogram(Mobilogram& m,
                           const String& mobilogram_type,
                           Int64 ms_level,
                           Int64 transition_id,
                           const String& transition_native_id,
                           double feature_rt,
                           Int64 feature_id)
    {
#ifndef WITH_PARQUET
      (void)m;
      (void)mobilogram_type;
      (void)ms_level;
      (void)transition_id;
      (void)transition_native_id;
      (void)feature_rt;
      (void)feature_id;
#else
      std::lock_guard<std::mutex> lock(write_mutex_);

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

      encodeMobilogram_(m);
      m.clear();
#endif
    }

    void setExpectedSize(Size expectedMobilograms)
    {
#ifdef WITH_PARQUET
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
#endif
    }

    void setExperimentalSettings(const ExperimentalSettings&)
    {
    }

    void finalize()
    {
#ifdef WITH_PARQUET
      std::lock_guard<std::mutex> lock(write_mutex_);
      if (wrote_) return;
      write_();
#endif
    }

  private:
    String filename_;
    UInt64 run_id_{0};
    String source_file_;
    bool wrote_{false};
    std::mutex write_mutex_;

#ifdef WITH_PARQUET
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

    void encodeMobilogram_(const Mobilogram& m)
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

      const bool use_lossy = true;
      const int64_t mobility_compression = use_lossy ? 5 : 1;
      const int64_t intensity_compression = use_lossy ? 6 : 1;

      String mobility_encoded;
      String intensity_encoded;
      if (mobility.empty())
      {
        appendBinary_(mobility_data_builder_, mobility_encoded, "MOBILITY_DATA");
        appendBinary_(intensity_data_builder_, intensity_encoded, "INTENSITY_DATA");
        appendOrThrow_(mobility_compression_builder_.Append(mobility_compression), "MOBILITY_COMPRESSION");
        appendOrThrow_(intensity_compression_builder_.Append(intensity_compression), "INTENSITY_COMPRESSION");
        return;
      }

      if (use_lossy)
      {
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
        ZlibCompression::compressString(mob_uncomp, mobility_encoded);

        String int_uncomp;
        MSNumpressCoder().encodeNPRaw(intensity, int_uncomp, npcfg_i);
        ZlibCompression::compressString(int_uncomp, intensity_encoded);
      }
      else
      {
        std::string mob_bytes(reinterpret_cast<const char*>(mobility.data()), mobility.size() * sizeof(double));
        ZlibCompression::compressString(mob_bytes, mobility_encoded);
        std::string int_bytes(reinterpret_cast<const char*>(intensity.data()), intensity.size() * sizeof(double));
        ZlibCompression::compressString(int_bytes, intensity_encoded);
      }

      appendBinary_(mobility_data_builder_, mobility_encoded, "MOBILITY_DATA");
      appendBinary_(intensity_data_builder_, intensity_encoded, "INTENSITY_DATA");
      appendOrThrow_(mobility_compression_builder_.Append(mobility_compression), "MOBILITY_COMPRESSION");
      appendOrThrow_(intensity_compression_builder_.Append(intensity_compression), "INTENSITY_COMPRESSION");
    }

    void write_()
    {
      auto schema = arrow::schema({
        arrow::field("RUN_ID", arrow::int64()),
        arrow::field("SOURCE_FILE", arrow::utf8()),
        arrow::field("MS_LEVEL", arrow::int64()),
        arrow::field("MOBILOGRAM_TYPE", arrow::utf8()),
        arrow::field("PRECURSOR_ID", arrow::int64()),
        arrow::field("TRANSITION_ID", arrow::int64()),
        arrow::field("FEATURE_ID", arrow::int64()),
        arrow::field("FEATURE_RT", arrow::float64()),
        arrow::field("MODIFIED_SEQUENCE", arrow::utf8()),
        arrow::field("PRECURSOR_CHARGE", arrow::int64()),
        arrow::field("PRODUCT_CHARGE", arrow::int64()),
        arrow::field("DETECTING_TRANSITION", arrow::int64()),
        arrow::field("PRECURSOR_DECOY", arrow::int64()),
        arrow::field("PRODUCT_DECOY", arrow::int64()),
        arrow::field("TRANSITION_ORDINAL", arrow::int64()),
        arrow::field("TRANSITION_TYPE", arrow::utf8()),
        arrow::field("ANNOTATION", arrow::utf8()),
        arrow::field("MOBILITY_DATA", arrow::binary()),
        arrow::field("INTENSITY_DATA", arrow::binary()),
        arrow::field("MOBILITY_COMPRESSION", arrow::int64()),
        arrow::field("INTENSITY_COMPRESSION", arrow::int64())
      });

      auto table = arrow::Table::Make(schema, {
        finishArray_(run_id_builder_, "RUN_ID"),
        finishArray_(source_file_builder_, "SOURCE_FILE"),
        finishArray_(ms_level_builder_, "MS_LEVEL"),
        finishArray_(mobilogram_type_builder_, "MOBILOGRAM_TYPE"),
        finishArray_(precursor_id_builder_, "PRECURSOR_ID"),
        finishArray_(transition_id_builder_, "TRANSITION_ID"),
        finishArray_(feature_id_builder_, "FEATURE_ID"),
        finishArray_(feature_rt_builder_, "FEATURE_RT"),
        finishArray_(modified_sequence_builder_, "MODIFIED_SEQUENCE"),
        finishArray_(precursor_charge_builder_, "PRECURSOR_CHARGE"),
        finishArray_(product_charge_builder_, "PRODUCT_CHARGE"),
        finishArray_(detecting_transition_builder_, "DETECTING_TRANSITION"),
        finishArray_(precursor_decoy_builder_, "PRECURSOR_DECOY"),
        finishArray_(product_decoy_builder_, "PRODUCT_DECOY"),
        finishArray_(transition_ordinal_builder_, "TRANSITION_ORDINAL"),
        finishArray_(transition_type_builder_, "TRANSITION_TYPE"),
        finishArray_(annotation_builder_, "ANNOTATION"),
        finishArray_(mobility_data_builder_, "MOBILITY_DATA"),
        finishArray_(intensity_data_builder_, "INTENSITY_DATA"),
        finishArray_(mobility_compression_builder_, "MOBILITY_COMPRESSION"),
        finishArray_(intensity_compression_builder_, "INTENSITY_COMPRESSION")
      });

      auto outfile_result = arrow::io::FileOutputStream::Open(std::string(filename_));
      if (!outfile_result.ok())
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename_);
      }
      auto outfile = outfile_result.ValueOrDie();
      parquet::WriterProperties::Builder builder;
      builder.compression(parquet::Compression::ZSTD);
      builder.compression_level(11);
      auto props = builder.build();
      auto status = parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 1024, props);
      if (!status.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to write mobilogram parquet table", status.ToString());
      }
      wrote_ = true;
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

    std::unordered_map<String, CompoundInfo> compound_info_;
    std::unordered_map<String, TransitionInfo> transition_info_;
    std::unordered_map<String, int64_t> transition_ids_;
#endif
  };

  MobilogramParquetConsumer::MobilogramParquetConsumer(const String& filename,
                                                       UInt64 run_id,
                                                       const String& source_file,
                                                       const OpenSwath::LightTargetedExperiment& transition_exp)
  {
    impl_ = std::make_unique<MobilogramParquetConsumerImpl>(filename, run_id, source_file, transition_exp);
  }

  MobilogramParquetConsumer::~MobilogramParquetConsumer() = default;

  void MobilogramParquetConsumer::consumeMobilogram(Mobilogram& m,
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

  void MobilogramParquetConsumer::setExperimentalSettings(const ExperimentalSettings& exp)
  {
    impl_->setExperimentalSettings(exp);
  }

} // namespace OpenMS
