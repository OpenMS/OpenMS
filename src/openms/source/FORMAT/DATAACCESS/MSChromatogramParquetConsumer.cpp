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
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/ZlibCompression.h>

#include <memory>
#ifdef WITH_PARQUET
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#endif

#include <unordered_map>
#include <vector>

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

    int64_t parseOrAssignId_(const String& text, int64_t& next_id)
    {
      try
      {
        int64_t value = text.toInt64();
        if (value >= next_id)
        {
          next_id = value + 1;
        }
        return value;
      }
      catch (Exception::ConversionError&)
      {
        return next_id++;
      }
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
#endif
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
#ifndef WITH_PARQUET
      (void)transition_exp;
      throw Exception::MissingFeature(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "OpenMS was built without Parquet support");
#else
      buildTransitionMaps_(transition_exp);
#endif
    }

    ~MSChromatogramParquetConsumerImpl()
    {
#ifdef WITH_PARQUET
      if (wrote_) return;
      write_();
#endif
    }

    void consumeSpectrum(MSChromatogramParquetConsumer::SpectrumType& s)
    {
      s.clear(false);
    }

    void consumeChromatogram(MSChromatogramParquetConsumer::ChromatogramType& c)
    {
#ifndef WITH_PARQUET
      (void)c;
#else
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
          appendOrThrow_(precursor_id_builder_.AppendNull(), "PRECURSOR_ID");
          appendOrThrow_(modified_sequence_builder_.AppendNull(), "MODIFIED_SEQUENCE");
          appendOrThrow_(precursor_charge_builder_.AppendNull(), "PRECURSOR_CHARGE");
          appendOrThrow_(precursor_decoy_builder_.AppendNull(), "PRECURSOR_DECOY");
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
          appendOrThrow_(precursor_id_builder_.AppendNull(), "PRECURSOR_ID");
          appendOrThrow_(transition_id_builder_.AppendNull(), "TRANSITION_ID");
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
      }

      encodeChromatogram_(c);
      c.clear(false);
#endif
    }

    void setExpectedSize(Size, Size expectedChromatograms)
    {
#ifdef WITH_PARQUET
      if (expectedChromatograms > 0)
      {
        run_id_builder_.Reserve(expectedChromatograms);
        source_file_builder_.Reserve(expectedChromatograms);
        ms_level_builder_.Reserve(expectedChromatograms);
        precursor_id_builder_.Reserve(expectedChromatograms);
        transition_id_builder_.Reserve(expectedChromatograms);
        modified_sequence_builder_.Reserve(expectedChromatograms);
        precursor_charge_builder_.Reserve(expectedChromatograms);
        product_charge_builder_.Reserve(expectedChromatograms);
        detecting_transition_builder_.Reserve(expectedChromatograms);
        precursor_decoy_builder_.Reserve(expectedChromatograms);
        product_decoy_builder_.Reserve(expectedChromatograms);
        transition_ordinal_builder_.Reserve(expectedChromatograms);
        transition_type_builder_.Reserve(expectedChromatograms);
        annotation_builder_.Reserve(expectedChromatograms);
        rt_data_builder_.Reserve(expectedChromatograms);
        intensity_data_builder_.Reserve(expectedChromatograms);
        rt_compression_builder_.Reserve(expectedChromatograms);
        intensity_compression_builder_.Reserve(expectedChromatograms);
      }
#endif
    }

    void setExperimentalSettings(const ExperimentalSettings&)
    {
    }

  private:
    String filename_;
    UInt64 run_id_{0};
    String source_file_;
    bool wrote_{false};
#ifdef WITH_PARQUET
    void buildTransitionMaps_(const OpenSwath::LightTargetedExperiment& transition_exp)
    {
      int64_t next_precursor_id = 1;
      int64_t next_transition_id = 1;

      std::unordered_map<String, int64_t> precursor_ids;
      precursor_ids.reserve(transition_exp.getCompounds().size());

      std::unordered_map<String, int64_t> precursor_decoy;
      precursor_decoy.reserve(transition_exp.getCompounds().size());

      for (const auto& compound : transition_exp.getCompounds())
      {
        const String compound_id = compound.id;
        const int64_t precursor_id = parseOrAssignId_(compound_id, next_precursor_id);
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
          transition_id = parseOrAssignId_(transition_name, next_transition_id);
          transition_ids_[transition_name] = transition_id;
        }
        else
        {
          transition_id = it->second;
        }

        auto precursor_it = precursor_ids.find(peptide_ref);
        if (precursor_it == precursor_ids.end())
        {
          const int64_t precursor_id = parseOrAssignId_(peptide_ref, next_precursor_id);
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
        std::string rt_bytes(reinterpret_cast<const char*>(rt_data.data()), rt_data.size() * sizeof(double));
        ZlibCompression::compressString(rt_bytes, rt_encoded);
        std::string int_bytes(reinterpret_cast<const char*>(intensity_data.data()), intensity_data.size() * sizeof(double));
        ZlibCompression::compressString(int_bytes, int_encoded);
      }

      appendBinary_(rt_data_builder_, rt_encoded, "RT_DATA");
      appendBinary_(intensity_data_builder_, int_encoded, "INTENSITY_DATA");
      appendOrThrow_(rt_compression_builder_.Append(rt_compression), "RT_COMPRESSION");
      appendOrThrow_(intensity_compression_builder_.Append(intensity_compression), "INTENSITY_COMPRESSION");
    }

    void write_()
    {
      auto schema = arrow::schema({
        arrow::field("RUN_ID", arrow::int64()),
        arrow::field("SOURCE_FILE", arrow::utf8()),
        arrow::field("MS_LEVEL", arrow::int64()),
        arrow::field("PRECURSOR_ID", arrow::int64()),
        arrow::field("TRANSITION_ID", arrow::int64()),
        arrow::field("MODIFIED_SEQUENCE", arrow::utf8()),
        arrow::field("PRECURSOR_CHARGE", arrow::int64()),
        arrow::field("PRODUCT_CHARGE", arrow::int64()),
        arrow::field("DETECTING_TRANSITION", arrow::int64()),
        arrow::field("PRECURSOR_DECOY", arrow::int64()),
        arrow::field("PRODUCT_DECOY", arrow::int64()),
        arrow::field("TRANSITION_ORDINAL", arrow::int64()),
        arrow::field("TRANSITION_TYPE", arrow::utf8()),
        arrow::field("ANNOTATION", arrow::utf8()),
        arrow::field("RT_DATA", arrow::binary()),
        arrow::field("INTENSITY_DATA", arrow::binary()),
        arrow::field("RT_COMPRESSION", arrow::int64()),
        arrow::field("INTENSITY_COMPRESSION", arrow::int64())
      });

      auto table = arrow::Table::Make(schema, {
        finishArray_(run_id_builder_, "RUN_ID"),
        finishArray_(source_file_builder_, "SOURCE_FILE"),
        finishArray_(ms_level_builder_, "MS_LEVEL"),
        finishArray_(precursor_id_builder_, "PRECURSOR_ID"),
        finishArray_(transition_id_builder_, "TRANSITION_ID"),
        finishArray_(modified_sequence_builder_, "MODIFIED_SEQUENCE"),
        finishArray_(precursor_charge_builder_, "PRECURSOR_CHARGE"),
        finishArray_(product_charge_builder_, "PRODUCT_CHARGE"),
        finishArray_(detecting_transition_builder_, "DETECTING_TRANSITION"),
        finishArray_(precursor_decoy_builder_, "PRECURSOR_DECOY"),
        finishArray_(product_decoy_builder_, "PRODUCT_DECOY"),
        finishArray_(transition_ordinal_builder_, "TRANSITION_ORDINAL"),
        finishArray_(transition_type_builder_, "TRANSITION_TYPE"),
        finishArray_(annotation_builder_, "ANNOTATION"),
        finishArray_(rt_data_builder_, "RT_DATA"),
        finishArray_(intensity_data_builder_, "INTENSITY_DATA"),
        finishArray_(rt_compression_builder_, "RT_COMPRESSION"),
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
                                      "Failed to write chromatogram parquet table", status.ToString());
      }
      wrote_ = true;
    }

    std::unordered_map<String, CompoundInfo> compound_info_;
    std::unordered_map<String, TransitionInfo> transition_info_;
    std::unordered_map<String, int64_t> transition_ids_;

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
#endif
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

} // namespace OpenMS
