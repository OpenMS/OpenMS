// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MobilogramParquetConsumer.h>

#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/ZlibCompression.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <memory>
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
#endif
  } // namespace

  class MobilogramParquetConsumerImpl
  {
  public:
    MobilogramParquetConsumerImpl(const String& filename,
                                  UInt64 run_id,
                                  const String& source_file,
                                  const OpenSwath::LightTargetedExperiment& /*transition_exp*/) :
      filename_(filename), run_id_(run_id), source_file_(source_file)
    {
#ifndef WITH_PARQUET
      (void)source_file_;
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
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

    void consumeMobilogram(Mobilogram& m)
    {
#ifndef WITH_PARQUET
      (void)m;
#else
      // Store RT (retention time) as context and precursor info if present
      appendOrThrow_(run_id_builder_.Append(static_cast<int64_t>(run_id_)), "RUN_ID");
      appendOptionalString_(source_file_builder_, source_file_, "SOURCE_FILE");

      // RT column (mobilogram RT)
      appendOrThrow_(rt_builder_.Append(m.getRT()), "RT");

      // Precursor info may be present in data arrays or metas; leave null for now
      appendOrThrow_(precursor_id_builder_.AppendNull(), "PRECURSOR_ID");
      appendOrThrow_(ident_id_builder_.AppendNull(), "IDENTIFICATION_ID");
      appendOrThrow_(modified_sequence_builder_.AppendNull(), "MODIFIED_SEQUENCE");
      appendOrThrow_(precursor_charge_builder_.AppendNull(), "PRECURSOR_CHARGE");

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
        reserveOrThrow_(rt_builder_.Reserve(expectedMobilograms), "RT");
        reserveOrThrow_(precursor_id_builder_.Reserve(expectedMobilograms), "PRECURSOR_ID");
        reserveOrThrow_(ident_id_builder_.Reserve(expectedMobilograms), "IDENTIFICATION_ID");
        reserveOrThrow_(modified_sequence_builder_.Reserve(expectedMobilograms), "MODIFIED_SEQUENCE");
        reserveOrThrow_(precursor_charge_builder_.Reserve(expectedMobilograms), "PRECURSOR_CHARGE");
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
      if (wrote_) return;
      write_();
#endif
    }

  private:
    String filename_;
    UInt64 run_id_{0};
    String source_file_;
    bool wrote_{false};

#ifdef WITH_PARQUET
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
        arrow::field("RT", arrow::float64()),
        arrow::field("PRECURSOR_ID", arrow::int64()),
        arrow::field("IDENTIFICATION_ID", arrow::int64()),
        arrow::field("MODIFIED_SEQUENCE", arrow::utf8()),
        arrow::field("PRECURSOR_CHARGE", arrow::int64()),
        arrow::field("MOBILITY_DATA", arrow::binary()),
        arrow::field("INTENSITY_DATA", arrow::binary()),
        arrow::field("MOBILITY_COMPRESSION", arrow::int64()),
        arrow::field("INTENSITY_COMPRESSION", arrow::int64())
      });

      auto table = arrow::Table::Make(schema, {
        finishArray_(run_id_builder_, "RUN_ID"),
        finishArray_(source_file_builder_, "SOURCE_FILE"),
        finishArray_(rt_builder_, "RT"),
        finishArray_(precursor_id_builder_, "PRECURSOR_ID"),
        finishArray_(ident_id_builder_, "IDENTIFICATION_ID"),
        finishArray_(modified_sequence_builder_, "MODIFIED_SEQUENCE"),
        finishArray_(precursor_charge_builder_, "PRECURSOR_CHARGE"),
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
    arrow::DoubleBuilder rt_builder_;
    arrow::Int64Builder precursor_id_builder_;
    arrow::Int64Builder ident_id_builder_;
    arrow::StringBuilder modified_sequence_builder_;
    arrow::Int64Builder precursor_charge_builder_;
    arrow::BinaryBuilder mobility_data_builder_;
    arrow::BinaryBuilder intensity_data_builder_;
    arrow::Int64Builder mobility_compression_builder_;
    arrow::Int64Builder intensity_compression_builder_;
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

  void MobilogramParquetConsumer::consumeMobilogram(Mobilogram& m)
  {
    impl_->consumeMobilogram(m);
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
