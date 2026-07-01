// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief Reader for extracted ion peak-map Parquet files (.xipm).

    Supports loading single or multiple files, either the full set of peak maps
    (load()) or a subset restricted by metadata columns such as precursor id,
    transition id, or MS level (getPeakMaps()).

    @note The .xipm schema is defined by XIPMSchema / XIPMParquetConsumer.
    @see OpenMS::XIPMParquetConsumer
  */
  class OPENMS_DLLAPI XIPMParquetFile
  {
  public:
    /**
      @brief Targeted raw peak map (mz/RT/ion-mobility/intensity arrays) for one analyte.
    */
    struct OPENMS_DLLAPI XIPMPeakMap
    {
      Int64 run_id{0};
      std::string source_file;
      Int64 ms_level{0};
      std::string peakmap_type;

      bool has_precursor_id{false};
      Int64 precursor_id{0};
      bool has_transition_id{false};
      Int64 transition_id{0};
      std::string modified_sequence;
      bool has_precursor_charge{false};
      Int64 precursor_charge{0};
      bool has_product_charge{false};
      Int64 product_charge{0};
      bool has_detecting_transition{false};
      Int64 detecting_transition{0};
      bool has_precursor_decoy{false};
      Int64 precursor_decoy{0};
      bool has_product_decoy{false};
      Int64 product_decoy{0};
      bool has_transition_ordinal{false};
      Int64 transition_ordinal{0};
      std::string transition_type;
      std::string annotation;

      double target_mz{0.0};
      bool has_target_rt{false};
      double target_rt{0.0};
      bool has_target_ion_mobility{false};
      double target_ion_mobility{0.0};
      bool has_rt_start{false};
      double rt_start{0.0};
      bool has_rt_end{false};
      double rt_end{0.0};

      std::vector<double> mz;
      std::vector<double> rt;
      std::vector<double> ion_mobility;
      std::vector<double> intensity;
    };

    /**
      @brief Unique run information (run_id, source_file).
    */
    struct OPENMS_DLLAPI XIPMRunInfo
    {
      Int64 run_id{0};
      std::string source_file;
    };

    /**
      @brief Construct from a single .xipm file.

      @param[in] filename Path to an OpenSWATH peak-map parquet file.
    */
    explicit XIPMParquetFile(const std::string& filename);

    /**
      @brief Construct from multiple .xipm files.

      @param[in] filenames Paths to OpenSWATH peak-map parquet files.
    */
    explicit XIPMParquetFile(const std::vector<std::string>& filenames);
    XIPMParquetFile(const XIPMParquetFile& rhs) = default;
    XIPMParquetFile& operator=(const XIPMParquetFile& rhs) = default;

    /**
      @brief Return the primary filename.

      For multi-file instances this is the first file in the list.

      @return Primary filename.
    */
    const std::string& getFilename() const;

    /**
      @brief Return all filenames associated with this instance.

      @return All filenames associated with this instance.
    */
    const std::vector<std::string>& getFilenames() const;

    /**
      @brief Load all peak maps from the file(s).

      @param[out] output Output peak maps.
    */
    void load(std::vector<XIPMPeakMap>& output) const;

    /**
      @brief Load peak maps with optional filtering.

      @param[out] output Output peak maps
      @param[in] precursor_id Optional precursor id filter (-1 to ignore)
      @param[in] transition_id Optional transition id filter (-1 to ignore)
      @param[in] modified_sequence Optional sequence filter (empty to ignore)
      @param[in] precursor_charge Optional precursor charge filter (-1 to ignore)
      @param[in] product_charge Optional product charge filter (-1 to ignore)
      @param[in] ms_level Optional MS level filter (-1 to ignore)
      @param[in] run_id Optional run_id filter (-1 to ignore)
      @param[in] peakmap_type Optional peak-map type filter (empty to ignore)
    */
    void getPeakMaps(std::vector<XIPMPeakMap>& output,
                     Int64 precursor_id = -1,
                     Int64 transition_id = -1,
                     const std::string& modified_sequence = "",
                     Int64 precursor_charge = -1,
                     Int64 product_charge = -1,
                     Int64 ms_level = -1,
                     Int64 run_id = -1,
                     const std::string& peakmap_type = "") const;

    /**
      @brief Return unique run metadata (run_id, source_file).
    */
    void getRuns(std::vector<XIPMRunInfo>& output) const;

    /**
      @brief Return the parquet schema column names.

      @param[out] output Column names.
    */
    void getColumns(std::vector<std::string>& output) const;

  private:
    std::string filename_;
    std::vector<std::string> filenames_;
  };
} // namespace OpenMS
