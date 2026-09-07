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
  */
  class OPENMS_DLLAPI XIPMParquetFile
  {
  public:
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

    struct OPENMS_DLLAPI XIPMRunInfo
    {
      Int64 run_id{0};
      std::string source_file;
    };

    explicit XIPMParquetFile(const std::string& filename);
    explicit XIPMParquetFile(const std::vector<std::string>& filenames);
    XIPMParquetFile(const XIPMParquetFile& rhs) = default;
    XIPMParquetFile& operator=(const XIPMParquetFile& rhs) = default;

    const std::string& getFilename() const;
    const std::vector<std::string>& getFilenames() const;

    void load(std::vector<XIPMPeakMap>& output) const;

    void getPeakMaps(std::vector<XIPMPeakMap>& output,
                     Int64 precursor_id = -1,
                     Int64 transition_id = -1,
                     const std::string& modified_sequence = "",
                     Int64 precursor_charge = -1,
                     Int64 product_charge = -1,
                     Int64 ms_level = -1,
                     Int64 run_id = -1,
                     const std::string& peakmap_type = "") const;

    void getRuns(std::vector<XIPMRunInfo>& output) const;
    void getColumns(std::vector<std::string>& output) const;

  private:
    std::string filename_;
    std::vector<std::string> filenames_;
  };
} // namespace OpenMS
