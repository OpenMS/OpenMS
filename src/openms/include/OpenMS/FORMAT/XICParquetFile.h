// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Reader for OpenSWATH chromatogram Parquet files (.xic).
  */
  class OPENMS_DLLAPI XICParquetFile
  {
  public:
    struct XICChromatogram
    {
      Int64 run_id{0};
      String source_file;
      Int64 ms_level{0};

      bool has_precursor_id{false};
      Int64 precursor_id{0};
      bool has_transition_id{false};
      Int64 transition_id{0};
      String modified_sequence;
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
      String transition_type;
      String annotation;

      std::vector<double> rt;
      std::vector<double> intensity;
    };

    explicit XICParquetFile(const String& filename);
    explicit XICParquetFile(const std::vector<String>& filenames);
    XICParquetFile(const XICParquetFile& rhs) = default;
    XICParquetFile& operator=(const XICParquetFile& rhs) = default;

    const String& getFilename() const;
    const std::vector<String>& getFilenames() const;

    /// Load all chromatograms from the file.
    void load(std::vector<XICChromatogram>& output) const;

    /**
      @brief Load chromatograms with optional filtering.

      @param[out] output Output chromatograms
      @param[in] precursor_id Optional precursor id (-1 to ignore)
      @param[in] transition_id Optional transition id (-1 to ignore)
      @param[in] modified_sequence Optional sequence filter (empty to ignore)
      @param[in] precursor_charge Optional charge filter (-1 to ignore)
      @param[in] product_charge Optional product charge filter (-1 to ignore)
      @param[in] ms_level Optional MS level filter (-1 to ignore)
      @param[in] run_id Optional run_id filter (-1 to ignore)
      @param[in] filter Optional filter expression on columns (e.g., "PRECURSOR_ID=1 OR TRANSITION_ID in [2,3]")
    */
    void getChromatograms(std::vector<XICChromatogram>& output,
                          Int64 precursor_id = -1,
                          Int64 transition_id = -1,
                          const String& modified_sequence = "",
                          Int64 precursor_charge = -1,
                          Int64 product_charge = -1,
                          Int64 ms_level = -1,
                          Int64 run_id = -1,
                          const String& filter = "") const;

  private:
    String filename_;
    std::vector<String> filenames_;
  };

  /// Convenience alias for the nested XIC chromatogram type.
  typedef XICParquetFile::XICChromatogram XICChromatogram;
} // namespace OpenMS
