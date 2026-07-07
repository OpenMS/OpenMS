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

    @note The .xipm schema is defined by XIPMParquetConsumer.
    @see OpenMS::XIPMParquetConsumer
  */
  class OPENMS_DLLAPI XIPMParquetFile
  {
  public:
    /**
      @brief Extracted targeted mz/RT/IM peak map, as stored in one row of a .xipm file.
    */
    struct OPENMS_DLLAPI XIPMPeakMap
    {
      Int64 run_id{0}; ///< Run identifier
      std::string source_file; ///< Input source filename
      Int64 ms_level{0}; ///< MS level (1 for precursor peak maps, 2 for fragment peak maps)
      std::string peakmap_type; ///< Peak-map type ("precursor" or "transition")

      bool has_precursor_id{false}; ///< Whether @ref precursor_id is set
      Int64 precursor_id{0}; ///< Precursor id
      bool has_transition_id{false}; ///< Whether @ref transition_id is set
      Int64 transition_id{0}; ///< Transition id
      std::string modified_sequence; ///< Modified peptide sequence
      bool has_precursor_charge{false}; ///< Whether @ref precursor_charge is set
      Int64 precursor_charge{0}; ///< Precursor charge
      bool has_product_charge{false}; ///< Whether @ref product_charge is set
      Int64 product_charge{0}; ///< Product charge
      bool has_detecting_transition{false}; ///< Whether @ref detecting_transition is set
      Int64 detecting_transition{0}; ///< Detecting transition flag
      bool has_precursor_decoy{false}; ///< Whether @ref precursor_decoy is set
      Int64 precursor_decoy{0}; ///< Precursor decoy flag
      bool has_product_decoy{false}; ///< Whether @ref product_decoy is set
      Int64 product_decoy{0}; ///< Product decoy flag
      bool has_transition_ordinal{false}; ///< Whether @ref transition_ordinal is set
      Int64 transition_ordinal{0}; ///< Transition ordinal
      std::string transition_type; ///< Transition type (e.g. y, b)
      std::string annotation; ///< Transition annotation (e.g. y3^1) or precursor annotation

      double target_mz{0.0}; ///< Extraction target m/z
      bool has_target_rt{false}; ///< Whether @ref target_rt is set
      double target_rt{0.0}; ///< Extraction target RT in the run coordinate system
      bool has_target_ion_mobility{false}; ///< Whether @ref target_ion_mobility is set
      double target_ion_mobility{0.0}; ///< Extraction target ion mobility
      bool has_rt_start{false}; ///< Whether @ref rt_start is set
      double rt_start{0.0}; ///< Extraction window start in RT
      bool has_rt_end{false}; ///< Whether @ref rt_end is set
      double rt_end{0.0}; ///< Extraction window end in RT

      std::vector<double> mz; ///< Extracted m/z array
      std::vector<double> rt; ///< Extracted RT array
      std::vector<double> ion_mobility; ///< Extracted ion mobility array
      std::vector<double> intensity; ///< Extracted intensity array
    };

    /**
      @brief Unique run information (run_id, source_file).
    */
    struct OPENMS_DLLAPI XIPMRunInfo
    {
      Int64 run_id{0}; ///< Run identifier
      std::string source_file; ///< Input source filename
    };

    /**
      @brief Construct from a single .xipm file.

      @param[in] filename Path to a peak-map parquet file.
    */
    explicit XIPMParquetFile(const std::string& filename);

    /**
      @brief Construct from multiple .xipm files.

      @param[in] filenames Paths to peak-map parquet files.
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

      This method never decodes mz/RT/ion mobility/intensity arrays and always
      returns distinct rows.

      @param[out] output Output run metadata.
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
