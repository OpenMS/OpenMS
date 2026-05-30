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
#include <OpenMS/FORMAT/ParquetFilter.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Reader for OpenSWATH mobilogram Parquet files (.xim).

    Supports loading single or multiple files and filtering on metadata
    columns (e.g., precursor id, transition id, annotations). Filters are
    applied before decoding mobility/intensity binary arrays.

    @section XIMParquetFile_FilterSyntax Filter syntax
    The @p filter argument in getMobilograms() accepts simple boolean
    expressions over column names. Supported operators are:
    - Comparison: =, ==, !=, <, <=, >, >=
    - Set membership: in [v1, v2, ...]
    - Boolean: AND/OR (also accepts &&, ||, &, |)

    Values can be integers or strings; strings may be unquoted if they contain
    no spaces or commas (e.g., annotation=y3^1), otherwise use quotes.

    Supported filter columns (case-insensitive):
    RUN_ID, SOURCE_FILE, MS_LEVEL, MOBILOGRAM_TYPE, PRECURSOR_ID, TRANSITION_ID,
    FEATURE_ID, FEATURE_RT,
    MODIFIED_SEQUENCE, PRECURSOR_CHARGE, PRODUCT_CHARGE, DETECTING_TRANSITION,
    PRECURSOR_DECOY, PRODUCT_DECOY, TRANSITION_ORDINAL, TRANSITION_TYPE,
    ANNOTATION. MOBILITY and INTENSITY are not filterable because they are stored
    as compressed binary arrays.

    @section XIMParquetFile_Internal Internal processing notes
    The implementation uses an Arrow-based pipeline:
    - If Arrow Dataset is available, filters are translated into Arrow
      expressions and pushed down via dataset scanning.
    - If dataset filtering is unavailable or fails, the same filter expression
      is evaluated in-memory using Arrow compute.
    - mobility/intensity binary arrays are decoded only after filtering.

    These steps are implemented in helper functions in the corresponding
    `.cpp` file (e.g., dataset scan vs. compute filter fallback and filter
    parsing). Keeping the helpers in the implementation file avoids exposing
    Arrow types in the public header.

    @note The .xim schema is defined by MobilogramParquetConsumer.
    @see OpenMS::MobilogramParquetConsumer
  */
  class OPENMS_DLLAPI XIMParquetFile
  {
  public:
    /**
      @brief Lightweight mobilogram container for XIM parquet rows.
    */
    struct OPENMS_DLLAPI XIMMobilogram
    {
      Int64 run_id{0};
      String source_file;
      Int64 ms_level{0};
      String mobilogram_type;

      bool has_precursor_id{false};
      Int64 precursor_id{0};
      bool has_transition_id{false};
      Int64 transition_id{0};
      bool has_feature_id{false};
      Int64 feature_id{0};
      bool has_feature_rt{false};
      double feature_rt{0.0};
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

      std::vector<double> mobility;
      std::vector<double> intensity;
    };

    /**
      @brief Unique run information (run_id, source_file).
    */
    struct OPENMS_DLLAPI XIMRunInfo
    {
      Int64 run_id{0};
      String source_file;
    };

    /**
      @brief Analyte metadata container.

      If @p nest_transitions is false in getAnalytes(), transition-level fields
      are stored in the scalar members (transition_id, product_charge, etc.).
      If @p nest_transitions is true, transition-level fields are stored in the
      vector members (transition_ids, product_charges, etc.), with one entry
      per unique transition belonging to the precursor.
    */
    struct OPENMS_DLLAPI XIMAnalyte
    {
      bool has_precursor_id{false};
      Int64 precursor_id{0};
      String modified_sequence;
      bool has_precursor_charge{false};
      Int64 precursor_charge{0};
      bool has_precursor_decoy{false};
      Int64 precursor_decoy{0};

      bool has_transition_id{false};
      Int64 transition_id{0};
      bool has_product_charge{false};
      Int64 product_charge{0};
      bool has_transition_ordinal{false};
      Int64 transition_ordinal{0};
      bool has_detecting_transition{false};
      Int64 detecting_transition{0};
      bool has_product_decoy{false};
      Int64 product_decoy{0};
      String transition_type;
      String annotation;

      std::vector<Int64> transition_ids;
      std::vector<Int64> product_charges;
      std::vector<Int64> transition_ordinals;
      std::vector<Int64> detecting_transitions;
      std::vector<Int64> product_decoys;
      std::vector<String> transition_types;
      std::vector<String> annotations;
    };

    /**
      @brief Construct from a single .xim file.

      @param[in] filename Path to an OpenSWATH mobilogram parquet file.
    */
    explicit XIMParquetFile(const String& filename);

    /**
      @brief Construct from multiple .xim files.

      @param[in] filenames Paths to OpenSWATH mobilogram parquet files.
    */
    explicit XIMParquetFile(const std::vector<String>& filenames);
    XIMParquetFile(const XIMParquetFile& rhs) = default;
    XIMParquetFile& operator=(const XIMParquetFile& rhs) = default;

    /**
      @brief Return the primary filename.

      For multi-file instances this is the first file in the list.

      @return Primary filename.
    */
    const String& getFilename() const;

    /**
      @brief Return all filenames associated with this instance.

      @return All filenames associated with this instance.
    */
    const std::vector<String>& getFilenames() const;

    /**
      @brief Load all mobilograms from the file(s).

      @param[out] output Output mobilograms.
    */
    void load(std::vector<XIMMobilogram>& output) const;

    /**
      @brief Load mobilograms with optional filtering.

      @param[out] output Output mobilograms
      @param[in] precursor_id Optional precursor id (-1 to ignore)
      @param[in] transition_id Optional transition id (-1 to ignore)
      @param[in] modified_sequence Optional sequence filter (empty to ignore)
      @param[in] precursor_charge Optional charge filter (-1 to ignore)
      @param[in] product_charge Optional product charge filter (-1 to ignore)
      @param[in] ms_level Optional MS level filter (-1 to ignore)
      @param[in] run_id Optional run_id filter (-1 to ignore)
      @param[in] mobilogram_type Optional mobilogram type filter (empty to ignore)
      @param[in] feature_id Optional feature id filter (-1 to ignore)
      @param[in] feature_rt Optional feature RT filter (< 0 to ignore)
      @param[in] filter Optional filter expression on columns
        (e.g., "PRECURSOR_ID=1 OR TRANSITION_ID in [2,3]")
    */
    void getMobilograms(std::vector<XIMMobilogram>& output,
                          Int64 precursor_id = -1,
                          Int64 transition_id = -1,
                          const String& modified_sequence = "",
                          Int64 precursor_charge = -1,
                          Int64 product_charge = -1,
                          Int64 ms_level = -1,
                          Int64 run_id = -1,
                          const String& mobilogram_type = "",
                          Int64 feature_id = -1,
                          double feature_rt = -1.0,
                          const String& filter = "") const;

    /**
      @brief Return mobilograms using a typed filter expression.

      @param[out] output Output mobilograms
      @param[in] filter Typed filter builder expression
    */
    void getMobilograms(std::vector<XIMMobilogram>& output,
                          const ParquetFilter& filter) const;

    /**
      @brief Return mobilograms using a typed filter builder.

      @param[out] output Output mobilograms
      @param[in] filter Typed filter builder
    */
    void getMobilograms(std::vector<XIMMobilogram>& output,
                          const ParquetFilterBuilder& filter) const;

    /**
      @brief Return unique run metadata (run_id, source_file).

      This method never decodes mobility/intensity arrays and always returns distinct
      rows.
    */
    void getRuns(std::vector<XIMRunInfo>& output) const;

    /**
      @brief Return unique analyte metadata.

      If @p nest_transitions is false, each row represents a unique
      precursor-transition pair. If @p nest_transitions is true, each row
      represents a unique precursor with transition-level fields aggregated
      into vectors.

      This method never decodes mobility/intensity arrays and always returns distinct
      entries.

      @param[out] output Output analyte metadata
      @param[in] columns Optional list of analyte columns to return (empty for defaults)
      @param[in] nest_transitions Aggregate transition fields per precursor
    */
    void getAnalytes(std::vector<XIMAnalyte>& output,
                     const std::vector<String>& columns = {},
                     bool nest_transitions = true) const;

    /**
      @brief Return the parquet schema column names.

      @param[out] output Column names.
    */
    void getColumns(std::vector<String>& output) const;

  private:
    void getMobilograms_(std::vector<XIMMobilogram>& output,
                           const FilterExpression& extra_filter,
                           Int64 precursor_id,
                           Int64 transition_id,
                           const String& modified_sequence,
                           Int64 precursor_charge,
                           Int64 product_charge,
                           Int64 ms_level,
                           Int64 run_id,
                           const String& mobilogram_type,
                           Int64 feature_id,
                           double feature_rt,
                           const String& filter) const;

    String filename_;
    std::vector<String> filenames_;
  };

  /// Convenience alias for the nested XIM mobilogram type.
  typedef XIMParquetFile::XIMMobilogram XIMMobilogram;
  /// Convenience alias for the nested run info type.
  typedef XIMParquetFile::XIMRunInfo XIMRunInfo;
  /// Convenience alias for the nested analyte type.
  typedef XIMParquetFile::XIMAnalyte XIMAnalyte;
} // namespace OpenMS
