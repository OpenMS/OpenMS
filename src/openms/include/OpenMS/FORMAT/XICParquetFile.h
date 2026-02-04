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
    @brief Reader for OpenSWATH chromatogram Parquet files (.xic).

    Supports loading single or multiple files and filtering on metadata
    columns (e.g., precursor id, transition id, annotations). Filters are
    applied before decoding RT/intensity binary arrays.

    @section XICParquetFile_FilterSyntax Filter syntax
    The @p filter argument in getChromatograms() accepts simple boolean
    expressions over column names. Supported operators are:
    - Comparison: =, ==, !=, <, <=, >, >=
    - Set membership: in [v1, v2, ...]
    - Boolean: AND/OR (also accepts &&, ||, &, |)

    Values can be integers or strings; strings may be unquoted if they contain
    no spaces or commas (e.g., annotation=y3^1), otherwise use quotes.

    Supported filter columns (case-insensitive):
    RUN_ID, SOURCE_FILE, MS_LEVEL, PRECURSOR_ID, TRANSITION_ID,
    MODIFIED_SEQUENCE, PRECURSOR_CHARGE, PRODUCT_CHARGE, DETECTING_TRANSITION,
    PRECURSOR_DECOY, PRODUCT_DECOY, TRANSITION_ORDINAL, TRANSITION_TYPE,
    ANNOTATION. RT and INTENSITY are not filterable because they are stored
    as compressed binary arrays.

    @section XICParquetFile_Internal Internal processing notes
    The implementation uses an Arrow-based pipeline:
    - If Arrow Dataset is available, filters are translated into Arrow
      expressions and pushed down via dataset scanning.
    - If dataset filtering is unavailable or fails, the same filter expression
      is evaluated in-memory using Arrow compute.
    - RT/intensity binary arrays are decoded only after filtering.

    These steps are implemented in helper functions in the corresponding
    `.cpp` file (e.g., dataset scan vs. compute filter fallback and filter
    parsing). Keeping the helpers in the implementation file avoids exposing
    Arrow types in the public header.

    @note The .xic schema is defined by MSChromatogramParquetConsumer.
    @see OpenMS::MSChromatogramParquetConsumer
  */
  class OPENMS_DLLAPI XICParquetFile
  {
  public:
    /**
      @brief Lightweight chromatogram container for XIC parquet rows.
    */
    struct OPENMS_DLLAPI XICChromatogram
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

    /**
      @brief Unique run information (run_id, source_file).
    */
    struct OPENMS_DLLAPI XICRunInfo
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
    struct OPENMS_DLLAPI XICAnalyte
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
      @brief Construct from a single .xic file.

      @param[in] filename Path to an OpenSWATH chromatogram parquet file.
    */
    explicit XICParquetFile(const String& filename);

    /**
      @brief Construct from multiple .xic files.

      @param[in] filenames Paths to OpenSWATH chromatogram parquet files.
    */
    explicit XICParquetFile(const std::vector<String>& filenames);
    XICParquetFile(const XICParquetFile& rhs) = default;
    XICParquetFile& operator=(const XICParquetFile& rhs) = default;

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
      @brief Load all chromatograms from the file(s).

      @param[out] output Output chromatograms.
    */
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
      @param[in] filter Optional filter expression on columns
        (e.g., "PRECURSOR_ID=1 OR TRANSITION_ID in [2,3]")
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

    /**
      @brief Return chromatograms using a typed filter expression.

      @param[out] output Output chromatograms
      @param[in] filter Typed filter builder expression
    */
    void getChromatograms(std::vector<XICChromatogram>& output,
                          const ParquetFilter& filter) const;

    /**
      @brief Return chromatograms using a typed filter builder.

      @param[out] output Output chromatograms
      @param[in] filter Typed filter builder
    */
    void getChromatograms(std::vector<XICChromatogram>& output,
                          const ParquetFilterBuilder& filter) const;

    /**
      @brief Return unique run metadata (run_id, source_file).

      This method never decodes RT/intensity arrays and always returns distinct
      rows.
    */
    void getRuns(std::vector<XICRunInfo>& output) const;

    /**
      @brief Return unique analyte metadata.

      If @p nest_transitions is false, each row represents a unique
      precursor-transition pair. If @p nest_transitions is true, each row
      represents a unique precursor with transition-level fields aggregated
      into vectors.

      This method never decodes RT/intensity arrays and always returns distinct
      entries.

      @param[out] output Output analyte metadata
      @param[in] columns Optional list of analyte columns to return (empty for defaults)
      @param[in] nest_transitions Aggregate transition fields per precursor
    */
    void getAnalytes(std::vector<XICAnalyte>& output,
                     const std::vector<String>& columns = {},
                     bool nest_transitions = true) const;

    /**
      @brief Return the parquet schema column names.

      @param[out] output Column names.
    */
    void getColumns(std::vector<String>& output) const;

  private:
    void getChromatograms_(std::vector<XICChromatogram>& output,
                           const FilterExpression& extra_filter,
                           Int64 precursor_id,
                           Int64 transition_id,
                           const String& modified_sequence,
                           Int64 precursor_charge,
                           Int64 product_charge,
                           Int64 ms_level,
                           Int64 run_id,
                           const String& filter) const;

    String filename_;
    std::vector<String> filenames_;
  };

  /// Convenience alias for the nested XIC chromatogram type.
  typedef XICParquetFile::XICChromatogram XICChromatogram;
  /// Convenience alias for the nested run info type.
  typedef XICParquetFile::XICRunInfo XICRunInfo;
  /// Convenience alias for the nested analyte type.
  typedef XICParquetFile::XICAnalyte XICAnalyte;
} // namespace OpenMS
