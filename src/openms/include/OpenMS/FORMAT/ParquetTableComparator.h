// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief Settings for ParquetTableComparator

    The tolerance pair mirrors FuzzyStringComparator: a numeric cell is accepted when
    <em>either</em> @p acceptable_absdiff or @p acceptable_ratio is satisfied. Use
    @p acceptable_absdiff for the "zero vs. epsilon" case, where no ratio is meaningful.
  */
  struct OPENMS_DLLAPI ParquetDiffSettings
  {
    /// Columns forming the primary key. Empty: derive from the file's QPX @c file_type metadata.
    std::vector<std::string> primary_key;
    /// Columns excluded from value comparison (they are still schema-checked).
    std::vector<std::string> ignore_columns;
    /// Maximum acceptable ratio max(|a|,|b|)/min(|a|,|b|) of two numbers. 1.0 = must be identical.
    double acceptable_ratio = 1.0;
    /// Maximum acceptable |a-b| of two numbers.
    double acceptable_absdiff = 0.0;
    /// Compare schemas only; skip key building and value comparison.
    bool schema_only = false;
    /// Compare list-valued cells as multisets rather than sequences.
    bool unordered_lists = false;
    /// Stop collecting messages of any one category after this many (0 = unlimited).
    Size max_reported = 25;
  };

  /**
    @brief Outcome of a table comparison

    The three message vectors are populated independently, so a run can report a schema
    mismatch and a value mismatch at once. @p equal is true only when all three are empty.
  */
  struct OPENMS_DLLAPI ParquetDiffResult
  {
    /// True when no schema, key or value difference was found.
    bool equal = false;
    /// Missing/extra columns, differing Arrow types, differing nullability.
    std::vector<std::string> schema_errors;
    /// Duplicate primary keys, and keys present in only one of the two tables.
    std::vector<std::string> key_errors;
    /// Per-cell differences, prefixed by the primary key of the row they occur in.
    std::vector<std::string> value_errors;
    /// Number of messages suppressed by ParquetDiffSettings::max_reported.
    Size suppressed = 0;

    Size rows_1 = 0;      ///< row count of the first table
    Size rows_2 = 0;      ///< row count of the second table
    Size rows_compared = 0;  ///< rows matched by primary key and compared cell-by-cell

    /// Largest ratio observed between two numeric cells (1.0 when all were identical).
    double max_ratio = 1.0;
    /// Largest absolute difference observed between two numeric cells.
    double max_absdiff = 0.0;

    /// The primary key actually used, whether supplied or auto-detected.
    std::vector<std::string> primary_key_used;
  };

  /**
    @brief Primary-key-aware, order-insensitive, tolerant comparison of Parquet tables

    Two Parquet files holding the same logical table may legitimately differ in row order and
    in the low-order bits of floating-point columns, so neither a byte comparison nor a text
    diff of a serialised form answers "are these the same result?". This class matches rows by
    a primary key, compares the matched rows cell-by-cell with a numeric tolerance, and reports
    schema drift separately from value drift.

    It is the table-shaped counterpart of @ref TOPP_FuzzyDiff, and is exposed as
    @ref TOPP_ParquetDiff.

    @note All Arrow/Parquet API use is confined to the implementation, so binaries linking
          libOpenMS do not need to import Arrow symbols (see ArrowIOHelpers.h).

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ParquetTableComparator
  {
  public:
    /**
      @brief Compare two Parquet files.

      Rows are matched on ParquetDiffSettings::primary_key. List-valued key columns are
      canonicalised by sorting their elements before the key is formed, matching how QPX
      treats set-valued key columns such as @c grouped_runs.

      @param[in] file_1 first Parquet file
      @param[in] file_2 second Parquet file
      @param[in] settings tolerance, key and reporting options
      @return the differences found; ParquetDiffResult::equal is true when there are none
    */
    static ParquetDiffResult compare(const std::string& file_1,
                                     const std::string& file_2,
                                     const ParquetDiffSettings& settings);

    /**
      @brief Check one Parquet file against a built-in QPX schema.

      Reports columns the view requires but the file lacks, columns whose Arrow type differs
      from the schema, and columns the file carries that the view does not define. Duplicate
      primary keys are reported as well, since a QPX primary key must be unique.

      @param[in] file the Parquet file to check
      @param[in] view one of @c psm, @c feature, @c pg
      @param[in] settings reporting options (tolerances are unused)
      @return the violations found; ParquetDiffResult::equal is true when there are none
      @exception Exception::IllegalArgument if @p view is not a known view name
    */
    static ParquetDiffResult validate(const std::string& file,
                                      const std::string& view,
                                      const ParquetDiffSettings& settings);

    /**
      @brief The primary key of a QPX view.

      @param[in] view one of @c psm, @c feature, @c pg
      @return the key columns, in order
      @exception Exception::IllegalArgument if @p view is not a known view name
    */
    static std::vector<std::string> qpxPrimaryKey(const std::string& view);

    /**
      @brief Map a QPX @c file_type metadata value to a view name.

      @param[in] file_type e.g. @c psm_file, @c feature_file, @c pg_file
      @return the view name, or an empty string when @p file_type is not recognised
    */
    static std::string viewFromFileType(const std::string& file_type);
  };

} // namespace OpenMS
