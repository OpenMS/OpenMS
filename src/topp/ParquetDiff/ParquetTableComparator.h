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
  struct ParquetDiffSettings
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
    /**
      Compare (and dump) the QPX identity and cross-reference columns as ordinary values.

      Off by default, and the default is the correct one for a reference comparison. The spec is
      explicit that "identity is meaningful within a file only; distinct QPX files must not be
      joined on @c feature_id alone" -- and comparing a file against a committed reference is
      exactly that join. The values are also hostile to it: @c feature_id hashes @c rt and
      @c observed_mz, so one ULP of platform drift does not shift a digit, it produces a wholly
      different int64 and turns a tolerable numeric difference into "row missing on both sides".

      What is asserted about the ids instead is stronger and portable: @c -schema mode checks that
      they are present, non-null and unique, and QPXIdentity_test pins the derivation itself
      against the reference implementation's own vectors.
    */
    bool compare_identity_columns = false;
  };

  /**
    @brief Outcome of a table comparison

    The three message vectors are populated independently, so a run can report a schema
    mismatch and a value mismatch at once. @p equal is true only when all three are empty.
  */
  struct ParquetDiffResult
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
  class ParquetTableComparator
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
      @brief Write a Parquet table out as TSV, one row per line, sorted by primary key.

      A Parquet file cannot be reviewed in a diff or patched by hand, so a committed binary
      reference can only ever be regenerated wholesale - which is exactly the operation that hides
      unrelated drift. Dumping to text puts Parquet output back on the same footing as every other
      reference in the suite: readable in a pull request, comparable with @ref TOPP_FuzzyDiff, and
      editable line by line.

      Rows are emitted in primary-key order rather than file order, so the dump does not depend on
      the order the producer happened to write - the same property that makes the comparison
      order-insensitive. List- and struct-valued cells are rendered in full; nulls print as
      @c null, which is distinct from an empty string or a zero.

      @param[in] file the Parquet file to read
      @param[in] out_file destination TSV path
      @param[in] settings only ParquetDiffSettings::primary_key is used; when empty the key is
                 derived from the file's QPX @c file_type metadata, and failing that rows are
                 emitted in file order
      @return true when the file was read and written successfully
    */
    static bool dumpToTsv(const std::string& file,
                          const std::string& out_file,
                          const ParquetDiffSettings& settings);

    /**
      @brief The primary key of a QPX view.

      @param[in] view one of @c psm, @c feature, @c pg
      @return the key columns, in order
      @exception Exception::IllegalArgument if @p view is not a known view name
    */
    static std::vector<std::string> qpxPrimaryKey(const std::string& view);

    /**
      @brief The opaque identity and cross-reference columns of a QPX view

      @c feature_id / @c psm_id / @c pg_id, plus the optional references between the views. They
      are derived from the other columns rather than measured, so a difference in one is always a
      difference in something else that is reported anyway -- see
      ParquetDiffSettings::compare_identity_columns for why comparing them across files is wrong
      rather than merely redundant.

      @param[in] view @c "psm", @c "feature" or @c "pg"
      @return The column names, or empty for an unknown view
    */
    static std::vector<std::string> qpxIdentityColumns(const std::string& view);

    /**
      @brief Map a QPX @c file_type metadata value to a view name.

      @param[in] file_type e.g. @c psm_file, @c feature_file, @c pg_file
      @return the view name, or an empty string when @p file_type is not recognised
    */
    static std::string viewFromFileType(const std::string& file_type);
  };

} // namespace OpenMS
