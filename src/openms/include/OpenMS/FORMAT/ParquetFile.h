// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/arrow/reader.h>

#include <memory>
#include <string>
#include <vector>

namespace OpenMS
{

  /**
      @brief Shared utilities for reading, writing, and packaging Parquet-based file formats.

      This class provides static helpers used by multiple OpenMS Parquet-backed
      I/O classes (e.g. TransitionParquetFile, OpenSwathOSWParquetWriter,
      XICParquetFile).

      Capabilities include:
      - Reading / writing Arrow Tables to individual Parquet files
      - Type-safe column accessors with flexible type coercion
      - Zip / unzip of Parquet directory bundles (e.g. .oswpq archives)
      - Common Arrow builder helpers (append, finish)

      All Parquet zip archives use store-only compression (`-0`) because Parquet
      files are already internally compressed; re-compressing with deflate wastes
      CPU for negligible size reduction.

      @ingroup FileIO
  */
  class OPENMS_DLLAPI ParquetFile
  {
  public:
    /** @name Arrow builder helpers
    */
    //@{

    /**
      @brief Append a value to an Arrow builder, throwing on failure.

      @param[in,out] status  Arrow status returned by the append call
      @param[in] column      Column name (used in error messages)

      @throws Exception::InvalidValue if the status is not OK
    */
    static void appendOrThrow(const arrow::Status& status, const std::string& column);

    /**
      @brief Finish an Arrow builder and return the resulting Array.

      @param[in,out] builder  Any Arrow ArrayBuilder subclass
      @param[in] name         Column name (used in error messages)

      @return Shared pointer to the finished Arrow Array

      @throws Exception::InvalidValue if finishing fails
    */
    template <typename BuilderT>
    static std::shared_ptr<arrow::Array> finishArray(BuilderT& builder, const std::string& name)
    {
      std::shared_ptr<arrow::Array> array;
      auto status = builder.Finish(&array);
      if (!status.ok())
      {
        throw_finish_error_(name, status.ToString());
      }
      return array;
    }

    //@}

    /** @name Parquet file I/O
    */
    //@{

    /**
      @brief Write an Arrow Table to a Parquet file.

      @param[in] table     The Arrow Table to write
      @param[in] filename  Output file path
      @param[in] row_group_size  Number of rows per row group (default: 262144)

      @throws Exception::FileNotWritable if the file cannot be opened
      @throws Exception::InvalidValue if writing fails
    */
    static void writeTable(const std::shared_ptr<arrow::Table>& table, const String& filename,
                int64_t row_group_size = 262144);

    /**
      @brief Read a Parquet file into an Arrow Table.

      The table is returned with chunks combined into single arrays per column.

      @param[in] filename  Input Parquet file path

      @return Shared pointer to the Arrow Table

      @throws Exception::InvalidValue if reading fails
    */
    static std::shared_ptr<arrow::Table> readTable(const String& filename);
    /**
      @brief Read a Parquet file from an Arrow RandomAccessFile into an Arrow Table.

      Allows reading Parquet data directly from an in-archive RandomAccessFile (e.g. libzip-backed).
    */
    static std::shared_ptr<arrow::Table> readTable(const std::shared_ptr<arrow::io::RandomAccessFile>& infile);

    //@}

    /** @name Column accessors
    */
    //@{

    /**
      @brief Get a required column from an Arrow Table by name.

      @param[in] table  The Arrow Table
      @param[in] name   Column name

      @return The first chunk of the column as an Arrow Array

      @throws Exception::MissingInformation if column not found
      @throws Exception::InvalidValue if column has no chunks
    */
    static std::shared_ptr<arrow::Array> getColumn(const std::shared_ptr<arrow::Table>& table,
                                                   const std::string& name);

    /**
      @brief Get an optional column from an Arrow Table by name.

      @param[in] table  The Arrow Table
      @param[in] name   Column name

      @return The first chunk of the column, or nullptr if not present

      @throws Exception::InvalidValue if column exists but has no chunks
    */
    static std::shared_ptr<arrow::Array> getOptionalColumn(const std::shared_ptr<arrow::Table>& table,
                                                           const std::string& name);

    //@}

    /** @name Type-safe value accessors
    */
    //@{

    /**
      @brief Read an integer value from an Arrow Array with type coercion.

      Supports Int8–Int64, UInt8–UInt64 types.

      @param[in] array         The Arrow Array (may be nullptr if allow_null)
      @param[in] row           Row index
      @param[in] default_value Value returned when array is nullptr or value is null
      @param[in] allow_null    If true, null values return default_value instead of throwing

      @return The integer value

      @throws Exception::MissingInformation if value is null and allow_null is false
      @throws Exception::InvalidValue for unsupported column types
    */
    static int64_t getInt64(const std::shared_ptr<arrow::Array>& array, int64_t row,
                            int64_t default_value, bool allow_null);

    /**
      @brief Read a floating-point value from an Arrow Array with type coercion.

      Supports Float, Double, and all integer types (coerced to double).

      @param[in] array         The Arrow Array (may be nullptr if allow_null)
      @param[in] row           Row index
      @param[in] default_value Value returned when array is nullptr or value is null
      @param[in] allow_null    If true, null values return default_value instead of throwing

      @return The double value

      @throws Exception::MissingInformation if value is null and allow_null is false
      @throws Exception::InvalidValue for unsupported column types
    */
    static double getDouble(const std::shared_ptr<arrow::Array>& array, int64_t row,
                            double default_value, bool allow_null);

    /**
      @brief Read a boolean value from an Arrow Array with type coercion.

      Supports Boolean and all integer types (non-zero = true).

      @param[in] array         The Arrow Array (may be nullptr if allow_null)
      @param[in] row           Row index
      @param[in] default_value Value returned when array is nullptr or value is null
      @param[in] allow_null    If true, null values return default_value instead of throwing

      @return The boolean value

      @throws Exception::MissingInformation if value is null and allow_null is false
      @throws Exception::InvalidValue for unsupported column types
    */
    static bool getBool(const std::shared_ptr<arrow::Array>& array, int64_t row,
                        bool default_value, bool allow_null);

    /**
      @brief Read a string value from an Arrow Array.

      Supports String and LargeString types.

      @param[in] array  The Arrow Array (may be nullptr)
      @param[in] row    Row index

      @return The string value, or empty string if array is nullptr or value is null

      @throws Exception::InvalidValue for unsupported column types
    */
    static std::string getString(const std::shared_ptr<arrow::Array>& array, int64_t row);

    /**
      @brief Read a list of strings from an Arrow Array.

      Supports String/LargeString (semicolon-delimited) and List/LargeList of strings.

      @param[in] array  The Arrow Array (may be nullptr)
      @param[in] row    Row index

      @return Vector of strings (empty if array is nullptr or value is null)

      @throws Exception::InvalidValue for unsupported column types
    */
    static std::vector<std::string> getStringList(const std::shared_ptr<arrow::Array>& array, int64_t row);

    //@}

    /** @name Misc helpers
    */
    //@{

    /**
      @brief Escape a string for safe embedding into JSON values.

      Mirrors the ad-hoc jsonEscape_ implementations used in several
      Parquet-writing sources so callers can reuse a single canonical
      implementation.
    */
    static std::string jsonEscape(const String& input);

    /**
      @brief Return the number of rows in a parquet file using the low-level
      parquet reader metadata. Returns 0 if the file does not exist.
    */
    static int64_t rowCount(const String& filename);

    //@}


  private:
    /// Internal helper to throw a consistent error from finishArray
    static void throw_finish_error_(const std::string& name, const std::string& error);
  };

} // namespace OpenMS
