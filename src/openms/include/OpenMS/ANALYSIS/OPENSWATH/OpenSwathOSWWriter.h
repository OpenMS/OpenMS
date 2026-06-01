// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#pragma once

// Interfaces
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/SqliteConnector.h>

#include <array>
#include <cstddef>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <mutex>

namespace OpenMS
{

  /**
    @brief Class to write out an OpenSwath OSW SQLite output (PyProphet input).

    The class can take a FeatureMap and create buffered rows from it suitable
    for output to OSW using the prepareRows function. The SQL data is directly
    linked to the PQP file format described in the TransitionPQPFile class. See
    also OpenSwathTSVWriter for another output format.

    The file format has the following tables:

      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>RUN</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ID</td> <td>INT</td> <td> Primary Key (run id)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">FILENAME</td> <td>TEXT</td> <td> Original filename associated with the run </td> </tr>
      </table>

      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>FEATURE</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ID</td> <td>INT</td> <td> Primary Key (feature id)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">RUN_ID</td> <td>INT</td> <td> Foreign Key (RUN.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">PRECURSOR_ID</td> <td>INT</td> <td> Foreign Key (TransitionPQPFile PRECURSOR.ID) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">EXP_RT</td> <td>REAL</td> <td>Experimental RT (retention time) of the feature </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">NORM_RT</td> <td>REAL</td> <td>Normalized RT (retention time) of the feature. The position of the peak group in the normalized retention time space (e.g. fx(RT) where fx describes the transformation fx)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">DELTA_RT</td> <td>REAL</td> <td>The difference in retention between expected retention time of the assay and the measured feature retention time (EXP_RT) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">LEFT_WIDTH</td> <td>REAL</td> <td>Retention time start of the peak (left width) in seconds</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">RIGHT_WIDTH</td> <td>REAL</td> <td>Retention time end of the peak (right width) in seconds</td> </tr>
      </table>

      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>FEATURE_MS1</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">FEATURE_ID</td> <td>INT</td> <td>Foreign Key (FEATURE.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">AREA_INTENSITY</td> <td>REAL</td> <td>%Precursor intensity (area) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">APEX_INTENSITY</td> <td>REAL</td> <td>%Precursor intensity (apex) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">VAR_...</td> <td>REAL</td> <td>%Precursor score used in pyProphet </td> </tr>
      </table>

      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>FEATURE_MS2</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">FEATURE_ID</td> <td>INT</td> <td> Foreign Key (FEATURE.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">AREA_INTENSITY</td> <td>REAL</td> <td>Summed fragment ion intensity (area) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TOTAL_AREA_INTENSITY </td> <td>REAL</td> <td>Summed total XIC of the whole chromatogram </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">APEX_INTENSITY</td> <td>REAL</td> <td>Summed fragment ion intensity (apex) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TOTAL_MI</td> <td>REAL</td> <td>Total mutual information (MI)  </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">VAR_...</td> <td>REAL</td> <td>Fragment ion score used in pyProphet </td> </tr>
      </table>

      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>FEATURE_PRECURSOR</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">FEATURE_ID</td> <td>INT</td> <td> Foreign Key (FEATURE.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">ISOTOPE</td> <td>INT</td> <td>Isotope identifier </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">AREA_INTENSITY</td> <td>REAL</td> <td>%Precursor isotope ion intensity (area) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">APEX_INTENSITY</td> <td>REAL</td> <td>%Precursor isotope ion intensity (apex) </td> </tr>
      </table>

      <table>
        <tr> <th BGCOLOR="#EBEBEB" colspan=3>FEATURE_TRANSITION</th> </tr>
        <tr> <td BGCOLOR="#EBEBEB">FEATURE_ID</td> <td>INT</td> <td> Foreign Key (FEATURE.ID)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TRANSITION_ID</td> <td>INT</td> <td> Foreign Key (transition identifier)</td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">AREA_INTENSITY</td> <td>REAL</td> <td>Fragment ion intensity (area) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TOTAL_AREA_INTENSITY </td> <td>REAL</td> <td>Total XIC of the whole chromatogram </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">APEX_INTENSITY</td> <td>REAL</td> <td>Fragment ion intensity (apex) </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">TOTAL_MI</td> <td>REAL</td> <td>Total mutual information (MI)  </td> </tr>
        <tr> <td BGCOLOR="#EBEBEB">VAR_...</td> <td>REAL</td> <td>Fragment ion score used in pyProphet  </td> </tr>
      </table>

   */
  class OPENMS_DLLAPI OpenSwathOSWWriter
  {
    String output_filename_;
    OpenMS::UInt64 run_id_ = 0;
    bool doWrite_;
    bool enable_uis_scoring_;
    // persistent sqlite connector to avoid reopening DB on each write
    mutable std::unique_ptr<SqliteConnector> conn_;
    mutable std::mutex conn_mutex_;

  public:
    /// Typed OSW cell value used for direct SQLite binding.
    struct OPENMS_DLLAPI OSWValue
    {
      enum class Type : unsigned char
      {
        Null,
        Int64,
        Double,
        Text
      };

      union Storage
      {
        Int64 int_value;
        double double_value;
        String text_value;

        Storage();
        ~Storage();
      };

      Type type;
      Storage storage;

      OSWValue();
      ~OSWValue();
      OSWValue(std::nullptr_t);
      OSWValue(const char* value);
      OSWValue(const String& value);
      OSWValue(const std::string& value);
      OSWValue(Int64 value);
      OSWValue(UInt64 value);
      OSWValue(int value);
      OSWValue(double value);
      OSWValue(const OSWValue& rhs);
      OSWValue(OSWValue&& rhs) noexcept;
      OSWValue& operator=(const OSWValue& rhs);
      OSWValue& operator=(OSWValue&& rhs) noexcept;

      /// Returns a NULL OSW cell.
      static OSWValue null();

      /// Returns a TEXT OSW cell without numeric/null parsing.
      static OSWValue text(const String& value);

      /// Returns true if the cell should be inserted as SQL NULL.
      bool isNull() const;

      /// Returns the integer payload. Only valid when type == Type::Int64.
      Int64 asInt() const;

      /// Returns the double payload. Only valid when type == Type::Double.
      double asDouble() const;

      /// Returns the text payload. Only valid when type == Type::Text.
      const String& asText() const;

    private:
      void destroyText_();
      void copyFrom_(const OSWValue& rhs);
      void moveFrom_(OSWValue&& rhs);
    };

    static constexpr Size FEATURE_COLUMN_COUNT = 11;
    static constexpr Size FEATURE_MS1_COLUMN_COUNT = 18;
    static constexpr Size FEATURE_MS2_COLUMN_COUNT = 38;
    static constexpr Size FEATURE_PRECURSOR_COLUMN_COUNT = 4;
    static constexpr Size FEATURE_TRANSITION_COLUMN_COUNT = 43;

    using FeatureRow = std::array<OSWValue, FEATURE_COLUMN_COUNT>;
    using FeatureMS1Row = std::array<OSWValue, FEATURE_MS1_COLUMN_COUNT>;
    using FeatureMS2Row = std::array<OSWValue, FEATURE_MS2_COLUMN_COUNT>;
    using FeaturePrecursorRow = std::array<OSWValue, FEATURE_PRECURSOR_COLUMN_COUNT>;
    using FeatureTransitionRow = std::array<OSWValue, FEATURE_TRANSITION_COLUMN_COUNT>;

    /// Buffered OSW table rows ready for insertion through prepared statements.
    struct OPENMS_DLLAPI OSWData
    {
      std::vector<FeatureRow> feature_rows;
      std::vector<FeatureMS1Row> feature_ms1_rows;
      std::vector<FeatureMS2Row> feature_ms2_rows;
      std::vector<FeaturePrecursorRow> feature_precursor_rows;
      std::vector<FeatureTransitionRow> feature_transition_rows;

      /// Reserve storage for roughly @p row_count rows per OSW table.
      void reserve(Size row_count);

      /// Reserve storage for expected feature-level and transition-level rows.
      void reserve(Size feature_row_count, Size transition_row_count);

      /// Append rows from @p rhs, moving row storage where possible.
      void append(OSWData&& rhs);

      /// Returns true if no rows were buffered.
      bool empty() const;

      /// Returns the total number of buffered rows across all OSW tables.
      Size rowCount() const;

      /// Estimate the memory retained by this buffer in bytes.
      UInt64 estimateMemoryUsage() const;
    };

    OpenSwathOSWWriter(const String& output_filename,
              bool uis_scores = false);

    bool isActive() const;

    /**
     * @brief Initializes file by generating SQLite tables
     *
     */
    void writeHeader();

    /// Add a RUN entry to the OSW file. Can be called multiple times to register multiple runs.
    void addRun(const UInt64 run_id, const String& input_filename);

    /// Set the current run id used when prepareLine generates FEATURE entries.
    void setRunId(const UInt64 run_id);

    /**
     * @brief Prepare scores for SQLite insertion
     *
     * Some scores might not be defined, those are reported as NULL
     *
     * @param[in] feature The feature being evaluated
     * @param[in] score_name The name of the queried score
     *
     * @returns A string with the queried score
     *
     */
    String getScore(const Feature& feature, const std::string& score_name) const;

    /**
     * @brief Prepare concatenated scores for SQLite insertion
     *
     * Some scores might not be defined, those are reported as NULL
     *
     * @param[in] feature The feature being evaluated
     * @param[in] score_name The name of the queried score
     *
     * @returns A vector of strings with the queried scores
     *
     */
    std::vector<String> getSeparateScore(const Feature& feature, const std::string& score_name) const;

    /**
     * @brief Prepare feature rows for SQLite insertion
     *
     * This is the preferred high-throughput path. It avoids constructing full
     * SQL statements for every feature and lets writeRows bind values into
     * reusable prepared statements.
     *
     * @param[in] pep The compound (peptide/metabolite) used for extraction
     * @param[in] transition The transition used for extraction
     * @param[in] output The feature map containing all features
     * @param[in] id The transition group identifier (peptide/metabolite id)
     *
     * @return Buffered OSW table rows
     */
    OSWData prepareRows(const OpenSwath::LightCompound& pep,
        const OpenSwath::LightTransition* transition,
        const FeatureMap& output, const String& id) const;

    /**
     * @brief Append feature rows for SQLite insertion into an existing buffer
     *
     * This avoids constructing small per-assay OSWData buffers in high-throughput
     * schedulers.
     *
     * @param[in,out] rows The row buffer to append to
     * @param[in] pep The compound (peptide/metabolite) used for extraction
     * @param[in] transition The transition used for extraction
     * @param[in] output The feature map containing all features
     * @param[in] id The transition group identifier (peptide/metabolite id)
     */
    void prepareRowsInto(OSWData& rows,
        const OpenSwath::LightCompound& pep,
        const OpenSwath::LightTransition* transition,
        const FeatureMap& output, const String& id) const;

    /**
     * @brief Prepare a single line (feature) for output
     *
     * The result can be flushed to disk using writeLines (either line by line
     * or after collecting several lines).
     *
     * @param[in] pep The compound (peptide/metabolite) used for extraction
     * @param[in] transition The transition used for extraction
     * @param[in] output The feature map containing all features (each feature will generate one entry in the output)
     * @param[in] id The transition group identifier (peptide/metabolite id)
     *
     * @returns A string to be written using writeLines
     *
     */
    String prepareLine(const OpenSwath::LightCompound& pep,
        const OpenSwath::LightTransition* transition,
        const FeatureMap& output, const String& id) const;

    /**
     * @brief Write buffered OSW rows to disk
     *
     * Uses one SQLite transaction and reusable prepared statements for each OSW
     * table. This is the preferred writer path for newly generated output.
     *
     * @param[in] osw_output Rows generated by prepareRows
     *
     * @note Try to call this function as little as possible (it opens a new
     * database connection each time)
     *
     * @note Only call inside an OpenMP critical section
     *
     */
    void writeRows(const OSWData& osw_output);

    /**
     * @brief Write data to disk
     *
     * Takes a set of pre-prepared data statements from prepareLine and flushes them to disk
     *
     * @param[in] to_osw_output Statements generated by prepareLine
     *
     * @note Try to call this function as little as possible (it opens a new
     * database connection each time)
     *
     * @note Only call inside an OpenMP critical section
     *
     */
    void writeLines(const std::vector<String>& to_osw_output);

  };


}
