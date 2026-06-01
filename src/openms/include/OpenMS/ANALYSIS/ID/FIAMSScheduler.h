// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Svetlana Kutuzova, Douglas McCloskey $
// $Authors: Svetlana Kutuzova, Douglas McCloskey $
// --------------------------------------------------------------------------
 
#pragma once

#include <OpenMS/FORMAT/CsvFile.h>
#include <map>

namespace OpenMS
{
class MzMLFile;
/**
    @brief Batch driver for FIA-MS analyses. Dispatches each sample of a CSV-described batch to @ref FIAMSDataProcessor.

    A CSV file describing the batch is parsed during construction. Each row defines one sample
    and is consumed by @ref run, which loads the corresponding mzML file, configures a
    @ref FIAMSDataProcessor instance, and runs the analysis once per requested time window.
    Samples are processed in parallel (OpenMP-parallelised loop).

    The CSV is comma-delimited; the first row is taken as the header. Each subsequent row is
    stored as a header-keyed @c std::map<String, String>. The implementation looks up the
    following keys per row (column names as read from the CSV):

    - @c filename — base filename of the sample's mzML (no path, no @c .mzML extension);
      the loader appends @c .mzML.
    - @c dir_input — sample subdirectory relative to @p base_dir.
    - @c dir_output — sample subdirectory relative to @p output_dir; created on demand.
    - @c resolution — instrument resolution; parsed with @c std::stof.
    - @c charge — instrument polarity; forwarded as the @c polarity parameter of
      @ref FIAMSDataProcessor (typical values @c "positive" / @c "negative").
    - @c db_mapping — accurate-mass database mapping file (relative to @p base_dir); forwarded
      as the @c db:mapping parameter (three tab-separated columns: mass, formula, identifier).
    - @c db_struct — accurate-mass database structure file (relative to @p base_dir);
      forwarded as the @c db:struct parameter (four tab-separated columns: identifier, name,
      SMILES, InChI).
    - @c positive_adducts / @c negative_adducts — adduct list files (relative to @p base_dir).
    - @c time — @c ";"-separated list of time windows in seconds (e.g. @c "30;60;90;180");
      @ref FIAMSDataProcessor::run is invoked once per window.

    The mzML path opened for each sample is therefore @c base_dir + dir_input + "/" + filename + ".mzML".

    No column validation is performed; missing or mis-named keys raise the usual
    @c std::map::at exception.

    @ingroup Analysis_ID
*/
class OPENMS_DLLAPI FIAMSScheduler
  {
public:
    /// Default-construct an empty scheduler (no samples loaded, paths set to defaults).
    FIAMSScheduler() = default;

    /**
      @brief Construct from a batch CSV and parse it immediately.

      The CSV at @p filename is read during construction (see class documentation for the
      expected columns) and the parsed rows are kept for subsequent @ref run calls.

      @param filename     Full path to the batch CSV file.
      @param base_dir     Base directory used to resolve the per-row @c dir_input, @c db_mapping,
                          @c db_struct, and adduct paths. Must end with a path separator.
      @param output_dir   Output directory under which each row's @c dir_output is created.
                          Must end with a path separator.
      @param load_cached_ Forwarded to @ref FIAMSDataProcessor::run as its @c load_cached
                          flag; if @c true, cached intermediate results are reused.
    */
    FIAMSScheduler(
      String filename,
      String base_dir = "/",
      String output_dir = "/",
      bool load_cached_ = true
    );

    /// Default destructor.
    ~FIAMSScheduler() = default;

    /// Copy constructor.
    FIAMSScheduler(const FIAMSScheduler& cp) = default;

    /// Copy assignment.
    FIAMSScheduler& operator=(const FIAMSScheduler& fdp) = default;

    /**
      @brief Run the batch: process every sample loaded from the CSV.

      For each sample (rows are processed by an @c \#pragma @c omp @c parallel @c for loop),
      this loads the mzML at @c base_dir/dir_input/filename.mzML, creates a fresh
      @ref FIAMSDataProcessor with the per-row parameters, creates the per-sample
      @c dir_output directory if missing, and runs the processor once per time window
      listed in the row's @c time cell. Progress messages are emitted via @c OPENMS_LOG_INFO
      at the start and end of every (sample, time) pair.
    */
    void run();

    /// Return the parsed batch as a vector of header-keyed @c std::map rows.
    const std::vector<std::map<String, String>>& getSamples();

    /// Return the base directory passed to the constructor.
    const String& getBaseDir();

    /// Return the output directory passed to the constructor.
    const String& getOutputDir();

private:
    /**
      @brief Parse the batch CSV at @c filename_ and populate @c samples_.

      Treats the first row as headers and every subsequent row as a sample; each row is
      stored as a header-keyed @c std::map<String, String>. Invoked by the constructor.
    */
    void loadSamples_();

    String filename_;                                       ///< CSV file describing the batch.
    String base_dir_;                                       ///< Base directory for per-sample input paths (must end with a path separator).
    String output_dir_;                                     ///< Base directory under which per-sample outputs are written (must end with a path separator).
    bool load_cached_;                                      ///< Forwarded to @ref FIAMSDataProcessor::run; reuse cached results when @c true.
    std::vector<std::map<String, String>> samples_;         ///< Parsed CSV rows; populated by @ref loadSamples_.
  };
} // namespace OpenMS
