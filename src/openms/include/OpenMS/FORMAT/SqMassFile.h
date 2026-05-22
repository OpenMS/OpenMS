// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

namespace OpenMS
{

  /**
    @brief Read and write mass-spectrometry data in the on-disk @c sqMass (SQLite) format.

    Persists spectra and chromatograms in a single SQLite database so that
    large collections can be queried, filtered, or streamed on demand
    without loading the whole experiment into memory.

    For spectra and chromatograms that carry precursor information, the
    @c "peptide_sequence" meta-value of the first precursor (if any) is
    round-tripped alongside the binary data.

    Compression and meta-data options are configured via
    @ref SqMassConfig (see @ref setConfig).

    @ingroup FileIO
  */
  class OPENMS_DLLAPI SqMassFile
  {
public:

    /**
      @brief Compression / metadata options used when reading or writing an @c sqMass file.
    */
    struct OPENMS_DLLAPI SqMassConfig
    {
      bool write_full_meta{true};      ///< Write full meta data (verbose, lossless).
      bool use_lossy_numpress{false};  ///< Use lossy numpress compression for the binary data arrays.
      double linear_fp_mass_acc{-1};   ///< Target mass accuracy for numpress linear encoding (@c -1 has no effect; use e.g. @c 0.0001 for 0.2 ppm @c @ 500 m/z).
    };

    /// Convenience alias used by @ref load / @ref store.
    typedef MSExperiment MapType;

    /** @name Constructors and Destructor
    */
    //@{
    /// Default constructor.
    SqMassFile();

    /// Destructor.
    ~SqMassFile();
    //@}

    /** @name Read / Write a complete mass spectrometric experiment
    */
    //@{

    /**
      @brief Read the contents of an @c sqMass file into an @c MSExperiment.

      Loads all spectra, chromatograms and experimental metadata from
      @p filename into @p map. The current @ref SqMassConfig is applied
      to the read.

      @param[in]  filename Path to the @c sqMass file to read.
      @param[out] map      Destination experiment; populated with the
                           file's spectra, chromatograms and metadata.

      @throws Exception::IllegalArgument    On malformed or inconsistent
                                            rows in the sqMass tables
                                            (mismatched native id,
                                            missing data type, ...) or
                                            when the file contains more
                                            than one run.
      @throws Exception::ConversionError    When a binary buffer's size
                                            does not match its declared
                                            data type.
      @throws Exception::SqlOperationFailed When the file's @c RUN table
                                            advertises an unsupported
                                            configuration.
    */
    void load(const String& filename, MapType& map) const;

    /**
      @brief Store an @c MSExperiment in @c sqMass format.

      Writes the spectra, chromatograms and experimental metadata of
      @p map to @p filename, creating the file (and the required SQLite
      tables) if necessary. The current @ref SqMassConfig is applied.

      The sqMass @c RUN::ID column is taken from
      @c MSExperiment::getSqlRunID; populate it via
      @c MSExperiment::setSqlRunID before calling @ref store if a
      specific value is required.

      @param[in] filename Path to the output @c sqMass file.
      @param[in] map      Experiment to serialise.

      @throws Exception::IllegalArgument    When SQL commands fail during
                                            table creation or data
                                            insertion.
      @throws Exception::SqlOperationFailed When the database file cannot
                                            be created or opened for
                                            writing.
    */
    void store(const String& filename, const MapType& map) const;

    /**
      @brief Stream the spectra and chromatograms of an @c sqMass file through @p consumer.

      Reads @p filename and feeds every spectrum and chromatogram to
      @p consumer in input order. The consumer is informed of the
      expected counts (via @c IMSDataConsumer::setExpectedSize) and of
      the experimental metadata (via @c setExperimentalSettings) before
      the per-element callbacks are invoked. The full experiment is not
      held in memory.

      @note The @p skip_full_count and @p skip_first_pass parameters are
            currently unused and have no effect; they are retained for
            API compatibility.

      @param[in] filename_in     Path to the input @c sqMass file.
      @param[in] consumer        Receives the streamed spectra and
                                 chromatograms. Must remain valid for
                                 the duration of the call.
      @param[in] skip_full_count Currently ignored.
      @param[in] skip_first_pass Currently ignored.

      @throws Exception::IllegalArgument    On malformed or inconsistent
                                            rows in the sqMass tables
                                            (mismatched native id,
                                            missing data type, ...) or
                                            when the file contains more
                                            than one run.
      @throws Exception::ConversionError    When a binary buffer's size
                                            does not match its declared
                                            data type.
      @throws Exception::SqlOperationFailed When the file's @c RUN table
                                            advertises an unsupported
                                            configuration.
    */
    void transform(const String& filename_in, Interfaces::IMSDataConsumer* consumer, bool skip_full_count = false, bool skip_first_pass = false) const;

    /**
      @brief Convert an sqMass file containing chromatogram data to an XIC Parquet file.

      This is a convenience function that will stream chromatograms from the
      input sqMass file into an OpenMS Parquet writer (MSChromatogramParquetConsumer)
      and write them to disk.

      Note: The transition/precursor metadata passed to the underlying Parquet
      consumer is optional from the API perspective, but the writer expects
      chromatogram native IDs to have matching transition/precursor metadata in
      many cases. If required metadata for a chromatogram (precursor or
      transition) is missing the writer will throw an exception and the
      conversion will fail. Therefore it is recommended to provide a
      populated transition experiment when possible.

      @param[in] filename_in Path to the input sqMass file.
      @param[in] xic_filename Path to the output .xic Parquet file.
      @param[in] run_id Run identifier to store with each chromatogram (default 0).
      @param[in] source_file Source filename to store in the parquet file (if empty, filename_in is used).
      @param[in] transition_exp Optional transition/precursor metadata (LightTargetedExperiment) used to annotate chromatograms. If left empty, no transition metadata will be available to the writer and conversion may fail for chromatograms that require such metadata. Prefer supplying a populated transition experiment when available.

      @throws Exception::InvalidValue If a chromatogram refers to a precursor or
      transition for which no matching metadata entry exists, or if Arrow/Parquet
      operations fail while assembling/writing the table.
      @throws Exception::FileNotWritable If the output parquet file cannot be
      opened for writing.
      @throws Exception::NotImplemented If Parquet support was not compiled in.
      @throws Exception::BaseException Other OpenMS exceptions propagated from
      the writer/encoding layers may also be thrown.
    */
    void convertToXICParquet(const String& filename_in,
                const String& xic_filename,
                UInt64 run_id = 0,
                const String& source_file = "",
                const OpenSwath::LightTargetedExperiment& transition_exp = OpenSwath::LightTargetedExperiment()) const;

    /**
      @brief Replace the compression / metadata configuration used by subsequent @ref load, @ref store and @ref transform calls.

      @param[in] config New configuration.
    */
    void setConfig(const SqMassConfig& config)
    {
      config_ = config;
    }

    // maybe later ...
    // static inline void readSpectrumFast(OpenSwath::BinaryDataArrayPtr data1,
    //                                     OpenSwath::BinaryDataArrayPtr data2, std::ifstream& ifs, int& ms_level,
    //                                     double& rt)

    // static inline void readChromatogramFast(OpenSwath::BinaryDataArrayPtr data1,
    //                                        OpenSwath::BinaryDataArrayPtr data2, std::ifstream& ifs)

protected:
      SqMassConfig config_;
  };
}


