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
    @brief An class that uses on-disk SQLite database to read and write spectra and chromatograms

    This class provides functions to read and write spectra and chromatograms
    to disk using a SQLite database and store them in sqMass format. This
    allows users to access, select and filter spectra and chromatograms
    on-demand even in a large collection of data.

    Spectra and chromatograms with precursor information will additionally load/store the metavalue
    'peptide_sequence' from the first precursor (if any).

  */
  class OPENMS_DLLAPI SqMassFile
  {
public:

  /**
    @brief Configuration class for SqMassFile

    Contains configuration options for SQLite file
  */
    struct OPENMS_DLLAPI SqMassConfig 
    {
      bool write_full_meta{true}; ///< write full meta data
      bool use_lossy_numpress{false}; ///< use lossy numpress compression
      double linear_fp_mass_acc{-1}; ///< desired mass accuracy for numpress linear encoding (-1 no effect, use 0.0001 for 0.2 ppm accuracy @ 500 m/z)
    };

    typedef MSExperiment MapType;

    /** @name Constructors and Destructor
    */
    //@{
    /// Default constructor
    SqMassFile();

    /// Default destructor
    ~SqMassFile();
    //@}

    /** @name Read / Write a complete mass spectrometric experiment
    */
    //@{

    void load(const String& filename, MapType& map) const;

    /**
     @brief Store an MSExperiment in sqMass format

     If you want a specific RUN::ID in the sqMass file,
     make sure to populate MSExperiment::setSqlRunID(UInt64 id) before.
    */
    void store(const String& filename, const MapType& map) const;

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


