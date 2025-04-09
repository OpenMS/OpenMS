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


