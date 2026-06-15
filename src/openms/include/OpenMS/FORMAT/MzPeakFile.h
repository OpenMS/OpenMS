// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

/**
  @brief File adapter for the mzPeak format (Apache Parquet payloads packaged in a ZIP container).

  The mzPeak format stores mass spectrometric data as one or more Apache Parquet
  tables bundled together inside a ZIP archive. This class is intended to provide
  functions to read and write spectra and chromatograms in this format and to
  stream them to a data consumer.

  @note This class is currently a stub. The load(), store() and transform()
  methods are not yet implemented and will throw Exception::NotImplemented when
  called. The actual mzPeak reading/writing logic will be added in a later phase.
*/
class OPENMS_DLLAPI MzPeakFile
{
public:
  typedef MSExperiment MapType;

  /** @name Constructors and Destructor
   */
  //@{
  /// Default constructor
  MzPeakFile();

  /// Default destructor
  ~MzPeakFile();
  //@}

  /** @name Read / Write a complete mass spectrometric experiment
   */
  //@{

  /**
    @brief Load an mzPeak file into an MSExperiment.

    @note Not yet implemented; throws Exception::NotImplemented.
  */
  void load(const String& filename, MapType& map) const;

  /**
    @brief Store an MSExperiment in mzPeak format.

    Writes the experiment as point-layout Parquet tables (profile spectra to
    @c spectra_data.parquet, centroid spectra to @c spectra_peaks.parquet) plus
    a per-spectrum @c spectra_metadata.parquet and an @c mzpeak_index.json,
    bundled (STORED) into a ZIP archive. The output round-trips through load()
    to an equivalent experiment (spectrum count, ms_level, type, RT, peaks).

    @note Run-level metadata and precursor facets are not yet emitted.
  */
  void store(const String& filename, const MapType& map) const;

  /**
    @brief Stream the contents of an mzPeak file into a data consumer.

    @note Not yet implemented; throws Exception::NotImplemented.
  */
  void transform(const String& filename_in, Interfaces::IMSDataConsumer* consumer, bool skip_full_count = false, bool skip_first_pass = false) const;
  //@}
};
} // namespace OpenMS
