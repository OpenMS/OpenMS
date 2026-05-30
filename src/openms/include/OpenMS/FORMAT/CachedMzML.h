// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>

#include <fstream>

namespace OpenMS
{

  /**
    @brief Random-access reader/writer for the OpenMS on-disk cache format.

    Spectra and chromatograms are split across two files: an @c .mzML file
    holding the experiment-level metadata and a side-car @c .mzML.cached
    binary file holding the per-spectrum / per-chromatogram payload. After
    @ref load (or constructing with a filename) the side-car's stream
    offsets are kept in an in-memory index, so single spectra or
    chromatograms can be fetched in any order from disk without scanning
    the full file.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI CachedmzML
  {

public:

    /** @name Constructors and Destructor
    */
    //@{
    /**
      @brief Default constructor.

      Produces an empty instance with no associated file.
    */
    CachedmzML();

    /**
      @brief Construct and load metadata + side-car index from @p filename in one step.

      Equivalent to default-constructing and then calling @ref load. See
      @ref load for the exact behaviour, file layout and exceptions.

      @param[in] filename Path of the @c .mzML metadata file; the side-car
                          is expected next to it as @p filename + @c ".cached".
    */
    CachedmzML(const String& filename);

    /**
      @brief Copy constructor.

      Each instance keeps its own file handle on the @c .mzML.cached
      side-car, so the copy opens a fresh stream on the same path -- the
      original's stream is not shared.

      @param[in] rhs Instance to copy.
    */
    CachedmzML(const CachedmzML & rhs);

    /**
      @brief Destructor.

      Closes the open file stream on the @c .mzML.cached side-car.
    */
    ~CachedmzML();
    //@}

    /**
      @brief Read one spectrum from the side-car file.

      Seeks the cached stream to the indexed offset for @p id, parses the
      binary payload and returns the resulting @c MSSpectrum. The metadata
      for the spectrum (instrument settings, native id, ...) is taken
      from the in-memory experiment loaded at @ref load time.

      @param[in] id Index of the spectrum to read; must be smaller than
                    @ref getNrSpectra. Out-of-range access is a programming
                    error (checked in debug builds).
      @return Decoded spectrum.

      @throws Exception::ParseError when @c seekg on the side-car stream
                                    fails (e.g. a 32-bit-system overflow
                                    on a > 2 GB cache file).
    */
    MSSpectrum getSpectrum(Size id);

    /**
      @brief Read one chromatogram from the side-car file.

      @param[in] id Index of the chromatogram to read; must be smaller than
                    @ref getNrChromatograms. Out-of-range access is a
                    programming error (checked in debug builds).
      @return Decoded chromatogram.

      @throws Exception::ParseError when @c seekg on the side-car stream
                                    fails.
    */
    MSChromatogram getChromatogram(Size id);

    /**
      @brief Number of spectra in the loaded experiment.

      @return Spectrum count.
    */
    size_t getNrSpectra() const;

    /**
      @brief Number of chromatograms in the loaded experiment.

      @return Chromatogram count.
    */
    size_t getNrChromatograms() const;

    /**
      @brief Experiment-level metadata loaded from the @c .mzML file.

      The returned experiment contains instrument and acquisition settings
      as well as the per-spectrum / per-chromatogram headers. Empty until
      @ref load has been called.

      @return Const reference to the loaded metadata experiment.
    */
    const MSExperiment& getMetaData() const
    {
      return meta_ms_experiment_;
    }

    /**
      @brief Write @p map to disk as a cached @c .mzML pair.

      Produces two files: a binary side-car at @p filename + @c ".cached"
      holding the spectrum / chromatogram payload, and an @c .mzML file at
      @p filename holding the metadata.

      @param[in] filename Path of the metadata file (should end in @c .mzML); the side-car is written next to it as @p filename + @c ".cached".
      @param[in] map      Experiment to write.

      @throws Exception::UnableToCreateFile when either output file cannot be created.
    */
    static void store(const String& filename, const PeakMap& map);

    /**
      @brief Load a cached @c .mzML pair into @p map.

      Reads the metadata from @p filename and builds the in-memory index
      for the side-car at @p filename + @c ".cached"; subsequent calls to
      @ref getSpectrum / @ref getChromatogram on @p map fetch the payload
      from the side-car on demand.

      @param[in]  filename Path of the @c .mzML metadata file (the side-car is expected next to it as @p filename + @c ".cached").
      @param[out] map      Destination instance; previous contents are replaced.

      @throws Exception::FileNotFound when either file cannot be opened.
      @throws Exception::ParseError   when parsing the metadata or building the side-car index fails.
    */
    static void load(const String& filename, CachedmzML& map);

protected:

    void load_(const String& filename);

    /// Meta data
    MSExperiment meta_ms_experiment_;

    /// Internal filestream 
    std::ifstream ifs_;

    /// Name of the mzML file
    String filename_;

    /// Name of the cached mzML file
    String filename_cached_;

    /// Indices
    std::vector<std::streampos> spectra_index_;
    std::vector<std::streampos> chrom_index_;

  };
}

