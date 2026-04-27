// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>
#include <memory>

namespace OpenMS
{
  class ExperimentalSettings;
  namespace Interfaces
  {
    class IMSDataConsumer;
  }
  

  /**
   * @brief File adapter for Swath files.
   *
   * This class can load SWATH files in different storage versions. The most
   * convenient file is a single MzML file which contains one experiment.
   * However, also the loading of a list of files is supported (loadSplit)
   * where it is assumed that each individual file only contains scans from one
   * precursor isolation window (one SWATH). Finally, experimental support for
   * mzXML is available but needs to be selected with a specific compile flag
   * (this is not for everyday use).
   *
   */
  class OPENMS_DLLAPI SwathFile :
    public ProgressLogger
  {
public:

    /// Loads a Swath run from a list of split mzML files
    std::vector<OpenSwath::SwathMap> loadSplit(StringList file_list,
                                               const String& tmp,
                                               std::shared_ptr<ExperimentalSettings>& exp_meta, 
                                               const String& readoptions = "normal");

    /**
      @brief Loads a Swath run from a single mzML file

      Using the @p plugin_consumer, you can provide a custom consumer which will be chained
      into the process of loading the data and making it available (depending on @p readoptions).
      This is useful if you want to modify the data a priori or extract some other information using
      MSDataTransformingConsumer (for example). Make sure it leaves the data intact, such that the 
      returned SwathMaps are actually useful.

      @param[in] file Input filename
      @param[in] tmp Temporary directory (for cached data)
      @param[out] exp_meta Experimental metadata from mzML file
      @param[in] readoptions How are spectra accessed after reading - tradeoff between memory usage and time (disk caching)
      @param[in] plugin_consumer An intermediate custom consumer
      @return Swath maps for MS2 and MS1 (unless readoptions == split, which returns no data)
    */
    std::vector<OpenSwath::SwathMap> loadMzML(const String& file, 
                                              const String& tmp,
                                              std::shared_ptr<ExperimentalSettings>& exp_meta,
                                              const String& readoptions = "normal",
                                              Interfaces::IMSDataConsumer* plugin_consumer = nullptr);

    /// Loads a Swath run from a single mzXML file
    std::vector<OpenSwath::SwathMap> loadMzXML(const String& file, 
                                               const String& tmp,
                                               std::shared_ptr<ExperimentalSettings>& exp_meta,
                                               const String& readoptions = "normal");

    /// Loads a Swath run from a single sqMass file
    std::vector<OpenSwath::SwathMap> loadSqMass(const String& file, std::shared_ptr<ExperimentalSettings>& /* exp_meta */);

#ifdef WITH_OPENTIMS
    /**
      @brief Loads a Swath run from a Bruker .d (TDF) directory

      @param[in] file Path to a Bruker .d (TDF) directory
      @param[in] tmp Temporary directory (for cached data)
      @param[in,out] exp_meta Will be filled with ExperimentalSettings metadata
      @param[in] readoptions How spectra are accessed: "normal" (in-memory) or "cache" (disk-cached)
      @return Vector of SwathMap structures representing the loaded Swath maps
    */
    std::vector<OpenSwath::SwathMap> loadBrukerTdf(const String& file,
                                                    const String& tmp,
                                                    std::shared_ptr<ExperimentalSettings>& exp_meta,
                                                    const String& readoptions);

    /// @brief Convenience overload: loads Bruker TDF in-memory (readoptions="normal")
    std::vector<OpenSwath::SwathMap> loadBrukerTdf(const String& file,
                                                    std::shared_ptr<ExperimentalSettings>& exp_meta);
#endif

protected:

    /// Cache a file to disk
    OpenSwath::SpectrumAccessPtr doCacheFile_(const String& in, const String& tmp, const String& tmp_fname,
                                              const std::shared_ptr<PeakMap >& experiment_metadata);

    /// Only read the meta data from a file and use it to populate exp_meta
    std::shared_ptr< PeakMap > populateMetaData_(const String& file);

    /// Counts the number of scans in a full Swath file (e.g. concatenated non-split file)
    void countScansInSwath_(const std::vector<MSSpectrum>& exp,
                            std::vector<int>& swath_counter, int& nr_ms1_spectra, 
                            std::vector<OpenSwath::SwathMap>& known_window_boundaries,
                            double TOLERANCE=1e-6);

  };
}
