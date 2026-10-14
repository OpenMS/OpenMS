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
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/DATASTRUCTURES/TypeAliases.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>
#include <memory>
#include <limits>

namespace OpenMS
{
  class ExperimentalSettings;
  class IMSExperimentRandomAccess;

  /**
    @brief One FAIMS compensation-voltage group represented as ordinary OpenSWATH maps.

    @p faims_cv is the compensation voltage in volts. For non-FAIMS input it is NaN,
    matching IMDataConverter::splitByFAIMSCV(). The individual SwathMap objects remain
    unchanged and retain their normal DIA isolation-window semantics.
  */
  struct OPENMS_DLLAPI FAIMSSwathMapGroup
  {
    double faims_cv{std::numeric_limits<double>::quiet_NaN()};
    std::vector<OpenSwath::SwathMap> swath_maps;
  };

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
                                               const std::string& tmp,
                                               std::shared_ptr<ExperimentalSettings>& exp_meta, 
                                               const std::string& readoptions = "normal");

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
    std::vector<OpenSwath::SwathMap> loadMzML(const std::string& file, 
                                              const std::string& tmp,
                                              std::shared_ptr<ExperimentalSettings>& exp_meta,
                                              const std::string& readoptions = "normal",
                                              Interfaces::IMSDataConsumer* plugin_consumer = nullptr);

    /**
      @brief Loads a Swath run from a pre-loaded in-memory MSExperiment

      Used when the input format does not have a streaming reader (e.g. Thermo .raw
      via openms-thermo-bridge), so the full experiment has already been materialized.
      Avoids the round-trip through a temporary mzML file.

      @param[in] exp The pre-loaded experiment (must contain all spectra and metadata)
      @param[in] tmp Temporary directory (used only for readoptions=="cache")
      @param[out] exp_meta ExperimentalSettings extracted from @p exp
      @param[in] readoptions "normal" (in-memory) or "cache" (disk-cached)
      @return Swath maps for MS2 and MS1
    */
    std::vector<OpenSwath::SwathMap> loadFromMSExperiment(const std::shared_ptr<PeakMap>& exp,
                                                          const std::string& tmp,
                                                          std::shared_ptr<ExperimentalSettings>& exp_meta,
                                                          const std::string& readoptions = "normal");

    /**
      @brief Partition a pre-loaded experiment once by FAIMS CV and convert each group to SWATH maps.

      This is the FAIMS-aware adapter for in-memory inputs such as Thermo RAW data loaded through
      ThermoRawFile. It reuses IMDataConverter::splitByFAIMSCV() for CV assignment and
      loadFromMSExperiment() for the existing SWATH-window conversion. No FAIMS-specific state is
      added to OpenSwath::SwathMap; the CV is carried by the enclosing FAIMSSwathMapGroup.

      The input experiment is consumed. For non-FAIMS data, one group with a NaN CV is returned.
      For single-CV FAIMS data, one group with that CV is returned. Multi-CV input produces one
      group per exact compensation voltage.

      @param[in,out] exp Pre-loaded experiment; spectra are moved into per-CV groups
      @param[in] tmp Temporary directory (used only for readoptions=="cache")
      @param[out] exp_meta Experimental metadata copied before @p exp is consumed
      @param[in] readoptions "normal" (in-memory) or "cache" (disk-cached)
      @return FAIMS CV groups containing ordinary OpenSWATH maps
    */
    std::vector<FAIMSSwathMapGroup> loadFromMSExperimentByFAIMSCV(
      PeakMap&& exp,
      const std::string& tmp,
      std::shared_ptr<ExperimentalSettings>& exp_meta,
      const std::string& readoptions = "normal");

    /// Loads a Swath run from a single mzXML file
    std::vector<OpenSwath::SwathMap> loadMzXML(const std::string& file,
                                               const std::string& tmp,
                                               std::shared_ptr<ExperimentalSettings>& exp_meta,
                                               const std::string& readoptions = "normal");

    /// Loads a Swath run from a single sqMass file
    std::vector<OpenSwath::SwathMap> loadSqMass(const std::string& file, std::shared_ptr<ExperimentalSettings>& /* exp_meta */);

#ifdef WITH_OPENTIMS
    /**
      @brief Loads a Swath run from a Bruker .d (TDF) directory

      @param[in] file Path to a Bruker .d (TDF) directory
      @param[in] tmp Temporary directory (for cached data)
      @param[in,out] exp_meta Will be filled with ExperimentalSettings metadata
      @param[in] readoptions How spectra are accessed: "normal" (in-memory) or "cache" (disk-cached)
      @return Vector of SwathMap structures representing the loaded Swath maps
    */
    std::vector<OpenSwath::SwathMap> loadBrukerTdf(const std::string& file,
                                                    const std::string& tmp,
                                                    std::shared_ptr<ExperimentalSettings>& exp_meta,
                                                    const std::string& readoptions);

    /// @brief Convenience overload: loads Bruker TDF in-memory (readoptions="normal")
    std::vector<OpenSwath::SwathMap> loadBrukerTdf(const std::string& file,
                                                    std::shared_ptr<ExperimentalSettings>& exp_meta);
#endif

#ifdef WITH_THERMO_RAW
    /**
      @brief Loads a Swath run from a Thermo .raw file by streaming it scan by scan

      Like loadMzML(): a first pass reads the spectrum metadata only (no peaks) to determine the
      SWATH windows, then ThermoRawFile::transform() streams the scans into the SWATH consumer
      (and @p plugin_consumer), so that the run is never held in memory as a whole (except with
      readoptions "normal", which keeps the SWATH maps in memory). Chromatograms are not read.

      Readoptions "onDisc" forwards to loadOnDisc(const std::string&, ...).

      @param[in] file Path to a Thermo .raw file
      @param[in] tmp Temporary directory (for readoptions "cache" and "split")
      @param[out] exp_meta Run metadata (and the spectra without peaks, except for "onDisc")
      @param[in] readoptions "normal" (in-memory), "cache" (disk-cached), "split" (mzML per window) or "onDisc" (read on access)
      @param[in] plugin_consumer Optional consumer that also receives every spectrum
      @return Swath maps for MS2 and MS1 (empty for "split")
    */
    std::vector<OpenSwath::SwathMap> loadThermoRaw(const std::string& file,
                                                   const std::string& tmp,
                                                   std::shared_ptr<ExperimentalSettings>& exp_meta,
                                                   const std::string& readoptions = "normal",
                                                   Interfaces::IMSDataConsumer* plugin_consumer = nullptr);

#endif

    /**
      @brief SWATH maps that read the spectra on access from a random-access experiment ("onDisc")

      The SWATH windows are determined from the spectrum metadata as for the other read options; each
      map is a SpectrumAccessOnDisc view with its own handle (clone()) on @p experiment.

      @param[in] experiment Open random-access experiment
      @param[out] exp_meta Run metadata
      @param[in] plugin_consumer Optional consumer that observes every spectrum once. The maps read the
                 spectra again, so it must not modify them.
      @return Swath maps for MS2 and MS1

      @throws Exception::NullPointer if @p experiment is null
      @throws Exception::IllegalArgument if the input has more than one FAIMS compensation voltage
              (use loadOnDiscByFAIMSCV()) or if @p plugin_consumer modifies a spectrum
    */
    std::vector<OpenSwath::SwathMap> loadOnDisc(const std::shared_ptr<IMSExperimentRandomAccess>& experiment,
                                                std::shared_ptr<ExperimentalSettings>& exp_meta,
                                                Interfaces::IMSDataConsumer* plugin_consumer = nullptr);

    /**
      @brief Loads a Swath run with readoptions "onDisc": the file is indexed once and the spectra are read on access

      The random-access backend is chosen by file type: OnDiscMSExperiment for indexed mzML,
      OnDiscThermoRawExperiment for Thermo .raw (requires WITH_THERMO_RAW; centroided, as with the
      other read options). Nothing is loaded or cached up front; the spectra and SWATH windows are the
      same as with the other read options. See loadOnDisc(const std::shared_ptr<IMSExperimentRandomAccess>&, ...).

      @param[in] file Indexed mzML or Thermo .raw file
      @param[out] exp_meta Run metadata
      @param[in] plugin_consumer Optional consumer that observes every spectrum once (it must not modify them)
      @return Swath maps for MS2 and MS1

      @throws Exception::InvalidValue if the file type has no random-access backend
      @throws Exception::ParseError if an mzML file is not indexed
      @throws Exception::IllegalArgument if the input has more than one FAIMS compensation voltage
              or if @p plugin_consumer modifies a spectrum
    */
    std::vector<OpenSwath::SwathMap> loadOnDisc(const std::string& file,
                                                std::shared_ptr<ExperimentalSettings>& exp_meta,
                                                Interfaces::IMSDataConsumer* plugin_consumer = nullptr);

    /**
      @brief Like loadOnDisc(), but partitions MS1 and MS2 spectra by FAIMS compensation voltage first

      The assignment follows IMDataConverter::splitByFAIMSCV(): a spectrum with drift time unit
      FAIMS_COMPENSATION_VOLTAGE goes to the CV in its drift time, an MS2+ spectrum without a CV to the
      most recent CV before it, other spectra are skipped. SWATH windows are then determined within
      each CV. Without FAIMS data, one group with CV NaN holds all spectra.

      @param[in] experiment Open random-access experiment
      @param[out] exp_meta Run metadata
      @param[in] plugin_consumer Optional consumer that observes every spectrum once (it must not modify them)
      @return SWATH maps per FAIMS compensation voltage, in ascending CV order

      @throws Exception::NullPointer if @p experiment is null
      @throws Exception::IllegalArgument if @p plugin_consumer modifies a spectrum
    */
    std::vector<FAIMSSwathMapGroup> loadOnDiscByFAIMSCV(const std::shared_ptr<IMSExperimentRandomAccess>& experiment,
                                                        std::shared_ptr<ExperimentalSettings>& exp_meta,
                                                        Interfaces::IMSDataConsumer* plugin_consumer = nullptr);

    /// loadOnDiscByFAIMSCV() on a file; the backend is chosen as in loadOnDisc(const std::string&, ...)
    std::vector<FAIMSSwathMapGroup> loadOnDiscByFAIMSCV(const std::string& file,
                                                        std::shared_ptr<ExperimentalSettings>& exp_meta,
                                                        Interfaces::IMSDataConsumer* plugin_consumer = nullptr);

protected:

    /// SpectrumAccessOnDisc maps for the spectra @p indices of @p experiment (@p metadata: all spectra without peaks)
    std::vector<OpenSwath::SwathMap> buildOnDiscMaps_(const std::shared_ptr<IMSExperimentRandomAccess>& experiment,
                                                      const PeakMap& metadata,
                                                      const std::vector<Size>& indices);

    /// Passes every spectrum to @p plugin_consumer; throws Exception::IllegalArgument if it modifies one
    void passOnDiscToPlugin_(IMSExperimentRandomAccess& experiment,
                             const PeakMap& metadata,
                             ExperimentalSettings& exp_meta,
                             Interfaces::IMSDataConsumer* plugin_consumer);

    /// Opens @p file with the random-access backend for its type (see loadOnDisc(const std::string&, ...))
    std::shared_ptr<IMSExperimentRandomAccess> openOnDisc_(const std::string& file);

    /// Cache a file to disk
    OpenSwath::SpectrumAccessPtr doCacheFile_(const std::string& in, const std::string& tmp, const std::string& tmp_fname,
                                              const std::shared_ptr<PeakMap >& experiment_metadata);

    /// Only read the meta data from a file and use it to populate exp_meta
    std::shared_ptr< PeakMap > populateMetaData_(const std::string& file);

    /// Counts the number of scans in a full Swath file (e.g. concatenated non-split file)
    void countScansInSwath_(const std::vector<MSSpectrum>& exp,
                            std::vector<int>& swath_counter, int& nr_ms1_spectra, 
                            std::vector<OpenSwath::SwathMap>& known_window_boundaries,
                            double TOLERANCE=1e-6);

  };
}
