// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest, Justin Sing$
// $Authors: Hannes Roest, Justin Sing$
// --------------------------------------------------------------------------

#pragma once

#include <memory>

// Consumers
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h>

// Files
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>

// Kernel and implementations
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// Helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

// Algorithms
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>

#include <cassert>
#include <limits>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

namespace OpenMS
{
  /**
    @brief Shared @ref TOPPBase scaffolding for OpenSWATH-family CLI tools (currently @ref TOPP_OpenSwathWorkflow).

    Bundles the input-loading and output-setup steps that any OpenSWATH DIA tool needs, so the
    individual TOPP front-ends can stay focused on configuration and workflow wiring. The
    protected methods cover the three I/O stages of a typical run:

      - loadSwathFiles() — read one or more mzML / mzXML / sqMass SWATH input files via
        @ref SwathFile (with the @c readoptions parameter selecting in-memory vs. cached
        access), optionally remap the per-map isolation windows from a @c swath_windows_file,
        and return both the per-map @ref OpenSwath::SwathMap pointers and a one-to-one list of
        source-file names.
      - prepareChromOutput() / prepareMobilogramOutput() — install the output @ref Interfaces::IMSDataConsumer
        for chromatograms (sqMass / mzML+numpress / Parquet @c .xic) and the optional Parquet
        mobilogram writer.
      - loadTransitionList() — read the spectral library from TraML / TSV / PQP into an
        @ref OpenSwath::LightTargetedExperiment.

    The @ref CalibrationResult inner struct is the return-channel for the per-run RT / m/z /
    IM calibration step run by the consuming tool (see @ref CalibrationWorkflow).

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI TOPPOpenSwathBase : public TOPPBase
  {

  public:
    /**
      @brief Per-run outputs of the RT / m/z / IM calibration step.

      Populated by the consuming tool (via @ref CalibrationWorkflow) and consumed later
      during the analytical extraction phase. Window fields default to @c -1.0 meaning
      "not computed" — callers must check before using them as extraction windows.
    */
    struct CalibrationResult
    {
      /// RT normalization transformation (fitted Trafo)
      OpenMS::TransformationDescription rt_trafo;

      /// MS2 m/z extraction window (full width, ppm). -1 if not computed.
      double ms2_mz_window_ppm{ -1.0 };

      /// MS2 ion mobility extraction window (full width, native IM units). -1 if not computed/applicable.
      double ms2_im_window{ -1.0 };

      /// MS1 m/z extraction window (full width, ppm). -1 if not computed.
      double ms1_mz_window_ppm{ -1.0 };

      /// MS1 ion mobility extraction window (full width, native IM units). -1 if not computed/applicable.
      double ms1_im_window{ -1.0 };
    };

    /**
      @brief Constructor

      Must match TOPPBase' Ctor!

      @param[in] name Tool name.
      @param[in] description Short description of the tool (one line).
      @param[in] official If this is an official TOPP tool contained in the OpenMS/TOPP release.
             If @em true the tool name is checked against the list of TOPP tools and a warning printed if missing.
      @param[in] citations Add one or more citations if they are associated specifically to this TOPP tool; they will be printed during `--help`
    */
    TOPPOpenSwathBase(String name, String description, bool official = true, const std::vector<Citation>& citations = {});

    /// Destructor
    ~TOPPOpenSwathBase() override;

  protected:
    /**
     * @brief Load the DIA files into internal data structures.
     *
     * Loads SWATH files into the provided OpenSwath::SwathMap data structures. It
     * uses the SwathFile class to load files from either mzML, mzXML or SqMass.
     * The files will be either loaded into memory or cached to disk (depending on
     * the readoptions parameter).
     *
     * @param[in] file_list The input file(s)
     * @param[out] exp_meta The output (meta data about experiment)
     * @param[out] swath_maps The output (ptr to raw data)
     * @param[out] swath_map_sources The source files corresponding to each swath map. This is used to track when multiple experiment input files are provided.
     * @param[in] split_file If loading a single file that contains a single SWATH window
     * @param[in] tmp Temporary directory
     * @param[in] readoptions Description on how to read the data ("normal", "cache")
     * @param[in] swath_windows_file Provided file containing the SWATH windows which will be mapped to the experimental windows
     * @param[in] min_upper_edge_dist Distance for each assay to the upper edge of the SWATH window
     * @param[in] force Whether to override the sanity check
     * @param[in] sort_swath_maps Whether to sort the provided windows first before mapping
     * @param[in] prm Whether data is in prm format; allows for overlap
     * @param[in,out] plugin_consumer Intermediate consumer for mzML input. See SwathFile::loadMzML() for details.
     *
     * @return Returns whether loading and sanity check was successful
     *
     */
    bool loadSwathFiles(const StringList& file_list,
                        std::shared_ptr<ExperimentalSettings >& exp_meta,
                        std::vector< OpenSwath::SwathMap >& swath_maps,
                        std::vector<String> & swath_map_sources,
                        const bool split_file,
                        const String& tmp,
                        const String& readoptions,
                        const String& swath_windows_file,
                        const double min_upper_edge_dist,
                        const bool force,
                        const bool sort_swath_maps,
                        const bool prm,
                        Interfaces::IMSDataConsumer* plugin_consumer = nullptr);

    /**
     * @brief Prepare chromatogram output
     *
     * Sets up the chromatogram output, either sqMass, mzML (using numpress
     * lossy compression), or xic (Parquet). This assumes that 0.05 accuracy
     * in RT is sufficient for all purposes.
     *
     * @param[out] chromatogramConsumer The consumer to process chromatograms
     * @param[in] exp_meta meta data about experiment
     * @param[in] transition_exp The spectral library
     * @param[in] out_chrom The output file for the chromatograms
     * @param[in] run_id Unique identifier which links the sqMass and OSW file
     * @param[in] source_file Source file name for chromatogram provenance
     */
    void prepareChromOutput(Interfaces::IMSDataConsumer ** chromatogramConsumer,
                            const std::shared_ptr<ExperimentalSettings>& exp_meta,
                            const OpenSwath::LightTargetedExperiment& transition_exp,
                            const String& out_chrom,
                            const UInt64 run_id,
                            const String& source_file);

    /**
     * @brief Prepare mobilogram output
     *
     * Sets up the mobilogram output, currently only supports Parquet via MobilogramParquetConsumer.
     * If no output requested, returns a null consumer by allocating a nullptr.
     *
     * @param[out] mobilogramConsumer The consumer to process mobilograms
     * @param[in] exp_meta meta data about experiment
     * @param[in] transition_exp The spectral library
     * @param[in] out_mobilogram The output file for the mobilograms
     * @param[in] run_id Unique identifier which links the mobilogram and OSW file
     * @param[in] source_file Source file name for provenance
     */
    void prepareMobilogramOutput(std::unique_ptr<class MobilogramParquetConsumer>& mobilogramConsumer,
                  const std::shared_ptr<ExperimentalSettings>& exp_meta,
                  const OpenSwath::LightTargetedExperiment& transition_exp,
                  const String& out_mobilogram,
                  const UInt64 run_id,
                  const String& source_file);

    /**
     * @brief Loads transition list from TraML / TSV or PQP
     *
     * @param[in] tr_type Input file type
     * @param[in] tr_file Input file name
     * @param[in] tsv_reader_param Parameters on how to interpret spectral data
     *
     */
    OpenSwath::LightTargetedExperiment loadTransitionList(const FileTypes::Type& tr_type,
                                                          const String& tr_file,
                                                          const Param& tsv_reader_param);

  private:
    void loadSwathFiles_(const StringList& file_list,
                         const bool split_file,
                         const String& tmp,
                         const String& readoptions,
                         std::shared_ptr<ExperimentalSettings > & exp_meta,
                         std::vector< OpenSwath::SwathMap > & swath_maps,
                         std::vector<String> & swath_map_sources,
                         Interfaces::IMSDataConsumer* plugin_consumer);
  }; // end TOPPOpenSwathBase
} //  end NS OpenMS
