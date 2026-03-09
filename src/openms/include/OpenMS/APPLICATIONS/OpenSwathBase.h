// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest, Justin Sing$
// $Authors: Hannes Roest, Justin Sing$
// --------------------------------------------------------------------------

#pragma once

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
  class OPENMS_DLLAPI TOPPOpenSwathBase : public TOPPBase
  {

  public:
    /// Outputs of RT, m/z, IM calibration
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

      @param name Tool name.
      @param description Short description of the tool (one line).
      @param official If this is an official TOPP tool contained in the OpenMS/TOPP release.
             If @em true the tool name is checked against the list of TOPP tools and a warning printed if missing.
      @param citations Add one or more citations if they are associated specifically to this TOPP tool; they will be printed during `--help`
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
     * @param file_list The input file(s)
     * @param exp_meta The output (meta data about experiment)
     * @param swath_maps The output (ptr to raw data)
     * @param split_file If loading a single file that contains a single SWATH window
     * @param tmp Temporary directory
     * @param readoptions Description on how to read the data ("normal", "cache")
     * @param swath_windows_file Provided file containing the SWATH windows which will be mapped to the experimental windows
     * @param min_upper_edge_dist Distance for each assay to the upper edge of the SWATH window
     * @param force Whether to override the sanity check
     * @param sort_swath_maps Whether to sort the provided windows first before mapping
     * @param prm Whether data is in prm format; allows for overlap
     * @param pasef Whether data is in PASEF format; allows for overlap
     * @param plugin_consumer Intermediate consumer for mzML input. See SwathFile::loadMzML() for details.
     *
     * @return Returns whether loading and sanity check was successful
     *
     */
    bool loadSwathFiles(const StringList& file_list,
                        std::shared_ptr<ExperimentalSettings >& exp_meta,
                        std::vector< OpenSwath::SwathMap >& swath_maps,
                        const bool split_file,
                        const String& tmp,
                        const String& readoptions,
                        const String& swath_windows_file,
                        const double min_upper_edge_dist,
                        const bool force,
                        const bool sort_swath_maps,
                        const bool prm,
                        const bool pasef,
                        Interfaces::IMSDataConsumer* plugin_consumer = nullptr);

    /**
     * @brief Prepare chromatogram output
     *
     * Sets up the chromatogram output, either sqMass or mzML (using numpress
     * lossy compression). This assumes that 0.05 accuracy in RT is sufficient
     * for all purposes.
     *
     * @param chromatogramConsumer The consumer to process chromatograms
     * @param exp_meta meta data about experiment
     * @param transition_exp The spectral library
     * @param out_chrom The output file for the chromatograms
     * @param run_id Unique identifier which links the sqMass and OSW file
     */
    void prepareChromOutput(Interfaces::IMSDataConsumer ** chromatogramConsumer,
                            const std::shared_ptr<ExperimentalSettings>& exp_meta,
                            const OpenSwath::LightTargetedExperiment& transition_exp,
                            const String& out_chrom,
                            const UInt64 run_id);

    /**
     * @brief Loads transition list from TraML / TSV or PQP
     *
     * @param tr_type Input file type
     * @param tr_file Input file name
     * @param tsv_reader_param Parameters on how to interpret spectral data
     *
     */
    OpenSwath::LightTargetedExperiment loadTransitionList(const FileTypes::Type& tr_type,
                                                          const String& tr_file,
                                                          const Param& tsv_reader_param);

    /**
     * @brief Perform retention time and m/z calibration
     *
     * This function will create the retention time transformation either by
     * loading a provided .trafoXML file or determine it from the data itself by
     * extracting the transitions specified in the irt_tr_file TraML file. It
     * will also perform the m/z calibration (when an irt_tr_file is provided).
     *
     * @note Internally, the retention time and @p m/z calibration are performed
     * by OpenMS::OpenSwathCalibrationWorkflow::performRTNormalization
     *
     * @param trafo_in Input trafoXML file (if not empty, transformation will be
     *                 loaded from this file)
     * @param irt_transitions  Input iRT transition experiment (if trafo_in
     *                     is empty, this will be used for iRT extraction)
     * @param swath_maps The raw data (swath maps)
     * @param min_rsq Minimal R^2 value that is expected for the RT regression
     * @param min_coverage Minimal coverage of the chromatographic space that needs to be achieved
     * @param feature_finder_param Parameter set for the feature finding in chromatographic dimension
     * @param cp_irt Parameter set for the chromatogram extraction
     * @param irt_detection_param Parameter set for the detection of the iRTs (outlier detection, peptides per bin etc)
     * @param calibration_param Parameter for the m/z and im calibration (see SwathMapMassCorrection)
     * @param debug_level Debug level (writes out the RT normalization chromatograms if larger than 1)
     * @param pasef whether the data is PASEF data with possible overlapping m/z windows (with different ion mobility). In this case, the "best" SWATH window (with precursor centered around IM) is chosen.
     * @param load_into_memory Whether to cache the current SWATH map in memory
     * @param irt_trafo_out Output trafoXML file (if not empty and no input trafoXML file is given,
     *        the transformation parameters will be stored in this file)
     * @param irt_mzml_out Output Chromatogram mzML containing the iRT peptides (if not empty,
     *        iRT chromatograms will be stored in this file)
     *
     * @return CalibrationResult with: \n
     *           - rt_trafo              : the RT normalization transformation \n
     *           - ms2_mz_window_ppm     : auto-estimated MS2 m/z window (full width, ppm) \n
     *           - ms2_im_window         : auto-estimated MS2 IM window (full width, native units) \n
     *           - ms1_mz_window_ppm     : auto-estimated MS1 m/z window (full width, ppm) \n
     *           - ms1_im_window         : auto-estimated MS1 IM window (full width, native units)
     */
    CalibrationResult performCalibration(String trafo_in,
                                         const OpenSwath::LightTargetedExperiment& irt_transitions,
                                         std::vector< OpenSwath::SwathMap > & swath_maps,
                                         double min_rsq,
                                         double min_coverage,
                                         const Param& feature_finder_param,
                                         const ChromExtractParams& cp_irt,
                                         const Param& irt_detection_param,
                                         const Param& calibration_param,
                                         Size debug_level,
                                         bool pasef,
                                         bool load_into_memory,
                                         const String& irt_trafo_out,
                                         const String& irt_mzml_out);

  private:
    void loadSwathFiles_(const StringList& file_list,
                         const bool split_file,
                         const String& tmp,
                         const String& readoptions,
                         std::shared_ptr<ExperimentalSettings > & exp_meta,
                         std::vector< OpenSwath::SwathMap > & swath_maps,
                         Interfaces::IMSDataConsumer* plugin_consumer);
  }; // end TOPPOpenSwathBase
} //  end NS OpenMS
