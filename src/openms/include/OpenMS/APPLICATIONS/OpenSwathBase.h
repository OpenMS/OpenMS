// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest$
// $Authors: Hannes Roest$
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

class TOPPOpenSwathBase :
  public TOPPBase
{

public:

  TOPPOpenSwathBase(String name, String description, bool official = true) :
    TOPPBase(name, description, official)
  {
  }

private:

  void loadSwathFiles_(const StringList& file_list,
                       const bool split_file,
                       const String& tmp,
                       const String& readoptions,
                       boost::shared_ptr<ExperimentalSettings > & exp_meta,
                       std::vector< OpenSwath::SwathMap > & swath_maps,
                       Interfaces::IMSDataConsumer* plugin_consumer)
  {
    SwathFile swath_file;
    swath_file.setLogType(log_type_);

    if (split_file || file_list.size() > 1)
    {
      // TODO cannot use data reduction here any more ...
      swath_maps = swath_file.loadSplit(file_list, tmp, exp_meta, readoptions);
    }
    else
    {
      FileTypes::Type in_file_type = FileHandler::getTypeByFileName(file_list[0]);
      if (in_file_type == FileTypes::MZML)
      {
        swath_maps = swath_file.loadMzML(file_list[0], tmp, exp_meta, readoptions, plugin_consumer);
      }
      else if (in_file_type == FileTypes::MZXML)
      {
        swath_maps = swath_file.loadMzXML(file_list[0], tmp, exp_meta, readoptions);
      }
      else if (in_file_type == FileTypes::SQMASS)
      {
        swath_maps = swath_file.loadSqMass(file_list[0], exp_meta);
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Input file needs to have ending mzML or mzXML");
      }
    }
  }

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
                      boost::shared_ptr<ExperimentalSettings >& exp_meta,
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
                      Interfaces::IMSDataConsumer* plugin_consumer = nullptr)
  {
    // (i) Load files
    loadSwathFiles_(file_list, split_file, tmp, readoptions, exp_meta, swath_maps, plugin_consumer);

    // (ii) Allow the user to specify the SWATH windows
    if (!swath_windows_file.empty())
    {
      SwathWindowLoader::annotateSwathMapsFromFile(swath_windows_file, swath_maps, sort_swath_maps, force);
    }

    for (Size i = 0; i < swath_maps.size(); i++)
    {
      OPENMS_LOG_DEBUG << "Found swath map " << i
        << " with lower " << swath_maps[i].lower
        << " and upper " << swath_maps[i].upper
        << " and im Lower bounds of " << swath_maps[i].imLower
        << " and im Upper bounds of " << swath_maps[i].imUpper
        << " and " << swath_maps[i].sptr->getNrSpectra()
        << " spectra." << std::endl;
    }

    // (iii) Sanity check: there should be no overlap between the windows:
    std::vector<std::pair<double, double>> sw_windows;
    for (Size i = 0; i < swath_maps.size(); i++)
    {
      if (!swath_maps[i].ms1)
      {
        sw_windows.push_back(std::make_pair(swath_maps[i].lower, swath_maps[i].upper));
      }
    }
    // sort by lower bound (first entry in pair)
    std::sort(sw_windows.begin(), sw_windows.end());

    for (Size i = 1; i < sw_windows.size(); i++)
    {
      double lower_map_end = sw_windows[i-1].second - min_upper_edge_dist;
      double upper_map_start = sw_windows[i].first;
      OPENMS_LOG_DEBUG << "Extraction will go up to " << lower_map_end << " and continue at " << upper_map_start << std::endl;

      if (prm) {continue;} // skip next step as expect them to overlap and have gaps...

      if (upper_map_start - lower_map_end > 0.01)
      {
        OPENMS_LOG_WARN << "Extraction will have a gap between " << lower_map_end << " and " << upper_map_start << std::endl;
        if (!force)
        {
          OPENMS_LOG_ERROR << "Extraction windows have a gap. Will abort (override with -force)" << std::endl;
          return false;
        }
      }

      if (pasef) {continue;} // skip this step, expect there to be overlap ...

      if (lower_map_end - upper_map_start > 0.01)
      {
        OPENMS_LOG_WARN << "Extraction will overlap between " << lower_map_end << " and " << upper_map_start << "!\n"
                 << "This will lead to multiple extraction of the transitions in the overlapping region "
                 << "which will lead to duplicated output. It is very unlikely that you want this." << "\n"
                 << "Please fix this by providing an appropriate extraction file with -swath_windows_file" << "\n"
                 << "Did you mean to set the -pasef Flag?" << std::endl;
        if (!force)
        {
          OPENMS_LOG_ERROR << "Extraction windows overlap. Will abort (override with -force)" << std::endl;
          return false;
        }
      }
    }
    return true;
  }

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
                          const boost::shared_ptr<ExperimentalSettings>& exp_meta,
                          const OpenSwath::LightTargetedExperiment& transition_exp,
                          const String& out_chrom,
                          const UInt64 run_id)
  {
    if (!out_chrom.empty())
    {
      String tmp = out_chrom;
      if (tmp.toLower().hasSuffix(".sqmass"))
      {
        bool full_meta = false; // can lead to very large files in memory
        bool lossy_compression = true;
        *chromatogramConsumer = new MSDataSqlConsumer(out_chrom, run_id, 500, full_meta, lossy_compression);
      }
      else
      {
        PlainMSDataWritingConsumer * chromConsumer = new PlainMSDataWritingConsumer(out_chrom);
        int expected_chromatograms = transition_exp.transitions.size();
        chromConsumer->setExpectedSize(0, expected_chromatograms);
        chromConsumer->setExperimentalSettings(*exp_meta);
        chromConsumer->getOptions().setWriteIndex(true);  // ensure that we write the index
        chromConsumer->addDataProcessing(getProcessingInfo_(DataProcessing::SMOOTHING));

        // prepare data structures for lossy compression
        MSNumpressCoder::NumpressConfig npconfig_mz;
        MSNumpressCoder::NumpressConfig npconfig_int;
        npconfig_mz.estimate_fixed_point = true; // critical
        npconfig_int.estimate_fixed_point = true; // critical
        npconfig_mz.numpressErrorTolerance = -1.0; // skip check, faster
        npconfig_int.numpressErrorTolerance = -1.0; // skip check, faster
        npconfig_mz.setCompression("linear");
        npconfig_int.setCompression("slof");
        npconfig_mz.linear_fp_mass_acc = 0.05; // set the desired RT accuracy in seconds

        chromConsumer->getOptions().setNumpressConfigurationMassTime(npconfig_mz);
        chromConsumer->getOptions().setNumpressConfigurationIntensity(npconfig_int);
        chromConsumer->getOptions().setCompression(true);

        *chromatogramConsumer = chromConsumer;
      }
    }
    else
    {
      *chromatogramConsumer = new NoopMSDataWritingConsumer("");
    }
  }

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
                                                        const Param& tsv_reader_param)
  {
    OpenSwath::LightTargetedExperiment transition_exp;
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    if (tr_type == FileTypes::TRAML)
    {
      progresslogger.startProgress(0, 1, "Load TraML file");
      TargetedExperiment targeted_exp;
      FileHandler().loadTransitions(tr_file, targeted_exp, {FileTypes::TRAML});
      OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, transition_exp);
      progresslogger.endProgress();
    }
    else if (tr_type == FileTypes::PQP)
    {
      progresslogger.startProgress(0, 1, "Load PQP file");
      TransitionPQPFile().convertPQPToTargetedExperiment(tr_file.c_str(), transition_exp);
      progresslogger.endProgress();
    }
    else if (tr_type == FileTypes::TSV)
    {
      progresslogger.startProgress(0, 1, "Load TSV file");
      TransitionTSVFile tsv_reader;
      tsv_reader.setParameters(tsv_reader_param);
      tsv_reader.convertTSVToTargetedExperiment(tr_file.c_str(), tr_type, transition_exp);
      progresslogger.endProgress();
    }
    else
    {
      OPENMS_LOG_ERROR << "Provide valid TraML, TSV or PQP transition file." << std::endl;
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Need to provide valid input file.");
    }
    return transition_exp;
  }

  /**
   * @brief Perform retention time and m/z calibration
   *
   * This function will create the retention time transformation either by
   * loading a provided .trafoXML file or determine it from the data itself by
   * extracting the transitions specified in the irt_tr_file TraML file. It
   * will also perform the m/z calibration.
   *
   * @note Internally, the retention time and @p m/z calibration are performed
   * by OpenMS::OpenSwathCalibrationWorkflow::performRTNormalization
   *
   * @param trafo_in Input trafoXML file (if not empty, transformation will be
   *                 loaded from this file)
   * @param irt_tr_file  Input TraML file containing transitions (if trafo_in
   *                     is empty, this file will be loaded and transitions
   *                     will be extracted)
   * @param swath_maps The raw data (swath maps)
   * @param min_rsq Minimal R^2 value that is expected for the RT regression
   * @param min_coverage Minimal coverage of the chromatographic space that needs to be achieved
   * @param feature_finder_param Parameter set for the feature finding in chromatographic dimension
   * @param cp_irt Parameter set for the chromatogram extraction
   * @param irt_detection_param Parameter set for the detection of the iRTs (outlier detection, peptides per bin etc)
   * @param calibration_param Parameter for the m/z and im calibration (see SwathMapMassCorrection)
   * @param debug_level Debug level (writes out the RT normalization chromatograms if larger than 1)
   * @param pasef whether the data is PASEF data with possible overlapping m/z windows (with different ion mobility). In this case, the "best" SWATH window (with precursor cetntered around IM) is chosen.
   * @param load_into_memory Whether to cache the current SWATH map in memory
   * @param irt_trafo_out Output trafoXML file (if not empty and no input trafoXML file is given,
   *        the transformation parameters will be stored in this file)
   * @param irt_mzml_out Output Chromatogram mzML containing the iRT peptides (if not empty,
   *        iRT chromatograms will be stored in this file)
   *
   */
  TransformationDescription performCalibration(String trafo_in,
        String irt_tr_file,
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
        const String& irt_mzml_out)
  {
    TransformationDescription trafo_rtnorm;

    if (!trafo_in.empty())
    {
      // get read RT normalization file
      FileHandler().loadTransformations(trafo_in, trafo_rtnorm, false, {FileTypes::TRANSFORMATIONXML});
      Param model_params = getParam_().copy("model:", true);
      model_params.setValue("symmetric_regression", "false");
      model_params.setValue("span", irt_detection_param.getValue("lowess:span"));
      model_params.setValue("num_nodes", irt_detection_param.getValue("b_spline:num_nodes"));
      String model_type = irt_detection_param.getValue("alignmentMethod").toString();
      trafo_rtnorm.fitModel(model_type, model_params);
    }
    else if (!irt_tr_file.empty())
    {
      // Loading iRT file
      std::cout << "Will load iRT transitions and try to find iRT peptides" << std::endl;
      FileTypes::Type tr_type = FileHandler::getType(irt_tr_file);
      Param tsv_reader_param = TransitionTSVFile().getDefaults();
      OpenSwath::LightTargetedExperiment irt_transitions = loadTransitionList(tr_type, irt_tr_file, tsv_reader_param);

      // If pasef flag is set, validate that IM is present
      if (pasef)
      {
        const auto& transitions = irt_transitions.getTransitions();

        for ( Size k=0; k < (Size)transitions.size(); k++ )
        {
          if (transitions[k].precursor_im == -1)
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: iRT Transition " + transitions[k].getNativeID() +  " does not have a valid IM value, this must be set to use the -pasef flag");
          }
        }
      }

      // perform extraction
      OpenSwathCalibrationWorkflow wf;
      wf.setLogType(log_type_);
      TransformationDescription im_trafo;
      trafo_rtnorm = wf.performRTNormalization(irt_transitions, swath_maps, im_trafo,
                                               min_rsq, min_coverage,
                                               feature_finder_param,
                                               cp_irt, irt_detection_param,
                                               calibration_param, irt_mzml_out, debug_level, pasef,
                                               load_into_memory);

      if (!irt_trafo_out.empty())
      {
        FileHandler().storeTransformations(irt_trafo_out, trafo_rtnorm, {FileTypes::TRANSFORMATIONXML});
      }
    }
    return trafo_rtnorm;
  }


};

}
