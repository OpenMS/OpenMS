// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest, Justin Sing$
// $Authors: Hannes Roest, Justin Sing$
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/OpenSwathBase.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h>

#include <algorithm>
#include <iostream>
#include <utility>

using namespace std;

namespace OpenMS
{

  TOPPOpenSwathBase::TOPPOpenSwathBase(String name, String description, bool official, const std::vector<Citation>& citations) :
    TOPPBase(name, description, official, citations)
  {
  }

  TOPPOpenSwathBase::~TOPPOpenSwathBase() = default;

  // Private

  void TOPPOpenSwathBase::loadSwathFiles_(const StringList& file_list,
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

  // Protected

  bool TOPPOpenSwathBase::loadSwathFiles(const StringList& file_list,
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
                      Interfaces::IMSDataConsumer* plugin_consumer)
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

  void TOPPOpenSwathBase::prepareChromOutput(Interfaces::IMSDataConsumer ** chromatogramConsumer,
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

  OpenSwath::LightTargetedExperiment TOPPOpenSwathBase::loadTransitionList(const FileTypes::Type& tr_type,
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

  TOPPOpenSwathBase::CalibrationResult TOPPOpenSwathBase::performCalibration(String trafo_in,
                                                                             const OpenSwath::LightTargetedExperiment& irt_transitions,
                                                                             std::vector<OpenSwath::SwathMap>& swath_maps,
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
    double auto_mz_w = 0;
    double auto_im_w = -1;
    double auto_ms1_mz_w = 0;
    double auto_ms1_im_w = -1;

    if (! trafo_in.empty())
    {
      // get read RT normalization file
      FileHandler().loadTransformations(trafo_in, trafo_rtnorm, false, {FileTypes::TRANSFORMATIONXML});
      Param model_params = getParam_().copy("model:", true);
      model_params.setValue("symmetric_regression", "false");
      model_params.setValue("span", irt_detection_param.getValue("lowess:span"));
      model_params.setValue("num_nodes", irt_detection_param.getValue("b_spline:num_nodes"));
      String model_type = irt_detection_param.getValue("alignmentMethod").toString();
      trafo_rtnorm.fitModel(model_type, model_params);
      // We don't calibrate for mz and IM if a user supplies a transformation function
      // TODO: Should we deprecate and remove the option of providing a transformation function?
      //  Not sure how often this is used in practice. @singjc, 2025-08-01
      auto_mz_w = -1.0;
      auto_im_w = -1.0;
      auto_ms1_mz_w = -1.0;
      auto_ms1_im_w = -1.0;
    }
    else if (! irt_transitions.getTransitions().empty())
    {
      // Loading iRT file
      std::cout << "Will load iRT transitions and try to find iRT peptides" << std::endl;

      // If pasef flag is set, validate that IM is present
      if (pasef)
      {
        const auto& transitions = irt_transitions.getTransitions();

        for (Size k = 0; k < (Size)transitions.size(); k++)
        {
          if (transitions[k].precursor_im == -1)
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                             "Error: iRT Transition " + transitions[k].getNativeID()
                                               + " does not have a valid IM value, this must be set to use the -pasef flag");
          }
        }
      }

      // perform extraction
      OpenSwathCalibrationWorkflow wf;
      wf.setLogType(log_type_);
      TransformationDescription im_trafo;
      trafo_rtnorm = wf.performRTNormalization(irt_transitions, swath_maps, im_trafo, min_rsq, min_coverage, feature_finder_param, cp_irt,
                                               irt_detection_param, calibration_param, irt_mzml_out, debug_level, pasef, load_into_memory);
      // Retrieve estimated mz and IM extraction windows
      auto_mz_w = wf.getEstimatedMzWindow();
      auto_im_w = wf.getEstimatedImWindow();
      // Retrieve estimate MS1 mz and IM extraction windows
      auto_ms1_mz_w = wf.getEstimatedMs1MzWindow();
      auto_ms1_im_w = wf.getEstimatedMs1ImWindow();

      if (! irt_trafo_out.empty()) { FileHandler().storeTransformations(irt_trafo_out, trafo_rtnorm, {FileTypes::TRANSFORMATIONXML}); }
    }

    CalibrationResult out;
    out.rt_trafo = std::move(trafo_rtnorm);
    out.ms2_mz_window_ppm = auto_mz_w;
    out.ms2_im_window = auto_im_w;
    out.ms1_mz_window_ppm = auto_ms1_mz_w;
    out.ms1_im_window = auto_ms1_im_w;
    return out;
  }
} //  end NS OpenMS