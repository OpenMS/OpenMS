// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>

// OpenSwathCalibrationWorkflow
namespace OpenMS
{

  OpenSwath::SpectrumAccessPtr loadMS1Map(const std::vector< OpenSwath::SwathMap > & swath_maps, bool load_into_memory)
  {
    OpenSwath::SpectrumAccessPtr ms1_map;
    // store reference to MS1 map for later -> note that this is *not* threadsafe!
    for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
    {
      // if (swath_maps[i].ms1 && use_ms1_traces_)
      if (swath_maps[i].ms1)
      {
        ms1_map = swath_maps[i].sptr;
      }
    }
    if (load_into_memory)
    {
      // This creates an InMemory object that keeps all data in memory
      // but provides the same access functionality to the raw data as
      // any object implementing ISpectrumAccess
      ms1_map = boost::shared_ptr<SpectrumAccessOpenMSInMemory>( new SpectrumAccessOpenMSInMemory(*ms1_map) );
    }
    return ms1_map;
  }

  TransformationDescription OpenSwathCalibrationWorkflow::performRTNormalization(
    const OpenSwath::LightTargetedExperiment& irt_transitions,
    std::vector< OpenSwath::SwathMap > & swath_maps,
    TransformationDescription& im_trafo,
    double min_rsq,
    double min_coverage,
    const Param& feature_finder_param,
    const ChromExtractParams& cp_irt,
    const Param& irt_detection_param,
    const Param& calibration_param,
    const String& irt_mzml_out,
    Size debug_level,
    bool pasef,
    bool load_into_memory)
  {
    OPENMS_LOG_DEBUG << "performRTNormalization method starting" << std::endl;
    std::vector< OpenMS::MSChromatogram > irt_chromatograms;
    TransformationDescription trafo; // dummy
    this->simpleExtractChromatograms_(swath_maps, irt_transitions, irt_chromatograms, trafo, cp_irt, pasef, load_into_memory);

    // debug output of the iRT chromatograms
    if (irt_mzml_out.empty() && debug_level > 1)
      {
        String irt_mzml_out = "debug_irts.mzML";
      }
    if (!irt_mzml_out.empty())
    {
      try
      {
        PeakMap exp;
        exp.setChromatograms(irt_chromatograms);
        FileHandler().storeExperiment(irt_mzml_out, exp, {FileTypes::MZML});
      }
      catch (OpenMS::Exception::UnableToCreateFile& /*e*/)
      {
        OPENMS_LOG_DEBUG << "Error creating file " + irt_mzml_out + ", not writing out iRT chromatogram file"  << std::endl;
      }
      catch (OpenMS::Exception::BaseException& /*e*/)
      {
        OPENMS_LOG_DEBUG << "Error writing to file " + irt_mzml_out + ", not writing out iRT chromatogram file"  << std::endl;
      }
    }
    OPENMS_LOG_DEBUG << "Extracted number of chromatograms from iRT files: " << irt_chromatograms.size() <<  std::endl;

    // perform RT and m/z correction on the data
    TransformationDescription tr = doDataNormalization_(irt_transitions,
        irt_chromatograms, im_trafo, swath_maps,
        min_rsq, min_coverage, feature_finder_param,
        irt_detection_param, calibration_param, pasef);
    return tr;
  }

  TransformationDescription OpenSwathCalibrationWorkflow::doDataNormalization_(
    const OpenSwath::LightTargetedExperiment& targeted_exp,
    const std::vector< OpenMS::MSChromatogram >& chromatograms,
    TransformationDescription& im_trafo,
    std::vector< OpenSwath::SwathMap > & swath_maps,
    double min_rsq,
    double min_coverage,
    const Param& default_ffparam,
    const Param& irt_detection_param,
    const Param& calibration_param,
    const bool pasef)
  {
    OPENMS_LOG_DEBUG << "Start of doDataNormalization_ method" << std::endl;
    this->startProgress(0, 1, "Retention time normalization");

    bool estimateBestPeptides = irt_detection_param.getValue("estimateBestPeptides").toBool();
    if (estimateBestPeptides)
    {
      OPENMS_LOG_DEBUG << "Activated the 'estimateBestPeptides' option." << std::endl;
    }

    // 1. Estimate the retention time range of the iRT peptides over all assays
    std::pair<double,double> RTRange = OpenSwathHelper::estimateRTRange(targeted_exp);
    OPENMS_LOG_DEBUG << "Detected retention time range from " << RTRange.first << " to " << RTRange.second << std::endl;

    // 2. Store the peptide retention times in an intermediate map
    std::map<OpenMS::String, double> PeptideRTMap;
    for (Size i = 0; i < targeted_exp.getCompounds().size(); i++)
    {
      PeptideRTMap[targeted_exp.getCompounds()[i].id] = targeted_exp.getCompounds()[i].rt;
    }

    // 3. Pick input chromatograms to identify RT pairs from the input data
    const OpenSwath::LightTargetedExperiment& transition_exp_used = targeted_exp;

    // Change the feature finding parameters:
    //  - no RT score (since we don't know the correct retention time)
    //  - no RT window
    //  - no elution model score
    //  - no peak quality (use all peaks)
    //  - if best peptides should be used, use peak quality
    MRMFeatureFinderScoring featureFinder;
    Param feature_finder_param(default_ffparam);
    feature_finder_param.setValue("Scores:use_rt_score", "false");
    feature_finder_param.setValue("Scores:use_elution_model_score", "false");
    feature_finder_param.setValue("rt_extraction_window", -1.0);
    feature_finder_param.setValue("stop_report_after_feature", 1);
    feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise", 1.0); // set to 1.0 in all cases
    feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "false"); // no peak quality -> take all peaks!
    if (estimateBestPeptides)
    {
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "true");
      feature_finder_param.setValue("TransitionGroupPicker:minimal_quality", irt_detection_param.getValue("InitialQualityCutoff"));
    }
    featureFinder.setParameters(feature_finder_param);

    FeatureMap featureFile; // for results
    OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map; // for results
    std::vector<OpenSwath::SwathMap> empty_swath_maps;
    TransformationDescription empty_trafo; // empty transformation

    // Prepare the data with the chromatograms
    boost::shared_ptr<PeakMap > xic_map(new PeakMap);
    xic_map->setChromatograms(chromatograms);
    OpenSwath::SpectrumAccessPtr chromatogram_ptr = OpenSwath::SpectrumAccessPtr(new OpenMS::SpectrumAccessOpenMS(xic_map));

    featureFinder.setStrictFlag(false); // TODO remove this, it should be strict (e.g. all transitions need to be present for RT norm)
    featureFinder.pickExperiment(chromatogram_ptr, featureFile, transition_exp_used, empty_trafo, empty_swath_maps, transition_group_map);

    // 4. Find most likely correct feature for each compound and add it to the
    // "pairs" vector by computing pairs of iRT and real RT.
    //
    // Note that the quality threshold will only be applied if
    // estimateBestPeptides is true
    std::vector<std::pair<double, double> > pairs; // store the RT pairs to write the output trafoXML
    std::map<std::string, double> best_features = OpenSwathHelper::simpleFindBestFeature(transition_group_map,
      estimateBestPeptides, irt_detection_param.getValue("OverallQualityCutoff"));
    OPENMS_LOG_DEBUG << "Extracted best features: " << best_features.size() << std::endl;

    // Create pairs vector and store peaks
    std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> trgrmap_allpeaks; // store all peaks above cutoff
    for (std::map<std::string, double>::iterator it = best_features.begin(); it != best_features.end(); ++it)
    {
      pairs.emplace_back(it->second, PeptideRTMap[it->first]); // pair<exp_rt, theor_rt>
      if (transition_group_map.find(it->first) != transition_group_map.end())
      {
        trgrmap_allpeaks[ it->first ] = &transition_group_map[ it->first];
      }
    }

    // 5. Perform the outlier detection
    std::vector<std::pair<double, double> > pairs_corrected;
    String outlier_method = irt_detection_param.getValue("outlierMethod").toString();
    if (outlier_method == "iter_residual" || outlier_method == "iter_jackknife")
    {
      pairs_corrected = MRMRTNormalizer::removeOutliersIterative(pairs, min_rsq, min_coverage,
      irt_detection_param.getValue("useIterativeChauvenet").toBool(), outlier_method);
    }
    else if (outlier_method == "ransac")
    {
      // First, estimate of the maximum deviation from RT that is tolerated:
      //   Because 120 min gradient can have around 4 min elution shift, we use
      //   a default value of 3 % of the gradient to find upper RT threshold (3.6 min).
      double pcnt_rt_threshold = irt_detection_param.getValue("RANSACMaxPercentRTThreshold");
      double max_rt_threshold = (RTRange.second - RTRange.first) * pcnt_rt_threshold / 100.0;

      pairs_corrected = MRMRTNormalizer::removeOutliersRANSAC(pairs, min_rsq, min_coverage,
        irt_detection_param.getValue("RANSACMaxIterations"), max_rt_threshold,
        irt_detection_param.getValue("RANSACSamplingSize"));
    }
    else if (outlier_method == "none")
    {
      pairs_corrected = pairs;
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        String("Illegal argument '") + outlier_method +
        "' used for outlierMethod (valid: 'iter_residual', 'iter_jackknife', 'ransac', 'none').");
    }
    OPENMS_LOG_DEBUG << "Performed outlier detection, left with features: " << pairs_corrected.size() << std::endl;

    // 6. Check whether the found peptides fulfill the binned coverage criteria
    // set by the user.
    if (estimateBestPeptides)
    {
      bool enoughPeptides = MRMRTNormalizer::computeBinnedCoverage(RTRange, pairs_corrected,
        irt_detection_param.getValue("NrRTBins"),
        irt_detection_param.getValue("MinPeptidesPerBin"),
        irt_detection_param.getValue("MinBinsFilled") );

      if (!enoughPeptides)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "There were not enough bins with the minimal number of peptides");
      }
    }
    if (pairs_corrected.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "There are less than 2 iRT normalization peptides, not enough for an RT correction.");
    }

    // 7. Select the "correct" peaks for m/z (and IM) correction (e.g. remove those not
    // part of the linear regression)
    std::map<String, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType *> trgrmap_final; // store all peaks above cutoff
    for (const auto& it : trgrmap_allpeaks)
    {
      if (it.second->getFeatures().empty() ) {continue;}
      const MRMFeature& feat = it.second->getBestFeature();

      // Check if the current feature is in the list of pairs used for the
      // linear RT regression (using other features may result in wrong
      // calibration values).
      // Matching only by RT is not perfect but should work for most cases.
      for (Size pit = 0; pit < pairs_corrected.size(); pit++)
      {
        if (fabs(feat.getRT() - pairs_corrected[pit].first ) < 1e-2)
        {
          trgrmap_final[ it.first ] = it.second;
          break;
        }
      }
    }

    // 8. Correct m/z (and IM) deviations using SwathMapMassCorrection
    // m/z correction is done with the -irt_im_extraction parameters
    SwathMapMassCorrection mc;
    mc.setParameters(calibration_param);

    mc.correctMZ(trgrmap_final, targeted_exp, swath_maps, pasef);
    mc.correctIM(trgrmap_final, targeted_exp, swath_maps, pasef, im_trafo);

    // 9. store RT transformation, using the selected model
    TransformationDescription trafo_out;
    trafo_out.setDataPoints(pairs_corrected);
    Param model_params;
    model_params.setValue("symmetric_regression", "false");
    model_params.setValue("span", irt_detection_param.getValue("lowess:span"));
    model_params.setValue("num_nodes", irt_detection_param.getValue("b_spline:num_nodes"));
    String model_type = irt_detection_param.getValue("alignmentMethod").toString();
    trafo_out.fitModel(model_type, model_params);

    OPENMS_LOG_DEBUG << "Final RT mapping:" << std::endl;
    for (Size i = 0; i < pairs_corrected.size(); i++)
    {
      OPENMS_LOG_DEBUG << pairs_corrected[i].first << " " <<  pairs_corrected[i].second << std::endl;
    }
    OPENMS_LOG_DEBUG << "End of doDataNormalization_ method" << std::endl;

    this->endProgress();
    return trafo_out;
  }

  void OpenSwathCalibrationWorkflow::simpleExtractChromatograms_(
    const std::vector< OpenSwath::SwathMap > & swath_maps,
    const OpenSwath::LightTargetedExperiment& irt_transitions,
    std::vector< OpenMS::MSChromatogram > & chromatograms,
    const TransformationDescription& trafo,
    const ChromExtractParams & cp,
    bool pasef,
    bool load_into_memory)
  {
    TransformationDescription trafo_inverse = trafo;
    trafo_inverse.invert();

    // If this is pasef data, do chromatogram extraction beforehand in unparallel workflow
    std::vector<int> tr_win_map; // maps transition k to dia map i from which it should be extracted, only used if pasef flag is on
    if (pasef)
    {
      // Before calling this function, check to ensure that precursors actually have IM data
      for (Size k = 0; k < irt_transitions.transitions.size(); k++)
      {
        const OpenSwath::LightTransition& tr = irt_transitions.transitions[k];
        if (tr.getPrecursorIM() == -1)
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Transition " + tr.getNativeID() +  " does not have a valid IM value, this must be set to use the -pasef flag");
        }
      }
      OpenSwathHelper::selectSwathTransitionsPasef(irt_transitions, tr_win_map, cp.min_upper_edge_dist, swath_maps);
    }

    this->startProgress(0, 1, "Extract iRT chromatograms");
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (SignedSize map_idx = 0; map_idx < boost::numeric_cast<SignedSize>(swath_maps.size()); ++map_idx)
    {
      std::vector< OpenMS::MSChromatogram > tmp_chromatograms;
      if (!swath_maps[map_idx].ms1) // skip MS1
      {

        OpenSwath::LightTargetedExperiment transition_exp_used;

        if (pasef)
        {
          // Step 1.2: select transitions based on matching PRM/PASEF window (best window)
          std::set<std::string> matching_compounds;
          for (Size k = 0; k < tr_win_map.size(); k++)
          {
            if (tr_win_map[k] == map_idx)
            {
               const OpenSwath::LightTransition& tr = irt_transitions.transitions[k];
               transition_exp_used.transitions.push_back(tr);
               matching_compounds.insert(tr.getPeptideRef());
               OPENMS_LOG_DEBUG << "Adding Precursor with m/z " << tr.getPrecursorMZ() << " and IM of " << tr.getPrecursorIM() <<  " to swath with mz lower of " << swath_maps[map_idx].lower << " m/z upper of " << swath_maps[map_idx].upper << " im lower of " << swath_maps[map_idx].imLower << " and im upper of " << swath_maps[map_idx].imUpper << std::endl;
            }
          }

          std::set<std::string> matching_proteins;
          for (Size i = 0; i < irt_transitions.compounds.size(); i++)
          {
            if (matching_compounds.find(irt_transitions.compounds[i].id) != matching_compounds.end())
            {
              transition_exp_used.compounds.push_back( irt_transitions.compounds[i] );
              for (Size j = 0; j < irt_transitions.compounds[i].protein_refs.size(); j++)
              {
                matching_proteins.insert(irt_transitions.compounds[i].protein_refs[j]);
              }
            }
          }
          for (Size i = 0; i < irt_transitions.proteins.size(); i++)
          {
            if (matching_proteins.find(irt_transitions.proteins[i].id) != matching_proteins.end())
            {
              transition_exp_used.proteins.push_back( irt_transitions.proteins[i] );
            }
          }
        }
        else
        {
          OpenSwathHelper::selectSwathTransitions(irt_transitions, transition_exp_used,
          cp.min_upper_edge_dist, swath_maps[map_idx].lower, swath_maps[map_idx].upper);
        }
        if (!transition_exp_used.getTransitions().empty()) // skip if no transitions found
        {

          std::vector< OpenSwath::ChromatogramPtr > tmp_out;
          std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
          ChromatogramExtractor extractor;

          OpenSwath::SpectrumAccessPtr current_swath_map = swath_maps[map_idx].sptr;
          if (load_into_memory)
          {
            // This creates an InMemory object that keeps all data in memory
            current_swath_map = boost::shared_ptr<SpectrumAccessOpenMSInMemory>( new SpectrumAccessOpenMSInMemory(*current_swath_map) );
          }

          prepareExtractionCoordinates_(tmp_out, coordinates, transition_exp_used, trafo_inverse, cp);
          extractor.extractChromatograms(current_swath_map, tmp_out, coordinates, cp.mz_extraction_window,
                cp.ppm, cp.im_extraction_window, cp.extraction_function);
          extractor.return_chromatogram(tmp_out, coordinates,
              transition_exp_used, SpectrumSettings(), tmp_chromatograms, false, cp.im_extraction_window);

#ifdef _OPENMP
#pragma omp critical (osw_write_chroms)
#endif
          {
            int nr_empty_chromatograms = 0;
            OPENMS_LOG_DEBUG << "[simple] Extracted "  << tmp_chromatograms.size() << " chromatograms from SWATH map " <<
              map_idx << " with m/z " << swath_maps[map_idx].lower << " to " << swath_maps[map_idx].upper << ":" << std::endl;
            for (Size chrom_idx = 0; chrom_idx < tmp_chromatograms.size(); chrom_idx++)
            {
              // Check TIC and remove empty chromatograms (can happen if the
              // extraction window is outside the mass spectrometric acquisition
              // window).
              double tic = std::accumulate(tmp_out[chrom_idx]->getIntensityArray()->data.begin(),
                                           tmp_out[chrom_idx]->getIntensityArray()->data.end(),0.0);
              OPENMS_LOG_DEBUG << "Chromatogram "  << coordinates[chrom_idx].id << " with size "
                << tmp_out[chrom_idx]->getIntensityArray()->data.size() << " and TIC " << tic  << std::endl;
              if (tic > 0.0)
              {
                // add the chromatogram to the output
                chromatograms.push_back(tmp_chromatograms[chrom_idx]);
              }
              else
              {
                OPENMS_LOG_DEBUG << " - Warning: Empty chromatogram " << coordinates[chrom_idx].id <<
                  " detected. Will skip it!" << std::endl;
                nr_empty_chromatograms++;
              }
            }

            if (nr_empty_chromatograms > 0)
            {
              std::cerr << " - Warning: Detected " << nr_empty_chromatograms << " empty chromatograms. Will skip them!" << std::endl;
            }
          }
        }
        else
        {
          OPENMS_LOG_DEBUG << "Extracted no transitions from SWATH map " << map_idx << " with m/z " <<
              swath_maps[map_idx].lower << " to " << swath_maps[map_idx].upper << std::endl;
        }
      }
    }

    this->endProgress();
  }

  void OpenSwathCalibrationWorkflow::addChromatograms(MSChromatogram& base_chrom, const MSChromatogram& newchrom)
  {
    if (base_chrom.empty())
    {
      base_chrom = newchrom;
    }

    LinearResamplerAlign ls;
    ls.raster(newchrom.begin(), newchrom.end(), base_chrom.begin(), base_chrom.end());
  }

}

// OpenSwathWorkflow
namespace OpenMS
{

  void OpenSwathWorkflow::performExtraction(
    const std::vector< OpenSwath::SwathMap > & swath_maps,
    const TransformationDescription& trafo,
    const ChromExtractParams & cp,
    const ChromExtractParams & cp_ms1,
    const Param & feature_finder_param,
    const OpenSwath::LightTargetedExperiment& transition_exp,
    FeatureMap& out_featureFile,
    bool store_features,
    OpenSwathOSWWriter & osw_writer,
    Interfaces::IMSDataConsumer * chromConsumer,
    int batchSize,
    int ms1_isotopes,
    bool load_into_memory)
  {
    osw_writer.writeHeader();

    bool ms1_only = (swath_maps.size() == 1 && swath_maps[0].ms1);

    // Compute inversion of the transformation
    TransformationDescription trafo_inverse = trafo;
    trafo_inverse.invert();

    std::cout << "Will analyze " << transition_exp.transitions.size() << " transitions in total." << std::endl;
    int progress = 0;
    this->startProgress(0, swath_maps.size(), "Extracting and scoring transitions");

    // (i) Obtain precursor chromatograms (MS1) if precursor extraction is enabled
    ChromExtractParams ms1_cp(cp_ms1);

    if (ms1_only && !use_ms1_traces_)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Error, you need to enable use_ms1_traces when run in MS1 mode." );
    }

    if (use_ms1_traces_) ms1_map_ = loadMS1Map(swath_maps, load_into_memory);

    // (ii) Precursor extraction only
    if (ms1_only)
    {
      std::vector< MSChromatogram > ms1_chromatograms;
      MS1Extraction_(ms1_map_, swath_maps, ms1_chromatograms, ms1_cp,
                     transition_exp, trafo_inverse, ms1_only, ms1_isotopes);

      FeatureMap featureFile;
      boost::shared_ptr<MSExperiment> empty_exp = boost::shared_ptr<MSExperiment>(new MSExperiment);

      const OpenSwath::LightTargetedExperiment& transition_exp_used = transition_exp;
      scoreAllChromatograms_(std::vector<MSChromatogram>(), ms1_chromatograms, swath_maps, transition_exp_used,
                            feature_finder_param, trafo,
                            cp.rt_extraction_window, featureFile, osw_writer, ms1_isotopes, true);

      // write features to output if so desired
      std::vector< OpenMS::MSChromatogram > chromatograms;
      writeOutFeaturesAndChroms_(chromatograms, ms1_chromatograms, featureFile, out_featureFile, store_features, chromConsumer);
    }

    // (iii) map transitions to individual DIA windows for cases where this is
    // non-trivial (e.g. when there is m/z overlap and a transition could be
    // extracted from more than one window
    std::vector<int> tr_win_map; // maps transition k to dia map i from which it should be extracted
    //
    // currently not supported to do PASEF and PRM
    if (prm_ & pasef_) {
      std::cerr << "Setting -pasef and -matching_window_only flags simultaneously is not currently supported." << std::endl;
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    else if (prm_)
    {
      // Here we deal with overlapping PRM / DIA windows: we only want to extract
      // each peptide from a single window and we assume that PRM windows are
      // centered around the target peptide. We therefore select for each peptide
      // the best-matching PRM / DIA window:
      tr_win_map.resize(transition_exp.transitions.size(), -1);
      for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
      {
        for (Size k = 0; k < transition_exp.transitions.size(); k++)
        {
          const OpenSwath::LightTransition& tr = transition_exp.transitions[k];

          // If the transition falls inside the current PRM / DIA window, check
          // if the window is potentially a better match for extraction than
          // the one previously stored in the map:
          if (swath_maps[i].lower < tr.getPrecursorMZ() && tr.getPrecursorMZ() < swath_maps[i].upper &&
              std::fabs(swath_maps[i].upper - tr.getPrecursorMZ()) >= cp.min_upper_edge_dist)
          {

            if (tr_win_map[k] == -1) tr_win_map[k] = i;
            if (
                std::fabs(swath_maps[ tr_win_map[k] ].center - tr.getPrecursorMZ() ) >
                std::fabs(swath_maps[ i ].center - tr.getPrecursorMZ() ) )
            {
              // current PRM / DIA window "i" is a better match
              tr_win_map[k] = i;
            }

          }
        }
      }
    }
    else if (pasef_)
    {
      // For PASEF experiments it is possible to have DIA windows with the same m/z however different IM.
      // Extract from the DIA window in which the precursor is more centered across its IM.

      tr_win_map.resize(transition_exp.transitions.size(), -1);
      for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
      {
        for (Size k = 0; k < transition_exp.transitions.size(); k++)
        {
          const OpenSwath::LightTransition& tr = transition_exp.transitions[k];

          // If the transition falls inside the current DIA window (both in IM and m/z axis), check
          // if the window is potentially a better match for extraction than
          // the one previously stored in the map:
          if (
             swath_maps[i].imLower < tr.getPrecursorIM() && tr.getPrecursorIM() < swath_maps[i].imUpper &&
             swath_maps[i].lower < tr.getPrecursorMZ() && tr.getPrecursorMZ() < swath_maps[i].upper &&
             std::fabs(swath_maps[i].upper - tr.getPrecursorMZ()) >= cp.min_upper_edge_dist )
          {
            if (tr_win_map[k] == -1) tr_win_map[k] = i;

            // Check if the current window is better than the previously assigned window (across IM)
            double imOld = std::fabs(((swath_maps[ tr_win_map[k] ].imLower + swath_maps [ tr_win_map[k] ].imUpper) / 2) - tr.getPrecursorIM() );
            double imNew = std::fabs(((swath_maps[ i ].imLower + swath_maps [ i ].imUpper) / 2) - tr.getPrecursorIM() );
            if (imOld > imNew)
            {
              // current DIA window "i" is a better match
              OPENMS_LOG_DEBUG << "For Precursor " << tr.getPrecursorIM() << "Replacing Swath Map with IM center of " <<
                imOld << " with swath map of im center " << imNew << std::endl;
              tr_win_map[k] = i;
            }

          }
        }
      }
    }
    else {
    };

    // (iv) Perform extraction and scoring of fragment ion chromatograms (MS2)
    // We set dynamic scheduling such that the maps are worked on in the order
    // in which they were given to the program / acquired. This gives much
    // better load balancing than static allocation.
#ifdef _OPENMP
#ifdef MT_ENABLE_NESTED_OPENMP
    int total_nr_threads = omp_get_max_threads(); // store total number of threads we are allowed to use
    if (threads_outer_loop_ > -1)
    {
      std::cout << "Setting up nested loop with " << std::min(threads_outer_loop_, omp_get_max_threads()) << " threads out of "<< omp_get_max_threads() << std::endl;
      omp_set_nested(1);
      omp_set_dynamic(0);
      omp_set_num_threads(std::min(threads_outer_loop_, omp_get_max_threads()) ); // use at most threads_outer_loop_ threads here
    }
    else
    {
      std::cout << "Use non-nested loop with " << total_nr_threads << " threads." << std::endl;
    }
#endif
#pragma omp parallel for schedule(dynamic,1)
#endif
    for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
    {
      if (!swath_maps[i].ms1) // skip MS1
      {

        // Step 1: select which transitions to extract (proceed in batches)
        OpenSwath::LightTargetedExperiment transition_exp_used_all;
        if (!(prm_ || pasef_))
        {
          // Step 1.1: select transitions matching the window
          OpenSwathHelper::selectSwathTransitions(transition_exp, transition_exp_used_all,
              cp.min_upper_edge_dist, swath_maps[i].lower, swath_maps[i].upper);
        }
        else
        {
          // Step 1.2: select transitions based on matching PRM/PASEF window (best window)
          std::set<std::string> matching_compounds;
          for (Size k = 0; k < tr_win_map.size(); k++)
          {
            if (tr_win_map[k] == i)
            {
               const OpenSwath::LightTransition& tr = transition_exp.transitions[k];
               transition_exp_used_all.transitions.push_back(tr);
               matching_compounds.insert(tr.getPeptideRef());
               OPENMS_LOG_DEBUG << "Adding Precursor with m/z " << tr.getPrecursorMZ() << " and IM of " << tr.getPrecursorIM() <<  " to swath with mz upper of " << swath_maps[i].upper << " im lower of " << swath_maps[i].imLower << " and im upper of " << swath_maps[i].imUpper << std::endl;
            }
          }

          std::set<std::string> matching_proteins;
          for (Size i = 0; i < transition_exp.compounds.size(); i++)
          {
            if (matching_compounds.find(transition_exp.compounds[i].id) != matching_compounds.end())
            {
              transition_exp_used_all.compounds.push_back( transition_exp.compounds[i] );
              for (Size j = 0; j < transition_exp.compounds[i].protein_refs.size(); j++)
              {
                matching_proteins.insert(transition_exp.compounds[i].protein_refs[j]);
              }
            }
          }
          for (Size i = 0; i < transition_exp.proteins.size(); i++)
          {
            if (matching_proteins.find(transition_exp.proteins[i].id) != matching_proteins.end())
            {
              transition_exp_used_all.proteins.push_back( transition_exp.proteins[i] );
            }
          }
        }

        if (!transition_exp_used_all.getTransitions().empty()) // skip if no transitions found
        {

          OpenSwath::SpectrumAccessPtr current_swath_map = swath_maps[i].sptr;
          if (load_into_memory)
          {
            // This creates an InMemory object that keeps all data in memory
            current_swath_map = boost::shared_ptr<SpectrumAccessOpenMSInMemory>( new SpectrumAccessOpenMSInMemory(*current_swath_map) );
          }

          int batch_size;
          if (batchSize <= 0 || batchSize >= (int)transition_exp_used_all.getCompounds().size())
          {
            batch_size = transition_exp_used_all.getCompounds().size();
          }
          else
          {
            batch_size = batchSize;
          }

          SignedSize nr_batches = (transition_exp_used_all.getCompounds().size() / batch_size);

#ifdef _OPENMP
#ifdef MT_ENABLE_NESTED_OPENMP
          // If we have a multiple of threads_outer_loop_ here, then use nested
          // parallelization here. E.g. if we use 8 threads for the outer loop,
          // but we have a total of 24 cores available, each of the 8 threads
          // will then create a team of 3 threads to work on the batches
          // individually.
          //
          // We should avoid oversubscribing the CPUs, therefore we use integer division.
          // -- see https://docs.oracle.com/cd/E19059-01/stud.10/819-0501/2_nested.html
          int outer_thread_nr = omp_get_thread_num();
          omp_set_num_threads(std::max(1, total_nr_threads / threads_outer_loop_) );
#pragma omp parallel for schedule(dynamic, 1)
#endif
#endif
          for (SignedSize pep_idx = 0; pep_idx <= nr_batches; pep_idx++)
          {
            OpenSwath::SpectrumAccessPtr current_swath_map_inner = current_swath_map;

#ifdef _OPENMP
#ifdef MT_ENABLE_NESTED_OPENMP
            // To ensure multi-threading safe access to the individual spectra, we
            // need to use a light clone of the spectrum access (if multiple threads
            // share a single filestream and call seek on it, chaos will ensue).
            if (total_nr_threads / threads_outer_loop_ > 1)
            {
              current_swath_map_inner = current_swath_map->lightClone();
            }
#endif
#pragma omp critical (osw_write_stdout)
#endif
            {
              std::cout << "Thread " <<
#ifdef _OPENMP
#ifdef MT_ENABLE_NESTED_OPENMP
              outer_thread_nr << "_" << omp_get_thread_num() << " " <<
#else
              omp_get_thread_num() << "_0 " <<
#endif
#else
              "0" <<
#endif
              "will analyze " << transition_exp_used_all.getCompounds().size() <<  " compounds and "
              << transition_exp_used_all.getTransitions().size() <<  " transitions "
              "from SWATH " << i << " (batch " << pep_idx << " out of " << nr_batches << ")" << std::endl;
            }

            // Create the new, batch-size transition experiment
            OpenSwath::LightTargetedExperiment transition_exp_used;
            selectCompoundsForBatch_(transition_exp_used_all, transition_exp_used, batch_size, pep_idx);

            // Extract MS1 chromatograms for this batch
            std::vector< MSChromatogram > ms1_chromatograms;
            if (ms1_map_ != nullptr)
            {
              OpenSwath::SpectrumAccessPtr threadsafe_ms1 = ms1_map_->lightClone();
              MS1Extraction_(threadsafe_ms1, swath_maps, ms1_chromatograms, ms1_cp,
                  transition_exp_used, trafo_inverse, ms1_only, ms1_isotopes);
            }

            // Step 2.1: extract these transitions
            ChromatogramExtractor extractor;
            std::vector< OpenSwath::ChromatogramPtr > chrom_list;
            std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;

            // Step 2.2: prepare the extraction coordinates and extract chromatograms
            // chrom_list contains one entry for each fragment ion (transition) in transition_exp_used
            prepareExtractionCoordinates_(chrom_list, coordinates, transition_exp_used, trafo_inverse, cp);
            extractor.extractChromatograms(current_swath_map_inner, chrom_list, coordinates, cp.mz_extraction_window,
                cp.ppm, cp.im_extraction_window, cp.extraction_function);

            // Step 2.3: convert chromatograms back to OpenMS::MSChromatogram and write to output
            PeakMap chrom_exp;
            extractor.return_chromatogram(chrom_list, coordinates, transition_exp_used,  SpectrumSettings(),
                                          chrom_exp.getChromatograms(), false, cp.im_extraction_window);


            // Step 3: score these extracted transitions
            FeatureMap featureFile;
            std::vector< OpenSwath::SwathMap > tmp = {swath_maps[i]};
            tmp.back().sptr = current_swath_map_inner;
            scoreAllChromatograms_(chrom_exp.getChromatograms(), ms1_chromatograms, tmp, transition_exp_used,
                feature_finder_param, trafo, cp.rt_extraction_window, featureFile, osw_writer, ms1_isotopes);

            // Step 4: write all chromatograms and features out into an output object / file
            // (this needs to be done in a critical section since we only have one
            // output file and one output map).
            #pragma omp critical (osw_write_out)
            {
              writeOutFeaturesAndChroms_(chrom_exp.getChromatograms(), ms1_chromatograms, featureFile, out_featureFile, store_features, chromConsumer);
            }
          }

        } // continue 2 (no continue due to OpenMP)
      } // continue 1 (no continue due to OpenMP)

      #pragma omp critical (progress)
      this->setProgress(++progress);

    }
    this->endProgress();

#ifdef _OPENMP
#ifdef MT_ENABLE_NESTED_OPENMP
    if (threads_outer_loop_ > -1)
    {
      omp_set_num_threads(total_nr_threads); // set number of available threads back to initial value
    }
#endif
#endif
  }

  void OpenSwathWorkflow::writeOutFeaturesAndChroms_(
    std::vector< OpenMS::MSChromatogram > & chromatograms,
    std::vector< MSChromatogram >& ms1_chromatograms,
    const FeatureMap & featureFile,
    FeatureMap& out_featureFile,
    bool store_features,
    Interfaces::IMSDataConsumer * chromConsumer)
  {
    // write out MS1 chromatograms to output if so desired
    for (Size j = 0; j < ms1_chromatograms.size(); j++)
    {
      if (ms1_chromatograms[j].empty()) continue; // skip empty chromatograms
      // write MS1 chromatograms to disk
      chromConsumer->consumeChromatogram( ms1_chromatograms[j] );
    }


    // write chromatograms to output if so desired
    for (Size chrom_idx = 0; chrom_idx < chromatograms.size(); ++chrom_idx)
    {
      if (!chromatograms[chrom_idx].empty())
      {
        chromConsumer->consumeChromatogram(chromatograms[chrom_idx]);
      }
    }

    // write features to output if so desired
    if (store_features)
    {
      for (FeatureMap::const_iterator feature_it = featureFile.begin();
           feature_it != featureFile.end(); ++feature_it)
      {
        out_featureFile.push_back(*feature_it);
      }
      for (std::vector<ProteinIdentification>::const_iterator protid_it =
             featureFile.getProteinIdentifications().begin();
           protid_it != featureFile.getProteinIdentifications().end();
           ++protid_it)
      {
        out_featureFile.getProteinIdentifications().push_back(*protid_it);
      }
    }
  }

  void OpenSwathWorkflowBase::MS1Extraction_(const OpenSwath::SpectrumAccessPtr& ms1_map,
                                             const std::vector< OpenSwath::SwathMap > & /* swath_maps */,
                                             std::vector< MSChromatogram >& ms1_chromatograms,
                                             const ChromExtractParams& cp,
                                             const OpenSwath::LightTargetedExperiment& transition_exp,
                                             const TransformationDescription& trafo_inverse,
                                             bool /* ms1_only */,
                                             int ms1_isotopes)
  {
    std::vector< OpenSwath::ChromatogramPtr > chrom_list;
    std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
    OpenSwath::LightTargetedExperiment transition_exp_used = transition_exp; // copy for const correctness
    ChromatogramExtractor extractor;

    // prepare the extraction coordinates and extract chromatogram
    prepareExtractionCoordinates_(chrom_list, coordinates, transition_exp_used, trafo_inverse, cp, true, ms1_isotopes);
    extractor.extractChromatograms(ms1_map, chrom_list, coordinates, cp.mz_extraction_window,
        cp.ppm, cp.im_extraction_window, cp.extraction_function);
    extractor.return_chromatogram(chrom_list, coordinates, transition_exp_used,
        SpectrumSettings(), ms1_chromatograms, true, cp.im_extraction_window);
  }

  void OpenSwathWorkflow::scoreAllChromatograms_(
    const std::vector< OpenMS::MSChromatogram > & ms2_chromatograms,
    const std::vector< OpenMS::MSChromatogram > & ms1_chromatograms,
    const std::vector< OpenSwath::SwathMap >& swath_maps,
    const OpenSwath::LightTargetedExperiment& transition_exp,
    const Param& feature_finder_param,
    const TransformationDescription& trafo,
    const double rt_extraction_window,
    FeatureMap& output,
    OpenSwathOSWWriter & osw_writer,
    int nr_ms1_isotopes,
    bool ms1only) const
  {
    TransformationDescription trafo_inv = trafo;
    trafo_inv.invert();

    MRMFeatureFinderScoring featureFinder;
    MRMTransitionGroupPicker trgroup_picker;

    // To ensure multi-threading safe access to the individual spectra, we
    // need to use a light clone of the spectrum access (if multiple threads
    // share a single filestream and call seek on it, chaos will ensue).
    if (use_ms1_traces_ && ms1_map_)
    {
      OpenSwath::SpectrumAccessPtr threadsafe_ms1 = ms1_map_->lightClone();
      featureFinder.setMS1Map( threadsafe_ms1 );
    }
    else if (use_ms1_traces_ && !ms1_map_)
    {
      OPENMS_LOG_WARN << "WARNING: Attempted to use MS1 traces but no MS1 map was provided: Will not use MS1 signal!" << std::endl;
    }

    // If use_total_mi_score is defined, we need to instruct MRMTransitionGroupPicker to compute the score
    Param trgroup_picker_param = feature_finder_param.copy("TransitionGroupPicker:", true);
    if ((bool)feature_finder_param.getValue("Scores:use_total_mi_score").toBool())
    {
      trgroup_picker_param.setValue("compute_total_mi", "true");
    }
    trgroup_picker.setParameters(trgroup_picker_param);

    featureFinder.setParameters(feature_finder_param);
    featureFinder.prepareProteinPeptideMaps_(transition_exp);

    // Map ms1 chromatogram id to sequence number
    std::map<String, int> ms1_chromatogram_map;
    for (Size i = 0; i < ms1_chromatograms.size(); i++)
    {
      ms1_chromatogram_map[ms1_chromatograms[i].getNativeID()] = boost::numeric_cast<int>(i);
    }

    // Map chromatogram id to sequence number
    std::map<String, int> chromatogram_map;
    for (Size i = 0; i < ms2_chromatograms.size(); i++)
    {
      chromatogram_map[ms2_chromatograms[i].getNativeID()] = boost::numeric_cast<int>(i);
    }
    // Map peptide id to sequence number
    std::map<String, int> assay_peptide_map;
    for (Size i = 0; i < transition_exp.getCompounds().size(); i++)
    {
      assay_peptide_map[transition_exp.getCompounds()[i].id] = boost::numeric_cast<int>(i);
    }

    // Map peptide id to corresponding transitions
    typedef std::map<String, std::vector< const TransitionType* > > AssayMapT;
    AssayMapT assay_map;
    // create an entry for each member (ensure there is one even if we don't
    // have any transitions for it, e.g. in the case of ms1 only)
    for (Size i = 0; i < transition_exp.getCompounds().size(); i++)
    {
      assay_map[transition_exp.getCompounds()[i].id] = std::vector< const TransitionType* >();
    }
    for (Size i = 0; i < transition_exp.getTransitions().size(); i++)
    {
      assay_map[transition_exp.getTransitions()[i].getPeptideRef()].push_back(&transition_exp.getTransitions()[i]);
    }

    std::vector<String> to_osw_output;
    ///////////////////////////////////
    // Start of main function
    // Iterating over all the assays
    ///////////////////////////////////
    for (AssayMapT::iterator assay_it = assay_map.begin(); assay_it != assay_map.end(); ++assay_it)
    {
      // Create new MRMTransitionGroup
      String id = assay_it->first;
      MRMTransitionGroupType transition_group;
      transition_group.setTransitionGroupID(id);
      double expected_rt = transition_exp.getCompounds()[ assay_peptide_map[id] ].rt;

      // 1. Go through all transitions, for each transition get
      // the chromatogram and the assay to the MRMTransitionGroup
      const TransitionType* detection_assay_it = nullptr; // store last detecting transition
      for (const TransitionType* transition : assay_it->second)
      {
        if (transition->isDetectingTransition())
        {
          detection_assay_it = transition;
        }

        // continue if we only have MS1 (we wont have any chromatograms for
        // the transitions)
        if (ms1only) {continue;}

        if (chromatogram_map.find(transition->getNativeID()) == chromatogram_map.end())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Error, did not find chromatogram for transition " + transition->getNativeID() );
        }

        // Convert chromatogram to MSChromatogram and filter
        auto chromatogram = ms2_chromatograms[ chromatogram_map[transition->getNativeID()] ];
        chromatogram.setNativeID(transition->getNativeID());
        if (rt_extraction_window > 0)
        {
          double de_normalized_experimental_rt = trafo_inv.apply(expected_rt);
          double rt_max = de_normalized_experimental_rt + rt_extraction_window;
          double rt_min = de_normalized_experimental_rt - rt_extraction_window;
          chromatogram.erase(std::remove_if(chromatogram.begin(), chromatogram.end(),
                                            [rt_min, rt_max](const ChromatogramPeak& chr)
                                            { return chr.getRT() > rt_max || chr.getRT() < rt_min; })
                             , chromatogram.end());
        }

        // Add the transition and the chromatogram to the MRMTransitionGroup
        transition_group.addTransition(*transition, transition->getNativeID());
        transition_group.addChromatogram(chromatogram, chromatogram.getNativeID());
      }

      // currently  .osw and .featureXML are mutually exclusive
      if (osw_writer.isActive()) { output.clear(); }

      // 2. Set the MS1 chromatograms for the different isotopes, if available
      // (note that for 3 isotopes, we include the monoisotopic peak plus three
      // isotopic traces)
      for (int iso = 0; iso <= nr_ms1_isotopes; iso++)
      {
        String prec_id = OpenSwathHelper::computePrecursorId(transition_group.getTransitionGroupID(), iso);
        if (!ms1_chromatograms.empty() && ms1_chromatogram_map.find(prec_id) != ms1_chromatogram_map.end())
        {
          MSChromatogram chromatogram = ms1_chromatograms[ ms1_chromatogram_map[prec_id] ];
          transition_group.addPrecursorChromatogram(chromatogram, chromatogram.getNativeID());
        }
      }

      // 3. / 4. Process the MRMTransitionGroup: find peakgroups and score them
      trgroup_picker.pickTransitionGroup(transition_group);
      featureFinder.scorePeakgroups(transition_group, trafo, swath_maps, output, ms1only);

      // Ensure that a detection transition is used to derive features for output
      if (detection_assay_it == nullptr && !output.empty())
      {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Error, did not find any detection transition for feature " + id );
      }

      // 5. Add to the output osw if given
      if (osw_writer.isActive() && !output.empty()) // implies that detection_assay_it was set
      {
        const OpenSwath::LightCompound pep;
        to_osw_output.push_back(osw_writer.prepareLine(OpenSwath::LightCompound(), // not used currently: transition_exp.getCompounds()[ assay_peptide_map[id] ],
                                                       nullptr, // not used currently: detection_assay_it,
                                                       output,
                                                       id));
      }
    }

    // Only write at the very end since this is a step that needs a barrier
    if (osw_writer.isActive())
    {
#ifdef _OPENMP
#pragma omp critical (osw_write_tsv)
#endif
      {
        osw_writer.writeLines(to_osw_output);
      }
    }
  }


  void OpenSwathWorkflow::selectCompoundsForBatch_(const OpenSwath::LightTargetedExperiment& transition_exp_used_all,
    OpenSwath::LightTargetedExperiment& transition_exp_used, int batch_size, size_t j)
  {
    // compute batch start/end
    size_t start = j * batch_size;
    size_t end = j * batch_size + batch_size;
    if (end > transition_exp_used_all.compounds.size())
    {
      end = transition_exp_used_all.compounds.size();
    }

    // Create the new, batch-size transition experiment
    transition_exp_used.proteins = transition_exp_used_all.proteins;
    transition_exp_used.compounds.insert(transition_exp_used.compounds.end(),
        transition_exp_used_all.compounds.begin() + start, transition_exp_used_all.compounds.begin() + end);
    copyBatchTransitions_(transition_exp_used.compounds, transition_exp_used_all.transitions, transition_exp_used.transitions);
  }

  void OpenSwathWorkflow::copyBatchTransitions_(const std::vector<OpenSwath::LightCompound>& used_compounds,
    const std::vector<OpenSwath::LightTransition>& all_transitions,
    std::vector<OpenSwath::LightTransition>& output)
  {
    std::set<std::string> selected_compounds;
    for (Size i = 0; i < used_compounds.size(); i++)
    {
      selected_compounds.insert(used_compounds[i].id);
    }

    for (Size i = 0; i < all_transitions.size(); i++)
    {
      if (selected_compounds.find(all_transitions[i].peptide_ref) != selected_compounds.end())
      {
        output.push_back(all_transitions[i]);
      }
    }
  }

  void OpenSwathWorkflowBase::prepareExtractionCoordinates_(std::vector< OpenSwath::ChromatogramPtr > & chrom_list,
                                                            std::vector< ChromatogramExtractorAlgorithm::ExtractionCoordinates > & coordinates,
                                                            const OpenSwath::LightTargetedExperiment & transition_exp_used,
                                                            const TransformationDescription& trafo_inverse,
                                                            const ChromExtractParams & cp,
                                                            const bool ms1,
                                                            const int ms1_isotopes) const
  {
    if (cp.rt_extraction_window < 0)
    {
      ChromatogramExtractor::prepare_coordinates(chrom_list, coordinates, transition_exp_used, cp.rt_extraction_window, ms1, ms1_isotopes);
    }
    else
    {
      // Use an rt extraction window of 0.0 which will just write the retention time in start / end positions
      // Then correct the start/end positions and add the extra_rt_extract parameter
      ChromatogramExtractor::prepare_coordinates(chrom_list, coordinates, transition_exp_used, 0.0, ms1, ms1_isotopes);
      for (std::vector< ChromatogramExtractor::ExtractionCoordinates >::iterator it = coordinates.begin(); it != coordinates.end(); ++it)
      {
        it->rt_start = trafo_inverse.apply(it->rt_start) - (cp.rt_extraction_window + cp.extra_rt_extract)/ 2.0;
        it->rt_end = trafo_inverse.apply(it->rt_end) + (cp.rt_extraction_window + cp.extra_rt_extract)/ 2.0;
      }
    }
  }
}