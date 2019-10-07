// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    bool sonar,
    bool load_into_memory)
  {
    OPENMS_LOG_DEBUG << "performRTNormalization method starting" << std::endl;
    std::vector< OpenMS::MSChromatogram > irt_chromatograms;
    TransformationDescription trafo; // dummy
    this->simpleExtractChromatograms_(swath_maps, irt_transitions, irt_chromatograms, trafo, cp_irt, sonar, load_into_memory);

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
        MzMLFile().store(irt_mzml_out, exp);
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
        irt_detection_param, calibration_param);
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
    const Param& calibration_param)
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
    OpenSwath::LightTargetedExperiment transition_exp_used = targeted_exp;

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
    feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise", 1.0); // set to 1.0 in all cases
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
      pairs.push_back(std::make_pair(it->second, PeptideRTMap[it->first])); // pair<exp_rt, theor_rt>
      if (transition_group_map.find(it->first) != transition_group_map.end())
      {
        trgrmap_allpeaks[ it->first ] = &transition_group_map[ it->first];
      }
    }

    // 5. Perform the outlier detection
    std::vector<std::pair<double, double> > pairs_corrected;
    String outlier_method = irt_detection_param.getValue("outlierMethod");
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

    // 7. Select the "correct" peaks for m/z correction (e.g. remove those not
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

    // 8. Correct m/z deviations using SwathMapMassCorrection
    SwathMapMassCorrection mc;
    mc.setParameters(calibration_param);
    mc.correctMZ(trgrmap_final, swath_maps, targeted_exp);
    mc.correctIM(trgrmap_final, swath_maps, im_trafo, targeted_exp);

    // 9. store RT transformation, using the selected model
    TransformationDescription trafo_out;
    trafo_out.setDataPoints(pairs_corrected);
    Param model_params;
    model_params.setValue("symmetric_regression", "false");
    model_params.setValue("span", irt_detection_param.getValue("lowess:span"));
    model_params.setValue("num_nodes", irt_detection_param.getValue("b_spline:num_nodes"));
    String model_type = irt_detection_param.getValue("alignmentMethod");
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
    bool sonar,
    bool load_into_memory)
  {
    TransformationDescription trafo_inverse = trafo;
    trafo_inverse.invert();

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
        OpenSwathHelper::selectSwathTransitions(irt_transitions, transition_exp_used,
            cp.min_upper_edge_dist, swath_maps[map_idx].lower, swath_maps[map_idx].upper);
        if (transition_exp_used.getTransitions().size() > 0) // skip if no transitions found
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

    if (sonar)
    {

      OPENMS_LOG_DEBUG << " got a total of " << chromatograms.size() << " chromatograms before SONAR addition " << std::endl;

      // for SONAR: group chromatograms together and then add them up (we will have one chromatogram for every single map)
      std::vector< OpenMS::MSChromatogram > chromatograms_new;
      std::map<std::string, std::vector<int> > chr_map;
      for (Size i = 0; i < chromatograms.size(); i++)
      {
        chr_map[ chromatograms[i].getNativeID() ].push_back(i);
      }

      for (std::map<std::string, std::vector<int> >::iterator it = chr_map.begin(); it != chr_map.end(); ++it)
      {
        MSChromatogram chrom_acc; // accumulator
        for (Size i = 0; i < it->second.size(); i++)
        {
          addChromatograms(chrom_acc, chromatograms[ it->second[i] ] );
        }
        chromatograms_new.push_back(chrom_acc);
      }
      chromatograms = chromatograms_new; // switch

      OPENMS_LOG_DEBUG << " got a total of " << chromatograms.size() << " chromatograms after SONAR addition " << std::endl;
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
    const TransformationDescription trafo,
    const ChromExtractParams & cp,
    const ChromExtractParams & cp_ms1,
    const Param & feature_finder_param,
    const OpenSwath::LightTargetedExperiment& transition_exp,
    FeatureMap& out_featureFile,
    bool store_features,
    OpenSwathTSVWriter & tsv_writer,
    OpenSwathOSWWriter & osw_writer,
    Interfaces::IMSDataConsumer * chromConsumer,
    int batchSize,
    int ms1_isotopes,
    bool load_into_memory)
  {
    tsv_writer.writeHeader();
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
    if (!use_ms1_ion_mobility_)
    {
      ms1_cp.im_extraction_window = -1;
    }

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
      MS1Extraction_(ms1_map_, swath_maps, ms1_chromatograms, chromConsumer, ms1_cp,
                     transition_exp, trafo_inverse, ms1_only, ms1_isotopes);

      FeatureMap featureFile;
      boost::shared_ptr<MSExperiment> empty_exp = boost::shared_ptr<MSExperiment>(new MSExperiment);

      OpenSwath::LightTargetedExperiment transition_exp_used = transition_exp;
      scoreAllChromatograms_(std::vector<MSChromatogram>(), ms1_chromatograms, swath_maps, transition_exp_used, 
                            feature_finder_param, trafo,
                            cp.rt_extraction_window, featureFile, tsv_writer, osw_writer, ms1_isotopes, true);

      // write features to output if so desired
      std::vector< OpenMS::MSChromatogram > chromatograms;
      writeOutFeaturesAndChroms_(chromatograms, featureFile, out_featureFile, store_features, chromConsumer);
    }

    std::vector<int> prm_map;
    if (prm_)
    {
      // Here we deal with overlapping PRM / DIA windows: we only want to extract
      // each peptide from a single window and we assume that PRM windows are
      // centered around the target peptide. We therefore select for each peptide
      // the best-matching PRM / DIA window:
      prm_map.resize(transition_exp.transitions.size(), -1);
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

            if (prm_map[k] == -1) prm_map[k] = i;
            if (
                std::fabs(swath_maps[ prm_map[k] ].center - tr.getPrecursorMZ() ) > 
                std::fabs(swath_maps[ i ].center - tr.getPrecursorMZ() ) )
            {
              // current PRM / DIA window "i" is a better match
              prm_map[k] = i;
            }

          }
        }
      }
    }

    // (iii) Perform extraction and scoring of fragment ion chromatograms (MS2)
    // We set dynamic scheduling such that the maps are worked on in the order
    // in which they were given to the program / acquired. This gives much
    // better load balancing than static allocation.
#ifdef _OPENMP
#ifdef MT_ENABLE_NESTED_OPENMP
    int total_nr_threads = omp_get_max_threads(); // store total number of threads we are allowed to use
    if (threads_outer_loop_ > -1)
    {
      omp_set_nested(1);
      omp_set_dynamic(0);
      omp_set_num_threads(std::min(threads_outer_loop_, omp_get_max_threads()) ); // use at most threads_outer_loop_ threads here
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
        if (!prm_)
        {
          // Step 1.1: select transitions matching the window
          OpenSwathHelper::selectSwathTransitions(transition_exp, transition_exp_used_all,
              cp.min_upper_edge_dist, swath_maps[i].lower, swath_maps[i].upper);
        }
        else
        {
          // Step 1.2: select transitions based on matching PRM window (best window)
          std::set<std::string> matching_compounds;
          for (Size k = 0; k < prm_map.size(); k++)
          {
            if (prm_map[k] == i)
            {
               const OpenSwath::LightTransition& tr = transition_exp.transitions[k];
               transition_exp_used_all.transitions.push_back(tr);
               matching_compounds.insert(tr.getPeptideRef());
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

        if (transition_exp_used_all.getTransitions().size() > 0) // skip if no transitions found
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
              MS1Extraction_(threadsafe_ms1, swath_maps, ms1_chromatograms, chromConsumer, ms1_cp,
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
                feature_finder_param, trafo, cp.rt_extraction_window, featureFile, tsv_writer, osw_writer, ms1_isotopes);

            // Step 4: write all chromatograms and features out into an output object / file
            // (this needs to be done in a critical section since we only have one
            // output file and one output map).
            #pragma omp critical (osw_write_out)
            {
              writeOutFeaturesAndChroms_(chrom_exp.getChromatograms(), featureFile, out_featureFile, store_features, chromConsumer);
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
    const FeatureMap & featureFile,
    FeatureMap& out_featureFile,
    bool store_features,
    Interfaces::IMSDataConsumer * chromConsumer)
  {
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

  void OpenSwathWorkflowBase::MS1Extraction_(const OpenSwath::SpectrumAccessPtr ms1_map,
                                             const std::vector< OpenSwath::SwathMap > & /* swath_maps */,
                                             std::vector< MSChromatogram >& ms1_chromatograms,
                                             Interfaces::IMSDataConsumer* chromConsumer,
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

    for (Size j = 0; j < coordinates.size(); j++)
    {
      if (ms1_chromatograms[j].empty()) continue; // skip empty chromatograms

#ifdef _OPENMP
#pragma omp critical (osw_write_out)
#endif
      {
        // write MS1 chromatograms to disk
        chromConsumer->consumeChromatogram( ms1_chromatograms[j] );
      }
    } // end of for coordinates

  }

  void OpenSwathWorkflow::scoreAllChromatograms_(
    const std::vector< OpenMS::MSChromatogram > & ms2_chromatograms,
    const std::vector< OpenMS::MSChromatogram > & ms1_chromatograms,
    const std::vector< OpenSwath::SwathMap >& swath_maps,
    const OpenSwath::LightTargetedExperiment& transition_exp,
    const Param& feature_finder_param,
    TransformationDescription trafo,
    const double rt_extraction_window,
    FeatureMap& output, 
    OpenSwathTSVWriter & tsv_writer,
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
    if (use_ms1_traces_)
    {
      if (ms1_map_ == nullptr) 
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Error, attempted to use MS1 traces, but no MS1 map was provided." );
      }
      OpenSwath::SpectrumAccessPtr threadsafe_ms1 = ms1_map_->lightClone();
      featureFinder.setMS1Map( threadsafe_ms1 );
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

    std::vector<String> to_tsv_output, to_osw_output;
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

      // 1. Go through all transitions, for each transition get chromatogram
      // and the chromatogram and the assay to the MRMTransitionGroup
      int detection_assay_it = -1; // store index for the last detection transition
      for (Size i = 0; i < assay_it->second.size(); i++)
      {
        const TransitionType* transition = assay_it->second[i];

        if (transition->isDetectingTransition())
        {
          detection_assay_it = i;
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
          auto new_end = std::remove_if(chromatogram.begin(), chromatogram.end(),
                                        [rt_min, rt_max](const ChromatogramPeak& chr)
                                        { return chr.getRT() > rt_max  || chr.getRT() < rt_min; });
          chromatogram.erase(new_end, chromatogram.end());
        }

        // Now add the transition and the chromatogram to the MRMTransitionGroup
        transition_group.addTransition(*transition, transition->getNativeID());
        transition_group.addChromatogram(chromatogram, chromatogram.getNativeID());
      }

      // currently .tsv, .osw and .featureXML are mutually exclusive
      if (tsv_writer.isActive() || osw_writer.isActive()) { output.clear(); }

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
      if (detection_assay_it < 0 && output.size() > 0)
      {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Error, did not find any detection transition for feature " + id );
      }

      // 5. Add to the output tsv if given
      if (tsv_writer.isActive() && output.size() > 0) // implies that detection_assay_it was set
      {
        const OpenSwath::LightCompound pep = transition_exp.getCompounds()[ assay_peptide_map[id] ];
        const TransitionType* transition = assay_it->second[detection_assay_it];
        to_tsv_output.push_back(tsv_writer.prepareLine(pep, transition, output, id));
      }

      // 6. Add to the output osw if given
      if (osw_writer.isActive() && output.size() > 0) // implies that detection_assay_it was set
      {
        const OpenSwath::LightCompound pep = transition_exp.getCompounds()[ assay_peptide_map[id] ];
        const TransitionType* transition = assay_it->second[detection_assay_it];
        to_osw_output.push_back(osw_writer.prepareLine(pep, transition, output, id));
      }
    }

    // Only write at the very end since this is a step that needs a barrier
    if (tsv_writer.isActive())
    {
#ifdef _OPENMP
#pragma omp critical (osw_write_tsv)
#endif
      {
        tsv_writer.writeLines(to_tsv_output);
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
                                                            const TransformationDescription trafo_inverse, 
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

// OpenSwathWorkflowSonar
namespace OpenMS
{

    void OpenSwathWorkflowSonar::performExtractionSonar(
           const std::vector< OpenSwath::SwathMap > & swath_maps,
           const TransformationDescription trafo,
           const ChromExtractParams & cp,
           const ChromExtractParams & cp_ms1,
           const Param & feature_finder_param,
           const OpenSwath::LightTargetedExperiment& transition_exp,
           FeatureMap& out_featureFile,
           bool store_features,
           OpenSwathTSVWriter & tsv_writer,
           OpenSwathOSWWriter & osw_writer,
           Interfaces::IMSDataConsumer * chromConsumer,
           int batchSize,
           bool load_into_memory)
    {
      tsv_writer.writeHeader();
      osw_writer.writeHeader();

      // Compute inversion of the transformation
      TransformationDescription trafo_inverse = trafo;
      trafo_inverse.invert();

      if (swath_maps.empty() )
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          String("No swath maps provided"));
      }

      if (use_ms1_traces_) ms1_map_ = loadMS1Map(swath_maps, load_into_memory);

      // (i) Obtain precursor chromatograms (MS1) if precursor extraction is enabled
      std::vector< MSChromatogram > ms1_chromatograms;
      if (ms1_map_ != nullptr)
      {
        MS1Extraction_(ms1_map_, swath_maps, ms1_chromatograms, chromConsumer, cp_ms1,
            transition_exp, trafo_inverse);
      }

      ///////////////////////////////////////////////////////////////////////////
      // (ii) Compute SONAR window sizes and upper/lower limit
      double sonar_winsize, sonar_start, sonar_end;
      int sonar_total_win;
      computeSonarWindows_(swath_maps, sonar_winsize, sonar_start, sonar_end, sonar_total_win);

      std::cout << "Will analyze " << transition_exp.transitions.size() << " transitions in total." << std::endl;
      int progress = 0;
      this->startProgress(0, sonar_total_win, "Extracting and scoring transitions");

      ///////////////////////////////////////////////////////////////////////////
      // Iterate through all SONAR windows
      // We set dynamic scheduling such that the SONAR windows are worked on in
      // the order in which they were given to the program / acquired. This
      // gives much better load balancing than static allocation.
      // TODO: this means that there is possibly some overlap between threads accessing sptr ... !!
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
      for (int sonar_idx = 0; sonar_idx < sonar_total_win; sonar_idx++)
      {
        double currwin_start = sonar_start + sonar_idx * sonar_winsize;
        double currwin_end = currwin_start + sonar_winsize;
        OPENMS_LOG_DEBUG << "   ====  sonar window " << sonar_idx << " from " << currwin_start << " to " << currwin_end << std::endl;

        // Step 1: select which transitions to extract with the current windows (proceed in batches)
        OpenSwath::LightTargetedExperiment transition_exp_used_all;
        OpenSwathHelper::selectSwathTransitions(transition_exp, transition_exp_used_all,
            0, currwin_start, currwin_end);

        if (transition_exp_used_all.getTransitions().size() > 0) // skip if no transitions found
        {


          ////////////////////////////////// 
          // Identify which SONAR windows to use for current set of transitions
          ////////////////////////////////// 
          std::vector< OpenSwath::SwathMap > used_maps;
          for (size_t i = 0; i < swath_maps.size(); ++i)
          {
            if (swath_maps[i].ms1) {continue;} // skip MS1
              
            // TODO: what if the swath map is smaller than the current window ??
            if (  (currwin_start >= swath_maps[i].lower && currwin_start <= swath_maps[i].upper  ) ||
                  (currwin_end >= swath_maps[i].lower && currwin_end <= swath_maps[i].upper  ) )
            {
#ifdef OPENSWATH_WORKFLOW_DEBUG
              std::cout << " will use curr window  " << i << " : " << swath_maps[i].lower << "-" <<
                                                                      swath_maps[i].upper << std::endl;
#endif
              used_maps.push_back(swath_maps[i]);
            }
          }

          ////////////////////////////////// 
          // Threadsafe loading of identified maps
          ////////////////////////////////// 
          for (Size i = 0; i < used_maps.size(); i++)
          {
#ifdef _OPENMP
#pragma omp critical (loadMemory)
#endif
            {
              // Loading the maps is not threadsafe if they overlap (e.g.
              // multiple threads could access the same maps) which often
              // happens in SONAR. Thus we either create a threadsafe light
              // clone or load them into memory if requested.
              if (load_into_memory)
              {
                used_maps[i].sptr = boost::shared_ptr<SpectrumAccessOpenMSInMemory>( new SpectrumAccessOpenMSInMemory(*used_maps[i].sptr) );
              }
              else
              {
                used_maps[i].sptr = used_maps[i].sptr->lightClone();
              }
            }
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

#ifdef _OPENMP
#pragma omp critical (osw_write_stdout)
#endif
          {
            std::cout << "Thread " <<
#ifdef _OPENMP
            omp_get_thread_num() << " " <<
#endif
            "will analyze " << transition_exp_used_all.getCompounds().size() <<  " compounds and "
            << transition_exp_used_all.getTransitions().size() <<  " transitions "
            "from SONAR SWATH " << sonar_idx << " in batches of " << batch_size << std::endl;
          }
          for (size_t pep_idx = 0; pep_idx <= (transition_exp_used_all.getCompounds().size() / batch_size); pep_idx++)
          {
            // Create the new, batch-size transition experiment
            OpenSwath::LightTargetedExperiment transition_exp_used;
            selectCompoundsForBatch_(transition_exp_used_all, transition_exp_used, batch_size, pep_idx);

            // Step 2.1: extract these transitions
            std::vector< OpenSwath::ChromatogramPtr > chrom_list;
            std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;

            // Step 2.2: prepare the extraction coordinates and extract chromatograms
            prepareExtractionCoordinates_(chrom_list, coordinates, transition_exp_used, trafo_inverse, cp);
            performSonarExtraction_(used_maps, coordinates, chrom_list, cp);

            // Step 2.3: convert chromatograms back to OpenMS::MSChromatogram and write to output
            PeakMap chrom_exp;
            ChromatogramExtractor().return_chromatogram(chrom_list, coordinates, transition_exp_used, SpectrumSettings(),
                                                        chrom_exp.getChromatograms(), false, cp.im_extraction_window);

            // Step 3: score these extracted transitions
            FeatureMap featureFile;
            scoreAllChromatograms_(chrom_exp.getChromatograms(), ms1_chromatograms, used_maps, transition_exp_used,
                                   feature_finder_param, trafo, cp.rt_extraction_window, featureFile, tsv_writer, osw_writer);

            // Step 4: write all chromatograms and features out into an output object / file
            // (this needs to be done in a critical section since we only have one
            // output file and one output map).
#ifdef _OPENMP
#pragma omp critical (osw_write_out)
#endif
            {
              writeOutFeaturesAndChroms_(chrom_exp.getChromatograms(), featureFile, out_featureFile, store_features, chromConsumer);
            }
          }
        }
#ifdef _OPENMP
#pragma omp critical (progress)
#endif
        this->setProgress(++progress);
      }
      this->endProgress();
    }


    void OpenSwathWorkflowSonar::computeSonarWindows_(const std::vector< OpenSwath::SwathMap > & swath_maps,
                                                      double & sonar_winsize,
                                                      double & sonar_start,
                                                      double & sonar_end,
                                                      int & sonar_total_win)
    {
      sonar_winsize = -1;
      sonar_start = std::numeric_limits<double>::max();
      sonar_end = -1;
      for (size_t i = 0; i < swath_maps.size(); ++i)
      {
        if (swath_maps[i].ms1) {continue;} // skip MS1

        // compute sonar window size (estimate)
        if (swath_maps[i].upper - swath_maps[i].lower > sonar_winsize)
        {
          sonar_winsize = swath_maps[i].upper - swath_maps[i].lower;
        }

        // compute start of SONAR range
        if (swath_maps[i].lower < sonar_start)
        {
          sonar_start = swath_maps[i].lower;
        }

        // compute end of SONAR range
        if (swath_maps[i].upper > sonar_end)
        {
          sonar_end = swath_maps[i].upper;
        }
      }

      // compute total number of windows
      sonar_total_win = int((sonar_end - sonar_start) / sonar_winsize) + 1;

#ifdef OPENSWATH_WORKFLOW_DEBUG
      std::cout << " will use  a total of " << sonar_total_win << " windows " << std::endl;
      for (int kk = 0; kk < sonar_total_win; kk++)
      {
        std::cout << " sonar window " << kk << " from " <<
          sonar_start + kk * sonar_winsize << " to " <<
          sonar_start + (kk+1) * sonar_winsize << std::endl;
      }
#endif

    }


    void OpenSwathWorkflowSonar::performSonarExtraction_(const std::vector< OpenSwath::SwathMap > & used_maps,
                                 const std::vector< ChromatogramExtractor::ExtractionCoordinates > & coordinates,
                                 std::vector< OpenSwath::ChromatogramPtr > & chrom_list,
                                 const ChromExtractParams & cp)
    {
      typedef std::vector< OpenSwath::ChromatogramPtr > chromatogramList;
      typedef std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinatesList;

      ChromatogramExtractor extractor;
      // Iterate over all SONAR maps we currently have and extract chromatograms from them
      for (size_t map_idx = 0; map_idx < used_maps.size(); map_idx++)
      {
        chromatogramList tmp_chromatogram_list;
        coordinatesList coordinates_used;

        for (size_t c_idx = 0; c_idx < coordinates.size(); c_idx++)
        {
          if (coordinates[c_idx].mz_precursor > used_maps[map_idx].lower &&
              coordinates[c_idx].mz_precursor < used_maps[map_idx].upper)
          {
            coordinates_used.push_back( coordinates[c_idx] );
            OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
            tmp_chromatogram_list.push_back(s);
          }
        }

#ifdef OPENSWATH_WORKFLOW_DEBUG
        std::cout << " in used maps, extract " << coordinates_used.size()
          << " coordinates from " << used_maps[map_idx].lower << "-" << used_maps[map_idx].upper << std::endl;
#endif

        extractor.extractChromatograms(used_maps[map_idx].sptr,
            tmp_chromatogram_list, coordinates_used,
            cp.mz_extraction_window, cp.ppm, cp.im_extraction_window, cp.extraction_function);

        // In order to reach maximal sensitivity and identify peaks in
        // the data, we will aggregate the data by adding all
        // chromatograms from different SONAR scans up
        size_t chrom_idx = 0;
        for (size_t c_idx = 0; c_idx < coordinates.size(); c_idx++)
        {
          if (coordinates[c_idx].mz_precursor > used_maps[map_idx].lower &&
              coordinates[c_idx].mz_precursor < used_maps[map_idx].upper)
          {

            OpenSwath::ChromatogramPtr s = tmp_chromatogram_list[chrom_idx];
            OpenSwath::ChromatogramPtr base_chrom = chrom_list[c_idx];

            /// add the new chromatogram to the one that we already have (the base chromatogram)
            chrom_list[c_idx] = addChromatograms(chrom_list[c_idx], tmp_chromatogram_list[chrom_idx]);

            chrom_idx++;
          }
        }
      }

#ifdef OPENSWATH_WORKFLOW_DEBUG
            // debug output ...
            std::cout << " done with extraction of all coordinates!!!" << std::endl;
            for (size_t c_idx = 0; c_idx < coordinates.size(); c_idx++)
            {
              {
                OpenSwath::ChromatogramPtr base_chrom = chrom_list[c_idx];

                std::cout << " coordinate  : " << coordinates[c_idx].id << " (" << coordinates[c_idx].mz << ")"<< std::endl;
                for (size_t kk = 0; kk < base_chrom->getIntensityArray()->data.size(); kk++)
                {
                  std::cout << " base chrom: " <<
                      base_chrom->getTimeArray()->data[kk] << " / "   <<
                      base_chrom->getIntensityArray()->data[kk] << std::endl;
                }
              }
            }
#endif


    }

    OpenSwath::ChromatogramPtr OpenSwathWorkflowSonar::addChromatograms(OpenSwath::ChromatogramPtr base_chrom, OpenSwath::ChromatogramPtr newchrom)
    {
      if (base_chrom->getTimeArray()->data.empty())
      {
        return newchrom;
      }

      LinearResamplerAlign ls;
      ls.raster(newchrom->getTimeArray()->data.begin(),
                newchrom->getTimeArray()->data.end(),
                newchrom->getIntensityArray()->data.begin(),
                newchrom->getIntensityArray()->data.end(),
                base_chrom->getTimeArray()->data.begin(),
                base_chrom->getTimeArray()->data.end(),
                base_chrom->getIntensityArray()->data.begin(),
                base_chrom->getIntensityArray()->data.end()
      );

      return base_chrom;
    }

}

