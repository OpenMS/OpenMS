// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/ANALYSIS/OPENSWATH/CalibrationWorkflow.h>
#include <OpenMS/ANALYSIS/TARGETED/IChromatogramHandler.h>
#include <OpenMS/ANALYSIS/TARGETED/ChromatogramProcessor.h>
#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <cmath>
#include <unordered_map>


// OpenSwathWorkflow
namespace OpenMS
{
  // Helper: load MS1 map (returns first swath map marked as ms1)
  OpenSwath::SpectrumAccessPtr OpenSwathWorkflow::loadMS1Map(const std::vector< OpenSwath::SwathMap > & swath_maps, bool load_into_memory)
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
    if (load_into_memory && ms1_map)
    {
      // This creates an InMemory object that keeps all data in memory
      // but provides the same access functionality to the raw data as
      // any object implementing ISpectrumAccess
      ms1_map = std::shared_ptr<SpectrumAccessOpenMSInMemory>( new SpectrumAccessOpenMSInMemory(*ms1_map) );
    }
    return ms1_map;
  }

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
    bool load_into_memory,
    const Param & mrm_mapping_param,
    MobilogramParquetConsumer * mobilogram_consumer)
  {
    bool ms1_only = (swath_maps.size() == 1 && swath_maps[0].ms1);

    if (mrm_)
    {
      this->startProgress(0, 1, "Extraction and Scoring");
      std::unique_ptr<IChromatogramHandler> provider = IChromatogramHandler::createDefault();
      std::vector<MSChromatogram> filtered_chroms = provider->extractAndMapChromatogramsForTransitions(swath_maps, transition_exp, cp, mrm_mapping_param);

      FeatureMap featureFile;
      scoreAllChromatograms_(filtered_chroms, std::vector<MSChromatogram>(), swath_maps, transition_exp,
                feature_finder_param, trafo, cp.rt_extraction_window, featureFile, osw_writer, ms1_isotopes, false, mobilogram_consumer);

      std::vector<MSChromatogram> empty_ms1_chromatograms;

      writeOutFeaturesAndChroms_(filtered_chroms, empty_ms1_chromatograms, featureFile, out_featureFile, store_features, chromConsumer);
      this->endProgress();
      return;
    }

    // Compute inversion of the transformation
    TransformationDescription trafo_inverse = trafo;
    trafo_inverse.invert();

    OPENMS_LOG_INFO << "Will analyze " << transition_exp.transitions.size() << " transitions in total." << std::endl;
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
      std::shared_ptr<MSExperiment> empty_exp = std::shared_ptr<MSExperiment>(new MSExperiment);

      const OpenSwath::LightTargetedExperiment& transition_exp_used = transition_exp;
      scoreAllChromatograms_(std::vector<MSChromatogram>(), ms1_chromatograms, swath_maps, transition_exp_used,
                feature_finder_param, trafo,
                cp.rt_extraction_window, featureFile, osw_writer, ms1_isotopes, true, mobilogram_consumer);

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
      OPENMS_LOG_ERROR << "Setting -pasef and -matching_window_only flags simultaneously is not currently supported." << std::endl;
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
                imOld << " with swath map of im center " << imNew << '\n';
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
  OPENMS_LOG_INFO << "Setting up nested loop with " << std::min(threads_outer_loop_, omp_get_max_threads()) << " threads out of "<< omp_get_max_threads() << std::endl;
      omp_set_nested(1);
      omp_set_dynamic(0);
      omp_set_num_threads(std::min(threads_outer_loop_, omp_get_max_threads()) ); // use at most threads_outer_loop_ threads here
    }
    else
    {
  OPENMS_LOG_INFO << "Use non-nested loop with " << total_nr_threads << " threads." << std::endl;
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
               OPENMS_LOG_DEBUG << "Adding Precursor with m/z " << tr.getPrecursorMZ() << " and IM of " << tr.getPrecursorIM() <<  " to swath with mz upper of " << swath_maps[i].upper << " im lower of " << swath_maps[i].imLower << " and im upper of " << swath_maps[i].imUpper << '\n';
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
            current_swath_map = std::shared_ptr<SpectrumAccessOpenMSInMemory>( new SpectrumAccessOpenMSInMemory(*current_swath_map) );
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

          const Size n_compounds = transition_exp_used_all.getCompounds().size();
          SignedSize nr_batches = 0;
          if (batch_size > 0)
          {
            nr_batches = static_cast<SignedSize>((n_compounds + batch_size - 1) / batch_size);
          }

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
          for (SignedSize pep_idx = 0; pep_idx < nr_batches; ++pep_idx)
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
              OPENMS_LOG_INFO << "Thread " <<
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
              feature_finder_param, trafo, cp.rt_extraction_window, featureFile, osw_writer, ms1_isotopes, false, mobilogram_consumer);

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
    bool ms1only,
    MobilogramParquetConsumer* mobilogram_consumer) const
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
      OPENMS_LOG_WARN << "WARNING: Attempted to use MS1 traces but no MS1 map was provided: Will not use MS1 signal!\n";
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

    std::unordered_map<String, int> ms1_chromatogram_map;
    ms1_chromatogram_map.reserve(ms1_chromatograms.size());
    for (Size i = 0; i < ms1_chromatograms.size(); i++)
    {
      ms1_chromatogram_map[ms1_chromatograms[i].getNativeID()] = boost::numeric_cast<int>(i);
    }

    std::unordered_map<String, int> chromatogram_map;
    chromatogram_map.reserve(ms2_chromatograms.size());
    for (Size i = 0; i < ms2_chromatograms.size(); i++)
    {
      const String cid = ms2_chromatograms[i].getNativeID();
      chromatogram_map[cid] = boost::numeric_cast<int>(i);
    }
    std::unordered_map<String, int> assay_peptide_map;
    assay_peptide_map.reserve(transition_exp.getCompounds().size());
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
          if (mrm_) // in SRM/MRM mode, we can skip missing chromatograms, because it's unlikely that we will map all transitions to the already targeted extracted chromatograms
          {
            OPENMS_LOG_DEBUG << "Did not find chromatogram for transition " << transition->getNativeID()
                             << "; skipping this transition." << std::endl;
            continue;
          }
          else
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "Error, did not find chromatogram for transition " + transition->getNativeID() );
          }
        }

        // Convert chromatogram to MSChromatogram and filter (properly handles DataArrays via select())
        auto chromatogram = ms2_chromatograms[ chromatogram_map[transition->getNativeID()] ];
        chromatogram.setNativeID(transition->getNativeID());
        if (rt_extraction_window > 0)
        {
          double de_normalized_experimental_rt = trafo_inv.apply(expected_rt);
          double rt_max = de_normalized_experimental_rt + rt_extraction_window;
          double rt_min = de_normalized_experimental_rt - rt_extraction_window;
          std::vector<Size> indices_to_keep;
          indices_to_keep.reserve(chromatogram.size());
          for (Size i = 0; i < chromatogram.size(); ++i)
          {
            double rt = chromatogram[i].getRT();
            if (rt >= rt_min && rt <= rt_max)
            {
              indices_to_keep.push_back(i);
            }
          }
          if (indices_to_keep.size() != chromatogram.size())
          {
            chromatogram.select(indices_to_keep);
          }
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
          const MSChromatogram& chromatogram = ms1_chromatograms[ ms1_chromatogram_map[prec_id] ];
          transition_group.addPrecursorChromatogram(chromatogram, chromatogram.getNativeID());
        }
      }

      // 3. / 4. Process the MRMTransitionGroup: find peakgroups and score them
      // For SRM/MRM, If there are no chromatograms added to this transition_group (e.g. all transitions were unmapped/skipped), skip processing to avoid range errors in the picker/feature finder.
      if (transition_group.getChromatograms().empty() && transition_group.getPrecursorChromatograms().empty())
      {
        OPENMS_LOG_DEBUG << "No chromatograms present for assay " << id << "; skipping scoring." << std::endl;
        continue;
      }

      // Ensure there is at least one chromatogram originating from a
      // detecting transition. MRMTransitionGroupPicker expects at least one
      // detecting chromatogram.
      bool has_detecting_chrom = false;
      // Check fragment/transition chromatograms first
      for (const auto & chrom : transition_group.getChromatograms())
      {
        String nid = chrom.getNativeID();
        if (transition_group.getTransitions().size() == 0)
        {
          has_detecting_chrom = true; break;
        }
        if (transition_group.hasTransition(nid) && transition_group.getTransition(nid).isDetectingTransition())
        {
          has_detecting_chrom = true; break;
        }
      }
      // If no fragment chromatograms were found, allow precursor (MS1) chromatograms
      if (!has_detecting_chrom && !transition_group.getPrecursorChromatograms().empty())
      {
        has_detecting_chrom = true;
      }
      if (!has_detecting_chrom)
      {
        if (mrm_)
        {
          OPENMS_LOG_DEBUG << "No detecting chromatograms for assay " << id
                           << "; skipping this assay." << std::endl;
          continue;
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Error, did not find any detecting chromatogram for assay " + id );
        }
      }

      try
      {
        trgroup_picker.pickTransitionGroup(transition_group);
      }
      catch (const Exception::InvalidRange & e)
      {
        // In SRM/MRM mode extraction/mapping failures are common for noisy or
        // partially-mapped chromatograms -> log and continue. In DIA mode an
        // InvalidRange usually indicates a real extraction/mapping bug and
        // should surface during CI, so re-throw to fail fast.
        if (!mrm_)
        {
          throw;
        }
        OPENMS_LOG_ERROR << "InvalidRange while picking transition group for assay " << id
                         << " - transitions=" << transition_group.getTransitions().size()
                         << ", chroms=" << transition_group.getChromatograms().size()
                         << ", prec_chroms=" << transition_group.getPrecursorChromatograms().size()
                         << ": " << e.getMessage() << " - skipping assay." << std::endl;
        continue;
      }
      catch (const Exception::IllegalArgument & e)
      {
        // IllegalArgument may be thrown when mapping fails (e.g. no matching
        // chromatograms). Treat mode-specifically similar to InvalidRange.
        if (!mrm_)
        {
          throw;
        }
        OPENMS_LOG_ERROR << "IllegalArgument while picking transition group for assay " << id
                         << ": " << e.getMessage() << " - skipping assay." << std::endl;
        continue;
      }
      catch (const std::exception & e)
      {
        OPENMS_LOG_ERROR << "Exception while picking transition group for assay " << id
                         << ": " << e.what() << " - skipping assay." << std::endl;
        continue;
      }

      try
      {
        featureFinder.scorePeakgroups(transition_group, trafo, swath_maps, output, ms1only, mobilogram_consumer);
      }
      catch (const Exception::InvalidRange & e)
      {
        if (!mrm_)
        {
          throw;
        }
        OPENMS_LOG_ERROR << "InvalidRange while scoring transition group for assay " << id
                         << " - transitions=" << transition_group.getTransitions().size()
                         << ", chroms=" << transition_group.getChromatograms().size()
                         << ", prec_chroms=" << transition_group.getPrecursorChromatograms().size()
                         << ": " << e.getMessage() << " - skipping assay." << std::endl;
        continue;
      }
      catch (const Exception::IllegalArgument & e)
      {
        if (!mrm_)
        {
          throw;
        }
        OPENMS_LOG_ERROR << "IllegalArgument while scoring transition group for assay " << id
                         << ": " << e.getMessage() << " - skipping assay." << std::endl;
        continue;
      }
      catch (const std::exception & e)
      {
        OPENMS_LOG_ERROR << "Exception while scoring transition group for assay " << id
                         << ": " << e.what() << " - skipping assay." << std::endl;
        continue;
      }

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
