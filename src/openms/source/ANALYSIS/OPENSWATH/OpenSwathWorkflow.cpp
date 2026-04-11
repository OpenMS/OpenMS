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
#include <OpenMS/SYSTEM/StopWatch.h>
#include <chrono>
#include <cmath>
#include <exception>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <utility>


// OpenSwathWorkflow
namespace OpenMS
{
  namespace
  {
    using ProfileClock = std::chrono::steady_clock;

    double elapsedProfileSeconds(const ProfileClock::time_point& start)
    {
      return std::chrono::duration<double>(ProfileClock::now() - start).count();
    }

    ProfileClock::time_point profileStart(bool enabled)
    {
      return enabled ? ProfileClock::now() : ProfileClock::time_point();
    }

    struct OpenSwathWorkflowPhaseTiming
    {
      double ms1_map_loading = 0.0;
      double window_mapping = 0.0;
      double transition_selection = 0.0;
      double swath_loading = 0.0;
      double batch_selection = 0.0;
      double ms1_extraction = 0.0;
      double coordinate_preparation = 0.0;
      double ms2_extraction = 0.0;
      double chromatogram_conversion = 0.0;
      double scoring = 0.0;
      double writing = 0.0;
      double total = 0.0;
      OpenSwathScoringPhaseTiming scoring_breakdown;
      Size swath_count = 0;
      Size batch_count = 0;
      Size compound_count = 0;
      Size transition_count = 0;
      Size ms1_chromatogram_count = 0;
      Size ms2_chromatogram_count = 0;
      Size feature_count = 0;

      void add(const OpenSwathWorkflowPhaseTiming& rhs)
      {
        ms1_map_loading += rhs.ms1_map_loading;
        window_mapping += rhs.window_mapping;
        transition_selection += rhs.transition_selection;
        swath_loading += rhs.swath_loading;
        batch_selection += rhs.batch_selection;
        ms1_extraction += rhs.ms1_extraction;
        coordinate_preparation += rhs.coordinate_preparation;
        ms2_extraction += rhs.ms2_extraction;
        chromatogram_conversion += rhs.chromatogram_conversion;
        scoring += rhs.scoring;
        writing += rhs.writing;
        total += rhs.total;
        scoring_breakdown.add(rhs.scoring_breakdown);
        swath_count += rhs.swath_count;
        batch_count += rhs.batch_count;
        compound_count += rhs.compound_count;
        transition_count += rhs.transition_count;
        ms1_chromatogram_count += rhs.ms1_chromatogram_count;
        ms2_chromatogram_count += rhs.ms2_chromatogram_count;
        feature_count += rhs.feature_count;
      }
    };

    String formatProfileSeconds(double seconds)
    {
      return StopWatch::toString(seconds);
    }

    String formatScoringBreakdown(const OpenSwathScoringPhaseTiming& timing)
    {
      std::ostringstream os;
      os << "score_detail=(map_setup=" << formatProfileSeconds(timing.map_setup)
         << ", assay_setup=" << formatProfileSeconds(timing.assay_setup)
         << ", picker=" << formatProfileSeconds(timing.transition_group_picker)
         << ", score_peakgroups=" << formatProfileSeconds(timing.score_peakgroups)
         << ", score_setup=" << formatProfileSeconds(timing.score_setup)
         << ", sn_setup=" << formatProfileSeconds(timing.signal_to_noise_setup)
         << ", chrom=" << formatProfileSeconds(timing.chrom_scores)
         << ", library=" << formatProfileSeconds(timing.library_scores)
         << ", dia=" << formatProfileSeconds(timing.dia_scores)
         << ", ms1=" << formatProfileSeconds(timing.ms1_scores)
         << ", uis=" << formatProfileSeconds(timing.uis_scores)
         << ", score_output=" << formatProfileSeconds(timing.score_output)
         << ", feature_sort_output=" << formatProfileSeconds(timing.feature_sort_output)
         << ", osw_prepare=" << formatProfileSeconds(timing.osw_prepare)
         << ", osw_write=" << formatProfileSeconds(timing.osw_write)
         << "; assays=" << timing.assay_count
         << ", scored_assays=" << timing.scored_assay_count
         << ", skipped_assays=" << timing.skipped_assay_count
         << ", scored_features=" << timing.scored_feature_count
         << ")";
      return String(os.str());
    }

    void logBatchProfile(SignedSize swath_index, SignedSize batch_index, SignedSize nr_batches, const OpenSwathWorkflowPhaseTiming& timing)
    {
      OPENMS_LOG_INFO << "[OpenSwathWorkflow profile] SWATH " << swath_index
                      << " batch " << batch_index << "/" << nr_batches
                      << ": total=" << formatProfileSeconds(timing.total)
                      << " (batch_select=" << formatProfileSeconds(timing.batch_selection)
                      << ", ms1_extract=" << formatProfileSeconds(timing.ms1_extraction)
                      << ", prepare=" << formatProfileSeconds(timing.coordinate_preparation)
                      << ", ms2_extract=" << formatProfileSeconds(timing.ms2_extraction)
                      << ", convert=" << formatProfileSeconds(timing.chromatogram_conversion)
                      << ", score=" << formatProfileSeconds(timing.scoring)
                      << ", write=" << formatProfileSeconds(timing.writing)
                      << "); compounds=" << timing.compound_count
                      << ", transitions=" << timing.transition_count
                      << ", ms1_chromatograms=" << timing.ms1_chromatogram_count
                      << ", ms2_chromatograms=" << timing.ms2_chromatogram_count
                      << ", features=" << timing.feature_count
                      << "; " << formatScoringBreakdown(timing.scoring_breakdown) << std::endl;
    }

    void logSwathProfile(SignedSize swath_index, const OpenSwathWorkflowPhaseTiming& timing)
    {
      OPENMS_LOG_INFO << "[OpenSwathWorkflow profile] SWATH " << swath_index
                      << " total=" << formatProfileSeconds(timing.total)
                      << " (transition_select=" << formatProfileSeconds(timing.transition_selection)
                      << ", swath_load=" << formatProfileSeconds(timing.swath_loading)
                      << ", batch_select=" << formatProfileSeconds(timing.batch_selection)
                      << ", ms1_extract=" << formatProfileSeconds(timing.ms1_extraction)
                      << ", prepare=" << formatProfileSeconds(timing.coordinate_preparation)
                      << ", ms2_extract=" << formatProfileSeconds(timing.ms2_extraction)
                      << ", convert=" << formatProfileSeconds(timing.chromatogram_conversion)
                      << ", score=" << formatProfileSeconds(timing.scoring)
                      << ", write=" << formatProfileSeconds(timing.writing)
                      << "); batches=" << timing.batch_count
                      << ", compounds=" << timing.compound_count
                      << ", transitions=" << timing.transition_count
                      << "; " << formatScoringBreakdown(timing.scoring_breakdown) << std::endl;
    }

    void logWorkflowProfileSummary(double workflow_wall_time, const OpenSwathWorkflowPhaseTiming& timing)
    {
      OPENMS_LOG_INFO << "[OpenSwathWorkflow profile] Workflow wall="
                      << formatProfileSeconds(workflow_wall_time)
                      << "; summed worker phases="
                      << "ms1_map_load=" << formatProfileSeconds(timing.ms1_map_loading)
                      << ", window_map=" << formatProfileSeconds(timing.window_mapping)
                      << ", transition_select=" << formatProfileSeconds(timing.transition_selection)
                      << ", swath_load=" << formatProfileSeconds(timing.swath_loading)
                      << ", batch_select=" << formatProfileSeconds(timing.batch_selection)
                      << ", ms1_extract=" << formatProfileSeconds(timing.ms1_extraction)
                      << ", prepare=" << formatProfileSeconds(timing.coordinate_preparation)
                      << ", ms2_extract=" << formatProfileSeconds(timing.ms2_extraction)
                      << ", convert=" << formatProfileSeconds(timing.chromatogram_conversion)
                      << ", score=" << formatProfileSeconds(timing.scoring)
                      << ", write=" << formatProfileSeconds(timing.writing)
                      << ", swath_total=" << formatProfileSeconds(timing.total)
                      << "; swaths=" << timing.swath_count
                      << ", batches=" << timing.batch_count
                      << ", compounds=" << timing.compound_count
                      << ", transitions=" << timing.transition_count
                      << ", ms1_chromatograms=" << timing.ms1_chromatogram_count
                      << ", ms2_chromatograms=" << timing.ms2_chromatogram_count
                      << ", features=" << timing.feature_count
                      << "; " << formatScoringBreakdown(timing.scoring_breakdown) << std::endl;
    }
  }

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

    const bool profile = phase_profiling_enabled_;
    const auto profile_workflow_start = profileStart(profile);
    OpenSwathWorkflowPhaseTiming profile_total;

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

    ProfileClock::time_point profile_phase_start = profileStart(profile);
    if (use_ms1_traces_) ms1_map_ = loadMS1Map(swath_maps, load_into_memory);
    if (profile)
    {
      profile_total.ms1_map_loading = elapsedProfileSeconds(profile_phase_start);
    }

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
    profile_phase_start = profileStart(profile);
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
    if (profile)
    {
      profile_total.window_mapping = elapsedProfileSeconds(profile_phase_start);
    }

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
      OpenSwathWorkflowPhaseTiming swath_timing;
      const auto swath_profile_start = profileStart(profile);
      bool profile_swath = false;

      if (!swath_maps[i].ms1) // skip MS1
      {
        profile_swath = true;

        // Step 1: select which transitions to extract (proceed in batches)
        auto phase_start = profileStart(profile);
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
        if (profile)
        {
          swath_timing.transition_selection = elapsedProfileSeconds(phase_start);
          if (transition_exp_used_all.getTransitions().empty())
          {
            swath_timing.compound_count = transition_exp_used_all.getCompounds().size();
            swath_timing.transition_count = transition_exp_used_all.getTransitions().size();
          }
        }

        if (!transition_exp_used_all.getTransitions().empty()) // skip if no transitions found
        {

          phase_start = profileStart(profile);
          OpenSwath::SpectrumAccessPtr current_swath_map = swath_maps[i].sptr;
          if (load_into_memory)
          {
            // This creates an InMemory object that keeps all data in memory
            current_swath_map = std::shared_ptr<SpectrumAccessOpenMSInMemory>( new SpectrumAccessOpenMSInMemory(*current_swath_map) );
          }
          if (profile)
          {
            swath_timing.swath_loading = elapsedProfileSeconds(phase_start);
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
            OpenSwathWorkflowPhaseTiming batch_timing;
            const auto batch_profile_start = profileStart(profile);
            batch_timing.batch_count = 1;
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
            auto batch_phase_start = profileStart(profile);
            OpenSwath::LightTargetedExperiment transition_exp_used;
            selectCompoundsForBatch_(transition_exp_used_all, transition_exp_used, batch_size, pep_idx);
            if (profile)
            {
              batch_timing.batch_selection = elapsedProfileSeconds(batch_phase_start);
              batch_timing.compound_count = transition_exp_used.getCompounds().size();
              batch_timing.transition_count = transition_exp_used.getTransitions().size();
            }

            // Extract MS1 chromatograms for this batch
            batch_phase_start = profileStart(profile);
            std::vector< MSChromatogram > ms1_chromatograms;
            if (ms1_map_ != nullptr)
            {
              OpenSwath::SpectrumAccessPtr threadsafe_ms1 = ms1_map_->lightClone();
              MS1Extraction_(threadsafe_ms1, swath_maps, ms1_chromatograms, ms1_cp,
                  transition_exp_used, trafo_inverse, ms1_only, ms1_isotopes);
            }
            if (profile)
            {
              batch_timing.ms1_extraction = elapsedProfileSeconds(batch_phase_start);
              batch_timing.ms1_chromatogram_count = ms1_chromatograms.size();
            }

            // Step 2.1: extract these transitions
            ChromatogramExtractor extractor;
            std::vector< OpenSwath::ChromatogramPtr > chrom_list;
            std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;

            // Step 2.2: prepare the extraction coordinates and extract chromatograms
            // chrom_list contains one entry for each fragment ion (transition) in transition_exp_used
            batch_phase_start = profileStart(profile);
            prepareExtractionCoordinates_(chrom_list, coordinates, transition_exp_used, trafo_inverse, cp);
            if (profile)
            {
              batch_timing.coordinate_preparation = elapsedProfileSeconds(batch_phase_start);
            }

            batch_phase_start = profileStart(profile);
            extractor.extractChromatograms(current_swath_map_inner, chrom_list, coordinates, cp.mz_extraction_window,
                cp.ppm, cp.im_extraction_window, cp.extraction_function);
            if (profile)
            {
              batch_timing.ms2_extraction = elapsedProfileSeconds(batch_phase_start);
            }

            // Step 2.3: convert chromatograms back to OpenMS::MSChromatogram and write to output
            batch_phase_start = profileStart(profile);
            PeakMap chrom_exp;
            extractor.return_chromatogram(chrom_list, coordinates, transition_exp_used,  SpectrumSettings(),
                                          chrom_exp.getChromatograms(), false, cp.im_extraction_window);
            if (profile)
            {
              batch_timing.chromatogram_conversion = elapsedProfileSeconds(batch_phase_start);
              batch_timing.ms2_chromatogram_count = chrom_exp.getChromatograms().size();
            }


            // Step 3: score these extracted transitions
            batch_phase_start = profileStart(profile);
            FeatureMap featureFile;
            std::vector< OpenSwath::SwathMap > tmp = {swath_maps[i]};
            tmp.back().sptr = current_swath_map_inner;
            OpenSwathScoringPhaseTiming scoring_profile;
            scoreAllChromatograms_(chrom_exp.getChromatograms(), ms1_chromatograms, tmp, transition_exp_used,
              feature_finder_param, trafo, cp.rt_extraction_window, featureFile, osw_writer, ms1_isotopes, false,
              mobilogram_consumer, profile ? &scoring_profile : nullptr);
            if (profile)
            {
              batch_timing.scoring = elapsedProfileSeconds(batch_phase_start);
              batch_timing.feature_count = featureFile.size();
              batch_timing.scoring_breakdown.add(scoring_profile);
            }

            // Step 4: write all chromatograms and features out into an output object / file
            // (this needs to be done in a critical section since we only have one
            // output file and one output map).
            #pragma omp critical (osw_write_out)
            {
              batch_phase_start = profileStart(profile);
              writeOutFeaturesAndChroms_(chrom_exp.getChromatograms(), ms1_chromatograms, featureFile, out_featureFile, store_features, chromConsumer);
              if (profile)
              {
                batch_timing.writing = elapsedProfileSeconds(batch_phase_start);
              }
            }

            if (profile)
            {
              batch_timing.total = elapsedProfileSeconds(batch_profile_start);
#ifdef _OPENMP
#pragma omp critical (osw_profile_stdout)
#endif
              {
                swath_timing.add(batch_timing);
                logBatchProfile(i, pep_idx, nr_batches, batch_timing);
              }
            }
          }

        } // continue 2 (no continue due to OpenMP)
      } // continue 1 (no continue due to OpenMP)

      if (profile && profile_swath)
      {
        swath_timing.total = elapsedProfileSeconds(swath_profile_start);
        swath_timing.swath_count = 1;
#ifdef _OPENMP
#pragma omp critical (osw_profile_stdout)
#endif
        {
          logSwathProfile(i, swath_timing);
          profile_total.add(swath_timing);
        }
      }

      #pragma omp critical (progress)
      this->setProgress(++progress);

    }
    this->endProgress();

    if (profile)
    {
      logWorkflowProfileSummary(elapsedProfileSeconds(profile_workflow_start), profile_total);
    }

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
    MobilogramParquetConsumer* mobilogram_consumer,
    OpenSwathScoringPhaseTiming* scoring_profile) const
  {
    const bool profile = scoring_profile != nullptr;
    if (profile)
    {
      *scoring_profile = OpenSwathScoringPhaseTiming();
    }

    auto profile_phase_start = profileStart(profile);
    TransformationDescription trafo_inv = trafo;
    trafo_inv.invert();

    if (use_ms1_traces_ && !ms1_map_)
    {
      OPENMS_LOG_WARN << "WARNING: Attempted to use MS1 traces but no MS1 map was provided: Will not use MS1 signal!\n";
    }

    // If use_total_mi_score is defined, we need to instruct MRMTransitionGroupPicker to compute the score
    Param trgroup_picker_param = feature_finder_param.copy("TransitionGroupPicker:", true);
    if ((bool)feature_finder_param.getValue("Scores:use_total_mi_score").toBool())
    {
      trgroup_picker_param.setValue("compute_total_mi", "true");
    }

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

    struct AssayWorkItem
    {
      String id;
      std::vector<const TransitionType*> transitions;
    };
    std::vector<AssayWorkItem> assays;
    assays.reserve(assay_map.size());
    for (const auto& assay : assay_map)
    {
      assays.push_back({assay.first, assay.second});
    }

#ifdef _OPENMP
    const int worker_count = omp_in_parallel() ? omp_get_num_threads() : 1;
#else
    const int worker_count = 1;
#endif

    std::vector<std::unique_ptr<MRMFeatureFinderScoring>> feature_finders;
    feature_finders.reserve(worker_count);
    std::vector<std::unique_ptr<MRMTransitionGroupPicker>> trgroup_pickers;
    trgroup_pickers.reserve(worker_count);
    std::vector<std::vector<OpenSwath::SwathMap>> worker_swath_maps;
    worker_swath_maps.reserve(worker_count);
    std::vector<OpenSwathScoringPhaseTiming> worker_profiles(worker_count);
    for (int worker_idx = 0; worker_idx < worker_count; ++worker_idx)
    {
      std::vector<OpenSwath::SwathMap> cloned_swath_maps = swath_maps;
      if (worker_count > 1)
      {
        for (auto& swath_map : cloned_swath_maps)
        {
          if (swath_map.sptr)
          {
            swath_map.sptr = swath_map.sptr->lightClone();
          }
        }
      }
      worker_swath_maps.push_back(std::move(cloned_swath_maps));

      auto feature_finder = std::make_unique<MRMFeatureFinderScoring>();
      feature_finder->setParameters(feature_finder_param);
      feature_finder->setScoringProfiling(profile);
      feature_finder->resetScoringProfile();
      feature_finder->prepareProteinPeptideMaps_(transition_exp);

      // SpectrumAccess implementations are not generally thread-safe. Each
      // OpenMP worker that may score DIA/MS1 spectra gets its own light clone.
      if (use_ms1_traces_ && ms1_map_)
      {
        feature_finder->setMS1Map(ms1_map_->lightClone());
      }
      feature_finders.push_back(std::move(feature_finder));

      auto trgroup_picker = std::make_unique<MRMTransitionGroupPicker>();
      trgroup_picker->setParameters(trgroup_picker_param);
      trgroup_pickers.push_back(std::move(trgroup_picker));
    }
    if (profile)
    {
      scoring_profile->map_setup += elapsedProfileSeconds(profile_phase_start);
    }

    const bool osw_active = osw_writer.isActive();
    struct AssayResult
    {
      FeatureMap output;
      String osw_line;
    };
    std::vector<AssayResult> assay_results(assays.size());
    FeatureMap last_osw_output;
    SignedSize last_osw_output_index = -1;
    std::exception_ptr first_exception;

    auto recordException = [&](std::exception_ptr exception)
    {
#ifdef _OPENMP
#pragma omp critical (openswath_score_exception)
#endif
      {
        if (!first_exception)
        {
          first_exception = exception;
        }
      }
    };

    auto updateLastOSWOutput = [&](Size assay_idx, FeatureMap& assay_output)
    {
#ifdef _OPENMP
#pragma omp critical (openswath_last_osw_output)
#endif
      {
        if (static_cast<SignedSize>(assay_idx) > last_osw_output_index)
        {
          last_osw_output = assay_output;
          last_osw_output_index = static_cast<SignedSize>(assay_idx);
        }
      }
    };

    auto processAssay = [&](Size assay_idx)
    {
      try
      {
#ifdef _OPENMP
        const int worker_idx = omp_in_parallel() ? omp_get_thread_num() : 0;
#else
        const int worker_idx = 0;
#endif
        MRMFeatureFinderScoring& feature_finder = *feature_finders[worker_idx];
        MRMTransitionGroupPicker& trgroup_picker = *trgroup_pickers[worker_idx];
        const std::vector<OpenSwath::SwathMap>& thread_swath_maps = worker_swath_maps[worker_idx];
        OpenSwathScoringPhaseTiming& worker_profile = worker_profiles[worker_idx];
        const AssayWorkItem& assay = assays[assay_idx];

      auto assay_setup_start = profileStart(profile);
      bool assay_setup_recorded = false;
      auto recordAssaySetup = [&]()
      {
        if (profile && !assay_setup_recorded)
        {
          worker_profile.assay_setup += elapsedProfileSeconds(assay_setup_start);
          assay_setup_recorded = true;
        }
      };
      if (profile)
      {
        ++worker_profile.assay_count;
      }

      // Create new MRMTransitionGroup
      String id = assay.id;
      MRMTransitionGroupType transition_group;
      transition_group.setTransitionGroupID(id);
      double expected_rt = transition_exp.getCompounds()[ assay_peptide_map[id] ].rt;

      // 1. Go through all transitions, for each transition get
      // the chromatogram and the assay to the MRMTransitionGroup
      const TransitionType* detection_assay_it = nullptr; // store last detecting transition
      for (const TransitionType* transition : assay.transitions)
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
        recordAssaySetup();
        if (profile)
        {
          ++worker_profile.skipped_assay_count;
        }
        return;
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
          recordAssaySetup();
          if (profile)
          {
            ++worker_profile.skipped_assay_count;
          }
          return;
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Error, did not find any detecting chromatogram for assay " + id );
        }
      }

      recordAssaySetup();
      auto picker_profile_start = profileStart(profile);
      bool picker_profile_recorded = false;
      auto recordPickerProfile = [&]()
      {
        if (profile && !picker_profile_recorded)
        {
          worker_profile.transition_group_picker += elapsedProfileSeconds(picker_profile_start);
          picker_profile_recorded = true;
        }
      };
      try
      {
        trgroup_picker.pickTransitionGroup(transition_group);
        recordPickerProfile();
      }
      catch (const Exception::InvalidRange & e)
      {
        recordPickerProfile();
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
        if (profile)
        {
          ++worker_profile.skipped_assay_count;
        }
        return;
      }
      catch (const Exception::IllegalArgument & e)
      {
        recordPickerProfile();
        // IllegalArgument may be thrown when mapping fails (e.g. no matching
        // chromatograms). Treat mode-specifically similar to InvalidRange.
        if (!mrm_)
        {
          throw;
        }
        OPENMS_LOG_ERROR << "IllegalArgument while picking transition group for assay " << id
                         << ": " << e.getMessage() << " - skipping assay." << std::endl;
        if (profile)
        {
          ++worker_profile.skipped_assay_count;
        }
        return;
      }
      catch (const std::exception & e)
      {
        recordPickerProfile();
        OPENMS_LOG_ERROR << "Exception while picking transition group for assay " << id
                         << ": " << e.what() << " - skipping assay." << std::endl;
        if (profile)
        {
          ++worker_profile.skipped_assay_count;
        }
        return;
      }

      auto score_profile_start = profileStart(profile);
      bool score_profile_recorded = false;
      auto recordScoreProfile = [&]()
      {
        if (profile && !score_profile_recorded)
        {
          worker_profile.score_peakgroups += elapsedProfileSeconds(score_profile_start);
          score_profile_recorded = true;
        }
      };
      FeatureMap assay_output;
      try
      {
        feature_finder.scorePeakgroups(transition_group, trafo, thread_swath_maps, assay_output, ms1only, mobilogram_consumer);
        recordScoreProfile();
        if (profile)
        {
          ++worker_profile.scored_assay_count;
        }
      }
      catch (const Exception::InvalidRange & e)
      {
        recordScoreProfile();
        if (!mrm_)
        {
          throw;
        }
        OPENMS_LOG_ERROR << "InvalidRange while scoring transition group for assay " << id
                         << " - transitions=" << transition_group.getTransitions().size()
                         << ", chroms=" << transition_group.getChromatograms().size()
                         << ", prec_chroms=" << transition_group.getPrecursorChromatograms().size()
                         << ": " << e.getMessage() << " - skipping assay." << std::endl;
        if (profile)
        {
          ++worker_profile.skipped_assay_count;
        }
        return;
      }
      catch (const Exception::IllegalArgument & e)
      {
        recordScoreProfile();
        if (!mrm_)
        {
          throw;
        }
        OPENMS_LOG_ERROR << "IllegalArgument while scoring transition group for assay " << id
                         << ": " << e.getMessage() << " - skipping assay." << std::endl;
        if (profile)
        {
          ++worker_profile.skipped_assay_count;
        }
        return;
      }
      catch (const std::exception & e)
      {
        recordScoreProfile();
        OPENMS_LOG_ERROR << "Exception while scoring transition group for assay " << id
                         << ": " << e.what() << " - skipping assay." << std::endl;
        if (profile)
        {
          ++worker_profile.skipped_assay_count;
        }
        return;
      }

      // Ensure that a detection transition is used to derive features for output
      if (detection_assay_it == nullptr && !assay_output.empty())
      {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Error, did not find any detection transition for feature " + id );
      }

      // 5. Add to the output osw if given
      if (osw_active && !assay_output.empty()) // implies that detection_assay_it was set
      {
        auto osw_prepare_start = profileStart(profile);
        assay_results[assay_idx].osw_line = osw_writer.prepareLine(OpenSwath::LightCompound(), // not used currently: transition_exp.getCompounds()[ assay_peptide_map[id] ],
                                                                   nullptr, // not used currently: detection_assay_it,
                                                                   assay_output,
                                                                   id);
        updateLastOSWOutput(assay_idx, assay_output);
        if (profile)
        {
          worker_profile.osw_prepare += elapsedProfileSeconds(osw_prepare_start);
        }
      }
      else if (!osw_active)
      {
        assay_results[assay_idx].output = assay_output;
      }
      }
      catch (...)
      {
        recordException(std::current_exception());
      }
    };

    ///////////////////////////////////
    // Start of main function
    // Iterating over all the assays
    ///////////////////////////////////
#ifdef _OPENMP
    const bool use_assay_tasks = omp_in_parallel() && worker_count > 1 && mobilogram_consumer == nullptr && assays.size() >= 64;
    if (use_assay_tasks)
    {
      constexpr Size assay_task_grain = 8;
#pragma omp taskgroup
      {
        for (Size task_begin = 0; task_begin < assays.size(); task_begin += assay_task_grain)
        {
          const Size task_end = std::min(task_begin + assay_task_grain, assays.size());
#pragma omp task firstprivate(task_begin, task_end)
          {
            for (Size assay_idx = task_begin; assay_idx < task_end; ++assay_idx)
            {
              processAssay(assay_idx);
            }
          }
        }
      }
    }
    else
#endif
    {
      for (Size assay_idx = 0; assay_idx < assays.size(); ++assay_idx)
      {
        processAssay(assay_idx);
      }
    }

    if (first_exception)
    {
      std::rethrow_exception(first_exception);
    }

    std::vector<String> to_osw_output;
    if (osw_active)
    {
      output = last_osw_output;
      to_osw_output.reserve(assay_results.size());
      for (const AssayResult& result : assay_results)
      {
        if (!result.osw_line.empty())
        {
          to_osw_output.push_back(result.osw_line);
        }
      }
    }
    else
    {
      for (const AssayResult& result : assay_results)
      {
        for (const Feature& feature : result.output)
        {
          output.push_back(feature);
        }
        for (const ProteinIdentification& protein : result.output.getProteinIdentifications())
        {
          output.getProteinIdentifications().push_back(protein);
        }
      }
    }

    if (profile)
    {
      for (int worker_idx = 0; worker_idx < worker_count; ++worker_idx)
      {
        scoring_profile->add(worker_profiles[worker_idx]);
        scoring_profile->add(feature_finders[worker_idx]->getScoringProfile());
      }
    }

    // Only write at the very end since this is a step that needs a barrier
    if (osw_writer.isActive())
    {
#ifdef _OPENMP
#pragma omp critical (osw_write_tsv)
#endif
      {
        auto osw_write_start = profileStart(profile);
        osw_writer.writeLines(to_osw_output);
        if (profile)
        {
          scoring_profile->osw_write += elapsedProfileSeconds(osw_write_start);
        }
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
