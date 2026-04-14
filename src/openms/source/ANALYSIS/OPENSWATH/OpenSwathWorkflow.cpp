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
#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <exception>
#include <iterator>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <thread>
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

    constexpr UInt64 OSW_MIN_WRITE_BUFFER_BYTES = 64ull * 1024ull * 1024ull;
    constexpr UInt64 OSW_FALLBACK_WRITE_BUFFER_BYTES = 512ull * 1024ull * 1024ull;
    constexpr UInt64 OSW_MAX_WRITE_BUFFER_BYTES = 2ull * 1024ull * 1024ull * 1024ull;
    constexpr UInt64 OSW_SYSTEM_MEMORY_RESERVE_BYTES = 4ull * 1024ull * 1024ull * 1024ull;
    constexpr Size OSW_MIN_SCORE_JOB_COMPOUNDS = 1000;
    constexpr Size OSW_MAX_SCORE_JOB_COMPOUNDS = 10000;

    UInt64 determineOSWWriteBufferBytes()
    {
      size_t available_kb = 0;
      if (!SysInfo::getFreeSystemMemory(available_kb))
      {
        return OSW_FALLBACK_WRITE_BUFFER_BYTES;
      }

      const UInt64 available_kb_64 = static_cast<UInt64>(available_kb);
      UInt64 available_bytes = 0;
      if (available_kb_64 > std::numeric_limits<UInt64>::max() / 1024ull)
      {
        available_bytes = std::numeric_limits<UInt64>::max();
      }
      else
      {
        available_bytes = available_kb_64 * 1024ull;
      }

      const UInt64 usable_bytes = available_bytes > OSW_SYSTEM_MEMORY_RESERVE_BYTES ?
        available_bytes - OSW_SYSTEM_MEMORY_RESERVE_BYTES :
        available_bytes / 8ull;
      const UInt64 budget_bytes = usable_bytes / 4ull;
      return std::clamp(budget_bytes, OSW_MIN_WRITE_BUFFER_BYTES, OSW_MAX_WRITE_BUFFER_BYTES);
    }

    Size estimateScoreJobCount(const std::vector<Size>& compound_counts, Size job_size)
    {
      Size job_count = 0;
      for (Size count : compound_counts)
      {
        job_count += (count + job_size - 1) / job_size;
      }
      return job_count;
    }

    Size determineScoreJobCompoundCount(const std::vector<Size>& compound_counts)
    {
      if (compound_counts.empty())
      {
        return OSW_MIN_SCORE_JOB_COMPOUNDS;
      }

#ifdef _OPENMP
      const Size thread_count = static_cast<Size>(std::max(1, omp_get_max_threads()));
#else
      const Size thread_count = 1;
#endif
      const Size active_swath_count = compound_counts.size();
      const Size target_job_count = thread_count > active_swath_count ?
        active_swath_count + 1 :
        active_swath_count;
      const Size max_compounds = *std::max_element(compound_counts.begin(), compound_counts.end());

      Size lower = OSW_MIN_SCORE_JOB_COMPOUNDS;
      Size upper = std::max(OSW_MIN_SCORE_JOB_COMPOUNDS, max_compounds);
      while (lower < upper)
      {
        const Size middle = lower + (upper - lower) / 2;
        if (estimateScoreJobCount(compound_counts, middle) <= target_job_count)
        {
          upper = middle;
        }
        else
        {
          lower = middle + 1;
        }
      }
      return std::clamp(lower, OSW_MIN_SCORE_JOB_COMPOUNDS, OSW_MAX_SCORE_JOB_COMPOUNDS);
    }

    // Queue OSW rows so scoring threads only block when the memory budget is exhausted.
    class OSWBufferedWriter
    {
    public:
      struct EnqueueProfile
      {
        double wait_seconds = 0.0;
        Size row_count = 0;
        UInt64 byte_count = 0;
      };

      struct Stats
      {
        double enqueue_wait = 0.0;
        double write_hold = 0.0;
        Size queued_jobs = 0;
        Size flush_count = 0;
        Size queued_rows = 0;
        Size written_rows = 0;
        UInt64 queued_bytes = 0;
        UInt64 written_bytes = 0;
        UInt64 max_buffered_bytes = 0;
        UInt64 buffer_budget_bytes = 0;
      };

      OSWBufferedWriter(OpenSwathOSWWriter& writer, UInt64 buffer_budget_bytes, bool profile) :
        writer_(writer),
        buffer_budget_bytes_(std::max(buffer_budget_bytes, OSW_MIN_WRITE_BUFFER_BYTES)),
        profile_(profile)
      {
        stats_.buffer_budget_bytes = buffer_budget_bytes_;
        worker_ = std::thread(&OSWBufferedWriter::writerLoop_, this);
      }

      OSWBufferedWriter(const OSWBufferedWriter&) = delete;
      OSWBufferedWriter& operator=(const OSWBufferedWriter&) = delete;

      ~OSWBufferedWriter()
      {
        try
        {
          finish();
        }
        catch (...)
        {
        }
      }

      EnqueueProfile enqueue(OpenSwathOSWWriter::OSWData&& rows)
      {
        EnqueueProfile profile;
        if (rows.empty())
        {
          return profile;
        }

        profile.row_count = rows.rowCount();
        profile.byte_count = std::max<UInt64>(rows.estimateMemoryUsage(), 1ull);
        const auto wait_start = profileStart(profile_);

        std::unique_lock<std::mutex> lock(mutex_);
        checkExceptionLocked_();
        cv_.wait(lock, [&]()
        {
          return exception_ != nullptr ||
                 finish_requested_ ||
                 (buffered_bytes_ <= buffer_budget_bytes_ &&
                  profile.byte_count <= buffer_budget_bytes_ - buffered_bytes_) ||
                 (buffered_bytes_ == 0 && queue_.empty());
        });
        checkExceptionLocked_();
        if (finish_requested_)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Cannot enqueue OSW rows after the buffered writer was finished.");
        }

        if (profile_)
        {
          profile.wait_seconds = elapsedProfileSeconds(wait_start);
        }

        queue_.push_back({std::move(rows), profile.byte_count, profile.row_count});
        buffered_bytes_ += profile.byte_count;
        stats_.enqueue_wait += profile.wait_seconds;
        stats_.queued_jobs += 1;
        stats_.queued_rows += profile.row_count;
        stats_.queued_bytes += profile.byte_count;
        stats_.max_buffered_bytes = std::max(stats_.max_buffered_bytes, buffered_bytes_);

        lock.unlock();
        cv_.notify_all();
        return profile;
      }

      void finish()
      {
        {
          std::lock_guard<std::mutex> lock(mutex_);
          finish_requested_ = true;
        }
        cv_.notify_all();

        if (worker_.joinable())
        {
          worker_.join();
        }

        std::lock_guard<std::mutex> lock(mutex_);
        checkExceptionLocked_();
      }

      Stats stats() const
      {
        std::lock_guard<std::mutex> lock(mutex_);
        return stats_;
      }

    private:
      struct QueueItem
      {
        OpenSwathOSWWriter::OSWData rows;
        UInt64 byte_count = 0;
        Size row_count = 0;
      };

      void checkExceptionLocked_() const
      {
        if (exception_ != nullptr)
        {
          std::rethrow_exception(exception_);
        }
      }

      void writerLoop_()
      {
        try
        {
          while (true)
          {
            OpenSwathOSWWriter::OSWData batch;
            UInt64 batch_bytes = 0;
            Size batch_rows = 0;

            {
              std::unique_lock<std::mutex> lock(mutex_);
              cv_.wait(lock, [&]()
              {
                return exception_ != nullptr ||
                       finish_requested_ ||
                       !queue_.empty();
              });

              if (exception_ != nullptr)
              {
                return;
              }

              if (queue_.empty() && finish_requested_)
              {
                return;
              }
              if (queue_.empty())
              {
                continue;
              }

              while (!queue_.empty())
              {
                QueueItem item = std::move(queue_.front());
                queue_.pop_front();
                batch_bytes += item.byte_count;
                batch_rows += item.row_count;
                batch.append(std::move(item.rows));
              }
            }

            if (!batch.empty())
            {
              const auto write_start = profileStart(profile_);
              writer_.writeRows(batch);
              const double write_hold = profile_ ? elapsedProfileSeconds(write_start) : 0.0;

              std::lock_guard<std::mutex> lock(mutex_);
              stats_.write_hold += write_hold;
              stats_.flush_count += 1;
              stats_.written_rows += batch_rows;
              stats_.written_bytes += batch_bytes;
              buffered_bytes_ -= batch_bytes;
              cv_.notify_all();
            }
          }
        }
        catch (...)
        {
          std::lock_guard<std::mutex> lock(mutex_);
          exception_ = std::current_exception();
          queue_.clear();
          buffered_bytes_ = 0;
          cv_.notify_all();
        }
      }

      OpenSwathOSWWriter& writer_;
      UInt64 buffer_budget_bytes_ = OSW_FALLBACK_WRITE_BUFFER_BYTES;
      bool profile_ = false;
      mutable std::mutex mutex_;
      std::condition_variable cv_;
      std::deque<QueueItem> queue_;
      UInt64 buffered_bytes_ = 0;
      bool finish_requested_ = false;
      std::exception_ptr exception_;
      Stats stats_;
      std::thread worker_;
    };

    void logOSWBufferedWriterProfile(const OSWBufferedWriter::Stats& stats)
    {
      OPENMS_LOG_INFO << "[OpenSwathWorkflow profile] OSW writer queue: "
                      << "wait=" << formatProfileSeconds(stats.enqueue_wait)
                      << ", write_hold=" << formatProfileSeconds(stats.write_hold)
                      << ", flushes=" << stats.flush_count
                      << ", jobs=" << stats.queued_jobs
                      << ", rows=" << stats.written_rows
                      << ", budget=" << bytesToHumanReadable(stats.buffer_budget_bytes)
                      << ", max_buffered=" << bytesToHumanReadable(stats.max_buffered_bytes)
                      << ", queued=" << bytesToHumanReadable(stats.queued_bytes)
                      << ", written=" << bytesToHumanReadable(stats.written_bytes)
                      << std::endl;
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
         << ", osw_write_wait=" << formatProfileSeconds(timing.osw_write_wait)
         << ", osw_write_hold=" << formatProfileSeconds(timing.osw_write_hold)
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

#ifdef _OPENMP
#ifdef MT_ENABLE_NESTED_OPENMP
    const bool nested_scheduler_requested = threads_outer_loop_ > -1;
#else
    const bool nested_scheduler_requested = false;
#endif
    const bool use_swath_range_scheduler = batchSize <= 0 && load_into_memory && !ms1_only &&
                                          mobilogram_consumer == nullptr && !nested_scheduler_requested;
#else
    const bool use_swath_range_scheduler = false;
#endif
    if (use_swath_range_scheduler)
    {
      struct SwathSchedulerContext
      {
        SignedSize swath_index = -1;
        bool active = false;
        int score_job_size = 0;
        SignedSize nr_score_jobs = 0;
        SignedSize remaining_score_jobs = 0;
        bool finalized = false;
        OpenSwathWorkflowPhaseTiming swath_timing;
        OpenSwath::LightTargetedExperiment transition_exp_used_all;
        OpenSwath::SpectrumAccessPtr current_swath_map;
        std::vector<MSChromatogram> ms1_chromatograms;
        PeakMap chrom_exp;
        FeatureMap feature_file;
      };

      struct SwathScoreJob
      {
        Size context_index = 0;
        SignedSize range_index = 0;
      };

      std::vector<SwathSchedulerContext> contexts(swath_maps.size());
      std::unique_ptr<OSWBufferedWriter> buffered_osw_writer;
      if (osw_writer.isActive())
      {
        const UInt64 osw_buffer_bytes = determineOSWWriteBufferBytes();
        OPENMS_LOG_INFO << "Use buffered OSW writer with "
                        << bytesToHumanReadable(osw_buffer_bytes)
                        << " queue budget and a single writer thread." << std::endl;
        buffered_osw_writer = std::make_unique<OSWBufferedWriter>(osw_writer, osw_buffer_bytes, profile);
      }

#ifdef _OPENMP
      OPENMS_LOG_INFO << "Use SWATH assay-range scheduler with " << omp_get_max_threads()
                      << " threads." << std::endl;
#pragma omp parallel for schedule(dynamic,1)
#endif
      for (SignedSize swath_index = 0; swath_index < boost::numeric_cast<SignedSize>(swath_maps.size()); ++swath_index)
      {
        SwathSchedulerContext& context = contexts[swath_index];
        context.swath_index = swath_index;

        if (swath_maps[swath_index].ms1)
        {
          continue;
        }

        auto phase_start = profileStart(profile);
        if (!(prm_ || pasef_))
        {
          OpenSwathHelper::selectSwathTransitions(transition_exp, context.transition_exp_used_all,
              cp.min_upper_edge_dist, swath_maps[swath_index].lower, swath_maps[swath_index].upper);
        }
        else
        {
          std::set<std::string> matching_compounds;
          for (Size k = 0; k < tr_win_map.size(); k++)
          {
            if (tr_win_map[k] == swath_index)
            {
              const OpenSwath::LightTransition& tr = transition_exp.transitions[k];
              context.transition_exp_used_all.transitions.push_back(tr);
              matching_compounds.insert(tr.getPeptideRef());
              OPENMS_LOG_DEBUG << "Adding Precursor with m/z " << tr.getPrecursorMZ()
                               << " and IM of " << tr.getPrecursorIM()
                               << " to swath with mz upper of " << swath_maps[swath_index].upper
                               << " im lower of " << swath_maps[swath_index].imLower
                               << " and im upper of " << swath_maps[swath_index].imUpper << '\n';
            }
          }

          std::set<std::string> matching_proteins;
          for (Size compound_idx = 0; compound_idx < transition_exp.compounds.size(); compound_idx++)
          {
            if (matching_compounds.find(transition_exp.compounds[compound_idx].id) != matching_compounds.end())
            {
              context.transition_exp_used_all.compounds.push_back(transition_exp.compounds[compound_idx]);
              for (Size protein_ref_idx = 0;
                   protein_ref_idx < transition_exp.compounds[compound_idx].protein_refs.size();
                   protein_ref_idx++)
              {
                matching_proteins.insert(transition_exp.compounds[compound_idx].protein_refs[protein_ref_idx]);
              }
            }
          }
          for (Size protein_idx = 0; protein_idx < transition_exp.proteins.size(); protein_idx++)
          {
            if (matching_proteins.find(transition_exp.proteins[protein_idx].id) != matching_proteins.end())
            {
              context.transition_exp_used_all.proteins.push_back(transition_exp.proteins[protein_idx]);
            }
          }
        }
        if (profile)
        {
          context.swath_timing.transition_selection = elapsedProfileSeconds(phase_start);
          if (context.transition_exp_used_all.getTransitions().empty())
          {
            context.swath_timing.compound_count = context.transition_exp_used_all.getCompounds().size();
            context.swath_timing.transition_count = context.transition_exp_used_all.getTransitions().size();
          }
        }

        if (context.transition_exp_used_all.getTransitions().empty())
        {
          continue;
        }

        phase_start = profileStart(profile);
        context.current_swath_map = swath_maps[swath_index].sptr;
        context.current_swath_map = std::shared_ptr<SpectrumAccessOpenMSInMemory>(
          new SpectrumAccessOpenMSInMemory(*context.current_swath_map));
        if (profile)
        {
          context.swath_timing.swath_loading = elapsedProfileSeconds(phase_start);
        }

#ifdef _OPENMP
#pragma omp critical (osw_write_stdout)
#endif
        {
          OPENMS_LOG_INFO << "Thread " <<
#ifdef _OPENMP
            omp_get_thread_num() << "_0 " <<
#else
            "0" <<
#endif
            "will extract " << context.transition_exp_used_all.getCompounds().size() << " compounds and "
                            << context.transition_exp_used_all.getTransitions().size() << " transitions "
                            << "from SWATH " << swath_index << std::endl;
        }

        phase_start = profileStart(profile);
        if (ms1_map_ != nullptr)
        {
          OpenSwath::SpectrumAccessPtr threadsafe_ms1 = ms1_map_->lightClone();
          MS1Extraction_(threadsafe_ms1, swath_maps, context.ms1_chromatograms, ms1_cp,
              context.transition_exp_used_all, trafo_inverse, ms1_only, ms1_isotopes);
        }
        if (profile)
        {
          context.swath_timing.ms1_extraction = elapsedProfileSeconds(phase_start);
          context.swath_timing.ms1_chromatogram_count = context.ms1_chromatograms.size();
        }

        ChromatogramExtractor extractor;
        std::vector<OpenSwath::ChromatogramPtr> chrom_list;
        std::vector<ChromatogramExtractor::ExtractionCoordinates> coordinates;

        phase_start = profileStart(profile);
        prepareExtractionCoordinates_(chrom_list, coordinates, context.transition_exp_used_all, trafo_inverse, cp);
        if (profile)
        {
          context.swath_timing.coordinate_preparation = elapsedProfileSeconds(phase_start);
        }

        phase_start = profileStart(profile);
        extractor.extractChromatograms(context.current_swath_map, chrom_list, coordinates, cp.mz_extraction_window,
            cp.ppm, cp.im_extraction_window, cp.extraction_function);
        if (profile)
        {
          context.swath_timing.ms2_extraction = elapsedProfileSeconds(phase_start);
        }

        phase_start = profileStart(profile);
        extractor.return_chromatogram(chrom_list, coordinates, context.transition_exp_used_all, SpectrumSettings(),
                                      context.chrom_exp.getChromatograms(), false, cp.im_extraction_window);
        if (profile)
        {
          context.swath_timing.chromatogram_conversion = elapsedProfileSeconds(phase_start);
          context.swath_timing.ms2_chromatogram_count = context.chrom_exp.getChromatograms().size();
        }

        context.active = true;
      }

      std::vector<SwathScoreJob> score_jobs;
      std::vector<Size> active_compound_counts;
      for (const SwathSchedulerContext& context : contexts)
      {
        if (context.active)
        {
          active_compound_counts.push_back(context.transition_exp_used_all.getCompounds().size());
        }
      }

      const Size score_job_compound_count = determineScoreJobCompoundCount(active_compound_counts);
      OPENMS_LOG_INFO << "Use " << score_job_compound_count
                      << " compounds per scoring job." << std::endl;

      for (Size context_index = 0; context_index < contexts.size(); ++context_index)
      {
        SwathSchedulerContext& context = contexts[context_index];
        if (!context.active)
        {
          continue;
        }

        const Size n_compounds = context.transition_exp_used_all.getCompounds().size();
        context.score_job_size = boost::numeric_cast<int>(std::min(score_job_compound_count, n_compounds));
        if (context.score_job_size > 0)
        {
          context.nr_score_jobs = static_cast<SignedSize>((n_compounds + context.score_job_size - 1) / context.score_job_size);
          context.remaining_score_jobs = context.nr_score_jobs;
        }

        OPENMS_LOG_INFO << "SWATH " << context.swath_index
                        << " will score " << context.transition_exp_used_all.getCompounds().size()
                        << " compounds and " << context.transition_exp_used_all.getTransitions().size()
                        << " transitions in " << context.nr_score_jobs
                        << " scoring jobs" << std::endl;

        for (SignedSize range_index = 0; range_index < context.nr_score_jobs; ++range_index)
        {
          score_jobs.push_back({context_index, range_index});
        }
      }

      auto finalize_context = [&](const Size context_index)
      {
        SwathSchedulerContext& context = contexts[context_index];
        if (!context.active || context.finalized)
        {
          return;
        }

        double writing_seconds = 0.0;
#ifdef _OPENMP
#pragma omp critical (osw_write_out)
#endif
        {
          auto write_start = profileStart(profile);
          writeOutFeaturesAndChroms_(context.chrom_exp.getChromatograms(), context.ms1_chromatograms,
              context.feature_file, out_featureFile, store_features, chromConsumer);
          if (profile)
          {
            writing_seconds = elapsedProfileSeconds(write_start);
          }
        }
        if (profile)
        {
          context.swath_timing.writing += writing_seconds;
          if (!swath_maps[context_index].ms1)
          {
            context.swath_timing.swath_count = 1;
            context.swath_timing.total = context.swath_timing.transition_selection +
                                         context.swath_timing.swath_loading +
                                         context.swath_timing.batch_selection +
                                         context.swath_timing.ms1_extraction +
                                         context.swath_timing.coordinate_preparation +
                                         context.swath_timing.ms2_extraction +
                                         context.swath_timing.chromatogram_conversion +
                                         context.swath_timing.scoring +
                                         context.swath_timing.writing;
#ifdef _OPENMP
#pragma omp critical (osw_scheduler_profile)
#endif
            {
              logSwathProfile(context.swath_index, context.swath_timing);
              profile_total.add(context.swath_timing);
            }
          }
        }

        context.finalized = true;
#ifdef _OPENMP
#pragma omp critical (progress)
#endif
        {
          this->setProgress(++progress);
        }
      };

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
      for (SignedSize job_index = 0; job_index < boost::numeric_cast<SignedSize>(score_jobs.size()); ++job_index)
      {
        const SwathScoreJob& job = score_jobs[job_index];
        SwathSchedulerContext& context = contexts[job.context_index];

        OpenSwathWorkflowPhaseTiming batch_timing;
        const auto batch_profile_start = profileStart(profile);
        batch_timing.batch_count = 1;

        auto batch_phase_start = profileStart(profile);
        OpenSwath::LightTargetedExperiment transition_exp_used;
        selectCompoundsForBatch_(context.transition_exp_used_all, transition_exp_used, context.score_job_size, job.range_index);
        if (profile)
        {
          batch_timing.batch_selection = elapsedProfileSeconds(batch_phase_start);
          batch_timing.compound_count = transition_exp_used.getCompounds().size();
          batch_timing.transition_count = transition_exp_used.getTransitions().size();
        }

        batch_phase_start = profileStart(profile);
        FeatureMap featureFile;
        std::vector<OpenSwath::SwathMap> tmp = {swath_maps[context.swath_index]};
        tmp.back().sptr = context.current_swath_map->lightClone();
        OpenSwathScoringPhaseTiming scoring_profile;
        OpenSwathOSWWriter::OSWData job_osw_output;
        scoreAllChromatograms_(context.chrom_exp.getChromatograms(), context.ms1_chromatograms, tmp, transition_exp_used,
          feature_finder_param, trafo, cp.rt_extraction_window, featureFile, osw_writer, ms1_isotopes, false,
          mobilogram_consumer, profile ? &scoring_profile : nullptr, &job_osw_output);
        if (profile)
        {
          batch_timing.scoring = elapsedProfileSeconds(batch_phase_start);
          batch_timing.feature_count = featureFile.size();
        }

        double job_osw_wait_seconds = 0.0;
        if (osw_writer.isActive() && !job_osw_output.empty())
        {
          OSWBufferedWriter::EnqueueProfile osw_enqueue_profile =
            buffered_osw_writer->enqueue(std::move(job_osw_output));
          if (profile)
          {
            job_osw_wait_seconds = osw_enqueue_profile.wait_seconds;
            scoring_profile.osw_write_wait += job_osw_wait_seconds;
            scoring_profile.osw_write += job_osw_wait_seconds;
          }
        }

        if (profile)
        {
          batch_timing.writing += job_osw_wait_seconds;
          batch_timing.scoring_breakdown.add(scoring_profile);
          batch_timing.total = elapsedProfileSeconds(batch_profile_start);
        }

        bool context_ready_to_finalize = false;
#ifdef _OPENMP
#pragma omp critical (osw_scheduler_context)
#endif
        {
          if (!osw_writer.isActive() && store_features)
          {
            for (FeatureMap::const_iterator feature_it = featureFile.begin(); feature_it != featureFile.end(); ++feature_it)
            {
              context.feature_file.push_back(*feature_it);
            }
            for (std::vector<ProteinIdentification>::const_iterator protid_it = featureFile.getProteinIdentifications().begin();
                 protid_it != featureFile.getProteinIdentifications().end(); ++protid_it)
            {
              context.feature_file.getProteinIdentifications().push_back(*protid_it);
            }
          }

          if (profile)
          {
            context.swath_timing.add(batch_timing);
            logBatchProfile(context.swath_index, job.range_index, context.nr_score_jobs, batch_timing);
          }

          --context.remaining_score_jobs;
          context_ready_to_finalize = context.remaining_score_jobs == 0;
        }

        if (context_ready_to_finalize)
        {
          finalize_context(job.context_index);
        }
      }

      for (Size context_index = 0; context_index < contexts.size(); ++context_index)
      {
        SwathSchedulerContext& context = contexts[context_index];
        if (context.active && !context.finalized)
        {
          finalize_context(context_index);
        }
        else if (!context.active)
        {
          this->setProgress(++progress);
        }
      }

      if (buffered_osw_writer != nullptr)
      {
        buffered_osw_writer->finish();
        if (profile)
        {
          const OSWBufferedWriter::Stats osw_writer_stats = buffered_osw_writer->stats();
          profile_total.scoring_breakdown.osw_write_hold += osw_writer_stats.write_hold;
          profile_total.scoring_breakdown.osw_write += osw_writer_stats.write_hold;
          logOSWBufferedWriterProfile(osw_writer_stats);
        }
      }
      this->endProgress();

      if (profile)
      {
        logWorkflowProfileSummary(elapsedProfileSeconds(profile_workflow_start), profile_total);
      }
      return;
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
    OpenSwathScoringPhaseTiming* scoring_profile,
    OpenSwathOSWWriter::OSWData* deferred_osw_output) const
  {
    const bool profile = scoring_profile != nullptr;
    if (profile)
    {
      *scoring_profile = OpenSwathScoringPhaseTiming();
    }

    auto profile_phase_start = profileStart(profile);
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
    featureFinder.setScoringProfiling(profile);
    featureFinder.resetScoringProfile();
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
    if (profile)
    {
      scoring_profile->map_setup += elapsedProfileSeconds(profile_phase_start);
    }

    OpenSwathOSWWriter::OSWData to_osw_output;
    if (osw_writer.isActive())
    {
      const int stop_report_after_feature =
        static_cast<int>(feature_finder_param.getValue("stop_report_after_feature"));
      const Size expected_features_per_assay = stop_report_after_feature > 0 ?
        static_cast<Size>(stop_report_after_feature) :
        5;
      to_osw_output.reserve(transition_exp.getCompounds().size() * expected_features_per_assay,
                            transition_exp.getTransitions().size() * expected_features_per_assay);
    }
    ///////////////////////////////////
    // Start of main function
    // Iterating over all the assays
    ///////////////////////////////////
    for (AssayMapT::iterator assay_it = assay_map.begin(); assay_it != assay_map.end(); ++assay_it)
    {
      auto assay_setup_start = profileStart(profile);
      bool assay_setup_recorded = false;
      auto recordAssaySetup = [&]()
      {
        if (profile && !assay_setup_recorded)
        {
          scoring_profile->assay_setup += elapsedProfileSeconds(assay_setup_start);
          assay_setup_recorded = true;
        }
      };
      if (profile)
      {
        ++scoring_profile->assay_count;
      }

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
        recordAssaySetup();
        if (profile)
        {
          ++scoring_profile->skipped_assay_count;
        }
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
          recordAssaySetup();
          if (profile)
          {
            ++scoring_profile->skipped_assay_count;
          }
          continue;
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
          scoring_profile->transition_group_picker += elapsedProfileSeconds(picker_profile_start);
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
          ++scoring_profile->skipped_assay_count;
        }
        continue;
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
          ++scoring_profile->skipped_assay_count;
        }
        continue;
      }
      catch (const std::exception & e)
      {
        recordPickerProfile();
        OPENMS_LOG_ERROR << "Exception while picking transition group for assay " << id
                         << ": " << e.what() << " - skipping assay." << std::endl;
        if (profile)
        {
          ++scoring_profile->skipped_assay_count;
        }
        continue;
      }

      auto score_profile_start = profileStart(profile);
      bool score_profile_recorded = false;
      auto recordScoreProfile = [&]()
      {
        if (profile && !score_profile_recorded)
        {
          scoring_profile->score_peakgroups += elapsedProfileSeconds(score_profile_start);
          score_profile_recorded = true;
        }
      };
      try
      {
        featureFinder.scorePeakgroups(transition_group, trafo, swath_maps, output, ms1only, mobilogram_consumer);
        recordScoreProfile();
        if (profile)
        {
          ++scoring_profile->scored_assay_count;
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
          ++scoring_profile->skipped_assay_count;
        }
        continue;
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
          ++scoring_profile->skipped_assay_count;
        }
        continue;
      }
      catch (const std::exception & e)
      {
        recordScoreProfile();
        OPENMS_LOG_ERROR << "Exception while scoring transition group for assay " << id
                         << ": " << e.what() << " - skipping assay." << std::endl;
        if (profile)
        {
          ++scoring_profile->skipped_assay_count;
        }
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
        auto osw_prepare_start = profileStart(profile);
        // Compound and transition are currently unused by the OSW row writer.
        osw_writer.prepareRowsInto(to_osw_output,
                                   OpenSwath::LightCompound(),
                                   nullptr,
                                   output,
                                   id);
        if (profile)
        {
          scoring_profile->osw_prepare += elapsedProfileSeconds(osw_prepare_start);
        }
      }
    }

    if (profile)
    {
      scoring_profile->add(featureFinder.getScoringProfile());
    }

    // Only write at the very end since this is a step that needs a barrier.
    // Schedulers can defer these rows to avoid blocking worker threads on the writer.
    if (osw_writer.isActive())
    {
      if (deferred_osw_output != nullptr)
      {
        deferred_osw_output->append(std::move(to_osw_output));
        return;
      }

      double osw_write_wait = 0.0;
      double osw_write_hold = 0.0;
      const auto osw_write_wait_start = profileStart(profile);
#ifdef _OPENMP
#pragma omp critical (osw_write_tsv)
#endif
      {
        if (profile)
        {
          osw_write_wait = elapsedProfileSeconds(osw_write_wait_start);
        }
        auto osw_write_start = profileStart(profile);
        osw_writer.writeRows(to_osw_output);
        if (profile)
        {
          osw_write_hold = elapsedProfileSeconds(osw_write_start);
        }
      }
      if (profile)
      {
        scoring_profile->osw_write_wait += osw_write_wait;
        scoring_profile->osw_write_hold += osw_write_hold;
        scoring_profile->osw_write += osw_write_wait + osw_write_hold;
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
