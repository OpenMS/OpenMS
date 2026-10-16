// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// The database search engine ProSE (ANALYSIS/ID/ProSEAlgorithm): a search of one spectrum
// file as the ProSE tool runs it, with its statistics, the modification analysis of open
// searches and progress reports; a search of several files pooled as the tool's -out_merged
// does, whose results (ProSEResultSet) JavaScript holds; the annotation of a spectrum, or of the
// stored peaks of an identified spectrum, with the fragment ions of a peptide, for spectrum viewers; and the fragment
// intensities that PeptDeep predicts for a peptide (modules with ONNX Runtime).

#include "bind_prose.h"
#include "binding_utils.h"

#include <OpenMS/ANALYSIS/ID/PeptDeepRescoring.h>
#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/config.h>

#ifdef WITH_ONNX
#include <OpenMS/ML/PEPTDEEP/PeptDeepMS2Inference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>
#endif

#include <pthread.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <atomic>
#include <charconv>
#include <cctype>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace emscripten;

namespace OpenMS::WASM
{
  using ProteinIdentificationList = std::vector<ProteinIdentification>;

  namespace
  {
    //------------------------------------------------------------------------------
    // Progress

    /// The position of the input file of a progress step "Searching <name> (<i> of <n>)...", the step of a
    /// multi-file search per input file (ProSEAlgorithm): i - 1, and whether it is the last file (i = n);
    /// nothing for other labels.
    std::optional<std::pair<int, bool>> fileOfProgressStep(const std::string& label)
    {
      static constexpr std::string_view prefix = "Searching ", suffix = ")...", of = " of ";
      const std::string_view l(label);
      if (l.size() <= prefix.size() + suffix.size() || l.substr(0, prefix.size()) != prefix
          || l.substr(l.size() - suffix.size()) != suffix)
      {
        return std::nullopt;
      }
      const auto open = l.rfind(" (");
      if (open == std::string_view::npos || open < prefix.size()) return std::nullopt;
      const std::string_view counts = l.substr(open + 2, l.size() - suffix.size() - open - 2); // "<i> of <n>"
      const auto of_pos = counts.find(of);
      if (of_pos == std::string_view::npos) return std::nullopt;
      const auto number = [](std::string_view s, int& n) {
        const auto [end, error] = std::from_chars(s.data(), s.data() + s.size(), n);
        return error == std::errc() && end == s.data() + s.size() && !s.empty();
      };
      int i = 0, n = 0;
      if (!number(counts.substr(0, of_pos), i) || !number(counts.substr(of_pos + of.size()), n) || i < 1 || i > n)
      {
        return std::nullopt;
      }
      return std::make_pair(i - 1, i == n);
    }

    class JSProgressLogger;

    /// The JSProgressLogger of each ProgressLogger that has one, so that a search can reset it (resetProgress()):
    /// ProgressLogger has no getter for its implementation. A JSProgressLogger enters itself when it is created and
    /// leaves when its ProgressLogger deletes it (for another callback or log type, or with the algorithm).
    std::unordered_map<const ProgressLogger*, JSProgressLogger*>& jsProgressLoggers()
    {
      static std::unordered_map<const ProgressLogger*, JSProgressLogger*> loggers;
      return loggers;
    }

    std::mutex& jsProgressLoggersMutex()
    {
      static std::mutex mutex;
      return mutex;
    }

    /// Reports the progress of an algorithm to a JavaScript function, called with
    /// {event: "start" | "progress" | "end", label, depth, begin, end, value, file}:
    /// - Every report is that of the innermost open step: its label, its range [begin, end] and
    ///   its nesting depth (0: outermost), as this logger counts the steps; value is in the range.
    /// - file: the position of the input file being searched, from the step "Searching <name>
    ///   (<i> of <n>)..." of a multi-file search: i - 1 from that step on, until the step of the
    ///   next file. -1 before the first file, after the step of the last one ends, and in steps
    ///   outside the files.
    /// The function is called only on the thread that set it (a JavaScript value belongs to the
    /// thread that created it), and only the steps of that thread count: the steps of other
    /// threads (e.g. a background reader) are not reported, and their values only count towards
    /// the open step (OpenMP threads). It runs while the search runs: in a Web Worker,
    /// postMessage() its argument to the page.
    /// A search resets it before it starts (reset()): the steps of the algorithms are not ended
    /// when they throw, so a search that threw may have left steps open.
    class JSProgressLogger : public ProgressLogger::ProgressLoggerImpl
    {
    public:
      /// Reports the progress of @p logger, which takes it (setLogger()), to @p callback
      JSProgressLogger(val callback, const ProgressLogger& logger) :
        callback_(std::move(callback)),
        owner_(pthread_self()),
        logger_(&logger)
      {
        const std::lock_guard<std::mutex> lock(jsProgressLoggersMutex());
        jsProgressLoggers()[logger_] = this;
      }

      JSProgressLogger(const JSProgressLogger&) = delete;
      JSProgressLogger& operator=(const JSProgressLogger&) = delete;

      ~JSProgressLogger() override
      {
        const std::lock_guard<std::mutex> lock(jsProgressLoggersMutex());
        const auto it = jsProgressLoggers().find(logger_);
        if (it != jsProgressLoggers().end() && it->second == this) jsProgressLoggers().erase(it);
      }

      /// Forgets the open steps and the file, before a search starts (on the thread of the search, with no step open)
      void reset() const
      {
        steps_.clear();
        file_ = -1;
        current_ = 0;
      }

      // ProgressLogger's depth is not used: it is static and shared by all threads.
      void startProgress(const SignedSize begin, const SignedSize end, const std::string& label, const int /*depth*/) const override
      {
        if (!isOwner_()) return;
        if (!steps_.empty()) steps_.back().value = current_;
        Step step{label, begin, end, begin, false};
        if (const auto file = fileOfProgressStep(label))
        {
          file_ = file->first;
          step.last_file = file->second;
        }
        else if (steps_.empty())
        {
          file_ = -1; // a step outside the files
        }
        steps_.push_back(std::move(step));
        current_ = begin;
        report_("start", steps_.back(), steps_.size() - 1, begin);
      }

      void setProgress(const SignedSize value, const int /*depth*/) const override
      {
        if (!isOwner_() || steps_.empty()) return;
        current_ = value;
        report_("progress", steps_.back(), steps_.size() - 1, value);
      }

      SignedSize nextProgress() const override
      {
        return ++current_;
      }

      void endProgress(const int /*depth*/, UInt64 /*bytes_processed*/) const override
      {
        if (!isOwner_() || steps_.empty()) return;
        const Step step = std::move(steps_.back());
        steps_.pop_back();
        if (!steps_.empty()) current_ = steps_.back().value;
        report_("end", step, steps_.size(), step.end);
        if (step.last_file) file_ = -1; // after the search of the files
      }

    private:
      /// An open step
      struct Step
      {
        std::string label;
        SignedSize begin, end;
        SignedSize value;   ///< the progress when a step inside it started
        bool last_file;     ///< the step of the last input file
      };

      bool isOwner_() const
      {
        return pthread_equal(pthread_self(), owner_);
      }

      /// Reports @p event of @p step, at nesting depth @p depth.
      void report_(const char* event, const Step& step, Size depth, SignedSize value) const
      {
        val o = val::object();
        o.set("event", std::string(event));
        o.set("label", step.label);
        o.set("depth", static_cast<int>(depth));
        o.set("begin", static_cast<double>(step.begin));
        o.set("end", static_cast<double>(step.end));
        o.set("value", static_cast<double>(value));
        o.set("file", file_);
        callback_(o);
      }

      val callback_;
      pthread_t owner_;
      const ProgressLogger* logger_; ///< the logger that has this one
      /// the progress of the innermost open step (all threads count)
      mutable std::atomic<SignedSize> current_{0};
      /// the open steps of the owner thread, innermost last (only that thread changes them)
      mutable std::vector<Step> steps_;
      mutable int file_ = -1;
    };

    /// Sets the function that receives the progress reports of @p algo (null or undefined: none).
    void setProgressCallback(ProgressLogger& logger, const JS::ProgressCallback& callback)
    {
      if (isNullish(callback))
      {
        logger.setLogType(ProgressLogger::NONE);
        return;
      }
      if (callback.typeOf().as<std::string>() != "function")
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "The progress callback has to be a function (or null).", "");
      }
      logger.setLogger(new JSProgressLogger(callback, logger)); // the logger owns it
    }

    /// Resets the progress reports of @p logger to JavaScript, if it has them (JSProgressLogger::reset()): at the
    /// start of every search, which has to report its steps from the outermost one.
    void resetProgress(const ProgressLogger& logger)
    {
      const std::lock_guard<std::mutex> lock(jsProgressLoggersMutex());
      const auto it = jsProgressLoggers().find(&logger);
      if (it != jsProgressLoggers().end()) it->second->reset();
    }

    //------------------------------------------------------------------------------
    // Search

    std::string exitCodeName(ProSEAlgorithm::ExitCodes code)
    {
      switch (code)
      {
        case ProSEAlgorithm::ExitCodes::EXECUTION_OK: return "EXECUTION_OK";
        case ProSEAlgorithm::ExitCodes::INPUT_FILE_EMPTY: return "INPUT_FILE_EMPTY";
        case ProSEAlgorithm::ExitCodes::UNEXPECTED_RESULT: return "UNEXPECTED_RESULT";
        case ProSEAlgorithm::ExitCodes::ILLEGAL_PARAMETERS: return "ILLEGAL_PARAMETERS";
        default: return "UNKNOWN_ERROR";
      }
    }

    /// A std::map<K, Size> as {key: count}
    template <typename Map>
    val histogramToJS(const Map& m)
    {
      val o = val::object();
      for (const auto& [key, count] : m)
      {
        o.set(std::to_string(key), static_cast<double>(count));
      }
      return o;
    }

    val runStatisticsToJS(const ProSEAlgorithm::RunStatistics& s)
    {
      val o = val::object();
      o.set("inputFile", s.input_file);
      o.set("ms2Spectra", static_cast<double>(s.ms2_spectra));
      o.set("matchedSpectra", static_cast<double>(s.matched_spectra));
      o.set("targetPSMs", static_cast<double>(s.target_psms));
      o.set("decoyPSMs", static_cast<double>(s.decoy_psms));
      o.set("fdrApplied", s.fdr_applied);
      o.set("achievedPsmFDR", s.achieved_psm_fdr < 0 ? val::null() : val(s.achieved_psm_fdr));
      o.set("uniquePeptides", static_cast<double>(s.unique_peptides));
      o.set("uniqueProteins", static_cast<double>(s.unique_proteins));
      o.set("chargeHistogram", histogramToJS(s.charge_histogram));
      o.set("missedCleavageHistogram", histogramToJS(s.missed_cleavage_histogram));
      if (s.score_stats_valid)
      {
        val h = val::object();
        h.set("min", s.hyperscore_min);
        h.set("median", s.hyperscore_median);
        h.set("max", s.hyperscore_max);
        o.set("hyperscore", h);
      }
      const auto error = [](bool valid, double median, double mad, double recommended) {
        if (!valid) return val::null();
        val e = val::object();
        e.set("median", median);
        e.set("mad", mad);
        e.set("recommended", recommended);
        return e;
      };
      o.set("precursorError", error(s.prec_tol_valid, s.prec_err_median, s.prec_err_mad, s.prec_err_recommended));
      o.set("fragmentError", error(s.frag_tol_valid, s.frag_err_median, s.frag_err_mad, s.frag_err_recommended));
      o.set("secondsSearch", s.seconds_search);
      o.set("secondsCalibration", s.seconds_calibration);
      o.set("secondsFDR", s.seconds_fdr);
      return o;
    }

    val sharedStatisticsToJS(const ProSEAlgorithm::SharedSearchStats& s)
    {
      val o = val::object();
      o.set("databaseFile", s.database_file);
      o.set("enzyme", s.enzyme);
      o.set("openSearch", s.open_search);
      o.set("calibration", s.calibration_enabled);
      o.set("chunked", s.chunked);
      o.set("decoyMode", s.decoy_mode);
      o.set("targetProteins", static_cast<double>(s.db_target_proteins));
      o.set("decoyProteins", static_cast<double>(s.db_decoy_proteins));
      o.set("indexedPeptides", static_cast<double>(s.indexed_peptides));
      o.set("indexedFragments", static_cast<double>(s.indexed_fragments));
      o.set("secondsIndexBuild", s.seconds_index_build);
      o.set("secondsTotal", s.seconds_total);
      return o;
    }

    val modificationAnalysisToJS(const OpenSearchModificationAnalysis::OpenSearchAnalysisResult& r)
    {
      const auto& d = r.delta_mass_stats;
      val delta = val::object();
      delta.set("totalPSMs", d.total_psms);
      delta.set("modifiedPSMs", d.modified_psms);
      delta.set("unmodifiedPSMs", d.unmodified_psms);
      delta.set("meanDeltaMass", d.mean_delta_mass);
      delta.set("medianDeltaMass", d.median_delta_mass);
      delta.set("entries", toArray(d.entries, [](const OpenSearchModificationAnalysis::DeltaMassEntry& e) {
        val o = val::object();
        o.set("deltaMass", e.delta_mass);
        o.set("count", e.count);
        o.set("uniquePeptides", e.unique_peptides);
        o.set("chargeStates", e.num_charge_states);
        o.set("percentage", e.percentage);
        o.set("modification", e.mapped_modification);
        o.set("known", e.is_known_modification);
        return o;
      }));

      const auto& p = r.ptm_stats;
      val ptms = val::object();
      ptms.set("totalModifiedPSMs", p.total_modified_psms);
      ptms.set("unknownModificationPSMs", p.unknown_modification_psms);
      ptms.set("uniqueModifications", p.num_unique_modifications);
      ptms.set("entries", toArray(p.entries, [](const OpenSearchModificationAnalysis::PTMEntry& e) {
        val o = val::object();
        o.set("name", e.name);
        o.set("theoreticalMass", e.theoretical_mass);
        o.set("observedMass", e.observed_mass);
        o.set("massDeviation", e.mass_deviation);
        o.set("count", e.count);
        o.set("uniquePeptides", e.unique_peptides);
        o.set("chargeStates", e.num_charge_states);
        o.set("percentage", e.percentage);
        o.set("targetResidues", e.target_residues);
        val residues = val::object();
        for (const auto& [residue, count] : e.residue_counts)
        {
          residues.set(std::string(1, residue), count);
        }
        o.set("residueCounts", residues);
        return o;
      }));

      val o = val::object();
      o.set("deltaMasses", delta);
      o.set("ptms", ptms);
      return o;
    }

    /// The option @p name of @p options (an object or nullish), or undefined.
    val option(const val& options, const char* name)
    {
      return isNullish(options) ? val::undefined() : options[name];
    }

    /// Sets the parameter 'rescore' of an algorithm to the option {rescore} of a call, if given, for as long as it
    /// exists; the parameters are restored afterwards, also when the call throws.
    class RescoreOverride
    {
    public:
      RescoreOverride(ProSEAlgorithm& algo, const val& options) :
        algo_(algo),
        saved_(algo.getParameters())
      {
        const val rescore = option(options, "rescore");
        if (isNullish(rescore)) return;
        Param p = saved_;
        p.setValue("rescore", rescore.as<bool>() ? "true" : "false");
        algo_.setParameters(p);
        active_ = true;
      }

      RescoreOverride(const RescoreOverride&) = delete;
      RescoreOverride& operator=(const RescoreOverride&) = delete;

      ~RescoreOverride()
      {
        if (!active_) return;
        try
        {
          algo_.setParameters(saved_);
        }
        catch (...) // the parameters were valid before
        {
        }
      }

      /// Whether the call rescores the PSMs (the parameter 'rescore' as set for it)
      bool rescore() const
      {
        return algo_.getParameters().getValue("rescore").toString() == "true";
      }

    private:
      ProSEAlgorithm& algo_;
      const Param saved_;
      bool active_ = false;
    };

    /// Searches the spectra of @p spectra_path (mzML) against the proteins of @p database_path
    /// (FASTA) as the ProSE tool does for one input file: the search, which by default (parameter
    /// 'rescore') rescores the PSMs with Percolator (in-process) before its PSM-level FDR
    /// (FDR:PSM), then, if FDR:protein is set and the database has decoys, protein inference with
    /// picked-protein FDR. Options {rescore} set the parameter 'rescore' for this call.
    JS::ProSESummary proseSearchWithOptions(ProSEAlgorithm& algo, const std::string& spectra_path, const std::string& database_path,
                                            ProteinIdentificationList& proteins, PeptideIdentificationList& peptides,
                                            const JS::ProSEOptions& options)
    {
      resetProgress(algo);
      bool rescore = false;
      const double protein_fdr = static_cast<double>(algo.getParameters().getValue("FDR:protein"));
      StopWatch watch;
      watch.start();
      ProSEAlgorithm::MultiFileSearchResult result;
      {
        // {rescore} overrides the parameter 'rescore' for the search
        const RescoreOverride rescore_override(algo, options);
        rescore = rescore_override.rescore();
        result = algo.searchWithModificationAnalysis(std::vector<std::string>{spectra_path}, database_path, {}, "", false);
      }
      if (result.per_file.size() != 1)
      {
        throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "ProSE returned no result for the input file.");
      }
      ProSEAlgorithm::SearchResult& run = result.per_file[0];
      const bool ok = run.exit_code == ProSEAlgorithm::ExitCodes::EXECUTION_OK;

      // the outcome of the search's Percolator rescoring
      val percolator = val::null();
      if (rescore && ok)
      {
        percolator = val::object();
        percolator.set("applied", run.percolator.applied);
        percolator.set("skipped", run.percolator.applied ? val::null() : val(run.percolator.skipped_reason));
      }

      bool protein_fdr_applied = false;
      if (ok && protein_fdr > 0.0 && result.have_decoys && !run.protein_ids.empty() && !run.peptide_ids.empty())
      {
        ProSEAlgorithm::applyCompleteSetProteinFDR(run.protein_ids, run.peptide_ids, result.decoy_string,
                                                   result.decoy_is_prefix, protein_fdr);
        protein_fdr_applied = true;
      }
      // The protein FDR changed the PSMs after the search collected its statistics: recount them
      // from the final PSMs, as the ProSE tool does.
      if (ok)
      {
        ProSEAlgorithm::updateFinalStats(run.stats, run.peptide_ids, result.shared.enzyme, run.stats.fdr_applied);
      }
      watch.stop();
      result.shared.seconds_total = watch.getClockTime();

      proteins = std::move(run.protein_ids);
      peptides = std::move(run.peptide_ids);

      val o = val::object();
      o.set("status", exitCodeName(run.exit_code));
      o.set("openSearch", run.is_open_search);
      o.set("haveDecoys", result.have_decoys);
      o.set("decoyString", result.decoy_string);
      o.set("decoyIsPrefix", result.decoy_is_prefix);
      o.set("proteinFDRApplied", protein_fdr_applied);
      o.set("percolator", percolator);
      o.set("run", runStatisticsToJS(run.stats));
      o.set("search", sharedStatisticsToJS(result.shared));
      o.set("modificationAnalysis", run.is_open_search ? modificationAnalysisToJS(run.modification_analysis) : val::null());
      return o;
    }

    JS::ProSESummary proseSearch(ProSEAlgorithm& algo, const std::string& spectra_path, const std::string& database_path,
                                 ProteinIdentificationList& proteins, PeptideIdentificationList& peptides)
    {
      return proseSearchWithOptions(algo, spectra_path, database_path, proteins, peptides, val::undefined());
    }

    void proseSetProgressCallback(ProSEAlgorithm& algo, const JS::ProgressCallback& callback)
    {
      setProgressCallback(algo, callback);
    }
  } // namespace

  //------------------------------------------------------------------------------
  // Several files: searchFiles() and its ProSEResultSet

  JS::ProSEMultiSummary ProSEResultSet::summary() const
  {
    const auto file_summary = [this](const ProSEFileResult& f) {
      const bool ok = f.exit_code == ProSEAlgorithm::ExitCodes::EXECUTION_OK;
      val percolator = val::null();
      if (rescore && ok)
      {
        percolator = val::object();
        percolator.set("applied", f.percolator.applied);
        percolator.set("skipped", f.percolator.applied ? val::null() : val(f.percolator.skipped_reason));
        percolator.set("warning", f.percolator.warning.empty() ? val::null() : val(f.percolator.warning));
      }
      val o = val::object();
      o.set("path", f.path);
      o.set("status", exitCodeName(f.exit_code));
      o.set("ms2Spectra", static_cast<double>(f.stats.ms2_spectra));
      o.set("targetPSMs", static_cast<double>(f.stats.target_psms));
      o.set("decoyPSMs", static_cast<double>(f.stats.decoy_psms));
      o.set("percolator", percolator);
      o.set("peptdeepNce", f.peptdeep_nce ? val(*f.peptdeep_nce) : val::null());
      o.set("secondsSearch", f.stats.seconds_search);
      return o;
    };
    val o = val::object();
    o.set("files", toArray(files, file_summary));
    o.set("proteinFdrApplied", merged.protein_fdr_applied);
    o.set("psmFdr", psm_fdr);
    o.set("proteinFdr", protein_fdr);
    o.set("decoyString", decoy_string);
    o.set("openSearch", open_search);
    o.set("modificationAnalysis", open_search ? modificationAnalysisToJS(modification_analysis) : val::null());
    return o;
  }

  void ProSEResultSet::storeIdXML(const std::string& path) const
  {
    IdXMLFile().store(path, merged.protein_ids, merged.peptide_ids);
  }

  namespace
  {
    /// Searches the spectra of @p files (mzML) against the proteins of @p database_path (FASTA) as the ProSE tool does
    /// with -out_merged: one search of all files, which by default (parameter 'rescore') rescores the PSMs of each file
    /// with Percolator before its FDR:PSM; then the PSMs of all files are pooled for the protein inference and, with
    /// FDR:protein and decoys, the picked-protein FDR (ProSEAlgorithm::mergeFileResults()). Options: {rescore} sets the
    /// parameter 'rescore' for the call; {removeDecoys} (default false) removes the decoys after the protein FDR, as
    /// the tool does; {keepSpectra} (default true) keeps the peaks of the identified spectra, which are read from the
    /// files at the end (ProSEResultSet::collectSpectra()).
    std::unique_ptr<ProSEResultSet> proseSearchFilesWithOptions(ProSEAlgorithm& algo, const JS::StringArray& files_js,
                                                                const std::string& database_path,
                                                                const JS::ProSEFilesOptions& options)
    {
      if (!val::global("Array").call<bool>("isArray", files_js))
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "searchFiles() takes the paths of the spectrum files as array.");
      }
      if (!isNullish(options) && options.typeOf().as<std::string>() != "object")
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The options have to be an object (or null).");
      }
      const std::vector<std::string> files = stringsFromJS(files_js);
      if (files.empty())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "searchFiles() needs at least one spectrum file.");
      }
      // The pooled PSMs refer to their file by its path ('id_merge_index'): the same path twice would merge two files.
      for (Size i = 1; i < files.size(); ++i)
      {
        if (std::find(files.begin(), files.begin() + i, files[i]) != files.begin() + i)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "searchFiles() needs distinct paths: '" + files[i] + "' is given twice.");
        }
      }
      const val remove_decoys_option = option(options, "removeDecoys");
      const val keep_spectra_option = option(options, "keepSpectra");
      const bool remove_decoys = !isNullish(remove_decoys_option) && remove_decoys_option.as<bool>();

      resetProgress(algo);
      auto set = std::make_unique<ProSEResultSet>();
      set->keep_spectra = isNullish(keep_spectra_option) || keep_spectra_option.as<bool>();
      set->psm_fdr = static_cast<double>(algo.getParameters().getValue("FDR:PSM"));
      set->protein_fdr = static_cast<double>(algo.getParameters().getValue("FDR:protein"));
      FASTAFile().load(database_path, set->database);

      // {rescore} overrides the parameter 'rescore' for the search and its evaluation
      const RescoreOverride rescore_override(algo, options);
      set->rescore = rescore_override.rescore();
      // No pooled copy of the PSMs (build_pooled_aggregate false): mergeFileResults() pools them by moving them out of
      // the results of the files. An open search still analyses the modifications of the pooled PSMs.
      ProSEAlgorithm::MultiFileSearchResult result = algo.searchWithModificationAnalysis(files, set->database, {}, "", false);
      if (result.per_file.size() != files.size())
      {
        throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "ProSE returned no result for some input files.");
      }

      std::string not_rescored; // the files Percolator did not rescore, with the reasons
      bool have_psms = false;   // whether a file has a PSM
      for (Size i = 0; i < files.size(); ++i)
      {
        ProSEAlgorithm::SearchResult& run = result.per_file[i];
        ProSEFileResult f;
        f.path = files[i];
        f.exit_code = run.exit_code;
        f.percolator = run.percolator;
        if (run.exit_code == ProSEAlgorithm::ExitCodes::EXECUTION_OK)
        {
          if (!run.protein_ids.empty())
          {
            ProteinIdentification& protein_id = run.protein_ids[0];
            protein_id.getSearchParameters().db = database_path; // as a search of the FASTA file records it
            if (protein_id.metaValueExists(PeptDeepRescoring::NCE_META_VALUE))
            {
              f.peptdeep_nce = static_cast<double>(protein_id.getMetaValue(PeptDeepRescoring::NCE_META_VALUE));
            }
            // the features of the run for tables(): the pooled run records those of one file only
            const ProteinIdentification::SearchParameters& sp = protein_id.getSearchParameters();
            if (sp.metaValueExists("extra_features"))
            {
              for (const std::string& name : ListUtils::create<std::string>(sp.getMetaValue("extra_features").toString()))
              {
                if (std::find(set->extra_features.begin(), set->extra_features.end(), name) == set->extra_features.end())
                {
                  set->extra_features.push_back(name);
                }
              }
            }
          }
          // the statistics of the PSMs of the file, as the ProSE tool reports them
          ProSEAlgorithm::updateFinalStats(run.stats, run.peptide_ids, result.shared.enzyme, run.stats.fdr_applied);
          have_psms = have_psms || std::any_of(run.peptide_ids.begin(), run.peptide_ids.end(),
                                               [](const PeptideIdentification& id) { return !id.getHits().empty(); });
          if (set->rescore && !run.percolator.applied)
          {
            not_rescored += (not_rescored.empty() ? "" : "; ") + files[i] + " (" + run.percolator.skipped_reason + ")";
          }
        }
        f.stats = run.stats;
        set->files.push_back(std::move(f));
      }
      set->open_search = result.aggregate.is_open_search;
      // a search of one file analyses the modifications of its PSMs, not of a pooled set
      set->modification_analysis = std::move(files.size() == 1 ? result.per_file[0].modification_analysis
                                                                : result.aggregate.modification_analysis);
      set->decoy_string = result.decoy_string;
      set->decoy_is_prefix = result.decoy_is_prefix;
      // how the search built its decoys, for the sequences of the decoy proteins (ProSEAlgorithm's decoy database)
      {
        const Param& p = algo.getParameters();
        set->decoys.generated = result.shared.decoy_mode == "generated";
        set->decoys.enzyme = p.getValue("enzyme").toString();
        set->decoys.unspecific = p.getValue("peptide:enzyme_specificity").toString() == "none";
        set->decoys.keep_initial_met = p.getValue("peptide:clip_nterm_methionine").toBool();
        set->min_charge = static_cast<int>(p.getValue("precursor:min_charge"));
        set->max_charge = static_cast<int>(p.getValue("precursor:max_charge"));
      }

      // Without PSMs, nothing is pooled: the results stay empty, as without a searched file. (mergeFileResults() would
      // throw: the protein inference needs PSMs.)
      if (!have_psms) return set;
      bool partly_rescored = false;
      try
      {
        set->merged = algo.mergeFileResults(result, files, remove_decoys);
      }
      catch (const Exception::InvalidValue&) // Percolator rescored only some of the files
      {
        partly_rescored = true;
      }
      if (partly_rescored)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Percolator rescored only some of the input files, and their PSMs cannot be pooled with the search engine "
          "scores of the others (not rescored: " + not_rescored + "). Search the files with { rescore: false }, or leave "
          "out those that could not be rescored.");
      }
      // the peaks of the identified spectra (spectra()), read now, while the input files are there
      set->collectSpectra(&algo);
      return set;
    }

    std::unique_ptr<ProSEResultSet> proseSearchFiles(ProSEAlgorithm& algo, const JS::StringArray& files, const std::string& database_path)
    {
      return proseSearchFilesWithOptions(algo, files, database_path, val::undefined());
    }

    //------------------------------------------------------------------------------
    // Annotation of a spectrum with the fragment ions of a peptide

    /// The fragment ions of a peptide as TheoreticalSpectrumGenerator makes them, sorted by m/z, with the data arrays of
    /// their names and charges (add_metainfo): what annotateSpectrum(), annotatePeaks() and predictMS2() share.
    struct TheoreticalFragments
    {
      MSSpectrum spectrum;       ///< the fragment ions sorted by m/z, with their names and charges
      Size names_array = 0;      ///< the data array of spectrum with the ion names
      Size charges_array = 0;    ///< the data array of spectrum with the ion charges

      const MSSpectrum::StringDataArray& names() const { return spectrum.getStringDataArrays()[names_array]; }
      const MSSpectrum::IntegerDataArray& charges() const { return spectrum.getIntegerDataArrays()[charges_array]; }
    };

    /// The fragment ions of @p peptide of the charges 1 to @p max_charge (at least 1) and the types @p ions ("a", "b", "c",
    /// "x", "y", "z"), with their neutral losses if @p losses and the precursor peaks (of the charge max_charge) if
    /// @p precursor. The first prefix ion is there (b1).
    TheoreticalFragments theoreticalFragments(const AASequence& peptide, const std::vector<std::string>& ions, int max_charge,
                                              bool losses, bool precursor)
    {
      TheoreticalSpectrumGenerator generator;
      Param p = generator.getParameters();
      for (const char* ion : {"a", "b", "c", "x", "y", "z"})
      {
        const bool wanted = std::find(ions.begin(), ions.end(), ion) != ions.end();
        p.setValue(std::string("add_") + ion + "_ions", wanted ? "true" : "false");
      }
      for (const auto& ion : ions)
      {
        if (ion.size() != 1 || std::string("abcxyz").find(ion) == std::string::npos)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown ion type (a, b, c, x, y, z)", ion);
        }
      }
      p.setValue("add_losses", losses ? "true" : "false");
      p.setValue("add_precursor_peaks", precursor ? "true" : "false");
      p.setValue("add_first_prefix_ion", "true");
      p.setValue("add_metainfo", "true");
      generator.setParameters(p);

      TheoreticalFragments result;
      generator.getSpectrum(result.spectrum, peptide, 1, std::max(1, max_charge));
      result.spectrum.sortByPosition(); // keeps the data arrays aligned

      // the arrays add_metainfo adds: ion names and charges
      MSSpectrum& spectrum = result.spectrum;
      const auto names_it = std::find_if(spectrum.getStringDataArrays().begin(), spectrum.getStringDataArrays().end(),
        [](const auto& a) { return a.getName() == Constants::UserParam::IonNames; });
      const auto charges_it = std::find_if(spectrum.getIntegerDataArrays().begin(), spectrum.getIntegerDataArrays().end(),
        [](const auto& a) { return a.getName() == "Charges"; });
      if (names_it == spectrum.getStringDataArrays().end() || charges_it == spectrum.getIntegerDataArrays().end())
      {
        throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The theoretical spectrum has no ion names or charges.");
      }
      result.names_array = static_cast<Size>(names_it - spectrum.getStringDataArrays().begin());
      result.charges_array = static_cast<Size>(charges_it - spectrum.getIntegerDataArrays().begin());
      return result;
    }

    /// The fragment ions of a peptide matched to the peaks of a spectrum: what annotateSpectrum() and annotatePeaks()
    /// share.
    struct FragmentAlignment
    {
      std::string unit;                             ///< "Da" or "ppm": the unit of the tolerance
      MSSpectrum theoretical;                       ///< the fragment ions sorted by m/z, with their names and charges
      MSSpectrum measured;                          ///< the peaks of the spectrum sorted by m/z
      std::vector<Size> measured_position;          ///< measured[k] is the peak measured_position[k] of the spectrum as given
      std::vector<std::pair<Size, Size>> alignment; ///< (theoretical, measured) index pairs by ascending m/z; see alignFragments()
      Size names_array = 0;                         ///< the data array of theoretical with the ion names
      Size charges_array = 0;                       ///< the data array of theoretical with the ion charges

      const MSSpectrum::StringDataArray& names() const { return theoretical.getStringDataArrays()[names_array]; }
      const MSSpectrum::IntegerDataArray& charges() const { return theoretical.getIntegerDataArrays()[charges_array]; }
    };

    /// Matches the fragment ions of @p peptide (charges 1 to maxCharge) to the peaks of @p spectrum, which need not be
    /// sorted. Options: tolerance (default 0.02), unit ("Da" or "ppm", default "Da"), maxCharge (default: precursor
    /// charge - 1, at least 1), ions (default ["b", "y"]; also "a", "c", "x", "z"), losses (neutral losses, default
    /// false), precursor (precursor peaks, of the charge maxCharge, default false).
    ///
    /// How the ions take peaks depends on the unit (SpectrumAlignment):
    /// - ppm: each ion takes the closest peak within the tolerance (MatchedIterator). Several ions can take the same
    ///   peak (e.g. z1++ and y1-NH3++, which have the same m/z), so a peak can be in several alignment pairs. Of peaks
    ///   of equal m/z the one first in the spectrum is taken (the sort is stable).
    /// - Da: a dynamic-programming alignment, one to one and in m/z order: a peak has at most one ion, and of two
    ///   ions of the same m/z only one has a peak (unless there are two peaks). Of peaks of equal m/z one is taken;
    ///   which one is unspecified (it depends on floating-point details of the alignment and on the neighbouring peaks).
    FragmentAlignment alignFragments(const MSSpectrum& spectrum, const AASequence& peptide, int precursor_charge,
                                     const JS::AnnotationOptions& options)
    {
      FragmentAlignment result;
      const bool has_options = !isNullish(options);
      const auto option = [&](const char* name) { return has_options ? options[name] : val::undefined(); };

      const double tolerance = isNullish(option("tolerance")) ? 0.02 : option("tolerance").as<double>();
      const std::string unit = isNullish(option("unit")) ? std::string("Da") : option("unit").as<std::string>();
      if (unit != "Da" && unit != "ppm")
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "unit has to be 'Da' or 'ppm'", unit);
      }
      result.unit = unit;
      const int max_charge = isNullish(option("maxCharge")) ? std::max(1, precursor_charge - 1) : option("maxCharge").as<int>();
      std::vector<std::string> ions{"b", "y"};
      if (!isNullish(option("ions"))) ions = stringsFromJS(option("ions"));
      const bool losses = !isNullish(option("losses")) && option("losses").as<bool>();
      const bool precursor = !isNullish(option("precursor")) && option("precursor").as<bool>();

      // the theoretical spectrum, with the names and charges of its ions
      TheoreticalFragments fragments = theoreticalFragments(peptide, ions, max_charge, losses, precursor);
      result.theoretical = std::move(fragments.spectrum);
      result.names_array = fragments.names_array;
      result.charges_array = fragments.charges_array;
      MSSpectrum& theoretical = result.theoretical;

      // the measured peaks, sorted by m/z (SpectrumAlignment needs that), and where they were in the spectrum
      result.measured_position.resize(spectrum.size());
      std::iota(result.measured_position.begin(), result.measured_position.end(), Size(0));
      std::stable_sort(result.measured_position.begin(), result.measured_position.end(),
                       [&spectrum](Size a, Size b) { return spectrum[a].getMZ() < spectrum[b].getMZ(); });
      MSSpectrum& measured = result.measured;
      measured.reserve(spectrum.size());
      for (const Size position : result.measured_position) measured.push_back(spectrum[position]);

      SpectrumAlignment aligner;
      Param ap = aligner.getParameters();
      ap.setValue("tolerance", tolerance);
      ap.setValue("is_relative_tolerance", unit == "ppm" ? "true" : "false");
      aligner.setParameters(ap);
      aligner.getSpectrumAlignment(result.alignment, theoretical, measured);

      return result;
    }

    /// Annotates @p spectrum with the fragment ions of @p peptide: see alignFragments() for the options.
    JS::SpectrumAnnotation annotateSpectrum(const MSSpectrum& spectrum, const AASequence& peptide, int precursor_charge,
                                            const JS::AnnotationOptions& options)
    {
      const FragmentAlignment fragments = alignFragments(spectrum, peptide, precursor_charge, options);
      const MSSpectrum& theoretical = fragments.theoretical;
      const MSSpectrum& measured = fragments.measured;
      const auto& names = fragments.names();
      const auto& charges = fragments.charges();
      const std::string& unit = fragments.unit;

      val matches = val::array();
      for (const auto& [t, m] : fragments.alignment)
      {
        const double theo_mz = theoretical[t].getMZ();
        const double mz = measured[m].getMZ();
        val o = val::object();
        o.set("mz", mz);
        o.set("intensity", measured[m].getIntensity());
        o.set("ion", names[t]);
        o.set("charge", charges[t]);
        o.set("theoreticalMZ", theo_mz);
        o.set("error", unit == "ppm" ? (mz - theo_mz) / theo_mz * 1e6 : mz - theo_mz);
        matches.call<void>("push", o);
      }

      val o = val::object();
      o.set("matches", matches);
      val theo = val::object();
      theo.set("mz", toTypedArray<double>(theoretical, [](const Peak1D& peak) { return peak.getMZ(); }));
      theo.set("ion", toArray(names));
      theo.set("charge", toArray(charges));
      o.set("theoretical", theo);
      o.set("unit", unit);
      return o;
    }

    /// A fragment ion as TheoreticalSpectrumGenerator names it, split up.
    struct IonLabel
    {
      std::string ion;  ///< the name without the '+' of charge 1: "b3", "y5++", "y7-H2O", "[M+2H]++"
      std::string type; ///< "a" to "z", or "precursor"
      int ordinal = 0;  ///< the ion number: 3 for "b3"; 0 for a precursor
      std::string loss; ///< the formula of the neutral loss: "H2O" for "y7-H2O"; "" for none
    };

    /// A neutral loss formula as written in ion names: without the count 1 of an element ("H2O" for the "H2O1" that
    /// TheoreticalSpectrumGenerator makes of an EmpiricalFormula), and the common losses as they are conventionally
    /// spelled, not in the formula's alphabetical order: "NH3" (ammonia; not "H3N"), "H3PO4" (not "H3O4P") and "HPO3"
    /// (not "HO3P"). Others stay in that order, e.g. "CH4OS" (the loss of oxidized methionine).
    std::string lossFormulaName(const std::string& formula)
    {
      std::string name;
      for (size_t i = 0; i < formula.size(); ++i)
      {
        const bool count_one = formula[i] == '1' && i > 0 && std::isalpha(static_cast<unsigned char>(formula[i - 1]))
                               && (i + 1 == formula.size() || !std::isdigit(static_cast<unsigned char>(formula[i + 1])));
        if (!count_one) name += formula[i];
      }
      static const std::unordered_map<std::string, std::string> conventional{{"H3N", "NH3"}, {"H3O4P", "H3PO4"}, {"HO3P", "HPO3"}};
      const auto it = conventional.find(name);
      return it == conventional.end() ? name : it->second;
    }

    /// Splits the name of a fragment ion that TheoreticalSpectrumGenerator gives a peak of @p charge (with add_metainfo).
    /// It names a fragment ion by its type, number, neutral loss if any (as the formula "H2O1", "H3N1", "CH4OS1"), and
    /// one '+' per charge: "b3+", "y5++", "y7-H2O1++"; a precursor peak by its adduct, with the loss inside: "[M+H]+",
    /// "[M+2H-NH3]++". The label has the loss as lossFormulaName() and keeps the '+' signs of charge 2 and more only:
    /// "b3", "y5++", "y7-H2O++".
    IonLabel parseIonName(const std::string& name, int charge)
    {
      const auto unexpected = [&name]() {
        return Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unexpected name of a fragment ion", name);
      };
      // without the charge: at most one '+' per charge (a name without charge, as of isotope peaks, has none)
      std::string base = name;
      for (int z = 0; z < charge && !base.empty() && base.back() == '+'; ++z) base.pop_back();

      IonLabel label;
      if (base.size() > 4 && base.compare(0, 3, "[M+") == 0 && base.back() == ']') // "[M+H]", "[M+2H-H2O]"
      {
        label.type = "precursor";
        const size_t minus = base.find('-');
        if (minus != std::string::npos)
        {
          label.loss = lossFormulaName(base.substr(minus + 1, base.size() - minus - 2)); // up to the ']'
          if (label.loss.empty()) throw unexpected();
          base = base.substr(0, minus + 1) + label.loss + "]";
        }
      }
      else
      {
        // the type (a letter, but z. and z' too), the number, and "-" with the formula of the loss
        size_t digits = 0;
        while (digits < base.size() && !std::isdigit(static_cast<unsigned char>(base[digits])) && base[digits] != '-') ++digits;
        size_t end = digits;
        while (end < base.size() && std::isdigit(static_cast<unsigned char>(base[end]))) ++end;
        if (digits == 0 || end == digits) throw unexpected();
        label.type = base.substr(0, digits);
        if (std::from_chars(base.data() + digits, base.data() + end, label.ordinal).ec != std::errc()) throw unexpected();
        if (end < base.size())
        {
          if (base[end] != '-' || end + 1 == base.size()) throw unexpected();
          label.loss = lossFormulaName(base.substr(end + 1));
          base = base.substr(0, end + 1) + label.loss;
        }
      }
      label.ion = charge > 1 ? base + std::string(static_cast<size_t>(charge), '+') : base;
      return label;
    }

    /// Annotates the peaks mz/intensity (arrays of one spectrum, in any order; same length) with the fragment ions of
    /// @p peptide_string (AASequence.toString(), e.g. "AM(Oxidation)K") of precursor charge @p charge: see
    /// alignFragments() for the options. The ions take peaks as in annotateSpectrum(): in ppm each ion takes the closest
    /// peak, so several ions can share one (a peak is then in several matches); in Da the matching is one to one. A
    /// match refers to its peak by its index in the arrays as passed; its error is always in ppm.
    JS::PeakAnnotationResult annotatePeaks(const JS::NumberInput& mz, const JS::NumberInput& intensity,
                                           const std::string& peptide_string, int charge, const JS::AnnotationOptions& options)
    {
      if (charge < 1)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The precursor charge has to be at least 1", std::to_string(charge));
      }
      const AASequence peptide = AASequence::fromString(peptide_string);
      if (peptide.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The peptide is empty", peptide_string);
      }
      // float m/z become doubles exactly
      const auto mzs = numbersFromJS<double>(mz);
      const auto intensities = numbersFromJS<float>(intensity);
      if (mzs.size() != intensities.size())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "m/z and intensity arrays differ in length (" + std::to_string(mzs.size()) + " vs. " + std::to_string(intensities.size()) + ")");
      }
      MSSpectrum spectrum;
      spectrum.reserve(mzs.size());
      for (Size i = 0; i < mzs.size(); ++i)
      {
        if (!std::isfinite(mzs[i]))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "m/z has to be a finite number", "peak " + std::to_string(i));
        }
        spectrum.emplace_back(mzs[i], intensities[i]);
      }

      const FragmentAlignment fragments = alignFragments(spectrum, peptide, charge, options);
      const auto& names = fragments.names();
      const auto& charges = fragments.charges();
      std::vector<IonLabel> labels;
      labels.reserve(fragments.theoretical.size());
      for (Size t = 0; t < fragments.theoretical.size(); ++t) labels.push_back(parseIonName(names[t], charges[t]));

      val matches = val::array();
      for (const auto& [t, m] : fragments.alignment)
      {
        const double theo_mz = fragments.theoretical[t].getMZ();
        const double observed = fragments.measured[m].getMZ();
        val o = val::object();
        o.set("peak", static_cast<double>(fragments.measured_position[m]));
        o.set("ion", labels[t].ion);
        o.set("type", labels[t].type);
        o.set("ordinal", labels[t].ordinal);
        o.set("charge", charges[t]);
        o.set("loss", labels[t].loss);
        o.set("theoreticalMz", theo_mz);
        o.set("errorPpm", (observed - theo_mz) / theo_mz * 1e6);
        matches.call<void>("push", o);
      }
      val theoretical = val::array();
      for (Size t = 0; t < fragments.theoretical.size(); ++t)
      {
        val o = val::object();
        o.set("ion", labels[t].ion);
        o.set("type", labels[t].type);
        o.set("ordinal", labels[t].ordinal);
        o.set("charge", charges[t]);
        o.set("mz", fragments.theoretical[t].getMZ());
        theoretical.call<void>("push", o);
      }
      val result = val::object();
      result.set("matches", matches);
      result.set("theoretical", theoretical);
      return result;
    }

    //------------------------------------------------------------------------------
    // Fragment intensities predicted by PeptDeep

    /// The MS2 model of predictMS2() unless another is named: the one of the package (models/ next to openms.mjs), as the
    /// application mounts it.
    constexpr const char* default_ms2_model = "/models/peptdeep_ms2_dynamic.onnx";

#ifdef WITH_ONNX
    /// The inference session of the MS2 model file @p model_path. The first call for a path loads the model (it takes a
    /// while); the session is then kept for the rest of the module's life and used by every later call, so the file is
    /// read once, whatever happens to it (or to its mount point) afterwards. A model that does not load throws and is not
    /// kept, so that a later call can load it (e.g. once the file is mounted).
    PeptDeepMS2Inference& ms2Inference(const std::string& model_path)
    {
      static std::mutex mutex;
      static std::unordered_map<std::string, std::unique_ptr<PeptDeepMS2Inference>> sessions;
      const std::lock_guard<std::mutex> lock(mutex);
      const auto it = sessions.find(model_path);
      if (it != sessions.end()) return *it->second;
      // The threads as ProSEAlgorithm gives them to PeptDeepRescoring (the threads OpenMP computes with), and the
      // default batch size (that of the 'batch_size' parameter, too).
#ifdef _OPENMP
      const int threads = std::max(1, omp_get_max_threads());
#else
      const int threads = 1;
#endif
      return *sessions.emplace(model_path, std::make_unique<PeptDeepMS2Inference>(model_path, threads)).first->second;
    }
#endif

    /// The fragment ions of the peptide @p peptide_string (AASequence.toString(), e.g. "AVM(Oxidation)DDFAAFVEK", with all
    /// modifications written out) of precursor charge @p charge, with the intensities that PeptDeep (the MS2 model at
    /// @p model_path) predicts for collision energy @p nce (PeptDeep's normalised collision energy, as the 'nce' of
    /// PeptDeepRescoring) on @p instrument ("QE", "Lumos", "timsTOF", "SciexTOF", "ThermoTOF", case-insensitively; any
    /// other name is peptdeep's catch-all instrument, ML::instrumentIndex()).
    ///
    /// - Only the b and y ions of charge 1 and 2 without neutral loss (not the modloss channels of the model); charge 2
    ///   only if the precursor has charge 2 or more.
    /// - Intensities relative to the most intense of these ions (which is 1), the ions of 0.01 or less left out; no ions
    ///   at all if the model predicts none.
    /// - m/z from TheoreticalSpectrumGenerator, as annotatePeaks() (b and y, charges 1 and 2).
    /// - By ascending m/z; a peptide of one residue has no ions.
    /// Without ONNX Runtime (WITH_ONNX) in the module it throws.
    JS::PredictedIons predictMS2(const std::string& peptide_string, int charge, double nce, const std::string& instrument,
                                 const std::string& model_path)
    {
#ifdef WITH_ONNX
      if (charge < 1)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The precursor charge has to be at least 1", std::to_string(charge));
      }
      if (!std::isfinite(nce))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The collision energy has to be a finite number", std::to_string(nce));
      }
      const AASequence peptide = AASequence::fromString(peptide_string);
      if (peptide.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The peptide is empty", peptide_string);
      }

      // As PeptDeepRescoring: the peptide as AASequence::toString() writes it, the collision energy as it is (not /100),
      // and the instrument as peptdeep's index.
      PeptDeepMS2Inference& model = ms2Inference(model_path);
      val ions = val::array();
      const Size length = peptide.size();
      if (length < 2) return ions; // no cleavage
      const std::vector<std::vector<float>> predicted =
        model.predictMS2({peptide.toString()}, {static_cast<float>(charge)}, {static_cast<float>(nce)}, {ML::instrumentIndex(instrument)});
      const std::vector<float>& flat = predicted.at(0);

      // The model gives eight channels for each of the length - 1 cleavages: [b_z1, b_z2, y_z1, y_z2, b_modloss_z1, ...],
      // of which the first four are the b and y ions. Cleavage `pos` (from the N-terminus) carries b(pos + 1) and
      // y(length - 1 - pos), as the model lists the fragments N to C: the y series is numbered from the other end.
      // Keep in sync with PeptDeepRescoring::predictedIons_() (PeptDeepRescoring.cpp: a helper in an anonymous namespace,
      // so not callable from here), which PeptDeepLibraryPredictor::predict() agrees with.
      constexpr int ion_b1 = 0, ion_b2 = 1, ion_y1 = 2, ion_y2 = 3;
      const Size channels = flat.size() / (length - 1);
      if (channels <= static_cast<Size>(ion_y2) || flat.size() != channels * (length - 1))
      {
        throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "PeptDeep gave " + std::to_string(flat.size()) + " intensities for the " + std::to_string(length - 1) + " cleavages of " + peptide_string);
      }
      struct Predicted
      {
        char type;
        int ordinal;
        int charge;
        double intensity;
      };
      std::vector<Predicted> candidates;
      candidates.reserve(4 * (length - 1));
      const bool doubly_charged = charge > 1; // a precursor of charge 1 has no doubly charged fragments
      for (Size pos = 0; pos + 1 < length; ++pos)
      {
        const Size base = pos * channels;
        const int b = static_cast<int>(pos + 1), y = static_cast<int>(length - 1 - pos);
        candidates.push_back({'b', b, 1, flat[base + ion_b1]});
        candidates.push_back({'y', y, 1, flat[base + ion_y1]});
        if (doubly_charged)
        {
          candidates.push_back({'b', b, 2, flat[base + ion_b2]});
          candidates.push_back({'y', y, 2, flat[base + ion_y2]});
        }
      }
      // relative to the most intense of these ions; the ions of 0.01 or less are left out
      double base_peak = 0;
      for (const Predicted& c : candidates) base_peak = std::max(base_peak, c.intensity);
      if (!(base_peak > 0)) return ions; // nothing predicted
      std::map<std::tuple<char, int, int>, double> intensity;
      for (const Predicted& c : candidates)
      {
        const double relative = c.intensity / base_peak;
        if (relative > 0.01) intensity[{c.type, c.ordinal, c.charge}] = relative;
      }

      // the ions with their m/z, by ascending m/z, named as annotatePeaks() names them
      const TheoreticalFragments theoretical = theoreticalFragments(peptide, {"b", "y"}, doubly_charged ? 2 : 1, false, false);
      const auto& names = theoretical.names();
      const auto& charges = theoretical.charges();
      Size found = 0;
      for (Size t = 0; t < theoretical.spectrum.size(); ++t)
      {
        const IonLabel label = parseIonName(names[t], charges[t]);
        const auto it = intensity.find({label.type.at(0), label.ordinal, charges[t]});
        if (it == intensity.end()) continue;
        ++found;
        val o = val::object();
        o.set("ion", label.ion);
        o.set("type", label.type);
        o.set("ordinal", label.ordinal);
        o.set("charge", charges[t]);
        o.set("mz", theoretical.spectrum[t].getMZ());
        o.set("intensity", it->second);
        ions.call<void>("push", o);
      }
      if (found != intensity.size())
      {
        throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "The theoretical spectrum of " + peptide_string + " lacks some of the predicted ions");
      }
      return ions;
#else
      (void)peptide_string; (void)charge; (void)nce; (void)instrument; (void)model_path;
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "predictMS2: this OpenMS module was built without ONNX Runtime");
#endif
    }

    /// predictMS2() with the default model
    JS::PredictedIons predictMS2WithDefaultModel(const std::string& peptide_string, int charge, double nce, const std::string& instrument)
    {
      return predictMS2(peptide_string, charge, nce, instrument, default_ms2_model);
    }
  } // namespace
} // namespace OpenMS::WASM

using namespace OpenMS;
using namespace OpenMS::WASM;

EMSCRIPTEN_BINDINGS(openms_prose)
{
  // The results of searchFiles(); delete() releases them.
  class_<ProSEResultSet>("ProSEResultSet")
    .function("summary", OMS_METHOD(ProSEResultSet, &ProSEResultSet::summary))
    .function("tables", OMS_METHOD(ProSEResultSet, &ProSEResultSet::tables))
    .function("spectra", OMS_METHOD(ProSEResultSet, &ProSEResultSet::spectra))
    .function("storeIdXML", OMS_METHOD(ProSEResultSet, &ProSEResultSet::storeIdXML))
    // for tests: (row, nativeId) replaces the native ID of a PSM row and reads the spectra from the input files again
    .function("_debugSetNativeId", OMS_METHOD(ProSEResultSet, &ProSEResultSet::debugSetNativeId))
    ;

  class_<ProSEAlgorithm, base<DefaultParamHandler>>("ProSEAlgorithm")
    .constructor<>()
    .function("search", OMS_FN(&proseSearch))
    .function("search", OMS_FN(&proseSearchWithOptions))
    .function("searchFiles", OMS_FN(&proseSearchFiles), nonnull<ret_val>())
    .function("searchFiles", OMS_FN(&proseSearchFilesWithOptions), nonnull<ret_val>())
    .function("setProgressCallback", OMS_FN(&proseSetProgressCallback))
    // q-values of the best hits by target-decoy competition (as FDR:PSM, without filtering)
    .function("annotatePsmQValues", OMS_METHOD(ProSEAlgorithm, &ProSEAlgorithm::annotatePsmQValues))
    ;

  function("annotateSpectrum", OMS_FN(&annotateSpectrum));
  function("annotatePeaks", OMS_FN(&annotatePeaks));
  function("predictMS2", OMS_FN(&predictMS2WithDefaultModel));
  function("predictMS2", OMS_FN(&predictMS2));
}
