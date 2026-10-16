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
// does, whose results (ProSEResultSet) JavaScript holds; and the annotation of a spectrum with
// the fragment ions of a peptide, for spectrum viewers.

#include "binding_utils.h"

#include <OpenMS/ANALYSIS/ID/PeptDeepRescoring.h>
#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/StopWatch.h>

#include <pthread.h>

#include <algorithm>
#include <atomic>
#include <charconv>
#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
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
    class JSProgressLogger : public ProgressLogger::ProgressLoggerImpl
    {
    public:
      explicit JSProgressLogger(val callback) :
        callback_(std::move(callback)),
        owner_(pthread_self())
      {
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
      logger.setLogger(new JSProgressLogger(callback)); // the logger owns it
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

  /// The search of one input file of ProSEAlgorithm.searchFiles()
  struct ProSEFileResult
  {
    std::string path;                                   ///< as given
    ProSEAlgorithm::ExitCodes exit_code = ProSEAlgorithm::ExitCodes::EXECUTION_OK;
    ProSEAlgorithm::RunStatistics stats;                ///< of the PSMs of the file (after its FDR:PSM)
    ProSEAlgorithm::PercolatorRescoring percolator;
    std::optional<double> peptdeep_nce;                 ///< the NCE of PeptDeep's predictions (PeptDeepRescoring::NCE_META_VALUE)
  };

  /// The results of ProSEAlgorithm.searchFiles(), which JavaScript holds until delete(): the search of each input
  /// file and the pooled evaluation of all of them (ProSEAlgorithm::mergeFileResults(), the ProSE tool's -out_merged).
  struct ProSEResultSet
  {
    std::vector<ProSEFileResult> files;      ///< in input order
    /// the pooled run and the PSMs of all files; the meta value 'id_merge_index' of a PSM is the position of its file
    ProSEAlgorithm::MergedResult merged;
    std::vector<FASTAFile::FASTAEntry> database; ///< the proteins of the FASTA file searched
    bool rescore = false;                    ///< whether the search rescored the PSMs with Percolator ('rescore')
    bool open_search = false;
    /// of the PSMs of all files (open search)
    OpenSearchModificationAnalysis::OpenSearchAnalysisResult modification_analysis;
    std::string decoy_string;                ///< the decoy marker of the accessions (empty: target-only)
    double psm_fdr = 0.0;                    ///< FDR:PSM
    double protein_fdr = 0.0;                ///< FDR:protein
    bool keep_spectra = true;                ///< ProSEFilesOptions.keepSpectra

    ProSEResultSet() = default;
    // no copies of the results by mistake (Embind copies a class returned by value)
    ProSEResultSet(const ProSEResultSet&) = delete;
    ProSEResultSet& operator=(const ProSEResultSet&) = delete;

    /// The summary of the search (ProSEMultiSummary)
    JS::ProSEMultiSummary summary() const;

    /// Stores the pooled run and PSMs as idXML (the ProSE tool's -out_merged)
    void storeIdXML(const std::string& path) const;
  };

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
    /// the tool does; {keepSpectra} (default true) keeps the peaks of the identified spectra.
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
      return set;
    }

    std::unique_ptr<ProSEResultSet> proseSearchFiles(ProSEAlgorithm& algo, const JS::StringArray& files, const std::string& database_path)
    {
      return proseSearchFilesWithOptions(algo, files, database_path, val::undefined());
    }

    //------------------------------------------------------------------------------
    // Annotation of a spectrum with the fragment ions of a peptide

    /// Matches the fragment ions of @p peptide (charges 1 to maxCharge) to the peaks of
    /// @p spectrum, which need not be sorted. Options: tolerance (default 0.02), unit ("Da"
    /// or "ppm", default "Da"), maxCharge (default: precursor charge - 1, at least 1), ions
    /// (default ["b", "y"]; also "a", "c", "x", "z"), losses (neutral losses, default false),
    /// precursor (precursor peaks, default false).
    JS::SpectrumAnnotation annotateSpectrum(const MSSpectrum& spectrum, const AASequence& peptide, int precursor_charge,
                                            const JS::AnnotationOptions& options)
    {
      const bool has_options = !isNullish(options);
      const auto option = [&](const char* name) { return has_options ? options[name] : val::undefined(); };

      const double tolerance = isNullish(option("tolerance")) ? 0.02 : option("tolerance").as<double>();
      const std::string unit = isNullish(option("unit")) ? std::string("Da") : option("unit").as<std::string>();
      if (unit != "Da" && unit != "ppm")
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "unit has to be 'Da' or 'ppm'", unit);
      }
      const int max_charge = isNullish(option("maxCharge")) ? std::max(1, precursor_charge - 1) : option("maxCharge").as<int>();
      std::vector<std::string> ions{"b", "y"};
      if (!isNullish(option("ions"))) ions = stringsFromJS(option("ions"));
      const bool losses = !isNullish(option("losses")) && option("losses").as<bool>();
      const bool precursor = !isNullish(option("precursor")) && option("precursor").as<bool>();

      // the theoretical spectrum, with the names and charges of its ions
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
      MSSpectrum theoretical;
      generator.getSpectrum(theoretical, peptide, 1, std::max(1, max_charge));
      theoretical.sortByPosition(); // keeps the data arrays aligned

      // the measured peaks, sorted by m/z (SpectrumAlignment needs that)
      MSSpectrum measured;
      measured.reserve(spectrum.size());
      for (const auto& peak : spectrum) measured.push_back(peak);
      measured.sortByPosition();

      SpectrumAlignment aligner;
      Param ap = aligner.getParameters();
      ap.setValue("tolerance", tolerance);
      ap.setValue("is_relative_tolerance", unit == "ppm" ? "true" : "false");
      aligner.setParameters(ap);
      std::vector<std::pair<Size, Size>> alignment;
      aligner.getSpectrumAlignment(alignment, theoretical, measured);

      // the arrays add_metainfo adds: ion names and charges
      const auto names_it = std::find_if(theoretical.getStringDataArrays().begin(), theoretical.getStringDataArrays().end(),
        [](const auto& a) { return a.getName() == Constants::UserParam::IonNames; });
      const auto charges_it = std::find_if(theoretical.getIntegerDataArrays().begin(), theoretical.getIntegerDataArrays().end(),
        [](const auto& a) { return a.getName() == "Charges"; });
      if (names_it == theoretical.getStringDataArrays().end() || charges_it == theoretical.getIntegerDataArrays().end())
      {
        throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The theoretical spectrum has no ion names or charges.");
      }
      const auto& names = *names_it;
      const auto& charges = *charges_it;
      val matches = val::array();
      for (const auto& [t, m] : alignment)
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
  } // namespace
} // namespace OpenMS::WASM

using namespace OpenMS;
using namespace OpenMS::WASM;

EMSCRIPTEN_BINDINGS(openms_prose)
{
  // The results of searchFiles(); delete() releases them.
  class_<ProSEResultSet>("ProSEResultSet")
    .function("summary", OMS_METHOD(ProSEResultSet, &ProSEResultSet::summary))
    .function("storeIdXML", OMS_METHOD(ProSEResultSet, &ProSEResultSet::storeIdXML))
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
}
