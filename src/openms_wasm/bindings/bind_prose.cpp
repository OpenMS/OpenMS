// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// The database search engine ProSE (ANALYSIS/ID/ProSEAlgorithm): a search of one spectrum
// file as the ProSE tool runs it, with its statistics, the modification analysis of open
// searches and progress reports; and the annotation of a spectrum with the fragment ions of
// a peptide, for spectrum viewers.

#include "binding_utils.h"

#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/SYSTEM/StopWatch.h>

#include <pthread.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <atomic>
#include <cctype>
#include <set>
#include <cmath>
#include <utility>

using namespace emscripten;

namespace OpenMS::WASM
{
  using ProteinIdentificationList = std::vector<ProteinIdentification>;

  namespace
  {
    //------------------------------------------------------------------------------
    // Progress

    /// Reports the progress of an algorithm to a JavaScript function, called with
    /// {event: "start" | "progress" | "end", label, depth, begin, end, value}.
    /// The function is called only on the thread that set it (a JavaScript value belongs
    /// to the thread that created it); reports of other threads (OpenMP) only count.
    /// It runs while the search runs: in a Web Worker, postMessage() its argument to the page.
    class JSProgressLogger : public ProgressLogger::ProgressLoggerImpl
    {
    public:
      explicit JSProgressLogger(val callback) :
        callback_(std::move(callback)),
        owner_(pthread_self())
      {
      }

      void startProgress(const SignedSize begin, const SignedSize end, const std::string& label, const int depth) const override
      {
        current_ = begin;
        if (depth == 0)
        {
          begin_ = begin;
          end_ = end;
        }
        report_("start", label, depth, begin, end, begin);
      }

      void setProgress(const SignedSize value, const int depth) const override
      {
        current_ = value;
        report_("progress", "", depth, begin_, end_, value);
      }

      SignedSize nextProgress() const override
      {
        return ++current_;
      }

      void endProgress(const int depth, UInt64 /*bytes_processed*/) const override
      {
        report_("end", "", depth, begin_, end_, end_);
      }

    private:
      void report_(const char* event, const std::string& label, int depth, SignedSize begin, SignedSize end, SignedSize value) const
      {
        if (!pthread_equal(pthread_self(), owner_))
        {
          return;
        }
        val o = val::object();
        o.set("event", std::string(event));
        o.set("label", label);
        o.set("depth", depth);
        o.set("begin", static_cast<double>(begin));
        o.set("end", static_cast<double>(end));
        o.set("value", static_cast<double>(value));
        callback_(o);
      }

      val callback_;
      pthread_t owner_;
      mutable std::atomic<SignedSize> current_{0};
      mutable SignedSize begin_ = 0;
      mutable SignedSize end_ = 0;
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

    //------------------------------------------------------------------------------
    // Percolator rescoring, in-process (Percolator::rescorePSMs). As the ProSE tool does it with
    // -rescore (src/topp/ProSE.cpp, rescoreWithPercolator_ and pinEnzymeName_): keep the two in sync.

    /// The enzyme names Percolator's PIN format knows
    std::string pinEnzymeName(const std::string& enzyme)
    {
      std::string name;
      for (const char c : enzyme)
      {
        if (c != ' ') name += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
      }
      if (name == "trypsin/p") return "trypsinp";
      if (name == "pepsina") return "pepsin";
      static const std::set<std::string> pin_enzymes = {"elastase", "pepsin", "proteinasek", "thermolysin", "chymotrypsin",
                                                        "lys-n", "lys-c", "arg-c", "asp-n", "glu-c", "trypsin", "trypsinp"};
      return pin_enzymes.contains(name) ? name : "no_enzyme";
    }

    /// Replaces the scores of the PSMs by Percolator q-values. Returns why it was skipped, or an
    /// empty string: at least 100 PSMs and decoys are needed.
    std::string rescoreWithPercolator(ProteinIdentificationList& protein_ids, PeptideIdentificationList& peptide_ids)
    {
      if (protein_ids.empty() || peptide_ids.size() < 100)
      {
        return "only " + std::to_string(peptide_ids.size()) + " PSMs, Percolator needs at least 100";
      }
      bool has_decoys = false;
      int min_charge = 10;
      int max_charge = 0;
      for (const PeptideIdentification& pid : peptide_ids)
      {
        for (const PeptideHit& hit : pid.getHits())
        {
          has_decoys = has_decoys || hit.isDecoy();
          min_charge = std::min(min_charge, hit.getCharge());
          max_charge = std::max(max_charge, hit.getCharge());
        }
      }
      if (!has_decoys)
      {
        return "no decoy PSMs for the target/decoy competition (set decoys to auto or generate)";
      }

      const ProteinIdentification::SearchParameters& sp = protein_ids.front().getSearchParameters();
      StringList feature_set = PercolatorInfile::getStandardFeatureSet(min_charge, max_charge);
      if (sp.metaValueExists("extra_features"))
      {
        StringList extra = ListUtils::create<std::string>(sp.getMetaValue("extra_features").toString());
        feature_set.insert(feature_set.end(), extra.begin(), extra.end());
      }
      else
      {
        feature_set.push_back("score");
      }
      feature_set.push_back("Peptide");
      feature_set.push_back("Proteins");

      Percolator perc;
      Param pp = perc.getDefaults();
      pp.setValue("train_best_positive", "true");
      pp.setValue("post_processing_tdc", "true");
      pp.setValue("subset_max_train", 0);
      pp.setValue("pep_method", "nonparametric");
#ifdef _OPENMP
      pp.setValue("num_threads", std::clamp(omp_get_max_threads(), 1, 3)); // one model per fold: more do not help
#else
      pp.setValue("num_threads", 1);
#endif
      perc.setParameters(pp);
      try
      {
        perc.rescorePSMs(peptide_ids, feature_set, pinEnzymeName(sp.digestion_enzyme.getName()), min_charge, max_charge, "q-value");
      }
      catch (const Exception::BaseException& e)
      {
        return std::string("Percolator failed: ") + e.what();
      }
      IDFilter::keepNBestHits(peptide_ids, 1);

      for (ProteinIdentification& run : protein_ids)
      {
        run.setSearchEngine("Percolator");
        run.setSearchEngineVersion("3.08-vendored");
        run.setMetaValue("percolator", "ProSE");
        ProteinIdentification::SearchParameters run_sp = run.getSearchParameters();
        run_sp.setMetaValue("Percolator:peptide_level_fdrs", "false");
        run_sp.setMetaValue("Percolator:protein_level_fdrs", "false");
        run_sp.setMetaValue("Percolator:testFDR", static_cast<double>(pp.getValue("test_fdr")));
        run_sp.setMetaValue("Percolator:trainFDR", static_cast<double>(pp.getValue("train_fdr")));
        run_sp.setMetaValue("Percolator:maxiter", static_cast<int>(pp.getValue("num_iterations")));
        run_sp.setMetaValue("Percolator:subset_max_train", static_cast<int>(pp.getValue("subset_max_train")));
        run_sp.setMetaValue("Percolator:seed", static_cast<int>(pp.getValue("seed")));
        run_sp.setMetaValue("Percolator:post_processing_tdc", "true");
        run_sp.setMetaValue("Percolator:train_best_positive", "true");
        run.setSearchParameters(run_sp);
      }
      return "";
    }

    /// Searches the spectra of @p spectra_path (mzML) against the proteins of @p database_path
    /// (FASTA) as the ProSE tool does for one input file: the search with PSM-level FDR
    /// (FDR:PSM), then, if FDR:protein is set and the database has decoys, protein inference
    /// with picked-protein FDR. By default (options {rescore: false} turn it off) it rescores the PSMs with Percolator
    /// (in-process) as the tool's -rescore does: the PSM and protein FDR then use its q-values.
    JS::ProSESummary proseSearchWithOptions(ProSEAlgorithm& algo, const std::string& spectra_path, const std::string& database_path,
                                            ProteinIdentificationList& proteins, PeptideIdentificationList& peptides,
                                            const JS::ProSEOptions& options)
    {
      // Percolator rescoring is on unless {rescore: false}
      const bool rescore = isNullish(options) || isNullish(options["rescore"]) || options["rescore"].as<bool>();
      const Param user_params = algo.getParameters();
      const double psm_fdr = static_cast<double>(user_params.getValue("FDR:PSM"));
      const double protein_fdr = static_cast<double>(user_params.getValue("FDR:protein"));
      StopWatch watch;
      watch.start();
      // With rescoring, the FDR filters wait for Percolator's q-values (as in the tool).
      if (rescore)
      {
        Param p = user_params;
        p.setValue("FDR:PSM", 0.0);
        p.setValue("FDR:protein", 0.0);
        algo.setParameters(p);
      }
      ProSEAlgorithm::MultiFileSearchResult result;
      try
      {
        result = algo.searchWithModificationAnalysis(std::vector<std::string>{spectra_path}, database_path, {}, "", false);
      }
      catch (...)
      {
        if (rescore) algo.setParameters(user_params);
        throw;
      }
      if (rescore) algo.setParameters(user_params);
      if (result.per_file.size() != 1)
      {
        throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "ProSE returned no result for the input file.");
      }
      ProSEAlgorithm::SearchResult& run = result.per_file[0];
      const bool ok = run.exit_code == ProSEAlgorithm::ExitCodes::EXECUTION_OK;

      val percolator = val::null();
      bool deferred_fdr_applied = false;
      if (rescore && ok)
      {
        StopWatch perc_watch;
        perc_watch.start();
        const std::string skipped = rescoreWithPercolator(run.protein_ids, run.peptide_ids);
        perc_watch.stop();
        percolator = val::object();
        percolator.set("applied", skipped.empty());
        percolator.set("skipped", skipped.empty() ? val::null() : val(skipped));
        percolator.set("seconds", perc_watch.getClockTime());
        // the deferred PSM FDR: on Percolator's q-values, or on ProSE's own if it was skipped
        if (psm_fdr > 0.0 && !run.peptide_ids.empty() && (skipped.empty() || result.have_decoys))
        {
          if (!skipped.empty()) algo.annotatePsmQValues(run.peptide_ids);
          IDFilter::filterHitsByScore(run.peptide_ids, psm_fdr);
          deferred_fdr_applied = true;
        }
      }

      bool protein_fdr_applied = false;
      if (ok && protein_fdr > 0.0 && result.have_decoys && !run.protein_ids.empty() && !run.peptide_ids.empty())
      {
        ProSEAlgorithm::applyCompleteSetProteinFDR(run.protein_ids, run.peptide_ids, result.decoy_string,
                                                   result.decoy_is_prefix, protein_fdr);
        protein_fdr_applied = true;
      }
      // Rescoring, the deferred PSM FDR and the protein FDR changed the PSMs after the search
      // collected its statistics: recount them from the final PSMs, as the ProSE tool does.
      if (ok)
      {
        ProSEAlgorithm::updateFinalStats(run.stats, run.peptide_ids, result.shared.enzyme,
                                         run.stats.fdr_applied || deferred_fdr_applied);
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
  class_<ProSEAlgorithm, base<DefaultParamHandler>>("ProSEAlgorithm")
    .constructor<>()
    .function("search", OMS_FN(&proseSearch))
    .function("search", OMS_FN(&proseSearchWithOptions))
    .function("setProgressCallback", OMS_FN(&proseSetProgressCallback))
    // q-values of the best hits by target-decoy competition (as FDR:PSM, without filtering)
    .function("annotatePsmQValues", OMS_METHOD(ProSEAlgorithm, &ProSEAlgorithm::annotatePsmQValues))
    ;

  function("annotateSpectrum", OMS_FN(&annotateSpectrum));
}
