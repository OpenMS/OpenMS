// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include "ProSEAlgorithmImpl.h"

#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>

#include <algorithm>
#include <cstring>
#include <exception>
#include <limits>
#include <iterator>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace OpenMS
{
  using namespace Internal::ProSEAlgorithmHelpers;

  // =====================================================================
  // Pooling of the searched files of a multi-file search (IDMergerAlgorithm), shared
  // by the search's pooled aggregate and mergeFileResults().
  // =====================================================================
  // static
  void ProSEAlgorithm::poolSearchedFiles_(std::vector<SearchResult>& per_file, bool move_psms,
                                          const StringList& in_spectra_files,
                                          std::vector<ProteinIdentification>& protein_ids,
                                          PeptideIdentificationList& peptide_ids)
  {
    // IDMergerAlgorithm merges only runs of one search engine (and version). With 'rescore', a run that Percolator did
    // not rescore keeps ProSE as its search engine, while the rescored ones report Percolator. The copies of the runs
    // that are pooled all get the engine of the files with PSMs: Percolator's if it rescored all of them (a file without
    // PSMs is not rescored), else ProSE's (the pooled evaluation then scores all files by ProSE's score, see
    // mergeFileResults()); the runs of the files stay as they are.
    // If no file has a PSM left (e.g. FDR:PSM removed all), the engine of the searched files by the same rule:
    // IDMergerAlgorithm requires one engine of runs without PSMs as well.
    const SearchResult* engine_of = nullptr; // the first file that Percolator did not rescore, else the first
    for (const bool require_psms : {true, false})
    {
      for (const SearchResult& pf : per_file)
      {
        if (!searched_(pf) || (require_psms && !hasPSMs_(pf.peptide_ids))) { continue; }
        if (engine_of == nullptr || (engine_of->percolator.applied && !pf.percolator.applied)) { engine_of = &pf; }
      }
      if (engine_of != nullptr) { break; }
    }
    std::string engine, engine_version;
    if (engine_of != nullptr)
    {
      engine = engine_of->protein_ids[0].getSearchEngine();
      engine_version = engine_of->protein_ids[0].getSearchEngineVersion();
    }

    // 'id_merge_index' of a PSM is the input position of its file, per occurrence of the input: IDMergerAlgorithm numbers
    // the primary MS run paths of the runs, so it would give two inputs of one path (the same file given twice) one
    // index, and after a failed file the index k of the k-th pooled file. It appends the PSMs of each insertRuns() call
    // in their order (those of the runs given, i.e. all of a searched file), so the input of each pooled PSM follows from
    // the number each file contributes.
    IDMergerAlgorithm merger;
    std::vector<std::pair<Size, Size>> contributed; // (input position, PSMs) per insertRuns() call
    for (Size i = 0; i < per_file.size(); ++i)
    {
      SearchResult& pf = per_file[i];
      if (pf.exit_code != ExitCodes::EXECUTION_OK) { continue; }
      // Moving: the merger takes (and clears) the PSMs and a copy of the runs, which stay with the files -- their search
      // parameters, meta values (e.g. PeptDeepRescoring's NCE) and run paths are still the caller's.
      std::vector<ProteinIdentification> runs = pf.protein_ids;
      if (!engine.empty())
      {
        for (ProteinIdentification& run : runs)
        {
          run.setSearchEngine(engine);
          run.setSearchEngineVersion(engine_version);
        }
      }
      PeptideIdentificationList psms = move_psms ? std::move(pf.peptide_ids) : PeptideIdentificationList(pf.peptide_ids);
      if (move_psms) { pf.peptide_ids = PeptideIdentificationList(); }
      // (the merger skips the PSMs of no given run, and all PSMs of a file without runs)
      std::unordered_set<std::string> run_ids;
      for (const ProteinIdentification& run : runs) { run_ids.insert(run.getIdentifier()); }
      contributed.emplace_back(i, static_cast<Size>(std::count_if(psms.begin(), psms.end(),
        [&run_ids](const PeptideIdentification& id) { return run_ids.contains(id.getIdentifier()); })));
      merger.insertRuns(std::move(runs), std::move(psms));
    }
    ProteinIdentification merged_proteins;
    merger.returnResultsAndClear(merged_proteins, peptide_ids);
    protein_ids.clear();
    protein_ids.emplace_back(std::move(merged_proteins));
    protein_ids[0].setPrimaryMSRunPath(in_spectra_files);
    Size k = 0;
    for (const auto& [i, n] : contributed)
    {
      for (const Size end = std::min(k + n, peptide_ids.size()); k < end; ++k)
      {
        peptide_ids[k].setMetaValue(Constants::UserParam::ID_MERGE_INDEX, static_cast<int>(i));
      }
    }
    OPENMS_POSTCONDITION(k == peptide_ids.size(), "IDMergerAlgorithm pooled other PSMs than those of the files")
  }

  // =====================================================================
  // 'rescore' of a multi-file search: the files that Percolator could not rescore on their own (e.g. fewer than 100
  // PSMs) get the scores of a model trained on the files that it could rescore, or, if it could rescore none, all
  // files are rescored together.
  // =====================================================================
  Size ProSEAlgorithm::transferQuota_() const
  {
    if (max_train_psms_ == 0) { return std::numeric_limits<Size>::max(); }
    return std::max(Size(1), max_train_psms_ / std::max(Size(1), pool_files_));
  }

  // static
  std::vector<Size> ProSEAlgorithm::transferSample_(const PeptideIdentificationList& ids, Size max_psms, Size hits_per_spectrum)
  {
    // place in the draw: FNV-1a of the spectrum reference, retention time and precursor m/z, mixed (splitmix64)
    const auto priority = [](const PeptideIdentification& id)
    {
      std::uint64_t h = 1469598103934665603ULL;
      const auto add = [&h](const unsigned char* p, Size n) { for (Size k = 0; k < n; ++k) { h ^= p[k]; h *= 1099511628211ULL; } };
      const std::string& ref = id.getSpectrumReference();
      add(reinterpret_cast<const unsigned char*>(ref.data()), ref.size());
      const double rt = id.getRT(), mz = id.getMZ();
      std::uint64_t bits;
      std::memcpy(&bits, &rt, sizeof bits);
      add(reinterpret_cast<const unsigned char*>(&bits), sizeof bits);
      std::memcpy(&bits, &mz, sizeof bits);
      add(reinterpret_cast<const unsigned char*>(&bits), sizeof bits);
      h += 0x9E3779B97F4A7C15ULL;
      h = (h ^ (h >> 30)) * 0xBF58476D1CE4E5B9ULL;
      h = (h ^ (h >> 27)) * 0x94D049BB133111EBULL;
      return h ^ (h >> 31);
    };
    std::vector<std::pair<std::uint64_t, Size>> order;
    Size total = 0;
    for (Size r = 0; r < ids.size(); ++r)
    {
      if (ids[r].getHits().empty()) { continue; }
      order.emplace_back(0, r);
      total += std::min(ids[r].getHits().size(), hits_per_spectrum);
    }
    std::vector<Size> rows;
    rows.reserve(order.size());
    if (total <= max_psms)
    {
      for (const auto& o : order) { rows.push_back(o.second); }
      return rows;
    }
    for (auto& o : order) { o.first = priority(ids[o.second]); }
    std::stable_sort(order.begin(), order.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    Size psms = 0;
    for (const auto& o : order)
    {
      const Size n = std::min(ids[o.second].getHits().size(), hits_per_spectrum);
      if (psms + n > max_psms) { continue; } // (a smaller spectrum later may still fit)
      rows.push_back(o.second);
      psms += n;
    }
    std::sort(rows.begin(), rows.end());
    return rows;
  }

  void ProSEAlgorithm::rescoreFilesTogether_(MultiFileSearchResult& mfres, const std::vector<std::string>& in_spectra_files) const
  {
    // Each file ends as a search without pooling would leave it: no copy of its PSMs before the rescoring unless asked
    // for (setKeepPreFdrPsms()); its proteins stay as the search produced them (see applyPsmFDR_()).
    const auto finish = [this](SearchResult& pf)
    {
      if (!keep_pre_fdr_psms_) { pf.pre_fdr_peptide_ids = PeptideIdentificationList(); }
      else if (rescore_top_hits_ > 0) { keepTopHits_(pf.pre_fdr_peptide_ids, report_top_hits_); } // (see applyPsmFDR_())
    };
    std::vector<Size> searched, own, alone; // files searched; rescored on their own; with PSMs, but not rescored
    StringList alone_names;
    for (Size i = 0; i < mfres.per_file.size(); ++i)
    {
      const SearchResult& pf = mfres.per_file[i];
      if (!searched_(pf)) { continue; }
      searched.push_back(i);
      if (pf.percolator.applied) { own.push_back(i); }
      else if (hasPSMs_(pf.pre_fdr_peptide_ids))
      {
        alone.push_back(i);
        alone_names.push_back(File::basename(in_spectra_files[i]) + " (" + pf.percolator.skipped_reason + ")");
      }
    }
    // in parallel over the files (27 HeLa runs: their proteins and copies one after another took 12 s)
    const auto finish_all = [&]()
    {
      std::vector<std::exception_ptr> errors(searched.size());
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize k = 0; k < static_cast<SignedSize>(searched.size()); ++k)
      {
        try { finish(mfres.per_file[searched[k]]); }
        catch (...) { errors[k] = std::current_exception(); }
      }
      for (const std::exception_ptr& error : errors)
      {
        if (error) { std::rethrow_exception(error); }
      }
    };
    if (alone.empty())
    {
      finish_all();
      return;
    }
#ifdef _OPENMP
    const int threads = omp_get_max_threads();
#else
    const int threads = 1;
#endif
    // Percolator as for the files rescored on their own (rescoring:..., the training cap percolator:max_train_psms)
    const Param pp = rescoringPercolatorParameters(Percolator().getDefaults(), threads, percolator_training_, Param());

    // A file back from the rescoring (its PSMs before the rescoring, now with Percolator's score, q-value and PEP):
    // q-values per file (FDR per run, as for a file rescored on its own), by target-decoy competition of the file's
    // PSMs on the SVM score (MS:1001492), (D + 1) / T as Percolator's own (tdcQValues_(), also 1/T for a file without
    // a decoy PSM); FDR:PSM filters on them.
    const std::string svm_score = "MS:1001492";
    const auto adopt = [&](Size i, PeptideIdentificationList&& ids, const std::string& warning)
    {
      SearchResult& pf = mfres.per_file[i];
      for (PeptideIdentification& id : ids)
      {
        for (PeptideHit& hit : id.getHits()) { hit.setScore(static_cast<double>(hit.getMetaValue(svm_score))); }
        id.setScoreType(svm_score);
        id.setHigherScoreBetter(true);
      }
      if (!ids.empty()) { tdcQValues_(ids, true); }
      for (PeptideIdentification& id : ids)
      {
        for (PeptideHit& hit : id.getHits())
        {
          hit.removeMetaValue(svm_score + "_score"); // (FalseDiscoveryRate's copy of the SVM score)
          hit.setMetaValue("MS:1001491", hit.getScore());
          hit.setMetaValue("percolator_q_value", hit.getScore());
        }
      }
      // rescoring:stratum_q_values: as for a file rescored on its own (applyRescoredPsmFDR_()), before FDR:PSM
      if (stratum_q_values_ && !ids.empty())
      {
        annotateStratumQValues(ids);
        for (ProteinIdentification& run : pf.protein_ids) run.getSearchParameters().setMetaValue("rescoring:stratum_q_values", "true");
      }
      pf.inference_only_peptide_ids = PeptideIdentificationList();
      bool filtered = false;
      if (fdr_psm_ > 0.0)
      {
        if (fdr_protein_ > 0.0) { collectInferenceOnlyPsms_(ids, fdr_psm_, pf.inference_only_peptide_ids); }
        IDFilter::filterHitsByScore(ids, fdr_psm_);
        filtered = true;
      }
      pf.peptide_ids = std::move(ids);
      markRescoredByPercolator_(pf.protein_ids, pp, percolator_training_);
      pf.percolator.applied = true;
      pf.percolator.skipped_reason.clear();
      pf.percolator.warning = warning;
      updateFinalStats(pf.stats, pf.peptide_ids, enzyme_, filtered);
      // (an open search keeps the modification analysis of the PSMs of the file's own scoring)
    };

    // Percolator gets the first rescoring:top_hits hits of a spectrum, as for a file rescored on its own (the copies
    // before the rescoring may hold more, for a .pin with report:top_hits above it; see applyPsmFDR_())
    const Size rescored_hits = rescore_top_hits_ > 0 ? rescoreTopHits_() : std::numeric_limits<Size>::max();
    const auto hits_rescored = [rescored_hits](const PeptideIdentification& id) { return std::min(id.getHits().size(), rescored_hits); };
    const auto cut_to_rescored = [&](PeptideIdentificationList& ids) { if (rescore_top_hits_ > 0) keepTopHits_(ids, rescored_hits); };

    // The precursor charges of the charge features: of all PSMs the model sees, for one feature set
    int min_charge = std::numeric_limits<int>::max(), max_charge = std::numeric_limits<int>::min();
    for (Size i : searched)
    {
      for (const PeptideIdentification& id : mfres.per_file[i].pre_fdr_peptide_ids)
      {
        for (Size h = 0; h < hits_rescored(id); ++h)
        {
          min_charge = std::min(min_charge, id.getHits()[h].getCharge());
          max_charge = std::max(max_charge, id.getHits()[h].getCharge());
        }
      }
    }
    const ProteinIdentification& first_run = mfres.per_file[searched[0]].protein_ids[0];
    const StringList feature_set = PercolatorInfile::getFeatureSet(first_run.getSearchParameters(), min_charge, max_charge);
    const std::string enzyme = rescoringEnzyme(first_run.getSearchParameters().digestion_enzyme.getName(), percolator_training_);

    // rescore_method 'reset' has no transfer of a model (Percolator::trainPSMs()/scorePSMs() are plain Percolator): the
    // files too small for a rescoring of their own keep their search engine scores, as before.
    if (!own.empty() && percolator_extra_.exists("algorithm"))
    {
      OPENMS_LOG_WARN << "[ProSE] Percolator could not rescore " << alone.size() << " of " << searched.size()
                      << " input files on their own: " << ListUtils::concatenate(alone_names, ", ")
                      << ". rescore_method 'reset' cannot score them with a model of the other files; they keep their "
                      << "search engine scores." << std::endl;
      finish_all();
      return;
    }
    if (!own.empty())
    {
      // Transfer: the files rescored on their own keep their results. A model trained on their PSMs (before the
      // rescoring; each file at most percolator:max_train_psms / (the input files), so that none dominates) scores the
      // others, on the scale of its training PSMs: a PSM gets the score it would get among them, and the PEP of the
      // training PSMs at that score (Percolator::scorePSMs()), not one estimated on a file too small for that.
      // 'id_merge_index' keeps the spectra of the training files apart (spectrum file numbers).
      // At most percolator:max_train_psms / (the input files) PSMs of each file (Percolator's rows: the hits it rescores),
      // reproducible random whole spectra (transferSample_(): a file rescored on its own kept only these, unless
      // setKeepPreFdrPsms(), see applyPsmFDR_())
      PeptideIdentificationList training;
      Size available = 0;
      for (Size i : own)
      {
        const PeptideIdentificationList& ids = mfres.per_file[i].pre_fdr_peptide_ids;
        for (const PeptideIdentification& id : ids) { available += hits_rescored(id); }
        for (Size r : transferSample_(ids, transferQuota_(), rescored_hits))
        {
          PeptideIdentification& copy = training.emplace_back(ids[r]);
          copy.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, static_cast<int>(i));
        }
      }
      cut_to_rescored(training);
      Size training_psms = 0;
      for (const PeptideIdentification& id : training) { training_psms += id.getHits().size(); }
      OPENMS_LOG_INFO << "[ProSE] Percolator could not rescore " << alone.size() << " of " << searched.size()
                      << " input files on their own: " << ListUtils::concatenate(alone_names, ", ")
                      << ". Scoring them with a model trained on " << training_psms << " PSMs of the " << own.size()
                      << " files it rescored on their own";
      if (training_psms < available)
      {
        OPENMS_LOG_INFO << " (sampled from " << available << ": at most " << transferQuota_() << " of each file, "
                        << "Search:percolator:max_train_psms " << max_train_psms_ << " / " << pool_files_ << " input files)";
      }
      OPENMS_LOG_INFO << "." << std::endl;
      try
      {
        Percolator perc;
        perc.setParameters(pp);
        const PercolatorModel model = perc.trainPSMs(training, feature_set, enzyme, min_charge, max_charge);
        const std::string transferred = "its PSMs are scored with a model trained on the PSMs of the "
          + StringUtils::toStr(own.size()) + " input(s) that Percolator rescored on their own, on the scale of that "
          "model's training PSMs (their PEP at the same score); the q-values are those of target-decoy competition of "
          "this input's PSMs on that score." + (perc.usedDefaultDirection() ? " Percolator's sanity check rejected "
          "the learned model, so the PSMs are scored with the best single feature." : "");
        // The PSMs of all these files in one call: it scores every PSM on its own (no target-decoy competition among
        // them), and prepares and scores the training PSMs (the reference) once.
        PeptideIdentificationList ids;
        std::vector<Size> file_of; // per element of ids
        for (Size i : alone)
        {
          PeptideIdentificationList& pre_fdr = mfres.per_file[i].pre_fdr_peptide_ids;
          for (PeptideIdentification& id : pre_fdr)
          {
            if (keep_pre_fdr_psms_) { ids.push_back(id); }
            else { ids.push_back(std::move(id)); }
            file_of.push_back(i);
          }
          if (!keep_pre_fdr_psms_) { pre_fdr = PeptideIdentificationList(); }
        }
        // the best hit of each spectrum by the SVM score (the q-values of the training PSMs at a score have plateaus,
        // which would leave ties to the search engine's order); target-decoy competition follows in adopt()
        cut_to_rescored(ids);
        perc.scorePSMs(ids, model, training, feature_set, enzyme, min_charge, max_charge, "svm");
        IDFilter::keepNBestHits(ids, 1);
        std::vector<PeptideIdentificationList> parts(mfres.per_file.size());
        for (Size k = 0; k < ids.size(); ++k) { parts[file_of[k]].push_back(std::move(ids[k])); }
        ids = PeptideIdentificationList();
        for (Size i : alone)
        {
          adopt(i, std::move(parts[i]), "Percolator could not rescore this input on its own ("
                                        + mfres.per_file[i].percolator.skipped_reason + "): " + transferred);
        }
      }
      catch (const Exception::BaseException& e)
      {
        OPENMS_LOG_WARN << "[ProSE] Percolator cannot score the files too small for a model of their own with a model "
                        << "of the other files either (" << e.what() << "); they keep their search engine scores."
                        << std::endl;
      }
      finish_all();
      return;
    }

    // No file could be rescored on its own: the PSMs of all files together, with one model, cross-validated as for one
    // file (Percolator samples at most percolator:max_train_psms PSMs to train on; every file has too few PSMs for a
    // model of its own, so none can dominate). 'id_merge_index' keeps the files apart (spectrum file numbers, for the
    // folds and the target-decoy competition of each spectrum).
    PeptideIdentificationList pooled;
    std::vector<Size> file_of; // per element of pooled
    std::vector<ProteinIdentification> runs = mfres.per_file[searched[0]].protein_ids; // the search parameters
    for (Size i : searched)
    {
      for (const PeptideIdentification& id : mfres.per_file[i].pre_fdr_peptide_ids)
      {
        PeptideIdentification& copy = pooled.emplace_back(id);
        copy.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, static_cast<int>(i));
        file_of.push_back(i);
      }
    }
    OPENMS_LOG_INFO << "[ProSE] Percolator could not rescore any of the " << searched.size()
                    << " input files on its own: " << ListUtils::concatenate(alone_names, ", ")
                    << ". Rescoring the " << pooled.size() << " PSMs of all files together, with one model." << std::endl;
    // (rescore_method 'reset': without a stem table, which names one run)
    Param extra = percolator_extra_;
    if (extra.exists("reset:stem_table")) { extra.setValue("reset:stem_table", ""); }
    cut_to_rescored(pooled);
    const PercolatorRescoring outcome = rescoreWithPercolator_(runs, pooled, threads, nullptr, percolator_training_, extra);
    if (!outcome.applied)
    {
      OPENMS_LOG_WARN << "[ProSE] Percolator cannot rescore the PSMs of all input files together either ("
                      << outcome.skipped_reason << "); they keep their search engine scores." << std::endl;
      finish_all();
      return;
    }
    std::vector<PeptideIdentificationList> parts(mfres.per_file.size());
    for (Size k = 0; k < pooled.size(); ++k)
    {
      pooled[k].removeMetaValue(Constants::UserParam::ID_MERGE_INDEX);
      parts[file_of[k]].push_back(std::move(pooled[k]));
    }
    pooled = PeptideIdentificationList();
    const std::string warning = "Percolator could not rescore any input on its own (too few PSMs or decoys); the PSMs of "
      "all " + StringUtils::toStr(searched.size()) + " inputs are rescored together, with one model; the q-values are "
      "those of target-decoy competition of this input's PSMs on that model's score." + (outcome.warning.empty() ? std::string() : " " + outcome.warning);
    for (Size i : searched) { adopt(i, std::move(parts[i]), warning); }
    finish_all();
  }

  // =====================================================================
  // The pooled evaluation of a multi-file search (the ProSE tool's -out_merged):
  // pooling, cross-file protein inference and the optional picked-protein FDR.
  // =====================================================================
  ProSEAlgorithm::MergedResult
  ProSEAlgorithm::mergeFileResults(MultiFileSearchResult& result, const StringList& in_spectra_files,
                                   bool remove_decoys) const
  {
    // The pooled run lists the input files, and 'id_merge_index' indexes them: one per file of the search.
    if (in_spectra_files.size() != result.per_file.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "mergeFileResults(): " + StringUtils::toStr(in_spectra_files.size()) + " input file(s) given for the result of a "
        "search of " + StringUtils::toStr(result.per_file.size()) + " file(s); pass the input files of the search, in its "
        "order.");
    }

    // Percolator ('rescore') rescores each file on its own; a file it could not rescore (e.g. one with fewer than 100
    // PSMs) keeps ProSE's scores. The PSMs of each file passed FDR:PSM on the q-values of their file, Percolator's or
    // target-decoy q-values of ProSE's score, which are valid per file and stay. Only the protein inference compares
    // PSMs across files, by the best PEP of a protein's PSMs: the PSMs that Percolator did not rescore have none. A
    // mixed set is therefore pooled on ProSE's own score, which every PSM has: the inference scores proteins by the
    // best 'ln(hyperscore)' (as without rescoring, see inferProteins(), which switches the scores of its copies of the
    // PSMs only); the PSMs themselves keep their scores.
    Size n_searched = 0, n_rescored = 0;
    std::vector<std::string> not_rescored;
    for (Size i = 0; i < result.per_file.size(); ++i)
    {
      const SearchResult& pf = result.per_file[i];
      if (pf.exit_code != ExitCodes::EXECUTION_OK || !hasPSMs_(pf.peptide_ids)) { continue; } // (no PSM to compare)
      ++n_searched;
      if (pf.percolator.applied) { ++n_rescored; }
      else { not_rescored.push_back(File::basename(in_spectra_files[i]) + ": " + pf.percolator.skipped_reason); }
    }
    const bool rescoring_mixed = n_rescored > 0 && n_rescored < n_searched;
    if (rescoring_mixed)
    {
      OPENMS_LOG_WARN << "[ProSE] Percolator rescored " << n_rescored << " of " << n_searched << " input files; not "
                      << ListUtils::concatenate(not_rescored, " ") << " The pooled evaluation scores the proteins of "
                      << "all files by ProSE's score ln(hyperscore) instead of Percolator's PEP, which only the rescored "
                      << "PSMs have; each PSM keeps the q-value of its file." << std::endl;
    }

    MergedResult merged;
    if (n_searched == 0) { return merged; } // no file was searched: nothing to pool

    // Use the effective decoy marker/position the search resolved (which may be an auto-detected external marker,
    // prefix or suffix), not the raw decoy_prefix parameter.
    const std::string& decoy_string = result.decoy_string;
    const bool decoy_is_prefix = result.decoy_is_prefix;

    // The multi-file search already produced exactly this merge: it pools the per-file results with
    // poolSearchedFiles_(), as below -- so re-merging here would only build a third copy of every PSM. Consume the
    // aggregate by move instead. Without an aggregate (build_pooled_aggregate false), pool the runs of the files here.
    const bool pooled_from_files = result.aggregate.protein_ids.empty();
    if (!pooled_from_files)
    {
      merged.protein_ids.emplace_back(std::move(result.aggregate.protein_ids[0]));
      merged.peptide_ids = std::move(result.aggregate.peptide_ids);
      result.aggregate.protein_ids.clear();
      result.aggregate.peptide_ids.clear();
    }
    else
    {
      // The PSMs are moved out of per_file, as documented; its runs are copied and stay, as do the statistics, exit
      // codes, Percolator outcomes and modification analyses.
      poolSearchedFiles_(result.per_file, /*move_psms=*/true, in_spectra_files, merged.protein_ids, merged.peptide_ids);
    }

    std::vector<ProteinIdentification>& merged_protein_ids = merged.protein_ids;
    PeptideIdentificationList& merged_peptides = merged.peptide_ids;
    merged_protein_ids[0].setPrimaryMSRunPath(in_spectra_files);
    // IDMergerAlgorithm carries over search engine/params but not a run date; set one so the
    // merged idXML is schema-valid (matches the single-file path).
    merged_protein_ids[0].setDateTime(DateTime::now());

    // If the evaluation fails, the files get the PSMs moved out of them back (as far as it changed them: at most the
    // best hit of each kept, or q-values as below), so that their own outputs can still be written.
    try
    {
      // One main score type for all PSMs of the pooled run: with rescoring_mixed and FDR:PSM 0, the files that Percolator
      // did not rescore still have ProSE's score as main score, the others Percolator's q-values; a reader of the run
      // (which takes the score type of its first PSM) would misread one of them. These files get target-decoy q-values
      // of ProSE's score per file instead, as with FDR:PSM > 0 (unfiltered).
      if (rescoring_mixed)
      {
        std::map<Size, std::vector<Size>> rows_of; // per file not rescored whose PSMs have no q-values: its pooled PSMs
        for (Size k = 0; k < merged_peptides.size(); ++k)
        {
          const PeptideIdentification& id = merged_peptides[k];
          if (!id.metaValueExists(Constants::UserParam::ID_MERGE_INDEX)) { continue; }
          const Size i = static_cast<Size>(static_cast<int>(id.getMetaValue(Constants::UserParam::ID_MERGE_INDEX)));
          if (i >= result.per_file.size() || result.per_file[i].percolator.applied || id.getScoreType() == "q-value") { continue; }
          rows_of[i].push_back(k);
        }
        for (const auto& [i, rows] : rows_of)
        {
          PeptideIdentificationList ids;
          for (Size k : rows) { ids.push_back(std::move(merged_peptides[k])); }
          annotatePsmQValues(ids);
          for (Size j = 0; j < rows.size(); ++j) { merged_peptides[rows[j]] = std::move(ids[j]); }
        }
      }

      // Cross-file protein inference (inferProteins(): the best PEP, or search engine score, per peptide per protein) and
      // the optional picked-protein FDR. Picked-protein FDR needs identified decoy proteins: gate on DB-level decoys here,
      // applyCompleteSetProteinFDR() checks the result of the inference (it could be target-only even with a decoy
      // database if no decoy survived in any file). It infers the proteins on the pooled PSMs together with the PSMs of
      // the files between FDR:PSM and FDR:protein_inference_PSM, and finalizes the result: with remove_decoys it strips
      // the decoys, and in any case the proteins without a pooled PSM and the references to removed proteins (peptide
      // evidences, indistinguishable proteins and protein groups), so that the result can still be stored as idXML.
      // Without the protein FDR the merged output is a pooled intermediate that retains the target+decoy evidence for
      // a downstream FDR.
      const bool protein_fdr = fdr_protein_ > 0.0 && result.have_decoys && !decoy_string.empty();
      startProgress(0, 1, protein_fdr ? "Pooled protein inference and protein FDR" : "Pooled protein inference");
      const OnScopeExit end_progress([this]() { endProgress(); }); // also when it fails
      if (protein_fdr)
      {
        PeptideIdentificationList inference_only;
        for (SearchResult& pf : result.per_file)
        {
          if (pf.exit_code != ExitCodes::EXECUTION_OK) { continue; }
          inference_only.insert(inference_only.end(), std::make_move_iterator(pf.inference_only_peptide_ids.begin()),
                                std::make_move_iterator(pf.inference_only_peptide_ids.end()));
          pf.inference_only_peptide_ids = PeptideIdentificationList();
        }
        merged.protein_fdr_applied = applyCompleteSetProteinFDR(merged_protein_ids, merged_peptides, decoy_string,
                                                                decoy_is_prefix, fdr_protein_, remove_decoys, inference_only);
      }
      else
      {
        IDFilter::keepNBestHits(merged_peptides, 1); // as the protein FDR keeps them
        inferProteins(merged_protein_ids, merged_peptides);
      }
    }
    catch (...)
    {
      if (pooled_from_files) { returnPooledPsms_(merged.peptide_ids, result.per_file); }
      throw;
    }

    return merged;
  }

  // static
  bool ProSEAlgorithm::searched_(const SearchResult& pf)
  {
    return pf.exit_code == ExitCodes::EXECUTION_OK && !pf.protein_ids.empty();
  }

  // static
  int ProSEAlgorithm::searchedFileOf_(const PeptideIdentification& id, const std::vector<SearchResult>& per_file)
  {
    if (!id.metaValueExists(Constants::UserParam::ID_MERGE_INDEX)) { return -1; }
    const int i = static_cast<int>(id.getMetaValue(Constants::UserParam::ID_MERGE_INDEX));
    return i >= 0 && static_cast<Size>(i) < per_file.size() && searched_(per_file[i]) ? i : -1;
  }

  // static
  void ProSEAlgorithm::returnPooledPsms_(PeptideIdentificationList& pooled, std::vector<SearchResult>& per_file)
  {
    for (PeptideIdentification& id : pooled)
    {
      const int i = searchedFileOf_(id, per_file);
      if (i < 0) { continue; }
      id.removeMetaValue(Constants::UserParam::ID_MERGE_INDEX);
      id.setIdentifier(per_file[i].protein_ids[0].getIdentifier());
      per_file[i].peptide_ids.push_back(std::move(id));
    }
    pooled = PeptideIdentificationList();
  }

  // static
  void ProSEAlgorithm::splitMergedPsms(MergedResult& merged, MultiFileSearchResult& result)
  {
    for (const PeptideIdentification& id : merged.peptide_ids)
    {
      if (searchedFileOf_(id, result.per_file) < 0)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "splitMergedPsms(): a pooled PSM has no 'id_merge_index' of a file of the search that was searched.");
      }
    }
    for (SearchResult& pf : result.per_file)
    {
      if (searched_(pf)) { pf.peptide_ids = PeptideIdentificationList(); }
    }
    returnPooledPsms_(merged.peptide_ids, result.per_file);
    if (merged.protein_ids.empty())
    {
      for (SearchResult& pf : result.per_file)
      {
        if (searched_(pf)) { removeUnreportedProteins(pf.protein_ids, pf.peptide_ids); }
      }
      return;
    }
    // The pooled run's scores (e.g. the protein q-values), and of its proteins and groups those that the file's PSMs
    // reference, the groups reduced to these proteins: indices of the pooled run, built once
    const ProteinIdentification& pooled = merged.protein_ids[0];
    std::unordered_map<std::string, Size> hit_index;
    for (Size h = 0; h < pooled.getHits().size(); ++h) { hit_index.emplace(pooled.getHits()[h].getAccession(), h); }
    const auto group_index = [](const std::vector<ProteinIdentification::ProteinGroup>& groups)
    {
      std::unordered_map<std::string, std::vector<Size>> index;
      for (Size g = 0; g < groups.size(); ++g)
      {
        for (const std::string& accession : groups[g].accessions) { index[accession].push_back(g); }
      }
      return index;
    };
    const std::vector<ProteinIdentification::ProteinGroup>& pooled_indistinguishable = pooled.getIndistinguishableProteins();
    const std::vector<ProteinIdentification::ProteinGroup>& pooled_groups = pooled.getProteinGroups();
    const auto indistinguishable_index = group_index(pooled_indistinguishable);
    const auto groups_index = group_index(pooled_groups);
    for (SearchResult& pf : result.per_file)
    {
      if (!searched_(pf)) { continue; }
      // the pooled proteins that the file's PSMs reference, in the pooled order
      std::unordered_set<std::string> referenced;
      for (const PeptideIdentification& id : pf.peptide_ids)
      {
        for (const PeptideHit& hit : id.getHits())
        {
          for (const PeptideEvidence& evidence : hit.getPeptideEvidences()) { referenced.insert(evidence.getProteinAccession()); }
        }
      }
      std::vector<Size> hit_rows;
      for (const std::string& accession : referenced)
      {
        const auto it = hit_index.find(accession);
        if (it != hit_index.end()) { hit_rows.push_back(it->second); }
      }
      std::sort(hit_rows.begin(), hit_rows.end());
      std::vector<ProteinHit> hits;
      hits.reserve(hit_rows.size());
      std::unordered_set<std::string> kept;
      for (Size h : hit_rows)
      {
        hits.push_back(pooled.getHits()[h]);
        kept.insert(hits.back().getAccession());
      }
      // the pooled groups with one of these proteins, reduced to them (as IDFilter::updateProteinGroups()), in order
      const auto select = [&kept](const std::vector<ProteinIdentification::ProteinGroup>& groups,
                                  const std::unordered_map<std::string, std::vector<Size>>& index)
      {
        std::vector<Size> rows;
        for (const std::string& accession : kept)
        {
          const auto it = index.find(accession);
          if (it != index.end()) { rows.insert(rows.end(), it->second.begin(), it->second.end()); }
        }
        std::sort(rows.begin(), rows.end());
        rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
        std::vector<ProteinIdentification::ProteinGroup> out;
        out.reserve(rows.size());
        for (Size g : rows)
        {
          const ProteinIdentification::ProteinGroup& group = groups[g];
          ProteinIdentification::ProteinGroup& reduced = out.emplace_back();
          for (const std::string& accession : group.accessions)
          {
            if (kept.contains(accession)) { reduced.accessions.push_back(accession); }
          }
          reduced.probability = group.probability;
          reduced.setFloatDataArrays(group.getFloatDataArrays());
          reduced.setStringDataArrays(group.getStringDataArrays());
          reduced.setIntegerDataArrays(group.getIntegerDataArrays());
        }
        return out;
      };
      std::vector<ProteinIdentification::ProteinGroup> indistinguishable = select(pooled_indistinguishable, indistinguishable_index);
      std::vector<ProteinIdentification::ProteinGroup> groups = select(pooled_groups, groups_index);
      for (ProteinIdentification& run : pf.protein_ids)
      {
        run.setScoreType(pooled.getScoreType());
        run.setHigherScoreBetter(pooled.isHigherScoreBetter());
        run.setInferenceEngine(pooled.getInferenceEngine());
        run.setInferenceEngineVersion(pooled.getInferenceEngineVersion());
        // (the file's PSMs are those of its first run; another run keeps no protein)
        run.getHits().clear();
        run.getIndistinguishableProteins().clear();
        run.getProteinGroups().clear();
      }
      pf.protein_ids[0].getHits() = std::move(hits);
      pf.protein_ids[0].getIndistinguishableProteins() = std::move(indistinguishable);
      pf.protein_ids[0].getProteinGroups() = std::move(groups);
      // (a PSM referencing a protein that the pooled run does not have: as for any run)
      if (kept.size() < referenced.size()) { IDFilter::removeDanglingProteinReferences(pf.peptide_ids, pf.protein_ids); }
    }
  }
} // namespace OpenMS
