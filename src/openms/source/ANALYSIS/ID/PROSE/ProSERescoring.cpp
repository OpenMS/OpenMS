// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PROSE/ProSEAlgorithm.h>
#include "ProSEImpl.h"

#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <exception>
#include <limits>
#include <numeric>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <thread>
#include <tuple>
#include <utility>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace OpenMS
{
  using namespace Internal::ProSEAlgorithmHelpers;

  namespace
  {
    // IDFilter::filterHitsByScore(ids, threshold), spectrum by spectrum in parallel: the same hits remain in the same
    // order, and the hits that fail the threshold are destroyed on the thread that removes them. On one thread,
    // destroying the hits of the ~140,000 spectra of a timsTOF run that fail FDR:PSM (each with its meta values and peak
    // annotations) took about a second per run while the other threads waited (5 of 44 stack samples of a 27-run search).
    void filterHitsByScoreInParallel(PeptideIdentificationList& ids, double threshold)
    {
      const SignedSize n = static_cast<SignedSize>(ids.size());
      std::exception_ptr error; // (an allocation failure in a move; must not leave the parallel region)
#pragma omp parallel for schedule(dynamic, 256) if (n >= 1024)
      for (SignedSize i = 0; i < n; ++i)
      {
        try
        {
          IDFilter::filterHitsByScore(ids[i], threshold);
        }
        catch (...)
        {
#pragma omp critical (ProSEAlgorithm_filter_hits_error)
          {
            if (!error) error = std::current_exception();
          }
        }
      }
      if (error) std::rethrow_exception(error);
    }
  }

  // static
  bool ProSEAlgorithm::hasPSMs_(const PeptideIdentificationList& ids)
  {
    return std::any_of(ids.begin(), ids.end(), [](const PeptideIdentification& id) { return !id.getHits().empty(); });
  }

  // static
  void ProSEAlgorithm::tdcQValues_(PeptideIdentificationList& peptide_ids, bool plus_one)
  {
    // With plus_one: the best hit of each PSM (its score, oriented higher-better, and whether it is a decoy), which
    // FalseDiscoveryRate keeps, before it replaces the scores.
    std::vector<std::tuple<double, bool, Size>> best; // score, decoy, PSM
    if (plus_one)
    {
      for (Size i = 0; i < peptide_ids.size(); ++i)
      {
        PeptideIdentification& id = peptide_ids[i];
        if (id.getHits().empty()) { continue; }
        id.sort();
        const PeptideHit& hit = id.getHits()[0];
        best.emplace_back(id.isHigherScoreBetter() ? hit.getScore() : -hit.getScore(), hit.isDecoy(), i);
      }
    }
    FalseDiscoveryRate fdr;
    Param fdr_params = fdr.getParameters();
    fdr_params.setValue("add_decoy_peptides", "true"); // keep decoys eligible (q-value filtered, but no decoy-specific stripping)
    fdr.setParameters(fdr_params);
    // Without a decoy PSM, FalseDiscoveryRate's D/T is 0 for every PSM, which would let any threshold pass PSMs
    // whose FDR nothing estimates. Use the +1-corrected estimate (D + 1) / T instead: 1/T for every PSM (the q-value
    // of the best of T targets, and of all below it), so that a threshold below 1/T reports no PSM, as an
    // unattainable protein FDR reports no protein.
    fdr.apply(peptide_ids); // keeps the best hit of each PSM, decoys included (add_decoy_peptides)
    Size targets = 0, decoys = 0;
    for (const PeptideIdentification& id : peptide_ids)
    {
      if (!id.getHits().empty()) { (id.getHits()[0].isDecoy() ? decoys : targets) += 1; }
    }
    if (decoys == 0 && targets > 0)
    {
      OPENMS_LOG_WARN << "[ProSE] FDR:PSM: none of the " << targets << " PSMs matched a decoy; their q-value is the "
                      << "conservative estimate 1/" << targets << " = " << 1.0 / static_cast<double>(targets)
                      << " (no PSM passes a lower threshold)." << std::endl;
      if (!plus_one)
      {
        for (PeptideIdentification& id : peptide_ids)
        {
          for (PeptideHit& hit : id.getHits()) { hit.setScore(1.0 / static_cast<double>(targets)); }
        }
      }
    }
    if (plus_one)
    {
      // (D + 1) / T of the PSMs scoring at least as well (ties together), the minimum of it at and below each score
      // (Percolator's PosteriorEstimator::getQValues() without pi0), at most 1. Without a decoy, 1/T for every PSM.
      std::stable_sort(best.begin(), best.end(), [](const auto& a, const auto& b) { return std::get<0>(a) > std::get<0>(b); });
      std::vector<double> q(best.size());
      Size d = 0, t = 0;
      for (Size k = 0; k < best.size();)
      {
        Size end = k; // past the PSMs tied with PSM k
        for (; end < best.size() && std::get<0>(best[end]) == std::get<0>(best[k]); ++end) { (std::get<1>(best[end]) ? d : t) += 1; }
        std::fill(q.begin() + k, q.begin() + end, std::min(1.0, static_cast<double>(d + 1) / static_cast<double>(std::max<Size>(1, t))));
        k = end;
      }
      for (Size k = best.size(); k-- > 1;) { q[k - 1] = std::min(q[k - 1], q[k]); }
      for (Size k = 0; k < best.size(); ++k)
      {
        for (PeptideHit& hit : peptide_ids[std::get<2>(best[k])].getHits()) { hit.setScore(q[k]); }
      }
    }
  }

  void ProSEAlgorithm::annotatePsmQValues(PeptideIdentificationList& peptide_ids) const
  {
    FalseDiscoveryRate fdr;
    Param fdr_params = fdr.getParameters();
    fdr_params.setValue("add_decoy_peptides", "true"); // keep decoys eligible (q-value filtered, but no decoy-specific stripping)
    fdr.setParameters(fdr_params);

    // Group the spectra by the number of fragment charges their best hit was scored with. FalseDiscoveryRate
    // competes the best hit of each spectrum (stable score order, as sort() below).
    std::map<int, std::vector<Size>> groups;
    bool separate = fdr_psm_by_scored_charges_ && scoring_multiple_charges_;
    if (separate)
    {
      std::map<int, std::pair<bool, bool>> has_target_decoy;
      std::vector<Size> without_hits;
      for (Size i = 0; i < peptide_ids.size(); ++i)
      {
        PeptideIdentification& id = peptide_ids[i];
        if (id.getHits().empty())
        {
          without_hits.push_back(i);
          continue;
        }
        id.sort();
        const PeptideHit& best = id.getHits()[0];
        const int group = scoringMaxCharge_(best.getCharge());
        groups[group].push_back(i);
        auto& [has_target, has_decoy] = has_target_decoy[group];
        (best.isDecoy() ? has_decoy : has_target) = true;
      }
      separate = groups.size() > 1
                 && std::all_of(has_target_decoy.begin(), has_target_decoy.end(),
                                [](const auto& entry) { return entry.second.first && entry.second.second; });
      if (groups.size() > 1 && ! separate)
      {
        OPENMS_LOG_WARN << "[ProSE] FDR:PSM_groups: a group of PSMs scored with the same number of fragment charges has no "
                        << "target or no decoy PSM; computing q-values over all PSMs together." << std::endl;
      }
      if (separate && ! without_hits.empty()) // they carry no score, but get the q-value score type with the others
      {
        auto& first = groups.begin()->second;
        first.insert(first.end(), without_hits.begin(), without_hits.end());
      }
    }
    if (! separate)
    {
      tdcQValues_(peptide_ids);
      return;
    }
    for (const auto& [group, indices] : groups)
    {
      PeptideIdentificationList part;
      part.reserve(indices.size());
      for (Size i : indices) { part.push_back(std::move(peptide_ids[i])); }
      fdr.apply(part);
      for (Size k = 0; k < indices.size(); ++k) { peptide_ids[indices[k]] = std::move(part[k]); }
    }
  }

  namespace
  {
    // A helper thread that destroys hits, joined when this goes (see applyPsmFDR_()).
    struct HitRelease
    {
      std::thread thread;

      ~HitRelease()
      {
        if (thread.joinable()) { thread.join(); }
      }
    };
  }

  // static
  bool ProSEAlgorithm::annotatePeptideQValues(PeptideIdentificationList& peptide_ids, Size* target_peptides_1pct)
  {
    if (target_peptides_1pct != nullptr) { *target_peptides_1pct = 0; }
    const std::string& name = Constants::UserParam::PEPTIDE_Q_VALUE;
    const auto remove_all = [&]()
    {
      for (PeptideIdentification& id : peptide_ids)
      {
        for (PeptideHit& hit : id.getHits()) { hit.removeMetaValue(name); }
      }
    };

    // FalseDiscoveryRate::applyBasicPeptideLevel() on a slim copy of the first hits: their sequence, label and PEP as
    // the score (lower is better; a PEP that is not a number ranks last). A copy, since it replaces the main scores and
    // removes the other hits; with add_decoy_peptides, the decoys keep their hit (and get the q-value of their peptide).
    constexpr Size NO_COPY = std::numeric_limits<Size>::max();
    std::vector<Size> copy_of_id(peptide_ids.size(), NO_COPY);
    PeptideIdentificationList copies;
    for (Size i = 0; i < peptide_ids.size(); ++i)
    {
      const std::vector<PeptideHit>& hits = peptide_ids[i].getHits();
      if (hits.empty()) { continue; }
      const PeptideHit& hit = hits.front();
      if (!hit.metaValueExists("percolator_pep") || !hit.metaValueExists("target_decoy"))
      {
        remove_all();
        return false;
      }
      double pep = static_cast<double>(hit.getMetaValue("percolator_pep"));
      if (std::isnan(pep)) { pep = std::numeric_limits<double>::infinity(); }
      PeptideHit copy_hit(pep, 1, hit.getCharge(), hit.getSequence());
      copy_hit.setMetaValue("target_decoy", hit.getMetaValue("target_decoy"));
      PeptideIdentification& copy = copies.emplace_back();
      copy.setScoreType("Posterior Error Probability");
      copy.setHigherScoreBetter(false);
      copy.insertHit(std::move(copy_hit));
      copy_of_id[i] = copies.size() - 1;
    }
    if (copies.empty()) { return true; } // (no hit to annotate)
    FalseDiscoveryRate fdr;
    Param fdr_param = fdr.getParameters();
    fdr_param.setValue("conservative", "true"); // (decoys + 1) / (targets + 1), as calculateFDRBasic_() computes it
    fdr_param.setValue("add_decoy_peptides", "true");
    fdr_param.setValue("no_qvalues", "false");
    fdr.setParameters(fdr_param);
    fdr.applyBasicPeptideLevel(copies);

    const UInt index = MetaInfoInterface::metaRegistry().registerName(name);
    std::unordered_set<std::string> target_peptides; // at q <= 0.01
    for (Size i = 0; i < peptide_ids.size(); ++i)
    {
      std::vector<PeptideHit>& hits = peptide_ids[i].getHits();
      for (Size h = 0; h < hits.size(); ++h)
      {
        const bool scored = h == 0 && copy_of_id[i] != NO_COPY && !copies[copy_of_id[i]].getHits().empty();
        if (!scored)
        {
          hits[h].removeMetaValue(index);
          continue;
        }
        const double q = copies[copy_of_id[i]].getHits().front().getScore();
        hits[h].setMetaValue(index, q);
        if (target_peptides_1pct != nullptr && !hits[h].isDecoy() && q <= 0.01)
        {
          target_peptides.insert(hits[h].getSequence().toUnmodifiedString());
        }
      }
    }
    if (target_peptides_1pct != nullptr) { *target_peptides_1pct = target_peptides.size(); }
    return true;
  }

  // static
  std::string ProSEAlgorithm::pinFeatureEnzyme(const std::string& enzyme, bool fixed_enzyme_termini)
  {
    // a ProteaseDB enzyme without a Percolator name keeps its name: PercolatorInfile then uses its cleavage rule
    const std::string name = pinEnzymeName(enzyme);
    return (name == "no_enzyme" ? enzyme : name) + (fixed_enzyme_termini ? "+termini" : "+keep_termini");
  }

  // static
  std::string ProSEAlgorithm::pinEnzymeName(const std::string& enzyme)
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

  namespace
  {
  // =====================================================================
  // rescore_reuse: cross-fitted Percolator models, trained once on PSMs of one or more files of a search and applied
  // to the files (static models: Fondrie & Noble, J Proteome Res 2020; joint models: mokapot, Fondrie & Noble, J
  // Proteome Res 2021; training subsampled over runs and scores relative to the decoys: PyProphet, Rosenberger et al.,
  // Nat Methods 2017).
  //
  // How the FDR stays valid: (1) No PSM is scored by a model that was trained on it, and no PSM's label enters the
  // transform of its own score. The spectra are split into two halves by a seeded hash of (file index, spectrum
  // index); one model is trained per half, and a PSM is scored by the model of the other half. All PSMs of a spectrum,
  // target and decoy, are in the same half, so Percolator's target-decoy competition (spectrum and experimental mass,
  // as without rescore_reuse) compares scores of the same model. A file that does not train ('static') is seen by
  // neither model. As in Percolator's cross-validation, the model that scores a PSM was trained on other spectra of
  // the same file (with 'joint' and 'own'): a PSM's label enters only models that score other PSMs. (2) The training
  // spectra are drawn by the hash, independently of labels and scores. (3) The two models' scores are put on one
  // scale, per model: relative to the mean and standard deviation of its scores of calibration decoys, decoys of its
  // own half that it was not trained on (5 % of the spectra of each half are kept out of the training for this);
  // these PSMs are scored by the other model. (4) q-values and PEPs are estimated per file by Percolator::score() on
  // the file's own PSMs; they are not pooled across files. (5) The domain-shift check of 'static' compares the feature
  // distributions of all PSMs of a file with those of the training PSMs, without their target/decoy labels, so the
  // choice between the static models and Percolator trained on the file does not depend on the labels it estimates
  // the FDR from.
  // =====================================================================
    UInt64 reuseMix(UInt64 x)
    {
      x += 0x9e3779b97f4a7c15ULL;
      x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
      x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
      return x ^ (x >> 31);
    }

    // stream 1: the half of a spectrum (top bit); stream 2: its uniform draw for the training sample
    UInt64 reuseHash(UInt64 seed, UInt64 stream, Size file, Size spectrum)
    {
      return reuseMix(reuseMix(reuseMix(seed * 4 + stream) ^ static_cast<UInt64>(file)) ^ static_cast<UInt64>(spectrum));
    }

    int reuseHalf(UInt64 seed, Size file, Size spectrum)
    {
      return static_cast<int>(reuseHash(seed, 1, file, spectrum) >> 63);
    }

    double reuseUniform(UInt64 seed, Size file, Size spectrum)
    {
      return static_cast<double>(reuseHash(seed, 2, file, spectrum) >> 11) * 0x1.0p-53;
    }

    // The raw SVM score of a feature row: weights x features + bias (the last weight)
    double reuseRawScore(const std::vector<double>& weights, const double* features, Size n_features)
    {
      double score = weights.back();
      for (Size j = 0; j < n_features; ++j) score += weights[j] * features[j];
      return score;
    }

    // AUC of @p a against @p b: the probability that a value of a exceeds one of b (ties count half)
    double reuseAUC(std::vector<double> a, std::vector<double> b)
    {
      if (a.empty() || b.empty()) return 0.5;
      std::sort(a.begin(), a.end());
      std::sort(b.begin(), b.end());
      double u = 0.0;
      Size below = 0, not_above = 0; // values of b below x, not above x
      for (const double x : a)
      {
        while (below < b.size() && b[below] < x) ++below;
        not_above = std::max(not_above, below);
        while (not_above < b.size() && b[not_above] <= x) ++not_above;
        u += static_cast<double>(below) + 0.5 * static_cast<double>(not_above - below);
      }
      return u / (static_cast<double>(a.size()) * static_cast<double>(b.size()));
    }
  }

  namespace Internal::ProSEAlgorithmHelpers
  {
    // The Percolator parameters of ProSE's rescoring: @p training (parameters 'rescoring:...', 'rescore_solver'), then
    // @p extra (rescore_method 'reset': the RESET parameters)
    Param rescoringPercolatorParameters(Param pp, int threads, const ProSEAlgorithm::PercolatorTraining& training,
                                        const Param& extra)
    {
      pp.setValue("train_best_positive", "true");
      pp.setValue("post_processing_tdc", "true");
      pp.setValue("subset_max_train", training.subset_max_train);  // 0 (the default): train on all PSMs
      pp.setValue("subset_mode", training.subset_mode);
      pp.setValue("subset_curve_start", training.curve_start);
      pp.setValue("subset_curve_min_gain", training.curve_min_gain);
      pp.setValue("num_iterations", training.iterations);
      pp.setValue("tron_epsilon", training.tron_epsilon);
      pp.setValue("warm_start", training.warm_start ? "true" : "false");
      pp.setValue("single_cneg", training.single_cneg ? "true" : "false");
      pp.setValue("pep_method", "nonparametric");  // as the percolator executable
      // cross-validation trains 9 models at a time (3 folds x 3 Cneg values) and sorts on more threads (up to 12); the
      // results do not depend on the threads
      pp.setValue("num_threads", threads <= 0 ? 3 : std::min(threads, 128));
      pp.setValue("solver", training.solver);
      pp.setValue("use_pi0", training.use_pi0);
      pp.setValue("calc_mass_feature", training.calc_mass_feature ? "true" : "false");
      if (training.peptide_level) pp.setValue("peptide_level", "true");
      for (auto it = extra.begin(); it != extra.end(); ++it)
      {
        pp.setValue(it.getName(), it->value);
      }
      return pp;
    }

    // The enzyme name of the PIN features of the rescoring: Percolator's name of the search's enzyme, with
    // PercolatorInfile's "+termini" modifier for 'rescoring:fixed_enzyme_termini'.
    std::string rescoringEnzyme(const std::string& enzyme, const ProSEAlgorithm::PercolatorTraining& training)
    {
      return ProSEAlgorithm::pinFeatureEnzyme(enzyme, training.fixed_enzyme_termini);
    }
  } // namespace Internal::ProSEAlgorithmHelpers

  namespace
  {
    // Training PSMs kept for the domain-shift check, and PSMs of a file compared with them
    constexpr Size REUSE_SHIFT_PSMS = 20000;
    // The spectra of a half whose uniform draw is at least this are calibration spectra: never trained on
    constexpr double REUSE_CALIBRATION_FROM = 0.95;
    // Calibration decoys used per model (every k-th)
    constexpr Size REUSE_CALIBRATION_DECOYS = 200000;
  }

  bool ProSEAlgorithm::reuseTrainingDraw(UInt64 seed, Size file, Size spectrum, double fraction)
  {
    return reuseUniform(seed, file, spectrum) < REUSE_CALIBRATION_FROM * fraction;
  }

  bool ProSEAlgorithm::reuseCalibrationDraw(UInt64 seed, Size file, Size spectrum, double fraction)
  {
    const double u = reuseUniform(seed, file, spectrum);
    return u >= REUSE_CALIBRATION_FROM && (u - REUSE_CALIBRATION_FROM) / (1.0 - REUSE_CALIBRATION_FROM) < fraction;
  }

  struct ProSEAlgorithm::ReuseModel_
  {
    std::string mode;                    ///< rescore_reuse:mode of the training
    UInt64 seed = 0;                     ///< seed of the random halves and training spectra (Percolator's)
    StringList features;                 ///< numeric features of the models
    int min_charge = 0;                  ///< charge range of the charge features
    int max_charge = 0;
    std::array<PercolatorModel, 2> models; ///< models[h]: trained on half h
    /// mean and standard deviation of models[h]'s scores of the calibration decoys of half h
    std::array<double, 2> center{0.0, 0.0};
    std::array<double, 2> scale{1.0, 1.0};
    std::array<Size, 2> training_psms{0, 0};
    std::vector<Size> training_files;    ///< indices of the files that trained
    std::vector<std::vector<double>> sample; ///< training PSMs (a sample, targets and decoys), for the domain-shift check
    std::vector<Size> used_features;     ///< features with a weight in a model, for the domain-shift check
    bool used_default_direction = false; ///< Percolator's sanity check replaced a learned model

    std::string describe() const
    {
      std::string files;
      for (Size f : training_files) files += (files.empty() ? "" : ",") + StringUtils::toStr(f + 1);
      return mode + " (cross-fitted models trained on " + StringUtils::toStr(training_psms[0]) + " + "
             + StringUtils::toStr(training_psms[1]) + " PSMs of file(s) " + files + ")";
    }
  };

  struct ProSEAlgorithm::ReuseScoring_
  {
    /// what the scorer did (written by it)
    struct Note
    {
      std::string model;
      bool used_default_direction = false;
      double shift = -1.0;
      std::string fallback; ///< why Percolator is to be trained on the run instead (then the scorer throws)
    };
    Percolator::PSMScorer scorer;
    int min_charge = 0;
    int max_charge = 0;
    std::shared_ptr<Note> note;
    std::string mode;
    UInt64 seed = 0;
    Size max_train = 0;
  };

  std::shared_ptr<const ProSEAlgorithm::ReuseModel_> ProSEAlgorithm::trainReuseModel_(const std::vector<ReuseRun_>& runs,
                                                                                       int threads) const
  {
    auto model = std::make_shared<ReuseModel_>();
    model->mode = reuse_mode_;
    {
      // The random halves and training spectra follow Percolator's seed (as its cross-validation folds do).
      Percolator perc;
      perc.setParameters(rescoringPercolatorParameters(perc.getDefaults(), threads, percolator_training_, Param()));
      model->seed = static_cast<UInt64>(static_cast<int>(perc.getParameters().getValue("seed")));
    }
    const UInt64 seed = model->seed;
    // The runs that rescoreWithPercolator_() would rescore: a run, at least 100 PSMs, decoys.
    std::vector<const ReuseRun_*> used;
    int min_charge = 10;
    int max_charge = 0;
    for (const ReuseRun_& run : runs)
    {
      if (run.protein_ids->empty() || run.peptide_ids->size() < 100) continue;
      bool has_decoys = false;
      int run_min_charge = 10;
      int run_max_charge = 0;
      for (const PeptideIdentification& pid : *run.peptide_ids)
      {
        for (const PeptideHit& hit : pid.getHits())
        {
          has_decoys = has_decoys || hit.isDecoy();
          run_min_charge = std::min(run_min_charge, hit.getCharge());
          run_max_charge = std::max(run_max_charge, hit.getCharge());
        }
      }
      if (!has_decoys) continue;
      min_charge = std::min(min_charge, run_min_charge);
      max_charge = std::max(max_charge, run_max_charge);
      used.push_back(&run);
    }
    if (used.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "No run with at least 100 PSMs and decoy PSMs to train the models on.");
    }
    model->min_charge = min_charge;
    model->max_charge = max_charge;
    const ProteinIdentification::SearchParameters& sp = used.front()->protein_ids->front().getSearchParameters();
    const StringList feature_set = PercolatorInfile::getFeatureSet(sp, min_charge, max_charge);
    const std::string enzyme = rescoringEnzyme(sp.digestion_enzyme.getName(), percolator_training_);

    // The training PSMs of a half: max_train shared equally by the files, a file with fewer PSMs gives all of them (the
    // hits of a spectrum count, as an estimate of its feature rows).
    const Size n_runs = used.size();
    std::vector<std::array<double, 2>> psms(n_runs, {0.0, 0.0});
    for (Size k = 0; k < n_runs; ++k)
    {
      const PeptideIdentificationList& ids = *used[k]->peptide_ids;
      for (Size s = 0; s < ids.size(); ++s)
      {
        psms[k][reuseHalf(seed, used[k]->file_index, s)] += static_cast<double>(ids[s].getHits().size());
      }
    }
    std::vector<std::array<double, 2>> fraction(n_runs, {1.0, 1.0});
    if (reuse_max_train_ > 0)
    {
      for (int h = 0; h < 2; ++h)
      {
        std::vector<Size> order(n_runs);
        std::iota(order.begin(), order.end(), Size(0));
        std::stable_sort(order.begin(), order.end(), [&](Size a, Size b) { return psms[a][h] < psms[b][h]; });
        double budget = static_cast<double>(reuse_max_train_);
        Size left = n_runs;
        for (const Size k : order)
        {
          const double quota = budget / static_cast<double>(left--);
          const double trainable = psms[k][h] * REUSE_CALIBRATION_FROM;
          if (trainable > quota) { fraction[k][h] = quota / trainable; }
          budget -= std::min(trainable, quota);
        }
      }
    }
    // The calibration spectra of a half: a share of those with u >= REUSE_CALIBRATION_FROM, for about
    // 2 x REUSE_CALIBRATION_DECOYS PSMs per half.
    std::array<double, 2> calibration_fraction{1.0, 1.0};
    for (int h = 0; h < 2; ++h)
    {
      double half_psms = 0.0;
      for (Size k = 0; k < n_runs; ++k) { half_psms += psms[k][h]; }
      const double calibration_psms = half_psms * (1.0 - REUSE_CALIBRATION_FROM);
      if (calibration_psms > 2.0 * REUSE_CALIBRATION_DECOYS) { calibration_fraction[h] = 2.0 * REUSE_CALIBRATION_DECOYS / calibration_psms; }
    }
    std::array<RescoreInput, 2> train;
    StopWatch sw_rows;
    sw_rows.start();
    std::array<std::vector<std::vector<double>>, 2> calibration; // calibration decoys of each half
    for (Size k = 0; k < n_runs; ++k)
    {
      const Size file = used[k]->file_index;
      // Only the drawn spectra get feature rows: a copy of them (spectra[i]: the index of the copy's spectrum i).
      const PeptideIdentificationList& ids = *used[k]->peptide_ids;
      std::vector<Size> spectra;
      for (Size s = 0; s < ids.size(); ++s)
      {
        const int h = reuseHalf(seed, file, s);
        // (disjoint: a calibration spectrum is never a training spectrum)
        const bool drawn = reuseCalibrationDraw(seed, file, s, calibration_fraction[h])
                           || reuseTrainingDraw(seed, file, s, fraction[k][h]);
        if (drawn && !ids[s].getHits().empty()) { spectra.push_back(s); }
      }
      PeptideIdentificationList drawn_ids;
      drawn_ids.resize(spectra.size());
#pragma omp parallel for schedule(dynamic, 1024) num_threads(std::max(1, threads))
      for (SignedSize i = 0; i < static_cast<SignedSize>(spectra.size()); ++i) { drawn_ids[i] = ids[spectra[i]]; }
      if (drawn_ids.empty()) { continue; }
      RescoreInput ri;
      Percolator::PSMRowLocations rows;
      Percolator::buildPSMRescoreInput(drawn_ids, feature_set, enzyme, min_charge, max_charge, ri, rows,
                                       percolator_training_.calc_mass_feature);
      if (model->features.empty()) { model->features = ri.feature_names; }
      else if (ri.feature_names != model->features)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "The files of the search have different Percolator features.", "");
      }
      // (the rows of buildPSMRescoreInput are dense: the training rows are appended to dense matrices as well)
      const Size n_feat = ri.featureCount();
      for (Size r = 0; r < ri.rowCount(); ++r)
      {
        const Size spectrum = spectra[rows[r].first];
        const int h = reuseHalf(seed, file, spectrum);
        if (reuseUniform(seed, file, spectrum) >= REUSE_CALIBRATION_FROM)
        {
          if (ri.is_decoy[r]) { calibration[h].emplace_back(ri.row(r), ri.row(r) + n_feat); }
          continue;
        }
        RescoreInput& t = train[h];
        t.feature_matrix.insert(t.feature_matrix.end(), ri.row(r), ri.row(r) + n_feat);
        t.is_decoy.push_back(ri.is_decoy[r]);
        // The spectrum (not its scan number, which may be missing): Percolator's cross-validation folds and its best
        // PSM per spectrum go by file and scan number.
        t.scan_numbers.push_back(static_cast<int>(spectrum));
        t.spec_file_numbers.push_back(static_cast<int>(k));
        t.exp_masses.push_back(ri.exp_masses[r]);
        t.calc_masses.push_back(ri.calc_masses[r]);
      }
      model->training_files.push_back(file);
    }
    for (int h = 0; h < 2; ++h)
    {
      train[h].feature_names = model->features;
      model->training_psms[h] = train[h].rowCount();
      if (train[h].rowCount() == 0)
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "No training PSMs in a half of the spectra.");
      }
      Percolator perc;
      perc.setParameters(rescoringPercolatorParameters(perc.getDefaults(), threads, percolator_training_, Param()));
      if (h == 0)
      {
        sw_rows.stop();
        OPENMS_LOG_INFO << "[ProSE] rescore_reuse: training PSMs of " << n_runs << " file(s) drawn in " << sw_rows.getClockTime()
                        << " s" << std::endl;
      }
      StopWatch sw_train;
      sw_train.start();
      model->models[h] = perc.train(train[h]);
      sw_train.stop();
      OPENMS_LOG_INFO << "[ProSE] rescore_reuse: model " << h + 1 << " trained on " << train[h].rowCount() << " PSMs in "
                      << sw_train.getClockTime() << " s" << std::endl;
      model->used_default_direction = model->used_default_direction || perc.usedDefaultDirection();
    }
    // The scale of each model: its scores of the calibration decoys of its half (not trained on; scored by the other model).
    for (int h = 0; h < 2; ++h)
    {
      const std::vector<std::vector<double>>& decoys = calibration[h];
      const Size step = std::max<Size>(1, (decoys.size() + REUSE_CALIBRATION_DECOYS - 1) / REUSE_CALIBRATION_DECOYS);
      double sum = 0.0, sum2 = 0.0;
      Size n = 0;
      for (Size r = 0; r < decoys.size(); r += step)
      {
        const double score = reuseRawScore(model->models[h].weights, decoys[r].data(), decoys[r].size());
        sum += score;
        sum2 += score * score;
        ++n;
      }
      if (n > 1)
      {
        const double mean = sum / static_cast<double>(n);
        const double sd = std::sqrt(std::max(0.0, (sum2 - sum * mean) / static_cast<double>(n - 1)));
        model->center[h] = mean;
        model->scale[h] = (sd > 0.0 && std::isfinite(sd)) ? sd : 1.0;
      }
    }
    // The domain-shift check: features with a weight, and a sample of the training PSMs (every k-th, without labels).
    for (Size j = 0; j < model->features.size(); ++j)
    {
      if (model->models[0].weights[j] != 0.0 || model->models[1].weights[j] != 0.0) { model->used_features.push_back(j); }
    }
    const Size n_train = train[0].rowCount() + train[1].rowCount();
    const Size step = std::max<Size>(1, (n_train + REUSE_SHIFT_PSMS - 1) / REUSE_SHIFT_PSMS);
    Size d = 0;
    for (const RescoreInput& t : train)
    {
      for (Size r = 0; r < t.rowCount(); ++r)
      {
        if (d++ % step == 0) { model->sample.emplace_back(t.row(r), t.row(r) + t.featureCount()); }
      }
    }
    return model;
  }

  ProSEAlgorithm::PercolatorRescoring ProSEAlgorithm::rescoreWithReuseModel_(
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids,
      Size file_index,
      const ReuseModel_& model,
      bool check_shift,
      int threads,
      std::vector<std::vector<PeptideHit>>* dropped) const
  {
    ReuseScoring_ reuse;
    reuse.min_charge = model.min_charge;
    reuse.max_charge = model.max_charge;
    reuse.note = std::make_shared<ReuseScoring_::Note>();
    reuse.mode = model.mode;
    reuse.seed = model.seed;
    reuse.max_train = reuse_max_train_;
    const UInt64 seed = model.seed;
    const double max_shift = reuse_max_shift_;
    std::shared_ptr<ReuseScoring_::Note> note = reuse.note;
    // the reused models are Percolator models with the training options (RESET's parameters do not apply)
    const PercolatorTraining training = percolator_training_;
    reuse.scorer = [&model, note, seed, file_index, check_shift, max_shift, threads, training](
      const RescoreInput& ri, const Percolator::PSMRowLocations& rows) -> RescoreOutput
    {
      Percolator perc;
      perc.setParameters(rescoringPercolatorParameters(perc.getDefaults(), threads, training, Param()));
      // Percolator trained on the run, as without rescore_reuse: rescoreWithReuseModel_() rescores the run again without
      // the models (with the run's own feature set)
      const auto trained_on_run = [&](const std::string& why) -> RescoreOutput
      {
        note->fallback = why;
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "rescore_reuse: Percolator is to be trained on the run (" + why + ")", "");
      };
      if (ri.feature_names != model.features)
      {
        return trained_on_run("its Percolator features differ from those of the models");
      }
      const Size n = ri.rowCount();
      const Size n_feat = ri.featureCount();
      if (check_shift)
      {
        // all PSMs, without their labels (see above)
        const Size step = std::max<Size>(1, (n + REUSE_SHIFT_PSMS - 1) / REUSE_SHIFT_PSMS);
        double shift = 0.0;
        std::string worst;
        for (const Size j : model.used_features)
        {
          std::vector<double> a, b;
          for (Size r = 0; r < n; r += step) { a.push_back(ri.row(r)[j]); }
          for (const std::vector<double>& row : model.sample) { b.push_back(row[j]); }
          const double s = std::fabs(2.0 * reuseAUC(std::move(a), std::move(b)) - 1.0);
          if (s > shift) { shift = s; worst = model.features[j]; }
        }
        note->shift = shift;
        OPENMS_LOG_INFO << "[ProSE] rescore_reuse: domain shift of the PSMs " << shift << " (largest in feature '" << worst
                        << "'; limit " << max_shift << ")" << std::endl;
        if (shift > max_shift)
        {
          return trained_on_run("domain shift " + StringUtils::toStr(shift) + " in feature '" + worst + "'");
        }
      }
      // The cross-fitted score of every row: the model of the other half, on its scale.
      RescoreInput one;
      one.feature_names = {"rescore_reuse_score"};
      one.features.assign(n, std::vector<double>(1, 0.0));
#pragma omp parallel for schedule(static) num_threads(std::max(1, threads))
      for (SignedSize r = 0; r < static_cast<SignedSize>(n); ++r)
      {
        const int m = 1 - reuseHalf(seed, file_index, rows[r].first);
        one.features[r][0] = (reuseRawScore(model.models[m].weights, ri.row(r), n_feat) - model.center[m]) / model.scale[m];
      }
      one.is_decoy = ri.is_decoy;
      {
        double sum_t = 0.0, sum_d = 0.0;
        Size n_d = 0;
        for (Size r = 0; r < n; ++r) { (ri.is_decoy[r] ? sum_d : sum_t) += one.features[r][0]; n_d += ri.is_decoy[r] ? 1 : 0; }
        OPENMS_LOG_DEBUG << "[ProSE] rescore_reuse: " << n << " rows (" << n_d << " decoys) of file " << file_index + 1
                         << ", mean cross-fitted score of the targets " << sum_t / std::max<double>(1.0, double(n - n_d))
                         << ", of the decoys " << sum_d / std::max<double>(1.0, double(n_d)) << std::endl;
      }
      // The spectrum: best PSM per spectrum (target-decoy competition) as in training
      one.scan_numbers.resize(n);
      for (Size r = 0; r < n; ++r) { one.scan_numbers[r] = static_cast<int>(rows[r].first); }
      one.exp_masses = ri.exp_masses;
      one.calc_masses = ri.calc_masses;
      PercolatorModel unit;
      unit.feature_names = one.feature_names;
      unit.weights = {1.0, 0.0};
      try
      {
        RescoreOutput out = perc.score(one, unit);
        note->model = model.describe();
        note->used_default_direction = model.used_default_direction;
        return out;
      }
      catch (const Exception::BaseException& e)
      {
        return trained_on_run(std::string("scoring with the models failed: ") + e.what());
      }
    };
    PercolatorRescoring outcome = rescoreWithPercolator_(protein_ids, peptide_ids, threads, dropped, percolator_training_,
                                                         Param(), &reuse);
    if (note->fallback.empty())
    {
      outcome.model = note->model;
      return outcome;
    }
    // The run is rescored as without rescore_reuse (rescorePSMs() left the PSMs unchanged when the scorer threw).
    OPENMS_LOG_INFO << "[ProSE] rescore_reuse: Percolator is trained on the file: " << note->fallback << std::endl;
    outcome = rescoreWithPercolator_(protein_ids, peptide_ids, threads, dropped, percolator_training_, Param());
    if (outcome.applied)
    {
      outcome.model = "trained on the file (" + note->fallback + ")";
      for (ProteinIdentification& run : protein_ids)
      {
        ProteinIdentification::SearchParameters run_sp = run.getSearchParameters();
        run_sp.setMetaValue("Percolator:rescore_reuse_mode", reuse.mode);
        run_sp.setMetaValue("Percolator:rescore_reuse_model", outcome.model);
        if (note->shift >= 0.0) { run_sp.setMetaValue("Percolator:rescore_reuse_shift", note->shift); }
        run.setSearchParameters(run_sp);
      }
    }
    return outcome;
  }

  // =====================================================================
  // In-process Percolator rescoring of one search run (parameter 'rescore'). Single
  // source of truth for the ProSE TOPP tool and other callers (e.g. the
  // WebAssembly bindings); it was a copy in the tool before.
  // =====================================================================
  // static
  ProSEAlgorithm::PercolatorRescoring ProSEAlgorithm::rescoreWithPercolator(
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids,
      int threads,
      Size max_train_psms)
  {
    PercolatorTraining training;
    training.subset_max_train = static_cast<int>(std::min<Size>(max_train_psms, std::numeric_limits<int>::max()));
    return rescoreWithPercolator_(protein_ids, peptide_ids, threads, nullptr, training, Param());
  }

  // static
  ProSEAlgorithm::PercolatorRescoring ProSEAlgorithm::rescoreWithPercolator_(
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids,
      int threads,
      std::vector<std::vector<PeptideHit>>* dropped,
      const PercolatorTraining& training,
      const Param& extra,
      const ReuseScoring_* reuse)
  {
    PercolatorRescoring outcome;
    if (protein_ids.empty())
    {
      outcome.skipped_reason = "No search run (protein identification) to rescore the PSMs of.";
      return outcome;
    }
    if (peptide_ids.size() < 100)
    {
      outcome.skipped_reason = "Only " + StringUtils::toStr(peptide_ids.size())
                               + " PSMs, Percolator rescoring needs at least 100.";
      return outcome;
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
      // The search's decoy mode is not known here: name both causes.
      outcome.skipped_reason = "No decoy PSMs for the target/decoy competition of Percolator. A search without "
                               "decoys needs decoys 'auto' or 'generate', or a FASTA with decoy proteins; with "
                               "decoys, no spectrum of this run matched one (e.g. a small run).";
      return outcome;
    }

    // The charge features of the models with rescore_reuse (a charge outside their range has none)
    if (reuse != nullptr)
    {
      min_charge = reuse->min_charge;
      max_charge = reuse->max_charge;
    }
    const ProteinIdentification::SearchParameters& sp = protein_ids.front().getSearchParameters();
    const StringList feature_set = PercolatorInfile::getFeatureSet(sp, min_charge, max_charge);

    Percolator perc;
    Param pp = rescoringPercolatorParameters(perc.getDefaults(), threads, training, extra);
    perc.setParameters(pp);
    if (training.subset_max_train > 0 && peptide_ids.size() > static_cast<Size>(training.subset_max_train))
    {
      OPENMS_LOG_INFO << "[ProSE] Percolator trains on " << training.subset_max_train << " of the " << peptide_ids.size()
                      << " PSMs (Search:percolator:max_train_psms) and scores all of them." << std::endl;
    }

    try
    {
      // (only the best hit of each spectrum is kept below: the others need no Percolator meta values)
      if (reuse == nullptr)
      {
        perc.rescorePSMs(peptide_ids, feature_set, rescoringEnzyme(sp.digestion_enzyme.getName(), training),
                         min_charge, max_charge, "q-value", /*best_hit_only=*/true);
      }
      else
      {
        perc.rescorePSMsWithScorer(peptide_ids, feature_set, rescoringEnzyme(sp.digestion_enzyme.getName(), training),
                                   min_charge, max_charge, "q-value", reuse->scorer, /*best_hit_only=*/true);
      }
    }
    catch (const Exception::BaseException& e)
    {
      outcome.skipped_reason = std::string("Percolator failed: ") + e.what();
      return outcome;
    }
    // The best hit per spectrum, as IDFilter::keepNBestHits(peptide_ids, 1), with the spectra in parallel: releasing
    // the other hits (and their meta values) takes most of the time. With the threads of the rescoring, as Percolator's
    // own loops (at most 9 threads, ProSE's former cap on Percolator's threads, released them 1.4 times slower at 16).
#ifdef _OPENMP
    const int max_threads = omp_get_max_threads();
#else
    const int max_threads = 1;
#endif
    [[maybe_unused]] const int trim_threads = std::max(1, std::min(threads <= 0 ? 3 : threads, max_threads));
    // With @p dropped, the other hits are moved there (to be destroyed later, see applyPsmFDR_()) if memory allows,
    // otherwise destroyed here.
    if (dropped != nullptr)
    {
      try { dropped->resize(peptide_ids.size()); }
      catch (...) { dropped = nullptr; }
    }
    std::vector<std::exception_ptr> errors(peptide_ids.size());
#pragma omp parallel for schedule(dynamic, 256) num_threads(trim_threads)
    for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
    {
      try
      {
        PeptideIdentification& pid = peptide_ids[i];
        pid.sort();
        std::vector<PeptideHit>& hits = pid.getHits();
        if (hits.size() > 1)
        {
          if (dropped != nullptr)
          {
            // moves no hit if the allocation fails (then they are destroyed below)
            try { (*dropped)[i].assign(std::make_move_iterator(hits.begin() + 1), std::make_move_iterator(hits.end())); }
            catch (...) {}
          }
          hits.resize(1);
        }
      }
      catch (...)
      {
        errors[i] = std::current_exception();
      }
    }
    for (const std::exception_ptr& error : errors)
    {
      if (error) std::rethrow_exception(error);
    }
    // The peptide-level q-values of the run, 'peptide q-value' (the name FalseDiscoveryRate writes and
    // EntrapmentEvaluator reads): annotatePeptideQValues() on all its PSMs, before FDR:PSM, as on the pooled PSMs of
    // several runs (mergeFileResults()). Not with RESET, which trained on part of the decoys (its own peptide-level
    // q-values are 'RESET:peptide_q_value'); none of an earlier rescoring stays.
    if (perc.getResetOutput() == nullptr) { annotatePeptideQValues(peptide_ids); }
    else
    {
      for (PeptideIdentification& pid : peptide_ids)
      {
        for (PeptideHit& hit : pid.getHits()) { hit.removeMetaValue(Constants::UserParam::PEPTIDE_Q_VALUE); }
      }
    }
    if (reuse != nullptr ? reuse->note->used_default_direction : perc.usedDefaultDirection())
    {
      outcome.warning = "Percolator's sanity check rejected the learned model (on the held-out folds it found no target "
                        "PSM at q < " + StringUtils::toStr(static_cast<double>(pp.getValue("test_fdr")))
                        + ", or fewer than the best single feature; typical for blank runs and runs with few confident "
                          "PSMs), so the PSMs are scored with that feature alone.";
    }

    markRescoredByPercolator_(protein_ids, pp, training, &perc, reuse);
    outcome.applied = true;
    return outcome;
  }

  // static
  void ProSEAlgorithm::markRescoredByPercolator_(std::vector<ProteinIdentification>& protein_ids, const Param& pp,
                                                 const PercolatorTraining& training, const Percolator* perc,
                                                 const ReuseScoring_* reuse)
  {
    // Mark the runs as post-processed by Percolator, as PercolatorAdapter does. Percolator's settings, in the order of
    // the existing output (flags as "true"/"false"); an empty value removes a setting that does not apply.
    std::vector<std::pair<std::string, DataValue>> settings = {
      {"peptide_level_fdrs", "false"},
      {"protein_level_fdrs", "false"},
      {"testFDR", static_cast<double>(pp.getValue("test_fdr"))},
      {"trainFDR", static_cast<double>(pp.getValue("train_fdr"))},
      {"maxiter", static_cast<int>(pp.getValue("num_iterations"))},
      {"subset_max_train", static_cast<int>(pp.getValue("subset_max_train"))}};
    // (only with a cap and an out-of-fold mode, or the learning curve: the output without a cap, or in_sample, stays
    // as it was)
    const std::string subset_mode = pp.getValue("subset_mode").toString();
    const bool record_subset_mode = subset_mode == "learning_curve" ||
                                    (static_cast<int>(pp.getValue("subset_max_train")) > 0 && subset_mode != "in_sample");
    settings.emplace_back("subset_mode", record_subset_mode ? DataValue(subset_mode) : DataValue::EMPTY);
    // the learning curve: "validation=<rows>; <sample rows>:<validation targets per fold>; ...; chosen=<rows per fold>"
    std::string curve;
    if (subset_mode == "learning_curve" && perc != nullptr && !perc->getCappedTrainingInfo().curve_sizes.empty())
    {
      const Percolator::CappedTrainingInfo& capped = perc->getCappedTrainingInfo();
      curve = "validation=" + std::to_string(capped.validation_rows);
      for (size_t k = 0; k < capped.curve_sizes.size(); ++k)
      {
        curve += "; " + std::to_string(capped.curve_sizes[k]) + ":";
        for (size_t f = 0; f < capped.curve_counts[k].size(); ++f)
        {
          curve += (f ? "/" : "") + std::to_string(capped.curve_counts[k][f]);
        }
      }
      curve += "; chosen=";
      for (size_t f = 0; f < capped.chosen_size.size(); ++f) curve += (f ? "/" : "") + std::to_string(capped.chosen_size[f]);
    }
    settings.emplace_back("subset_curve", curve.empty() ? DataValue::EMPTY : DataValue(curve));
    settings.emplace_back("seed", static_cast<int>(pp.getValue("seed")));
    settings.emplace_back("post_processing_tdc", "true");
    settings.emplace_back("train_best_positive", "true");
    settings.emplace_back("solver", training.solver != "tron" ? DataValue(training.solver) : DataValue::EMPTY);
    const auto rr = perc != nullptr ? perc->getResetOutput() : nullptr;
    const bool peptide_q_values = training.peptide_level && reuse == nullptr && perc != nullptr && !rr;
    settings.emplace_back("peptide_q_values", peptide_q_values ? DataValue("percolator_peptide_q_value, percolator_peptide_pep")
                                                               : DataValue::EMPTY);
    if (rr) { settings.emplace_back("algorithm", "reset"); }
    Percolator::markRescoredRuns(protein_ids, "ProSE", settings);

    // what RESET and rescore_reuse did, after them
    if (!rr && reuse == nullptr) { return; }
    for (ProteinIdentification& run : protein_ids)
    {
      ProteinIdentification::SearchParameters& run_sp = run.getSearchParameters();
      if (rr)
      {
        run_sp.setMetaValue("RESET:salt", StringUtils::toStr(rr->salt));
        run_sp.setMetaValue("RESET:protocol", "PSM-only");
        run_sp.setMetaValue("RESET:learner", rr->learner);
        if (!rr->learner_settings.empty()) run_sp.setMetaValue("RESET:learner_settings", rr->learner_settings);
        run_sp.setMetaValue("RESET:fdr_level", "peptide (stem): RESET:peptide_q_value; PSM q-values heuristic");
        run_sp.setMetaValue("RESET:target_stems_1pct", static_cast<int>(rr->discoveries_1pct));
        if (pp.getValue("reset:decoys").toString() == "two")
        {
          run_sp.setMetaValue("RESET:decoys", "two (second set: " + pp.getValue("reset:second_decoy_prefix").toString() + ")");
          run_sp.setMetaValue("RESET:decoy_stems_second_set", static_cast<int>(rr->decoy_stems_second_set));
          run_sp.setMetaValue("RESET:decoy_stems", static_cast<int>(rr->decoy_stems));
          run_sp.setMetaValue("RESET:decoy_stems_both_sets", static_cast<int>(rr->decoy_stems_both_sets));
        }
      }
      if (reuse != nullptr)
      {
        run_sp.setMetaValue("Percolator:rescore_reuse_mode", reuse->mode);
        run_sp.setMetaValue("Percolator:rescore_reuse_model", reuse->note->model);
        run_sp.setMetaValue("Percolator:rescore_reuse_max_train", static_cast<int>(reuse->max_train));
        run_sp.setMetaValue("Percolator:rescore_reuse_seed", static_cast<int>(reuse->seed));
        if (reuse->note->shift >= 0.0) { run_sp.setMetaValue("Percolator:rescore_reuse_shift", reuse->note->shift); }
      }
    }
  }

  void ProSEAlgorithm::collectInferenceOnlyPsms_(const PeptideIdentificationList& peptide_ids, double psm_fdr,
                                                 PeptideIdentificationList& inference_only_psms) const
  {
    const double cutoff = fdr_protein_inference_psm_;
    if (cutoff <= psm_fdr) { return; }
    for (const PeptideIdentification& id : peptide_ids)
    {
      std::vector<PeptideHit> hits;
      for (const PeptideHit& hit : id.getHits())
      {
        if (hit.getScore() > psm_fdr && hit.getScore() <= cutoff) { hits.push_back(hit); }
      }
      if (hits.empty()) { continue; }
      PeptideIdentification& copy = inference_only_psms.emplace_back(id);
      copy.setHits(std::move(hits));
    }
  }

  bool ProSEAlgorithm::applyDeferredPsmFDR(const std::vector<ProteinIdentification>& protein_ids,
                                           PeptideIdentificationList& peptide_ids,
                                           bool rescored,
                                           double psm_fdr,
                                           PeptideIdentificationList* inference_only_psms) const
  {
    if (inference_only_psms != nullptr) { inference_only_psms->clear(); }
    if (psm_fdr <= 0.0 || protein_ids.empty() || peptide_ids.empty())
    {
      return false;
    }
    // Rescored PSMs carry Percolator q-values. The others need q-values from their search engine
    // scores, and so decoys to estimate them with.
    if (!rescored)
    {
      const auto& hits = protein_ids[0].getHits();
      const bool has_decoys = std::any_of(hits.begin(), hits.end(), [](const ProteinHit& ph)
      {
        return ph.metaValueExists("target_decoy") && ph.getMetaValue("target_decoy").toString() == "decoy";
      });
      if (!has_decoys)
      {
        return false;
      }
      annotatePsmQValues(peptide_ids); // FDR:PSM_groups as in the search's own FDR:PSM
    }
    // The PSMs the filter removes but the protein inference uses (FDR:protein_inference_PSM), on the same q-values.
    if (inference_only_psms != nullptr) { collectInferenceOnlyPsms_(peptide_ids, psm_fdr, *inference_only_psms); }
    // Filters target and decoy PSMs alike by q-value; decoy removal is the protein FDR's finalization.
    IDFilter::filterHitsByScore(peptide_ids, psm_fdr);
    return true;
  }

  namespace Internal::ProSEAlgorithmHelpers
  {
  // rescoring:stratum_q_values (ProSEAlgorithm::stratumQValues() in the header)
    // rescoring:stratum_q_values: grouped target-decoy competition (as mokapot's grouped confidence): the rescored best
    // hit of every spectrum belongs to the stratum 'protein_nterm' (an evidence at a protein N-terminus, or only
    // Met-clipped evidences) or 'other'; within each stratum the hits are ordered by percolator_score (ties: decoys
    // first), q = (decoys + 1) / max(targets, 1) from the top, made monotone from the bottom. Meta values
    // 'stratum' and 'stratum_q_value'; the main score stays the pooled q-value.
    void annotateStratumQValues(PeptideIdentificationList& peptide_ids)
    {
      // In parallel (serially ~1 s for a timsTOF run, by the meta value lookups by name): per spectrum, the stratum of
      // its best hit (-1: no hit, 2: not rescored), its score and whether it is a decoy.
      const SignedSize n = static_cast<SignedSize>(peptide_ids.size());
      std::vector<int8_t> stratum_of(n, -1);
      std::vector<double> score_of(n, 0.0);
      std::vector<uint8_t> decoy_of(n, 0);
      MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
      const UInt score_index = registry.getIndex("percolator_score");
      const bool score_registered = score_index != UInt(-1); // (else no hit has the meta value)
      FirstError collect_error;
#pragma omp parallel for schedule(dynamic, 1024)
      for (SignedSize i = 0; i < n; ++i)
      {
        if (peptide_ids[i].getHits().empty()) continue;
        const PeptideHit& hit = peptide_ids[i].getHits()[0];
        if (!score_registered || !hit.metaValueExists(score_index))
        {
          stratum_of[i] = 2;
          continue;
        }
        try
        {
          bool nterm = false, all_clipped = !hit.getPeptideEvidences().empty();
          for (const PeptideEvidence& ev : hit.getPeptideEvidences())
          {
            nterm = nterm || ev.getAABefore() == PeptideEvidence::N_TERMINAL_AA;
            all_clipped = all_clipped && ev.getAABefore() == 'M' && ev.getStart() == 1;
          }
          score_of[i] = static_cast<double>(hit.getMetaValue(score_index));
          decoy_of[i] = hit.isDecoy() ? 1 : 0;
          stratum_of[i] = (nterm || all_clipped) ? 1 : 0;
        }
        catch (...)
        {
          collect_error.set(i);
        }
      }
      if (collect_error.error)
      {
        // as the serial loop: the best hits before the failing one that are not rescored have their q-value
        for (SignedSize i = 0; i < collect_error.index; ++i)
        {
          if (stratum_of[i] == 2) peptide_ids[i].getHits()[0].setMetaValue("stratum_q_value", 1.0);
        }
        collect_error.rethrow();
      }
      // The rows of each stratum in spectrum order. The meta value names are registered in the order in which the
      // serial annotation set them first (MetaInfo orders the values of a hit by registry index): 'stratum_q_value'
      // first if a best hit is not rescored (its q-value was set in a first pass), else 'stratum'.
      std::vector<SignedSize> rows[2];
      bool not_rescored = false;
      for (SignedSize i = 0; i < n; ++i)
      {
        if (stratum_of[i] == 2) not_rescored = true;
        else if (stratum_of[i] >= 0) rows[stratum_of[i]].push_back(i);
      }
      UInt stratum_index = 0, q_index = 0;
      if (not_rescored) q_index = registry.registerName("stratum_q_value");
      if (!rows[0].empty() || !rows[1].empty())
      {
        stratum_index = registry.registerName("stratum");
        q_index = registry.registerName("stratum_q_value");
      }
      FirstError write_error;
      if (not_rescored)
      {
        const DataValue one(1.0); // not rescored (q = 1 in the pooled list as well)
#pragma omp parallel for schedule(static)
        for (SignedSize i = 0; i < n; ++i)
        {
          try
          {
            if (stratum_of[i] == 2) peptide_ids[i].getHits()[0].setMetaValue(q_index, one);
          }
          catch (...)
          {
            write_error.set(i);
          }
        }
        write_error.rethrow();
      }
      static const char* names[2] = {"other", "protein_nterm"};
      for (int g = 0; g < 2; ++g)
      {
        // the rows (in spectrum order) sorted as before: the same stable sort with the same comparison (also for NaN)
        std::vector<SignedSize>& r = rows[g];
        if (std::any_of(r.begin(), r.end(), [&score_of](const SignedSize i) { return std::isnan(score_of[i]); }))
        {
          // With NaN the comparison is no strict weak order, and the permutation depends on the algorithm's path (e.g.
          // on the temporary buffer it gets, counted in elements): the former sort of rows of the former size, exactly
          // (the former Row held a PeptideHit* where this one holds the spectrum index). (Only the buffer that the heap
          // grants can differ, as it could for the former code between runs.)
          struct Row { double score; bool decoy; SignedSize spectrum; };
          static_assert(sizeof(SignedSize) == sizeof(PeptideHit*) && alignof(SignedSize) == alignof(PeptideHit*), "rows of the former size");
          std::vector<Row> former(r.size());
          for (Size i = 0; i < r.size(); ++i) former[i] = {score_of[r[i]], decoy_of[r[i]] != 0, r[i]};
          std::stable_sort(former.begin(), former.end(), [](const Row& a, const Row& b)
                           { return a.score > b.score || (a.score == b.score && a.decoy && !b.decoy); });
          for (Size i = 0; i < r.size(); ++i) r[i] = former[i].spectrum;
        }
        else
        {
          // a strict weak order: every stable sort gives the same (unique) result
          const auto before = [&score_of, &decoy_of](const SignedSize a, const SignedSize b)
                              { return score_of[a] > score_of[b] || (score_of[a] == score_of[b] && decoy_of[a] && !decoy_of[b]); };
#ifdef PROSE_GNU_PARALLEL
          if (r.size() >= 65536 && omp_get_max_threads() > 1) { __gnu_parallel::stable_sort(r.begin(), r.end(), before); }
          else
#endif
          { std::stable_sort(r.begin(), r.end(), before); }
        }
        std::vector<double> q(r.size());
        Size targets = 0, decoys = 0;
        for (Size i = 0; i < r.size(); ++i)
        {
          (decoy_of[r[i]] ? decoys : targets) += 1;
          q[i] = std::min(1.0, static_cast<double>(decoys + 1) / static_cast<double>(std::max<Size>(targets, 1)));
        }
        for (Size i = r.size(); i-- > 1;) q[i - 1] = std::min(q[i - 1], q[i]);
        const DataValue name(names[g]);
#pragma omp parallel for schedule(static)
        for (SignedSize i = 0; i < static_cast<SignedSize>(r.size()); ++i)
        {
          try
          {
            PeptideHit& hit = peptide_ids[r[i]].getHits()[0];
            hit.setMetaValue(stratum_index, name);
            hit.setMetaValue(q_index, q[i]);
          }
          catch (...)
          {
            write_error.set(i);
          }
        }
        write_error.rethrow();
      }
    }
  } // namespace Internal::ProSEAlgorithmHelpers

  void ProSEAlgorithm::stratumQValues(PeptideIdentificationList& peptide_ids)
  {
    annotateStratumQValues(peptide_ids);
  }

  void ProSEAlgorithm::applyRescoredPsmFDR_(std::vector<ProteinIdentification>& protein_ids,
                                            PeptideIdentificationList& peptide_ids,
                                            bool has_decoys,
                                            RunStatistics& stats,
                                            bool rescored,
                                            PeptideIdentificationList* inference_only_psms) const
  {
    if (rescored && stratum_q_values_)
    {
      annotateStratumQValues(peptide_ids);
      for (ProteinIdentification& run : protein_ids) run.getSearchParameters().setMetaValue("rescoring:stratum_q_values", "true");
    }
    // A run that Percolator did not rescore reports at most report:top_hits hits per spectrum, which can be fewer than
    // the rescoring:top_hits that the search kept for Percolator.
    if (!rescored && rescore_top_hits_ > 0) keepTopHits_(peptide_ids, report_top_hits_);
    // The hits per spectrum that Percolator rescored; not recorded with rescoring:top_hits 0 (= report:top_hits, as
    // earlier ProSE versions, which recorded nothing).
    if (rescored && rescore_top_hits_ > 0)
    {
      for (ProteinIdentification& run : protein_ids)
      {
        run.getSearchParameters().setMetaValue("Percolator:top_hits", static_cast<int>(rescoreTopHits_()));
      }
    }
    // PSM-level FDR: Percolator's q-values, or target-decoy q-values of the search engine scores (for a run that
    // Percolator did not rescore, also without a decoy PSM, see annotatePsmQValues()), as applyPsmFDR_() without
    // rescoring. The filter is never lifted: a threshold that the run cannot reach reports no PSM.
    if (inference_only_psms != nullptr) { inference_only_psms->clear(); }
    if (fdr_psm_ > 0.0 && (rescored || has_decoys))
    {
      StopWatch sw_fdr; sw_fdr.start();
      if (!rescored && !peptide_ids.empty()) { annotatePsmQValues(peptide_ids); }
      // The PSMs the filter removes but the protein inference uses (FDR:protein_inference_PSM), on the same q-values.
      if (inference_only_psms != nullptr) { collectInferenceOnlyPsms_(peptide_ids, fdr_psm_, *inference_only_psms); }
      // Filters target and decoy PSMs alike by q-value; decoy removal is the protein FDR's finalization.
      filterHitsByScoreInParallel(peptide_ids, fdr_psm_);
      stats.fdr_applied = true;
      stats.achieved_psm_fdr = maxRetainedScore_(peptide_ids);
      sw_fdr.stop();
      stats.seconds_fdr = sw_fdr.getClockTime();
    }
    else if (fdr_psm_ > 0.0)
    {
      OPENMS_LOG_WARN << "FDR:PSM is set but the search has no decoys (decoys=ignore). "
                         "Use decoys=auto to search/generate decoys and enable FDR. Skipping FDR filtering." << std::endl;
    }
  }

  void ProSEAlgorithm::finishDeferredRescoring_(MultiFileSearchResult& mfres, const std::vector<Size>& deferred) const
  {
    reuse_defer_ = false;
    if (deferred.empty()) { return; }
#ifdef _OPENMP
    const int threads = omp_get_max_threads();
#else
    const int threads = 1;
#endif
    std::vector<ReuseRun_> runs;
    for (const Size d : deferred)
    {
      SearchResult& result = mfres.per_file[d];
      if (result.exit_code == ExitCodes::EXECUTION_OK) { runs.push_back({d, &result.protein_ids, &result.peptide_ids}); }
    }
    std::shared_ptr<const ReuseModel_> model;
    StopWatch sw; sw.start();
    OPENMS_LOG_INFO << "[ProSE] rescore_reuse " << reuse_mode_ << ": training the Percolator models on " << runs.size()
                    << " file(s)..." << std::endl;
    try
    {
      model = trainReuseModel_(runs, threads);
      OPENMS_LOG_INFO << "[ProSE] rescore_reuse: " << model->describe() << " in " << sw.getClockTime() << " s" << std::endl;
    }
    catch (const Exception::BaseException& e)
    {
      OPENMS_LOG_WARN << "[ProSE] rescore_reuse: cannot train the models (" << e.what()
                      << "); every file is rescored by Percolator trained on it, as with 'per_file'." << std::endl;
    }
    // 'static': the files after the training files are rescored with these models after their search (applyPsmFDR_())
    reuse_model_ = reuse_mode_ == "static" ? model : nullptr;
    for (const Size d : deferred)
    {
      SearchResult& result = mfres.per_file[d];
      if (result.exit_code != ExitCodes::EXECUTION_OK) { continue; }
      OPENMS_LOG_INFO << "[ProSE] Rescoring the PSMs of " << result.stats.input_file << " with Percolator..." << std::endl;
      StopWatch sw_file;
      sw_file.start();
      result.percolator = model != nullptr
                          ? rescoreWithReuseModel_(result.protein_ids, result.peptide_ids, d, *model, false, threads, nullptr)
                          : rescoreWithPercolator_(result.protein_ids, result.peptide_ids, threads, nullptr, percolator_training_,
                                                   Param());
      // (the PSMs only the protein inference uses: with FDR:protein, as applyPsmFDR_() collects them)
      applyRescoredPsmFDR_(result.protein_ids, result.peptide_ids, mfres.have_decoys, result.stats, result.percolator.applied,
                           fdr_protein_ > 0.0 ? &result.inference_only_peptide_ids : nullptr);
      // The statistics of the final PSMs (the search collected them before its rescoring); the spectra are gone.
      const Size ms2_spectra = result.stats.ms2_spectra;
      result.stats.charge_histogram.clear();
      result.stats.missed_cleavage_histogram.clear();
      collectRunStatistics_(PeakMap(), result.protein_ids, result.peptide_ids, result.stats);
      result.stats.ms2_spectra = ms2_spectra;
      sw_file.stop();
      OPENMS_LOG_INFO << "[ProSE] rescore_reuse: " << result.stats.input_file << " rescored in " << sw_file.getClockTime() << " s"
                      << std::endl;
    }
    sw.stop();
    OPENMS_LOG_INFO << "[ProSE] rescore_reuse: " << deferred.size() << " file(s) rescored in " << sw.getClockTime() << " s" << std::endl;
  }

  // static
  void ProSEAlgorithm::keepTopHits_(PeptideIdentificationList& peptide_ids, Size n)
  {
    // The hits of a spectrum are in the order of postProcessHits_() (best search engine score first): the first n
    // are its top n. Nothing changes if no spectrum has more (e.g. always with rescoring:top_hits 0).
    if (std::none_of(peptide_ids.begin(), peptide_ids.end(),
                     [n](const PeptideIdentification& pid) { return pid.getHits().size() > n; }))
    {
      return;
    }
#pragma omp parallel for schedule(dynamic, 256)
    for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
    {
      std::vector<PeptideHit>& hits = peptide_ids[i].getHits();
      if (hits.size() > n) hits.erase(hits.begin() + n, hits.end()); // (destroys hits only: cannot throw)
    }
  }

  namespace
  {
    // rescoring:terminal_features: 'protein_nterm' = 1 if an evidence of the hit is at a protein N-terminus; 'met_clipped'
    // = 1 if every evidence is the peptide after an initial Met (aa_before M, start 1: in the search space only through
    // peptide:clip_nterm_methionine). Every hit (targets and decoys); the names join the rescoring features.
    void annotateTerminalFeatures(std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids)
    {
      static const std::string nterm_name = "protein_nterm", clip_name = "met_clipped";
      // In parallel (serially ~0.8 s for a timsTOF run, by the meta value lookups by name), with the names registered
      // as by the first hit of the serial loop (MetaInfo orders the values of a hit by registry index)
      if (std::any_of(peptide_ids.begin(), peptide_ids.end(), [](const PeptideIdentification& pid) { return !pid.getHits().empty(); }))
      {
        MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
        const UInt nterm_index = registry.registerName(nterm_name);
        const UInt clip_index = registry.registerName(clip_name);
        const DataValue one(1.0), zero(0.0);
        FirstError error; // (an allocation failure)
#pragma omp parallel for schedule(dynamic, 1024)
        for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
        {
          try
          {
            for (PeptideHit& hit : peptide_ids[i].getHits())
            {
              bool nterm = false, all_clipped = !hit.getPeptideEvidences().empty();
              for (const PeptideEvidence& ev : hit.getPeptideEvidences())
              {
                nterm = nterm || ev.getAABefore() == PeptideEvidence::N_TERMINAL_AA;
                all_clipped = all_clipped && ev.getAABefore() == 'M' && ev.getStart() == 1;
              }
              hit.setMetaValue(nterm_index, nterm ? one : zero);
              hit.setMetaValue(clip_index, (!nterm && all_clipped) ? one : zero);
            }
          }
          catch (...)
          {
            error.set(i);
          }
        }
        error.rethrow();
      }
      for (ProteinIdentification& prot : protein_ids)
      {
        ProteinIdentification::SearchParameters sp = prot.getSearchParameters();
        StringList features;
        if (sp.metaValueExists("extra_features"))
        {
          const std::string existing = sp.getMetaValue("extra_features").toString();
          if (!existing.empty()) features = ListUtils::create<std::string>(existing);
        }
        for (const std::string& name : {nterm_name, clip_name})
        {
          if (std::find(features.begin(), features.end(), name) == features.end()) features.push_back(name);
        }
        sp.setMetaValue("extra_features", ListUtils::concatenate(features, ","));
        sp.setMetaValue("rescoring:terminal_features", "true");
        prot.setSearchParameters(std::move(sp));
      }
    }
  }

  ProSEAlgorithm::PercolatorRescoring ProSEAlgorithm::applyPsmFDR_(std::vector<ProteinIdentification>& protein_ids,
                                                                   PeptideIdentificationList& peptide_ids,
                                                                   bool has_decoys,
                                                                   RunStatistics& stats,
                                                                   PeptideIdentificationList* pre_fdr_psms,
                                                                   PeptideIdentificationList* inference_only_psms) const
  {
    PercolatorRescoring rescoring;
    // The PSMs between FDR:PSM and FDR:protein_inference_PSM are only needed for the protein FDR.
    if (inference_only_psms != nullptr) { inference_only_psms->clear(); }
    PeptideIdentificationList* inference_only = fdr_protein_ > 0.0 ? inference_only_psms : nullptr;
    if (terminal_features_) annotateTerminalFeatures(protein_ids, peptide_ids); // also for a .pin without rescoring
    // The copy for a Percolator .pin (setKeepPreFdrPsms()) holds the search engine scores, every reported hit and the
    // decoys: it is the input of a later rescoring, so it is taken before the rescoring as well as before FDR:PSM.
    if (pre_fdr_psms != nullptr)
    {
      // (with rescoring:top_hits above report:top_hits the search kept more hits, for Percolator). A multi-file search
      // that may rescore its files together keeps the hits Percolator rescores (as well, for a .pin);
      // rescoreFilesTogether_() uses the first rescoring:top_hits of them and cuts the .pin copy to report:top_hits.
      const Size kept_hits = rescore_top_hits_ == 0 ? std::numeric_limits<Size>::max()
                             : !pool_rescoring_ ? report_top_hits_
                             : keep_pre_fdr_psms_ ? std::max(report_top_hits_, rescoreTopHits_()) : rescoreTopHits_();
      // in parallel (a HeLa run: 540,000 identifications)
      pre_fdr_psms->clear();
      pre_fdr_psms->resize(peptide_ids.size());
      std::vector<std::exception_ptr> errors(peptide_ids.size());
#pragma omp parallel for schedule(dynamic, 1024)
      for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
      {
        try
        {
          PeptideIdentification& copy = (*pre_fdr_psms)[i];
          copy = peptide_ids[i];
          std::vector<PeptideHit>& hits = copy.getHits();
          if (hits.size() > kept_hits) { hits.erase(hits.begin() + kept_hits, hits.end()); }
        }
        catch (...)
        {
          errors[i] = std::current_exception();
        }
      }
      for (const std::exception_ptr& error : errors)
      {
        if (error) { std::rethrow_exception(error); }
      }
    }
    if (rescore_)
    {
      // Percolator rescores the first rescoring:top_hits hits of each spectrum (the search kept more for the .pin copy).
      // Not with decoys 'ignore': the search kept report:top_hits, and Percolator skips the run (no decoys). Nothing is
      // cut with rescoring:top_hits 0 (identifications that the caller's list held before the search included).
      if (rescore_top_hits_ > 0 && decoy_mode_ != DecoyMode_::IGNORE) keepTopHits_(peptide_ids, rescoreTopHits_());
      // Percolator gets all target and decoy PSMs; the PSM FDR follows on its q-values (or, if it
      // skipped the run, on the search engine scores).
      startDeferredRead(); // the next file of a multi-file search, with threads Percolator leaves idle
      if (reuse_defer_)
      {
        // rescore_reuse joint/static: rescored (and FDR:PSM applied) by finishDeferredRescoring_() once the models are trained
        rescoring.skipped_reason = "The rescoring waits for the models of the multi-file search (rescore_reuse).";
        return rescoring;
      }
#ifdef _OPENMP
      const int threads = omp_get_max_threads();
#else
      const int threads = 1;
#endif
      OPENMS_LOG_INFO << "[ProSE] Rescoring the PSMs with Percolator..." << std::endl;
      // The hits that the rescoring drops (all but the best per spectrum) are destroyed on a helper thread while the
      // search and its output continue: destroying the 65,000 dropped hits of an Astral run took the rescoring's 16
      // threads 42-60 ms (allocator contention; 150 ms on one thread) while the search waited. The phases after the
      // rescoring keep 11-14 of 16 threads busy. The helper is joined when the holder goes: at the next rescoring of
      // this object, or with this object (its last copy). Not on Windows, where a holder destroyed while a DLL is
      // unloaded would join under the loader lock, which the ending helper thread may need.
      released_hits_.reset(); // the release of the search before (multi-file search)
#ifdef OPENMS_WINDOWSPLATFORM
      constexpr bool release_later = false;
#else
      const bool release_later = threads > 1 && !helpersOversubscribe(); // (inline, in parallel, when a helper's team would exceed the CPUs)
#endif
      std::vector<std::vector<PeptideHit>> dropped;
      // rescore_method 'reset' with a directory as reset:stem_table: one table per run, named after its spectra file
      Param extra = percolator_extra_;
      if (extra.exists("reset:stem_table"))
      {
        std::string path = extra.getValue("reset:stem_table").toString();
        if (!path.empty() && (path.back() == '/' || File::isDirectory(path)))
        {
          std::string base = reset_run_name_.empty() ? std::string("run") : File::basename(reset_run_name_);
          if (const auto dot = base.rfind('.'); dot != std::string::npos && dot > 0) base.erase(dot);
          extra.setValue("reset:stem_table", (path.back() == '/' ? path : path + "/") + base + ".reset_stems.tsv");
        }
      }
      if (reuse_model_ != nullptr && reuse_mode_ == "static")
      {
        // rescore_reuse static: a file after the training files
        rescoring = rescoreWithReuseModel_(protein_ids, peptide_ids, reuse_file_index_, *reuse_model_, true, threads,
                                           release_later ? &dropped : nullptr);
      }
      else if (reuse_mode_ == "own")
      {
        std::shared_ptr<const ReuseModel_> own;
        try
        {
          own = trainReuseModel_({ReuseRun_{reuse_file_index_, &protein_ids, &peptide_ids}}, threads);
        }
        catch (const Exception::BaseException& e)
        {
          OPENMS_LOG_WARN << "[ProSE] rescore_reuse: cannot train the models (" << e.what()
                          << "); Percolator is trained on the run as with 'per_file'." << std::endl;
        }
        rescoring = own != nullptr
                    ? rescoreWithReuseModel_(protein_ids, peptide_ids, reuse_file_index_, *own, false, threads,
                                             release_later ? &dropped : nullptr)
                    : rescoreWithPercolator_(protein_ids, peptide_ids, threads, release_later ? &dropped : nullptr,
                                             percolator_training_, Param());
      }
      else
      {
        rescoring = rescoreWithPercolator_(protein_ids, peptide_ids, threads, release_later ? &dropped : nullptr,
                                           percolator_training_, extra);
      }
      if (!dropped.empty())
      {
        try
        {
          auto release = std::make_shared<HitRelease>();
          // With a quarter of the threads (on one thread, the 2 million hits that the rescoring of a timsTOF run drops took
          // 2.5 s, longer than the rest of the search: the holder then waited for it)
          release->thread = std::thread([hits = std::move(dropped), release_threads = std::max(1, threads / 4)]() mutable
          {
#ifdef _OPENMP
            omp_set_num_threads(release_threads);
#endif
            (void)release_threads;
            freeInParallel(hits);
          });
          released_hits_ = std::move(release);
        }
        catch (...) {} // if the thread cannot be started, the hits are destroyed in this function
      }
      takeBackThreadsOfEndedRead();
      applyRescoredPsmFDR_(protein_ids, peptide_ids, has_decoys, stats, rescoring.applied, inference_only);
      // A file rescored on its own needs of its copy only what the model for the files too small for one of their own
      // may train on (rescoreFilesTogether_()): its sample, not the whole copy until the search ends.
      if (pre_fdr_psms != nullptr && pool_rescoring_ && !keep_pre_fdr_psms_ && rescoring.applied)
      {
        const std::vector<Size> rows = transferSample_(*pre_fdr_psms, transferQuota_(),
                                                       rescore_top_hits_ > 0 ? rescoreTopHits_() : std::numeric_limits<Size>::max());
        if (rows.size() < pre_fdr_psms->size())
        {
          PeptideIdentificationList sample;
          sample.reserve(rows.size());
          for (Size r : rows) { sample.push_back(std::move((*pre_fdr_psms)[r])); }
          // the others released in parallel (destroys only: cannot throw)
#pragma omp parallel for schedule(dynamic, 1024)
          for (SignedSize i = 0; i < static_cast<SignedSize>(pre_fdr_psms->size()); ++i)
          {
            PeptideIdentification().swap((*pre_fdr_psms)[i]);
          }
          *pre_fdr_psms = std::move(sample);
        }
      }
      return rescoring;
    }

    if (fdr_psm_ > 0.0 && (rescoring.applied || has_decoys))
    {
      // PSM-level FDR: Percolator's q-values, or target-decoy q-values of the search engine scores (also for a run
      // that Percolator did not rescore, and for one without a decoy PSM, see annotatePsmQValues()); target and
      // decoy PSMs are filtered alike by the threshold. No decoy-specific stripping happens here (decoupled from
      // decoy removal); the decoys that pass are kept because downstream protein-level FDR and cross-file merging
      // need them. Categorical decoy removal happens only at protein-FDR finalization (applyCompleteSetProteinFDR()).
      // The filter is never lifted: a threshold that the run cannot reach reports no PSM.
      StopWatch sw_fdr; sw_fdr.start();
      if (!rescoring.applied && !peptide_ids.empty()) { annotatePsmQValues(peptide_ids); }
      if (inference_only != nullptr) { collectInferenceOnlyPsms_(peptide_ids, fdr_psm_, *inference_only); }
      filterHitsByScoreInParallel(peptide_ids, fdr_psm_);
      stats.fdr_applied = true;
      stats.achieved_psm_fdr = maxRetainedScore_(peptide_ids);
      sw_fdr.stop();
      stats.seconds_fdr = sw_fdr.getClockTime();
    }
    else if (fdr_psm_ > 0.0)
    {
      OPENMS_LOG_WARN << "FDR:PSM is set but the search has no decoys (decoys=ignore). "
                         "Use decoys=auto to search/generate decoys and enable FDR. Skipping FDR filtering." << std::endl;
    }
    // The run keeps the proteins PeptideIndexing added (those of all hits, also of the PSMs the filters removed): the
    // protein FDR reduces them to those of the reported PSMs (applyCompleteSetProteinFDR()); without it, an output
    // holds the run's proteins as the search produced them.
    return rescoring;
  }
} // namespace OpenMS
