// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/METADATA/AnnotatedMSRun.h>

#include <algorithm>
#include <numeric>

using namespace std;

namespace OpenMS
{

  MapAlignmentAlgorithmIdentification::MapAlignmentAlgorithmIdentification() :
    DefaultParamHandler("MapAlignmentAlgorithmIdentification"),
    ProgressLogger(), reference_index_(-1), reference_(), min_run_occur_(0), min_score_(0.)
  {
    defaults_.setValue("score_type", "", "Name of the score type to use for ranking and filtering (.oms input only). If left empty, a score type is picked automatically.");

    defaults_.setValue("score_cutoff", "false", "Use only IDs above a score cut-off (parameter 'min_score') for alignment?");
    defaults_.setValidStrings("score_cutoff", {"true", "false"});

    defaults_.setValue("min_score", 0.05, "If 'score_cutoff' is 'true': Minimum score for an ID to be considered.\nUnless you have very few runs or identifications, increase this value to focus on more informative peptides.");

    defaults_.setValue("min_run_occur", 2, "Minimum number of runs (incl. reference, if any) in which a peptide must occur to be used for the alignment.\nUnless you have very few runs or identifications, increase this value to focus on more informative peptides.");
    defaults_.setMinInt("min_run_occur", 2);

    defaults_.setValue("max_rt_shift", 0.5, "Maximum realistic RT difference for a peptide (median per run vs. reference). Peptides with higher shifts (outliers) are not used to compute the alignment.\nIf 0, no limit (disable filter); if > 1, the final value in seconds; if <= 1, taken as a fraction of the range of the reference RT scale.");
    defaults_.setMinFloat("max_rt_shift", 0.0);

    defaults_.setValue("use_unassigned_peptides", "true", "Should unassigned peptide identifications be used when computing an alignment of feature or consensus maps? If 'false', only peptide IDs assigned to features will be used.");
    defaults_.setValidStrings("use_unassigned_peptides", {"true", "false"});

    defaults_.setValue("use_feature_rt", "false", "When aligning feature or consensus maps, don't use the retention time of a peptide identification directly; instead, use the retention time of the centroid of the feature (apex of the elution profile) that the peptide was matched to. If different identifications are matched to one feature, only the peptide closest to the centroid in RT is used.\nPrecludes 'use_unassigned_peptides'.");
    defaults_.setValidStrings("use_feature_rt", {"true", "false"});

    defaults_.setValue("use_adducts", "true", "If IDs contain adducts, treat differently adducted variants of the same molecule as different.");
    defaults_.setValidStrings("use_adducts", {"true", "false"});

    defaults_.setValue("auto_reference", "best_run", "Reference to align to if none is given (neither a reference file nor an input index): 'best_run' - the input that shares the most identified sequences with every other input (on ties, the one with the most identified sequences). A consensus is used instead if no input shares at least two sequences with every other input. If the chosen input leaves other inputs with too few alignment points, other inputs and a consensus are tried as well (see 'auto_reference_min_points'). 'consensus' - median RTs per sequence over all inputs. A consensus favors none of the inputs, but only partly corrects larger RT shifts, because every input contributes to the consensus it is aligned to.");
    defaults_.setValidStrings("auto_reference", {"best_run", "consensus"});

    defaults_.setValue("auto_reference_min_points", 11, "If 'auto_reference' is 'best_run': number of alignment points (after removing outliers, see 'max_rt_shift') that the reference should provide for every other input. If the chosen input leaves inputs with fewer points, and one of them shares at least this many sequences with other inputs, every input and a consensus of all inputs are tried as the reference. The choice that gives the most inputs at least this many points is used (the reference counts); on ties, the first choice is kept, and an input is preferred over a consensus. The default is the smallest number of points to which ProteomicsLFQ and MS1LabeledWorkflow fit an RT model. 0 disables the check.");
    defaults_.setMinInt("auto_reference_min_points", 0);

    defaultsToParam_();
  }

  MapAlignmentAlgorithmIdentification::~MapAlignmentAlgorithmIdentification() = default;

  void MapAlignmentAlgorithmIdentification::checkParameters_(Size runs)
  {
    min_run_occur_ = (int)param_.getValue("min_run_occur");

    // reference is not counted as a regular run:
    if (!reference_.empty()) runs++;

    use_feature_rt_ = param_.getValue("use_feature_rt").toBool();
    if (min_run_occur_ > runs)
    {
      std::string msg = "Warning: Value of parameter 'min_run_occur' (here: " +
        StringUtils::toStr(min_run_occur_) + ") is higher than the number of runs incl. "
        "reference (here: " + StringUtils::toStr(runs) + "). Using " + StringUtils::toStr(runs) +
        " instead.";
      OPENMS_LOG_WARN << msg << endl;
      min_run_occur_ = runs;
    }
    score_cutoff_ = param_.getValue("score_cutoff").toBool();
    // score type may have been set by reference already - don't overwrite it:
    if (score_cutoff_ && score_type_.empty())
    {
      score_type_ = StringUtils::toStr(param_.getValue("score_type"));
    }
    min_score_ = param_.getValue("min_score");
    use_adducts_ = param_.getValue("use_adducts").toBool();
    consensus_reference_ = (param_.getValue("auto_reference").toString() == "consensus");
    auto_reference_min_points_ = Size(int(param_.getValue("auto_reference_min_points")));
}

  Int MapAlignmentAlgorithmIdentification::selectReference_(const vector<SeqToList>& rt_data) const
  {
    if (rt_data.empty()) return -1;

    Size best = 0;
    if (rt_data.size() > 1)
    {
      // only sequences that occur in at least "min_run_occur" inputs are used for the alignment;
      // number them, so that the inputs can be compared quickly:
      map<std::string, Size> n_inputs;
      for (const SeqToList& input : rt_data)
      {
        for (const auto& entry : input) ++n_inputs[entry.first];
      }
      map<std::string, Size> seq_index;
      for (const auto& entry : n_inputs)
      {
        if (entry.second >= max(min_run_occur_, Size(2)))
        {
          seq_index.emplace_hint(seq_index.end(), entry.first, seq_index.size());
        }
      }
      vector<vector<Size>> usable(rt_data.size()); // sorted, like "seq_index"
      for (Size i = 0; i < rt_data.size(); ++i)
      {
        for (const auto& entry : rt_data[i])
        {
          auto pos = seq_index.find(entry.first);
          if (pos != seq_index.end()) usable[i].push_back(pos->second);
        }
      }
      auto n_shared = [&usable](Size i, Size j)
      {
        Size count = 0;
        auto it_i = usable[i].begin(), it_j = usable[j].begin();
        while ((it_i != usable[i].end()) && (it_j != usable[j].end()))
        {
          if (*it_i < *it_j) ++it_i;
          else if (*it_j < *it_i) ++it_j;
          else { ++count; ++it_i; ++it_j; }
        }
        return count;
      };

      // Every other input is aligned using the sequences it shares with the reference, so the
      // reference is the input whose smallest number of shared sequences with any other input
      // is largest. On ties, the input with the most identified sequences wins (then the first).
      vector<Size> candidates(rt_data.size());
      iota(candidates.begin(), candidates.end(), 0);
      stable_sort(candidates.begin(), candidates.end(),
                  [&rt_data](Size i, Size j) { return rt_data[i].size() > rt_data[j].size(); });
      // inputs with few usable sequences limit the overlap most, so compare with those first:
      vector<Size> others(rt_data.size());
      iota(others.begin(), others.end(), 0);
      stable_sort(others.begin(), others.end(),
                  [&usable](Size i, Size j) { return usable[i].size() < usable[j].size(); });
      Size best_overlap = 0;
      bool found = false;
      for (Size candidate : candidates)
      {
        if (found && (usable[candidate].size() <= best_overlap)) continue; // can't do better
        Size overlap = numeric_limits<Size>::max();
        for (Size other : others)
        {
          if (other == candidate) continue;
          overlap = min(overlap, n_shared(candidate, other));
          if (found && (overlap <= best_overlap)) break;
        }
        if (!found || (overlap > best_overlap))
        {
          best = candidate;
          best_overlap = overlap;
          found = true;
        }
      }

      if (best_overlap < 2) // too few for any transformation model
      {
        OPENMS_LOG_WARN << "No reference given, and no input shares at least two identified "
                        << "sequences with every other input - aligning to a consensus of all "
                        << "inputs instead." << endl;
        return -1;
      }
      OPENMS_LOG_INFO << "No reference given - aligning to input " << best + 1
                      << ", which shares at least " << best_overlap
                      << " identified sequences with every other input." << endl;
    }

    return Int(best);
  }

  void MapAlignmentAlgorithmIdentification::alignToInput_(
    vector<SeqToList>& rt_data, Size index, vector<TransformationDescription>& transforms,
    bool sorted, bool verbose)
  {
    SeqToList ref_data;
    ref_data.swap(rt_data[index]);
    rt_data.erase(rt_data.begin() + index);
    reference_index_ = Int(index);
    computeMedians_(ref_data, reference_, sorted);
    computeTransformations_(rt_data, transforms, sorted, verbose);
    reference_.clear(); // taken from the inputs, so it must not carry over into the next call
    rt_data.insert(rt_data.begin() + index, SeqToList());
    rt_data[index].swap(ref_data);
  }

  void MapAlignmentAlgorithmIdentification::alignToAutoReference_(
    vector<SeqToList>& rt_data, vector<TransformationDescription>& transforms, bool sorted)
  {
    Int first = consensus_reference_ ? -1 : selectReference_(rt_data);
    if (first < 0) // align to a consensus of all inputs
    {
      computeTransformations_(rt_data, transforms, sorted);
      return;
    }
    alignToInput_(rt_data, first, transforms, sorted, true); // (all RT lists are sorted now)

    const Size min_points = auto_reference_min_points_;
    Size n_short = 0;
    for (Size i = 0; i < transforms.size(); ++i)
    {
      if ((Int(i) != first) && (transforms[i].getDataPoints().size() < min_points)) ++n_short;
    }
    if (n_short == 0) return;

    // an input can't get more alignment points than it has sequences that occur in other inputs
    // (at least "min_run_occur" inputs in total) - is there an input that could get enough?
    map<std::string, Size> n_inputs;
    for (const SeqToList& input : rt_data)
    {
      for (const auto& entry : input) ++n_inputs[entry.first];
    }
    const Size min_occur = max(min_run_occur_, Size(2));
    bool recoverable = false;
    for (Size i = 0; (i < transforms.size()) && !recoverable; ++i)
    {
      if ((Int(i) == first) || (transforms[i].getDataPoints().size() >= min_points)) continue;
      Size usable = count_if(rt_data[i].begin(), rt_data[i].end(), [&](const auto& entry)
                             { return n_inputs[entry.first] >= min_occur; });
      recoverable = (usable >= min_points);
    }
    OPENMS_LOG_INFO << n_short << " input(s) get fewer than 'auto_reference_min_points' ("
                    << min_points << ") alignment points from input " << first + 1;
    if (!recoverable)
    {
      OPENMS_LOG_INFO << ", but share too few identified sequences with the other inputs for any "
                      << "other choice to give them that many." << endl;
      return;
    }
    OPENMS_LOG_INFO << " - trying every input as reference, and a consensus of all inputs."
                    << endl;

    // (number of inputs with at least "min_points" alignment points - the reference counts -,
    // smallest number of points of the other inputs among them):
    auto assess = [min_points](const vector<TransformationDescription>& trafos, Int ref)
    {
      pair<Size, Size> score(0, numeric_limits<Size>::max());
      for (Size i = 0; i < trafos.size(); ++i)
      {
        Size n_points = trafos[i].getDataPoints().size();
        if (Int(i) == ref) ++score.first;
        else if (n_points >= min_points)
        {
          ++score.first;
          score.second = min(score.second, n_points);
        }
      }
      return score;
    };
    // the first choice is only replaced by one that gives more inputs enough points;
    // among those, other inputs are preferred over a consensus:
    const pair<Size, Size> first_score = assess(transforms, first);
    Int choice = first; // -1 for a consensus
    pair<Size, Size> choice_score = first_score;
    vector<TransformationDescription> trial;
    for (Size i = 0; i < rt_data.size(); ++i)
    {
      if (Int(i) == first) continue;
      alignToInput_(rt_data, i, trial, true, false);
      pair<Size, Size> score = assess(trial, i);
      if ((score.first > first_score.first) && ((choice == first) || (score > choice_score)))
      {
        choice = Int(i);
        choice_score = score;
      }
    }
    reference_index_ = -1;
    computeTransformations_(rt_data, trial, true, false);
    pair<Size, Size> consensus_score = assess(trial, -1);
    if ((consensus_score.first > first_score.first) &&
        ((choice == first) || (consensus_score.first > choice_score.first)))
    {
      choice = -1;
      choice_score = consensus_score;
    }

    if (choice == first)
    {
      OPENMS_LOG_INFO << "Keeping input " << first + 1 << " as reference: no other choice gives "
                      << "more inputs at least " << min_points << " alignment points." << endl;
      reference_index_ = first;
      return;
    }
    if (choice >= 0)
    {
      OPENMS_LOG_INFO << "Aligning to input " << choice + 1 << " instead of input " << first + 1
                      << ": it gives " << choice_score.first << " of " << rt_data.size()
                      << " inputs at least " << min_points << " alignment points (input "
                      << first + 1 << ": " << first_score.first << ")." << endl;
      alignToInput_(rt_data, choice, transforms, true, true);
    }
    else
    {
      OPENMS_LOG_WARN << "Aligning to a consensus of all inputs instead of input " << first + 1
                      << ": it gives " << choice_score.first << " of " << rt_data.size()
                      << " inputs at least " << min_points << " alignment points (input "
                      << first + 1 << ": " << first_score.first << ")." << endl;
      computeTransformations_(rt_data, transforms, true, true);
    }
  }

  // RT lists in "rt_data" will be sorted (unless "sorted" is true)
  void MapAlignmentAlgorithmIdentification::computeMedians_(SeqToList& rt_data,
                                                            SeqToValue& medians,
                                                            bool sorted)
  {
    medians.clear();
    for (SeqToList::iterator rt_it = rt_data.begin();
         rt_it != rt_data.end(); ++rt_it)
    {
      double median = Math::median(rt_it->second.begin(),
                                   rt_it->second.end(), sorted);
      medians.insert(medians.end(), make_pair(rt_it->first, median));
    }
  }

  // lists of peptide hits in "peptides" will be sorted
  bool MapAlignmentAlgorithmIdentification::getRetentionTimes_(
      const PeptideIdentificationList& peptides, SeqToList& rt_data)
  {
    for (auto pep_it = peptides.cbegin(); pep_it != peptides.cend(); ++pep_it)
    {
      if (!pep_it->getHits().empty())
      {
        const PeptideHit* best_hit = getBestScoringHit(pep_it->getHits(), pep_it->isHigherScoreBetter());
        if (better_(best_hit->getScore(), min_score_))
        {
          const std::string& seq = best_hit->getSequence().toString();
          rt_data[seq].push_back(pep_it->getRT());
        }
      }
    }
    return false;
  }

  LegacyIdentificationData::ScoreTypeRef
  MapAlignmentAlgorithmIdentification::handleIdDataScoreType_(const LegacyIdentificationData& id_data)
  {
    LegacyIdentificationData::ScoreTypeRef score_ref;
    if (score_type_.empty()) // choose a score type
    {
      score_ref = id_data.pickScoreType(id_data.getObservationMatches());
      if (score_ref == id_data.getScoreTypes().end())
      {
        std::string msg = "no scores found";
        throw Exception::MissingInformation(__FILE__, __LINE__,
                                            OPENMS_PRETTY_FUNCTION, msg);
      }
      score_type_ = score_ref->cv_term.getName();
      OPENMS_LOG_INFO << "Using score type: " << score_type_ << endl;
    }
    else
    {
      score_ref = id_data.findScoreType(score_type_);
      if (score_ref == id_data.getScoreTypes().end())
      {
        std::string msg = "score type '" + score_type_ + "' not found";
        throw Exception::MissingInformation(__FILE__, __LINE__,
                                            OPENMS_PRETTY_FUNCTION, msg);
      }
    }
    return score_ref;
  }


  bool MapAlignmentAlgorithmIdentification::getRetentionTimes_(
    const LegacyIdentificationData& id_data, SeqToList& rt_data)
  {
    // @TODO: should this get handled as an error?
    if (id_data.getObservationMatches().empty()) return true;

    LegacyIdentificationData::ScoreTypeRef score_ref =
      handleIdDataScoreType_(id_data);

    vector<LegacyIdentificationData::ObservationMatchRef> top_hits =
      id_data.getBestMatchPerObservation(score_ref);

    for (const auto& hit : top_hits)
    {
      bool include = true;
      if (score_cutoff_)
      {
        pair<double, bool> result = hit->getScore(score_ref);
        if (!result.second ||
            score_ref->isBetterScore(min_score_, result.first))
        {
          include = false;
        }
      }
      if (include)
      {
        std::string molecule = hit->identified_molecule_var.toString();
        if (use_adducts_ && hit->adduct_opt)
        {
          molecule += "+[" + (*hit->adduct_opt)->getName() + "]";
        }
        rt_data[molecule].push_back(hit->observation_ref->rt);
      }
    }
    return false;
  }

  void MapAlignmentAlgorithmIdentification::computeTransformations_(
    vector<SeqToList>& rt_data, vector<TransformationDescription>& transforms,
    bool sorted, bool verbose)
  {
    Int size = rt_data.size(); // not Size because we compare to Ints later
    transforms.clear();

    // filter RT data (remove peptides that elute in several fractions):
    // TODO

    // compute RT medians:
    OPENMS_LOG_DEBUG << "Computing RT medians..." << endl;
    vector<SeqToValue> medians_per_run(size);
    for (Int i = 0; i < size; ++i)
    {
      computeMedians_(rt_data[i], medians_per_run[i], sorted);
    }
    SeqToList medians_per_seq;
    for (vector<SeqToValue>::iterator run_it = medians_per_run.begin();
         run_it != medians_per_run.end(); ++run_it)
    {
      for (SeqToValue::iterator med_it = run_it->begin();
           med_it != run_it->end(); ++med_it)
      {
        medians_per_seq[med_it->first].push_back(med_it->second);
      }
    }

    // get reference retention time scale: either directly from reference file,
    // or compute consensus time scale
    bool reference_given = !reference_.empty(); // reference file given
    if (reference_given)
    {
      // remove peptides that don't occur in enough runs:
      OPENMS_LOG_DEBUG << "Removing peptides that occur in too few runs..." << endl;
      SeqToValue temp;
      for (SeqToValue::iterator ref_it = reference_.begin();
           ref_it != reference_.end(); ++ref_it)
      {
        SeqToList::iterator med_it = medians_per_seq.find(ref_it->first);
        if ((med_it != medians_per_seq.end()) &&
            (med_it->second.size() + 1 >= min_run_occur_))
        {
          temp.insert(temp.end(), *ref_it); // new items should go at the end
        }
      }
      OPENMS_LOG_DEBUG << "Removed " << reference_.size() - temp.size() << " of "
                << reference_.size() << " peptides." << endl;
      temp.swap(reference_);
    }
    else // compute overall RT median per sequence (median of medians per run)
    {
      OPENMS_LOG_DEBUG << "Computing overall RT medians per sequence..." << endl;

      // remove peptides that don't occur in enough runs (at least two):
      OPENMS_LOG_DEBUG << "Removing peptides that occur in too few runs..." << endl;
      SeqToList temp;
      for (SeqToList::iterator med_it = medians_per_seq.begin();
           med_it != medians_per_seq.end(); ++med_it)
      {
        if (med_it->second.size() >= min_run_occur_)
        {
          temp.insert(temp.end(), *med_it);
        }
      }
      OPENMS_LOG_DEBUG << "Removed " << medians_per_seq.size() - temp.size() << " of "
                << medians_per_seq.size() << " peptides." << endl;
      temp.swap(medians_per_seq);
      computeMedians_(medians_per_seq, reference_);
    }

    if (verbose && reference_.empty())
    {
      OPENMS_LOG_WARN << "No reference RT information left after filtering!" << endl;
    }

    double max_rt_shift = (double)param_.getValue("max_rt_shift");
    if (max_rt_shift <= 1)
    {
      // compute max. allowed shift from overall retention time range:
      double rt_min = numeric_limits<double>::infinity(), rt_max = -rt_min;
      for (SeqToValue::iterator it = reference_.begin(); it != reference_.end();
           ++it)
      {
        rt_min = min(rt_min, it->second);
        rt_max = max(rt_max, it->second);
      }
      double rt_range = rt_max - rt_min;
      max_rt_shift *= rt_range;
      // in the degenerate case of only one reference point, "max_rt_shift"
      // should be zero (because "rt_range" is zero) - this is covered below
    }
    if (max_rt_shift == 0)
    {
      max_rt_shift = numeric_limits<double>::max();
    }
    OPENMS_LOG_DEBUG << "Max. allowed RT shift (in seconds): " << max_rt_shift << endl;

    // generate RT transformations:
    OPENMS_LOG_DEBUG << "Generating RT transformations..." << endl;
    if (verbose) OPENMS_LOG_INFO << "\nAlignment based on:" << endl; // diagnostic output
    Size offset = 0; // offset in case of internal reference
    for (Int i = 0; i < size + 1; ++i)
    {
      if (i == reference_index_)
      {
        // if one of the input maps was used as reference, it has been skipped
        // so far - now we have to consider it again:
        TransformationDescription trafo;
        trafo.fitModel("identity");
        transforms.push_back(trafo);
        if (verbose)
        {
          OPENMS_LOG_INFO << "- " << reference_.size() << " data points for sample "
                          << i + 1 << " (reference)\n";
        }
        offset = 1;
      }

      if (i >= size) break;

      if (reference_.empty())
      {
        TransformationDescription trafo;
        trafo.fitModel("identity");
        transforms.push_back(trafo);
        continue;
      }
                
      // to be useful for the alignment, a peptide sequence has to occur in the
      // current run ("medians_per_run[i]"), but also in at least one other run
      // ("medians_overall"):
      TransformationDescription::DataPoints data;
      Size n_outliers = 0;
      for (SeqToValue::iterator med_it = medians_per_run[i].begin();
           med_it != medians_per_run[i].end(); ++med_it)
      {
        SeqToValue::const_iterator pos = reference_.find(med_it->first);
        if (pos != reference_.end())
        {
          if (abs(med_it->second - pos->second) <= max_rt_shift)
          { // found, and satisfies "max_rt_shift" condition!
            TransformationDescription::DataPoint point(med_it->second,
                                                       pos->second, pos->first);
            data.push_back(point);
          }
          else
          {
            n_outliers++;
          }
        }
      }
      transforms.emplace_back(data);
      if (verbose)
      {
        OPENMS_LOG_INFO << "- " << data.size() << " data points for sample "
                        << i + offset + 1;
        if (n_outliers) OPENMS_LOG_INFO << " (" << n_outliers << " outliers removed)";
        OPENMS_LOG_INFO << "\n";
      }
    }
    if (verbose) OPENMS_LOG_INFO << endl;

    // delete temporary reference
    if (!reference_given) reference_.clear();
  }

  // explicit template instantiation for Windows DLL:
  template bool OPENMS_DLLAPI MapAlignmentAlgorithmIdentification::getRetentionTimes_<>(const ConsensusMap& features, SeqToList& rt_data);

  // explicit template instantiation for Windows DLL:
  template bool OPENMS_DLLAPI MapAlignmentAlgorithmIdentification::getRetentionTimes_<>(const FeatureMap& features, SeqToList& rt_data);

  const PeptideHit* MapAlignmentAlgorithmIdentification::getBestScoringHit(const std::vector<PeptideHit>& hits, const bool is_higher_score_better)
  {
    auto scoreComparator = PeptideIdentification::getScoreComparator(is_higher_score_better);
    const PeptideHit* best_hit = nullptr;
    for (const auto& hit : hits)
    {
      if (!best_hit || scoreComparator(hit, *best_hit))
      {
        best_hit = &hit;
      }
    }
    return best_hit;
  }

} //namespace
