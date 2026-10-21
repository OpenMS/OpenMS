// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/IDBoostGraph.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <unordered_map>
#include <utility>

namespace OpenMS
{
  using Internal::IDBoostGraph;

  BasicProteinInferenceAlgorithm::BasicProteinInferenceAlgorithm():
      DefaultParamHandler("BasicProteinInferenceAlgorithm"),
      ProgressLogger()
  {
    //TODO allow min_unique_peptides_per_protein (not the same as "use_shared = F" if you want to score the shared ones)
    defaults_.setValue("min_peptides_per_protein", 1,
        "Minimal number of peptides needed for a protein identification."
        " If set to zero, unmatched proteins get a score of -Infinity."
        " If bigger than zero, proteins with less peptides are filtered and evidences removed from the PSMs."
        " PSMs that do not reference any proteins anymore are removed but the spectrum info is kept.");
    defaults_.setMinInt("min_peptides_per_protein", 0);
    defaults_.setValue("score_aggregation_method",
                       "best",
                       "How to aggregate scores of peptides matching to the same protein?");
    defaults_.setValidStrings("score_aggregation_method", {"best","product","sum","maximum"});

    defaults_.setValue("treat_charge_variants_separately", "true",
                       "If this is true, different charge variants of the same peptide sequence count as individual evidences.");
    defaults_.setValidStrings("treat_charge_variants_separately", {"true","false"});

    defaults_.setValue("treat_modification_variants_separately", "true",
                       "If this is true, different modification variants of the same peptide sequence count as individual evidences.");
    defaults_.setValidStrings("treat_modification_variants_separately", {"true","false"});

    defaults_.setValue("use_shared_peptides", "true", "If this is true, shared peptides are used as evidences. Note: shared_peptides are not deleted and potentially resolved in postprocessing as well.");
    defaults_.setValidStrings("use_shared_peptides", {"true","false"});

    defaults_.setValue("skip_count_annotation", "false", "If this is set, peptide counts won't be annotated at the proteins.");
    defaults_.setValidStrings("skip_count_annotation", {"true","false"});

    defaults_.setValue("annotate_indistinguishable_groups", "true", "If this is true, calculates and annotates indistinguishable protein groups.");
    defaults_.setValidStrings("annotate_indistinguishable_groups", {"true","false"});

    defaults_.setValue("greedy_group_resolution", "false", "If this is true, shared peptides will be associated to best proteins only (i.e. become potentially quantifiable razor peptides).");
    defaults_.setValidStrings("greedy_group_resolution", {"true","false"});

    defaults_.setValue("score_type", "", "PSM score type to use for inference. (default: empty = main score)");
    defaults_.setValidStrings("score_type", {"", "PEP", "q-value", "RAW"});

    defaultsToParam_();
  }

  void BasicProteinInferenceAlgorithm::run(PeptideIdentificationList &pep_ids,
                                           ProteinIdentification &prot_id) const
  {
    Size min_peptides_per_protein = static_cast<Size>(param_.getValue("min_peptides_per_protein"));

    SequenceToChargeToPSM best_pep;
    std::unordered_map<std::string, std::pair<ProteinHit*, Size>> acc_to_protein_hitP_and_count;

    std::string requested_score_type_as_string = param_.getValue("score_type").toString();
    auto isr = IDScoreSwitcherAlgorithm::switchToScoreType(pep_ids, requested_score_type_as_string);

    processRun_(
        acc_to_protein_hitP_and_count,
        best_pep,
        prot_id,
        pep_ids
    );

    if (min_peptides_per_protein > 0) //potentially sth was filtered
    {
      std::vector<ProteinIdentification> tmp(1);
      std::swap(tmp[0], prot_id);
      IDFilter::removeDanglingProteinReferences(pep_ids, tmp, true); //TODO allow keeping PSMs without evidence?
      std::swap(tmp[0], prot_id);
    }

    IDScoreSwitcherAlgorithm::switchBackScoreType(pep_ids, isr);
  }

  void BasicProteinInferenceAlgorithm::run(ConsensusMap& cmap, ProteinIdentification& prot_run, bool include_unassigned) const
  {
    runNative_(cmap, include_unassigned);
    if (const auto result = IdentificationDataAdapter::pooledInferenceResult(cmap.getIdentificationData()))
    {
      prot_run = result->proteins;
    }
  }

  void BasicProteinInferenceAlgorithm::run(ConsensusMap& cmap, bool include_unassigned) const
  {
    runNative_(cmap, include_unassigned);
  }

  void BasicProteinInferenceAlgorithm::runNative_(ConsensusMap& cmap, bool include_unassigned) const
  {
    bool group(param_.getValue("annotate_indistinguishable_groups").toBool());
    bool resolve(param_.getValue("greedy_group_resolution").toBool());
    Size min_peptides_per_protein = static_cast<Size>(param_.getValue("min_peptides_per_protein"));
    bool treat_charge_variants_separately(param_.getValue("treat_charge_variants_separately").toBool());
    bool treat_modification_variants_separately(param_.getValue("treat_modification_variants_separately").toBool());
    bool use_shared_peptides(param_.getValue("use_shared_peptides").toBool());

    std::string agg_method_string(param_.getValue("score_aggregation_method").toString());
    AggregationMethod aggregation_method = aggFromString_(agg_method_string);

    // The proteins of the protein run of the identifications (as export writes it) are scored: those of the inference
    // result that pools the runs, else of the run's database sequences. The result replaces that inference result.
    auto& data = cmap.getIdentificationData();
    auto pooled = IdentificationDataAdapter::pooledInferenceResult(data);
    if (!pooled)
    {
      return; // no peptide identifications
    }
    IdentificationData::InferenceResult result = std::move(*pooled);
    ProteinIdentification& prot_run = result.proteins;

    prot_run.setInferenceEngine("TOPPProteinInference");
    prot_run.setInferenceEngineVersion(VersionInfo::getVersion());
    ProteinIdentification::SearchParameters sp = prot_run.getSearchParameters();
    sp.setMetaValue("TOPPProteinInference:aggregation_method", agg_method_string);
    sp.setMetaValue("TOPPProteinInference:use_shared_peptides", use_shared_peptides);
    sp.setMetaValue("TOPPProteinInference:treat_charge_variants_separately", treat_charge_variants_separately);
    sp.setMetaValue("TOPPProteinInference:treat_modification_variants_separately", treat_modification_variants_separately);
    prot_run.setSearchParameters(sp);
    auto& prot_hits = prot_run.getHits();

    IDFilter::keepNBestPeptideHits(cmap, 1); // we should filter for best psm per spec only, since those will be the psms used, also filterUnreferencedProteins depends on it (e.g. after resolution)
    IdentificationDataAdapter::storeInferenceResult(data, result);

    // determine requested score type. This can be a search engine score name or a broader score category (e.g. PEP)
    std::string requested_score_type_as_string = param_.getValue("score_type").toString();

    // The PSMs are scored with the requested score type (e.g., "RAW" = a search engine score, "PEP", "q-value"): the
    // main score, else that of a copy of the map whose main score is switched to it. The result also contains name,
    // score orientation and type before and after switching.
    std::optional<ConsensusMap> switched;
    {
      std::string main_score_name;
      bool main_higher_better = true;
      IDScoreSwitcherAlgorithm::ScoreType main_score_type;
      IDScoreSwitcherAlgorithm().determineScoreNameOrientationAndType(cmap, main_score_name, main_higher_better, main_score_type, include_unassigned);
      if (!requested_score_type_as_string.empty() && IDScoreSwitcherAlgorithm::toScoreTypeEnum(requested_score_type_as_string) != main_score_type)
      {
        switched = cmap;
      }
    }
    const ConsensusMap& scored = switched ? *switched : cmap;
    const auto isr = IDScoreSwitcherAlgorithm::switchToScoreType(switched ? *switched : cmap, requested_score_type_as_string, include_unassigned);

    // The best match of each identification (of the features, then the unassigned ones if requested), the first of the
    // best score as sorting peptide hits makes it
    std::vector<PeptideCandidate> candidates;
    const auto addCandidate = [&](const IdentificationData::QueryMatches& identification)
    {
      const auto& run = *identification.run;
      if (identification.matches.empty() || run.getMoleculeKind() != IdentificationData::MoleculeKind::PEPTIDE) return;
      const auto primary = *run.getPrimaryScore();
      const IdentificationData::Match* best = nullptr;
      double best_score = 0.0;
      for (const auto* match : identification.matches)
      {
        const double score = *run.getScore(match->getId(), primary);
        if (best == nullptr || (isr.requested_score_higher_better ? score > best_score : score < best_score))
        {
          best = match;
          best_score = score;
        }
      }
      PeptideCandidate candidate;
      candidate.score = best_score;
      const AASequence sequence = AASequence::fromString(best->representation);
      candidate.sequence = sequence.toString();
      candidate.unmodified_sequence = sequence.toUnmodifiedString();
      candidate.charge = best->charge;
      if (best->metaValueExists("protein_references")) candidate.protein_references = best->getMetaValue("protein_references").toString();
      for (const auto& evidence : best->sequence_evidence) candidate.accessions.insert(evidence.accession);
      candidates.push_back(std::move(candidate));
    };
    for (const auto& feature : scored)
    {
      for (const auto& identification : feature.getLinkedIdentifications(scored.getIdentificationData()))
      {
        addCandidate(identification);
      }
    }
    if (include_unassigned)
    {
      for (const auto& identification : scored.getUnassignedIdentifications())
      {
        addCandidate(identification);
      }
    }

    // determine inital values for protein scores based on score type
    bool pep_scores = (isr.requested_score_type == IDScoreSwitcherAlgorithm::ScoreType::PEP);
    double initScore = getInitScoreForAggMethod_(aggregation_method, pep_scores || isr.requested_score_higher_better); // if we have pep scores, we will complement to pp during aggregation

    SequenceToChargeToPSM best_pep;
    std::unordered_map<std::string, std::pair<ProteinHit*, Size>> acc_to_protein_hitP_and_count;
    for (auto& prothit : prot_hits)
    {
      prothit.setScore(initScore);
      acc_to_protein_hitP_and_count[prothit.getAccession()] = std::make_pair<ProteinHit*, Size>(&prothit, 0);
    }

    checkCompat_(isr.requested_score_name, aggregation_method);

    aggregatePeptideScores_(best_pep, candidates, isr.requested_score_higher_better);

    updateProteinScores_(
        acc_to_protein_hitP_and_count,
        best_pep,
        pep_scores,
        isr.requested_score_higher_better
    );

    if (pep_scores)
    {
      prot_run.setScoreType("Posterior Probability");
      prot_run.setHigherScoreBetter(true);
    }
    else
    {
      prot_run.setScoreType(isr.requested_score_name); // score name after switch (e.g. "PEP", "q-value")
      prot_run.setHigherScoreBetter(isr.requested_score_higher_better);
    }

    if (min_peptides_per_protein > 0)
    {
      IDFilter::removeMatchingItems<std::vector<ProteinHit>>(prot_run.getHits(),
          IDFilter::HasMaxMetaValue<ProteinHit>("nr_found_peptides", static_cast<int>(min_peptides_per_protein) - 1));
      IdentificationDataAdapter::storeInferenceResult(data, result);
      IDFilter::removeDanglingProteinReferences(cmap, prot_run, true);
    }

    if (group || resolve)
    {
      //TODO you could actually also do the aggregation/inference as well as the resolution on the Graph structure.
      // Groups would be clustered already. Saving some time.
      // But it is quite fast right now already.
      {
        IDBoostGraph ibg{prot_run, cmap, 1, false, include_unassigned, false};
        ibg.computeConnectedComponents();
        // No PSM maps to one of the run's proteins (e.g. none passed an FDR filter): the graph is empty, there is
        // nothing to group or resolve, and IDBoostGraph refuses to work on an empty graph.
        const bool empty_graph = ibg.getNrConnectedComponents() == 0;
        if (resolve && !empty_graph)
        {
          ibg.clusterIndistProteinsAndPeptides(); //TODO check in resolve or do it there if not done yet!
          //Note: the above does not add singleton groups to graph
          ibg.resolveGraphPeptideCentric(true);
          if (group) ibg.annotateIndistProteins(true); // this does not really add singletons since they are not in the graph
          // resolution removed protein references of the PSMs
          ibg.applyToIdentificationData(data, false);
        }
        else if (!resolve && !empty_graph)
        {
          ibg.calculateAndAnnotateIndistProteins(true);
        }
      }
      if (resolve)
      {
        IdentificationDataAdapter::storeInferenceResult(data, result);
        IDFilter::removeUnreferencedProteins(cmap, include_unassigned);
        result = *IdentificationDataAdapter::pooledInferenceResult(data);
        IDFilter::updateProteinGroups(prot_run.getIndistinguishableProteins(), prot_run.getHits());
        if (group) prot_run.fillIndistinguishableGroupsWithSingletons();
      }
      if (group)
      {
        auto & ipg = prot_run.getIndistinguishableProteins();
        std::sort(std::begin(ipg), std::end(ipg));
      }
    }

    prot_run.sort();
    IdentificationDataAdapter::storeInferenceResult(data, result);
  }

  void BasicProteinInferenceAlgorithm::annotateIndistinguishableGroups(ProteinIdentification& proteins,
                                                                       const PeptideIdentificationList& peptides,
                                                                       Size use_top_psms,
                                                                       bool add_singletons)
  {
    // The graph keeps pointers to the hits it is built from and needs non-const ones. Building
    // it and grouping do not modify them; the copy is what lets the interface promise that.
    PeptideIdentificationList peptides_copy{peptides};
    IDBoostGraph ibg{proteins, peptides_copy, use_top_psms, false, false};
    // The graph appends its groups. Start from none, so that they replace the run's groups
    // (also those of an earlier call) instead of adding to them, and restore those on failure.
    auto previous = std::exchange(proteins.getIndistinguishableProteins(), {});
    try
    {
      ibg.calculateAndAnnotateIndistProteins(add_singletons);
    }
    catch (...)
    {
      proteins.getIndistinguishableProteins() = std::move(previous);
      throw;
    }
  }

  void BasicProteinInferenceAlgorithm::run(PeptideIdentificationList &pep_ids,
                                           std::vector<ProteinIdentification> &prot_ids) const
  {
    Size min_peptides_per_protein = static_cast<Size>(param_.getValue("min_peptides_per_protein"));
    IDFilter::keepNBestHits(pep_ids,1); // we should filter for best psm per spec only, since those will be the psms used, also filterUnreferencedProteins depends on it (e.g. after resolution)
    SequenceToChargeToPSM best_pep;
    std::unordered_map<std::string, std::pair<ProteinHit*, Size>> acc_to_protein_hitP_and_count;

    // determine requested score type. This can be a search engine score name or a broader score category (e.g. PEP)
    std::string requested_score_type_as_string = param_.getValue("score_type").toString();
    auto isr = IDScoreSwitcherAlgorithm::switchToScoreType(pep_ids, requested_score_type_as_string);

    for (auto &prot_run : prot_ids)
    {
      processRun_(
          acc_to_protein_hitP_and_count,
          best_pep,
          prot_run,
          pep_ids
          );
    }

    if (min_peptides_per_protein > 0) //potentially sth was filtered
    {
      IDFilter::removeDanglingProteinReferences(pep_ids, prot_ids, true); //TODO allow keeping PSMs without evidence?
    }

    IDScoreSwitcherAlgorithm::switchBackScoreType(pep_ids, isr); // NOP if no switch was performed
  }

  std::vector<BasicProteinInferenceAlgorithm::PeptideCandidate> BasicProteinInferenceAlgorithm::candidates_(
      PeptideIdentificationList& pep_ids,
      const std::string& overall_score_type,
      const std::string& run_id)
  {
    std::vector<PeptideCandidate> candidates;
    for (auto &pep : pep_ids)
    {
      if (pep.getScoreType() != overall_score_type)
      {
        throw OpenMS::Exception::InvalidParameter(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "Differing score_types in the PeptideHits. Aborting...");
      }
      //skip if it does not belong to run
      if (!run_id.empty() && pep.getIdentifier() != run_id)
        continue;
      //skip if no hits (which almost could be considered and error or warning)
      if (pep.getHits().empty())
        continue;
      //make sure that first = best hit
      pep.sort();

      //TODO think about if using any but the best PSM per spectrum makes sense in such a simple aggregation scheme
      const PeptideHit &hit = pep.getHits()[0];
      PeptideCandidate candidate;
      candidate.score = hit.getScore();
      candidate.sequence = hit.getSequence().toString();
      candidate.unmodified_sequence = hit.getSequence().toUnmodifiedString();
      candidate.charge = hit.getCharge();
      if (hit.metaValueExists("protein_references")) candidate.protein_references = hit.getMetaValue("protein_references").toString();
      candidate.accessions = hit.extractProteinAccessionsSet();
      candidates.push_back(std::move(candidate));
    }
    return candidates;
  }

  void BasicProteinInferenceAlgorithm::aggregatePeptideScores_(
      SequenceToChargeToPSM& best_pep,
      const std::vector<PeptideCandidate>& candidates,
      bool higher_better) const
  {
    bool treat_charge_variants_separately(param_.getValue("treat_charge_variants_separately").toBool());
    bool treat_modification_variants_separately(param_.getValue("treat_modification_variants_separately").toBool());
    bool use_shared_peptides(param_.getValue("use_shared_peptides").toBool());

    for (const auto &hit : candidates)
    {
      //skip if shared and option not enabled
      //TODO warn if not present but requested?
      //TODO use nr of evidences to re-calculate sharedness?
      if (!use_shared_peptides &&
      (hit.protein_references.empty() || (hit.protein_references == "non-unique")))
        continue;

      //TODO refactor: this is very similar to IDFilter best per peptide functionality
      const std::string& lookup_seq = treat_modification_variants_separately ? hit.sequence : hit.unmodified_sequence;

      int lookup_charge = 0;
      if (treat_charge_variants_separately)
      {
        lookup_charge = hit.charge;
      }

      if (auto current_best_pep_it = best_pep.find(lookup_seq); current_best_pep_it == best_pep.end())
      { // no entry exist for sequence? initialize seq->charge->&hit
        best_pep[lookup_seq][lookup_charge] = &hit;
      }
      else
      { // a peptide hit for the current sequence exists
        if (auto current_best_pep_charge_it = current_best_pep_it->second.find(lookup_charge); current_best_pep_charge_it == current_best_pep_it->second.end())
        { // no entry for charge? add hit
          current_best_pep_it->second[lookup_charge] = &hit;
        }
        else if (
            (higher_better && (hit.score > current_best_pep_charge_it->second->score)) ||
            (!higher_better && (hit.score < current_best_pep_charge_it->second->score)))
        { // seq with charge already exists? replace if new value has better score
          current_best_pep_charge_it->second = &hit;
        }
      }
    }
  }


  void BasicProteinInferenceAlgorithm::updateProteinScores_(
      std::unordered_map<std::string, std::pair<ProteinHit*, Size>>& acc_to_protein_hitP_and_count,
      const SequenceToChargeToPSM& best_pep,
      bool pep_scores,
      bool higher_better) const
  {
    //TODO Allow count as aggregation method -> i.e. set as protein score?

    if (!higher_better && pep_scores)
    {
      higher_better = true; // We will convert the scores to PPs for multiplication
    }

    bool skip_count_annotation(param_.getValue("skip_count_annotation").toBool());

    std::string agg_method_string(param_.getValue("score_aggregation_method").toString());

    AggregationMethod aggregation_method = aggFromString_(agg_method_string);
    const auto& aggregation_fun = aggFunFromEnum_(aggregation_method, higher_better);

    // update protein scores
    for (const auto &seq_to_map_from_charge_to_pep_hit : best_pep)
    {
      // The next line assumes that PeptideHits of different charge states necessarily share the same
      // protein accessions (might not be the case if we ever allow modified protein databases e.g. in PEF format)
      // TODO this could be done for mods, too (first hashing AASeq, then the mods)
      const std::map<Int, const PeptideCandidate*>& charge_to_peptide_hit = seq_to_map_from_charge_to_pep_hit.second;
      const PeptideCandidate& first_peptide_hit = *charge_to_peptide_hit.begin()->second;
      for (const auto &acc : first_peptide_hit.accessions)
      {
        for (const auto &charge_pep_hit_pair : charge_to_peptide_hit)
        {
          auto prot_count_pair_it = acc_to_protein_hitP_and_count.find(acc);
          if (prot_count_pair_it == acc_to_protein_hitP_and_count.end())
          {
            OPENMS_LOG_WARN << "Warning, skipping pep that maps to a non existent protein accession. "
            << first_peptide_hit.unmodified_sequence << std::endl;
            continue; // very weird, has an accession that was not in the proteins loaded in the beginning
            //TODO error? Suppress log?
          }

          ProteinHit *protein = prot_count_pair_it->second.first;
          prot_count_pair_it->second.second++;

          double new_score = charge_pep_hit_pair.second->score;

          if (pep_scores) // convert PEP to PP
            new_score = 1. - new_score;

          protein->setScore(aggregation_fun(protein->getScore(), new_score));
        }
      }
    }
    if (!skip_count_annotation)
    {
      for (auto& entry : acc_to_protein_hitP_and_count)
      {
        entry.second.first->setMetaValue("nr_found_peptides", entry.second.second);
      }
    }

    //normalize in case of SUM
    if (aggregation_method == AggregationMethod::SUM)
    {
      for (auto& entry : acc_to_protein_hitP_and_count)
      {
        ProteinHit* phitp = entry.second.first;
        phitp->setScore(phitp->getScore() / entry.second.second);
      }
    }
  }

  void BasicProteinInferenceAlgorithm::checkCompat_(
        const std::string& score_name,
        const AggregationMethod& aggregation_method
  ) const
  {
    //TODO do something smart about the scores, e.g. let the user specify a general score type
    // he wants to use and then switch all of them
    if (!IDScoreSwitcherAlgorithm().isScoreType(score_name, IDScoreSwitcherAlgorithm::ScoreType::PEP) &&
        !IDScoreSwitcherAlgorithm().isScoreType(score_name, IDScoreSwitcherAlgorithm::ScoreType::PP) &&
        aggregation_method == AggregationMethod::PROD)
    {
      OPENMS_LOG_WARN << "ProteinInference with multiplicative aggregation "
                         " should probably use Posterior (Error) Probabilities in the Peptide Hits."
                         " Use Percolator with PEP score or run IDPosteriorErrorProbability first.\n";
    }
  }

  void BasicProteinInferenceAlgorithm::checkCompat_(
        const IDScoreSwitcherAlgorithm::ScoreType& score_type,
        const AggregationMethod& aggregation_method
  ) const
  {
    //TODO do something smart about the scores, e.g. let the user specify a general score type
    // he wants to use and then switch all of them
    if ((score_type != IDScoreSwitcherAlgorithm::ScoreType::PEP) &&
        (score_type != IDScoreSwitcherAlgorithm::ScoreType::PP) &&
        aggregation_method == AggregationMethod::PROD)
    {
      OPENMS_LOG_WARN << "ProteinInference with multiplicative aggregation "
                         " should probably use Posterior (Error) Probabilities in the Peptide Hits."
                         " Use Percolator with PEP score or run IDPosteriorErrorProbability first.\n";
    }
  }  

  BasicProteinInferenceAlgorithm::AggregationMethod BasicProteinInferenceAlgorithm::aggFromString_(const std::string& agg_method_string) const
  {
    if (agg_method_string == "best")
    {
      return AggregationMethod::BEST;
    }
    else if (agg_method_string == "product")
    {
      return AggregationMethod::PROD;
    }
    else if (agg_method_string == "sum")
    {
      return AggregationMethod::SUM;
    }
    else if (agg_method_string == "maximum")
    {
      return AggregationMethod::BEST;
    }
    else
    {
      return AggregationMethod::BEST;
    }
  }

  BasicProteinInferenceAlgorithm::fptr BasicProteinInferenceAlgorithm::aggFunFromEnum_(const BasicProteinInferenceAlgorithm::AggregationMethod& agg_method, bool higher_better) const
  {
    switch (agg_method)
    {
      case AggregationMethod::PROD :
        return [](double old_score, double new_score){
          if (new_score > 0.0) //TODO for 0 probability peptides we could also multiply a minimum value
          {
            return old_score * new_score;
          }
          else
          {
            return old_score;
          }
        };
      case AggregationMethod::BEST :
        if (higher_better)
        {
          return [](double old_score, double new_score){return std::fmax(old_score, new_score);};
        }
        else
        {
          return [](double old_score, double new_score){return std::fmin(old_score, new_score);};
        }
      case AggregationMethod::SUM :
        return [](double old_score, double new_score){return old_score + new_score;};
      default:
        throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }
  }

  double BasicProteinInferenceAlgorithm::getInitScoreForAggMethod_(const AggregationMethod& aggregation_method, bool higher_better) const
  {
    switch (aggregation_method)
    {
      //TODO for 0 probability peptides we could also multiply a minimum value
      case AggregationMethod::PROD :
        return 1.0;
      case AggregationMethod::BEST :
        return higher_better ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
      case AggregationMethod::SUM :
        return 0.0;
      default:
        throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
    }
  }

  void BasicProteinInferenceAlgorithm::processRun_(
      std::unordered_map<std::string, std::pair<ProteinHit*, Size>>& acc_to_protein_hitP_and_count,
      SequenceToChargeToPSM& best_pep,
      ProteinIdentification& prot_run,
      PeptideIdentificationList& pep_ids) const
  {
    // TODO actually clearing the scores should be enough, since this algorithm does not change the grouping
    prot_run.getProteinGroups().clear();
    prot_run.getIndistinguishableProteins().clear();

    bool group(param_.getValue("annotate_indistinguishable_groups").toBool());
    bool resolve(param_.getValue("greedy_group_resolution").toBool());
    Size min_peptides_per_protein = static_cast<Size>(param_.getValue("min_peptides_per_protein"));
    bool treat_charge_variants_separately(param_.getValue("treat_charge_variants_separately").toBool());
    bool treat_modification_variants_separately(param_.getValue("treat_modification_variants_separately").toBool());
    bool use_shared_peptides(param_.getValue("use_shared_peptides").toBool());

    std::string agg_method_string(param_.getValue("score_aggregation_method").toString());
    AggregationMethod aggregation_method = aggFromString_(agg_method_string);

    //TODO think about only clearing values or using a big map for all runs together
    acc_to_protein_hitP_and_count.clear();
    best_pep.clear();

    prot_run.setInferenceEngine("TOPPProteinInference");
    prot_run.setInferenceEngineVersion(VersionInfo::getVersion());
    ProteinIdentification::SearchParameters sp = prot_run.getSearchParameters();
    sp.setMetaValue("TOPPProteinInference:aggregation_method", agg_method_string);
    sp.setMetaValue("TOPPProteinInference:use_shared_peptides", use_shared_peptides);
    sp.setMetaValue("TOPPProteinInference:treat_charge_variants_separately", treat_charge_variants_separately);
    sp.setMetaValue("TOPPProteinInference:treat_modification_variants_separately", treat_modification_variants_separately);
    prot_run.setSearchParameters(sp);

    std::string main_score_name;
    bool main_higher_better = true;
    IDScoreSwitcherAlgorithm::ScoreType main_score_type;
    IDScoreSwitcherAlgorithm().determineScoreNameOrientationAndType(pep_ids, main_score_name, main_higher_better, main_score_type);

    bool pep_scores = (main_score_type == IDScoreSwitcherAlgorithm::ScoreType::PEP);
    double initScore = getInitScoreForAggMethod_(aggregation_method, pep_scores || main_higher_better); // if we have pep scores, we will complement to pp during aggregation

    //create Accession to ProteinHit and peptide count map. To have quick access later.
    //If a protein occurs in multiple runs, it picks the last
    for (auto &phit : prot_run.getHits())
    {
      acc_to_protein_hitP_and_count[phit.getAccession()] = std::make_pair<ProteinHit*, Size>(&phit, 0);
      phit.setScore(initScore);
    }

    checkCompat_(main_score_name, aggregation_method);

    const auto candidates = candidates_(pep_ids, main_score_name, prot_run.getIdentifier());
    aggregatePeptideScores_(best_pep, candidates, main_higher_better);

    updateProteinScores_(acc_to_protein_hitP_and_count, best_pep, pep_scores, main_higher_better);

    if (pep_scores) // we converted/ will convert
    {
      prot_run.setScoreType("Posterior Probability");
      prot_run.setHigherScoreBetter(true);
    }
    else
    {
      prot_run.setScoreType(main_score_name);
      prot_run.setHigherScoreBetter(main_higher_better);
    }

    if (min_peptides_per_protein > 0)
    {
      IDFilter::removeMatchingItems<std::vector<ProteinHit>>(prot_run.getHits(),
          IDFilter::HasMaxMetaValue<ProteinHit>("nr_found_peptides", static_cast<int>(min_peptides_per_protein) - 1));
    }

    if (group)
    {
      //TODO you could actually also do the aggregation/inference as well as the resolution on the Graph structure.
      // Groups would be clustered already. Saving some time.
      // But it is quite fast right now already.
      IDBoostGraph ibg{prot_run, pep_ids, 1, false, false};

      ibg.computeConnectedComponents();
      // No PSM of the run maps to one of its proteins (e.g. none passed an FDR filter): the graph is empty, there is
      // nothing to group or resolve, and IDBoostGraph refuses to work on an empty graph.
      const bool empty_graph = ibg.getNrConnectedComponents() == 0;
      if (resolve)
      {
        if (!empty_graph)
        {
          ibg.clusterIndistProteinsAndPeptides(); //TODO check in resolve or do it there if not done yet!
          //Note: the above does not add singleton groups to graph
          ibg.resolveGraphPeptideCentric(true);
          ibg.annotateIndistProteins(true); // this does not really add singletons since they are not in the graph
        }
        IDFilter::removeUnreferencedProteins(prot_run, pep_ids);
        IDFilter::updateProteinGroups(prot_run.getIndistinguishableProteins(), prot_run.getHits());
        prot_run.fillIndistinguishableGroupsWithSingletons();
      }
      else if (!empty_graph)
      {
        ibg.calculateAndAnnotateIndistProteins(true);
      }

      auto & ipg = prot_run.getIndistinguishableProteins();
      std::sort(std::begin(ipg), std::end(ipg));
    }
    else
    {
      if (resolve) // resolution needs groups anyway, so this is very similar to above, except that we remove them in the end.
      {
        IDBoostGraph ibg{prot_run, pep_ids, 1, false, false};

        ibg.computeConnectedComponents();
        if (ibg.getNrConnectedComponents() > 0) // an empty graph has nothing to resolve (see above)
        {
          ibg.clusterIndistProteinsAndPeptides(); //TODO check in resolve or do it there if not done yet!
          //Note: the above does not add singleton groups to graph
          ibg.resolveGraphPeptideCentric(true);
        }
        prot_run.getIndistinguishableProteins().clear();
        IDFilter::removeUnreferencedProteins(prot_run, pep_ids);
      }
    }
  }
} //namespace OpenMS
