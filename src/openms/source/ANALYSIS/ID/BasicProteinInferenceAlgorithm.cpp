// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

#include <algorithm>
#include <iostream>
#include <map>
#include <unordered_map>

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

  void BasicProteinInferenceAlgorithm::run(std::vector<PeptideIdentification> &pep_ids,
                                           ProteinIdentification &prot_id) const
  {
    Size min_peptides_per_protein = static_cast<Size>(param_.getValue("min_peptides_per_protein"));

    std::unordered_map<std::string, std::map<Int, PeptideHit*>> best_pep;
    std::unordered_map<std::string, std::pair<ProteinHit*, Size>> acc_to_protein_hitP_and_count;

    String requested_score_type_as_string = param_.getValue("score_type").toString();
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
      IDFilter::updateProteinReferences(pep_ids, tmp, true); //TODO allow keeping PSMs without evidence?
      std::swap(tmp[0], prot_id);
    }

    IDScoreSwitcherAlgorithm::switchBackScoreType(pep_ids, isr);
  }

  void BasicProteinInferenceAlgorithm::run(ConsensusMap& cmap, ProteinIdentification& prot_run, bool include_unassigned) const
  {
    bool group(param_.getValue("annotate_indistinguishable_groups").toBool());
    bool resolve(param_.getValue("greedy_group_resolution").toBool());
    Size min_peptides_per_protein = static_cast<Size>(param_.getValue("min_peptides_per_protein"));
    bool treat_charge_variants_separately(param_.getValue("treat_charge_variants_separately").toBool());
    bool treat_modification_variants_separately(param_.getValue("treat_modification_variants_separately").toBool());
    bool use_shared_peptides(param_.getValue("use_shared_peptides").toBool());

    std::unordered_map<std::string, std::map<Int, PeptideHit*>> best_pep;
    std::unordered_map<std::string, std::pair<ProteinHit*, Size>> acc_to_protein_hitP_and_count;

    String agg_method_string(param_.getValue("score_aggregation_method").toString());
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
    auto& prot_hits = prot_run.getHits();

    IDFilter::keepNBestPeptideHits(cmap, 1); // we should filter for best psm per spec only, since those will be the psms used, also filterUnreferencedProteins depends on it (e.g. after resolution)

    // determine requested score type. This can be a search engine score name or a broader score category (e.g. PEP)
    String requested_score_type_as_string = param_.getValue("score_type").toString();

    // switch scores to requested score type (e.g., "RAW" = a search engine score, "PEP", "q-value")
    // for convenience, the results als contain name, score orientation and type before and after switching
    // as well as a flag that indicates if a swtich was performed (or e.g., score was already the main score)
    auto isr = IDScoreSwitcherAlgorithm::switchToScoreType(cmap, requested_score_type_as_string, include_unassigned);
    // Here: we can be sure that we have the requested score type that as the main score now    
 
    // determine inital values for protein scores based on score type
    bool pep_scores = (isr.requested_score_type == IDScoreSwitcherAlgorithm::ScoreType::PEP);
    double initScore = getInitScoreForAggMethod_(aggregation_method, pep_scores || isr.requested_score_higher_better); // if we have pep scores, we will complement to pp during aggregation

    for (auto& prothit : prot_hits)
    {
      prothit.setScore(initScore);
      acc_to_protein_hitP_and_count[prothit.getAccession()] = std::make_pair<ProteinHit*, Size>(&prothit, 0);
    }

    checkCompat_(isr.requested_score_name, aggregation_method);

    for (auto& cf : cmap)
    {
      aggregatePeptideScores_(best_pep, 
        cf.getPeptideIdentifications(), 
        isr.requested_score_name, 
        isr.requested_score_higher_better, 
        "");
    }

    if (include_unassigned)
    {
      aggregatePeptideScores_(best_pep, 
        cmap.getUnassignedPeptideIdentifications(),
        isr.requested_score_name,
        isr.requested_score_higher_better,
        "");
    }

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

      IDFilter::updateProteinReferences(cmap, prot_run, true);
    }

    if (group)
    {
      //TODO you could actually also do the aggregation/inference as well as the resolution on the Graph structure.
      // Groups would be clustered already. Saving some time.
      // But it is quite fast right now already.
      IDBoostGraph ibg{prot_run, cmap, 1, false, include_unassigned, false};

      ibg.computeConnectedComponents();
      if (resolve)
      {
        ibg.clusterIndistProteinsAndPeptides(); //TODO check in resolve or do it there if not done yet!
        //Note: the above does not add singleton groups to graph
        ibg.resolveGraphPeptideCentric(true);
        ibg.annotateIndistProteins(true); // this does not really add singletons since they are not in the graph
        IDFilter::removeUnreferencedProteins(cmap, include_unassigned);
        IDFilter::updateProteinGroups(prot_run.getIndistinguishableProteins(), prot_run.getHits());
        prot_run.fillIndistinguishableGroupsWithSingletons();
      }
      else
      {
        ibg.calculateAndAnnotateIndistProteins(true);
      }

      auto & ipg = prot_run.getIndistinguishableProteins();
      std::sort(std::begin(ipg), std::end(ipg));
    }
    else
    {
      if (resolve)
      {
        IDBoostGraph ibg{prot_run, cmap, 1, false, include_unassigned, false};

        ibg.computeConnectedComponents();
        ibg.clusterIndistProteinsAndPeptides(); //TODO check in resolve or do it there if not done yet!
        //Note: the above does not add singleton groups to graph
        ibg.resolveGraphPeptideCentric(true);
        IDFilter::removeUnreferencedProteins(cmap, include_unassigned);
        IDFilter::updateProteinGroups(prot_run.getIndistinguishableProteins(), prot_run.getHits());
      }
    }

    prot_run.sort();
    
    IDScoreSwitcherAlgorithm::switchBackScoreType(cmap, isr, include_unassigned); // NOP if no switch was performed
  }

  void BasicProteinInferenceAlgorithm::run(std::vector<PeptideIdentification> &pep_ids,
                                           std::vector<ProteinIdentification> &prot_ids) const
  {
    Size min_peptides_per_protein = static_cast<Size>(param_.getValue("min_peptides_per_protein"));
    IDFilter::keepNBestHits(pep_ids,1); // we should filter for best psm per spec only, since those will be the psms used, also filterUnreferencedProteins depends on it (e.g. after resolution)
    std::unordered_map<std::string, std::map<Int, PeptideHit*>> best_pep;
    std::unordered_map<std::string, std::pair<ProteinHit*, Size>> acc_to_protein_hitP_and_count;

    // determine requested score type. This can be a search engine score name or a broader score category (e.g. PEP)
    String requested_score_type_as_string = param_.getValue("score_type").toString();
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
      IDFilter::updateProteinReferences(pep_ids, prot_ids, true); //TODO allow keeping PSMs without evidence?
    }

    IDScoreSwitcherAlgorithm::switchBackScoreType(pep_ids, isr); // NOP if no switch was performed
  }

  void BasicProteinInferenceAlgorithm::aggregatePeptideScores_(
      std::unordered_map<std::string, std::map<Int, PeptideHit*>>& best_pep,
      std::vector<PeptideIdentification>& pep_ids,
      const String& overall_score_type,
      bool higher_better,
      const std::string& run_id) const
  {
    bool treat_charge_variants_separately(param_.getValue("treat_charge_variants_separately").toBool());
    bool treat_modification_variants_separately(param_.getValue("treat_modification_variants_separately").toBool());
    bool use_shared_peptides(param_.getValue("use_shared_peptides").toBool());

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
      //for (auto& hit : pep.getHits())
      //{
      PeptideHit &hit = pep.getHits()[0];
      //skip if shared and option not enabled
      //TODO warn if not present but requested?
      //TODO use nr of evidences to re-calculate sharedness?
      if (!use_shared_peptides &&
      (!hit.metaValueExists("protein_references") || (hit.getMetaValue("protein_references") == "non-unique")))
        continue;

      //TODO refactor: this is very similar to IDFilter best per peptide functionality
      String lookup_seq;
      if (!treat_modification_variants_separately)
      {
        lookup_seq = hit.getSequence().toUnmodifiedString();
      }
      else
      {
        lookup_seq = hit.getSequence().toString();
      }

      int lookup_charge = 0;
      if (treat_charge_variants_separately)
      {
        lookup_charge = hit.getCharge();
      }

      auto current_best_pep_it = best_pep.find(lookup_seq);
      if (current_best_pep_it == best_pep.end())
      { // no entry exist for sequence? initialize seq->charge->&hit
        best_pep[lookup_seq][lookup_charge] = &hit;
      }
      else
      { // a peptide hit for the current sequence exists
        auto current_best_pep_charge_it = current_best_pep_it->second.find(lookup_charge);
        if (current_best_pep_charge_it == current_best_pep_it->second.end())
        { // no entry for charge? add hit
          current_best_pep_it->second[lookup_charge] = &hit;
        }
        else if (
            (higher_better && (hit.getScore() > current_best_pep_charge_it->second->getScore())) ||
            (!higher_better && (hit.getScore() < current_best_pep_charge_it->second->getScore())))
        { // seq with charge already exists? replace if new value has better score
          current_best_pep_charge_it->second = &hit;
        }
      }
      //}
    }
  }


  void BasicProteinInferenceAlgorithm::updateProteinScores_(
      std::unordered_map<std::string, std::pair<ProteinHit*, Size>>& acc_to_protein_hitP_and_count,
      const std::unordered_map<std::string, std::map<Int, PeptideHit*>>& best_pep,
      bool pep_scores,
      bool higher_better) const
  {
    //TODO Allow count as aggregation method -> i.e. set as protein score?

    if (!higher_better && pep_scores)
    {
      higher_better = true; // We will convert the scores to PPs for multiplication
    }

    bool skip_count_annotation(param_.getValue("skip_count_annotation").toBool());

    String agg_method_string(param_.getValue("score_aggregation_method").toString());

    AggregationMethod aggregation_method = aggFromString_(agg_method_string);
    const auto& aggregation_fun = aggFunFromEnum_(aggregation_method, higher_better);

    // update protein scores
    for (const auto &seq_to_map_from_charge_to_pep_hit : best_pep)
    {
      // The next line assumes that PeptideHits of different charge states necessarily share the same
      // protein accessions (might not be the case if we ever allow modified protein databases e.g. in PEF format)
      // TODO this could be done for mods, too (first hashing AASeq, then the mods)
      const std::map<Int, PeptideHit*>& charge_to_peptide_hit = seq_to_map_from_charge_to_pep_hit.second;
      const PeptideHit& first_peptide_hit = *charge_to_peptide_hit.begin()->second;
      for (const auto &acc : first_peptide_hit.extractProteinAccessionsSet())
      {
        for (const auto &charge_pep_hit_pair : charge_to_peptide_hit)
        {
          auto prot_count_pair_it = acc_to_protein_hitP_and_count.find(std::string(acc));
          if (prot_count_pair_it == acc_to_protein_hitP_and_count.end())
          {
            OPENMS_LOG_WARN << "Warning, skipping pep that maps to a non existent protein accession. "
            << first_peptide_hit.getSequence().toUnmodifiedString() << std::endl;
            continue; // very weird, has an accession that was not in the proteins loaded in the beginning
            //TODO error? Suppress log?
          }

          ProteinHit *protein = prot_count_pair_it->second.first;
          prot_count_pair_it->second.second++;

          const PeptideHit& pep_hit = *charge_pep_hit_pair.second;
          double new_score = pep_hit.getScore();

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
        const String& score_name,
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
      std::unordered_map<std::string, std::map<Int, PeptideHit*>>& best_pep,
      ProteinIdentification& prot_run,
      std::vector<PeptideIdentification>& pep_ids) const
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

    String agg_method_string(param_.getValue("score_aggregation_method").toString());
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

    String main_score_name;
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

    aggregatePeptideScores_(best_pep, pep_ids, main_score_name, main_higher_better, prot_run.getIdentifier());

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
      if (resolve)
      {
        ibg.clusterIndistProteinsAndPeptides(); //TODO check in resolve or do it there if not done yet!
        //Note: the above does not add singleton groups to graph
        ibg.resolveGraphPeptideCentric(true);
        ibg.annotateIndistProteins(true); // this does not really add singletons since they are not in the graph
        IDFilter::removeUnreferencedProteins(prot_run, pep_ids);
        IDFilter::updateProteinGroups(prot_run.getIndistinguishableProteins(), prot_run.getHits());
        prot_run.fillIndistinguishableGroupsWithSingletons();
      }
      else
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
        ibg.clusterIndistProteinsAndPeptides(); //TODO check in resolve or do it there if not done yet!
        //Note: the above does not add singleton groups to graph
        ibg.resolveGraphPeptideCentric(true);
        prot_run.getIndistinguishableProteins().clear();
        IDFilter::removeUnreferencedProteins(prot_run, pep_ids);
      }
    }
  }
} //namespace OpenMS
