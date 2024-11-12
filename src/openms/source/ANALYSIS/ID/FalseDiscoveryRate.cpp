// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <boost/foreach.hpp> // must be first, otherwise Q_FOREACH macro will wreak havoc
#include <boost/regex.hpp>

#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDScoreGetterSetter.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <algorithm>

// #define FALSE_DISCOVERY_RATE_DEBUG
// #undef  FALSE_DISCOVERY_RATE_DEBUG

using namespace std;

namespace OpenMS
{
  FalseDiscoveryRate::FalseDiscoveryRate() :
    DefaultParamHandler("FalseDiscoveryRate")
  {
    defaults_.setValue("no_qvalues", "false", "If 'true' strict FDRs will be calculated instead of q-values (the default)");
    defaults_.setValidStrings("no_qvalues", {"true","false"});
    defaults_.setValue("use_all_hits", "false", "If 'true' not only the first hit, but all are used (peptides only)");
    defaults_.setValidStrings("use_all_hits", {"true","false"});
    defaults_.setValue("split_charge_variants", "false", "If 'true' charge variants are treated separately (for peptides of combined target/decoy searches only).");
    defaults_.setValidStrings("split_charge_variants", {"true","false"});
    defaults_.setValue("treat_runs_separately", "false", "If 'true' different search runs are treated separately (for peptides of combined target/decoy searches only).");
    defaults_.setValidStrings("treat_runs_separately", {"true","false"});
    defaults_.setValue("add_decoy_peptides", "false", "If 'true' decoy peptides will be written to output file, too. The q-value is set to the closest target score.");
    defaults_.setValidStrings("add_decoy_peptides", {"true","false"});
    defaults_.setValue("add_decoy_proteins", "false", "If 'true' decoy proteins will be written to output file, too. The q-value is set to the closest target score.");
    defaults_.setValidStrings("add_decoy_proteins", {"true","false"});
    defaults_.setValue("conservative", "true", "If 'true' (D+1)/T instead of (D+1)/(T+D) is used as a formula.");
    defaults_.setValidStrings("conservative", {"true","false"});
    //defaults_.setValue("subprotein_level", "PSM", "Choose PSM or peptide or PSM+peptide");
    //defaults_.setValidStrings("subprotein_level", {"PSM","peptide","PSM+peptide"});
    defaultsToParam_();
  }

  bool isFirstBetterScore(double first, double second, bool isHigherBetter)
  {
    if (isHigherBetter) return first > second; else return first < second;
  }

  void FalseDiscoveryRate::apply(vector<PeptideIdentification>& ids, bool annotate_peptide_fdr) const
  {
    bool q_value = !param_.getValue("no_qvalues").toBool();
    bool use_all_hits = param_.getValue("use_all_hits").toBool();
    bool treat_runs_separately = param_.getValue("treat_runs_separately").toBool();
    bool split_charge_variants = param_.getValue("split_charge_variants").toBool();
    bool add_decoy_peptides = param_.getValue("add_decoy_peptides").toBool();

    if (ids.empty())
    {
      OPENMS_LOG_WARN << "No peptide identifications given to FalseDiscoveryRate! No calculation performed.\n";
      return;
    }

    bool higher_score_better = ids.begin()->isHigherScoreBetter();

    // first search for all identifiers and charge variants
    set<String> identifiers;
    set<SignedSize> charge_variants;
    for (auto it = ids.begin(); it != ids.end(); ++it)
    {
      identifiers.insert(it->getIdentifier());
      it->sort();

      if (!use_all_hits && it->getHits().size() > 1)
      {
        it->getHits().resize(1);
      }

      for (auto pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        charge_variants.insert(pit->getCharge());
      }
    }

#ifdef FALSE_DISCOVERY_RATE_DEBUG
    cerr << "#id-runs: " << identifiers.size() << " ";
    for (auto it = identifiers.begin(); it != identifiers.end(); ++it)
    {
      cerr << "," << *it;
    }
    cerr << endl;


    cerr << "#of charge states: " << charge_variants.size() << " ";
    for (auto it = charge_variants.begin(); it != charge_variants.end(); ++it)
    {
      cerr << "," << *it;
    }
    cerr << endl;
#endif

    for (auto zit = charge_variants.begin(); zit != charge_variants.end(); ++zit)
    {
#ifdef FALSE_DISCOVERY_RATE_DEBUG
      cerr << "Charge variant=" << *zit << endl;
#endif

      // for all identifiers
      for (auto iit = identifiers.begin(); iit != identifiers.end(); ++iit)
      {
        if (!treat_runs_separately && iit != identifiers.begin())
        {
          continue; //only take the first run
        }

#ifdef FALSE_DISCOVERY_RATE_DEBUG
        cerr << "Id-run: " << *iit << endl;
#endif
        // get the scores of all peptide hits
        vector<double> target_scores, decoy_scores;
        map<String, double> peptide_to_best_decoy_score, peptide_to_best_target_score;
        for (auto it = ids.begin(); it != ids.end(); ++it)
        {
          // if runs should be treated separately, the identifiers must be the same
          if (treat_runs_separately && it->getIdentifier() != *iit)
          {
            continue;
          }

          for (Size i = 0; i < it->getHits().size(); ++i)
          {
            if (split_charge_variants && it->getHits()[i].getCharge() != *zit)
            {
              continue;
            }

            if (!it->getHits()[i].metaValueExists("target_decoy"))
            {
              OPENMS_LOG_FATAL_ERROR << "Meta value 'target_decoy' does not exists, reindex the idXML file with 'PeptideIndexer' first (run-id='" << it->getIdentifier() << ", rank=" << i + 1 << " of " << it->getHits().size() << ")!" << endl;
              throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Meta value 'target_decoy' does not exist!");
            }

            String target_decoy(it->getHits()[i].getMetaValue("target_decoy"));
            const String peptide_sequence = it->getHits()[i].getSequence().toUnmodifiedString();
            const double score = it->getHits()[i].getScore();

            if (target_decoy == "target" || target_decoy == "target+decoy")
            {
              target_scores.push_back(score);

              if (annotate_peptide_fdr)
              {
                // store best score for peptide (unmodified sequence)              
                auto [entry_it, success] = peptide_to_best_target_score.emplace(peptide_sequence, score); // try to construct in place (performance)

                if (!success && // emplace failed because key was already present -> replace if current score is better?
                    isFirstBetterScore(score, entry_it->second, higher_score_better)) 
                {
                  entry_it->second = score;
                }                           
              }
            }
            else
            {
              if (target_decoy == "decoy")
              {
                decoy_scores.push_back(score);

                if (annotate_peptide_fdr)
                {
                  auto [entry_it, success] = peptide_to_best_decoy_score.emplace(peptide_sequence, score); // try to construct in place (performance)

                  if (!success && // emplace failed because key was already present -> replace if current score is better?
                      isFirstBetterScore(score, entry_it->second, higher_score_better)) 
                  {
                    entry_it->second = score;
                  }
                }
              }
              else
              {
                if (!target_decoy.empty())
                {
                  throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown value of meta value 'target_decoy'", target_decoy);
                }
              }
            }
          }
        }

#ifdef FALSE_DISCOVERY_RATE_DEBUG
        cerr << "#target-scores=" << target_scores.size() << ", #decoy-scores=" << decoy_scores.size() << endl;
#endif

        // check decoy scores
        if (decoy_scores.empty())
        {
          String error_string = "FalseDiscoveryRate: #decoy sequences is zero! Setting all target sequences to q-value/FDR 0! ";
          if (split_charge_variants || treat_runs_separately)
          {
            error_string += "(";
            if (split_charge_variants)
            {
              error_string += "charge_variant=" + String(*zit) + " ";
            }
            if (treat_runs_separately)
            {
              error_string += "run-id=" + *iit;
            }
            error_string += ")";
          }
          OPENMS_LOG_ERROR << error_string << std::endl;
        }

        // check target scores
        if (target_scores.empty())
        {
          String error_string = "FalseDiscoveryRate: #target sequences is zero! Ignoring. ";
          if (split_charge_variants || treat_runs_separately)
          {
            error_string += "(";
            if (split_charge_variants)
            {
              error_string += "charge_variant=" + String(*zit) + " ";
            }
            if (treat_runs_separately)
            {
              error_string += "run-id=" + *iit;
            }
            error_string += ")";
          }
          OPENMS_LOG_ERROR << error_string << std::endl;
        }

        if (target_scores.empty() || decoy_scores.empty())
        {
          // now remove the relevant entries, or put 'pseudo-scores' in
          for (auto it = ids.begin(); it != ids.end(); ++it)
          {
            // if runs should be treated separately, the identifiers must be the same
            if (treat_runs_separately && it->getIdentifier() != *iit)
            {
              continue;
            }

            vector<PeptideHit> hits(it->getHits()), new_hits;
            for (Size i = 0; i < hits.size(); ++i)
            {
              if (split_charge_variants && hits[i].getCharge() != *zit)
              {
                new_hits.push_back(hits[i]);
                continue;
              }

              if (!hits[i].metaValueExists("target_decoy"))
              {
                OPENMS_LOG_FATAL_ERROR << "Meta value 'target_decoy' does not exists, reindex the idXML file with 'PeptideIndexer' (run-id='" << it->getIdentifier() << ", rank=" << i + 1 << " of " << hits.size() << ")!" << endl;
                throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Meta value 'target_decoy' does not exist!");
              }

              String target_decoy(hits[i].getMetaValue("target_decoy"));
              if (target_decoy == "target" || target_decoy == "target+decoy")
              {
                // if it is a target hit, there are no decoys, fdr/q-value should be zero then
                new_hits.push_back(hits[i]);
                String score_type = it->getScoreType() + "_score";
                new_hits.back().setMetaValue(score_type, new_hits.back().getScore());
                new_hits.back().setScore(0);
              }
              else
              {
                if (target_decoy != "decoy")
                {
                  throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown value of meta value 'target_decoy'", target_decoy);
                }
              }
            }
            it->setHits(new_hits);
          }
          continue;
        }

        // calculate fdr for the forward scores
        map<double, double> score_to_fdr;
        calculateFDRs_(score_to_fdr, target_scores, decoy_scores, q_value, higher_score_better);

        // calculate peptide FDR
        if (annotate_peptide_fdr)
        {
          vector<double> decoy_peptide_scores, target_peptide_scores;
          for (const auto& ps : peptide_to_best_decoy_score)
          {
            decoy_peptide_scores.push_back(ps.second);
          }
          for (const auto& ps : peptide_to_best_target_score)
          {
            target_peptide_scores.push_back(ps.second);
          }      
          map<double, double> score_to_peptide_fdr;
          calculateFDRs_(score_to_peptide_fdr, target_peptide_scores, decoy_peptide_scores, q_value, higher_score_better);
          // overwrite best peptide score with peptide q-value
          for (auto& ps : peptide_to_best_decoy_score)
          {
            ps.second = score_to_peptide_fdr[ps.second];
          }
          for (auto& ps : peptide_to_best_target_score)
          {
            ps.second = score_to_peptide_fdr[ps.second];
          }
        }

        // annotate fdr
        for (auto it = ids.begin(); it != ids.end(); ++it)
        {
          // if runs should be treated separately, the identifiers must be the same
          if (treat_runs_separately && it->getIdentifier() != *iit)
          {
            continue;
          }

          String score_type = it->getScoreType() + "_score";
          vector<PeptideHit> hits;
          for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
          {
            PeptideHit hit = *pit;

            if (split_charge_variants && pit->getCharge() != *zit)
            {
              hits.push_back(*pit);
              continue;
            }
            if (hit.metaValueExists("target_decoy"))
            {
              String meta_value = (String)hit.getMetaValue("target_decoy");
              if (meta_value == "decoy" && !add_decoy_peptides)
              {
                continue;
              }

              if (annotate_peptide_fdr)
              {
                const String peptide_sequence = hit.getSequence().toUnmodifiedString();
                double peptide_fdr;
                if (meta_value == "decoy")
                {
                  peptide_fdr = peptide_to_best_decoy_score[peptide_sequence];
                }
                else
                {
                  peptide_fdr = peptide_to_best_target_score[peptide_sequence];
                }
                if (q_value)
                {
                  hit.setMetaValue("peptide q-value", peptide_fdr);
                }
                else
                {
                  hit.setMetaValue("peptide FDR", peptide_fdr);
                }
              }              
            }
            hit.setMetaValue(score_type, pit->getScore());
            hit.setScore(score_to_fdr[pit->getScore()]);
            hits.push_back(hit);
          }
          it->getHits().swap(hits);
        }
      }
      if (!split_charge_variants)
      {
        break;
      }
    }

    // higher-score-better can be set now, calculations are finished
    for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
    {
      if (q_value)
      {
        if (it->getScoreType() != "q-value")
        {
          it->setScoreType("q-value");
        }
      }
      else
      {
        if (it->getScoreType() != "FDR")
        {
          it->setScoreType("FDR");
        }
      }
      it->setHigherScoreBetter(false);
      it->assignRanks();
    }

    return;
  }

  void FalseDiscoveryRate::apply(vector<PeptideIdentification>& fwd_ids, vector<PeptideIdentification>& rev_ids) const
  {
    if (fwd_ids.empty() || rev_ids.empty())
    {
      return;
    }
    vector<double> target_scores, decoy_scores;
    // get the scores of all peptide hits
    for (vector<PeptideIdentification>::const_iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
    {
      for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        target_scores.push_back(pit->getScore());
      }
    }

    for (vector<PeptideIdentification>::const_iterator it = rev_ids.begin(); it != rev_ids.end(); ++it)
    {
      for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        decoy_scores.push_back(pit->getScore());
      }
    }

    bool q_value = !param_.getValue("no_qvalues").toBool();
    bool higher_score_better = fwd_ids.begin()->isHigherScoreBetter();
    bool add_decoy_peptides = param_.getValue("add_decoy_peptides").toBool();
    // calculate fdr for the forward scores
    map<double, double> score_to_fdr;
    calculateFDRs_(score_to_fdr, target_scores, decoy_scores, q_value, higher_score_better);

    // annotate fdr
    String score_type = fwd_ids.begin()->getScoreType() + "_score";
    for (vector<PeptideIdentification>::iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
    {
      if (q_value)
      {
        it->setScoreType("q-value");
      }
      else
      {
        it->setScoreType("FDR");
      }

      it->setHigherScoreBetter(false);
      vector<PeptideHit> hits = it->getHits();
      for (vector<PeptideHit>::iterator pit = hits.begin(); pit != hits.end(); ++pit)
      {
#ifdef FALSE_DISCOVERY_RATE_DEBUG
        cerr << pit->getScore() << " " << score_to_fdr[pit->getScore()] << endl;
#endif
        pit->setMetaValue(score_type, pit->getScore());
        pit->setScore(score_to_fdr[pit->getScore()]);
      }
      it->setHits(hits);
    }
    //write as well decoy peptides
    if (add_decoy_peptides)
    {
      score_type = rev_ids.begin()->getScoreType() + "_score";
      for (vector<PeptideIdentification>::iterator it = rev_ids.begin(); it != rev_ids.end(); ++it)
      {
        if (q_value)
        {
          it->setScoreType("q-value");
        }
        else
        {
          it->setScoreType("FDR");
        }

        it->setHigherScoreBetter(false);
        vector<PeptideHit> hits = it->getHits();
        for (vector<PeptideHit>::iterator pit = hits.begin(); pit != hits.end(); ++pit)
        {
#ifdef FALSE_DISCOVERY_RATE_DEBUG
          cerr << pit->getScore() << " " << score_to_fdr[pit->getScore()] << endl;
#endif
          pit->setMetaValue(score_type, pit->getScore());
          pit->setScore(score_to_fdr[pit->getScore()]);
        }
        it->setHits(hits);
      }
    }

    return;
  }

  void FalseDiscoveryRate::apply(vector<ProteinIdentification>& ids) const
  {
    bool q_value = !param_.getValue("no_qvalues").toBool();
    bool higher_score_better = ids.begin()->isHigherScoreBetter();
    bool add_decoy_proteins = param_.getValue("add_decoy_proteins").toBool();

    if (ids.empty())
    {
      OPENMS_LOG_WARN << "No protein identifications given to FalseDiscoveryRate! No calculation performed.\n";
      return;
    }

    vector<double> target_scores, decoy_scores;
    for (auto it = ids.begin(); it != ids.end(); ++it)
    {
      for (auto pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        if (!pit->metaValueExists("target_decoy"))
        {
          OPENMS_LOG_FATAL_ERROR << "Meta value 'target_decoy' does not exists, reindex the idXML file with 'PeptideIndexer' (run-id='" << it->getIdentifier() << ", accession=" << pit->getAccession() << ")!" << endl;
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Meta value 'target_decoy' does not exist!");
        }

        String target_decoy = pit->getMetaValue("target_decoy");
        if (target_decoy == "decoy")
        {
          decoy_scores.push_back(pit->getScore());
        }
        else if (target_decoy == "target")
        {
          target_scores.push_back(pit->getScore());
        }
        else
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown value of meta value 'target_decoy'", target_decoy);
        }
      }
    }


    // calculate fdr for the forward scores
    map<double, double> score_to_fdr;
    calculateFDRs_(score_to_fdr, target_scores, decoy_scores, q_value, higher_score_better);

    // annotate fdr
    String score_type = ids.begin()->getScoreType() + "_score";
    for (auto it = ids.begin(); it != ids.end(); ++it)
    {
      if (q_value)
      {
        it->setScoreType("q-value");
      }
      else
      {
        it->setScoreType("FDR");
      }
      it->setHigherScoreBetter(false);
      const vector<ProteinHit>& old_hits = it->getHits();
      vector<ProteinHit> new_hits;
      for (auto hit : old_hits) // NOTE: performs copy
      {
        // Add decoy proteins only if add_decoy_proteins is set
        if (add_decoy_proteins || hit.getMetaValue("target_decoy") != "decoy")
        {
          hit.setMetaValue(score_type, hit.getScore());
          hit.setScore(score_to_fdr[hit.getScore()]);
          new_hits.push_back(std::move(hit));
        }
      }
      it->setHits(std::move(new_hits));
    }
  }

  void FalseDiscoveryRate::apply(vector<ProteinIdentification>& fwd_ids, vector<ProteinIdentification>& rev_ids) const
  {
    if (fwd_ids.empty() || rev_ids.empty())
    {
      return;
    }
    vector<double> target_scores, decoy_scores;
    // get the scores of all peptide hits
    for (vector<ProteinIdentification>::const_iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
    {
      for (vector<ProteinHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        target_scores.push_back(pit->getScore());
      }
    }
    for (vector<ProteinIdentification>::const_iterator it = rev_ids.begin(); it != rev_ids.end(); ++it)
    {
      for (vector<ProteinHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        decoy_scores.push_back(pit->getScore());
      }
    }

    bool q_value = !param_.getValue("no_qvalues").toBool();
    bool higher_score_better = fwd_ids.begin()->isHigherScoreBetter();
    // calculate fdr for the forward scores
    map<double, double> score_to_fdr;
    calculateFDRs_(score_to_fdr, target_scores, decoy_scores, q_value, higher_score_better);

    // annotate fdr
    String score_type = fwd_ids.begin()->getScoreType() + "_score";
    for (vector<ProteinIdentification>::iterator it = fwd_ids.begin(); it != fwd_ids.end(); ++it)
    {
      if (q_value)
      {
        it->setScoreType("q-value");
      }
      else
      {
        it->setScoreType("FDR");
      }
      it->setHigherScoreBetter(false);
      vector<ProteinHit> hits = it->getHits();
      for (vector<ProteinHit>::iterator pit = hits.begin(); pit != hits.end(); ++pit)
      {
        pit->setMetaValue(score_type, pit->getScore());
        pit->setScore(score_to_fdr[pit->getScore()]);
      }
      it->setHits(hits);
    }
  }

  IdentificationData::ScoreTypeRef FalseDiscoveryRate::applyToObservationMatches(
      IdentificationData& id_data, IdentificationData::ScoreTypeRef score_ref)
  const
  {
    bool use_all_hits = param_.getValue("use_all_hits").toBool();
    bool include_decoys = param_.getValue("add_decoy_peptides").toBool();
    vector<double> target_scores, decoy_scores;
    map<IdentificationData::IdentifiedMolecule, bool> molecule_to_decoy;
    map<IdentificationData::ObservationMatchRef, double> match_to_score;
    if (use_all_hits)
    {
      for (auto it = id_data.getObservationMatches().begin();
           it != id_data.getObservationMatches().end(); ++it)
      {
        handleObservationMatch_(it, score_ref, target_scores, decoy_scores,
                          molecule_to_decoy, match_to_score);
      }
    }
    else
    {
      vector<IdentificationData::ObservationMatchRef> best_matches =
          id_data.getBestMatchPerObservation(score_ref);
      for (const auto& match_ref : best_matches)
      {
        handleObservationMatch_(match_ref, score_ref, target_scores, decoy_scores,
                          molecule_to_decoy, match_to_score);
      }
    }

    map<double, double> score_to_fdr;
    bool higher_better = score_ref->higher_better;
    bool use_qvalue = !param_.getValue("no_qvalues").toBool();
    calculateFDRs_(score_to_fdr, target_scores, decoy_scores, use_qvalue,
                   higher_better);

    IdentificationData::ScoreType fdr_score;
    fdr_score.higher_better = false;
    if (use_qvalue)
    {
      fdr_score.cv_term = CVTerm("MS:1002354", "PSM-level q-value", "MS");
    }
    else
    {
      fdr_score.cv_term = CVTerm("MS:1002355", "PSM-level FDRScore", "MS");
    }
    IdentificationData::ScoreTypeRef fdr_ref =
        id_data.registerScoreType(fdr_score);
    for (IdentificationData::ObservationMatches::iterator it =
           id_data.getObservationMatches().begin(); it !=
           id_data.getObservationMatches().end(); ++it)
    {
      if (!include_decoys)
      {
        auto pos = molecule_to_decoy.find(it->identified_molecule_var);
        if ((pos != molecule_to_decoy.end()) && pos->second) continue;
      }
      auto pos = match_to_score.find(it);
      if (pos == match_to_score.end()) continue;
      double fdr = score_to_fdr.at(pos->second);
      id_data.addScore(it, fdr_ref, fdr);
    }
    return fdr_ref;
  }


  void FalseDiscoveryRate::handleObservationMatch_(
    IdentificationData::ObservationMatchRef match_ref,
    IdentificationData::ScoreTypeRef score_ref,
    vector<double>& target_scores, vector<double>& decoy_scores,
    map<IdentificationData::IdentifiedMolecule, bool>& molecule_to_decoy,
    map<IdentificationData::ObservationMatchRef, double>& match_to_score) const
  {
    const IdentificationData::IdentifiedMolecule& molecule_var =
      match_ref->identified_molecule_var;
    IdentificationData::MoleculeType molecule_type =
      molecule_var.getMoleculeType();
    if (molecule_type == IdentificationData::MoleculeType::COMPOUND)
    {
      return; // compounds don't have parents with target/decoy status
    }
    pair<double, bool> score = match_ref->getScore(score_ref);
    if (!score.second) return; // no score of this type
    match_to_score[match_ref] = score.first;
    auto pos = molecule_to_decoy.find(molecule_var);
    bool is_decoy;
    if (pos == molecule_to_decoy.end()) // new molecule
    {
      if (molecule_type == IdentificationData::MoleculeType::PROTEIN)
      {
        is_decoy = molecule_var.getIdentifiedPeptideRef()->allParentsAreDecoys();
      }
      else // if (molecule_type == IdentificationData::MoleculeType::RNA)
      {
        is_decoy = molecule_var.getIdentifiedOligoRef()->allParentsAreDecoys();
      }
      molecule_to_decoy[molecule_var] = is_decoy;
    }
    else
    {
      is_decoy = pos->second;
    }
    if (is_decoy)
    {
      decoy_scores.push_back(score.first);
    }
    else
    {
      target_scores.push_back(score.first);
    }
  }


  void FalseDiscoveryRate::calculateFDRs_(map<double, double>& score_to_fdr, vector<double>& target_scores, vector<double>& decoy_scores, bool q_value, bool higher_score_better) const
  {
    Size number_of_target_scores = target_scores.size();
    // sort the scores
    if (higher_score_better && !q_value)
    {
      sort(target_scores.rbegin(), target_scores.rend());
      sort(decoy_scores.rbegin(), decoy_scores.rend());
    }
    else if (!higher_score_better && !q_value)
    {
      sort(target_scores.begin(), target_scores.end());
      sort(decoy_scores.begin(), decoy_scores.end());
    }
    else if (higher_score_better)
    {
      sort(target_scores.begin(), target_scores.end());
      sort(decoy_scores.rbegin(), decoy_scores.rend());
    }
    else
    {
      sort(target_scores.rbegin(), target_scores.rend());
      sort(decoy_scores.begin(), decoy_scores.end());
    }

    Size j = 0;

    if (q_value)
    {
      double minimal_fdr = 1.;
      for (Size i = 0; i != target_scores.size(); ++i)
      {
        if (decoy_scores.empty())
        {
          // set FDR to 0 (done below automatically)
        }
        else if (i == 0 && j == 0)
        {
          while (j != decoy_scores.size()
                && ((target_scores[i] <= decoy_scores[j] && higher_score_better) ||
                    (target_scores[i] >= decoy_scores[j] && !higher_score_better)))
          {
            ++j;
          }
        }
        else
        {
          if (j == decoy_scores.size())
          {
            j--;
          }
          while (j != 0
                && ((target_scores[i] > decoy_scores[j] && higher_score_better) ||
                    (target_scores[i] < decoy_scores[j] && !higher_score_better)))
          {
            --j;
          }
          // Since j has to be equal to the number of fps above the threshold we add one
          if ((target_scores[i] <= decoy_scores[j] && higher_score_better)
             || (target_scores[i] >= decoy_scores[j] && !higher_score_better))
          {
            ++j;
          }
        }

#ifdef FALSE_DISCOVERY_RATE_DEBUG
        cerr << target_scores[i] << " " << decoy_scores[j] << " " << i << " " << j << " ";
#endif

        double fdr = 0.;

        if (minimal_fdr >= (double)j / (number_of_target_scores - i))
        {
          minimal_fdr = (double)j / (number_of_target_scores - i);
        }
        fdr = minimal_fdr;

#ifdef FALSE_DISCOVERY_RATE_DEBUG
        cerr << fdr << endl;
#endif
        score_to_fdr[target_scores[i]] = fdr;

      }
    }
    else
    {
      for (Size i = 0; i != target_scores.size(); ++i)
      {
        while (j != decoy_scores.size() &&
               ((target_scores[i] <= decoy_scores[j] && higher_score_better) ||
                (target_scores[i] >= decoy_scores[j] && !higher_score_better)))
        {
          ++j;
        }

#ifdef FALSE_DISCOVERY_RATE_DEBUG
        cerr << target_scores[i] << " " << decoy_scores[j] << " " << i << " " << j << " ";
#endif
        double fdr(0);

        fdr = (double)j / (double)(i + 1);

#ifdef FALSE_DISCOVERY_RATE_DEBUG
        cerr << fdr << endl;
#endif
        score_to_fdr[target_scores[i]] = fdr;
      }
    }

    // assign q-value of decoy_score to closest target_score
    for (Size i = 0; i != decoy_scores.size(); ++i)
    {
      const double& ds = decoy_scores[i];

      // advance target index until score is better than decoy score
      size_t k{0};
      while (k != target_scores.size() &&
             ((target_scores[k] <= ds && higher_score_better) ||
              (target_scores[k] >= ds && !higher_score_better)))
      {
        ++k;
      }

      // corner cases
      if (k == 0)
      {
        if (!target_scores.empty())
        {
          score_to_fdr[ds] = score_to_fdr[target_scores[0]];
          continue;
        }
        else
        {
          score_to_fdr[ds] = 1.0;
          continue;
        }
      }

      if (k == target_scores.size()) { score_to_fdr[ds] = score_to_fdr[target_scores.back()]; continue; }

      if (fabs(target_scores[k] - ds) < fabs(target_scores[k - 1] - ds))
      {
        score_to_fdr[ds] = score_to_fdr[target_scores[k]];
      }
      else
      {
        score_to_fdr[ds] = score_to_fdr[target_scores[k - 1]];
      }
    }
  }

  //TODO does not support "by run" and/or "by charge"
  //TODO could be done for a percentage of FalsePos instead of a number
  //TODO can be templated for proteins
  double FalseDiscoveryRate::rocN(const vector<PeptideIdentification>& ids, Size fp_cutoff) const
  {
    bool higher_score_better(ids.begin()->isHigherScoreBetter());
    bool use_all_hits = param_.getValue("use_all_hits").toBool();

    ScoreToTgtDecLabelPairs scores_labels;
    IDScoreGetterSetter::getScores_(scores_labels, ids, use_all_hits);

    if (scores_labels.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No scores could be extracted!");
    }

    if (higher_score_better)
    { // decreasing
      std::sort(scores_labels.rbegin(), scores_labels.rend());
    }
    else
    { // increasing
      std::sort(scores_labels.begin(), scores_labels.end());
    }
    // if fp_cutoff is zero do the full AUC.
    return rocN(scores_labels, fp_cutoff == 0 ? scores_labels.size() : fp_cutoff);
  }

  double FalseDiscoveryRate::rocN(const vector<PeptideIdentification>& ids, Size fp_cutoff, const String& identifier) const
  {
    bool higher_score_better(ids.begin()->isHigherScoreBetter());
    bool use_all_hits = param_.getValue("use_all_hits").toBool();

    ScoreToTgtDecLabelPairs scores_labels;
    IDScoreGetterSetter::getScores_(scores_labels, ids, [&identifier](const PeptideIdentification& id){return identifier == id.getIdentifier();}, use_all_hits);

    if (scores_labels.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No scores could be extracted!");
    }

    if (higher_score_better)
    { // decreasing
      std::sort(scores_labels.rbegin(), scores_labels.rend());
    }
    else
    { // increasing
      std::sort(scores_labels.begin(), scores_labels.end());
    }
    // if fp_cutoff is zero do the full AUC.
    return rocN(scores_labels, fp_cutoff == 0 ? scores_labels.size() : fp_cutoff);
  }

  double FalseDiscoveryRate::rocN(const ConsensusMap& ids, Size fp_cutoff, bool include_unassigned_peptides) const
  {
    bool higher_score_better(false);
    // Check first ID in a feature for the score orientation.
    for (const auto& f : ids)
    {
      const auto& pepids = f.getPeptideIdentifications();
      if (!pepids.empty())
      {
        higher_score_better = pepids[0].isHigherScoreBetter();
        break;
      }
    }
    bool use_all_hits = param_.getValue("use_all_hits").toBool();

    ScoreToTgtDecLabelPairs scores_labels;
    IDScoreGetterSetter::getPeptideScoresFromMap_(scores_labels, ids, include_unassigned_peptides, use_all_hits);

    if (scores_labels.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No scores could be extracted!");
    }

    if (higher_score_better)
    { // decreasing
      std::sort(scores_labels.rbegin(), scores_labels.rend());
    }
    else
    { // increasing
      std::sort(scores_labels.begin(), scores_labels.end());
    }
    // if fp_cutoff is zero do the full AUC.
    return rocN(scores_labels, fp_cutoff == 0 ? scores_labels.size() : fp_cutoff);
  }

  double FalseDiscoveryRate::rocN(const ConsensusMap& ids, Size fp_cutoff, const String& identifier, bool include_unassigned_peptides) const
  {
    bool higher_score_better(ids[0].getPeptideIdentifications().begin()->isHigherScoreBetter());
    bool use_all_hits = param_.getValue("use_all_hits").toBool();

    ScoreToTgtDecLabelPairs scores_labels;
    IDScoreGetterSetter::getPeptideScoresFromMap_(scores_labels, ids, include_unassigned_peptides, [&identifier](const PeptideIdentification& id){return identifier == id.getIdentifier();}, use_all_hits);

    if (scores_labels.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No scores could be extracted!");
    }

    if (higher_score_better)
    { // decreasing
      std::sort(scores_labels.rbegin(), scores_labels.rend());
    }
    else
    { // increasing
      std::sort(scores_labels.begin(), scores_labels.end());
    }
    // if fp_cutoff is zero do the full AUC.
    return rocN(scores_labels, fp_cutoff == 0 ? scores_labels.size() : fp_cutoff);
  }

  void FalseDiscoveryRate::applyBasic(ConsensusMap & cmap, bool include_unassigned_peptides)
  {
    bool q_value = !param_.getValue("no_qvalues").toBool();
    const string& score_type = q_value ? "q-value" : "FDR";
    bool all_hits = param_.getValue("use_all_hits").toBool();

    bool treat_runs_separately = param_.getValue("treat_runs_separately").toBool();
    bool split_charge_variants = param_.getValue("split_charge_variants").toBool();

    //TODO this assumes all used search engine scores have the same score orientation
    // include the determination of orientation in the getScores methods instead
    bool higher_score_better = false;
    for (const auto& cf : cmap)
    {
      const auto& pep_ids = cf.getPeptideIdentifications();
      if (!pep_ids.empty())
      {
        higher_score_better = pep_ids[0].isHigherScoreBetter();
        break;
      }
    }
    if (cmap.empty())
    {
      for (const auto& id : cmap.getUnassignedPeptideIdentifications())
      {
        higher_score_better = id.isHigherScoreBetter();
        break;
      }
    }

    bool add_decoy_peptides = param_.getValue("add_decoy_peptides").toBool();
    ScoreToTgtDecLabelPairs scores_labels;

    //Warning: this assumes that there are no dangling identifier references in the PeptideIDs
    // because this disables checking
    if (cmap.getProteinIdentifications().size() == 1)
    {
      treat_runs_separately = false;
    }

    if (treat_runs_separately)
    {
      for (const auto& protID : cmap.getProteinIdentifications())
      {
        if (split_charge_variants)
        {
          pair<int, int> chargeRange = protID.getSearchParameters().getChargeRange();
          for (int c = chargeRange.first; c <= chargeRange.second; ++c)
          {
            if (c == 0) continue;
            IDScoreGetterSetter::getPeptideScoresFromMap_(scores_labels, cmap, include_unassigned_peptides,
             [&protID](const PeptideIdentification& id){return protID.getIdentifier() == id.getIdentifier();}, all_hits,
             [&c](const PeptideHit& hit){return c == hit.getCharge();});
            map<double, double> scores_to_fdr;
            calculateFDRBasic_(scores_to_fdr, scores_labels, q_value, higher_score_better);
            IDScoreGetterSetter::setPeptideScoresForMap_(scores_to_fdr, cmap, include_unassigned_peptides, score_type, higher_score_better, add_decoy_peptides, c, protID.getIdentifier());
          }
        }
        else
        {
          IDScoreGetterSetter::getPeptideScoresFromMap_(scores_labels, cmap, include_unassigned_peptides, [&protID](const PeptideIdentification& id){return protID.getIdentifier() == id.getIdentifier();}, all_hits);
          map<double, double> scores_to_fdr;
          calculateFDRBasic_(scores_to_fdr, scores_labels, q_value, higher_score_better);
          IDScoreGetterSetter::setPeptideScoresForMap_(scores_to_fdr, cmap, include_unassigned_peptides, score_type, higher_score_better, add_decoy_peptides, protID.getIdentifier());
        }
      }
    }
    else
    {
      IDScoreGetterSetter::getPeptideScoresFromMap_(scores_labels, cmap, include_unassigned_peptides, all_hits);
      map<double, double> scores_to_fdr;
      calculateFDRBasic_(scores_to_fdr, scores_labels, q_value, higher_score_better);
      IDScoreGetterSetter::setPeptideScoresForMap_(scores_to_fdr, cmap, include_unassigned_peptides, score_type, higher_score_better, add_decoy_peptides);
    }
  }

  //TODO Add another overload that iterates over a vector. to be consistent with old interface
  //TODO Make it return a double for the AUC or max FDR (i.e., without cutoff)
  void FalseDiscoveryRate::applyBasic(ProteinIdentification & id, bool groups_too)
  {
    bool add_decoy_proteins = param_.getValue("add_decoy_proteins").toBool();

    bool q_value = !param_.getValue("no_qvalues").toBool();
    //TODO Check naming conventions. Ontology? Make class member?
    const string& score_type = q_value ? "q-value" : "FDR";
    bool higher_score_better(id.isHigherScoreBetter());

    ScoreToTgtDecLabelPairs scores_labels;
    scores_labels.reserve(id.getHits().size());
    std::map<double,double> scores_to_FDR;

    // TODO this could be a separate function.. And it could actually be sped up.
    //  We could store the number of decoys/targets in the group, or we only update the
    //  scores of proteins that are actually in groups (rest stays the same)
    // do groups first, if keep_decoy is false, we would otherwise miss those proteins
    if (groups_too)
    {
      // Prepare lookup map for decoy proteins (since there is no direct way back from group to protein)
      // TODO we could also require the decoy affix to be specified
      unordered_set<string> decoy_accs;
      for (const auto& prot : id.getHits())
      {
        if (!prot.metaValueExists("target_decoy") || prot.getMetaValue("target_decoy") == "decoy")
        {
          decoy_accs.insert(prot.getAccession());
        }
      }
      IDScoreGetterSetter::getScores_(scores_labels, id.getIndistinguishableProteins(), decoy_accs);
      calculateFDRBasic_(scores_to_FDR, scores_labels, q_value, higher_score_better);
      if (!scores_labels.empty())
        IDScoreGetterSetter::setScores_(scores_to_FDR, id.getIndistinguishableProteins(), score_type, false);
    }

    scores_to_FDR.clear();
    scores_labels.clear();
    scores_labels.reserve(id.getHits().size());

    IDScoreGetterSetter::getScores_(scores_labels, id);
    if (scores_labels.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No scores could be extracted!");
    }
    calculateFDRBasic_(scores_to_FDR, scores_labels, q_value, higher_score_better);
    if (!scores_labels.empty())
    {
      IDScoreGetterSetter::setScores_(scores_to_FDR, id, score_type, false, add_decoy_proteins);
    }
    else
    {
     OPENMS_LOG_WARN << "Warning: No scores could be extracted for proteins. No FDR calculation performed.";
    }

    scores_to_FDR.clear();
  }

  void FalseDiscoveryRate::applyBasic(const std::vector<ProteinIdentification> & run_info, std::vector<PeptideIdentification> & ids)
  {
    if (ids.empty()) return;
    bool treat_runs_separately = param_.getValue("treat_runs_separately").toBool();
    bool split_charge_variants = param_.getValue("split_charge_variants").toBool();
    // TODO decide on param interface. Consolidate with FDR tool parameters.
    //  Then pass the bool to the applyBasic calls.
    //bool best_per_pep = param_.getValue("pepFDR").toBool();

    String identifier = "";
    bool higher_score_better = true;
    if (treat_runs_separately)
    {
      for (const auto& run : run_info)
      {
        identifier = run.getIdentifier();
        for (const auto& pepid : ids)
        {
          if (pepid.getIdentifier() == identifier)
          {
            higher_score_better = pepid.isHigherScoreBetter();
            break;
          }
        }
        if (split_charge_variants)
        {
          pair<int, int> chargeRange = run.getSearchParameters().getChargeRange();
          for (int c = chargeRange.first; c <= chargeRange.second; ++c)
          {
            if (c == 0) continue;
            applyBasic(ids, higher_score_better, c, identifier);
          }
        }
        else
        {
          applyBasic(ids, higher_score_better, 0, identifier);
        }
        
      }
    }
    else if (split_charge_variants)
    {
      pair<int, int> chargeRange = {10000,-10000};
      for (const auto& run : run_info)
      {
        chargeRange.first = std::min(run.getSearchParameters().getChargeRange().first, chargeRange.first);
        chargeRange.second = std::max(run.getSearchParameters().getChargeRange().first, chargeRange.second);
      }
      higher_score_better = ids[0].isHigherScoreBetter();
      for (int c = chargeRange.first; c <= chargeRange.second; ++c)
      {
        if (c == 0) continue;
        applyBasic(ids, higher_score_better, c);
      }
    }
    else // altogether
    {
      higher_score_better = ids[0].isHigherScoreBetter();
      applyBasic(ids, higher_score_better);
    }
  }

  void FalseDiscoveryRate::applyBasicPeptideLevel(ConsensusMap & map, bool include_unassigned)
  {
    bool q_value = !param_.getValue("no_qvalues").toBool();
    //TODO Check naming conventions. Ontology?
    const string& score_type = q_value ? "peptide q-value" : "peptide FDR";
    bool add_decoy_peptides = param_.getValue("add_decoy_peptides").toBool();
    // since we do not support multiple runs here yet, we take the orientation of the first ID
    bool higher_better = true;
    for (const auto& f : map)
    {
      if (!f.getPeptideIdentifications().empty())
      {
        higher_better = f.getPeptideIdentifications()[0].isHigherScoreBetter();
      }
    }

    unordered_map<String, ScoreToTgtDecLabelPair> seq_to_score_labels;
    IDScoreGetterSetter::fillPeptideScoreMap_(seq_to_score_labels, map, include_unassigned);

    ScoreToTgtDecLabelPairs pairs;
    for (auto const & seq_to_score_label : seq_to_score_labels)
    {
      pairs.push_back(seq_to_score_label.second);
    }
    std::map<double,double> score_to_fdr;
    calculateFDRBasic_(score_to_fdr, pairs, q_value, higher_better);
    // convert scores in unordered map to FDR/qvalues
    for (auto & seq_to_score_label : seq_to_score_labels)
    {
      if (higher_better)
      {
        auto ub = score_to_fdr.upper_bound(seq_to_score_label.second.first);
        if (ub != score_to_fdr.begin()) ub--;
        seq_to_score_label.second.first = ub->second;
      }
      else
      {
        seq_to_score_label.second.first = score_to_fdr.lower_bound(seq_to_score_label.second.first)->second;
      }
    }

    IDScoreGetterSetter::setPeptideScoresFromMap_(seq_to_score_labels, map, score_type, add_decoy_peptides, include_unassigned);
  }

  void FalseDiscoveryRate::applyBasicPeptideLevel(std::vector<PeptideIdentification> & ids)
  {
    bool q_value = !param_.getValue("no_qvalues").toBool();
    //TODO Check naming conventions. Ontology?
    const string& score_type = q_value ? "peptide q-value" : "peptide FDR";
    bool add_decoy_peptides = param_.getValue("add_decoy_peptides").toBool();
    // since we do not support multiple runs here yet, we take the orientation of the first ID
    bool higher_better = ids[0].isHigherScoreBetter();

    unordered_map<String, ScoreToTgtDecLabelPair> seq_to_score_labels;
    IDScoreGetterSetter::fillPeptideScoreMap_(seq_to_score_labels, ids);

    ScoreToTgtDecLabelPairs pairs;
    for (auto const & seq_to_score_label : seq_to_score_labels)
    {
      pairs.push_back(seq_to_score_label.second);
    }
    map<double,double> score_to_fdr;
    calculateFDRBasic_(score_to_fdr, pairs, q_value, higher_better);
    // convert scores in unordered map to FDR/qvalues
    for (auto & seq_to_score_label : seq_to_score_labels)
    {
      if (higher_better)
      {
        auto ub = score_to_fdr.upper_bound(seq_to_score_label.second.first);
        if (ub != score_to_fdr.begin()) ub--;
        seq_to_score_label.second.first = ub->second;
      }
      else
      {
        seq_to_score_label.second.first = score_to_fdr.lower_bound(seq_to_score_label.second.first)->second;
      }
    }

    IDScoreGetterSetter::setPeptideScoresFromMap_(seq_to_score_labels, ids, score_type, add_decoy_peptides);
  }

  // TODO why again do we need higher_score_better here?
  void FalseDiscoveryRate::applyBasic(std::vector<PeptideIdentification> & ids, bool higher_score_better, int charge, String identifier, bool only_best_per_pep)
  {
    bool q_value = !param_.getValue("no_qvalues").toBool();
    //TODO Check naming conventions. Ontology?
    const string& score_type = q_value ? "q-value" : "FDR";

    bool use_all_hits = param_.getValue("use_all_hits").toBool();

    bool add_decoy_peptides = param_.getValue("add_decoy_peptides").toBool();

    ScoreToTgtDecLabelPairs scores_labels;
    std::map<double,double> scores_to_FDR;
    auto idcheck = [&identifier](const PeptideIdentification& id){return identifier == id.getIdentifier();};
    
    if (charge == 0 && !only_best_per_pep)
    {
      if (identifier.empty())
      {
        IDScoreGetterSetter::getScores_(scores_labels, ids, use_all_hits);
      }
      else
      {
        IDScoreGetterSetter::getScores_(scores_labels, ids, idcheck, use_all_hits);
      }
    }
    else
    {
      if (!only_best_per_pep /*&& charge != 0*/)
      {
        if (identifier.empty())
        {
          IDScoreGetterSetter::getScores_(scores_labels, ids, use_all_hits, [&charge](const PeptideHit& hit){return charge == hit.getCharge();});
        }
        else
        {
          IDScoreGetterSetter::getScores_(scores_labels, ids, idcheck, use_all_hits, [&charge](const PeptideHit& hit){return charge == hit.getCharge();});
        }
      }
      else if (charge == 0 /* && only_best_per_pep */)
      {
        if (identifier.empty())
        {
          IDScoreGetterSetter::getScores_(scores_labels, ids, use_all_hits, [](const PeptideHit& hit){return hit.metaValueExists("best_per_peptide") && int(hit.getMetaValue("best_per_peptide")) == 1;});
        }
        else
        {
          IDScoreGetterSetter::getScores_(scores_labels, ids, idcheck, use_all_hits, [](const PeptideHit& hit){return hit.metaValueExists("best_per_peptide") && int(hit.getMetaValue("best_per_peptide")) == 1;});
        }
      }
      else /*charge != 0 && only_best_per_pep */
      {
        if (identifier.empty())
        {
          IDScoreGetterSetter::getScores_(scores_labels, ids, use_all_hits, [&charge](const PeptideHit& hit){return (charge == hit.getCharge()) && hit.metaValueExists("best_per_peptide") && int(hit.getMetaValue("best_per_peptide")) == 1;});
        }
        else
        {
          IDScoreGetterSetter::getScores_(scores_labels, ids, idcheck, use_all_hits, [&charge](const PeptideHit& hit){return (charge == hit.getCharge()) && hit.metaValueExists("best_per_peptide") && int(hit.getMetaValue("best_per_peptide")) == 1;});
        }
      }

    }
    
    if (scores_labels.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No scores could be extracted for FDR!");
    }
    calculateFDRBasic_(scores_to_FDR, scores_labels, q_value, higher_score_better);
    if (!scores_labels.empty())
      IDScoreGetterSetter::setScores_<PeptideIdentification>(scores_to_FDR, ids, score_type, false, add_decoy_peptides);
    scores_to_FDR.clear();
  }

  //TODO could be implemented for PeptideIDs, too
  //TODO iterate over the vector. to be consistent with old interface
  void FalseDiscoveryRate::applyEstimated(std::vector<ProteinIdentification> &ids) const
  {
    //Note: this is actually unused because I think with that approach you will always get q-values.
    //bool q_value = !param_.getValue("no_qvalues").toBool();
    bool higher_score_better(ids.begin()->isHigherScoreBetter());

    bool add_decoy_proteins = param_.getValue("add_decoy_proteins").toBool();

    //TODO not yet supported (if ever)
    //bool treat_runs_separately = param_.getValue("treat_runs_separately").toBool();
    if (ids.size() > 1)
    {
     OPENMS_LOG_WARN << "More than one set of ProteinIdentifications found. Only using the first one for FDR calculation.\n";
    }

    if (ids[0].getScoreType() != "Posterior Probability" && ids[0].getScoreType() != "Posterior Error Probability")
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Proteins in ProteinIdentification do not have a posterior (error) probability assigned. Please run an inference first.", ids[0].getScoreType());
    }

    ScoreToTgtDecLabelPairs scores_labels;
    std::map<double,double> scores_to_FDR;
    //TODO actually we do not need the labels for estimated FDR and it currently fails if we do not have TD annotations
    //TODO maybe separate getScores and getScoresAndLabels
    IDScoreGetterSetter::getScores_(scores_labels, ids[0]);
    calculateEstimatedQVal_(scores_to_FDR, scores_labels, higher_score_better);
    if (!scores_labels.empty())
      IDScoreGetterSetter::setScores_(scores_to_FDR, ids[0], "Estimated Q-Values", false, add_decoy_proteins);
  }


  //TODO remove?
  double FalseDiscoveryRate::applyEvaluateProteinIDs(const std::vector<ProteinIdentification>& ids, double pepCutoff, UInt fpCutoff, double diffWeight) const
  {
    //TODO not yet supported (if ever)
    //bool treat_runs_separately = param_.getValue("treat_runs_separately").toBool();
    if (ids.size() > 1)
    {
     OPENMS_LOG_WARN << "More than one set of ProteinIdentifications found. Only using the first one for calculation.\n";
    }

    if (ids[0].getScoreType() != "Posterior Probability")
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Proteins in ProteinIdentification do not have a posterior probability assigned. Please run an inference first.", ids[0].getScoreType());
    }

    ScoreToTgtDecLabelPairs scores_labels;
    IDScoreGetterSetter::getScores_(scores_labels, ids[0]);
    std::sort(scores_labels.rbegin(), scores_labels.rend());
    return diffEstimatedEmpirical(scores_labels, pepCutoff) * diffWeight +
        rocN(scores_labels, fpCutoff) * (1 - diffWeight);
  }

  double FalseDiscoveryRate::applyEvaluateProteinIDs(const ProteinIdentification& ids, double pepCutoff, UInt fpCutoff, double diffWeight) const
  {
    if (ids.getScoreType() != "Posterior Probability")
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Proteins in ProteinIdentification do not have a posterior probability assigned. Please run an inference first.", ids.getScoreType());
    }

    ScoreToTgtDecLabelPairs scores_labels;
    IDScoreGetterSetter::getScores_(scores_labels, ids);
    std::sort(scores_labels.rbegin(), scores_labels.rend());
    double diff = diffEstimatedEmpirical(scores_labels, pepCutoff);
    double auc = rocN(scores_labels, fpCutoff);
    OPENMS_LOG_INFO << "Evaluation of protein probabilities: Difference estimated vs. T-D FDR = " << diff << " and roc" << fpCutoff << " = " << auc << std::endl;
    // we want the score to get higher the lesser the difference. Subtract from one.
    // Then convex combination with the AUC.
    return (1.0 - diff) * (1.0 - diffWeight) + auc * diffWeight;
  }

  double FalseDiscoveryRate::applyEvaluateProteinIDs(ScoreToTgtDecLabelPairs& scores_labels, double pepCutoff, UInt fpCutoff, double diffWeight) const
  {
    std::sort(scores_labels.rbegin(), scores_labels.rend());
    double diff = diffEstimatedEmpirical(scores_labels, pepCutoff);
    double auc = rocN(scores_labels, fpCutoff);
    OPENMS_LOG_INFO << "Evaluation of protein probabilities: Difference estimated vs. T-D FDR = " << diff << " and roc" << fpCutoff << " = " << auc << std::endl;
    // we want the score to get higher the lesser the difference. Subtract from one.
    // Then convex combination with the AUC.
    return (1.0 - diff) * (1.0 - diffWeight) + auc * diffWeight;
  }

  //TODO this probably could work on group level, too, but only if peptide-level decoys were used, such that
  // decoys are indistinguishable iff targets are indistinguishable
  void FalseDiscoveryRate::applyPickedProteinFDR(ProteinIdentification & id, String decoy_string, bool prefix, bool groups_too)
  {
    bool add_decoy_proteins = param_.getValue("add_decoy_proteins").toBool();
    bool q_value = !param_.getValue("no_qvalues").toBool();
    //TODO Check naming conventions. Ontology?
    const string& score_type = q_value ? "q-value" : "FDR";

    //TODO this assumes all runs have the same ordering! Otherwise do it per identifier.
    bool higher_score_better(id.isHigherScoreBetter());

    if (decoy_string.empty())
    {
      auto r = FalseDiscoveryRate::DecoyStringHelper::findDecoyString(id);
      if (!r.success)
      {
        r.is_prefix = true;
        r.name = "DECOY_";
        OPENMS_LOG_WARN << "Unable to determine decoy string automatically (not enough decoys were detected)! Using default " << (r.is_prefix ? "prefix" : "suffix") << " decoy string '" << r.name << "'\n"
        << "If you think that this is incorrect, please provide a decoy_string and its position manually!" << std::endl;
      }
      prefix = r.is_prefix;
      decoy_string = r.name;
      // decoy string and position was extracted successfully
      OPENMS_LOG_INFO << "Using " << (prefix ? "prefix" : "suffix") << " decoy string '" << decoy_string << "'" << std::endl;
    }

    ScoreToTgtDecLabelPairs scores_labels;
    std::map<double,double> scores_to_FDR;
    std::unordered_map<String, ScoreToTgtDecLabelPair> picked_scores;
    IDScoreGetterSetter::getPickedProteinScores_(picked_scores, id, decoy_string, prefix);
    scores_labels.reserve(picked_scores.size());

    if (groups_too)
    {
      IDScoreGetterSetter::getPickedProteinGroupScores_(picked_scores, scores_labels, id.getIndistinguishableProteins(), decoy_string, prefix);
      calculateFDRBasic_(scores_to_FDR, scores_labels, q_value, higher_score_better);
      IDScoreGetterSetter::setScores_(scores_to_FDR, id.getIndistinguishableProteins(), score_type, false);
      scores_to_FDR.clear();
      scores_labels.clear();
    }

    // for single proteins just take all scores
    for (auto& kv : picked_scores)
    {
      scores_labels.emplace_back(std::move(kv.second)); // move all. We do not need them anymore
    }

    if (scores_labels.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No scores could be extracted for FDR calculation!");
    }
    calculateFDRBasic_(scores_to_FDR, scores_labels, q_value, higher_score_better);
    IDScoreGetterSetter::setScores_(scores_to_FDR, id, score_type, false, add_decoy_proteins);
    scores_to_FDR.clear();
    scores_labels.clear();
  }

  //TODO the following two methods assume sortedness. Add precondition and/or doxygen comment
  double FalseDiscoveryRate::diffEstimatedEmpirical(const ScoreToTgtDecLabelPairs& scores_labels, double pepCutoff) const
  {
    bool conservative = param_.getValue("conservative").toBool();
    if (scores_labels.empty())
    {
     OPENMS_LOG_WARN << "Warning: No scores extracted for FDR calculation. Skipping. Do you have target-decoy annotated Hits?" << std::endl;
      return 1.0;
    }

    double diffArea = 0.0;
    double est = 0.0, estPrev = 0.0, emp = 0.0, empPrev = 0.0;
    double pepSum = 0.0;
    UInt truePos = 0u, falsePos = 0u; //, truePosPrev = 0u, falsePosPrev = 0u;

    auto pit = scores_labels.cbegin();
    for (; pit != scores_labels.end()-1; ++pit)
    {
      pit->second ? truePos++ : falsePos++;
      pepSum += (1 - pit->first);

      //Look ahead. Running variables have already been incremented
      if ((pit+1)->first != pit->first)
      {
        est = pepSum / (truePos + falsePos);
        if (conservative)
        {
          if (truePos == 0.)
          {
            emp = 1.;
          }
          else
          {
            emp = static_cast<double>(falsePos) / (truePos);
          }
        }
        else
        {
          emp = static_cast<double>(falsePos) / (truePos+falsePos);
        }

        diffArea += trapezoidal_area_xEqy(estPrev, est, empPrev, emp);

        //truePosPrev = truePos;
        //falsePosPrev = falsePos;
        estPrev = est;
        empPrev = emp;
      }
    }

    //Last item. Always add areas there
    pit->second ? truePos++ : falsePos++;
    pepSum += (1 - pit->first);
    est = pepSum / (truePos + falsePos);
    emp = static_cast<double>(falsePos) / (truePos + falsePos);
    diffArea += trapezoidal_area_xEqy(estPrev, est, empPrev, emp);

    //scale by max PEP value achievable (= width); height = empFDR can be 1.0
    diffArea /= std::min(est, pepCutoff);

    return diffArea;
  }

  double FalseDiscoveryRate::rocN(const ScoreToTgtDecLabelPairs& scores_labels, Size fpCutoff) const
  {
    if (scores_labels.empty())
    {
     OPENMS_LOG_WARN << "Warning: No scores extracted for FDR calculation. Skipping. Do you have target-decoy annotated Hits?" << std::endl;
      return 0.0;
    }

    double rocN = 0.0;
    UInt truePos = 0u, falsePos = 0u, truePosPrev = 0u, falsePosPrev = 0u;

    auto pit = scores_labels.cbegin();
    for (; pit != scores_labels.cend()-1; ++pit)
    {
      pit->second ? truePos++ : falsePos++;

      //Look ahead. Running variables have already been incremented
      if ((pit+1)->first != pit->first)
      {
        rocN += trapezoidal_area(falsePos, falsePosPrev, truePos, truePosPrev);
        if (falsePos >= fpCutoff) return rocN / (falsePos * truePos); //if with the last batch, you have >= N FPs, return scaled
        truePosPrev = truePos;
        falsePosPrev = falsePos;
      }
    }

    //Last item if not returned. Always add areas there
    pit->second ? truePos++ : falsePos++;
    rocN += trapezoidal_area(falsePos, falsePosPrev, truePos, truePosPrev);

    if (falsePos == 0) return 1;
    return rocN / (falsePos * truePos);
  }

  /// x2 has to be bigger than x1,
  /// handles possible intersections
  double FalseDiscoveryRate::trapezoidal_area_xEqy(double x1, double x2, double y1, double y2) const
  {
    double height = x2 - x1;
    double b1 = y1 - x1;
    double b2 = y2 - x2;

    if (std::signbit(b1) == std::signbit(b2))
    {
      return (std::fabs(b1) + std::fabs(b2)) * height / 2.0;
    }
    else
    {
      // it is intersecting the x=y line. Add the area of the resulting triangles
      return (b1*b1 + b2*b2) * height / (2.0 * (std::fabs(b1)+std::fabs(b2)));
    }
  }

  /// assumes a flat base
  double FalseDiscoveryRate::trapezoidal_area(double x1, double x2, double y1, double y2) const
  {
    double base = fabs(x1 - x2);
    double avgHeight = (y1+y2)/2.0;
    return base * avgHeight;
  }


  // Actually this does not need the bool entries in the scores_labels, but leads to less code
  // Assumes P(E)Probabilities as scores
  void FalseDiscoveryRate::calculateEstimatedQVal_(std::map<double, double> &scores_to_FDR,
                                                   ScoreToTgtDecLabelPairs &scores_labels,
                                                   bool higher_score_better) const
  {
    if (scores_labels.empty())
    {
     OPENMS_LOG_WARN << "Warning: No scores extracted for FDR calculation. Skipping. Do you have target-decoy annotated Hits?" << std::endl;
      return;
    }

    if (higher_score_better)
    { // decreasing
      std::sort(scores_labels.rbegin(), scores_labels.rend());
    }
    else
    { // increasing
      std::sort(scores_labels.begin(), scores_labels.end());
    }

    //TODO I think we can just do it "in-place" to save space
    std::vector<double> estimatedFDR;
    estimatedFDR.reserve(scores_labels.size());

    // Basically a running average
    double sum = 0.0;

    for (size_t j = 0; j < scores_labels.size(); ++j)
    {
      sum += scores_labels[j].first;
      estimatedFDR[j] = sum / (double(j)+1.0);
    }

    if (higher_score_better) // Transform to PEP
    {
      std::transform(estimatedFDR.begin(), estimatedFDR.end(), estimatedFDR.begin(), [&](double d) { return 1 - d; });
    }

    // In case of multiple equal scores, this will add the _last_ fdr that it finds. Since fdrs are decreasing
    // in either way in estimatedFDR, this adds the minimal FDR to the map for this score.
    std::transform(scores_labels.begin(), scores_labels.end(), estimatedFDR.begin(), std::inserter(scores_to_FDR, scores_to_FDR.begin()), [&](std::pair<double,bool> sl, double fdr){return make_pair(sl.first, fdr);});
  }

  /*
   void FalseDiscoveryRate::calculateFDRBasic_(
    std::map<double,double>& scores_to_FDR,
    ScoreToTgtDecLabelPairs& scores_labels,
    std::vector<size_t>& ordering,
    bool qvalue,
    bool higher_score_better) const
  {

  }

  void FalseDiscoveryRate::calculateFDRBasicOnSorted_(
  std::map<double,double>& scores_to_FDR,
  ScoreToTgtDecLabelPairs& scores_labels,
  bool qvalue,
  bool higher_score_better) const
  {

  }*/

  void FalseDiscoveryRate::calculateFDRBasic_(
      std::map<double,double>& scores_to_FDR,
      ScoreToTgtDecLabelPairs& scores_labels,
      bool qvalue,
      bool higher_score_better) const
  {
    //TODO put in separate function to avoid ifs in iteration
    bool conservative = param_.getValue("conservative").toBool();
    if (scores_labels.empty())
    {
      OPENMS_LOG_WARN << "Warning: No scores extracted for FDR calculation. Skipping. Do you have target-decoy annotated Hits?" << std::endl;
      return;
    }

    if (higher_score_better)
    { // decreasing
      std::sort(scores_labels.rbegin(), scores_labels.rend());
    }
    else
    { // increasing
      std::sort(scores_labels.begin(), scores_labels.end());
    }

    //uniquify scores and add decoy proportions
    double decoys = 0.; // double to account for "partial" decoys
    double last_score = scores_labels[0].first;

    size_t j = 0;
    for (; j < scores_labels.size(); ++j)
    {
      //Although we do not really care about equality we compare with tolerance to make it (more?) compiler independent.
      if (std::abs(scores_labels[j].first - last_score) > 1e-12)
      {
        #ifdef FALSE_DISCOVERY_RATE_DEBUG
        std::cerr << "Recording score: " << last_score << " with " << decoys << " decoys at index+1 = " << (j+1) << " -> fdr: " << decoys/(j+1.0) << std::endl;
        #endif
        //we are using the conservative formula (Decoy + 1) / (Tgts)
        if (conservative)
        {
          scores_to_FDR[last_score] = (decoys+1.0)/(double(j)+1.0-decoys);
        }
        else
        {
          scores_to_FDR[last_score] = (decoys+1.0)/(double(j)+1.0);
        }

        last_score = scores_labels[j].first;
      }

      decoys += 1. - scores_labels[j].second;

      /* The following was for the binary interpretation. Now we allow for partial decoy contributions as above
      if (!scores_labels[j].second)
      {
        decoys++;
      }
      */
    }

    // in case there is only one score and generally to include the last score, I guess we need to do this
    if (conservative)
    {
      scores_to_FDR[last_score] = (decoys+1.0)/(double(j)+1.0-decoys);
    }
    else
    {
      scores_to_FDR[last_score] = (decoys+1.0)/(double(j)+1.0);
    }

    if (qvalue) //apply a cumulative minimum on the map (from low to high fdrs)
    {
      double cummin = 1.0;

      if (higher_score_better)
      {
        for (auto&& rit = scores_to_FDR.begin(); rit != scores_to_FDR.end(); ++rit)
        {
        #ifdef FALSE_DISCOVERY_RATE_DEBUG
          std::cerr << "Comparing " << rit->second << " to " << cummin << std::endl;
        #endif
          cummin = std::min(rit->second, cummin);
          rit->second = cummin;
        }
      }
      else
      {
        for (auto&& rit = scores_to_FDR.rbegin(); rit != scores_to_FDR.rend(); ++rit)
        {
        #ifdef FALSE_DISCOVERY_RATE_DEBUG
          std::cerr << "Comparing " << rit->second << " to " << cummin << std::endl;
        #endif
          cummin = std::min(rit->second, cummin);
          rit->second = cummin;
        }
      }
    }
  }

  using DecoyStringToAffixCount = std::unordered_map<std::string, std::pair<Size, Size>>;
  using CaseInsensitiveToCaseSensitiveDecoy = std::unordered_map<std::string, std::string>;
  /**
    @brief Heuristic to determine the decoy string given a set of protein names

    Tested decoy strings are "decoy", "dec", "reverse", "rev", "__id_decoy", "xxx", "shuffled", "shuffle", "pseudo" and "random".
    Both prefix and suffix is tested and if one of the candidates above is found in at least 40% of all proteins,
    it is returned as the winner (see DecoyHelper::Result).
  */
   FalseDiscoveryRate::DecoyStringHelper::Result FalseDiscoveryRate::DecoyStringHelper::findDecoyString(const ProteinIdentification& proteins)
  {
    // common decoy strings in FASTA files
    // note: decoy prefixes/suffices must be provided in lower case
    static const std::vector<std::string> affixes{ "decoy", "dec", "reverse", "rev", "reversed", "__id_decoy", "xxx", "shuffled", "shuffle", "pseudo", "random" };

    // map decoys to counts of occurrences as prefix/suffix
    DecoyStringToAffixCount decoy_count;
    // map case-insensitive strings back to original case (as used in fasta)
    CaseInsensitiveToCaseSensitiveDecoy decoy_case_sensitive;

    // setup prefix- and suffix regex strings
    // TODO extend regex to allow skipping the underscore? i.e. with "?"
    const std::string regexstr_prefix = std::string("^(") + ListUtils::concatenate<std::string>(affixes, "_*|") + "_*)";
    const std::string regexstr_suffix = std::string("(_") + ListUtils::concatenate<std::string>(affixes, "*|_") + ")$";

    // setup regexes
    const boost::regex pattern_prefix(regexstr_prefix);
    const boost::regex pattern_suffix(regexstr_suffix);

    Size all_prefix_occur(0), all_suffix_occur(0), all_proteins_count(0);

    for (const auto& prot : proteins.getHits())
    {
      all_proteins_count += 1;

      boost::smatch sm;
      const String& seq = prot.getAccession();

      String seq_lower = seq;
      seq_lower.toLower();

      // search for prefix
      bool found_prefix = boost::regex_search(seq_lower, sm, pattern_prefix);
      if (found_prefix)
      {
        std::string match = sm[0];
        all_prefix_occur++;

        // increase count of observed prefix
        decoy_count[match].first++;

        // store observed (case-sensitive and with special characters)
        std::string seq_decoy = StringUtils::prefix(seq, match.length());
        decoy_case_sensitive[match] = seq_decoy;
      }

      // search for suffix
      bool found_suffix = boost::regex_search(seq_lower, sm, pattern_suffix);
      if (found_suffix)
      {
        std::string match = sm[0];
        all_suffix_occur++;

        // increase count of observed suffix
        decoy_count[match].second++;

        // store observed (case-sensitive and with special characters)
        std::string seq_decoy = StringUtils::suffix(seq, match.length());
        decoy_case_sensitive[match] = seq_decoy;
      }
    }

    // DEBUG ONLY: print counts of found decoys
    for (auto &a : decoy_count)
    {
      OPENMS_LOG_DEBUG << a.first << "\t" << a.second.first << "\t" << a.second.second << std::endl;
    }

    // less than 30% of proteins are decoys -> won't be able to determine a decoy string and its position
    // return default values
    if (static_cast<double>(all_prefix_occur + all_suffix_occur) < 0.3 * static_cast<double>(all_proteins_count))
    {
      OPENMS_LOG_ERROR << "Unable to determine decoy string (not enough occurrences; <30%)!" << std::endl;
      return {false, "?", true};
    }

    if (all_prefix_occur == all_suffix_occur)
    {
      OPENMS_LOG_ERROR << "Unable to determine decoy string (prefix and suffix occur equally often)!" << std::endl;
      return {false, "?", true};
    }

    // Decoy prefix occurred at least 80% of all prefixes + observed in at least 30% of all proteins -> set it as prefix decoy
    for (const auto& pair : decoy_count)
    {
      const std::string & case_insensitive_decoy_string = pair.first;
      const std::pair<Size, Size>& prefix_suffix_counts = pair.second;
      double freq_prefix = static_cast<double>(prefix_suffix_counts.first) / static_cast<double>(all_prefix_occur);
      double freq_prefix_in_proteins = static_cast<double>(prefix_suffix_counts.first) / static_cast<double>(all_proteins_count);

      if (freq_prefix >= 0.8 && freq_prefix_in_proteins >= 0.3)
      {
        if (prefix_suffix_counts.first != all_prefix_occur)
        {
          OPENMS_LOG_WARN << "More than one decoy prefix observed!" << std::endl;
          OPENMS_LOG_WARN << "Using most frequent decoy prefix (" << (int)(freq_prefix * 100) << "%)" << std::endl;
        }

        return { true, decoy_case_sensitive[case_insensitive_decoy_string], true};
      }
    }

    // Decoy suffix occurred at least 80% of all suffixes + observed in at least 30% of all proteins -> set it as suffix decoy
    for (const auto& pair : decoy_count)
    {
      const std::string& case_insensitive_decoy_string = pair.first;
      const std::pair<Size, Size>& prefix_suffix_counts = pair.second;
      double freq_suffix = static_cast<double>(prefix_suffix_counts.second) / static_cast<double>(all_suffix_occur);
      double freq_suffix_in_proteins = static_cast<double>(prefix_suffix_counts.second) / static_cast<double>(all_proteins_count);

      if (freq_suffix >= 0.8 && freq_suffix_in_proteins >= 0.3)
      {
        if (prefix_suffix_counts.second != all_suffix_occur)
        {
          OPENMS_LOG_WARN << "More than one decoy suffix observed!" << std::endl;
          OPENMS_LOG_WARN << "Using most frequent decoy suffix (" << (int)(freq_suffix * 100) << "%)" << std::endl;
        }

        return { true, decoy_case_sensitive[case_insensitive_decoy_string], false};
      }
    }

    OPENMS_LOG_ERROR << "Unable to determine decoy string and its position. Please provide a decoy string and its position as parameters." << std::endl;
    return {false, "?", true};
  }

} // namespace OpenMS
