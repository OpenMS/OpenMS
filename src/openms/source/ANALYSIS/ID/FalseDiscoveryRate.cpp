// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/LogStream.h>

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
    defaults_.setValidStrings("no_qvalues", ListUtils::create<String>("true,false"));
    defaults_.setValue("use_all_hits", "false", "If 'true' not only the first hit, but all are used (peptides only)");
    defaults_.setValidStrings("use_all_hits", ListUtils::create<String>("true,false"));
    defaults_.setValue("split_charge_variants", "false", "If 'true' charge variants are treated separately (for peptides of combined target/decoy searches only).");
    defaults_.setValidStrings("split_charge_variants", ListUtils::create<String>("true,false"));
    defaults_.setValue("treat_runs_separately", "false", "If 'true' different search runs are treated separately (for peptides of combined target/decoy searches only).");
    defaults_.setValidStrings("treat_runs_separately", ListUtils::create<String>("true,false"));
    defaults_.setValue("add_decoy_peptides", "false", "If 'true' decoy peptides will be written to output file, too. The q-value is set to the closest target score.");
    defaults_.setValidStrings("add_decoy_peptides", ListUtils::create<String>("true,false"));
    defaultsToParam_();
  }

  void FalseDiscoveryRate::apply(vector<PeptideIdentification>& ids)
  {
    bool q_value = !param_.getValue("no_qvalues").toBool();
    bool use_all_hits = param_.getValue("use_all_hits").toBool();
    bool treat_runs_separately = param_.getValue("treat_runs_separately").toBool();
    bool split_charge_variants = param_.getValue("split_charge_variants").toBool();
    bool add_decoy_peptides = param_.getValue("add_decoy_peptides").toBool();
#ifdef FALSE_DISCOVERY_RATE_DEBUG
    cerr << "Parameters: no_qvalues=" << !q_value << ", use_all_hits=" << use_all_hits << ", treat_runs_separately=" << treat_runs_separately << ", split_charge_variants=" << split_charge_variants << endl;
#endif


    if (ids.empty())
    {
      LOG_WARN << "No peptide identifications given to FalseDiscoveryRate! No calculation performed.\n";
      return;
    }

    // first search for all identifiers and charge variants
    set<String> identifiers;
    set<SignedSize> charge_variants;
    for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
    {
      identifiers.insert(it->getIdentifier());
      it->sort();

      if (!use_all_hits)
      {
        it->getHits().resize(1);
      }

      for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        charge_variants.insert(pit->getCharge());
      }
    }

#ifdef FALSE_DISCOVERY_RATE_DEBUG
    cerr << "#id-runs: " << identifiers.size() << " ";
    for (set<String>::const_iterator it = identifiers.begin(); it != identifiers.end(); ++it)
    {
      cerr << "," << *it;
    }
    cerr << endl;


    cerr << "#of charge states: " << charge_variants.size() << " ";
    for (set<SignedSize>::const_iterator it = charge_variants.begin(); it != charge_variants.end(); ++it)
    {
      cerr << "," << *it;
    }
    cerr << endl;
#endif

    for (set<SignedSize>::const_iterator zit = charge_variants.begin(); zit != charge_variants.end(); ++zit)
    {
#ifdef FALSE_DISCOVERY_RATE_DEBUG
      cerr << "Charge variant=" << *zit << endl;
#endif

      // for all identifiers
      for (set<String>::const_iterator iit = identifiers.begin(); iit != identifiers.end(); ++iit)
      {
        if (!treat_runs_separately && iit != identifiers.begin())
        {
          continue;
        }

#ifdef FALSE_DISCOVERY_RATE_DEBUG
        cerr << "Id-run: " << *iit << endl;
#endif
        // get the scores of all peptide hits
        vector<double> target_scores, decoy_scores;
        for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
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
              LOG_FATAL_ERROR << "Meta value 'target_decoy' does not exists, reindex the idXML file with 'PeptideIndexer' first (run-id='" << it->getIdentifier() << ", rank=" << i + 1 << " of " << it->getHits().size() << ")!" << endl;
              throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Meta value 'target_decoy' does not exist!");
            }

            String target_decoy(it->getHits()[i].getMetaValue("target_decoy"));
            if (target_decoy == "target" || target_decoy == "target+decoy")
            {
              target_scores.push_back(it->getHits()[i].getScore());
            }
            else
            {
              if (target_decoy == "decoy")
              {
                decoy_scores.push_back(it->getHits()[i].getScore());
              }
              else
              {
                if (target_decoy != "")
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
          LOG_ERROR << error_string << std::endl;
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
          LOG_ERROR << error_string << std::endl;
        }

        if (target_scores.empty() || decoy_scores.empty())
        {
          // no remove the the relevant entries, or put 'pseudo-scores' in
          for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
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
                LOG_FATAL_ERROR << "Meta value 'target_decoy' does not exists, reindex the idXML file with 'PeptideIndexer' (run-id='" << it->getIdentifier() << ", rank=" << i + 1 << " of " << hits.size() << ")!" << endl;
                throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Meta value 'target_decoy' does not exist!");
              }

              String target_decoy(hits[i].getMetaValue("target_decoy"));
              if (target_decoy == "target" || target_decoy == "target+decoy")
              {
                // if it is a target hit, there are now decoys, fdr/q-value should be zero then
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
        bool higher_score_better(ids.begin()->isHigherScoreBetter());
        Map<double, double> score_to_fdr;
        calculateFDRs_(score_to_fdr, target_scores, decoy_scores, q_value, higher_score_better);

        // annotate fdr
        for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
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

  void FalseDiscoveryRate::apply(vector<PeptideIdentification>& fwd_ids, vector<PeptideIdentification>& rev_ids)
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
    Map<double, double> score_to_fdr;
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

  void FalseDiscoveryRate::apply(vector<ProteinIdentification>& ids)
  {
    if (ids.empty())
    {
      LOG_WARN << "No protein identifications given to FalseDiscoveryRate! No calculation performed.\n";
      return;
    }

    vector<double> target_scores, decoy_scores;
    for (vector<ProteinIdentification>::const_iterator it = ids.begin(); it != ids.end(); ++it)
    {
      for (vector<ProteinHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        if (!pit->metaValueExists("target_decoy"))
        {
          LOG_FATAL_ERROR << "Meta value 'target_decoy' does not exists, reindex the idXML file with 'PeptideIndexer' (run-id='" << it->getIdentifier() << ", accession=" << pit->getAccession() << ")!" << endl;
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

    bool q_value = !param_.getValue("no_qvalues").toBool();
    bool higher_score_better = ids.begin()->isHigherScoreBetter();

    // calculate fdr for the forward scores
    Map<double, double> score_to_fdr;
    calculateFDRs_(score_to_fdr, target_scores, decoy_scores, q_value, higher_score_better);

    // annotate fdr
    String score_type = ids.begin()->getScoreType() + "_score";
    for (vector<ProteinIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
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

    return;
  }

  void FalseDiscoveryRate::apply(vector<ProteinIdentification>& fwd_ids, vector<ProteinIdentification>& rev_ids)
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
    Map<double, double> score_to_fdr;
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

    return;
  }

  void FalseDiscoveryRate::calculateFDRs_(Map<double, double>& score_to_fdr, vector<double>& target_scores, vector<double>& decoy_scores, bool q_value, bool higher_score_better)
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
      if (k == 0) { score_to_fdr[ds] = score_to_fdr[target_scores[0]]; continue; }

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

} // namespace OpenMS
