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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>
#include <map>
#include <unordered_map>

namespace OpenMS
{
  BasicProteinInferenceAlgorithm::BasicProteinInferenceAlgorithm():
      DefaultParamHandler("BasicProteinInferenceAlgorithm"),
      ProgressLogger()
  {
    //registerIntOption_("min_peptides_per_protein", "<num>", 2, "Minimal number of peptides needed for a protein identification", false);
    //setMinInt_("min_peptides_per_protein", 1);
    defaults_.setValue("score_aggregation_method",
                       "maximum",
                       "How to aggregate scores of PSM matching to the same protein?");
    defaults_.setValidStrings("score_aggregation_method", ListUtils::create<String>("maximum,product,sum"));
    defaults_.setValue("treat_charge_variants_separately", "true",
                       "If this is true, different charge variants of the same peptide sequence count as inidividual evidences.");
    defaults_.setValue("treat_modification_variants_separately", "true",
                       "If this is true, different modification variants of the same peptide sequence count as individual evidences.");
    defaults_.setValue("use_shared_peptides", "true", "If this is true, shared peptides are used as evidences.");
    defaults_.setValue("skip_count_annotation", "false", "If this is true, peptide counts won't be annotated at the proteins.");
    defaultsToParam_();
  }

  void BasicProteinInferenceAlgorithm::run(std::vector<PeptideIdentification> &pep_ids,
                                           std::vector<ProteinIdentification> &prot_ids) const
  {
    //Size min_peptides_per_protein = getIntOption_("min_peptides_per_protein");
    bool treat_charge_variants_separately(param_.getValue("treat_charge_variants_separately").toBool());
    bool treat_modification_variants_separately(param_.getValue("treat_modification_variants_separately").toBool());
    bool use_shared_peptides(param_.getValue("use_shared_peptides").toBool());
    bool skip_count_annotation(param_.getValue("skip_count_annotation").toBool());

    String aggMethodString(param_.getValue("score_aggregation_method").toString());
    AggregationMethod aggregation_method = AggregationMethod::MAXIMUM;

    if (aggMethodString == "maximum")
    {
      aggregation_method = AggregationMethod::MAXIMUM;
    }
    else if (aggMethodString == "product")
    {
      aggregation_method = AggregationMethod::PROD;
    }
    else if (aggMethodString == "sum")
    {
      aggregation_method = AggregationMethod::SUM;
    }

    std::unordered_map<std::string, std::map<Int, PeptideHit*>> best_pep{};

    // iterate over runs
    std::unordered_map<std::string, std::pair<ProteinHit *, Size>> acc_to_protein_hitP_and_count;

    for (auto &prot_run : prot_ids)
    {
      acc_to_protein_hitP_and_count.clear();
      best_pep.clear();

      ProteinIdentification::SearchParameters sp = prot_run.getSearchParameters();
      prot_run.setSearchEngine("TOPPProteinInference_" + aggMethodString);
      sp.setMetaValue("use_shared_peptides", use_shared_peptides);
      sp.setMetaValue("treat_charge_variants_separately", treat_charge_variants_separately);
      sp.setMetaValue("treat_modification_variants_separately", treat_modification_variants_separately);
      prot_run.setSearchParameters(sp);

      //create Accession to ProteinHit and peptide count map. To have quick access later.
      for (auto &phit : prot_run.getHits())
      {
        acc_to_protein_hitP_and_count[phit.getAccession()] = std::make_pair<ProteinHit *, Size>(&phit, 0);
      }

      for (auto &pep : pep_ids)
      {
        //skip if it does not belong to run
        if (pep.getIdentifier() != prot_run.getIdentifier())
          continue;
        //make sure that first = best hit
        pep.sort();

        //TODO think about if using any but the best PSM makes sense in such a simple aggregation scheme
        //for (auto& hit : pep.getHits())
        //{
        PeptideHit &hit = pep.getHits()[0];
        //skip if shared and option not enabled
        if (!use_shared_peptides &&
            (!hit.metaValueExists("protein_references") || (hit.getMetaValue("protein_references") == "non-unique")))
          continue;

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
        {
          auto &new_entry = best_pep[lookup_seq];
          new_entry = std::map<Int, PeptideHit *>();
          new_entry[lookup_charge] = &hit;
        }
        else
        {
          auto current_best_pep_charge_it = current_best_pep_it->second.find(lookup_charge);
          if (current_best_pep_charge_it == current_best_pep_it->second.end())
          {
            current_best_pep_it->second[lookup_charge] = &hit;
          }
          else if (
              (prot_run.isHigherScoreBetter() && (hit.getScore() > current_best_pep_charge_it->second->getScore())) ||
              (!prot_run.isHigherScoreBetter() && (hit.getScore() < current_best_pep_charge_it->second->getScore())))
          {
            current_best_pep_charge_it->second = &hit;
          }
        }
        //}
      }

      // update protein scores
      for (const auto &charge_to_pep_hit_map : best_pep)
      {
        // The next line assumes that PeptideHits of different charge states necessarily share the same
        // protein accessions
        // TODO this could be done for mods, too (first hashing AASeq, then the mods)
        for (const auto &acc : charge_to_pep_hit_map.second.begin()->second->extractProteinAccessionsSet())
        {
          for (const auto &pep_hit : charge_to_pep_hit_map.second)
          {
            auto& prot_count_pair = acc_to_protein_hitP_and_count[acc];
            ProteinHit *protein = prot_count_pair.first;
            prot_count_pair.second++;

            double new_score = pep_hit.second->getScore();

            // Note: This requires/works only with Posterior (Error) Probabilities
            if (!prot_run.isHigherScoreBetter())
              new_score = 1. - new_score;
            switch (aggregation_method)
            {
              //TODO for 0 probability peptides we could also multiply a minimum value
              case AggregationMethod::PROD :
                if (new_score > 0.0)
                  protein->setScore(protein->getScore() * new_score);
                break;
              case AggregationMethod::SUM :
                protein->setScore(protein->getScore() + new_score);
                break;
              case AggregationMethod::MAXIMUM :
                protein->setScore(std::fmax(double(protein->getScore()), new_score));
                break;
            }
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
      //TODO Allow count as aggregation method -> i.e. set as protein score?
    }

    //TODO Filtering? I think this should be done separate afterwards with IDFilter
    //for all protein hits for the id run, only accept proteins that have at least 'min_peptides_per_protein' peptides
  }
} //namespace OpenMS
