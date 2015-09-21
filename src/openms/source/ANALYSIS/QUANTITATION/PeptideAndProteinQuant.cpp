// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <algorithm> // for "equal"

using namespace std;

namespace OpenMS
{

  PeptideAndProteinQuant::PeptideAndProteinQuant() :
    DefaultParamHandler("PeptideAndProteinQuant"), stats_(), pep_quant_(),
    prot_quant_()
  {
    defaults_.setValue("top", 3, "Calculate protein abundance from this number of proteotypic peptides (most abundant first; '0' for all)");
    defaults_.setMinInt("top", 0);

    defaults_.setValue("average", "median", "Averaging method used to compute protein abundances from peptide abundances");
    defaults_.setValidStrings("average", ListUtils::create<String>("median,mean,weighted_mean,sum"));

    StringList true_false = ListUtils::create<String>("true,false");

    defaults_.setValue("include_all", "false", "Include results for proteins with fewer proteotypic peptides than indicated by 'top' (no effect if 'top' is 0 or 1)");
    defaults_.setValidStrings("include_all", true_false);

    defaults_.setValue("filter_charge", "false", "Distinguish between charge states of a peptide. For peptides, abundances will be reported separately for each charge;\nfor proteins, abundances will be computed based only on the most prevalent charge of each peptide.\nBy default, abundances are summed over all charge states.");
    defaults_.setValidStrings("filter_charge", true_false);

    defaults_.setValue("consensus:normalize", "false", "Scale peptide abundances so that medians of all samples are equal");
    defaults_.setValidStrings("consensus:normalize", true_false);

    defaults_.setValue("consensus:fix_peptides", "false", "Use the same peptides for protein quantification across all samples.\nWith 'top 0', all peptides that occur in every sample are considered.\nOtherwise ('top N'), the N peptides that occur in the most samples (independently of each other) are selected,\nbreaking ties by total abundance (there is no guarantee that the best co-ocurring peptides are chosen!).");
    defaults_.setValidStrings("consensus:fix_peptides", true_false);

    defaults_.setSectionDescription("consensus", "Additional options for consensus maps (and identification results comprising multiple runs)");

    defaultsToParam_();
  }


  void PeptideAndProteinQuant::countPeptides_(vector<PeptideIdentification>&
                                              peptides)
  {
    for (vector<PeptideIdentification>::iterator pep_it =
           peptides.begin(); pep_it != peptides.end(); ++pep_it)
    {
      if (!pep_it->getHits().empty())
      {
        pep_it->sort();
        const PeptideHit& hit = pep_it->getHits()[0];
        PeptideData& data = pep_quant_[hit.getSequence()];
        data.id_count++;
        data.abundances[hit.getCharge()]; // insert empty element for charge
        // add protein accessions:
        set<String> protein_accessions = hit.extractProteinAccessions();
        data.accessions.insert(protein_accessions.begin(), protein_accessions.end());
      }
    }
  }


  PeptideHit PeptideAndProteinQuant::getAnnotation_(
    vector<PeptideIdentification>& peptides)
  {
    // hits in IDs must already be sorted by score! (done in "countPeptides_")
    if (peptides.empty()) return PeptideHit();

    const PeptideHit& hit = peptides[0].getHits()[0];
    for (vector<PeptideIdentification>::iterator pep_it = ++peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      const PeptideHit& current = pep_it->getHits()[0];
      if (current.getSequence() != hit.getSequence())
      {
        return PeptideHit();
      }
    }
    return hit;
  }


  void PeptideAndProteinQuant::quantifyFeature_(const FeatureHandle& feature,
                                                const PeptideHit& hit)
  {
    if (hit == PeptideHit())
    {
      return; // annotation for the feature is ambiguous or missing
    }
    stats_.quant_features++;
    const AASequence& seq = hit.getSequence();
    pep_quant_[seq].abundances[hit.getCharge()][feature.getMapIndex()] +=
      feature.getIntensity(); // new map element is initialized with 0
  }


  void PeptideAndProteinQuant::quantifyPeptides(
    const vector<PeptideIdentification>& peptides)
  {
    // first, use peptide-level results from protein inference:
    // - remove peptides not supported by inference results
    // - update protein accessions according to inference results

    // mapping: peptide seq. (unmodified) -> protein accessions
    // (in protXML, only unmodified peptides are listed)
    map<String, set<String> > pep_info;
    for (vector<PeptideIdentification>::const_iterator pep_it = 
           peptides.begin(); pep_it != peptides.end(); ++pep_it)
    {
      for (vector<PeptideHit>::const_iterator hit_it =
             pep_it->getHits().begin(); hit_it != pep_it->getHits().end();
           ++hit_it)
      {
        String seq = hit_it->getSequence().toUnmodifiedString();
        set<String> accessions = hit_it->extractProteinAccessions();
        // If a peptide is seen multiple times, the protein accessions should
        // always be the same, so only the first time it should be necessary to
        // insert them. However, just in case there a differences in the
        // accessions, we accumulate them all (probably unnecessary work):
        pep_info[seq].insert(accessions.begin(), accessions.end());
      }
    }
    // if inference results are given, filter quant. data accordingly:
    if (!pep_info.empty())
    {
      PeptideQuant filtered;
      for (PeptideQuant::iterator q_it = pep_quant_.begin(); 
           q_it != pep_quant_.end(); ++q_it)
      {
        String seq = q_it->first.toUnmodifiedString();
        map<String, set<String> >::iterator pos = pep_info.find(seq);
        if (pos != pep_info.end()) // sequence found in protein inference data
        {
          q_it->second.accessions = pos->second; // replace accessions
          filtered.insert(*q_it);
        }
      }
      pep_quant_ = filtered;
    }

    // now perform the actual peptide quantification:
    for (PeptideQuant::iterator q_it = pep_quant_.begin();
         q_it != pep_quant_.end(); ++q_it)
    {
      if (param_.getValue("filter_charge") == "true")
      {
        // find charge state with abundances for highest number of samples
        // (break ties by total abundance):
        IntList charges; // sorted charge states (best first)
        orderBest_(q_it->second.abundances, charges);
        if (charges.empty()) continue; // only identified, not quantified
        Int best_charge = charges[0];

        // quantify according to the best charge state only:
        for (SampleAbundances::iterator samp_it =
               q_it->second.abundances[best_charge].begin(); samp_it !=
             q_it->second.abundances[best_charge].end(); ++samp_it)
        {
          q_it->second.total_abundances[samp_it->first] = samp_it->second;
        }
      }
      else
      {
        // sum up abundances over all charge states:
        for (map<Int, SampleAbundances>::iterator ab_it =
               q_it->second.abundances.begin(); ab_it !=
             q_it->second.abundances.end(); ++ab_it)
        {
          for (SampleAbundances::iterator samp_it = ab_it->second.begin();
               samp_it != ab_it->second.end(); ++samp_it)
          {
            q_it->second.total_abundances[samp_it->first] += samp_it->second;
          }
        }
      }
      if (!q_it->second.total_abundances.empty())
        stats_.quant_peptides++;
    }

    if ((stats_.n_samples > 1) &&
        (param_.getValue("consensus:normalize") == "true"))
    {
      normalizePeptides_();
    }
  }


  void PeptideAndProteinQuant::normalizePeptides_()
  {
    // gather data:
    map<UInt64, DoubleList> abundances; // all peptide abundances by sample
    for (PeptideQuant::iterator q_it = pep_quant_.begin();
         q_it != pep_quant_.end(); ++q_it)
    {
      // maybe TODO: treat missing abundance values as zero
      for (SampleAbundances::iterator samp_it =
             q_it->second.total_abundances.begin(); samp_it !=
           q_it->second.total_abundances.end(); ++samp_it)
      {
        abundances[samp_it->first].push_back(samp_it->second);
      }
    }
    if (abundances.size() <= 1) return;

    // compute scale factors for all samples:
    SampleAbundances medians; // median abundance by sample
    for (map<UInt64, DoubleList>::iterator ab_it = abundances.begin();
         ab_it != abundances.end(); ++ab_it)
    {
      medians[ab_it->first] = Math::median(ab_it->second.begin(),
                                           ab_it->second.end());
    }
    DoubleList all_medians;
    for (SampleAbundances::iterator med_it = medians.begin();
         med_it != medians.end(); ++med_it)
    {
      all_medians.push_back(med_it->second);
    }
    double overall_median = Math::median(all_medians.begin(),
                                         all_medians.end());
    SampleAbundances scale_factors;
    for (SampleAbundances::iterator med_it = medians.begin();
         med_it != medians.end(); ++med_it)
    {
      scale_factors[med_it->first] = overall_median / med_it->second;
    }

    // scale all abundance values:
    for (PeptideQuant::iterator q_it = pep_quant_.begin();
         q_it != pep_quant_.end(); ++q_it)
    {
      for (SampleAbundances::iterator tot_it =
             q_it->second.total_abundances.begin(); tot_it !=
           q_it->second.total_abundances.end(); ++tot_it)
      {
        tot_it->second *= scale_factors[tot_it->first];
      }
      for (map<Int, SampleAbundances>::iterator ab_it =
             q_it->second.abundances.begin(); ab_it !=
           q_it->second.abundances.end(); ++ab_it)
      {
        for (SampleAbundances::iterator samp_it = ab_it->second.begin();
             samp_it != ab_it->second.end(); ++samp_it)
        {
          samp_it->second *= scale_factors[samp_it->first];
        }
      }
    }
  }


  String PeptideAndProteinQuant::getAccession_(
    const set<String>& pep_accessions, map<String, String>& accession_to_leader)
  {
    if (accession_to_leader.empty())
    {
      // no info about indistinguishable proteins available
      if (pep_accessions.size() == 1) return *pep_accessions.begin();
    }
    else
    {
      // if all accessions belong to the same group of indistinguishable
      // proteins, return accession of the group leader
      StringList leaders;
      for (set<String>::const_iterator it = pep_accessions.begin();
           it != pep_accessions.end(); ++it)
      {
        map<String, String>::const_iterator pos = accession_to_leader.find(*it);
        if (pos != accession_to_leader.end()) leaders.push_back(pos->second);
        // if the protein accession was not found, this is not an error:
        // if there's not enough evidence for a protein, it won't occur in
        // the protXML - so we won't quantify it
      }
      if (leaders.empty()) return "";

      bool all_equal = equal(leaders.begin(), --leaders.end(),
                             ++leaders.begin());
      if (all_equal) return leaders[0];
    }
    return "";
  }


  void PeptideAndProteinQuant::quantifyProteins(const ProteinIdentification&
                                                proteins)
  {
    // if information about (indistinguishable) protein groups is available, map
    // each accession to the accession of the leader of its group of proteins:
    map<String, String> accession_to_leader;
    if (!proteins.getIndistinguishableProteins().empty())
    {
      for (vector<ProteinIdentification::ProteinGroup>::const_iterator pg_it =
             proteins.getIndistinguishableProteins().begin(); pg_it !=
             proteins.getIndistinguishableProteins().end(); ++pg_it)
      {
        for (StringList::const_iterator acc_it = pg_it->accessions.begin();
             acc_it != pg_it->accessions.end(); ++acc_it)
        {
          // each accession should only occur once, but we don't check...
          accession_to_leader[*acc_it] = pg_it->accessions[0];
        }
      }
    }

    for (PeptideQuant::const_iterator pep_it = pep_quant_.begin();
         pep_it != pep_quant_.end(); ++pep_it)
    {
      String accession = getAccession_(pep_it->second.accessions,
                                       accession_to_leader);
      if (!accession.empty()) // proteotypic peptide
      {
        prot_quant_[accession].id_count += pep_it->second.id_count;
        for (SampleAbundances::const_iterator tot_it =
               pep_it->second.total_abundances.begin(); tot_it !=
             pep_it->second.total_abundances.end(); ++tot_it)
        {
          // add up contributions of same peptide with different mods:
          String raw_peptide = pep_it->first.toUnmodifiedString();
          prot_quant_[accession].abundances[raw_peptide][tot_it->first] +=
            tot_it->second;
        }
      }
    }

    Size top = param_.getValue("top");
    String average = param_.getValue("average");
    bool include_all = param_.getValue("include_all") == "true";
    bool fix_peptides = param_.getValue("consensus:fix_peptides") == "true";

    for (ProteinQuant::iterator prot_it = prot_quant_.begin();
         prot_it != prot_quant_.end(); ++prot_it)
    {
      if ((top > 0) && (prot_it->second.abundances.size() < top))
      {
        stats_.too_few_peptides++;
        if (!include_all)
          continue; // not enough proteotypic peptides
      }

      vector<String> peptides; // peptides selected for quantification
      if (fix_peptides && (top == 0))
      {
        // consider all peptides that occur in every sample:
        for (map<String, SampleAbundances>::iterator ab_it =
               prot_it->second.abundances.begin(); ab_it !=
             prot_it->second.abundances.end(); ++ab_it)
        {
          if (ab_it->second.size() == stats_.n_samples)
          {
            peptides.push_back(ab_it->first);
          }
        }
      }
      else if (fix_peptides && (top > 0) &&
               (prot_it->second.abundances.size() > top))
      {
        orderBest_(prot_it->second.abundances, peptides);
        peptides.resize(top);
      }
      else
      {
        // consider all peptides:
        for (map<String, SampleAbundances>::iterator ab_it =
               prot_it->second.abundances.begin(); ab_it !=
             prot_it->second.abundances.end(); ++ab_it)
        {
          peptides.push_back(ab_it->first);
        }
      }

      map<UInt64, DoubleList> abundances; // all peptide abundances by sample
      // consider only the peptides selected above for quantification:
      for (vector<String>::iterator pep_it = peptides.begin();
           pep_it != peptides.end(); ++pep_it)
      {
        SampleAbundances& current_ab = prot_it->second.abundances[*pep_it];
        for (SampleAbundances::iterator samp_it = current_ab.begin();
             samp_it != current_ab.end(); ++samp_it)
        {
          abundances[samp_it->first].push_back(samp_it->second);
        }
      }

      for (map<UInt64, DoubleList>::iterator ab_it = abundances.begin();
           ab_it != abundances.end(); ++ab_it)
      {
        if (!include_all && (top > 0) && (ab_it->second.size() < top))
        {
          continue; // not enough peptide abundances for this sample
        }
        if ((top > 0) && (ab_it->second.size() > top))
        {
          // sort descending:
          sort(ab_it->second.begin(), ab_it->second.end(), greater<double>());
          ab_it->second.resize(top); // remove all but best "top" values
        }

        double result;
        if (average == "median")
        {
          result = Math::median(ab_it->second.begin(), ab_it->second.end());
        }
        else if (average == "mean")
        {
          result = Math::mean(ab_it->second.begin(), ab_it->second.end());
        }
        else if (average == "weighted_mean")
        {
          double sum_intensities = 0;
          double sum_intensities_squared = 0;
          for (DoubleList::const_iterator it_intensities =
                 ab_it->second.begin(); it_intensities != ab_it->second.end();
               ++it_intensities)
          {
            sum_intensities += (*it_intensities);
            sum_intensities_squared += (*it_intensities) * (*it_intensities);
          }
          result = sum_intensities_squared / sum_intensities;
        }
        else // "sum"
        {
          result = Math::sum(ab_it->second.begin(), ab_it->second.end());
        }
        prot_it->second.total_abundances[ab_it->first] = result;
      }

      // update statistics:
      if (prot_it->second.total_abundances.empty()) stats_.too_few_peptides++;
      else stats_.quant_proteins++;
    }
  }


  void PeptideAndProteinQuant::readQuantData(FeatureMap& features)
  {
    updateMembers_(); // clear data
    stats_.n_samples = 1;
    stats_.total_features = features.size();

    for (FeatureMap::Iterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it)
    {
      if (feat_it->getPeptideIdentifications().empty())
      {
        stats_.blank_features++;
        continue;
      }
      countPeptides_(feat_it->getPeptideIdentifications());
      PeptideHit hit = getAnnotation_(feat_it->getPeptideIdentifications());
      FeatureHandle handle(0, *feat_it);
      quantifyFeature_(handle, hit); // updates "stats_.quant_features"
    }
    countPeptides_(features.getUnassignedPeptideIdentifications());
    stats_.total_peptides = pep_quant_.size();
    stats_.ambig_features = stats_.total_features - stats_.blank_features -
                            stats_.quant_features;
  }


  void PeptideAndProteinQuant::readQuantData(ConsensusMap& consensus)
  {
    updateMembers_(); // clear data
    stats_.n_samples = consensus.getFileDescriptions().size();

    for (ConsensusMap::Iterator cons_it = consensus.begin();
         cons_it != consensus.end(); ++cons_it)
    {
      stats_.total_features += cons_it->getFeatures().size();
      if (cons_it->getPeptideIdentifications().empty())
      {
        stats_.blank_features += cons_it->getFeatures().size();
        continue;
      }
      countPeptides_(cons_it->getPeptideIdentifications());
      PeptideHit hit = getAnnotation_(cons_it->getPeptideIdentifications());
      for (ConsensusFeature::HandleSetType::const_iterator feat_it =
             cons_it->getFeatures().begin(); feat_it !=
           cons_it->getFeatures().end(); ++feat_it)
      {
        quantifyFeature_(*feat_it, hit); // updates "stats_.quant_features"
      }
    }
    countPeptides_(consensus.getUnassignedPeptideIdentifications());
    stats_.total_peptides = pep_quant_.size();
    stats_.ambig_features = stats_.total_features - stats_.blank_features -
                            stats_.quant_features;
  }


  void PeptideAndProteinQuant::readQuantData(
    vector<ProteinIdentification>& proteins,
    vector<PeptideIdentification>& peptides)
  {
    updateMembers_(); // clear data
    stats_.n_samples = proteins.size();
    stats_.total_features = peptides.size();

    countPeptides_(peptides);

    // treat identification runs as different samples - otherwise we could just
    // use the "id_count" element of PeptideData (filled by "countPeptides_") to
    // get the total spectral count for each peptide:
    map<String, Size> identifiers;
    for (Size i = 0; i < proteins.size(); ++i)
    {
      identifiers[proteins[i].getIdentifier()] = i;
    }

    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      if (pep_it->getHits().empty()) continue;
      const PeptideHit& hit = pep_it->getHits()[0];
      stats_.quant_features++;
      const AASequence& seq = hit.getSequence();
      Size sample = identifiers[pep_it->getIdentifier()];
      pep_quant_[seq].abundances[hit.getCharge()][sample] += 1;
    }
    stats_.total_peptides = pep_quant_.size();
  }


  void PeptideAndProteinQuant::updateMembers_()
  {
    // reset everything:
    stats_ = Statistics();
    pep_quant_.clear();
    prot_quant_.clear();
  }


  const PeptideAndProteinQuant::Statistics&
  PeptideAndProteinQuant::getStatistics()
  {
    return stats_;
  }


  const PeptideAndProteinQuant::PeptideQuant&
  PeptideAndProteinQuant::getPeptideResults()
  {
    return pep_quant_;
  }


  const PeptideAndProteinQuant::ProteinQuant&
  PeptideAndProteinQuant::getProteinResults()
  {
    return prot_quant_;
  }

}
