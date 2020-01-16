// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

    defaults_.setValue("best_charge_and_fraction", "false", "Distinguish between fraction and charge states of a peptide. For peptides, abundances will be reported separately for each fraction and charge;\nfor proteins, abundances will be computed based only on the most prevalent charge observed of each peptide (over all fractions).\nBy default, abundances are summed over all charge states.");
    defaults_.setValidStrings("best_charge_and_fraction", true_false);

    defaults_.setValue("consensus:normalize", "false", "Scale peptide abundances so that medians of all samples are equal");
    defaults_.setValidStrings("consensus:normalize", true_false);

    defaults_.setValue("consensus:fix_peptides", "false", "Use the same peptides for protein quantification across all samples.\nWith 'top 0', all peptides that occur in every sample are considered.\nOtherwise ('top N'), the N peptides that occur in the most samples (independently of each other) are selected,\nbreaking ties by total abundance (there is no guarantee that the best co-ocurring peptides are chosen!).");
    defaults_.setValidStrings("consensus:fix_peptides", true_false);

    defaults_.setSectionDescription("consensus", "Additional options for consensus maps (and identification results comprising multiple runs)");

    defaultsToParam_();
  }


  void PeptideAndProteinQuant::countPeptides_(
    vector<PeptideIdentification>& peptides)
  {
    // TODO FRACTION: map ids to fractions
    const int fraction = 1; //TODO: determine from ID?
    for (auto & pep : peptides)
    {
      if (!pep.getHits().empty())
      {
        pep.sort();
        const PeptideHit& hit = pep.getHits()[0]; // get best hit
        PeptideData& data = pep_quant_[hit.getSequence()];
        data.id_count++;
        data.abundances[fraction][hit.getCharge()]; // insert empty element for charge

        // add protein accessions:
        set<String> protein_accessions = hit.extractProteinAccessionsSet();
        data.accessions.insert(protein_accessions.begin(), protein_accessions.end());
      }
    }
  }


  PeptideHit PeptideAndProteinQuant::getAnnotation_(
    vector<PeptideIdentification>& peptides)
  {
    // hits in IDs must already be sorted by score! (done in "countPeptides_")
    if (peptides.empty()) return PeptideHit();

    // get best hit
    const PeptideHit& hit = peptides[0].getHits()[0];

    for (auto pep_it = ++peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      const PeptideHit& current = pep_it->getHits()[0];
      if (current.getSequence() != hit.getSequence())
      {
        // TODO?: warn/error that ambiguous sequences are annotated. check if this can happen
        return PeptideHit();
      }
    }
    return hit;
  }


  void PeptideAndProteinQuant::quantifyFeature_(const FeatureHandle& feature,
                                                const size_t fraction,
                                                const size_t sample,
                                                const PeptideHit& hit)
  {
    // return if annotation for the feature is ambiguous or missing
    if (hit == PeptideHit()) { return; }

    stats_.quant_features++;
    const AASequence& seq = hit.getSequence();
    pep_quant_[seq].abundances[fraction][hit.getCharge()][sample] +=
      feature.getIntensity(); // new map element is initialized with 0
  }


  void PeptideAndProteinQuant::quantifyPeptides(
    const vector<PeptideIdentification>& peptides)
  {
    //////////////////////////////////////////////////////
    // first, use peptide-level results from protein inference:
    // - remove peptides not supported by inference results
    // - update protein accessions according to inference results

    // mapping: peptide seq. (unmodified) -> protein accessions
    // (in protXML, only unmodified peptides are listed)
    map<String, set<String> > pep_info;
    for (auto const & pep : peptides)
    {
      for (auto const & hit : pep.getHits())
      {
        String seq = hit.getSequence().toUnmodifiedString();
        set<String> accessions = hit.extractProteinAccessionsSet();

        // cout << "Sequence: " << seq << " size: " << accessions.size() << " " << *(accessions.begin()) << endl;

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
      if (pep_quant_.empty())
      {
        OPENMS_LOG_ERROR << "No peptides quantified (pep_quant_ is empty)!" << endl;
      }

      PeptideQuant filtered;

      for (auto & pep_q : pep_quant_)  // for all quantified peptides
      {
        String seq = pep_q.first.toUnmodifiedString();
        OPENMS_LOG_DEBUG << "Sequence: " << seq << endl;
        map<String, set<String> >::iterator pos = pep_info.find(seq);
        if (pos != pep_info.end()) // sequence found in protein inference data
        {
          OPENMS_LOG_DEBUG << "Accessions: ";
          for (auto & a : pos->second) {
            OPENMS_LOG_DEBUG << a << "\t";
          }
          OPENMS_LOG_DEBUG << "\n";
          pep_q.second.accessions = pos->second; // replace accessions
          filtered.insert(pep_q);
        }
        else
        {
          OPENMS_LOG_DEBUG << "not found in inference data." << endl;
        }
      }
      pep_quant_ = filtered;
    }

    //////////////////////////////////////////////////////
    // second, perform the actual peptide quantification:
    for (auto & pep_q : pep_quant_)
    {
      if (param_.getValue("best_charge_and_fraction") == "true")
      { // quantify according to the best charge state only:

        // determine which fraction and charge state yields the maximum number of abundances 
        // (break ties by total abundance)
        std::pair<size_t, size_t> best_fraction_and_charge;

        // return false: only identified, not quantified
        if (!getBest_(pep_q.second.abundances, best_fraction_and_charge)) 
        { 
          continue;
        }
        
        // quantify according to the best fraction and charge state only:
        for (auto & sa : pep_q.second.abundances[best_fraction_and_charge.first][best_fraction_and_charge.second])
        {
          pep_q.second.total_abundances[sa.first] = sa.second;
        }
      }
      else
      { // sum up sample abundances over all fractions and charge states:

        for (auto & fa : pep_q.second.abundances)  // for all fractions 
        {
          for (auto & ca : fa.second) // for all charge states
          {  
            for (auto & sa : ca.second) // loop over abundances
            {
              const UInt64 & sample_id = sa.first;
              const double & sample_abundance = sa.second;
              pep_q.second.total_abundances[sample_id] += sample_abundance;
            }
          }
        }
      }

      // count quantified peptide
      if (!pep_q.second.total_abundances.empty()) { stats_.quant_peptides++; }
    }

    //////////////////////////////////////////////////////
    // normalize (optional):
    if ((stats_.n_samples > 1) &&
       (param_.getValue("consensus:normalize") == "true"))
    {      
      normalizePeptides_();
    }
  }

  void PeptideAndProteinQuant::normalizePeptides_()
  {
    /////////////////////////////////////////////////////
    // calculate total peptide abundances 
    // depending on earlier options, these include:
    // - all charges or only the best charge state
    // - all fractions (if multiple fractions are analyzed)
    map<UInt64, DoubleList> abundances; // all peptide abundances by sample
    for (auto & pq : pep_quant_)
    {
      // maybe TODO: treat missing abundance values as zero
      for (auto & sa : pq.second.total_abundances)
      {
        abundances[sa.first].push_back(sa.second);
      }
    }
    if (abundances.size() <= 1) { return; }

    /////////////////////////////////////////////////////
    // compute scale factors on the sample level:
    SampleAbundances medians; // median abundance by sample
    for (auto & ab : abundances)
    {
      medians[ab.first] = Math::median(ab.second.begin(), ab.second.end());
    }

    DoubleList all_medians;
    for (auto & sa : medians)
    {
      all_medians.push_back(sa.second);
    }
    double overall_median = Math::median(all_medians.begin(),
                                         all_medians.end());
    SampleAbundances scale_factors;
    for (auto & med : medians)
    {
      scale_factors[med.first] = overall_median / med.second;
    }

    /////////////////////////////////////////////////////
    // scale all abundance values:
    for (auto & pep_q : pep_quant_)
    {
      // scale total abundances
      for (auto & sta : pep_q.second.total_abundances)
      {
        sta.second *= scale_factors[sta.first];
      }

      // scale individual abundances
      for (auto & fa : pep_q.second.abundances) // for all fractions
      {
        for (auto & ca : fa.second) // for all charge states
        {
          for (auto & sa : ca.second) // loop over abundances
          {
            sa.second *= scale_factors[sa.first];
          }
        }
      }
    }
  }

  String PeptideAndProteinQuant::getAccession_(
    const set<String>& pep_accessions, 
    map<String, String>& accession_to_leader)
  {
    if (accession_to_leader.empty())
    {
      // no info about indistinguishable proteins available
      if (pep_accessions.size() == 1) { return *pep_accessions.begin(); }
    }
    else
    {
      // if all accessions belong to the same group of indistinguishable
      // proteins, return accession of the group leader
      StringList leaders;
      for (auto const & acc : pep_accessions)
      {
        map<String, String>::const_iterator pos = accession_to_leader.find(acc);
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
    OPENMS_LOG_DEBUG << "LEADERS EMPTY: " << endl;
    for (auto const & acc : pep_accessions)
    {
      OPENMS_LOG_DEBUG << acc << endl;
    } 
    return "";
  }


  void PeptideAndProteinQuant::quantifyProteins(const ProteinIdentification&
                                                proteins)
  {
    if (pep_quant_.empty())
    {
      OPENMS_LOG_WARN << "Warning: No peptides quantified." << endl;
    }

    // if information about (indistinguishable) protein groups is available, map
    // each accession to the accession of the leader of its group of proteins:
    map<String, String> accession_to_leader;
    if (!proteins.getIndistinguishableProteins().empty())
    {
      for (auto const & pg : proteins.getIndistinguishableProteins())
      {
        for (auto const & acc : pg.accessions)
        {
          // each accession should only occur once, but we don't check...
          accession_to_leader[acc] = pg.accessions[0];
        }
      }
    }

    // for (auto & a : accession_to_leader) { std::cout << a.first << "\tis led by:\t" << a.second << endl; }

    for (auto const& pep_q : pep_quant_)
    {
      String accession = getAccession_(pep_q.second.accessions,
                                       accession_to_leader);
      OPENMS_LOG_DEBUG << "Peptide id mapped to leader: " << accession << endl;
      if (!accession.empty()) // proteotypic peptide
      {
        prot_quant_[accession].id_count += pep_q.second.id_count;
        for (auto const & sta : pep_q.second.total_abundances)
        {
          // add up contributions of same peptide with different mods:
          String raw_peptide = pep_q.first.toUnmodifiedString();
          prot_quant_[accession].abundances[raw_peptide][sta.first] +=
            sta.second;
        }
      }
    }

    Size top = param_.getValue("top");
    String average = param_.getValue("average");
    bool include_all = param_.getValue("include_all") == "true";
    bool fix_peptides = param_.getValue("consensus:fix_peptides") == "true";

    for (auto & prot_q : prot_quant_)
    {
      if ((top > 0) && (prot_q.second.abundances.size() < top))
      {
        stats_.too_few_peptides++;
        if (!include_all) { continue; } // not enough proteotypic peptides
      }

      vector<String> peptides; // peptides selected for quantification
      if (fix_peptides && (top == 0))
      {
        // consider all peptides that occur in every sample:
        for (auto const & ab : prot_q.second.abundances)
        {
          if (ab.second.size() == stats_.n_samples)
          {
            peptides.push_back(ab.first);
          }
        }
      }
      else if (fix_peptides && (top > 0) &&
               (prot_q.second.abundances.size() > top))
      {
        orderBest_(prot_q.second.abundances, peptides);
        peptides.resize(top);
      }
      else
      {
        // consider all peptides of the protein:
        for (auto const & ab : prot_q.second.abundances)
        {
          peptides.push_back(ab.first);
        }
      }

      map<UInt64, DoubleList> abundances; // all peptide abundances by sample

      // consider only the selected peptides for quantification:
      for (auto & pep : peptides)
      {       
        for (auto & sa : prot_q.second.abundances[pep])
        {
          abundances[sa.first].push_back(sa.second);
        }
      }

      for (auto & ab : abundances)
      {
        // check if the protein has enough peptides in this sample
        if (!include_all && (top > 0) && (ab.second.size() < top))
        {
          continue;
        }

        // if we have more than "top", reduce to the top ones
        if ((top > 0) && (ab.second.size() > top))
        {
          // sort descending:
          sort(ab.second.begin(), ab.second.end(), greater<double>());
          ab.second.resize(top); // remove all but best "top" values
        }

        double result;
        if (average == "median")
        {
          result = Math::median(ab.second.begin(), ab.second.end());
        }
        else if (average == "mean")
        {
          result = Math::mean(ab.second.begin(), ab.second.end());
        }
        else if (average == "weighted_mean")
        {
          double sum_intensities = 0;
          double sum_intensities_squared = 0;
          for (auto const & in : ab.second)
          {
            sum_intensities += in;
            sum_intensities_squared += in * in;
          }
          result = sum_intensities_squared / sum_intensities;
        }
        else // "sum"
        {
          result = Math::sum(ab.second.begin(), ab.second.end());
        }
        prot_q.second.total_abundances[ab.first] = result;
      }

      // update statistics:
      if (prot_q.second.total_abundances.empty()) 
      { 
        stats_.too_few_peptides++; 
      }
      else 
      {
        stats_.quant_proteins++;
      }
    }
  }


  // FRACTIONS: DONE
  void PeptideAndProteinQuant::readQuantData(
    FeatureMap& features,
    const ExperimentalDesign& ed)
  {
    updateMembers_(); // clear data

    stats_.n_samples = ed.getNumberOfSamples();
    stats_.n_fractions = 1;
    stats_.n_ms_files = ed.getNumberOfMSFiles();

    stats_.total_features = features.size();

    for (auto & f : features)
    {
      if (f.getPeptideIdentifications().empty())
      {
        stats_.blank_features++;
        continue;
      }
      countPeptides_(f.getPeptideIdentifications());
      PeptideHit hit = getAnnotation_(f.getPeptideIdentifications());
      FeatureHandle handle(0, f);
      const size_t fraction(1), sample(1);
      quantifyFeature_(handle, fraction, sample, hit); // updates "stats_.quant_features"
    }
    countPeptides_(features.getUnassignedPeptideIdentifications());
    stats_.total_peptides = pep_quant_.size();
    stats_.ambig_features = stats_.total_features - stats_.blank_features -
                            stats_.quant_features;
  }


  void PeptideAndProteinQuant::readQuantData(
    ConsensusMap& consensus, 
    const ExperimentalDesign& ed)
  {
    updateMembers_(); // clear data

    if (consensus.empty())
    {
      OPENMS_LOG_ERROR << "Empty consensus map passed to readQuantData." << endl;
      return;
    }

    stats_.n_samples = ed.getNumberOfSamples();
    stats_.n_fractions = ed.getNumberOfFractions();
    stats_.n_ms_files = ed.getNumberOfMSFiles();

    OPENMS_LOG_DEBUG << "Reading quant data: " << endl;
    OPENMS_LOG_DEBUG << "  MS files        : " << stats_.n_ms_files << endl;
    OPENMS_LOG_DEBUG << "  Fractions       : " << stats_.n_fractions << endl;
    OPENMS_LOG_DEBUG << "  Samples (Assays): " << stats_.n_samples << endl;

    for (auto & c : consensus)
    {
      stats_.total_features += c.getFeatures().size();

      // count features without id
      if (c.getPeptideIdentifications().empty())
      {
        stats_.blank_features += c.getFeatures().size();
        continue;
      }

      countPeptides_(c.getPeptideIdentifications());
      PeptideHit hit = getAnnotation_(c.getPeptideIdentifications());
      for (auto const & f : c.getFeatures())
      {
        // indices in experimental design are 1-based (as in text file)
        // so we need to convert between them
        //TODO MULTIPLEXED: needs to be adapted for multiplexed experiments
        size_t row = f.getMapIndex();
        size_t fraction = ed.getMSFileSection()[row].fraction;
        size_t sample = ed.getMSFileSection()[row].sample;
        quantifyFeature_(f, fraction, sample, hit); // updates "stats_.quant_features"
      }
    }
    countPeptides_(consensus.getUnassignedPeptideIdentifications());
    stats_.total_peptides = pep_quant_.size();
    stats_.ambig_features = stats_.total_features - stats_.blank_features -
                            stats_.quant_features;
  }


  void PeptideAndProteinQuant::readQuantData(
    vector<ProteinIdentification>& proteins,
    vector<PeptideIdentification>& peptides,
    const ExperimentalDesign& ed)
  {
    updateMembers_(); // clear data

    stats_.n_samples = ed.getNumberOfSamples();
    stats_.n_fractions = ed.getNumberOfFractions();
    stats_.n_ms_files = ed.getNumberOfMSFiles();

    stats_.total_features = peptides.size();

    countPeptides_(peptides);

    map<String, String> identifier_to_ms_file;
    for (Size i = 0; i < proteins.size(); ++i)
    {
      StringList ms_files;
      proteins[i].getPrimaryMSRunPath(ms_files);

      if (ms_files.empty()) 
      {
        throw Exception::MissingInformation(
          __FILE__, 
          __LINE__, 
          OPENMS_PRETTY_FUNCTION, 
          "No MS file annotated in protein identification.");
      }
      if (ms_files.size() >= 2) 
      {
        throw Exception::MissingInformation(
          __FILE__, 
          __LINE__, 
          OPENMS_PRETTY_FUNCTION, 
          "More than one ms file annotated in protein identification.");
      }
      identifier_to_ms_file[proteins[i].getIdentifier()] = ms_files[0];
    }

    for (auto & p : peptides)
    {
      if (p.getHits().empty()) { continue; }

      const PeptideHit& hit = p.getHits()[0];
      stats_.quant_features++;
      const AASequence& seq = hit.getSequence();
      const String& ms_file_path = identifier_to_ms_file[p.getIdentifier()];

      // determine sample and fraction by MS file name (stored in protein identification)
      const ExperimentalDesign::MSFileSection& run_section = ed.getMSFileSection();
      auto row = find_if(begin(run_section), end(run_section), 
        [&ms_file_path](const ExperimentalDesign::MSFileSectionEntry& r)
          { 
            return r.path == ms_file_path; 
          });
      size_t sample = row->sample;
      size_t fraction = row->fraction;

      // TODO MULTIPLEXING: think about how id-based quant is done for SILAC, TMT, etc.
      // count peptides in the different fractions, charge states, and samples
      pep_quant_[seq].abundances[fraction][hit.getCharge()][sample] += 1;
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

  // static
  void PeptideAndProteinQuant::annotateQuantificationsToProteins(
    const ProteinQuant& protein_quants, 
    ProteinIdentification& proteins,
    const UInt n_samples)
  {
    auto & id_groups = proteins.getIndistinguishableProteins();

    OPENMS_LOG_DEBUG << "Quantified ind. group with accessions:\n";
    for (const auto& q : protein_quants)
    {
      // accession of quantified protein(group)
      const String & acc = q.first;

      if (q.second.total_abundances.empty()) 
      {
        //TODO maybe just count the number of unquantifiable proteins and report that?
        OPENMS_LOG_DEBUG << "Protein: " << acc << " not quantified." << endl;
        continue; 
      } // not quantified
 
      // lambda to check if a ProteinGroup has accession "acc"
      auto hasProteinInGroup = [&acc] (const ProteinIdentification::ProteinGroup& g)->bool 
      { 
        return find(g.accessions.begin(), g.accessions.end(), acc) != g.accessions.end(); 
      }; 

      // retrieve protein group with accession "acc"
      auto id_group = std::find_if(id_groups.begin(), id_groups.end(), hasProteinInGroup);  

      if (id_group != id_groups.end())
      {
        // copy abundances to float data array
        SampleAbundances total_abundances = q.second.total_abundances;
        // TODO: OPENMS_ASSERT(id_group->float_data_arrays.empty(), "Protein group float data array not empty!.");
        id_group->getFloatDataArrays().resize(1);
        ProteinIdentification::ProteinGroup::FloatDataArray & abundances = id_group->getFloatDataArrays()[0];
        abundances.setName("abundances");        
        abundances.resize(n_samples);

        for (auto const& a : id_group->accessions)
        {
          OPENMS_LOG_DEBUG << a << "\t";
        }
        OPENMS_LOG_DEBUG << "\n";

        for (auto const & s : total_abundances)
        {
          // Note: sample indices are one-based
          abundances[s.first - 1] = s.second;
          OPENMS_LOG_DEBUG << s.second << "\t";
        }
        OPENMS_LOG_DEBUG << endl;
      }
      else
      {
        throw Exception::MissingInformation(
          __FILE__, 
          __LINE__, 
          OPENMS_PRETTY_FUNCTION, 
          "Protein group quantified that is not present in inference data.");
      } 
    }

   // remove all protein groups that have not been quantified
   auto notQuantified = [] (const ProteinIdentification::ProteinGroup& g)->bool { return g.getFloatDataArrays().empty(); }; 
   id_groups.erase(
     remove_if(id_groups.begin(), id_groups.end(), notQuantified), 
     id_groups.end());
  } 

}
