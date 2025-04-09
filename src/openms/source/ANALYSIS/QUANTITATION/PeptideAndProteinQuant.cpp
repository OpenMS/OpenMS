// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>

using namespace std;

namespace OpenMS
{

  PeptideAndProteinQuant::PeptideAndProteinQuant() :
    DefaultParamHandler("PeptideAndProteinQuant"), stats_(), pep_quant_(),
    prot_quant_()
  {
    std::vector<std::string> true_false = {"true","false"};

    defaults_.setValue("method", "top", "- top - quantify based on three most abundant peptides (number can be changed in 'top').\n- iBAQ (intensity based absolute quantification), calculate the sum of all peptide peak intensities divided by the number of theoretically observable tryptic peptides (https://rdcu.be/cND1J). Warning: only consensusXML or featureXML input is allowed!");
    defaults_.setValidStrings("method", {"top","iBAQ"});

    defaults_.setValue("top:N", 3, "Calculate protein abundance from this number of proteotypic peptides (most abundant first; '0' for all)");
    defaults_.setMinInt("top:N", 0);

    defaults_.setValue("top:aggregate", "median", "Aggregation method used to compute protein abundances from peptide abundances");
    defaults_.setValidStrings("top:aggregate", {"median","mean","weighted_mean","sum"});

    defaults_.setValue("top:include_all", "false", "Include results for proteins with fewer proteotypic peptides than indicated by 'N' (no effect if 'N' is 0 or 1)");
    defaults_.setValidStrings("top:include_all", true_false);

    defaults_.setSectionDescription("top", "Additional options for custom quantification using top N peptides.");


    defaults_.setValue("best_charge_and_fraction", "false", "Distinguish between fraction and charge states of a peptide. For peptides, abundances will be reported separately for each fraction and charge;\nfor proteins, abundances will be computed based only on the most prevalent charge observed of each peptide (over all fractions).\nBy default, abundances are summed over all charge states.");
    defaults_.setValidStrings("best_charge_and_fraction", true_false);

    defaults_.setValue("consensus:normalize", "false", "Scale peptide abundances so that medians of all samples are equal");
    defaults_.setValidStrings("consensus:normalize", true_false);

    defaults_.setValue("consensus:fix_peptides", "false", "Use the same peptides for protein quantification across all samples.\nWith 'N 0',"
     "all peptides that occur in every sample are considered.\nOtherwise ('N'), the N peptides that occur in the most samples (independently of each other) are selected,\nbreaking ties by total abundance (there is no guarantee that the best co-ocurring peptides are chosen!).");
    defaults_.setValidStrings("consensus:fix_peptides", true_false);

    defaults_.setSectionDescription("consensus", "Additional options for consensus maps (and identification results comprising multiple runs)");

    defaultsToParam_();
  }

  // doesn't only count but also some initialization TODO: rename
  void PeptideAndProteinQuant::countPeptides_(
    vector<PeptideIdentification>& peptides)
  {
    for (auto & pep : peptides)
    {
      if (pep.getHits().empty()) continue;
      pep.sort(); // TODO: move this out of count peptides
      const PeptideHit& hit = pep.getHits()[0]; // get best hit
      PeptideData& data = pep_quant_[hit.getSequence()];
      data.psm_count++;

      // add protein accessions:
      set<String> protein_accessions = hit.extractProteinAccessionsSet();
      data.accessions.insert(protein_accessions.begin(), protein_accessions.end());
    }
  }


  PeptideHit PeptideAndProteinQuant::getAnnotation_(
    vector<PeptideIdentification>& peptides)
  {
    // hits in IDs must already be sorted by score! (done in "countPeptides_")
    if (peptides.empty() || peptides[0].getHits().empty()) return {};

    // get best hit
    const PeptideHit& hit = peptides[0].getHits()[0];

    // check for ambiguities
    for (auto pep_it = ++peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      const PeptideHit& current = pep_it->getHits()[0];
      if (current.getSequence() != hit.getSequence())
      {
        // TODO?: warn/error that ambiguous sequences are annotated. check if this can happen
        return {};
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
    //TODO The practice of inserting elements with the [] should be forbidden.
    // It is a debugging nightmare because if you try to access it and it is
    // not there, you are adding another element. In a next iteration this whole
    // class should be rewritten to use insert/emplace and find or better yet,
    // since we have "normal" 0-based values for samples now, vectors.
    pep_quant_[seq].abundances[fraction][hit.getCharge()][sample] +=
      feature.getIntensity(); // new map element is initialized with 0
  }

  bool PeptideAndProteinQuant::getBest_(const std::map<Int, std::map<Int, SampleAbundances>>& peptide_abundances, std::pair<size_t, size_t>& best)
  {
    size_t best_n_quant(0);
    double best_abundance(0);
    best = std::make_pair(0,0);

    for (auto & fa : peptide_abundances) // for all fractions 
    {
      for (auto & ca : fa.second) // for all charge states
      {
        const Int & fraction = fa.first;
        const Int & charge = ca.first;

        double current_abundance = std::accumulate(
          std::begin(ca.second),
          std::end(ca.second),
          0.0,
          [] (int value, const SampleAbundances::value_type& p)
          { return value + p.second; }
          ); // loop over all samples and sum abundances

        if (current_abundance <= 0) { continue; }

        const size_t current_n_quant = ca.second.size();
        if (current_n_quant > best_n_quant)
        {           
          best_abundance = current_abundance;
          best_n_quant = current_n_quant;
          best = std::make_pair(fraction, charge);
        }
        else if (current_n_quant == best_n_quant 
                 && current_abundance > best_abundance) // resolve tie by abundance
        {
          best_abundance = current_abundance;
          best = std::make_pair(fraction, charge);
        }
      }
    }
    return best_abundance > 0.;
  }


  void PeptideAndProteinQuant::quantifyPeptides(
    const vector<PeptideIdentification>& peptides)
  {
    OPENMS_LOG_INFO << "Quantifying peptides..." << std::endl;
    
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
          for (auto & a : pos->second) { OPENMS_LOG_DEBUG << a << "\t"; }
          OPENMS_LOG_DEBUG << "\n";
          pep_q.second.accessions = pos->second; // replace accessions
          filtered.insert(pep_q);
        }
        else
        {
          OPENMS_LOG_DEBUG << "not found in inference data." << endl;
        }
      }
      pep_quant_ = std::move(filtered);
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
            for (auto & sa : ca.second) // loop over all sample abundances
            {
              const UInt64 & sample_id = sa.first;
              const double & sample_abundance = sa.second;
              pep_q.second.total_abundances[sample_id] += sample_abundance;
            }
          }
        }
      }

      // for PSM counts we cover all fractions and charge states
      for (auto & fa : pep_q.second.psm_counts) // for all fractions 
      {
        for (auto & ca : fa.second) // for all charge states
        {  
          for (auto & sa : ca.second) // loop over all psm counts
          {
            const UInt64 & sample_id = sa.first;
            const double & sample_counts = sa.second;
            pep_q.second.total_psm_counts[sample_id] += sample_counts;
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
          for (auto & sa : ca.second) // loop over all sample abundances
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
      for (auto const& pg : proteins.getIndistinguishableProteins())
      {
        for (auto const& acc : pg.accessions)
        {
          // each accession should only occur once, but we don't check...
          accession_to_leader[acc] = pg.accessions[0];
        }
      }
    }

    // for (auto & a : accession_to_leader) { std::cout << a.first << "\tis led by:\t" << a.second << endl; }

    bool contains_accessions {false};

    for (auto const& pep_q : pep_quant_)
    {
      String accession = getAccession_(pep_q.second.accessions, accession_to_leader);
      OPENMS_LOG_DEBUG << "Peptide id mapped to leader: " << accession << endl;

      // not enough evidence or mapping to multiple groups
      if (accession.empty())
        continue;

      contains_accessions = true;
      // proteotypic peptide
      const String peptide = pep_q.first.toUnmodifiedString();

      prot_quant_[accession].psm_count += pep_q.second.psm_count;

      // transfer abundances and counts from peptides->protein
      // summarize abundances and counts between different peptidoforms
      for (auto const& sta : pep_q.second.total_abundances)
      {
        prot_quant_[accession].abundances[peptide][sta.first] += sta.second;
      }

      for (auto const& sta : pep_q.second.total_psm_counts)
      {
        prot_quant_[accession].psm_counts[peptide][sta.first] += sta.second;
      }
    }

    if (!contains_accessions)
    {
      OPENMS_LOG_FATAL_ERROR << "No protein matches found, cannot quantify proteins." << endl;
      throw Exception::MissingInformation(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        "No protein matches found, cannot quantify proteins.");
    }

    std::string method = param_.getValue("method");
    Size top_n = param_.getValue("top:N");
    std::string aggregate = param_.getValue("top:aggregate");
    bool include_all = param_.getValue("top:include_all") == "true";
    bool fix_peptides = param_.getValue("consensus:fix_peptides") == "true";

    if (method == "iBAQ")
    {
      top_n = 0;
      aggregate = "sum";
    }

    for (auto& prot_q : prot_quant_)
    {
      const ProteinData& pd = prot_q.second;

      // calculate PSM counts based on all (!) peptides of a protein (group)
      for (auto const& pep2sa : pd.psm_counts)
      { // for all peptides of this protein (group)
        const SampleAbundances& sas = pep2sa.second;
        for (auto const& sa : sas)
        {
          const Size& sample_id = sa.first;
          const Size& psms = sa.second;
          if (psms > 0)
            prot_q.second.total_distinct_peptides[sample_id]++; // count this peptide sequence once if observed in sample
          prot_q.second.total_psm_counts[sample_id] += psms;    // count all PSMs of this protein in this sample
        }
      }

      // select which peptides of the current protein (group) are quantified
      if ((top_n > 0) && (prot_q.second.abundances.size() < top_n))
      { // not enough proteotypic peptides? skip protein (except if user chose to include the nevertheless)
        stats_.too_few_peptides++;
        if (!include_all)
        {
          continue;
        }
      }

      vector<String> peptides; // peptides selected for quantification
      if (fix_peptides && (top_n == 0))
      {
        // consider all peptides that occur in every sample:
        for (auto const& ab : prot_q.second.abundances)
        {
          if (ab.second.size() == stats_.n_samples)
          {
            peptides.push_back(ab.first);
          }
        }
      }
      else if (fix_peptides && (top_n > 0) && (prot_q.second.abundances.size() > top_n))
      {
        orderBest_(prot_q.second.abundances, peptides);
        peptides.resize(top_n);
      }
      else
      {
        // consider all peptides of the protein:
        for (auto const& ab : prot_q.second.abundances)
        {
          peptides.push_back(ab.first);
        }
      }
      // done selecting peptides for quantification

      // consider only the selected peptides for quantification:
      map<UInt64, DoubleList> abundances; // all peptide abundances by sample
      for (const auto& pep : peptides)    // for all selected peptides
      {
        for (auto& sa : prot_q.second.abundances[pep]) // copy over all abundances
        {
          abundances[sa.first].push_back(sa.second);
        }
      }

      for (auto& ab : abundances)
      {
        // check if the protein has enough peptides in this sample
        if (!include_all && (top_n > 0) && (ab.second.size() < top_n))
        {
          continue;
        }

        // if we have more than "top", reduce to the top ones
        if ((top_n > 0) && (ab.second.size() > top_n))
        {
          // sort descending:
          sort(ab.second.begin(), ab.second.end(), greater<double>());
          ab.second.resize(top_n); // remove all but best N values
        }

        double abundance_result;
        if (aggregate == "median")
        {
          abundance_result = Math::median(ab.second.begin(), ab.second.end());
        }
        else if (aggregate == "mean")
        {
          abundance_result = Math::mean(ab.second.begin(), ab.second.end());
        }
        else if (aggregate == "weighted_mean")
        {
          double sum_intensities = 0;
          double sum_intensities_squared = 0;
          for (auto const& in : ab.second)
          {
            sum_intensities += in;
            sum_intensities_squared += in * in;
          }
          abundance_result = sum_intensities_squared / sum_intensities;
        }
        else // "sum"
        {
          abundance_result = Math::sum(ab.second.begin(), ab.second.end());
        }

        prot_q.second.total_abundances[ab.first] = abundance_result;
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
    if (method == "iBAQ")
    {
      EnzymaticDigestion digest{};
      for (auto & hit : proteins.getHits())
      {
        const OpenMS::String & hit_accession = hit.getAccession();
        const OpenMS::String & hit_sequence = hit.getSequence();

        if (prot_quant_.find(hit_accession) != prot_quant_.end())
        {
          if (hit_sequence.empty())
          {
            prot_quant_.erase(hit_accession);
            OPENMS_LOG_WARN << "Removed " << hit_accession <<  ", no protein sequence found!" << endl;
          }
          else
          {
            std::vector<StringView> peptides {};
            digest.digestUnmodified(StringView(hit_sequence), peptides);
            for (auto& total_abundance : prot_quant_[hit_accession].total_abundances)
            {
              total_abundance.second /= double(peptides.size());
            }
          }
        }
      }
    }
  }


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
      const size_t fraction(1), sample(0);
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
    // TODO check that the file section of the experimental design is compatible with what can be parsed from the consensus map.
    updateMembers_(); // clear data

    if (consensus.empty())
    {
      OPENMS_LOG_ERROR << "Empty consensus map passed to readQuantData." << endl;
      return;
    }

    // n_fractions are also used to initialize enough
    stats_.n_fractions = ed.getNumberOfFractions();
    stats_.n_ms_files = ed.getNumberOfMSFiles();
    stats_.n_samples = ed.getNumberOfSamples();

    OPENMS_LOG_DEBUG << "Reading quant data: " << endl;
    OPENMS_LOG_DEBUG << "  MS files        : " << stats_.n_ms_files << endl;
    OPENMS_LOG_DEBUG << "  Fractions       : " << stats_.n_fractions << endl;
    OPENMS_LOG_DEBUG << "  Samples (Assays): " << stats_.n_samples << endl;

   
    // map filename and label of experimental design to the full experimental design entry for faster lookup
    const auto& ms_section = ed.getMSFileSection();
    std::unordered_map<String, ExperimentalDesign::MSFileSectionEntry> fileAndLabel2MSFileSectionEntry;
    for (const auto& e : ms_section)
    {
      String ed_filename = FileHandler::stripExtension(File::basename(e.path));
      String ed_label = e.label;
      fileAndLabel2MSFileSectionEntry[ed_filename + ed_label] = e;
    }

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
        //TODO MULTIPLEXED: needs to be adapted for multiplexed experiments
        size_t row = f.getMapIndex();
        const auto& h = consensus.getColumnHeaders().at(row);
        const String c_fn = FileHandler::stripExtension(File::basename(h.filename)); // filename according to experimental design in consensus map
        const size_t c_lab = h.getLabelAsUInt(consensus.getExperimentType());

        // find entry in experimental design (ignore extension and folder) that corresponds to current column header entry
        if (auto it = fileAndLabel2MSFileSectionEntry.find(c_fn + String(c_lab)); it != fileAndLabel2MSFileSectionEntry.end())
        {
          const size_t fraction = it->second.fraction;
          const size_t sample = it->second.sample;
          quantifyFeature_(f, fraction, sample, hit); // updates "stats_.quant_features"          
        }
        else
        {
          OPENMS_LOG_FATAL_ERROR << "File+Label referenced in consensus header not found in experimental design.\n"  
                                 << "File+Label:" << c_fn << "\t" << c_lab << std::endl;
        }
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

    OPENMS_LOG_DEBUG << "Reading quant data: " << endl;
    OPENMS_LOG_DEBUG << "  MS files        : " << stats_.n_ms_files << endl;
    OPENMS_LOG_DEBUG << "  Fractions       : " << stats_.n_fractions << endl;
    OPENMS_LOG_DEBUG << "  Samples (Assays): " << stats_.n_samples << endl;

    stats_.total_features = peptides.size();
    
    countPeptides_(peptides);

    map<pair<String,Size>, String> identifier_idmergeidx_to_ms_file;
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
      for (Size s = 0; s < ms_files.size(); ++s)
      {
        identifier_idmergeidx_to_ms_file[{proteins[i].getIdentifier(), s}] = ms_files[s];
      }

      OPENMS_LOG_DEBUG << "  run index : MS file " << i << " : " << ListUtils::concatenate(ms_files, ", ") << endl;
    }

    for (auto & p : peptides)
    {
      if (p.getHits().empty()) { continue; }
      Size id_merge_idx = p.getMetaValue("id_merge_idx",0);
      const PeptideHit& hit = p.getHits()[0];

      // don't quantify decoys
      if ((std::string)hit.getMetaValue("target_decoy", DataValue("target")) == "decoy") continue;

      stats_.quant_features++;
      const AASequence& seq = hit.getSequence();
      const String& ms_file_path = identifier_idmergeidx_to_ms_file[{p.getIdentifier(),id_merge_idx}];

      // determine sample and fraction by MS file name (stored in protein identification)
      const ExperimentalDesign::MSFileSection& run_section = ed.getMSFileSection();
      auto row = find_if(begin(run_section), end(run_section), 
        [&ms_file_path](const ExperimentalDesign::MSFileSectionEntry& r)
          { 
            return File::basename(r.path) == File::basename(ms_file_path); 
          });

      if (row == end(run_section))
      {
        OPENMS_LOG_ERROR << "MS file: " << ms_file_path << " not found in experimental design." << endl;
        for (const auto& r : run_section)
        {
          OPENMS_LOG_ERROR << r.path << endl;
        }
        throw Exception::MissingInformation(
          __FILE__, 
          __LINE__, 
          OPENMS_PRETTY_FUNCTION, 
          "MS file annotated in protein identification doesn't match any in the experimental design.");
      }

      size_t sample = row->sample;
      size_t fraction = row->fraction;

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

  void PeptideAndProteinQuant::annotateQuantificationsToProteins(
    const ProteinQuant& protein_quants,
    ProteinIdentification& proteins,
    bool remove_unquantified)
  {
    auto & id_groups = proteins.getIndistinguishableProteins();

    for (const auto& q : protein_quants)
    {
      // accession of quantified protein(group)
      const String & acc = q.first;

      if (q.second.total_abundances.empty()) 
      {
        //TODO maybe just count the number of unquantifiable proteins and report that?
        OPENMS_LOG_DEBUG << "Protein " << acc << " not quantified." << endl;
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
        const SampleAbundances& total_abundances = q.second.total_abundances;
        const SampleAbundances& total_psm_counts = q.second.total_psm_counts;
        const SampleAbundances& total_distinct_peptides = q.second.total_distinct_peptides;

        // TODO: OPENMS_ASSERT(id_group->float_data_arrays.empty(), "Protein group float data array not empty!.");
        id_group->getFloatDataArrays().resize(3);
        ProteinIdentification::ProteinGroup::FloatDataArray & abundances = id_group->getFloatDataArrays()[0];
        Size n_samples = getStatistics().n_samples;
        abundances.setName("abundances");
        abundances.resize(n_samples);

        auto & psm_counts = id_group->getFloatDataArrays()[1];
        psm_counts.setName("psm_count");
        psm_counts.resize(n_samples);

        auto & peptide_counts = id_group->getFloatDataArrays()[2];
        peptide_counts.setName("distinct_peptides");
        peptide_counts.resize(n_samples);

        for (auto const & s : total_abundances)
        {
          abundances[s.first] = (float) s.second;
        }
        for (auto const & s : total_psm_counts)
        {
          psm_counts[s.first] = (float) s.second;
        }
        for (auto const & s : total_distinct_peptides)
        {
          peptide_counts[s.first] = (float) s.second;
        }
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

    if (remove_unquantified)
    {
      // remove all protein groups that have not been quantified
      auto notQuantified = [] (const ProteinIdentification::ProteinGroup& g)->bool { return g.getFloatDataArrays().empty(); };
      id_groups.erase(
          remove_if(id_groups.begin(), id_groups.end(), notQuantified),
          id_groups.end());
    }
  } 

}
