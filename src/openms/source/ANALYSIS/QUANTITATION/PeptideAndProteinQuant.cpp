// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------
//


#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <algorithm>

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
    PeptideIdentificationList& peptides)
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
    PeptideIdentificationList& peptides)
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
                                                const String& filename,
                                                const PeptideHit& hit,
                                                UInt channel_or_label)
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
    pep_quant_[seq].abundances[fraction][filename][hit.getCharge()][channel_or_label] +=
      feature.getIntensity(); // new map element is initialized with 0
  }

  bool PeptideAndProteinQuant::getBest_(const std::map<Int, std::map<String, std::map<Int, std::map<UInt, double>>>>& peptide_abundances, std::tuple<size_t, String, size_t, UInt>& best)
  {
    size_t best_n_quant(0);
    double best_abundance(0);
    best = std::make_tuple(0, "", 0, 0);

    for (auto & fa : peptide_abundances) // for all fractions
    {
      for (auto & fna : fa.second) // for all filenames
      {
        for (auto & ca : fna.second) // for all charge states
        {
          for (auto & cha : ca.second) // for all channels
          {
            const Int & fraction = fa.first;
            const String & filename = fna.first;
            const Int & charge = ca.first;
            const UInt & channel = cha.first;

            double current_abundance = cha.second;

            if (current_abundance <= 0) { continue; }

            const size_t current_n_quant = 1; // Each entry represents one quantification
            if (current_n_quant > best_n_quant)
            {
              best_abundance = current_abundance;
              best_n_quant = current_n_quant;
              best = std::make_tuple(fraction, filename, charge, channel);
            }
            else if (current_n_quant == best_n_quant
                     && current_abundance > best_abundance) // resolve tie by abundance
            {
              best_abundance = current_abundance;
              best = std::make_tuple(fraction, filename, charge, channel);
            }
          }
        }
      }
    }
    
    return best_n_quant > 0; // Return true if at least one abundance was found
  }

  size_t PeptideAndProteinQuant::getSampleIDFromFilenameAndChannel_(const String& filename,
                                                                 UInt channel_or_label,
                                                                 const ExperimentalDesign& ed) const
  {
    // Map filename and label to sample using experimental design
    const auto& ms_section = ed.getMSFileSection();
    for (const auto& entry : ms_section)
    {
      String ed_filename = FileHandler::stripExtension(File::basename(entry.path));
      if (ed_filename == filename && entry.label == channel_or_label)
      {
        return entry.sample;
      }
    }
    
    // If not found, throw an exception with detailed information
    throw Exception::MissingInformation(
      __FILE__, 
      __LINE__, 
      OPENMS_PRETTY_FUNCTION, 
      "Could not find sample mapping for filename '" + filename + "' and channel '" + String(channel_or_label) + "' in experimental design.");
  }

  void PeptideAndProteinQuant::quantifyPeptides(
    const PeptideIdentificationList& peptides)
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

        // determine which fraction, filename, charge state, and channel yields the maximum abundance
        // (break ties by total abundance)
        std::tuple<size_t, String, size_t, UInt> best_combination;

        // return false: only identified, not quantified
        if (!getBest_(pep_q.second.abundances, best_combination))
        {
          continue;
        }
        
        // quantify according to the best combination only:
        size_t best_fraction = std::get<0>(best_combination);
        String best_filename = std::get<1>(best_combination);
        size_t best_charge = std::get<2>(best_combination);
        UInt best_channel = std::get<3>(best_combination);
        
        double abundance = pep_q.second.abundances[best_fraction][best_filename][best_charge][best_channel];
        size_t sample_id = getSampleIDFromFilenameAndChannel_(best_filename, best_channel, experimental_design_);
        pep_q.second.total_abundances[sample_id] = abundance;
      }
      else
      { // sum up sample abundances over all fractions, filenames, charge states, and channels:
        for (auto & fa : pep_q.second.abundances)  // for all fractions
        {
          for (auto & fna : fa.second) // for all filenames
          {
            for (auto & ca : fna.second) // for all charge states
            {
              for (auto & cha : ca.second) // for all channels
              {
                const String & filename = fna.first;
                const UInt & channel = cha.first;
                const double & abundance = cha.second;
                
                // Map (filename, channel) to sample using ExperimentalDesign
                size_t sample_id = getSampleIDFromFilenameAndChannel_(filename, channel, experimental_design_);
                pep_q.second.total_abundances[sample_id] += abundance;
              }
            }
          }
        }
      }

      // for PSM counts we cover all fractions, filenames, charge states.
      for (auto & fa : pep_q.second.psm_counts) // for all fractions
      {
        for (auto & fna : fa.second) // for all filenames
        {
          for (auto & ca : fna.second) // for all charge states
          {
            const String & filename = fna.first;
            const double & psm_counts = ca.second;
            
            // In multiplexed design, e.g. TMT, a signle PSM is associated with all samples measured in the different channels/labels 
            for (Size channel = 1; channel <= experimental_design_.getNumberOfLabels(); ++channel)
            {
              size_t sample_id = getSampleIDFromFilenameAndChannel_(filename, channel, experimental_design_);              
              pep_q.second.total_psm_counts[sample_id] += psm_counts; // accumulate PSM counts for spectral counting
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
        for (auto & fna : fa.second) // for all filenames
        {
          for (auto & ca : fna.second) // for all charge states
          {
            for (auto & cha : ca.second) // for all channels
            {
              const String & filename = fna.first;
              const UInt & channel = cha.first;
              size_t sample_id = getSampleIDFromFilenameAndChannel_(filename, channel, experimental_design_);
              cha.second *= scale_factors[sample_id];
            }
          }
        }
      }
    }
  }

  String PeptideAndProteinQuant::getAccession_(
    const set<String>& pep_accessions, 
    const map<String, String>& accession_to_leader) const
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


  void PeptideAndProteinQuant::quantifyProteins(const ProteinIdentification& proteins)
  {
    if (pep_quant_.empty())
    {
      OPENMS_LOG_WARN << "Warning: No peptides quantified." << endl;
      return;
    }

    // Phase 1: Transfer peptide data to protein structures
    transferPeptideDataToProteins_(proteins);

    // Phase 2: Extract and validate parameters
    std::string method = param_.getValue("method");
    Size top_n = param_.getValue("top:N");
    std::string aggregate = param_.getValue("top:aggregate");
    bool include_all = param_.getValue("top:include_all") == "true";
    bool fix_peptides = param_.getValue("consensus:fix_peptides") == "true";

    // Handle iBAQ parameter overrides
    if (method == "iBAQ")
    {
      top_n = 0;
      aggregate = "sum";
    }

    // Phase 3: Process each protein
    for (auto& prot_q : prot_quant_)
    {
      const String& accession = prot_q.first;
      const ProteinData& pd = prot_q.second;

      // Calculate PSM counts based on all peptides of a protein (group)
      for (auto const& pep2sa : pd.peptide_psm_counts)
      {
        const SampleAbundances& sas = pep2sa.second;
        for (auto const& sa : sas)
        {
          const Size& sample_id = sa.first;
          const Size& psms = sa.second;
          if (psms > 0)
            prot_q.second.total_distinct_peptides[sample_id]++;
          prot_q.second.total_psm_counts[sample_id] += psms;
        }
      }

      // Check if protein has enough peptides (for statistics)
      if ((top_n > 0) && (prot_q.second.peptide_abundances.size() < top_n))
      {
        stats_.too_few_peptides++;
        if (!include_all)
        {
          continue;
        }
      }

      // Select peptides for quantification
      std::vector<String> selected_peptides = selectPeptidesForQuantification_(
          accession, top_n, fix_peptides);

      // Calculate protein abundances
      calculateProteinAbundances_(accession, selected_peptides, aggregate, top_n, include_all);

      // if information about (indistinguishable) protein groups is available, map
      // each accession to the accession of the leader of its group of proteins:
      auto accession_to_leader = mapAccessionToLeader(proteins);
      
      calculateFileAndChannelLevelProteinAbundances_(accession, selected_peptides, aggregate,
                                          top_n, include_all, accession_to_leader);

      // Update statistics
      if (prot_q.second.total_abundances.empty())
      {
        stats_.too_few_peptides++;
      }
      else
      {
        stats_.quant_proteins++;
      }
    }

    // Phase 4: Post-processing
    if (method == "iBAQ")
    {
      performIbaqNormalization_(proteins);
    }
  }


  std::map<OpenMS::String, OpenMS::String> PeptideAndProteinQuant::mapAccessionToLeader(const OpenMS::ProteinIdentification& proteins) const
  {
    std::map<OpenMS::String, OpenMS::String> accession_to_leader;
    if (! proteins.getIndistinguishableProteins().empty())
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
    return accession_to_leader;
  }
  void PeptideAndProteinQuant::readQuantData(FeatureMap& features, const ExperimentalDesign& ed)
  {
    updateMembers_(); // clear data
    experimental_design_ = ed; // store experimental design for aggregation

    stats_.n_samples = ed.getNumberOfSamples();
    stats_.n_fractions = 1;
    stats_.n_ms_files = ed.getNumberOfMSFiles();

    stats_.total_features = features.size();

    // For FeatureMap, extract filename from metadata or use default
    String filename = "default";
    if (features.metaValueExists("filename"))
    {
      filename = FileHandler::stripExtension(File::basename(features.getMetaValue("filename")));
    }
    else if (!ed.getMSFileSection().empty())
    {
      // Use first MS file from experimental design as fallback
      filename = FileHandler::stripExtension(File::basename(ed.getMSFileSection()[0].path));
    }

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
      const size_t fraction(1);
      const Int label(1); // Default label for LFQ data
      quantifyFeature_(handle, fraction, filename, hit, label); // updates "stats_.quant_features"
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
    experimental_design_ = ed; // store experimental design for aggregation

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
          quantifyFeature_(f, fraction, c_fn, hit, c_lab); // updates "stats_.quant_features"
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
    std::vector<ProteinIdentification>& proteins,
    PeptideIdentificationList& peptides,
    const ExperimentalDesign& ed)
  {
    updateMembers_(); // clear data
    experimental_design_ = ed; // store experimental design for aggregation

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
      if (hit.isDecoy()) continue;

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

      size_t fraction = row->fraction;
      String filename = FileHandler::stripExtension(File::basename(ms_file_path));
      Int label = row->label; // Use label from experimental design

      // count peptides in the different fractions, filenames, charge states, and channels
      pep_quant_[seq].abundances[fraction][filename][hit.getCharge()][label] += 1;
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
    // read experimental design as it is needed to annotate quantities in the correct order
    ExperimentalDesign::MSFileSection msfile_section = experimental_design_.getMSFileSection();
    
    // Extract the Spectra Filepath column from the design
    map<UInt64, map<UInt64, String>> design_group_fraction_filename;
    UInt64 n_files = 0;
    for (ExperimentalDesign::MSFileSectionEntry const& f : msfile_section)
    {
      const String fn = FileHandler::stripExtension(File::basename(f.path));
      design_group_fraction_filename[f.fraction_group][f.fraction] = fn;
      n_files++;
    }

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
        const auto& file_level_psm_counts = q.second.file_level_psm_counts;
      
        // TODO: OPENMS_ASSERT(id_group->float_data_arrays.empty(), "Protein group float data array not empty!.");
        id_group->getFloatDataArrays().resize(4);
        id_group->getStringDataArrays().resize(2);
        id_group->getIntegerDataArrays().resize(2);
        
        // Sample-level arrays (indices 0-2)
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

        // Add file/channel level abundances
        auto& file_channel_level_abundance = id_group->getFloatDataArrays()[3];
        file_channel_level_abundance.setName("file_channel_level_abundance");
        auto& file_channel_level_filename = id_group->getStringDataArrays()[0];
        file_channel_level_filename.setName("file_channel_level_filename");
        auto& file_channel_level_channel = id_group->getIntegerDataArrays()[0];
        file_channel_level_channel.setName("file_channel_level_channel");
     

        // We loop over the filenames in the design file, as this is the order we expect in the output.
        for (const auto& [group_id, fraction_to_filename_map] : design_group_fraction_filename)
        {
          for (auto [fraction, design_filename] : fraction_to_filename_map)
          {
            // Process each filename within the fraction group
            // important: strip file extension and path to find the entry
            design_filename = FileHandler::stripExtension(File::basename(design_filename));
            
            #ifdef DEBUG_PROTEINQUANTIFIER
            std::cout 
              << "Experimental design: fraction group: " << group_id 
              << ", filename: '" << design_filename
              << "', fraction: " << fraction
              << " of the experimental design." << std::endl;
            #endif

            // for each file in the design, fill the channels quantity
            for (Size c = 1; c <= experimental_design_.getNumberOfLabels(); ++c) // label/channel numbers are 1-based
            {
              double channel_abundance{};

              const auto& filename_to_channel_map = q.second.channel_level_abundances;

              if (auto file_level_it = filename_to_channel_map.find(design_filename); 
                file_level_it != filename_to_channel_map.end())
              {
                if (file_level_it->second.find(0) != file_level_it->second.end()) throw Exception::MissingInformation(
                  __FILE__, 
                  __LINE__, 
                  OPENMS_PRETTY_FUNCTION, 
                  "Channel found that should not exist.");

                // Found the file, now search for the channel
                if (auto channel_it = file_level_it->second.find(c);
                    channel_it != file_level_it->second.end())
                {
                  channel_abundance = channel_it->second;
                }
              }

              file_channel_level_abundance.push_back(channel_abundance);
              file_channel_level_filename.push_back(design_filename);
              file_channel_level_channel.push_back(c);
         
              #if DEBUG_PEPTIDEANDPROTEINQUANT
              std::cout << "DEBUG: Adding abundance for protein to meta value " << acc
                        << " filename " << design_filename
                        << " channel " << c
                        << ": " << channel_abundance << endl;
              #endif

            }    
          }
        }

        // Add file level PSM counts
        auto& file_level_psm_count = id_group->getIntegerDataArrays()[1];
        file_level_psm_count.setName("file_level_psm_count");
        auto& file_level_filename = id_group->getStringDataArrays()[1];
        file_level_filename.setName("file_level_filename");

        for (const auto& filename : file_level_psm_counts)
        {
          file_level_psm_count.push_back((int)filename.second);
          file_level_filename.push_back(filename.first);
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


  void PeptideAndProteinQuant::transferPeptideDataToProteins_(const ProteinIdentification& proteins)
  {
    // if information about (indistinguishable) protein groups is available, map
    // each accession to the accession of the leader of its group of proteins:
    map<String, String> accession_to_leader = mapAccessionToLeader(proteins);

    bool contains_accessions{ false}; // flag to check if any accessions were found

    for (auto const& pep_q : pep_quant_)
    {
      String leader_accession = getAccession_(pep_q.second.accessions, accession_to_leader);
      OPENMS_LOG_DEBUG << "Peptide id mapped to leader: " << leader_accession << endl;

      // not enough evidence or mapping to multiple groups
      if (leader_accession.empty())
        continue;

      contains_accessions = true;
      // proteotypic peptide
      const String peptide = pep_q.first.toUnmodifiedString();

      prot_quant_[leader_accession].psm_count += pep_q.second.psm_count; // total PSM count for this group of proteins (represented by the leader accession)

      // transfer abundances and counts from peptides->protein
      // summarize abundances and counts between different peptidoforms
      for (auto const& sta : pep_q.second.total_abundances)
      {
        prot_quant_[leader_accession].peptide_abundances[peptide][sta.first] += sta.second;
      }

      for (auto const& sta : pep_q.second.total_psm_counts)
      {
        prot_quant_[leader_accession].peptide_psm_counts[peptide][sta.first] += sta.second;
      }

      // transfer detailed abundances from peptide to protein
      for (auto const& fraction : pep_q.second.abundances)
      {
        for (auto const& filename : fraction.second)
        {
          for (auto const& charge : filename.second)
          {
            for (auto const& channel : charge.second)
            {
              prot_quant_[leader_accession].channel_level_abundances[filename.first][channel.first] += channel.second;
#ifdef DEBUG_PROTEINQUANTIFIER                
              std::cout << "DEBUG: Adding abundance for protein " << accession
                        << " fraction " << fraction.first
                        << " filename " << filename.first
                        << " charge " << charge.first
                        << " channel " << channel.first
                        << ": " << channel.second << endl; 
#endif
            }
          }
        }
      }

      // transfer detailed PSM counts from peptide to protein
      for (auto const& fraction : pep_q.second.psm_counts)
      {
        for (auto const& filename : fraction.second)
        {
          for (auto const& charge : filename.second)
          {
            prot_quant_[leader_accession].file_level_psm_counts[filename.first] += (UInt)charge.second;
          }
        }
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
  }

  std::vector<String> PeptideAndProteinQuant::selectPeptidesForQuantification_(const String& protein_accession,
                                                                              Size top_n,
                                                                              bool fix_peptides)
  {
    std::vector<String> peptides;
    
    auto prot_it = prot_quant_.find(protein_accession);
    if (prot_it == prot_quant_.end())
    {
      return peptides; // empty vector
    }
    
    const ProteinData& pd = prot_it->second;

    if (fix_peptides && (top_n == 0))
    {
      // consider all peptides that occur in every sample:
      for (auto const& ab : pd.peptide_abundances)
      {
        if (ab.second.size() == stats_.n_samples)
        {
          peptides.push_back(ab.first);
        }
      }
    }
    else if (fix_peptides && (top_n > 0) && (pd.peptide_abundances.size() > top_n))
    {
      orderBest_(pd.peptide_abundances, peptides);
      peptides.resize(top_n);
    }
    else
    {
      // consider all peptides of the protein:
      for (auto const& ab : pd.peptide_abundances)
      {
        peptides.push_back(ab.first);
      }
    }

    return peptides;
  }

  double PeptideAndProteinQuant::aggregateAbundances_(const std::vector<double>& abundances,
                                                     const String& method) const
  {
    if (abundances.empty())
    {
      return 0.0;
    }

    if (method == "median")
    {
      std::vector<double> sorted_abundances = abundances; // make a copy for sorting
      return Math::median(sorted_abundances.begin(), sorted_abundances.end());
    }
    else if (method == "mean")
    {
      return Math::mean(abundances.begin(), abundances.end());
    }
    else if (method == "weighted_mean")
    {
      double sum_intensities = 0;
      double sum_intensities_squared = 0;
      for (auto const& intensity : abundances)
      {
        sum_intensities += intensity;
        sum_intensities_squared += intensity * intensity;
      }
      return sum_intensities_squared / sum_intensities;
    }
    else // "sum"
    {
      return Math::sum(abundances.begin(), abundances.end());
    }
  }

  void PeptideAndProteinQuant::calculateProteinAbundances_(const String& protein_accession,
                                                          const std::vector<String>& selected_peptides,
                                                          const String& aggregate_method,
                                                          Size top_n,
                                                          bool include_all)
  {
    auto prot_it = prot_quant_.find(protein_accession);
    if (prot_it == prot_quant_.end())
    {
      return;
    }

    ProteinData& pd = prot_it->second;

    // consider only the selected peptides for quantification:
    map<UInt64, DoubleList> abundances; // all peptide abundances by sample
    for (const auto& pep : selected_peptides)    // for all selected peptides
    {
      auto pep_it = pd.peptide_abundances.find(pep);
      if (pep_it != pd.peptide_abundances.end())
      {
        for (auto& sa : pep_it->second) // copy over all abundances
        {
          abundances[sa.first].push_back(sa.second);
        }
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

      double abundance_result = aggregateAbundances_(ab.second, aggregate_method);
      pd.total_abundances[ab.first] = abundance_result;
    }
  }

  void PeptideAndProteinQuant::calculateFileAndChannelLevelProteinAbundances_(const String& protein_accession,
                                                                  const std::vector<String>& selected_peptides,
                                                                  const String& aggregate_method,
                                                                  Size top_n,
                                                                  bool include_all,
                                                                  const std::map<String, String>& accession_to_leader)
  {
    auto prot_it = prot_quant_.find(protein_accession);
    if ( prot_it == prot_quant_.end())
    {
      return;
    }

    ProteinData& pd = prot_it->second;

    // organize detailed abundances by (fraction, filename, channel) combinations
    map<tuple<Int, String, UInt>, DoubleList> channel_level_abundances_for_selected_peptides;
    
    // collect detailed abundances from selected peptides
    for (const auto& pep : selected_peptides)    // for all selected peptides
    {
      // find the original peptide data to get detailed abundances
      for (auto const& pep_q_check : pep_quant_)
      {
        if (pep_q_check.first.toUnmodifiedString() == pep)
        {
          String check_accession = getAccession_(pep_q_check.second.accessions, accession_to_leader);
          if (check_accession == protein_accession) // this peptide belongs to current protein
          {
            // collect detailed abundances from this peptide
            for (auto const& fraction : pep_q_check.second.abundances)
            {
              for (auto const& filename : fraction.second)
              {
                for (auto const& charge : filename.second)
                {
                  for (auto const& channel : charge.second)
                  {
                    auto peptide = make_tuple(fraction.first, filename.first, channel.first);
                    channel_level_abundances_for_selected_peptides[peptide].push_back(channel.second);

                    #ifdef DEBUG_PEPTIDEANDPROTEINQUANT
                    std::cout << "DEBUG: Adding abundance for leader " <<
                              getAccession_(pep_q_check.second.accessions, const_cast<std::map<String, String>&>(accession_to_leader))
                              << pep
                              << " fraction " << fraction.first
                              << " filename " << filename.first
                              << " charge " << charge.first
                              << " channel " << channel.first
                              << ": " << channel.second << endl;
                    #endif
                  }
                }
              }
            }
            break; // found the peptide, no need to continue searching
          }
        }
      }
    }
    
    // now aggregate using the same aggregation method
    for (auto& detailed_ab : channel_level_abundances_for_selected_peptides)
    {
      const auto& selected_peptide = detailed_ab.first;
      String filename = get<1>(selected_peptide);
      UInt channel = get<2>(selected_peptide);
      
      DoubleList& all_abundances = detailed_ab.second;
      
      if (all_abundances.empty()) continue;
      
      // check if we have enough peptides for this detailed key
      if (!include_all && (top_n > 0) && (all_abundances.size() < top_n))
      {
        continue;
      }
      
      // if we have more than "top", reduce to the top ones
      if ((top_n > 0) && (all_abundances.size() > top_n))
      {
        // sort descending:
        sort(all_abundances.begin(), all_abundances.end(), greater<double>());
        all_abundances.resize(top_n); // remove all but best N values
      }
      
      double abundance_result = aggregateAbundances_(all_abundances, aggregate_method);

      // store the aggregated result in channel_level_abundances
      pd.channel_level_abundances[filename][channel] = abundance_result;
      #ifdef DEBUG_PEPTIDEANDPROTEINQUANT
      Int fraction = get<0>(selected_peptide);
      std::cout << "DEBUG: Protein " << protein_accession
                << " leader " << getAccession_(protein_accession, const_cast<std::map<String, String>&>(accession_to_leader))
                << " fraction " << fraction
                << " filename " << filename
                << " channel " << channel
                << ": " << abundance_result << endl;
      #endif
    }
  }

  void PeptideAndProteinQuant::performIbaqNormalization_(const ProteinIdentification& proteins)
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
