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
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <set>
#include <string_view>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
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
    best_charge_by_peptidoform_(), prot_quant_()
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


    defaults_.setValue("best_charge_and_fraction", "false", "Distinguish between fraction and charge states in detailed peptide output. For protein quantification, select one charge per modified peptide globally: maximize the number of samples with a positive abundance, then break ties by total abundance; retain that charge's values in every sample and sum them over fractions.\nBy default, protein abundances are summed over all charge states.");
    defaults_.setValidStrings("best_charge_and_fraction", true_false);

    defaults_.setValue("consensus:normalize", "false", "Scale peptide abundances so that the median of each sample matches the overall median.\nAbundances of zero count as 'not detected' and are left out of the medians; a sample without any positive abundance takes no part in the normalization.");
    defaults_.setValidStrings("consensus:normalize", true_false);

    defaults_.setValue("consensus:fix_peptides", "false", "Use the same peptides for protein quantification across all samples.\nWith 'N 0',"
     "all peptides that occur in every sample are considered.\nOtherwise ('N'), the N peptides that occur in the most samples (independently of each other) are selected,\nbreaking ties by total abundance (there is no guarantee that the best co-ocurring peptides are chosen!).\nA peptide counts as occurring in a sample only where its abundance is positive: an abundance stored as zero means 'not detected' (e.g. an isobaric reporter below 'min_reporter_intensity'), not a measurement of absence.");
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
      set<std::string> protein_accessions = hit.extractProteinAccessionsSet();
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
                                                const std::string& filename,
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

  bool PeptideAndProteinQuant::getBestCharge_(
    const std::map<Int, std::map<std::string, std::map<Int, std::map<UInt, double>>>>& peptide_abundances,
    Int& best_charge) const
  {
    // Collapse every charge to sample-level abundances before ranking it. Multiple fractions,
    // files, or channels that resolve to the same sample count only once for prevalence.
    std::map<Int, SampleAbundances> charge_abundances;
    for (const auto& fa : peptide_abundances) // for all fractions
    {
      for (const auto& fna : fa.second) // for all filenames
      {
        for (const auto& ca : fna.second) // for all charge states
        {
          for (const auto& cha : ca.second) // for all channels
          {
            if (cha.second <= 0.0) { continue; }
            const Size sample = getDesignCellFromFilenameAndChannel_(fna.first, cha.first).sample;
            charge_abundances[ca.first][sample] += cha.second;
          }
        }
      }
    }

    bool found = false;
    Size best_sample_count = 0;
    double best_total_abundance = 0.0;
    for (const auto& [charge, sample_abundances] : charge_abundances)
    {
      double total_abundance = 0.0;
      for (const auto& sample_abundance : sample_abundances)
      {
        total_abundance += sample_abundance.second;
      }

      const Size sample_count = sample_abundances.size();
      if (!found || sample_count > best_sample_count ||
          (sample_count == best_sample_count && total_abundance > best_total_abundance))
      {
        found = true;
        best_charge = charge;
        best_sample_count = sample_count;
        best_total_abundance = total_abundance;
      }
      // Iteration is charge-sorted, so an exact tie deliberately keeps the lower charge.
    }

    return found;
  }

  void PeptideAndProteinQuant::buildSampleIDLookup_()
  {
    // Build the design-cell lookups once, so that per-peptide aggregation does not have to
    // linearly scan the MS file section (and recompute File::stemName) for every value.
    design_cell_lookup_.clear();
    fraction_group_label_to_sample_.clear();
    quantification_fraction_group_labels_.clear();
    std::map<std::string, UInt> run_to_fraction_group;
    const Size n_samples = experimental_design_.getNumberOfSamples();
    for (const auto& entry : experimental_design_.getMSFileSection())
    {
      const std::string run = File::stemName(entry.path);
      const auto file_label = std::make_pair(run, entry.label);
      const DesignCell cell{entry.sample, entry.fraction_group};

      if (entry.sample >= n_samples)
      {
        throw Exception::MissingInformation(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Experimental-design cell ('" + run + "', label " + StringUtils::toStr(entry.label)
          + ") refers to sample " + StringUtils::toStr(entry.sample) + ", but the sample section contains "
          + StringUtils::toStr(n_samples) + " sample(s).");
      }

      auto [cell_it, cell_inserted] = design_cell_lookup_.emplace(file_label, cell);
      if (!cell_inserted &&
          (cell_it->second.sample != cell.sample || cell_it->second.fraction_group != cell.fraction_group))
      {
        throw Exception::MissingInformation(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Experimental design assigns conflicting sample or fraction-group coordinates to ('"
          + run + "', label " + StringUtils::toStr(entry.label) + ").");
      }

      auto [run_it, run_inserted] = run_to_fraction_group.emplace(run, entry.fraction_group);
      if (!run_inserted && run_it->second != entry.fraction_group)
      {
        throw Exception::MissingInformation(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Experimental design assigns run '" + run + "' to fraction groups "
          + StringUtils::toStr(run_it->second) + " and " + StringUtils::toStr(entry.fraction_group) + ".");
      }

      const auto unit_label = std::make_pair(entry.fraction_group, entry.label);
      quantification_fraction_group_labels_.insert(unit_label);
      auto [unit_it, unit_inserted] = fraction_group_label_to_sample_.emplace(unit_label, entry.sample);
      if (!unit_inserted && unit_it->second != entry.sample)
      {
        throw Exception::MissingInformation(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Label " + StringUtils::toStr(entry.label) + " of fraction group "
          + StringUtils::toStr(entry.fraction_group) + " resolves to samples "
          + StringUtils::toStr(unit_it->second) + " and " + StringUtils::toStr(entry.sample) + ".");
      }
    }
  }

  const PeptideAndProteinQuant::DesignCell& PeptideAndProteinQuant::getDesignCellFromFilenameAndChannel_(
    const std::string& filename, UInt channel_or_label) const
  {
    if (auto it = design_cell_lookup_.find({filename, channel_or_label}); it != design_cell_lookup_.end())
    {
      return it->second;
    }

    // If not found, throw an exception with detailed information
    throw Exception::MissingInformation(
      __FILE__,
      __LINE__,
      OPENMS_PRETTY_FUNCTION,
      "Could not find sample mapping for filename '" + filename + "' and channel '" + StringUtils::toStr(channel_or_label) + "' in experimental design.");
  }

  Size PeptideAndProteinQuant::getSampleIDFromFilenameAndChannel_(const std::string& filename,
                                                                  UInt channel_or_label) const
  {
    return getDesignCellFromFilenameAndChannel_(filename, channel_or_label).sample;
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
    map<std::string, set<std::string> > pep_info;
    for (auto const & pep : peptides)
    {
      for (auto const & hit : pep.getHits())
      {
        std::string seq = hit.getSequence().toUnmodifiedString();
        set<std::string> accessions = hit.extractProteinAccessionsSet();

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
        std::string seq = pep_q.first.toUnmodifiedString();
        OPENMS_LOG_DEBUG << "Sequence: " << seq << endl;
        map<std::string, set<std::string> >::iterator pos = pep_info.find(seq);
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
    // number of labels/channels is constant across the design; compute once
    // (getNumberOfLabels() scans the whole MS file section on every call).
    const Size n_labels = experimental_design_.getNumberOfLabels();
    for (auto & pep_q : pep_quant_)
    {
      if (param_.getValue("best_charge_and_fraction") == "true")
      { // quantify according to the best charge state only:
        Int best_charge = 0;
        // A peptidoform without a single positive observation has no charge to select, so it stays
        // unquantified. Do not skip the rest of the loop body for it: the PSM-count aggregation
        // below is identification evidence and must not depend on the charge-selection flag - an
        // isobaric peptide whose reporters are all 0.0 (undetected or below
        // 'min_reporter_intensity') would otherwise silently lose its spectral counts, and with
        // them the protein's 'psm_count'/'distinct_peptides' in mzTab and the QPX pg row.
        if (getBestCharge_(pep_q.second.abundances, best_charge))
        {
          best_charge_by_peptidoform_[pep_q.first] = best_charge;

          // Retain every observation of the selected charge. Sample abundances span all fractions;
          // fraction-group/label abundances retain the QPX quantification grain. Explicit zeros are
          // kept in both maps but did not count as evidence when the charge was selected above.
          for (const auto& fa : pep_q.second.abundances)
          {
            for (const auto& fna : fa.second)
            {
              auto charge_it = fna.second.find(best_charge);
              if (charge_it == fna.second.end()) { continue; }
              for (const auto& [channel, abundance] : charge_it->second)
              {
                const DesignCell& cell = getDesignCellFromFilenameAndChannel_(fna.first, channel);
                pep_q.second.total_abundances[cell.sample] += abundance;
                pep_q.second.fraction_group_abundances[cell.fraction_group][channel] += abundance;
              }
            }
          }
        }
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
                const std::string & filename = fna.first;
                const UInt & channel = cha.first;
                const double & abundance = cha.second;
                
                // Retain both grains. The sample value remains the established mzTab/API result;
                // the fraction-group/label value is the QPX 1.1 protein-group quantification unit.
                const DesignCell& cell = getDesignCellFromFilenameAndChannel_(filename, channel);
                pep_q.second.total_abundances[cell.sample] += abundance;
                pep_q.second.fraction_group_abundances[cell.fraction_group][channel] += abundance;
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
            const std::string & filename = fna.first;
            const double & psm_counts = ca.second;
            
            // In multiplexed design, e.g. TMT, a signle PSM is associated with all samples measured in the different channels/labels 
            for (Size channel = 1; channel <= n_labels; ++channel)
            {
              size_t sample_id = getSampleIDFromFilenameAndChannel_(filename, channel);              
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
    // Zero abundances are deliberately not collected. A zero is not a measurement of "no protein":
    // IsobaricChannelExtractor stores a reporter it could not find - or one below
    // 'min_reporter_intensity' - as 0.0, so a zero means "not detected". Counting those would give a
    // channel in which more than half the peptides went undetected a median of exactly 0, and the
    // division below then scaled every positive abundance in it to infinity (and the zeros to NaN,
    // since 0 * inf is NaN). Excluding them instead normalizes such a channel on the data it does
    // have. Keeping zeros out of the statistic is the convention elsewhere: MSstatsTMT turns them
    // into NA before taking median(log2Intensity, na.rm = TRUE), isobar nulls sub-1 values after
    // impurity correction, and scp calls zeroIsNA() first, "to avoid artefacts in downstream steps".
    // Note this is about the zeros, not about log versus linear space - isobar and scp both divide
    // in linear space, as here, and are safe purely because no zero ever reaches the median.
    map<UInt64, DoubleList> abundances; // positive peptide abundances by sample
    set<UInt64> samples; // every sample seen, including ones without a single positive abundance
    for (auto & pq : pep_quant_)
    {
      for (auto & sa : pq.second.total_abundances)
      {
        samples.insert(sa.first);
        if (sa.second > 0.0) { abundances[sa.first].push_back(sa.second); }
      }
    }
    if (abundances.size() <= 1)
    { // no reference to scale to - say so, rather than returning as if normalization had happened
      OPENMS_LOG_WARN << "Warning: fewer than two samples contain a peptide with a positive abundance, "
                      << "so there is nothing to normalize against. Abundances are left unchanged."
                      << endl;
      return;
    }

    /////////////////////////////////////////////////////
    // compute scale factors on the sample level:
    SampleAbundances medians; // median abundance by sample
    for (auto & ab : abundances)
    {
      medians[ab.first] = Math::median(ab.second.begin(), ab.second.end());
    }

    // Only samples that carry signal define the reference. Letting a dead channel contribute its
    // median of 0 would drag the overall median down and mis-scale every other channel - quietly,
    // with no infinity and no warning to show for it.
    DoubleList all_medians;
    for (auto & sa : medians)
    {
      all_medians.push_back(sa.second);
    }
    double overall_median = Math::median(all_medians.begin(),
                                         all_medians.end());
    SampleAbundances scale_factors;
    for (UInt64 sample : samples)
    {
      auto med_it = medians.find(sample);
      if (med_it == medians.end()) // not a single positive abundance in this sample
      {
        // Name the (file, channel) pairs behind the sample: the sample id is a design-internal index
        // that does not appear in any input the user wrote.
        std::string where;
        for (const auto& lookup : design_cell_lookup_)
        {
          if (lookup.second.sample != sample) { continue; }
          if (!where.empty()) { where += ", "; }
          where += "'" + lookup.first.first + "' channel " + StringUtils::toStr(lookup.first.second);
        }
        OPENMS_LOG_WARN << "Warning: sample " << sample
                        << (where.empty() ? "" : " (" + where + ")")
                        << " has no peptide with a positive abundance - for isobaric data, no reporter "
                        << "ion was detected for it in any quantified peptide. It does not take part in "
                        << "the normalization: it defines no reference and its values are left as they "
                        << "are (they are all zero)." << endl;
        continue; // no factor at all - see scaleFactorFor below
      }
      // a median over strictly positive values is itself positive, so this cannot divide by zero
      scale_factors[sample] = overall_median / med_it->second;
    }

    // Never look a factor up with operator[]: that default-constructs 0.0 for an unknown sample and
    // would silently wipe its abundances. A sample without any positive signal has no factor but can
    // still reach the loops below through explicit zero observations. Leave those exactly as they are.
    auto scaleFactorFor = [&scale_factors](size_t sample_id)
    {
      auto it = scale_factors.find(sample_id);
      return (it == scale_factors.end()) ? 1.0 : it->second;
    };

    /////////////////////////////////////////////////////
    // scale all abundance values:
    for (auto & pep_q : pep_quant_)
    {
      // scale total abundances
      for (auto & sta : pep_q.second.total_abundances)
      {
        sta.second *= scaleFactorFor(sta.first);
      }

      // Scale each fraction-group/label quantity with the factor of the sample to which that
      // design cell resolves. This keeps a sample split across technical-replicate fraction groups
      // on the same normalization scale without merging those independently quantified units.
      for (auto& [fraction_group, label_abundances] : pep_q.second.fraction_group_abundances)
      {
        for (auto& [label, abundance] : label_abundances)
        {
          const auto sample_it = fraction_group_label_to_sample_.find({fraction_group, label});
          if (sample_it == fraction_group_label_to_sample_.end())
          {
            throw Exception::MissingInformation(
              __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Could not resolve fraction group " + StringUtils::toStr(fraction_group) + ", label "
              + StringUtils::toStr(label) + " to a sample for normalization.");
          }
          abundance *= scaleFactorFor(sample_it->second);
        }
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
              const std::string & filename = fna.first;
              const UInt & channel = cha.first;
              size_t sample_id = getSampleIDFromFilenameAndChannel_(filename, channel);
              cha.second *= scaleFactorFor(sample_id);
            }
          }
        }
      }
    }
  }

  std::string PeptideAndProteinQuant::getAccession_(
    const set<std::string>& pep_accessions, 
    const map<std::string, std::string>& accession_to_leader) const
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
        map<std::string, std::string>::const_iterator pos = accession_to_leader.find(acc);
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

    // The accession -> group-leader map depends only on @p proteins and is loop-invariant;
    // build it once here instead of rebuilding it inside the per-protein loop (which made
    // this phase O(P^2) in the number of proteins).
    const std::map<std::string, std::string> accession_to_leader = mapAccessionToLeader(proteins);

    // Reverse index: unmodified peptide -> all pep_quant_ entries (modified peptidoforms)
    // sharing that unmodified sequence. Built once by a single forward pass over pep_quant_.
    // This replaces the previous full rescan of pep_quant_ (with a toUnmodifiedString() per
    // entry) that calculateFileAndChannelLevelProteinAbundances_ performed for every
    // (protein, selected peptide) pair - the dominant O(N_pep^2) cost on large inputs.
    UnmodifiedToEntriesIndex unmod_to_entries;
    for (const auto& pep_q : pep_quant_)
    {
      unmod_to_entries[pep_q.first.toUnmodifiedString()].push_back(&pep_q);
    }

    for (auto& prot_q : prot_quant_)
    {
      const std::string& accession = prot_q.first;
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
      std::vector<std::string> selected_peptides = selectPeptidesForQuantification_(
          accession, top_n, fix_peptides);

      // Calculate protein abundances
      calculateProteinAbundances_(accession, selected_peptides, aggregate, top_n, include_all);
      calculateFractionGroupLevelProteinAbundances_(accession, selected_peptides, aggregate,
                                                    top_n, include_all);

      // accession_to_leader and unmod_to_entries are loop-invariant and computed once before
      // the loop (see above).
      calculateFileAndChannelLevelProteinAbundances_(accession, selected_peptides, aggregate,
                                          top_n, include_all, accession_to_leader, unmod_to_entries);

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


  std::map<std::string, std::string> PeptideAndProteinQuant::mapAccessionToLeader(const OpenMS::ProteinIdentification& proteins) const
  {
    std::map<std::string, std::string> accession_to_leader;
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
    buildSampleIDLookup_(); // precompute (basename, channel) -> sample lookup

    stats_.n_samples = ed.getNumberOfSamples();
    stats_.n_fractions = 1;
    stats_.n_ms_files = ed.getNumberOfMSFiles();

    stats_.total_features = features.size();

    // For FeatureMap, extract filename from metadata or use default
    std::string filename = "default";
    if (features.metaValueExists("filename"))
    {
      filename = File::stemName(features.getMetaValue("filename").toString());
    }
    else if (!ed.getMSFileSection().empty())
    {
      // Use first MS file from experimental design as fallback
      filename = File::stemName(ed.getMSFileSection()[0].path);
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
    buildSampleIDLookup_(); // precompute (basename, channel) -> sample lookup

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
    using FileAndLabel = std::pair<std::string, UInt>;
    // Pure lookup (only emplace + find below, never iterated), so an unordered map gives
    // O(1) access per consensus feature without affecting output order.
    std::unordered_map<FileAndLabel, ExperimentalDesign::MSFileSectionEntry, FileLabelHash> file_and_label_to_msfile_entry;
    for (const auto& e : ms_section)
    {
      const std::string ed_filename = File::stemName(e.path);
      const FileAndLabel key(ed_filename, e.label);
      const auto [it, inserted] = file_and_label_to_msfile_entry.emplace(key, e);
      if (!inserted &&
          ((it->second.fraction != e.fraction) ||
           (it->second.fraction_group != e.fraction_group) ||
           (it->second.sample != e.sample)))
      {
        throw Exception::MissingInformation(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Ambiguous basename+label mapping in experimental design for '" + ed_filename +
          "' and label '" + StringUtils::toStr(e.label) + "'.");
      }
    }

    // The design may list files that are intentionally absent from this ConsensusMap. QPX builds
    // its quantification units from the column headers, so remember that same set of cells for the
    // protein annotations. Otherwise every absent design file would add an invented zero-valued
    // quantity that the exporter correctly refuses as an extra key.
    quantification_fraction_group_labels_.clear();
    for (const auto& [map_index, header] : consensus.getColumnHeaders())
    {
      (void)map_index;
      const std::string filename = File::stemName(header.filename);
      const UInt label = header.getLabelAsUInt(consensus.getExperimentType());
      auto it = file_and_label_to_msfile_entry.find({filename, label});
      if (it == file_and_label_to_msfile_entry.end())
      {
        throw Exception::MissingInformation(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "File+Label referenced in consensus header not found in experimental design: "
          + filename + "\t" + StringUtils::toStr(label));
      }
      quantification_fraction_group_labels_.emplace(it->second.fraction_group, label);
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
        const std::string c_fn = File::stemName(h.filename); // filename according to experimental design in consensus map
        const UInt c_lab = h.getLabelAsUInt(consensus.getExperimentType());

        // find entry in experimental design (ignore extension and folder) that corresponds to current column header entry
        if (auto it = file_and_label_to_msfile_entry.find({c_fn, c_lab}); it != file_and_label_to_msfile_entry.end())
        {
          const size_t fraction = it->second.fraction;
          quantifyFeature_(f, fraction, c_fn, hit, c_lab); // updates "stats_.quant_features"
        }
        else
        {
          throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "File+Label referenced in consensus header not found in experimental design: " +
            c_fn + "\t" + StringUtils::toStr(c_lab));
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
    buildSampleIDLookup_(); // precompute (basename, channel) -> sample lookup

    stats_.n_samples = ed.getNumberOfSamples();
    stats_.n_fractions = ed.getNumberOfFractions();
    stats_.n_ms_files = ed.getNumberOfMSFiles();

    OPENMS_LOG_DEBUG << "Reading quant data: " << endl;
    OPENMS_LOG_DEBUG << "  MS files        : " << stats_.n_ms_files << endl;
    OPENMS_LOG_DEBUG << "  Fractions       : " << stats_.n_fractions << endl;
    OPENMS_LOG_DEBUG << "  Samples (Assays): " << stats_.n_samples << endl;

    stats_.total_features = peptides.size();
    
    countPeptides_(peptides);

    map<pair<std::string,Size>, std::string> identifier_idmergeidx_to_ms_file;
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
      // Origin file index within a merged ID run, annotated by IDMerger as
      // Constants::UserParam::ID_MERGE_INDEX ("id_merge_index"). Falls back to 0
      // for unmerged input (single primary MS run path).
      Size id_merge_idx = (Size)p.getMetaValue(Constants::UserParam::ID_MERGE_INDEX, 0);
      const PeptideHit& hit = p.getHits()[0];

      // don't quantify decoys
      if (hit.isDecoy()) continue;

      stats_.quant_features++;
      const AASequence& seq = hit.getSequence();
      const std::string& ms_file_path = identifier_idmergeidx_to_ms_file[{p.getIdentifier(),id_merge_idx}];

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
      std::string filename = File::stemName(ms_file_path);
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
    best_charge_by_peptidoform_.clear();
    prot_quant_.clear();
    design_cell_lookup_.clear();
    fraction_group_label_to_sample_.clear();
    quantification_fraction_group_labels_.clear();
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

    // number of labels/channels is constant across the design; compute once
    // (getNumberOfLabels() scans the whole MS file section on every call).
    const Size n_labels = experimental_design_.getNumberOfLabels();

    // Extract the Spectra Filepath column from the design
    map<UInt64, map<UInt64, std::string>> design_group_fraction_filename;
    for (ExperimentalDesign::MSFileSectionEntry const& f : msfile_section)
    {
      const std::string fn = File::stemName(f.path);
      design_group_fraction_filename[f.fraction_group][f.fraction] = fn;
    }

    auto & id_groups = proteins.getIndistinguishableProteins();

    // Build an accession -> protein-group index once. Previously the group for each
    // quantified protein was located with a linear std::find_if over all groups, making
    // annotation O(#quantified_proteins * #groups). For each accession we keep the FIRST
    // group (in id_groups order) that contains it, matching the previous find_if semantics.
    std::map<std::string, ProteinIdentification::ProteinGroup*> accession_to_group;
    for (auto & g : id_groups)
    {
      for (const auto & a : g.accessions)
      {
        accession_to_group.emplace(a, &g);
      }
    }

    for (const auto& q : protein_quants)
    {
      // accession of quantified protein(group)
      const std::string & acc = q.first;

      if (q.second.total_abundances.empty()) 
      {
        //TODO maybe just count the number of unquantifiable proteins and report that?
        OPENMS_LOG_DEBUG << "Protein " << acc << " not quantified." << endl;
        continue;
      } // not quantified
 
      // retrieve protein group with accession "acc" via the precomputed index
      auto group_it = accession_to_group.find(acc);

      if (group_it != accession_to_group.end())
      {
        ProteinIdentification::ProteinGroup* id_group = group_it->second;
        // copy abundances to float data array
        const SampleAbundances& total_abundances = q.second.total_abundances;
        const SampleAbundances& total_psm_counts = q.second.total_psm_counts;
        const SampleAbundances& total_distinct_peptides = q.second.total_distinct_peptides;
        const auto& file_level_psm_counts = q.second.file_level_psm_counts;
      
        // TODO: OPENMS_ASSERT(id_group->float_data_arrays.empty(), "Protein group float data array not empty!.");
        // Keep the existing positional arrays unchanged for mzTab and older consumers. QPX unit
        // quantities are appended and consumed by name, never by these indices.
        id_group->getFloatDataArrays().resize(5);
        id_group->getStringDataArrays().resize(2);
        id_group->getIntegerDataArrays().resize(4);
        
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

        // QPX 1.1 protein-group quantities are keyed by the exact OpenMS fraction_group and
        // label. Materialize every cell represented by the quantification input in deterministic
        // order, using zero for a quantified protein without sufficient evidence in that cell just
        // as the sample-level abundance array above does for an absent sample result.
        auto& fraction_group_level_abundance = id_group->getFloatDataArrays()[4];
        fraction_group_level_abundance.setName("fraction_group_level_abundance");
        auto& fraction_group_level_fraction_group = id_group->getIntegerDataArrays()[2];
        fraction_group_level_fraction_group.setName("fraction_group_level_fraction_group");
        auto& fraction_group_level_label = id_group->getIntegerDataArrays()[3];
        fraction_group_level_label.setName("fraction_group_level_label");
        for (const auto& [fraction_group, label] : quantification_fraction_group_labels_)
        {
          double abundance = 0.0;
          if (auto group_it = q.second.fraction_group_abundances.find(fraction_group);
              group_it != q.second.fraction_group_abundances.end())
          {
            if (auto label_it = group_it->second.find(label); label_it != group_it->second.end())
            {
              abundance = label_it->second;
            }
          }
          fraction_group_level_abundance.push_back(static_cast<float>(abundance));
          fraction_group_level_fraction_group.push_back(static_cast<Int>(fraction_group));
          fraction_group_level_label.push_back(static_cast<Int>(label));
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
            design_filename = File::stemName(design_filename);
            
            #ifdef DEBUG_PROTEINQUANTIFIER
            std::cout 
              << "Experimental design: fraction group: " << group_id 
              << ", filename: '" << design_filename
              << "', fraction: " << fraction
              << " of the experimental design." << std::endl;
            #endif

            // for each file in the design, fill the channels quantity
            for (Size c = 1; c <= n_labels; ++c) // label/channel numbers are 1-based
            {
              double channel_abundance{};

              const auto& filename_to_channel_map = q.second.channel_level_abundances;

              if (auto file_level_it = filename_to_channel_map.find(design_filename); 
                file_level_it != filename_to_channel_map.end())
              {
                if (file_level_it->second.contains(0)) throw Exception::MissingInformation(
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
    map<std::string, std::string> accession_to_leader = mapAccessionToLeader(proteins);

    bool contains_accessions{ false}; // flag to check if any accessions were found

    for (auto const& pep_q : pep_quant_)
    {
      std::string leader_accession = getAccession_(pep_q.second.accessions, accession_to_leader);
      OPENMS_LOG_DEBUG << "Peptide id mapped to leader: " << leader_accession << endl;

      // not enough evidence or mapping to multiple groups
      if (leader_accession.empty())
        continue;

      contains_accessions = true;
      // proteotypic peptide
      const std::string peptide = pep_q.first.toUnmodifiedString();

      prot_quant_[leader_accession].psm_count += pep_q.second.psm_count; // total PSM count for this group of proteins (represented by the leader accession)

      // transfer abundances and counts from peptides->protein
      // summarize abundances and counts between different peptidoforms
      for (auto const& sta : pep_q.second.total_abundances)
      {
        prot_quant_[leader_accession].peptide_abundances[peptide][sta.first] += sta.second;
      }

      for (const auto& [fraction_group, label_abundances] : pep_q.second.fraction_group_abundances)
      {
        for (const auto& [label, abundance] : label_abundances)
        {
          prot_quant_[leader_accession]
            .peptide_fraction_group_abundances[peptide][fraction_group][label] += abundance;
        }
      }

      for (auto const& sta : pep_q.second.total_psm_counts)
      {
        prot_quant_[leader_accession].peptide_psm_counts[peptide][sta.first] += sta.second;
      }

      // NOTE: "channel_level_abundances" is deliberately NOT filled here. It is the same
      // quantity as "total_abundances" at a finer granularity (per (file, channel) rather
      // than per sample, a 1:N relationship under fractionation) and is therefore owned
      // exclusively by calculateFileAndChannelLevelProteinAbundances_(), just like
      // "total_abundances" is owned by calculateProteinAbundances_(). Pre-filling it with
      // a raw sum over all peptides/peptidoforms/charges left that un-aggregated sum in
      // the output whenever the aggregation step later declined to quantify a cell (too
      // few peptides for 'top:N'), i.e. exactly where the value should have been absent.

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

  Size PeptideAndProteinQuant::countQuantifiedSamples_(const SampleAbundances& abundances)
  {
    Size n_quantified = 0;
    for (const auto& sample_abundance : abundances)
    {
      if (sample_abundance.second > 0.0)
      {
        ++n_quantified;
      }
    }
    return n_quantified;
  }

  std::vector<std::string> PeptideAndProteinQuant::selectPeptidesForQuantification_(const std::string& protein_accession,
                                                                              Size top_n,
                                                                              bool fix_peptides)
  {
    std::vector<std::string> peptides;
    
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
        if (countQuantifiedSamples_(ab.second) == stats_.n_samples)
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
                                                     const std::string& method) const
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

  void PeptideAndProteinQuant::calculateProteinAbundances_(const std::string& protein_accession,
                                                          const std::vector<std::string>& selected_peptides,
                                                          const std::string& aggregate_method,
                                                          Size top_n,
                                                          bool include_all)
  {
    auto prot_it = prot_quant_.find(protein_accession);
    if (prot_it == prot_quant_.end())
    {
      return;
    }

    ProteinData& pd = prot_it->second;

    // Consider only the selected peptides, and within them only the samples they were actually
    // detected in. An abundance stored as 0.0 is "not detected", not a measured absence:
    // IsobaricChannelExtractor writes a reporter it could not find (or one below
    // 'min_reporter_intensity') as 0.0 and inserts the handle regardless, so every peptide of an
    // isobaric run carries an entry for every channel. Aggregating those is what this tool's own
    // documentation says must not happen - "mean and median ignore missing cases, averaging only
    // present values" - and it is the artefact normalizePeptides_() excludes zeros to avoid: a
    // channel in which most of a protein's peptides went undetected gets a median of exactly 0,
    // i.e. the protein is reported absent from a sample it was measured in.
    // The samples are still tracked separately, so a sample whose peptides were all undetected
    // keeps its (zero) entry under 'include_all' rather than disappearing from the output.
    map<UInt64, DoubleList> abundances; // detected peptide abundances by sample
    set<UInt64> samples;                // every sample the selected peptides reach
    for (const auto& pep : selected_peptides)    // for all selected peptides
    {
      auto pep_it = pd.peptide_abundances.find(pep);
      if (pep_it != pd.peptide_abundances.end())
      {
        for (auto& sa : pep_it->second) // copy over the abundances that are measurements
        {
          samples.insert(sa.first);
          if (sa.second > 0.0) { abundances[sa.first].push_back(sa.second); }
        }
      }
    }

    for (UInt64 sample : samples)
    {
      DoubleList& values = abundances[sample]; // empty if nothing was detected in this sample

      // check if the protein has enough detected peptides in this sample
      if (!include_all && (top_n > 0) && (values.size() < top_n))
      {
        continue;
      }

      // if we have more than "top", reduce to the top ones
      if ((top_n > 0) && (values.size() > top_n))
      {
        // sort descending:
        sort(values.begin(), values.end(), greater<double>());
        values.resize(top_n); // remove all but best N values
      }

      double abundance_result = aggregateAbundances_(values, aggregate_method);
      pd.total_abundances[sample] = abundance_result;
    }
  }

  void PeptideAndProteinQuant::calculateFractionGroupLevelProteinAbundances_(
    const std::string& protein_accession,
    const std::vector<std::string>& selected_peptides,
    const std::string& aggregate_method,
    Size top_n,
    bool include_all)
  {
    auto prot_it = prot_quant_.find(protein_accession);
    if (prot_it == prot_quant_.end())
    {
      return;
    }

    ProteinData& pd = prot_it->second;
    std::map<std::pair<UInt, UInt>, DoubleList> abundances;
    for (const auto& peptide : selected_peptides)
    {
      auto peptide_it = pd.peptide_fraction_group_abundances.find(peptide);
      if (peptide_it == pd.peptide_fraction_group_abundances.end())
      {
        continue;
      }
      for (const auto& [fraction_group, label_abundances] : peptide_it->second)
      {
        for (const auto& [label, abundance] : label_abundances)
        {
          // Touch the cell even where the peptide was not detected, so a cell whose peptides were
          // all undetected still reports zero rather than dropping out. Only measurements become
          // values: an abundance stored as 0.0 is "not detected" - IsobaricChannelExtractor writes
          // a reporter it could not find (or one below 'min_reporter_intensity') as 0.0 and
          // inserts the handle anyway - so counting it here would let undetected peptides satisfy
          // the 'top_n' gate and drag a median or mean of the cell towards zero. Same rule as
          // calculateProteinAbundances_(), calculateFileAndChannelLevelProteinAbundances_() and
          // normalizePeptides_().
          DoubleList& cell = abundances[{fraction_group, label}];
          if (abundance > 0.0) { cell.push_back(abundance); }
        }
      }
    }

    for (auto& [unit_label, values] : abundances)
    {
      // Detected abundances only; empty means every peptide of this cell was undetected. That is
      // not the same as "no peptide reached this cell" - such a cell has no key at all - so it is
      // deliberately not skipped: the gate below drops it when 'top_n' peptides are required, and
      // aggregateAbundances_() renders it as 0.0 otherwise.
      if (!include_all && top_n > 0 && values.size() < top_n)
      {
        continue;
      }
      if (top_n > 0 && values.size() > top_n)
      {
        std::sort(values.begin(), values.end(), std::greater<double>());
        values.resize(top_n);
      }
      pd.fraction_group_abundances[unit_label.first][unit_label.second] =
        aggregateAbundances_(values, aggregate_method);
    }
  }

  void PeptideAndProteinQuant::calculateFileAndChannelLevelProteinAbundances_(const std::string& protein_accession,
                                                                  const std::vector<std::string>& selected_peptides,
                                                                  const std::string& aggregate_method,
                                                                  Size top_n,
                                                                  bool include_all,
                                                                  const std::map<std::string, std::string>& accession_to_leader,
                                                                  const UnmodifiedToEntriesIndex& unmod_to_entries)
  {
    auto prot_it = prot_quant_.find(protein_accession);
    if ( prot_it == prot_quant_.end())
    {
      return;
    }

    ProteinData& pd = prot_it->second;
    const bool select_best_charge = param_.getValue("best_charge_and_fraction") == "true";

    // Granularity note: a sample is a (file, label) pair of the experimental design, so in a
    // fractionated experiment ONE sample spans SEVERAL files (all fractions of its fraction
    // group at that label) and total_abundances aggregates over them. The cells computed here
    // are per (file, channel), i.e. one sample maps to N cells, and both the peptide selection
    // and the aggregation below are applied PER FILE. That is deliberate - these are per-file
    // values - but it means:
    //   - the cells only decompose the sample value for top_n == 0 with aggregate == "sum";
    //     top-N/median/mean/weighted_mean do not commute with aggregation across fractions;
    //   - the "at least top_n peptides" requirement is enforced per file, which is much
    //     stricter than the sample-level rule for fractionated data: a protein quantified at
    //     the sample level can have ALL of its cells dropped because no single fraction
    //     contains top_n peptides.
    // organize detailed abundances by (fraction, filename, channel) combinations
    map<tuple<Int, std::string, UInt>, DoubleList> channel_level_abundances_for_selected_peptides;
    
    // collect detailed abundances from selected peptides
    for (const auto& pep : selected_peptides)    // for all selected peptides
    {
      // Look up all modified peptidoforms sharing this unmodified sequence via the
      // precomputed index.
      auto idx_it = unmod_to_entries.find(pep);
      if (idx_it == unmod_to_entries.end()) { continue; }

      // Contribute exactly ONE value per selected peptide and (fraction, file, channel)
      // cell, summed over all matching peptidoforms. By default all charge states contribute;
      // with 'best_charge_and_fraction', each peptidoform contributes only the charge selected
      // before normalization. This mirrors its sample-level abundance and makes 'top_n' count
      // peptides rather than letting peptidoforms or charge states occupy separate top-N slots.
      map<tuple<Int, std::string, UInt>, double> abundances_of_peptide;

      for (const auto* pep_q_ptr : idx_it->second)
      {
        const auto& pep_q_check = *pep_q_ptr;
        std::string check_accession = getAccession_(pep_q_check.second.accessions, accession_to_leader);
        if (check_accession != protein_accession) { continue; } // peptidoform of another protein

        const Int* selected_charge = nullptr;
        if (select_best_charge)
        {
          auto selected_it = best_charge_by_peptidoform_.find(pep_q_check.first);
          if (selected_it == best_charge_by_peptidoform_.end()) { continue; }
          selected_charge = &selected_it->second;
        }

        for (auto const& fraction : pep_q_check.second.abundances)
        {
          for (auto const& filename : fraction.second)
          {
            for (auto const& charge : filename.second)
            {
              if (selected_charge != nullptr && charge.first != *selected_charge) { continue; }
              for (auto const& channel : charge.second)
              {
                abundances_of_peptide[make_tuple(fraction.first, filename.first, channel.first)] += channel.second;

                #ifdef DEBUG_PEPTIDEANDPROTEINQUANT
                std::cout << "DEBUG: Adding abundance for leader " << check_accession
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
      }

      for (auto const& [key, abundance] : abundances_of_peptide)
      {
        // Touch the cell even where the peptide was not detected, so a cell whose peptides were
        // all undetected still reports zero rather than dropping out of the parallel
        // (filename, channel, abundance) arrays. Only measurements become values: an abundance
        // stored as 0.0 is "not detected" - IsobaricChannelExtractor writes a reporter it could
        // not find (or one below 'min_reporter_intensity') as 0.0 and inserts the handle anyway -
        // so counting it here would let undetected peptides satisfy the per-file 'top_n' gate and
        // drag a median or mean of the cell towards zero. Same rule as
        // calculateProteinAbundances_() and normalizePeptides_().
        DoubleList& cell = channel_level_abundances_for_selected_peptides[key];
        if (abundance > 0.0) { cell.push_back(abundance); }
      }
    }
    
    // now aggregate using the same aggregation method
    for (auto& detailed_ab : channel_level_abundances_for_selected_peptides)
    {
      const auto& selected_peptide = detailed_ab.first;
      std::string filename = get<1>(selected_peptide);
      UInt channel = get<2>(selected_peptide);
      
      // Detected abundances only; empty means every peptide of this cell was undetected. That is
      // not the same as "no peptide reached this cell" - such a cell has no key at all - so it is
      // deliberately not skipped here: the gate below drops it when 'top_n' peptides are
      // required, and aggregateAbundances_() renders it as 0.0 otherwise.
      DoubleList& all_abundances = detailed_ab.second;

      // check if we have enough detected peptides for this detailed key
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
      auto leader_it = accession_to_leader.find(protein_accession);
      const std::string& leader_accession = (leader_it != accession_to_leader.end()) ? leader_it->second : protein_accession;
      std::cout << "DEBUG: Protein " << protein_accession
                << " leader " << leader_accession
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
      const std::string & hit_accession = hit.getAccession();
      const std::string & hit_sequence = hit.getSequence();

      if (prot_quant_.contains(hit_accession))
      {
        if (hit_sequence.empty())
        {
          prot_quant_.erase(hit_accession);
          OPENMS_LOG_WARN << "Removed " << hit_accession <<  ", no protein sequence found!" << endl;
        }
        else
        {
          std::vector<std::string_view> peptides {};
          digest.digestUnmodified(hit_sequence, peptides);
          const double n_theoretical_peptides = double(peptides.size());
          ProteinData& pd = prot_quant_[hit_accession];
          for (auto& total_abundance : pd.total_abundances)
          {
            total_abundance.second /= n_theoretical_peptides;
          }
          // the file+channel level abundances are the same quantity at a finer
          // granularity, so they need the same normalization
          for (auto& file_entry : pd.channel_level_abundances)
          {
            for (auto& channel_entry : file_entry.second)
            {
              channel_entry.second /= n_theoretical_peptides;
            }
          }
          for (auto& [fraction_group, label_abundances] :
               prot_quant_[hit_accession].fraction_group_abundances)
          {
            (void)fraction_group;
            for (auto& [label, abundance] : label_abundances)
            {
              (void)label;
              abundance /= double(peptides.size());
            }
          }
        }
      }
    }
  }
}
