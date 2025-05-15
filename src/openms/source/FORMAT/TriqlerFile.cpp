// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Lukas Heumos $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/TriqlerFile.h>

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>

#include <tuple>

using namespace std;
using namespace OpenMS;


const String TriqlerFile::na_string_ = "NA";

void TriqlerFile::checkConditionLFQ_(const ExperimentalDesign::SampleSection& sampleSection,
                                             const String& condition)
{
  // Sample Section must contain the column that contains the condition used for Triqler
  if (!sampleSection.hasFactor(condition))
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Sample Section of the experimental design does not contain condition column: " + condition);
  }
}

//TODO why do we need this method and store everything three times??? (Once in the CMap, once in the feature
// of aggregatedConsensusInfo, and once in the other fields of aggregatedConsensusInfo)
// Cant we just get this stuff on the fly?
// We go through the features anyway again.
TriqlerFile::AggregatedConsensusInfo TriqlerFile::aggregateInfo_(const ConsensusMap& consensus_map,
                                                                                 const std::vector<String>& spectra_paths)
{
  TriqlerFile::AggregatedConsensusInfo aggregatedInfo; //results
  const auto &column_headers = consensus_map.getColumnHeaders(); // needed for label_id

  for (const ConsensusFeature &consensus_feature : consensus_map)
  {
    vector<String> filenames;
    vector<TriqlerFile::Intensity> intensities;
    vector<TriqlerFile::Coordinate> retention_times;
    vector<unsigned> cf_labels;

    // Store the file names and the run intensities of this feature
    const ConsensusFeature::HandleSetType& fs(consensus_feature.getFeatures());
    for (const auto& feat : fs)
    {
      filenames.push_back(spectra_paths[feat.getMapIndex()]);
      intensities.push_back(feat.getIntensity());
      retention_times.push_back(feat.getRT());

      // Get the label_id from the file description MetaValue
      auto &column = column_headers.at(feat.getMapIndex());
      if (column.metaValueExists("channel_id"))
      {
        cf_labels.push_back(Int(column.getMetaValue("channel_id")));
      }
      else
      {
        // label id 1 is used in case the experimental design specifies a LFQ experiment
        //TODO Not really, according to the if-case it only cares about the metavalue.
        // which could be missing due to other reasons
        cf_labels.push_back(1u);
      }
    }
    aggregatedInfo.consensus_feature_labels.push_back(cf_labels);
    aggregatedInfo.consensus_feature_filenames.push_back(filenames);
    aggregatedInfo.consensus_feature_intensities.push_back(intensities);
    aggregatedInfo.consensus_feature_retention_times.push_back(retention_times);
    aggregatedInfo.features.push_back(consensus_feature);
  }
  return aggregatedInfo;
}

//@todo LineType should be a template only for the line, not for the whole
// mapping structure. More exact type matching/info then.
void TriqlerFile::constructFile_(TextFile& csv_out,
                                         const std::set<String>& peptideseq_quantifyable,
                                         const MapSequenceToLines_& peptideseq_to_line) const

{
  for (const auto& peptideseq : peptideseq_quantifyable)
  {
    for (const auto& line : peptideseq_to_line.at(peptideseq)) { csv_out.addLine(line.toString()); } 
  }
}

void TriqlerFile::storeLFQ(const String& filename,
                                   const ConsensusMap& consensus_map,
                                   const ExperimentalDesign& design,
                                   const StringList& reannotate_filenames,
                                   const String& condition)
{
  // Experimental Design file
  const ExperimentalDesign::SampleSection& sampleSection = design.getSampleSection();

  if (design.getNumberOfLabels() != 1)
  {
     throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Too many labels for a label-free quantitation experiments. Please select the appropriate method, or validate the experimental design.");
  }

  checkConditionLFQ_(sampleSection, condition);

  // assemble lookup table for run (each combination of pathname and fraction is a run)
  std::map< pair< String, unsigned>, unsigned > run_map{};
  assembleRunMap_(run_map, design);

  // Maps run in Triqler input to run for OpenMS
  map< unsigned, unsigned > Triqler_run_to_openms_fractiongroup;

  // Mapping of filepath and label to sample and fraction
  map< pair< String, unsigned >, unsigned> path_label_to_sample = design.getPathLabelToSampleMapping(true);
  map< pair< String, unsigned >, unsigned> path_label_to_fraction = design.getPathLabelToFractionMapping(true);
  map< pair< String, unsigned >, unsigned> path_label_to_fractiongroup = design.getPathLabelToFractionGroupMapping(true);

  ExperimentalDesign::MSFileSection msfile_section = design.getMSFileSection();

  // Extract the Spectra Filepath column from the design
  std::vector<String> design_filenames{};
  for (ExperimentalDesign::MSFileSectionEntry const& f : msfile_section)
  {
    const String fn = File::basename(f.path);
    design_filenames.push_back(fn);
  }

  //vector< BaseFeature> features{};
  vector< String > spectra_paths{};

  //features.reserve(consensus_map.size());

  if (reannotate_filenames.empty())
  {
    consensus_map.getPrimaryMSRunPath(spectra_paths);
  }
  else
  {
    spectra_paths = reannotate_filenames;
  }

  // Reduce spectra path to the basename of the files
  for (Size i = 0; i < spectra_paths.size(); ++i)
  {
    spectra_paths[i] = File::basename(spectra_paths[i]);
  }

  if (!checkUnorderedContent_(spectra_paths, design_filenames))
  {
    OPENMS_LOG_FATAL_ERROR << "The filenames (extension ignored) in the consensusXML file are not the same as in the experimental design" << endl;
    OPENMS_LOG_FATAL_ERROR << "Spectra files (consensus map): \n";
    for (auto const & s : spectra_paths)
    {
      OPENMS_LOG_FATAL_ERROR << s << endl;
    }
    OPENMS_LOG_FATAL_ERROR << "Spectra files (design): \n";
    for (auto const & s : design_filenames)
    {
      OPENMS_LOG_FATAL_ERROR << s << endl;
    }
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The filenames (extension ignored) in the consensusXML file are not the same as in the experimental design");
  }

  // Extract information from the consensus features.
  TriqlerFile::AggregatedConsensusInfo aggregatedInfo = TriqlerFile::aggregateInfo_(consensus_map, spectra_paths);

  // The output file of the Triqler converter
  TextFile csv_out;

  // output the header
  csv_out.addLine("run\tcondition\tcharge\tsearchScore\tintensity\tpeptide\tproteins");

  if (consensus_map.getProteinIdentifications().empty())
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
	  "No protein information found in the ConsensusXML.");
  }

  // warn if we have more than one protein ID run
  //TODO actually allow having more than one inference run e.g. for different conditions
  if (consensus_map.getProteinIdentifications().size() > 1)
  {
    OPENMS_LOG_WARN << "Found " +
    String(consensus_map.getProteinIdentifications().size()) +
    " protein runs in consensusXML. Using first one only to parse inference data for now." << std::endl;
  }

  if (!consensus_map.getProteinIdentifications()[0].hasInferenceData())
  {
    OPENMS_LOG_WARN << "No inference was performed on the first run, defaulting to one-peptide-rule." << std::endl;
  }

  // We quantify indistinguishable groups with one (corner case) or multiple proteins.
  // If indistinguishable groups are not annotated (no inference or only trivial inference has been performed) we assume
  // that all proteins can be independently quantified (each forming an indistinguishable group).
  // TODO currently we always create the mapping. If groups are missing we create it based on singletons which is
  //  quite unnecessary. Think about skipping if no groups are present

  //consensus_map.getProteinIdentifications()[0].fillIndistinguishableGroupsWithSingletons();
  const IndProtGrps& ind_prots = consensus_map.getProteinIdentifications()[0].getIndistinguishableProteins();

  // Map protein accession to its indistinguishable group
  std::unordered_map< String, const IndProtGrp* > accession_to_group = getAccessionToGroupMap_(ind_prots);

  // To aggregate/uniquify on peptide sequence-level and save if a peptide is quantifyable
  std::set<String> peptideseq_quantifyable; //set for deterministic ordering

  // Stores all the lines that will be present in the final Triqler output
  MapSequenceToLines_ peptideseq_to_line;
  IDScoreSwitcherAlgorithm scores;
  for (Size i = 0; i < aggregatedInfo.features.size(); ++i)
  {
    const BaseFeature &base_feature = aggregatedInfo.features[i];

    for (const PeptideIdentification &pep_id : base_feature.getPeptideIdentifications())
    {
      if (!scores.isScoreType(pep_id.getScoreType(), IDScoreSwitcherAlgorithm::ScoreType::PEP))
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "TriqlerFile export expects Posterior Error Probabilities in the IDs of all features"
          " to convert them to Posterior Probabilities.");
      }
      for (const PeptideHit & pep_hit : pep_id.getHits())
      {
        const String & sequence = pep_hit.getSequence().toString(); // to modified string

        const double & search_score = pep_hit.getScore();

        // check if all referenced protein accessions are part of the same indistinguishable group
        // if so, we mark the sequence as quantifiable
        std::set<String> accs = pep_hit.extractProteinAccessionsSet();

        //Note: In general as long as we only support merged proteins across conditions,
        // we check if the map is already set at this sequence since
        // it cannot happen that to peptides with the same sequence map to different proteins unless something is wrong.
        // Also I think Triqler cannot handle different associations to proteins across conditions.
        if (isQuantifyable_(accs, accession_to_group))
        {
          peptideseq_quantifyable.emplace(sequence);
        }
        else
        {
          continue; // we don't need the rest of the loop
        }

        // Variables of the peptide hit
        const Int precursor_charge = pep_hit.getCharge();

        String accession  = ListUtils::concatenate(accs, accdelim_);
        if (accession.empty())
        {
          accession = na_string_; // shouldn't really matter since we skip unquantifiable peptides
        }
        // Write new line for each run
        for (Size j = 0; j < aggregatedInfo.consensus_feature_filenames[i].size(); j++)
        {
          const String &current_filename = aggregatedInfo.consensus_feature_filenames[i][j];
          const Intensity intensity(aggregatedInfo.consensus_feature_intensities[i][j]);
          // const Coordinate retention_time(aggregatedInfo.consensus_feature_retention_times[i][j]);
          const unsigned label(aggregatedInfo.consensus_feature_labels[i][j]);

          const pair< String, unsigned> tpl1 = make_pair(current_filename, label);
          const unsigned sample = path_label_to_sample[tpl1];
          const unsigned fraction = path_label_to_fraction[tpl1];

          const pair< String, unsigned> tpl2 = make_pair(current_filename, fraction);

          // Resolve run
          const unsigned run = run_map[tpl2];  // Triqler run according to the file table
          const unsigned openms_fractiongroup = path_label_to_fractiongroup[tpl1];
          Triqler_run_to_openms_fractiongroup[run] = openms_fractiongroup;

          // Store Triqler line          
          peptideseq_to_line[sequence].insert(TriqlerLine_(
                  String(run),
                  sampleSection.getFactorValue(sample, condition),
                  String(precursor_charge),
                  String(1. - search_score),
                  String(intensity), 
                  sequence,
                  accession));
        }
      }
    }
  }

  // Print the run mapping between Triqler and OpenMS
  for (const auto& run_mapping : Triqler_run_to_openms_fractiongroup)
  {
    cout << "Triqler run " << String(run_mapping.first)
         << " corresponds to OpenMS fraction group " << String(run_mapping.second) << endl;
  }

  constructFile_(csv_out,
                 peptideseq_quantifyable,
                 peptideseq_to_line);

  // Store the final assembled CSV file
  csv_out.store(filename);
}

bool TriqlerFile::checkUnorderedContent_(const std::vector<String> &first, const std::vector<String> &second)
{
  const std::set< String > lhs(first.begin(), first.end());
  const std::set< String > rhs(second.begin(), second.end());
  return lhs == rhs
         && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

void TriqlerFile::assembleRunMap_(
    std::map< std::pair<String, unsigned>, unsigned> &run_map,
    const ExperimentalDesign &design)
{
  run_map.clear();
  const ExperimentalDesign::MSFileSection& msfile_section = design.getMSFileSection();
  unsigned run_counter = 1;

  for (ExperimentalDesign::MSFileSectionEntry const& r : msfile_section)
  {
    std::pair< String, unsigned> tpl = std::make_pair(File::basename(r.path), r.fraction);
    if (run_map.find(tpl) == run_map.end())
    {
      run_map[tpl] = run_counter++;
    }
  }
}

std::unordered_map<String, const IndProtGrp* > TriqlerFile::getAccessionToGroupMap_(const IndProtGrps& ind_prots)
{
  std::unordered_map<String, const IndProtGrp* > res{};
  for (const IndProtGrp& pgrp : ind_prots)
  {
    for (const String& a : pgrp.accessions)
    {
      res[a] = &(pgrp);
    }
  }
  return res;
}

bool TriqlerFile::isQuantifyable_(
    const std::set<String>& accs,
    const std::unordered_map<String, const IndProtGrp*>& accession_to_group) const
{
  if (accs.empty())
  {
    return false;
  }
  if (accs.size() == 1)
  {
    return true;
  }
  auto git = accession_to_group.find(*accs.begin());
  if (git == accession_to_group.end())
  {
    return false;
  }
  const IndProtGrp* grp = git->second;

  // every prot accession in the set needs to belong to the same indist. group to make this peptide
  // eligible for quantification
  auto accit = ++accs.begin();
  for (; accit != accs.end(); ++accit)
  {
    const auto it = accession_to_group.find(*accit);

    // we assume that it is a singleton. Cannot be quantifiable anymore.
    // Set makes them unique. Non-membership in groups means that there is at least one other
    // non-agreeing protein in the set.
    if (it == accession_to_group.end())
    {
      return false;
    }
    // check if two different groups
    if (it->second != grp)
    {
      return false;
    }
  }
  
  return true;
}

String TriqlerFile::TriqlerLine_::toString() const
{
  const String delim("\t");
  return  run_
          + delim + condition_
          + delim + precursor_charge_
          + delim + search_score_
          + delim + intensity_
          + delim + sequence_
          + delim + accession_;
}
