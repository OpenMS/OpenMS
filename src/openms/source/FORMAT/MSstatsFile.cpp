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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Lukas Heumos $
// --------------------------------------------------------------------------

#include <include/OpenMS/ANALYSIS/ID/ConsensusMapMergerAlgorithm.h>
#include "OpenMS/FORMAT/MSstatsFile.h"

using namespace std;


OpenMS::MSstatsFile::MSstatsFile() = default;
OpenMS::MSstatsFile::~MSstatsFile() = default;

const OpenMS::String OpenMS::MSstatsFile::na_string_ = "NA";

void OpenMS::MSstatsFile::checkConditionLFQ_(const ExperimentalDesign::SampleSection& sampleSection,
                                             const String& bioreplicate,
                                             const String& condition)
{
  // Sample Section must contain the column that contains the condition used for MSstats
  if (!sampleSection.hasFactor(condition))
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Sample Section of the experimental design does not contain MSstats_Condition");
  }

  // Sample Section must contain column for the Bioreplicate
  if (!sampleSection.hasFactor(bioreplicate))
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Sample Section of the experimental design does not contain MSstats_BioReplicate");
  }
}

void OpenMS::MSstatsFile::checkConditionISO_(const ExperimentalDesign::SampleSection& sampleSection,
                                             const String& bioreplicate,
                                             const String& condition,
                                             const String& mixture)
{
  checkConditionLFQ_(sampleSection, bioreplicate, condition);

  // Sample Section must contain column for Mixture
  if (!sampleSection.hasFactor(mixture))
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Sample Section of the experimental design does not contain MSstats_Mixture");
  }
}


//TODO why do we need this method and store everything three times??? (Once in the CMap, once in the feature
// of aggregatedConsensusInfo, and once in the other fields of aggregatedConsensusInfo)
// Cant we just get this stuff on the fly?
// We go through the features anyway again.
OpenMS::MSstatsFile::AggregatedConsensusInfo OpenMS::MSstatsFile::aggregateInfo_(const ConsensusMap& consensus_map,
                                                                                 const std::vector<String>& spectra_paths)
{
  OpenMS::MSstatsFile::AggregatedConsensusInfo aggregatedInfo; //results
  const auto &column_headers = consensus_map.getColumnHeaders(); // needed for label_id

  for (const OpenMS::ConsensusFeature &consensus_feature : consensus_map)
  {

    vector<OpenMS::String> filenames;
    vector<OpenMS::MSstatsFile::Intensity> intensities;
    vector<OpenMS::MSstatsFile::Coordinate> retention_times;
    vector<unsigned> cf_labels;

    // Store the file names and the run intensities of this feature
    const OpenMS::ConsensusFeature::HandleSetType& fs(consensus_feature.getFeatures());
    for (const auto& feat : fs)
    {
      filenames.push_back(spectra_paths[feat.getMapIndex()]);
      intensities.push_back(feat.getIntensity());
      retention_times.push_back(feat.getRT());

      // Get the label_id from the file description MetaValue
      auto &column = column_headers.at(feat.getMapIndex());
      if (column.metaValueExists("channel_id"))
      {
        cf_labels.push_back(OpenMS::Int(column.getMetaValue("channel_id")));
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
template <class LineType>
void OpenMS::MSstatsFile::constructFile_(const String& retention_time_summarization_method,
                                         const bool rt_summarization_manual,
                                         TextFile& csv_out,
                                         const std::set<String>& peptideseq_quantifyable,
                                         LineType& peptideseq_to_prefix_to_intensities) const

{
  // sanity check that the triples (peptide_sequence, precursor_charge, run) only appears once
  set<tuple<OpenMS::String, OpenMS::String, OpenMS::String> > peptideseq_precursor_charge_run;

  int count_similar = 0;
  for (const auto &peptideseq : peptideseq_quantifyable)
  {
    for (const auto &line :
            peptideseq_to_prefix_to_intensities[peptideseq])
    {
      // First, we collect all retention times and intensities
      set<OpenMS::MSstatsFile::Coordinate> retention_times{};
      set<OpenMS::MSstatsFile::Intensity> intensities{};
      for (const auto &p : line.second)
      {
        if (retention_times.find(get<1>(p)) != retention_times.end())
        {
          OPENMS_LOG_WARN << "Peptide ion appears multiple times at the same retention time."
                             " This is not expected."
                          << endl;
        }
        else
        {
          retention_times.insert(get<1>(p));
          intensities.insert(get<0>(p));
        }
      }

      tuple<OpenMS::String, OpenMS::String, OpenMS::String> tpl = make_tuple(
          line.first.sequence(), line.first.precursor_charge(), line.first.run());

      if (peptideseq_precursor_charge_run.find(tpl) != peptideseq_precursor_charge_run.end())
      {
        //TODO What is this doing here??
        count_similar += 1;
      }
      peptideseq_precursor_charge_run.insert(tpl);

      // If the rt summarization method is set to manual, we simply output all it,rt pairs
      if (rt_summarization_manual)
      {
        for (const auto &ity_rt_file : line.second)
        {
          //RT, common prefix items, intensity, "unique ID (file+spectrumID)"
          csv_out.addLine(
              String(get<1>(ity_rt_file)) + ',' + line.first.toString() + ',' + String(get<0>(ity_rt_file)) + ','
              + quote_ + get<2>(ity_rt_file) + quote_);
        }
      }
      // Otherwise, the intensities are resolved over the retention times
      else
      {
        OpenMS::MSstatsFile::Intensity intensity(0);
        if (retention_time_summarization_method == "max")
        {
          intensity = *(max_element(intensities.begin(), intensities.end()));
        }
        else if (retention_time_summarization_method == "min")
        {
          intensity = *(min_element(intensities.begin(), intensities.end()));
        }
        else if (retention_time_summarization_method == "mean")
        {
          intensity = meanIntensity_(intensities);
        }
        else if (retention_time_summarization_method == "sum")
        {
          intensity = sumIntensity_(intensities);
        }
        //common prefix items, aggregated intensity, "unique ID (file of first spectrum in the set of 'same')"
        //@todo we could collect all spectrum references contributing to this intensity instead
        csv_out.addLine(
            line.first.toString() + delim_ + OpenMS::String(intensity) + delim_ + quote_ +
            get<2>(*line.second.begin()) + quote_);
      }
    }
  }
}

void OpenMS::MSstatsFile::storeLFQ(const String& filename,
                                   const ConsensusMap& consensus_map,
                                   const ExperimentalDesign& design,
                                   const StringList& reannotate_filenames,
                                   const bool is_isotope_label_type,
                                   const String& bioreplicate,
                                   const String& condition,
                                   const String& retention_time_summarization_method)
{
  // Experimental Design file
  const ExperimentalDesign::SampleSection& sampleSection = design.getSampleSection();

  if (design.getNumberOfLabels() != 1)
  {
     throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Too many labels for a label-free quantitation experiments. Please select the appropriate method, or validate the experimental desing.");
  }

  checkConditionLFQ_(sampleSection, bioreplicate, condition);

  // assemble lookup table for run (each combination of pathname and fraction is a run)
  std::map< pair< String, unsigned>, unsigned > run_map{};
  assembleRunMap_(run_map, design);

  // Maps run in MSstats input to run for OpenMS
  map< unsigned, unsigned > msstats_run_to_openms_fractiongroup;

  // Mapping of filepath and label to sample and fraction
  map< pair< String, unsigned >, unsigned> path_label_to_sample = design.getPathLabelToSampleMapping(true);
  map< pair< String, unsigned >, unsigned> path_label_to_fraction = design.getPathLabelToFractionMapping(true);
  map< pair< String, unsigned >, unsigned> path_label_to_fractiongroup = design.getPathLabelToFractionGroupMapping(true);

  // The Retention Time is additionally written to the output as soon as the user wants to resolve multiple peptides manually
  const bool rt_summarization_manual(retention_time_summarization_method == "manual");

  if (rt_summarization_manual)
  {
    OPENMS_LOG_WARN << "WARNING: rt_summarization set to manual."
                       " One feature might appear at multiple retention times in the output file."
                       " This is invalid input for standard MSstats."
                       " Combining of features over retention times is recommended!" << endl;
  }

  ExperimentalDesign::MSFileSection msfile_section = design.getMSFileSection();

  // Extract the Spectra Filepath column from the design
  std::vector<String> design_filenames{};
  for (ExperimentalDesign::MSFileSectionEntry const& f : msfile_section)
  {
    const String fn = File::basename(f.path);
    design_filenames.push_back(fn);
  }

  // Determine if the experiment has fractions
  const bool has_fraction = design.isFractionated();

  //vector< OpenMS::BaseFeature> features{};
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
  OpenMS::MSstatsFile::AggregatedConsensusInfo aggregatedInfo = OpenMS::MSstatsFile::aggregateInfo_(consensus_map, spectra_paths);

  // The output file of the MSstats converter
  TextFile csv_out;
  csv_out.addLine(
    String(rt_summarization_manual ? "RetentionTime,": "") +
    "ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,"
    "ProductCharge,IsotopeLabelType,Condition,BioReplicate,Run," +
    String(has_fraction ? "Fraction,": "") + "Intensity,Reference");

  // From the MSstats user guide: endogenous peptides (use "L") or labeled reference peptides (use "H").
  String isotope_label_type = "L";
  if (is_isotope_label_type) //@todo remove? not sure if this is correct. I think DDA LFQ is always "L"
  {
    // use the channel_id information (?)
    isotope_label_type = "H";
  }

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


  // Stores all the lines that will be present in the final MSstats output
  // Several things needs to be considered:
  // - We need to map peptide sequences to full features, because then we can ignore peptides
  //   that are mapped to multiple proteins.
  // - We also need to map to the intensities, such that we combine intensities over multiple retention times.
  map< String, map< MSstatsLine_, set< tuple<Intensity, Coordinate, String> > > > peptideseq_to_prefix_to_intensities;

  for (Size i = 0; i < aggregatedInfo.features.size(); ++i)
  {
    const OpenMS::BaseFeature &base_feature = aggregatedInfo.features[i];

    for (const OpenMS::PeptideIdentification &pep_id : base_feature.getPeptideIdentifications())
    {
      for (const OpenMS::PeptideHit & pep_hit : pep_id.getHits())
      {
        //TODO Really double check with Meena Choi (MSStats author) or make it an option! I can't find any info
        // on what is correct. For TMT we include them (since it is necessary) (see occurrence above as well when map is built!)
        const String & sequence = pep_hit.getSequence().toString(); // to modified string

        // check if all referenced protein accessions are part of the same indistinguishable group
        // if so, we mark the sequence as quantifiable
        std::set<String> accs = pep_hit.extractProteinAccessionsSet();

        //Note: In general as long as we only support merged proteins across conditions,
        // we check if the map is already set at this sequence since
        // it cannot happen that to peptides with the same sequence map to different proteins unless something is wrong.
        // Also I think MSstats cannot handle different associations to proteins across conditions.
        if (isQuantifyable_(accs, accession_to_group))
        {
          peptideseq_quantifyable.emplace(sequence);
        }
        else
        {
          continue; // we dont need the rest of the loop
        }

        // Variables of the peptide hit
        // MSstats User manual 3.7.3: Unknown precursor charge should be set to 0
        const Int precursor_charge = pep_hit.getCharge();

        // Unused for DDA data anyway
        String fragment_ion = na_string_;
        String frag_charge = "0";

        String accession  = ListUtils::concatenate(accs,accdelim_);
        if (accession.empty()) accession = na_string_; //shouldn't really matter since we skip unquantifyable peptides

        // Write new line for each run
        for (Size j = 0; j < aggregatedInfo.consensus_feature_filenames[i].size(); j++)
        {
          const String &current_filename = aggregatedInfo.consensus_feature_filenames[i][j];
          const Intensity intensity(aggregatedInfo.consensus_feature_intensities[i][j]);
          const Coordinate retention_time(aggregatedInfo.consensus_feature_retention_times[i][j]);
          const unsigned label(aggregatedInfo.consensus_feature_labels[i][j]);

          const pair< String, unsigned> tpl1 = make_pair(current_filename, label);
          const unsigned sample = path_label_to_sample[tpl1];
          const unsigned fraction = path_label_to_fraction[tpl1];

          const pair< String, unsigned> tpl2 = make_pair(current_filename, fraction);

          // Resolve run
          const unsigned run = run_map[tpl2];  // MSstats run according to the file table
          const unsigned openms_fractiongroup = path_label_to_fractiongroup[tpl1];
          msstats_run_to_openms_fractiongroup[run] = openms_fractiongroup;

          // Assemble MSstats line
          //TODO since a lot of cols are constant in DDA LFQ, we could reduce the prefix and add the constant
          // cols on-the-fly during constructFile_ (so we save during checking duplicates)
          MSstatsLine_ prefix(
                  has_fraction,
                  accession,
                  sequence,
                  precursor_charge,
                  fragment_ion,
                  frag_charge,
                  isotope_label_type,
                  sampleSection.getFactorValue(sample, condition),
                  sampleSection.getFactorValue(sample, bioreplicate),
                  String(run),
                  (has_fraction ? String(fraction) : "")
          );
          tuple<Intensity, Coordinate, String> intensity_retention_time = make_tuple(intensity, retention_time, current_filename);
          peptideseq_to_prefix_to_intensities[sequence][prefix].insert(intensity_retention_time);
        }
      }
    }
  }

  // Print the run mapping between MSstats and OpenMS
  for (const auto& run_mapping : msstats_run_to_openms_fractiongroup)
  {
    cout << "MSstats run " << String(run_mapping.first)
         << " corresponds to OpenMS fraction group " << String(run_mapping.second) << endl;
  }

  constructFile_(retention_time_summarization_method,
                 rt_summarization_manual,
                 csv_out,
                 peptideseq_quantifyable,
                 peptideseq_to_prefix_to_intensities);

  // Store the final assembled CSV file
  csv_out.store(filename);
}

void OpenMS::MSstatsFile::storeISO(const String& filename,
                                   const ConsensusMap& consensus_map,
                                   const ExperimentalDesign& design,
                                   const StringList& reannotate_filenames,
                                   const String& bioreplicate,
                                   const String& condition,
                                   const String& mixture,
                                   const String& retention_time_summarization_method)
{
  // Experimental Design file
  const ExperimentalDesign::SampleSection& sampleSection = design.getSampleSection();

  checkConditionISO_(sampleSection, bioreplicate, condition, mixture);

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

  // Maps run in MSstats input to run for OpenMS
  map< unsigned, unsigned > msstats_run_to_openms_fractiongroup;

  // Mapping of filepath and label to sample and fraction
  map< pair< String, unsigned >, unsigned> path_label_to_sample = design.getPathLabelToSampleMapping(true);
  map< pair< String, unsigned >, unsigned> path_label_to_fraction = design.getPathLabelToFractionMapping(true);
  map< pair< String, unsigned >, unsigned> path_label_to_fractiongroup = design.getPathLabelToFractionGroupMapping(true);

  // The Retention Time is additionally written to the output as soon as the user wants to resolve multiple peptides manually
  bool rt_summarization_manual(retention_time_summarization_method == "manual");

  if (!rt_summarization_manual)
  {
    OPENMS_LOG_WARN << "WARNING: rt_summarization set to something else than 'manual' but MSstatsTMT does aggregation of"
            " intensities of peptide-chargestate combinations in the same file itself."
            " Reverting to 'manual'" << endl;
    rt_summarization_manual = true;
  }

  ExperimentalDesign::MSFileSection msfile_section = design.getMSFileSection();

  // Extract the Spectra Filepath column from the design
  std::vector<String> design_filenames;
  for (ExperimentalDesign::MSFileSectionEntry const& f : msfile_section)
  {
    const String fn = File::basename(f.path);
    design_filenames.push_back(fn);
  }

  vector< OpenMS::BaseFeature> features;
  vector< String > spectra_paths;

  features.reserve(consensus_map.size());

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
  OpenMS::MSstatsFile::AggregatedConsensusInfo AggregatedInfo = OpenMS::MSstatsFile::aggregateInfo_(consensus_map, spectra_paths);

  // The output file of the MSstatsConverter
  TextFile csv_out;
  csv_out.addLine(String(rt_summarization_manual ? "RetentionTime,": "") +
    "ProteinName,PeptideSequence,Charge,Channel,Condition,BioReplicate,Run,Mixture,TechRepMixture,Fraction,Intensity,Reference");

  // We quantify indistinguishable groups with one (corner case) or multiple proteins.
  // If indistinguishable groups are not annotated (no inference or only trivial inference has been performed) we assume
  // that all proteins can be independently quantified (each forming an indistinguishable group).
  //TODO refactor since shared with LFQ and ISO
  const IndProtGrps& ind_prots = consensus_map.getProteinIdentifications()[0].getIndistinguishableProteins();

  // Map protein accession to its indistinguishable group
  std::unordered_map< String, const IndProtGrp* > accession_to_group = getAccessionToGroupMap_(ind_prots);

  std::set<String> peptideseq_quantifyable; //set for deterministic ordering

  // Stores all the lines that will be present in the final MSstats output,
  // We need to map peptide sequences to full features, because then we can ignore peptides
  // that are mapped to multiple proteins. We also need to map to the
  // intensities, such that we combine intensities over multiple retention times.
  map< String, map< MSstatsTMTLine_, set< tuple<Intensity, Coordinate, String> > > > peptideseq_to_prefix_to_intensities;

  for (Size i = 0; i < AggregatedInfo.features.size(); ++i)
  {
    const OpenMS::BaseFeature &base_feature = AggregatedInfo.features[i];

    for (const OpenMS::PeptideIdentification &pep_id : base_feature.getPeptideIdentifications())
    {
      String nativeID = "NONATIVEID";
      if (pep_id.metaValueExists("spectrum_reference"))
      {
        nativeID = pep_id.getMetaValue("spectrum_reference");
      }

      for (const OpenMS::PeptideHit & pep_hit : pep_id.getHits())
      {
        // Variables of the peptide hit
        // MSstats User manual 3.7.3: Unknown precursor charge should be set to 0
        const Int precursor_charge = (std::max)(pep_hit.getCharge(), 0);
        const String & sequence = pep_hit.getSequence().toString();

        // check if all referenced protein accessions are part of the same indistinguishable group
        // if so, we mark the sequence as quantifiable
        std::set<String> accs = pep_hit.extractProteinAccessionsSet();

        // When using extractProteinAccessionSet, we do not really need to loop over Evidences
        // anymore since MSStats does not care about anything else but the Protein accessions

        if (isQuantifyable_(accs, accession_to_group))
        {
          peptideseq_quantifyable.emplace(sequence);
        }
        else
        {
          continue; // we dont need the rest of the loop
        }

        String accession = ListUtils::concatenate(accs,accdelim_);
        if (accession.empty()) accession = na_string_; //shouldn't really matter since we skip unquantifyable peptides

        // Write new line for each run
        for (Size j = 0; j < AggregatedInfo.consensus_feature_filenames[i].size(); j++)
        {
          const String &current_filename = AggregatedInfo.consensus_feature_filenames[i][j];

          const Intensity intensity(AggregatedInfo.consensus_feature_intensities[i][j]);
          const Coordinate retention_time(AggregatedInfo.consensus_feature_retention_times[i][j]);
          const unsigned channel(AggregatedInfo.consensus_feature_labels[i][j] + 1);

          const pair< String, unsigned> tpl1 = make_pair(current_filename, channel);
          const unsigned sample = path_label_to_sample[tpl1];
          const unsigned fraction = path_label_to_fraction[tpl1];

          // Resolve techrepmixture, run
          const unsigned openms_fractiongroup = path_label_to_fractiongroup[tpl1];
          String techrepmixture = String(sampleSection.getFactorValue(sample, mixture)) + "_" + String(openms_fractiongroup);
          String run = techrepmixture + "_" + String(fraction);

          // Assemble MSstats line
          MSstatsTMTLine_ prefix(
              accession,
              sequence,
              precursor_charge,
              String(channel),
              sampleSection.getFactorValue(sample, condition),
              sampleSection.getFactorValue(sample, bioreplicate),
              String(run),
              sampleSection.getFactorValue(sample, mixture),
              String(techrepmixture),
              String(fraction)
          );

          String identifier = current_filename;
          if (rt_summarization_manual)
          {
            identifier += "_" + nativeID;
          }
          tuple<Intensity, Coordinate, String> intensity_retention_time = make_tuple(intensity, retention_time, identifier);
          peptideseq_to_prefix_to_intensities[sequence][prefix].insert(intensity_retention_time);
        }
      }
    }
  }

  // Print the run mapping between MSstats and OpenMS
  for (const auto& run_mapping : msstats_run_to_openms_fractiongroup)
  {
    cout << "MSstats run " << String(run_mapping.first)
         << " corresponds to OpenMS TechRepMixture " << String(run_mapping.second) << endl;
  }

  constructFile_(retention_time_summarization_method,
                 rt_summarization_manual,
                 csv_out,
                 peptideseq_quantifyable,
                 peptideseq_to_prefix_to_intensities);

  // Store the final assembled CSV file
  csv_out.store(filename);
}

bool OpenMS::MSstatsFile::checkUnorderedContent_(const std::vector<String> &first, const std::vector<String> &second)
{
  const std::set< String > lhs(first.begin(), first.end());
  const std::set< String > rhs(second.begin(), second.end());
  return lhs == rhs
         && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

void OpenMS::MSstatsFile::assembleRunMap_(
    std::map< std::pair<OpenMS::String, unsigned>, unsigned> &run_map,
    const OpenMS::ExperimentalDesign &design)
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

std::unordered_map<OpenMS::String, const OpenMS::IndProtGrp* > OpenMS::MSstatsFile::getAccessionToGroupMap_(const OpenMS::IndProtGrps& ind_prots)
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

bool OpenMS::MSstatsFile::isQuantifyable_(
    const std::set<OpenMS::String>& accs,
    const std::unordered_map<OpenMS::String, const OpenMS::IndProtGrp*>& accession_to_group) const
{
  bool quantifyable = true;

  if (accs.empty())
  {
    quantifyable = false;
  }
  else if (accs.size() > 1)
      // every prot accession in the set needs to belong to the same indist. group to make this peptide
      // eligible for quantification
  {
    std::set<const IndProtGrp*> maps_to_indgrps;

    const IndProtGrp* grp = nullptr;
    auto git = accession_to_group.find(*accs.begin());
    if (git != accession_to_group.end()) grp = git->second;

    auto accit = ++accs.begin();
    for (; accit != accs.end(); ++accit)
    {
      const auto it = accession_to_group.find(*accit);
      if (it != accession_to_group.end())
      {
        maps_to_indgrps.insert(it->second);
        if (it->second == grp)
        {
          continue;
        }
        else
        {
          // two different groups
          quantifyable = false;
          break;
        }
      }
      else
      {
        // we assume that it is a singleton. Cannot be quantifiable anymore.
        // Set makes them unique. Non-membership in groups means that there is at least one other
        // non-agreeing protein in the set.
        quantifyable = false;
        break;
      }
    }
  }
  return quantifyable;
}
