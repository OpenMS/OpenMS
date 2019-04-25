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

#include "OpenMS/FORMAT/MSstatsFile.h"

using namespace std;

OpenMS::MSstatsFile::MSstatsFile()
{

}

OpenMS::MSstatsFile::~MSstatsFile()
{

}

void OpenMS::MSstatsFile::checkConditionLFQ_(const ExperimentalDesign::SampleSection& sampleSection, const String& bioreplicate, const String& condition)
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

void OpenMS::MSstatsFile::checkConditionISO_(const ExperimentalDesign::SampleSection sampleSection, const String& bioreplicate, const String& condition, const String& mixture)
{
  checkConditionLFQ_(sampleSection, bioreplicate, condition);
  
  // Sample Section must contain column for Mixture
  if (!sampleSection.hasFactor(mixture))
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Sample Section of the experimental design does not contain MSstats_Mixture");
  } 
}

void OpenMS::MSstatsFile::storeLFQ(const OpenMS::String &filename, const ConsensusMap &consensus_map,
                                const OpenMS::ExperimentalDesign& design, const StringList& reannotate_filenames,
                                const bool is_isotope_label_type, const String& bioreplicate, const String& condition,
                                const String& retention_time_summarization_method)
{
  // Experimental Design file
  ExperimentalDesign::SampleSection sampleSection = design.getSampleSection();

  if (design.getNumberOfLabels() != 1)
  {
     throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Too many lables for a label-free quantitation experiments. Please select the appropriate method, or validate the experimental desing.");
  }

  checkConditionLFQ_(sampleSection, bioreplicate, condition);

  // assemble lookup table for run (each combination of pathname and fraction is a run)
  std::map< pair< String, unsigned>, unsigned > run_map;
  assembleRunMap(run_map, design);

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
    cout << "WARNING: One feature might appear at multiple retention times in the output file. This is invalid input for MSstats. Combining of features over retention times is needed!" << endl;
  }

  typedef OpenMS::Peak2D::IntensityType Intensity;
  typedef OpenMS::Peak2D::CoordinateType Coordinate;

  ExperimentalDesign::MSFileSection msfile_section = design.getMSFileSection();

  // Extract the Spectra Filepath column from the design
  std::vector<String> design_filenames;
  for (ExperimentalDesign::MSFileSectionEntry const& f : msfile_section)
  {
    const String fn = File::basename(f.path);
    design_filenames.push_back(fn);
  }

  // Determine if the experiment has fractions
  const bool has_fraction = design.isFractionated();

  // label id 1 is used in case the experimental design specifies a LFQ experiment
  const unsigned label_lfq(1);

  vector< OpenMS::BaseFeature> features;
  vector< String > spectra_paths;

  // For each ConsensusFeature, store several attributes
  vector< vector< String > > consensus_feature_filenames;           // Filenames of ConsensusFeature
  vector< vector< Intensity > > consensus_feature_intensites;       // Intensites of ConsensusFeature
  vector< vector< Coordinate > > consensus_feature_retention_times; // Retention times of ConsensusFeature
  vector< vector< unsigned > > consensus_feature_labels;          // Labels of ConsensusFeature

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
    LOG_FATAL_ERROR << "The filenames (extension ignored) in the consensusXML file are not the same as in the experimental design" << endl;
    LOG_FATAL_ERROR << "Spectra files (consensus map): \n";
    for (auto const & s : spectra_paths)
    {
      LOG_FATAL_ERROR << s << endl;
    }
    LOG_FATAL_ERROR << "Spectra files (design): \n";
    for (auto const & s : design_filenames)
    {
      LOG_FATAL_ERROR << s << endl;
    }
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The filenames (extension ignored) in the consensusXML file are not the same as in the experimental design");
  };

  // Extract information from the consensus features.
  for (const ConsensusFeature &consensus_feature : consensus_map)
  {
    features.push_back(consensus_feature);

    vector< String > filenames;
    vector< Intensity > intensities;
    vector< Coordinate > retention_times;
    vector< unsigned > cf_labels;

    // Store the file names and the run intensities of this feature
    const ConsensusFeature::HandleSetType fs(consensus_feature.getFeatures());
    for (ConsensusFeature::HandleSetType::const_iterator fit = fs.begin(); fit != fs.end(); ++fit)
    {
      filenames.push_back(spectra_paths[fit->getMapIndex()]);
      intensities.push_back(fit->getIntensity());
      retention_times.push_back(fit->getRT());
      cf_labels.push_back(label_lfq);
    }
    consensus_feature_labels.push_back(cf_labels);
    consensus_feature_filenames.push_back(filenames);
    consensus_feature_intensites.push_back(intensities);
    consensus_feature_retention_times.push_back(retention_times);
  }

  // The output file of the MSstats converter (TODO Change to CSV file once store for CSV files has been implemented)
  TextFile csv_out;
  csv_out.addLine(String(rt_summarization_manual ? "RetentionTime,": "") + "ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType," + "Condition,BioReplicate,Run," + String(has_fraction ? "Fraction,": "") + "Intensity");
  
  // Regex definition for fragment ions
  boost::regex regex_msstats_FragmentIon("[abcxyz][0-9]+");

  // These are placeholder fragment annotations and peptide evidences in case the original ones are empty

  // Placeholder fragment annotation
  PeptideHit::PeakAnnotation new_peak_ann;
  new_peak_ann.annotation = na_string;
  new_peak_ann.charge = -1;
  std::vector< PeptideHit::PeakAnnotation > placeholder_fragment_annotations = {new_peak_ann};

  // Placeholder peptide evidence
  PeptideEvidence new_pep_ev;
  new_pep_ev.setProteinAccession(na_string);
  std::vector< PeptideEvidence > placeholder_peptide_evidences = {new_pep_ev};

  // From the MSstats user guide: endogenous peptides (use "L") or labeled reference peptides (use "H").
  String isotope_label_type = "L";
  if (is_isotope_label_type)
  {
    // use the channel_id information (?)
    isotope_label_type = "H";
  }
  const String delim(",");

  // Keeps track of unique peptides (Size of value set is 1)
  std::map< String, std::set<String > > peptideseq_to_accessions;

  // Stores all the lines that will be present in the final MSstats output,
  // We need to map peptide sequences to full features, because then we can ignore peptides
  // that are mapped to multiple proteins. We also need to map to the
  // intensities, such that we combine intensities over multiple retention times.
  map< String, map< MSstatsLine, set< pair<Intensity, Coordinate> > > > peptideseq_to_prefix_to_intensities;

  for (Size i = 0; i < features.size(); ++i)
  {
    const OpenMS::BaseFeature &base_feature = features[i];

    for (const OpenMS::PeptideIdentification &pep_id : base_feature.getPeptideIdentifications())
    {
      for (const OpenMS::PeptideHit & pep_hit : pep_id.getHits())
      {
        const std::vector< PeptideHit::PeakAnnotation > & original_fragment_annotations = pep_hit.getPeakAnnotations();
        const std::vector< PeptideEvidence > & original_peptide_evidences = pep_hit.getPeptideEvidences();

        // Decide whether to use original or placeholder iterator
        const std::vector< PeptideHit::PeakAnnotation > & fragment_annotations = (original_fragment_annotations.size() == 0) ? placeholder_fragment_annotations : original_fragment_annotations;
        const std::vector< PeptideEvidence> & peptide_evidences = (original_peptide_evidences.size() == 0) ? placeholder_peptide_evidences : original_peptide_evidences;

        // Variables of the peptide hit
        // MSstats User manual 3.7.3: Unknown precursor charge should be set to 0
        const Int precursor_charge = (std::max)(pep_hit.getCharge(), 0);
        const String & sequence = pep_hit.getSequence().toUnmodifiedString();

        // Have to combine all fragment annotations with all peptide evidences
        for (const OpenMS::PeptideHit::PeakAnnotation & frag_ann : fragment_annotations)
        {
          String fragment_ion = na_string;

          // Determine if the FragmentIon field can be assigned
          if (frag_ann.annotation != na_string)
          {
            std::set< std::string > frag_ions;
            boost::smatch sm;
            boost::regex_search(frag_ann.annotation, sm, regex_msstats_FragmentIon);
            frag_ions.insert(sm.begin(), sm.end());
            if (frag_ions.size() == 1)
            {
              for (const auto& frag_ions_elem : frag_ions)
              {
                fragment_ion = frag_ions_elem;
              }
            }
          }
          const Int frag_charge = (std::max)(frag_ann.charge, 0);

          for (const OpenMS::PeptideEvidence &pep_ev : peptide_evidences)
          {
            // Write new line for each run
            for (Size j = 0; j < consensus_feature_filenames[i].size(); j++)
            {
              const String &current_filename = consensus_feature_filenames[i][j];
              const Intensity intensity(consensus_feature_intensites[i][j]);
              const Coordinate retention_time(consensus_feature_retention_times[i][j]);
              const unsigned label(consensus_feature_labels[i][j]);

              const String & accession = pep_ev.getProteinAccession();
              peptideseq_to_accessions[sequence].insert(accession);

              const pair< String, unsigned> tpl1 = make_pair(current_filename, label);
              const unsigned sample = path_label_to_sample[tpl1];
              const unsigned fraction = path_label_to_fraction[tpl1];

              const pair< String, unsigned> tpl2 = make_pair(current_filename, fraction);

              // Resolve run
              const unsigned run = run_map[tpl2];  // MSstats run according to the file table
              const unsigned openms_fractiongroup = path_label_to_fractiongroup[tpl1];
              msstats_run_to_openms_fractiongroup[run] = openms_fractiongroup;

              // Assemble MSstats line
              MSstatsLine prefix(
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
                      //(has_fraction ? delim + String(fraction) : "") //not sure if this delim is correct
              );
              pair<Intensity, Coordinate> intensity_retention_time = make_pair(intensity, retention_time);
              peptideseq_to_prefix_to_intensities[sequence][prefix].insert(intensity_retention_time);
            }
          }
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

  // sanity check that the triples (peptide_sequence, precursor_charge, run) only appears once
  set< tuple<String, String, String> > peptideseq_precursor_charge_run;

  // test
  int count_similar = 0;

  for (const pair< String, set< String> > &peptideseq_accessions : peptideseq_to_accessions)
  {
    // Only write if unique peptide
    if (peptideseq_accessions.second.size() == 1)
    {
      for (const pair< MSstatsLine, set< pair< Intensity, Coordinate > > > &line :
              peptideseq_to_prefix_to_intensities[peptideseq_accessions.first])
      {
        // First, we collect all retention times and intensities
        set< Coordinate > retention_times;
        set< Intensity > intensities;
        for (const pair< Intensity, Coordinate >& p : line.second)
        {
          if (retention_times.find(p.second) != retention_times.end())
          {
            LOG_WARN <<  "Peptide ion appears multiple times at the same retention time. This is not expected" << endl;
          }
          else
          {
            retention_times.insert(p.second);
            intensities.insert(p.first);
          }
        }

        tuple<String, String, String > tpl = make_tuple(
                line.first.sequence(), line.first.precursor_charge(), line.first.run());
 
        if (peptideseq_precursor_charge_run.find(tpl) != peptideseq_precursor_charge_run.end())
        {
         count_similar += 1;
         // throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide ion appears multiple times for the same run!");
        }
        peptideseq_precursor_charge_run.insert(tpl);

        // If the rt summarization method is set to manual, we simply output all it,rt pairs
        if (rt_summarization_manual)
        {
          for (const pair< Intensity, Coordinate > &intensity : line.second)
          {
            csv_out.addLine(String(intensity.second) + ',' + line.first.toString() + ',' + String(intensity.first));
          }
        }
          // Otherwise, the intensities are resolved over the retention times
        else
        {
          Intensity intensity(0);
          if (retention_time_summarization_method == "max")
          {
            intensity = *(std::max_element(intensities.begin(), intensities.end()));
          }
          else if (retention_time_summarization_method == "min")
          {
            intensity = *(std::min_element(intensities.begin(), intensities.end()));
          }
          else if (retention_time_summarization_method == "mean")
          {
            intensity = meanIntensity(intensities);
          }
          else if (retention_time_summarization_method == "sum")
          {
            intensity = sumIntensity(intensities);
          }
          csv_out.addLine(line.first.toString() + delim + String(intensity));
        }

      }
    }
  }

  // Store the final assembled CSV file
  csv_out.store(filename);
}


void OpenMS::MSstatsFile::storeISO(const OpenMS::String &filename, const ConsensusMap &consensus_map,
                                   const OpenMS::ExperimentalDesign& design, const StringList& reannotate_filenames,
                                   const String& bioreplicate, const String& condition,
                                   const String& mixture, const String& retention_time_summarization_method)
{
  // Experimental Design file
  ExperimentalDesign::SampleSection sampleSection = design.getSampleSection();

  checkConditionISO_(sampleSection, bioreplicate, condition, mixture);

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
    cout << "WARNING: One feature might appear at multiple retention times in the output file. This is invalid input for MSstats. Combining of features over retention times is needed!" << endl;
  }

  typedef OpenMS::Peak2D::IntensityType Intensity;
  typedef OpenMS::Peak2D::CoordinateType Coordinate;

  ExperimentalDesign::MSFileSection msfile_section = design.getMSFileSection();

  // Extract the Spectra Filepath column from the design
  std::vector<String> design_filenames;
  for (ExperimentalDesign::MSFileSectionEntry const& f : msfile_section)
  {
    const String fn = File::basename(f.path);
    design_filenames.push_back(fn);
  }

  // Determine if the experiment has fractions
  const bool has_fraction = design.isFractionated();

  vector< OpenMS::BaseFeature> features;
  vector< String > spectra_paths;

  // For each ConsensusFeature, store several attributes
  vector< vector< String > > consensus_feature_filenames;           // Filenames of ConsensusFeature
  vector< vector< Intensity > > consensus_feature_intensites;       // Intensities of ConsensusFeature
  vector< vector< Coordinate > > consensus_feature_retention_times; // Retention times of ConsensusFeature
  vector< vector< unsigned > > consensus_feature_labels;            // Labels of ConsensusFeature

  features.reserve(consensus_map.size());

  if (reannotate_filenames.empty())
  {
    consensus_map.getPrimaryMSRunPath(spectra_paths);
  }
  else
  {
    spectra_paths = reannotate_filenames;
  }
  const auto& column_headers = consensus_map.getColumnHeaders(); // needed for label_id

  // Reduce spectra path to the basename of the files
  for (Size i = 0; i < spectra_paths.size(); ++i)
  {
    spectra_paths[i] = File::basename(spectra_paths[i]);
  }

  if (!checkUnorderedContent_(spectra_paths, design_filenames))
  {
    LOG_FATAL_ERROR << "The filenames (extension ignored) in the consensusXML file are not the same as in the experimental design" << endl;
    LOG_FATAL_ERROR << "Spectra files (consensus map): \n";
    for (auto const & s : spectra_paths)
    {
      LOG_FATAL_ERROR << s << endl;
    }
    LOG_FATAL_ERROR << "Spectra files (design): \n";
    for (auto const & s : design_filenames)
    {
      LOG_FATAL_ERROR << s << endl;
    }
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The filenames (extension ignored) in the consensusXML file are not the same as in the experimental design");
  }

  // Extract information from the consensus features.
  for (const ConsensusFeature &consensus_feature : consensus_map)
  {
    features.push_back(consensus_feature);

    vector< String > filenames;
    vector< Intensity > intensities;
    vector< Coordinate > retention_times;
    vector< unsigned > cf_labels;

    // Store the file names and the run intensities of this feature
    const ConsensusFeature::HandleSetType fs(consensus_feature.getFeatures());
    for (ConsensusFeature::HandleSetType::const_iterator fit = fs.begin(); fit != fs.end(); ++fit)
    {
      filenames.push_back(spectra_paths[fit->getMapIndex()]);
      intensities.push_back(fit->getIntensity());
      retention_times.push_back(fit->getRT());

      // Get the label_id form the file description MetaValue
      auto & column = column_headers.at(fit->getMapIndex());
      if (column.metaValueExists("channel_id"))
      {
        cf_labels.push_back(Int(column.getMetaValue("channel_id")));
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No channel_id metavalue assigned to mapindex in consensusMap.");
      }
    }
    consensus_feature_labels.push_back(cf_labels);
    consensus_feature_filenames.push_back(filenames);
    consensus_feature_intensites.push_back(intensities);
    consensus_feature_retention_times.push_back(retention_times);
  }

  // The output file of the MSstats converter (TODO Change to CSV file once store for CSV files has been implemented)
  TextFile csv_out;
  csv_out.addLine(String(rt_summarization_manual ? "RetentionTime,": "") + "ProteinName,PeptideSequence,Charge,Channel,Condition,BioReplicate,Run,Mixture,TechMixture," + String(has_fraction ? "Fraction,": "") + "Intensity");

  // These are placeholder peptide evidences in case the original ones are empty
  PeptideEvidence new_pep_ev;
  new_pep_ev.setProteinAccession(na_string);
  std::vector< PeptideEvidence > placeholder_peptide_evidences = {new_pep_ev};

  const String delim(",");

  // Keeps track of unique peptides (Size of value set is 1)
  std::map< String, std::set<String > > peptideseq_to_accessions;

  // Stores all the lines that will be present in the final MSstats output,
  // We need to map peptide sequences to full features, because then we can ignore peptides
  // that are mapped to multiple proteins. We also need to map to the
  // intensities, such that we combine intensities over multiple retention times.
  map< String, map< MSstatsTMTLine, set< pair<Intensity, Coordinate> > > > peptideseq_to_prefix_to_intensities;

  for (Size i = 0; i < features.size(); ++i)
  {
    const OpenMS::BaseFeature &base_feature = features[i];

    for (const OpenMS::PeptideIdentification &pep_id : base_feature.getPeptideIdentifications())
    {
      for (const OpenMS::PeptideHit & pep_hit : pep_id.getHits())
      {
        const std::vector< PeptideEvidence > & original_peptide_evidences = pep_hit.getPeptideEvidences();

        // Decide whether to use original or placeholder iterator
        const std::vector< PeptideEvidence> & peptide_evidences = (original_peptide_evidences.size() == 0) ? placeholder_peptide_evidences : original_peptide_evidences;

        // Variables of the peptide hit
        // MSstats User manual 3.7.3: Unknown precursor charge should be set to 0
        const Int precursor_charge = (std::max)(pep_hit.getCharge(), 0);
        const String & sequence = pep_hit.getSequence().toString();
          
          for (const OpenMS::PeptideEvidence &pep_ev : peptide_evidences)
          {
            // Write new line for each run
            for (Size j = 0; j < consensus_feature_filenames[i].size(); j++)
            {
              const String &current_filename = consensus_feature_filenames[i][j];

              const Intensity intensity(consensus_feature_intensites[i][j]);
              const Coordinate retention_time(consensus_feature_retention_times[i][j]);
              const unsigned channel(consensus_feature_labels[i][j]);

              const String & accession = pep_ev.getProteinAccession();
              peptideseq_to_accessions[sequence].insert(accession);
              
              const pair< String, unsigned> tpl1 = make_pair(current_filename, channel);
              const unsigned sample = path_label_to_sample[tpl1];
              const unsigned fraction = path_label_to_fraction[tpl1];

              // Resolve techmixture, run
              const unsigned openms_fractiongroup = path_label_to_fractiongroup[tpl1];
              String techmixture = String(sampleSection.getFactorValue(sample, mixture)) + "_" + String(openms_fractiongroup);
              String run = techmixture + (has_fraction ? String("_" + String(fraction)) : "");  
                          
              // Assemble MSstats line
              MSstatsTMTLine prefix(
                      has_fraction,
                      accession,
                      sequence,
                      precursor_charge,
                      String(channel),
                      sampleSection.getFactorValue(sample, condition),
                      sampleSection.getFactorValue(sample, bioreplicate),
                      String(run),
                      sampleSection.getFactorValue(sample, mixture),
                      String(techmixture),
                      (has_fraction ? String(fraction) : "")
              );
              pair<Intensity, Coordinate> intensity_retention_time = make_pair(intensity, retention_time);
              peptideseq_to_prefix_to_intensities[sequence][prefix].insert(intensity_retention_time);
            }
          }
        }
      }
    }

  // Print the run mapping between MSstats and OpenMS
  for (const auto& run_mapping : msstats_run_to_openms_fractiongroup)
  {
    cout << "MSstats run " << String(run_mapping.first)
         << " corresponds to OpenMS TechMixture " << String(run_mapping.second) << endl;
  }

  // sanity check that the triples (peptide_sequence, precursor_charge, run) only appears once
  set< tuple<String, String, String> > peptideseq_precursor_charge_run;

  // test
  int count_similar = 0;

  for (const pair< String, set< String> > &peptideseq_accessions : peptideseq_to_accessions)
  {
    // Only write if unique peptide
    if (peptideseq_accessions.second.size() == 1)
    {
      for (const pair< MSstatsTMTLine, set< pair< Intensity, Coordinate > > > &line :
              peptideseq_to_prefix_to_intensities[peptideseq_accessions.first])
      {
        // First, we collect all retention times and intensities
        set< Coordinate > retention_times;
        set< Intensity > intensities;
        for (const pair< Intensity, Coordinate >& p : line.second)
        {
          if (retention_times.find(p.second) != retention_times.end())
          {
            LOG_WARN <<  "Peptide ion appears multiple times at the same retention time. This is not expected" << endl;
          }
          else
          {
            retention_times.insert(p.second);
            intensities.insert(p.first);
          }
        }

        tuple<String, String, String > tpl = make_tuple(
                line.first.sequence(), line.first.precursor_charge(), line.first.run());
 
        if (peptideseq_precursor_charge_run.find(tpl) != peptideseq_precursor_charge_run.end())
        {
         count_similar += 1;
         // throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide ion appears multiple times for the same run!");
        }
        peptideseq_precursor_charge_run.insert(tpl);

        // If the rt summarization method is set to manual, we simply output all it,rt pairs
        if (rt_summarization_manual)
        {
          for (const pair< Intensity, Coordinate > &intensity : line.second)
          {
            csv_out.addLine(String(intensity.second) + ',' + line.first.toString() + ',' + String(intensity.first));
          }
        }
          // Otherwise, the intensities are resolved over the retention times
        else
        {
          Intensity intensity(0);
          if (retention_time_summarization_method == "max")
          {
            intensity = *(std::max_element(intensities.begin(), intensities.end()));
          }
          else if (retention_time_summarization_method == "min")
          {
            intensity = *(std::min_element(intensities.begin(), intensities.end()));
          }
          else if (retention_time_summarization_method == "mean")
          {
            intensity = meanIntensity(intensities);
          }
          else if (retention_time_summarization_method == "sum")
          {
            intensity = sumIntensity(intensities);
          }
          csv_out.addLine(line.first.toString() + delim + String(intensity));
        }

      }
    }
  }

  // Store the final assembled CSV file
  csv_out.store(filename);
}

