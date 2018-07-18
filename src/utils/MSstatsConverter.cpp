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
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <boost/regex.hpp>


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_MSstatsConverter

    @brief Converter to input for MSstats

    This util consumes an ID-mapped consensusXML file and OpenMS experimental design in TSV format to create a CSV file which can subsequently be used as input for the R package MSstats [1].

    [1] M. Choi et al. MSstats: an R package for statistical analysis for quantitative mass spectrometry-based proteomic experiments. Bioinformatics (2014), 30 (17): 2524-2526

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_MSstats.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_MSstats.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMSstatsConverter final :
public TOPPBase
{
public:

  TOPPMSstatsConverter() :
    TOPPBase("MSstatsConverter", "Converter to input for MSstats", false)
{
}

protected:

  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() final override
  {
    // Input consensusXML
    registerInputFile_(TOPPMSstatsConverter::param_in, "<in>", "", "Input consensusXML with peptide intensities", true, false);
    setValidFormats_(TOPPMSstatsConverter::param_in, ListUtils::create<String>("consensusXML"), true);

    registerInputFile_(TOPPMSstatsConverter::param_in_design, "<in_design>", "", "Experimental Design file", true, false);
    setValidFormats_(TOPPMSstatsConverter::param_in_design, ListUtils::create<String>("tsv"), true);

    registerStringOption_(TOPPMSstatsConverter::param_msstats_bioreplicate, "<msstats_bioreplicate>", "MSstats_BioReplicate", "Which column in the condition table should be used for MSstats 'BioReplicate'", false, false);
    registerStringOption_(TOPPMSstatsConverter::param_msstats_condition, "<msstats_condition>", "MSstats_Condition", "Which column in the condition table should be used for MSstats 'Condition'", false, false);

    // advanced option to overwrite MS file annotations in consensusXML
    registerInputFileList_("reannotate_filenames", "<file(s)>", StringList(), "Overwrite MS file names in consensusXML", false, true);

    // Isotope label type
    registerFlag_(TOPPMSstatsConverter::param_labeled_reference_peptides, "If set, IsotopeLabelType is 'H', else 'L'");

    // Specifies how peptide ions eluding at different retention times should be resolved
    registerStringOption_(TOPPMSstatsConverter::param_retention_time_summarization_method, "<retention_time_summarization_method>", "max", "How undistinguishable peptides at different retention times should be treated", false, true);
    setValidStrings_(TOPPMSstatsConverter::param_retention_time_summarization_method, ListUtils::create<String>("manual,max,min,mean,sum"));

    // Output CSV file
    registerOutputFile_(TOPPMSstatsConverter::param_out, "<out>", "", "Input CSV file for MSstats.", true, false);
    setValidFormats_(TOPPMSstatsConverter::param_out, ListUtils::create<String>("csv"));
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) final override
  {
    try
    {
      // Input file, must be consensusXML
      const String arg_in(getStringOption_(TOPPMSstatsConverter::param_in));
      const FileTypes::Type in_type(FileHandler::getType(arg_in));

      fatalErrorIf_(
              in_type != FileTypes::CONSENSUSXML,
              "Input type is not consensusXML!",
              ILLEGAL_PARAMETERS);

      // Tool arguments
      const String arg_out = getStringOption_(TOPPMSstatsConverter::param_out);
      const String arg_msstats_condition = getStringOption_(TOPPMSstatsConverter::param_msstats_condition);
      const String arg_msstats_bioreplicate = getStringOption_(TOPPMSstatsConverter::param_msstats_bioreplicate);
      const String arg_retention_time_summarization_method = getStringOption_(TOPPMSstatsConverter::param_retention_time_summarization_method);

      // Experimental Design file
      const String arg_in_design = getStringOption_(TOPPMSstatsConverter::param_in_design);
      const ExperimentalDesign design = ExperimentalDesignFile::load(arg_in_design, false);
      ExperimentalDesign::SampleSection sampleSection = design.getSampleSection();

      // Sample Section must contain the column that contains the condition used for MSstats
      fatalErrorIf_(
              sampleSection.hasFactor(arg_msstats_condition) == false,
              "Sample Section of experimental design does not contain MSstats_Condition",
              ILLEGAL_PARAMETERS
      );

      // Sample Section must contain column for the Bioreplicate
      fatalErrorIf_(
              sampleSection.hasFactor(arg_msstats_bioreplicate) == false,
              "Sample Section does not contain column for biological replicate",
              ILLEGAL_PARAMETERS
      );

      // assemble lookup table for run (each combination of pathname and fraction is a run)
      map< pair< String, unsigned>, unsigned > run_map;
      assembleRunMap(run_map, design);

      // Maps run in MSstats input to run for OpenMS
      map< unsigned, unsigned > msstats_run_to_openms_fractiongroup;

      // Mapping of filepath and label to sample and fraction
      map< pair< String, unsigned >, unsigned> path_label_to_sample = design.getPathLabelToSampleMapping(true);
      map< pair< String, unsigned >, unsigned> path_label_to_fraction = design.getPathLabelToFractionMapping(true);
      map< pair< String, unsigned >, unsigned> path_label_to_fractiongroup = design.getPathLabelToFractionGroupMapping(true);

      // The Retention Time is additionally written to the output as soon as the user wants to resolve multiple peptides manually
      const bool rt_summarization_manual(arg_retention_time_summarization_method == "manual");

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

      // Determine if the experimental design is LFQ (one label token in design file)
      const bool isLabelFree = (design.getNumberOfLabels() == 1);

      // Currently, we cannot support multiple labels for MSstats
      fatalErrorIf_(
        isLabelFree == false,
        "MSstatsConverter can only support label-free quantitation experiments",
        ILLEGAL_PARAMETERS
      );

      // label id 1 is used in case the experimental design specifies a LFQ experiment
      const unsigned label_lfq(1);

      vector< OpenMS::BaseFeature> features;
      vector< String > spectra_paths;

      // For each ConsensusFeature, store several attributes
      vector< vector< String > > consensus_feature_filenames;           // Filenames of ConsensusFeature
      vector< vector< Intensity > > consensus_feature_intensites;       // Intensites of ConsensusFeature
      vector< vector< Coordinate > > consensus_feature_retention_times; // Retention times of ConsensusFeature
      vector< vector< unsigned > > consensus_feature_labels;          // Labels of ConsensusFeature

      ConsensusMap consensus_map;
      features.reserve(consensus_map.size());
      ConsensusXMLFile().load(arg_in, consensus_map);

      StringList reannotate_filenames = getStringList_("reannotate_filenames");
      if (reannotate_filenames.empty())
      {
        consensus_map.getPrimaryMSRunPath(spectra_paths);
      }
      else
      {
        spectra_paths = reannotate_filenames;
      }
      ConsensusMap::ColumnHeaders& column_headers = consensus_map.getColumnHeaders(); // needed for label_id

      // Reduce spectra path to the basename of the files
      for (Size i = 0; i < spectra_paths.size(); ++i)
      {
        spectra_paths[i] = File::basename(spectra_paths[i]);
      }

      if (checkUnorderedContent_(spectra_paths, design_filenames) == false)
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
        return ILLEGAL_PARAMETERS;
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

          // If experiment is labelfree, us placeholder value for label (usually '1' in experimental design,
          // else get the label_id form the file description MetaValue
          cf_labels.push_back(
            isLabelFree ? label_lfq
                        : Int(column_headers[fit->getMapIndex()].getMetaValue("label_id"))
          );
        }
        consensus_feature_labels.push_back(cf_labels);
        consensus_feature_filenames.push_back(filenames);
        consensus_feature_intensites.push_back(intensities);
        consensus_feature_retention_times.push_back(retention_times);
      }

      // The output file of the MSstats converter (TODO Change to CSV file once store for CSV files has been implemented)
      TextFile csv_out;
      csv_out.addLine(String(rt_summarization_manual ? "RetentionTime,": "") + "ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,Condition,BioReplicate,Run," + String(has_fraction ? "Fraction,": "") + "Intensity");

      // Regex definition for fragment ions
      boost::regex regex_msstats_FragmentIon("[abcxyz][0-9]+");

      // These are placeholder fragment annotations and peptide evidences in case the original ones are empty

      // Placeholder fragment annotation
      PeptideHit::PeakAnnotation new_peak_ann;
      new_peak_ann.annotation = TOPPMSstatsConverter::na_string;
      new_peak_ann.charge = -1;
      std::vector< PeptideHit::PeakAnnotation > placeholder_fragment_annotations = {new_peak_ann};

      // Placeholder peptide evidence
      PeptideEvidence new_pep_ev;
      new_pep_ev.setProteinAccession(TOPPMSstatsConverter::na_string);
      std::vector< PeptideEvidence > placeholder_peptide_evidences = {new_pep_ev};

      // From the MSstats user guide: endogenous peptides (use "L") or labeled reference peptides (use "H").
      const String isotope_label_type = getFlag_(TOPPMSstatsConverter::param_labeled_reference_peptides) ? "H" : "L";
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
              String fragment_ion = TOPPMSstatsConverter::na_string;

              // Determine if the FragmentIon field can be assigned
              if (frag_ann.annotation != TOPPMSstatsConverter::na_string)
              {
                std::set< std::string > frag_ions;
                boost::smatch sm;
                boost::regex_search(frag_ann.annotation, sm, regex_msstats_FragmentIon);
                frag_ions.insert(sm.begin(), sm.end());
                if (frag_ions.size() == 1)
                {
                  for (auto frag_ions_elem : frag_ions)
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
                  const String &filename = consensus_feature_filenames[i][j];
                  const Intensity intensity(consensus_feature_intensites[i][j]);
                  const Coordinate retention_time(consensus_feature_retention_times[i][j]);
                  const unsigned label(consensus_feature_labels[i][j]);

                  const String & accession = pep_ev.getProteinAccession();
                  peptideseq_to_accessions[sequence].insert(accession);

                  const pair< String, unsigned> tpl1 = make_pair(filename, label);
                  const unsigned sample = path_label_to_sample[tpl1];
                  const unsigned fraction = path_label_to_fraction[tpl1];

                  const pair< String, unsigned> tpl2 = make_pair(filename, fraction);

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
                          sampleSection.getFactorValue(sample, arg_msstats_condition),
                          sampleSection.getFactorValue(sample, arg_msstats_bioreplicate),
                          String(run),
                          (has_fraction ? delim + String(fraction) : "")
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
      for (const auto run_mapping : msstats_run_to_openms_fractiongroup)
      {
        cout << "MSstats run " << String(run_mapping.first)
             << " corresponds to OpenMS fraction group " << String(run_mapping.second) << endl;
      }

      // sanity check that the triples (peptide_sequence, precursor_charge, run) only appears once
      set< tuple<String, String, String> > peptideseq_precursor_charge_run;

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
            for (const pair< Intensity, Coordinate > p : line.second)
            {
              if (retention_times.find(p.second) != retention_times.end())
              {
                LOG_WARN << "Peptide ion appears multiple times at the same retention time. This is not expected." << endl;
              }
              else
              {
                retention_times.insert(p.second);
                intensities.insert(p.first);
              }
            }

            tuple<String, String, String > tpl = make_tuple(
                    line.first.sequence(), line.first.precursor_charge(), line.first.run());
            fatalErrorIf_(
                    peptideseq_precursor_charge_run.find(tpl) != peptideseq_precursor_charge_run.end(),
                    "Peptide ion appears multiple times for the same run!",
                    ILLEGAL_PARAMETERS
            );
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
              if (arg_retention_time_summarization_method == "max")
              {
                intensity = *(std::max_element(intensities.begin(), intensities.end()));
              }
              else if (arg_retention_time_summarization_method == "min")
              {
                intensity = *(std::min_element(intensities.begin(), intensities.end()));
              }
              else if (arg_retention_time_summarization_method == "mean")
              {
                intensity = meanIntensity(intensities);
              }
              else if (arg_retention_time_summarization_method == "sum")
              {
                intensity = sumIntensity(intensities);
              }
              csv_out.addLine(line.first.toString() + delim + String(intensity));
            }

          }
        }
      }

      // Store the final assembled CSV file
      csv_out.store(arg_out);
      return EXECUTION_OK;
    }
    catch(const ExitCodes &exit_code)
    {
      return exit_code;
    }
      }

private:

  class MSstatsLine
  {
  public :
      MSstatsLine(
              bool _has_fraction,
              const String& _accession,
              const String& _sequence,
              const String& _precursor_charge,
              const String& _fragment_ion,
              const String& _frag_charge,
              const String& _isotope_label_type,
              const String& _condition,
              const String& _bioreplicate,
              const String& _run,
              const String& _fraction
      ): has_fraction_(_has_fraction),
         accession_(_accession),
         sequence_(_sequence),
         precursor_charge_(_precursor_charge),
         fragment_ion_(_fragment_ion),
         frag_charge_(_frag_charge),
         isotope_label_type_(_isotope_label_type),
         condition_(_condition),
         bioreplicate_(_bioreplicate),
         run_(_run),
         fraction_(_fraction) {}

      const String& accession() const {return this->accession_;}
      const String& sequence() const {return this->sequence_;}
      const String& precursor_charge() const {return this->precursor_charge_;}
      const String& run() const {return this->run_;}

      String toString() const
      {
        const String delim(",");
        return  accession_
                + delim + sequence_
                + delim + precursor_charge_
                + delim + fragment_ion_
                + delim + frag_charge_
                + delim + isotope_label_type_
                + delim + condition_
                + delim + bioreplicate_
                + delim + run_
                + (this->has_fraction_ ? delim + String(fraction_) : "");
      }

      friend bool operator<(const MSstatsLine &l,
                            const MSstatsLine &r) {

        return std::tie(l.accession_, l.run_, l.condition_, l.bioreplicate_, l.precursor_charge_, l.sequence_) <
               std::tie(r.accession_, r.run_, r.condition_, r.bioreplicate_, r.precursor_charge_, r.sequence_);
      }


  private:
      bool has_fraction_;
      String accession_;
      String sequence_;
      String precursor_charge_;
      String fragment_ion_;
      String frag_charge_;
      String isotope_label_type_;
      String condition_;
      String bioreplicate_;
      String run_;
      String fraction_;
  };

  static const String param_in;
  static const String param_in_design;
  static const String param_msstats_bioreplicate;
  static const String param_msstats_condition;
  static const String param_out;
  static const String param_labeled_reference_peptides;
  static const String param_retention_time_summarization_method;

  static const String na_string;

  // The meta value of the peptide identification which is going to be used for the experimental design link
  static const String meta_value_exp_design_key;


  static void fatalErrorIf_(const bool error_condition, const String &message, const int exit_code)
  {
    if (error_condition)
    {
      LOG_FATAL_ERROR << "FATAL: " << message << std::endl;
      throw exit_code;
    }
  }

  /*
   *  MSstats treats runs differently than OpenMS. In MSstats, runs are an enumeration of (SpectraFilePath, Fraction)
   *  In OpenMS, a run is split into multiple fractions.
   *
   */
  static void assembleRunMap(
          std::map< std::pair< String, unsigned>, unsigned> &run_map,
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

  bool checkUnorderedContent_(const std::vector< String> &first, const std::vector< String > &second)
  {
    const std::set< String > lhs(first.begin(), first.end());
    const std::set< String > rhs(second.begin(), second.end());
    return lhs.size() == rhs.size()
        && std::equal(lhs.begin(), lhs.end(), rhs.begin());
  }

  OpenMS::Peak2D::IntensityType sumIntensity(const set< OpenMS::Peak2D::IntensityType > &intensities)
  {
    OpenMS::Peak2D::IntensityType result = 0;
    for (const OpenMS::Peak2D::IntensityType &intensity : intensities)
    {
      result += intensity;
    }
    return result;
  }

  OpenMS::Peak2D::IntensityType meanIntensity(const set< OpenMS::Peak2D::IntensityType > &intensities)
  {
    return sumIntensity(intensities) / intensities.size();
  }
};

const String TOPPMSstatsConverter::param_in = "in";
const String TOPPMSstatsConverter::param_in_design = "in_design";
const String TOPPMSstatsConverter::param_msstats_bioreplicate = "msstats_bioreplicate";
const String TOPPMSstatsConverter::param_msstats_condition = "msstats_condition";
const String TOPPMSstatsConverter::param_out = "out";
const String TOPPMSstatsConverter::na_string = "NA";
const String TOPPMSstatsConverter::param_labeled_reference_peptides = "labeled_reference_peptides";
const String TOPPMSstatsConverter::meta_value_exp_design_key = "spectra_data";
const String TOPPMSstatsConverter::param_retention_time_summarization_method = "retention_time_summarization_method";


// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPMSstatsConverter tool;
  return tool.main(argc, argv);
}
/// @endcond
