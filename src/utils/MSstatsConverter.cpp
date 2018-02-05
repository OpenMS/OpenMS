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

#include <regex>


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_MSstatsConverter

    @brief Converter to input for MSstats

    This util consumes an ID-mapped consensusXML file and OpenMS experimental design in TSV format to create a file which can subsequently be used as input for the R package MSstats [1].

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
    this->registerInputFile_(TOPPMSstatsConverter::param_in, "<in>", "", "Input consensusXML or featureXML with peptide intensities", true, false);
    this->setValidFormats_(TOPPMSstatsConverter::param_in, ListUtils::create<String>("consensusXML"), true);

    this->registerExperimentalDesignInputFile_(TOPPMSstatsConverter::param_in_design_run, "<in_design_run>", "TSV file containing the run description (run table)");
    this->registerExperimentalDesignInputFile_(TOPPMSstatsConverter::param_in_design_condition, "<in_design_condition>", "TSV file containing the condition description (condition table)");

    this->registerStringOption_(TOPPMSstatsConverter::param_msstats_bioreplicate, "<msstats_bioreplicate>", "Biological Replicate", "Which column in the condition table should be used for MSstats 'BioReplicate'", false, false);
    this->registerStringOption_(TOPPMSstatsConverter::param_msstats_condition, "<msstats_condition>", "", "Which column in the condition table should be used for MSstats 'Condition'", true, false);

    // Isotope label type
    this->registerFlag_(TOPPMSstatsConverter::param_labeled_reference_peptides, "If set, IsotopeLabelType is 'H', else 'L'");

    // Non-unique Peptides
    this->registerFlag_(TOPPMSstatsConverter::param_ambiguous_peptides, "If set, the output CSV file can contain peptides that have been assigned to multiple protein ids. Attention: you normally do not want to do this for MSstats", true);

    // Specifies how peptide ions eluding at different retention times should be resolved
    this->registerStringOption_(TOPPMSstatsConverter::param_retention_time_resolution_method, "<retention_time_resolution_method>", "", "How undistinguishable peptides at different retention times should be treated", true, false);
    this->setValidStrings_(TOPPMSstatsConverter::param_retention_time_resolution_method, ListUtils::create<String>("manual"));

    // Output CSV file
    this->registerOutputFile_(TOPPMSstatsConverter::param_out, "<out>", "", "Input CSV file for MSstats.", true, false);
    this->setValidFormats_(TOPPMSstatsConverter::param_out, ListUtils::create<String>("csv"));
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) final override
      {
    try
    {
      // Tool arguments
      const String arg_out(this->getStringOption_(TOPPMSstatsConverter::param_out));
      const String &arg_msstats_condition = this->getStringOption_(TOPPMSstatsConverter::param_msstats_condition);
      const String &arg_msstats_bioreplicate = this->getStringOption_(TOPPMSstatsConverter::param_msstats_bioreplicate);
      const String &arg_retention_time_resolution_method = this->getStringOption_(TOPPMSstatsConverter::param_retention_time_resolution_method);

      // The Retention Time is additionally written to the output as soon as the user wants to resolve multiple peptides manually
      const bool write_retention_time(arg_retention_time_resolution_method == "manual");

      // Load the experimental design
      DesignFile file_run(this->getStringOption_(TOPPMSstatsConverter::param_in_design_run), ListUtils::create<String>("Run,Condition"), "Spectra File");
      DesignFile file_condition(this->getStringOption_(TOPPMSstatsConverter::param_in_design_condition), ListUtils::create<String>("Biological Replicate"), "Condition");
      conditionalFatalError_(
          "At least one of the specified column names is not part of condition table!",
          file_condition.isColumnName(arg_msstats_condition) == false || file_condition.isColumnName(arg_msstats_bioreplicate) ==  false,
          ILLEGAL_PARAMETERS);

      typedef OpenMS::Peak2D::IntensityType Intensity;
      typedef OpenMS::Peak2D::CoordinateType Coordinate;

      // Load the filenames of the design run file
      vector< String > design_run_filenames;
      file_run.getRowNames(design_run_filenames);

      conditionalFatalError_(
          "File Names in Design Run File are not unique. Cannot continue.",
          this->checkUnique_(design_run_filenames) == false,
          ILLEGAL_PARAMETERS);

      // Add the header line (With fraction if that has been specified in the experimental design)
      const bool has_fraction(file_run.isColumnName("Fraction"));

      // Input file, must be featureXML or consensusXML
      const String arg_in(this->getStringOption_(TOPPMSstatsConverter::param_in));
      const FileTypes::Type in_type(FileHandler::getType(arg_in));

      conditionalFatalError_(
          "Input type is not consensusXML!",
          (in_type != FileTypes::CONSENSUSXML),
          ILLEGAL_PARAMETERS);

      std::vector< OpenMS::BaseFeature> features;
      std::vector< OpenMS::String > spectra_paths;

      // For each ConsensusFeature, store the filenames
      std::vector< std::vector< String > > consensus_feature_filenames;
      // For each ConsensusFeature, store the intensities
      std::vector< std::vector< OpenMS::Peak2D::IntensityType > > consensus_feature_intensites;
      // For each ConsensusFeature, store the retention_times
      std::vector< std::vector< OpenMS::Peak2D::CoordinateType > > consensus_feature_retention_times;

      ConsensusMap consensus_map;
      features.reserve(consensus_map.size());
      ConsensusXMLFile().load(arg_in, consensus_map);
      consensus_map.getPrimaryMSRunPath(spectra_paths);

      // Reduce spectra path to the basename of the files
      for (Size i = 0; i < spectra_paths.size(); ++i)
      {
        spectra_paths[i] = File::basename(spectra_paths[i]);
      }
      conditionalFatalError_(
          "The number of files in the consensusXML does not coincide with the number of files in the experimental design!",
          spectra_paths.size() != design_run_filenames.size(),
          ILLEGAL_PARAMETERS);

      conditionalFatalError_(
          "The spectra file names in the consensusXML are not unique!",
          this->checkUnique_(spectra_paths) == false,
          ILLEGAL_PARAMETERS);

      conditionalFatalError_(
          "The filenames in the consensusXML file are not the same as in the experimental design",
          this->checkUnorderedContent_(spectra_paths, design_run_filenames) == false,
          ILLEGAL_PARAMETERS);

      for (const ConsensusFeature &consensus_feature : consensus_map)
      {
        features.push_back(consensus_feature);

        vector< String > filenames;
        vector< Intensity > intensities;
        vector< Coordinate > retention_times;

        // Store the file names and the run intensities of this feature
        const ConsensusFeature::HandleSetType fs(consensus_feature.getFeatures());
        for (ConsensusFeature::HandleSetType::const_iterator fit = fs.begin(); fit != fs.end(); ++fit)
        {
          filenames.push_back(spectra_paths[fit->getMapIndex()]);
          intensities.push_back(fit->getIntensity());
          retention_times.push_back(fit->getRT());
        }
        consensus_feature_filenames.push_back(filenames);
        consensus_feature_intensites.push_back(intensities);
        consensus_feature_retention_times.push_back(retention_times);
      }

      // The output file of the MSstats converter (TODO Change to CSV file once store for CSV files has been implemented)
      TextFile csv_out;
      csv_out.addLine(String(write_retention_time ? "RetentionTime,": "") + "ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,Condition,BioReplicate,Run," + String(has_fraction ? "Fraction,": "") + "Intensity");

      // Regex definition for fragment ions
      std::regex regex_msstats_FragmentIon("[abcxyz][0-9]+");

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
      const String isotope_label_type = this->getFlag_(TOPPMSstatsConverter::param_labeled_reference_peptides) ? "H" : "L";
      const char delim(',');

      // Keeps track of unique peptides (Size of value set is 1)
      std::map< String, std::set<String > > peptideseq_to_accessions;

      // Stores all the lines that will be present in the final MSstats output,
      map< String, map< String, set< pair<Intensity, Coordinate> > > > peptideseq_to_prefix_to_intensities;

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
                std::smatch sm;
                std::regex_search(frag_ann.annotation, sm, regex_msstats_FragmentIon);
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
                  const OpenMS::Peak2D::IntensityType intensity(consensus_feature_intensites[i][j]);
                  const OpenMS::Peak2D::CoordinateType retention_time(consensus_feature_retention_times[i][j]);
                  const String & accession = pep_ev.getProteinAccession();
                  peptideseq_to_accessions[sequence].insert(accession);

                  const String & condition = file_run.get(filename, "Condition");
                  const String & bioreplicate = file_condition.get(condition, arg_msstats_bioreplicate);
                  const String & fraction = String(has_fraction ? (delim + file_run.get(filename, "Fraction") + delim) : "");

                  const String prefix(
                      accession
                      + delim + sequence
                      + delim + precursor_charge
                      + delim + fragment_ion
                      + delim + frag_charge
                      + delim + isotope_label_type
                      + delim + file_condition.get(condition, arg_msstats_condition)
                      + delim + bioreplicate
                      + delim + file_run.get(filename, "Run")
                      + fraction);
                  pair<Intensity, Coordinate> intensity_retention_time;
                  intensity_retention_time.first = intensity;
                  intensity_retention_time.second = retention_time;

                  peptideseq_to_prefix_to_intensities[sequence][prefix].insert(intensity_retention_time);
                }
              }
            }
          }
        }
      }

      const bool write_ambigous_peptides(this->getFlag_(TOPPMSstatsConverter::param_ambiguous_peptides));

      for (const pair< String, set< String> > &peptideseq_accessions : peptideseq_to_accessions)
      {
        // Only write if unique peptide
        if (write_ambigous_peptides || peptideseq_accessions.second.size() == 1)
        {
          for (const pair< String, set< pair< Intensity, Coordinate > > > &line : peptideseq_to_prefix_to_intensities[peptideseq_accessions.first])
          {

        	set< Coordinate > retention_times;

            for (const pair< Intensity, Coordinate > &intensity : line.second)
            {
              // We do not expect that we encounter the same retention time twice
              conditionalFatalError_(
                "Peptide ion appears multiple times at the same retention time. This is not expected",
                retention_times.find(intensity.second) != retention_times.end(),
				ILLEGAL_PARAMETERS);

              retention_times.insert(intensity.second);
              csv_out.addLine( (write_retention_time ? String(intensity.second) + ',' : "") + line.first + ',' + String(intensity.first));
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

  static const String param_in;
  static const String param_in_design_run;
  static const String param_in_design_condition;
  static const String param_msstats_bioreplicate;
  static const String param_msstats_condition;
  static const String param_out;
  static const String param_labeled_reference_peptides;
  static const String param_ambiguous_peptides;
  static const String param_retention_time_resolution_method;

  static const String na_string;

  // The meta value of the peptide identification which is going to be used for the experimental design link
  static const String meta_value_exp_design_key;


  static void conditionalFatalError_(const String & message, bool error_condition, int exit_code)
  {
    if (error_condition)
    {
      LOG_FATAL_ERROR << "FATAL: " << message << std::endl;
      throw exit_code;
    }
  }

  class DesignFile final
  {
  public:
    DesignFile(const String & filename, const std::vector< String > & required_headers, const String & index_column)
  : _n_columns(0), _entries(), _rowname_to_rowindex(), _columnname_to_columnindex()
  {
      TextFile input_file;
      input_file.load(filename, true, -1, true);

      std::set< String > row_names;
      Size index_row = -1;

      for (TextFile::ConstIterator it = input_file.begin(); it != input_file.end(); ++it)
      {
        ++index_row;
        std::vector< String > line;

        // Because we assume that the experimental design file is tab separated
        it->split("\t", line);

        Size n_entries(line.size());

        // Trim all entries
        for (Size i = 0; i < n_entries; ++i)
        {
          line[i] = line[i].trim();
          conditionalFatalError_("Entry in design table is not allowed to be empty!", line[i].empty(), INPUT_FILE_CORRUPT);
        }

        // Line is record
        if (this->_n_columns > 0)
        {
          conditionalFatalError_("Conflicting number of entries in design table: " + String(n_entries) + " vs. " + String(this->_n_columns), n_entries != this->_n_columns, INPUT_FILE_CORRUPT);
          const String & row_name = line[this->_columnname_to_columnindex[index_column]];
          conditionalFatalError_("Row name " + row_name + " appears multiple times!", row_names.find(row_name) != row_names.end(), INPUT_FILE_CORRUPT);
          this->_entries.push_back(line);
          this->_rowname_to_rowindex[row_name] = row_names.size();
          row_names.insert(row_name);
        }
        // Line is header
        else
        {
          for (Size i = 0; i < n_entries; ++i)
          {
            // Ensure that the header lines are unique
            for (Size j = 0; j < i; ++j)
            {
              conditionalFatalError_("Header names in design table must be unique, but " + line[i] + " appears several times!", line[i] == line[j], INPUT_FILE_CORRUPT);
            }
            // Remember the index at which this header appears
            this->_columnname_to_columnindex[line[i]] = i;
          }

          // Make sure that all required headers exist
          for (const String & header : required_headers)
          {
            conditionalFatalError_("Header '" + header + "' does not exist in input design file!", std::find(line.begin(), line.end(), header) == line.end(), INPUT_FILE_CORRUPT);
          }
          // Make sure that the index column appears in the header
          conditionalFatalError_("Index column is not a header!", std::find(line.begin(), line.end(), index_column) == line.end(), INPUT_FILE_CORRUPT);
          this->_n_columns = n_entries;
        }
      }
  }

    inline String get(const String & row_name, const String & column_name)
    {
      conditionalFatalError_("Tried to access invalid row or column name!", this->isRowName(row_name) == false || this->isColumnName(column_name) == false, ExitCodes::INCOMPATIBLE_INPUT_DATA);
      return (this->_entries[this->_rowname_to_rowindex[row_name]])[this->_columnname_to_columnindex[column_name]];
    }

    inline bool isRowName(const String & row_name) const
    {
      return this->_rowname_to_rowindex.find(row_name) != this->_rowname_to_rowindex.end();
    }

    inline bool isColumnName(const String & column_name) const
    {
      return this->_columnname_to_columnindex.find(column_name) != this->_columnname_to_columnindex.end();
    }

    void getRowNames(std::vector< String > &rownames) const
    {
      rownames.clear();
      rownames.resize(this->_rowname_to_rowindex.size());
      for (const std::pair< String, Size >  &pair : this->_rowname_to_rowindex)
      {
        rownames[pair.second] = pair.first;
      }
    }


  private:

    // Number of columns in the file
    Size _n_columns;

    // The entries of the file
    std::vector< std::vector < String > > _entries;

    // Maps the row_name (column can be specified in constructor) to row index
    std::map< String, Size > _rowname_to_rowindex;

    // Maps the column name to the index
    std::map< String, Size > _columnname_to_columnindex;
  };

  void registerExperimentalDesignInputFile_(const String & param_name, const String & argument, const String & description)
  {
    static const StringList valid_formats = ListUtils::create<String>("tsv");
    this->registerInputFile_(param_name, argument, "", description, true, false);
    this->setValidFormats_(param_name, valid_formats, true);
  }

  bool checkUnorderedContent_(const std::vector< String> &first, const std::vector< String > &second)
  {
    const std::set< String > lhs(first.begin(), first.end());
    const std::set< String > rhs(second.begin(), second.end());
    return lhs.size() == rhs.size()
        && std::equal(lhs.begin(), lhs.end(), rhs.begin());
  }

  /**
   * Checks whether the elements in a vector of OpenMS::String are all unique
   */
  bool checkUnique_(const OpenMS::StringList &vec) const
  {
    const std::set< String > vec_set(vec.begin(), vec.end());
    return vec_set.size() == vec.size();
  }
};

const String TOPPMSstatsConverter::param_in = "in";
const String TOPPMSstatsConverter::param_in_design_run = "in_design_run";
const String TOPPMSstatsConverter::param_in_design_condition = "in_design_condition";
const String TOPPMSstatsConverter::param_msstats_bioreplicate = "msstats_bioreplicate";
const String TOPPMSstatsConverter::param_msstats_condition = "msstats_condition";
const String TOPPMSstatsConverter::param_out = "out";
const String TOPPMSstatsConverter::na_string = "NA";
const String TOPPMSstatsConverter::param_labeled_reference_peptides = "labeled_reference_peptides";
const String TOPPMSstatsConverter::meta_value_exp_design_key = "spectra_data";
const String TOPPMSstatsConverter::param_ambiguous_peptides = "ambiguous_peptides";
const String TOPPMSstatsConverter::param_retention_time_resolution_method = "retention_time_resolution_method";


// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPMSstatsConverter tool;
  return tool.main(argc, argv);
}
/// @endcond
