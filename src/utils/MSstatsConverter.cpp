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
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
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

    [1] M. Choi et al. “MSstats: an R package for statistical analysis for quantitative mass spectrometry-based proteomic experiments.” Bioinformatics (2014), 30 (17): 2524-2526

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

  static const String param_in_consensusxml;
  static const String param_in_design_run;
  static const String param_in_design_condition;

  static const String param_msstats_bioreplicate;
  static const String param_msstats_condition;

  static const String param_out;
  static const String param_labeled_reference_peptides;

  static const String na_string;

  // The meta value of the peptide identification which is going to be used for the experimental design link
  static const String meta_value_exp_design_key;

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
    this->registerInputFile_(TOPPMSstatsConverter::param_in_consensusxml, "<in_consensusxml>", "", "Input consensusXML with peptide intensities", true, false);
    this->setValidFormats_(TOPPMSstatsConverter::param_in_consensusxml, ListUtils::create<String>("consensusXML"));

    this->_registerExperimentalDesignInputFile(TOPPMSstatsConverter::param_in_design_run, "<in_design_run>", "TSV file containing the run description (run table)");
    this->_registerExperimentalDesignInputFile(TOPPMSstatsConverter::param_in_design_condition, "<in_design_condition>", "TSV file containing the condition description (condition table)");

    this->registerStringOption_(TOPPMSstatsConverter::param_msstats_bioreplicate, "<msstats_bioreplicate>", "Biological Replicate", "Which column in the condition table should be used for MSstats 'BioReplicate'", false, false);
    this->registerStringOption_(TOPPMSstatsConverter::param_msstats_condition, "<msstats_condition>", "", "Which column in the condition table should be used for MSstats 'Condition'", true, false);

    // Isotope label type
    this->registerFlag_(TOPPMSstatsConverter::param_labeled_reference_peptides, "If set, IsotopeLabelType is 'H', else 'L'");

    // Output CSV file
    this->registerOutputFile_(TOPPMSstatsConverter::param_out, "<out>", "", "Input CSV file for MSstats.", true, false);
    this->setValidFormats_(TOPPMSstatsConverter::param_out, ListUtils::create<String>("csv"));
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) final override
  {
    try
    {
      DesignFile file_run(this->getStringOption_(TOPPMSstatsConverter::param_in_design_run), ListUtils::create<String>("Run,Condition"), "Spectra File");
      DesignFile file_condition(this->getStringOption_(TOPPMSstatsConverter::param_in_design_condition), ListUtils::create<String>("Biological Replicate"), "Condition");

      const String & arg_msstats_condition = this->getStringOption_(TOPPMSstatsConverter::param_msstats_condition);
      const String & arg_msstats_bioreplicate = this->getStringOption_(TOPPMSstatsConverter::param_msstats_bioreplicate);

      if (file_condition.isColumnName(arg_msstats_condition) == false || file_condition.isColumnName(arg_msstats_bioreplicate) ==  false)
      {
        LOG_FATAL_ERROR << "FATAL: At least one of the specified column names is not part of condition table!" << std::endl;
        return ILLEGAL_PARAMETERS;
      }


      // Read the input files
      ConsensusMap consensus_map;
      ConsensusXMLFile().load(this->getStringOption_(TOPPMSstatsConverter::param_in_consensusxml), consensus_map);

      // The output file of the MSstats converter (TODO Change to CSV file once store for CSV files has been implemented)
      TextFile csv_out;

      // Add the header line (With fraction if that has been specified in the experimental design)
      const bool has_fraction = file_run.isColumnName("Fraction");

      csv_out.addLine("ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,Condition,BioReplicate,Run,Intensity" + String(has_fraction ? ",Fraction" : ""));

      // Regex definition for fragment ions
      std::regex regex_msstats_FragmentIon("[abcxyz][0-9]+");

      // Iterate protein identifications and collect spectra_data metavalue and ID
      std::map< String, std::set< String> > protid_to_filenames;
      for (const OpenMS::ProteinIdentification & protein_identification : consensus_map.getProteinIdentifications())
      {
        const String & identifier = protein_identification.getIdentifier();

        if (protein_identification.metaValueExists(TOPPMSstatsConverter::meta_value_exp_design_key) == false)
        {
          LOG_FATAL_ERROR << "FATAL: ProteinIdentification does not have meta value for original file. Cannot continue!" << std::endl;
          return ILLEGAL_PARAMETERS;
        }
        const StringList & exp_design_key = protein_identification.getMetaValue(TOPPMSstatsConverter::meta_value_exp_design_key).toStringList();
        for (const String & meta_value : exp_design_key)
        {
          protid_to_filenames[identifier].insert(File::basename(meta_value));
        }
      }

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

      // Keeps track of precursor charges and accessions of a peptide sequence to avoid duplicate lines in output (logic currently not used)
      //String previous_sequence = "";
      //std::set< Int > previous_precursor_charges;
      //std::set< String > previous_prot_accs;

      for (const OpenMS::ConsensusFeature & consensus_feature : consensus_map)
      {
        Peak2D::IntensityType intensity = consensus_feature.getIntensity();
        assert(intensity > 0);

        for (const OpenMS::PeptideIdentification & pep_id : consensus_feature.getPeptideIdentifications())
        {
          // Get runs that belong to this identifier
          std::set< String > filenames;
          assert(protid_to_filenames.find(pep_id.getIdentifier()) != protid_to_filenames.end());

          for (const String & filename : protid_to_filenames[pep_id.getIdentifier()])
          {
            // Test whether the experimental design specifies the encountered filename
            if (file_run.isRowName(filename) == false)
            {
              LOG_FATAL_ERROR << "FATAL: Experimental design does not contain information on file " << filename << ". Cannot continue!" << std::endl;
              return ILLEGAL_PARAMETERS;
            }
            filenames.insert(filename);
          }

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

            // Have to combine all fragment annotations with all peptide evidences
            for (const String & frag_ann : fragment_annotations)
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

              for (const OpenMS::PeptideEvidence & pep_ev : peptide_evidences)
              {
                // Write new line for each protein accession and for each run
                for (const String & filename : filenames)
                {
                  const String & accession = pep_ev.getProteinAccession();
                  const String & sequence = pep_hit.getSequence().toUnmodifiedString();

                  // Find the Biological replicate by filename
                  const String & condition = file_run.get(filename, "Condition");


                  // The peptide sequence has changed, the charges and accessions of the previous peptide can be removed
                  //if (sequence != previous_sequence)
                  //{
                  //  previous_precursor_charges.clear();
                  //  previous_prot_accs.clear();
                  //  previous_sequence = sequence;
                  //}

                  // Determine whether we need to write the current protein hit (precursor charge or accession changes)
                  //if (   previous_prot_accs.find(accession) == previous_prot_accs.end()
                  //    || previous_precursor_charges.find(precursor_charge) == previous_precursor_charges.end())
                  //{
                  csv_out.addLine(  accession
                                    + ',' + sequence
                                    + ',' + precursor_charge
                                    + ',' + fragment_ion
                                    + ',' + frag_charge
                                    + ',' + isotope_label_type
                                    + ',' + file_condition.get(condition, arg_msstats_condition)
                                    + ',' + file_condition.get(condition, arg_msstats_bioreplicate)
                                    + ',' + file_run.get(filename, "Run")
                                    + ',' + intensity
                                    + String(has_fraction ? ("," + file_run.get(filename, "Fraction")) : ""));
                  //previous_prot_accs.insert(accession);
                  //previous_precursor_charges.insert(precursor_charge);
                  //}
                }
              }
            }
          }
        }
      }
      // Store the final assembled CSV file
      csv_out.store(this->getStringOption_(TOPPMSstatsConverter::param_out));
      return EXECUTION_OK;
    }
    catch(int)
    {
      return INPUT_FILE_CORRUPT;
    }
      }

private:

  static void _conditionalFatalError(const String & message, bool error_condition)
  {
    if (error_condition)
    {
      LOG_FATAL_ERROR << "FATAL: " << message << std::endl;
      throw 1;
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
        it->split("\t", line);

        Size n_entries(line.size());
        // Trim all entries
        for (Size i = 0; i < n_entries; ++i)
        {
          line[i] = line[i].trim();
          _conditionalFatalError("Entry in design table is not allowed to be empty!", line[i].empty());
        }

        // Line is record
        if (this->_n_columns > 0)
        {
          _conditionalFatalError("Conflicting number of entries in design table: " + String(n_entries) + " vs. " + String(this->_n_columns), n_entries != this->_n_columns);
          const String & row_name = line[this->_columnname_to_columnindex[index_column]];
          _conditionalFatalError("Row name " + row_name + " appears multiple times!", row_names.find(row_name) != row_names.end());
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
              _conditionalFatalError("Header names in design table must be unique, but " + line[i] + " appears several times!", line[i] == line[j]);
            }
            // Remember the index at which this header appears
            this->_columnname_to_columnindex[line[i]] = i;
          }

          // Make sure that all required headers exist
          for (const String & header : required_headers)
          {
            _conditionalFatalError("Header '" + header + "' does not exist in input design file!", std::find(line.begin(), line.end(), header) == line.end());
          }
          // Make sure that the index column appears in the header
          _conditionalFatalError("Index column is not a header!", std::find(line.begin(), line.end(), index_column) == line.end());
          this->_n_columns = n_entries;
        }
      }
  }

    inline String get(const String & row_name, const String & column_name)
    {
      _conditionalFatalError("Tried to access invalid row or column name!", this->isRowName(row_name) == false || this->isColumnName(column_name) == false);
      return (this->_entries[this->_rowname_to_rowindex[row_name]])[this->_columnname_to_columnindex[column_name]];
    }

    inline bool isRowName(const String & row_name)
    {
      return this->_rowname_to_rowindex.find(row_name) != this->_rowname_to_rowindex.end();
    }

    inline bool isColumnName(const String & column_name)
    {
      return this->_columnname_to_columnindex.find(column_name) != this->_columnname_to_columnindex.end();
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


  // Advances
  //this->registerInputFile_(TOPPMSstatsConverter::param_in_experimental_design, "<in_experimental_design>", "",
  //                           "Experimental design as CSV file. The required columns are FileName,Condition,BioReplicate,Run", true, false

  void _registerExperimentalDesignInputFile(const String & param_name, const String & argument, const String & description)
  {
    static const StringList valid_formats = ListUtils::create<String>("tsv");
    this->registerInputFile_(param_name, argument, "", description, true, false);
    this->setValidFormats_(param_name, valid_formats, true);
  }
};

const String TOPPMSstatsConverter::param_in_consensusxml = "in";
const String TOPPMSstatsConverter::param_in_design_run = "in_design_run";
const String TOPPMSstatsConverter::param_in_design_condition = "in_design_condition";

const String TOPPMSstatsConverter::param_msstats_bioreplicate = "msstats_bioreplicate";
const String TOPPMSstatsConverter::param_msstats_condition = "msstats_condition";

const String TOPPMSstatsConverter::param_out = "out";
const String TOPPMSstatsConverter::na_string = "NA";
const String TOPPMSstatsConverter::param_labeled_reference_peptides = "labeled_reference_peptides";
const String TOPPMSstatsConverter::meta_value_exp_design_key = "spectra_data";




// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPMSstatsConverter tool;
  return tool.main(argc, argv);
}
/// @endcond
