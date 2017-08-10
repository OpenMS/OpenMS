// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
#include <OpenMS/FORMAT/CsvFile.h>
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

    This util consumes an ID-mapped consensusXML file and an experimental design in CSV format to create a file which can subsequently be used as input for the R package MSstats [1].
    The input CSV file for the experimental design must consist of exactly four columns: FileName, BioReplicate, Run, and Condition.
    The order of the columns does not matter.

    [1] M. Choi et al. “MSstats: an R package for statistical analysis for quantitative mass spectrometry-based proteomic experiments.” Bioinformatics (2014), 30 (17): 2524-2526

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_MSstats.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_MSstats.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMSstatsConverter :
    public TOPPBase
{
public:

  static const String param_in_consensusxml;
  static const String param_in_experimental_design;
  static const String param_out;
  static const String param_labeled_reference_peptides;

  static const String msstats_header_filename;
  static const String msstats_header_bioreplicate;
  static const String msstats_header_run;
  static const String msstats_header_condition;

  static const String na_string;

  // The meta value of the peptide identification which is gonna used for the exp design association
  static const String meta_value_exp_design_key;

  TOPPMSstatsConverter() :
    TOPPBase("MSstatsConverter", "Converter to input for MSstats", false)
  {

  }

protected:

  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    // Input consensusXML
    this->registerInputFile_(TOPPMSstatsConverter::param_in_consensusxml, "<in_consensusxml>", "", "Input consensusXML with peptide intensities", true, false);
    this->setValidFormats_(TOPPMSstatsConverter::param_in_consensusxml, ListUtils::create<String>("consensusXML"));

    // Input file for the experimental design
    this->registerInputFile_(TOPPMSstatsConverter::param_in_experimental_design, "<in_experimental_design>", "",
                             "Experimental design as CSV file. The required columns are FileName,Condition,BioReplicate,Run", true, false);
    this->setValidFormats_(TOPPMSstatsConverter::param_in_experimental_design, ListUtils::create<String>("csv"));

    // Isotope label type
    this->registerFlag_(TOPPMSstatsConverter::param_labeled_reference_peptides, "If set, IsotopeLabelType is 'H', else 'L'");

    // Output CSV file
    this->registerOutputFile_(TOPPMSstatsConverter::param_out, "<out>", "", "Input CSV file for MSstats.", true, false);
    this->setValidFormats_(TOPPMSstatsConverter::param_out, ListUtils::create<String>("csv"));
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) override
  {
    // Read the experimental design file and validate the format
    CsvFile file_experimental_design;
    file_experimental_design.fload(this->getStringOption_(TOPPMSstatsConverter::param_in_experimental_design));

    // Read the experimental design file, validate the format, and map the exp_design_key (edkey) to the index where it can
    // be fond in the CSVfile (the edkey normally is the filename with the raw data of the experiment)
    std::map< String, std::set< Size > > edkey_to_rowindex;
    std::map< String, Size > columnname_to_columnindex;

    {
      const Size n_lines = file_experimental_design.rowCount();
      std::set < String > headers = {
          TOPPMSstatsConverter::msstats_header_filename,
          TOPPMSstatsConverter::msstats_header_bioreplicate,
          TOPPMSstatsConverter::msstats_header_run,
          TOPPMSstatsConverter::msstats_header_condition
      };
      const Size headers_size = headers.size();

      // Go to the header line in the read file
      Size i = 0;
      for(; i < n_lines; ++i)
      {
        std::vector< String > line;
        file_experimental_design.getRow(i, line);

        // Skip empty lines
        if (line.empty())
        {
          continue;
        }

        std::set< String > col_set;
        col_set.insert(line.begin(), line.end());
        const Size n_cols = col_set.size();
        assert(n_cols > 0);

        // Compare the encountered column names with the expected ones
        std::set< String > diff;
        bool headers_valid = false;
        if (n_cols <= headers_size)
        {
          std::set_difference(headers.begin(), headers.end(), col_set.begin(), col_set.end(),
                              std::inserter(diff, diff.begin()));
          if (diff.size() != 0)
          {
            LOG_ERROR << "ERROR: Columns in experimental design file are missing. The following columns could not be found:";
          }
          else
          {
            // All required columns are present
            headers_valid = true;
          }
        }
        else
        {
          // n_cols > headers_size
          std::set_difference(col_set.begin(), col_set.end(), headers.begin(), headers.end(),
                              std::inserter(diff, diff.begin()));
          assert(diff.size() > 0);
          LOG_ERROR << "ERROR: Too many columns in experimental design input file. The following columns are unrecognized:";
        }

        if (headers_valid)
        {
          // Map the column name to the index
          for (Size j = 0; j < n_cols; ++j)
          {
            columnname_to_columnindex[line[j]] = j;
          }
          break;
        }
        else
        {
          for (const auto & entry: diff)
          {
            LOG_ERROR << " " << entry;
          }
          LOG_ERROR << std::endl;
          return ILLEGAL_PARAMETERS;
        }
      }

      // Iterate all remaining lines in the experimental design file
      for (++i; i < n_lines; ++i)
      {
        std::vector< String > line;
        file_experimental_design.getRow(i, line);

        // Skip empty lines
        if (line.empty())
        {
          continue;
        }
        // Check whether the number of entries in this line is as expected
        if (line.size() != headers_size)
        {
          LOG_ERROR << "ERROR: Wrong number of entries in line "
              << (i + 1) << ". Have: "  << line.size() << ". Expected: " << headers_size << std::endl;
          return ILLEGAL_PARAMETERS;
        }
        const String & filename = line[columnname_to_columnindex[TOPPMSstatsConverter::msstats_header_filename]];
        edkey_to_rowindex[filename].insert(i);
      }
    }

    // Print column name to index mapping
    if (this->debug_level_ > 0)
    {
      this->writeDebug_("Experimental Design Columns:", 1);
      for (const auto & columnname : columnname_to_columnindex)
      {
        this->writeDebug_(columnname.first + " : " + columnname.second, 1);
      }
    }
    // Read the input files
    ConsensusMap consensus_map;
    ConsensusXMLFile().load(this->getStringOption_(TOPPMSstatsConverter::param_in_consensusxml), consensus_map);

    // The output file of the MSstats converter (TODO Change to CSV file once store for CSV files has been implemented)
    TextFile csv_out;

    // Add the header line
    csv_out.addLine("ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,Condition,BioReplicate,Run,Intensity");

    // Regex definition for fragment ions
    std::regex regex_msstats_FragmentIon("[abcxyz][0-9]+");

    // Iterate protein identifications and collect spectra_data metavalue and ID
    std::map< String, std::set< String> > protid_to_filenames;
    for (const auto & protein_identification : consensus_map.getProteinIdentifications())
    {
      const String & identifier = protein_identification.getIdentifier();

      if (protein_identification.metaValueExists(TOPPMSstatsConverter::meta_value_exp_design_key) == false)
      {
        LOG_ERROR << "ERROR: ProteinIdentification does not have meta value for original file. Cannot continue!" << std::endl;
        return ILLEGAL_PARAMETERS;
      }
      const StringList & exp_design_key = protein_identification.getMetaValue(TOPPMSstatsConverter::meta_value_exp_design_key).toStringList();
      for (const auto & meta_value : exp_design_key)
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
    String isotope_label_type = this->getFlag_(TOPPMSstatsConverter::param_labeled_reference_peptides) ? "H" : "L";

    // Keeps track of precursor charges and accessions of a peptide sequence to avoid duplicate lines in output
    String previous_sequence = "";
    std::set< Int > previous_precursor_charges;
    std::set< String > previous_prot_accs;

    for (const auto & consensus_feature : consensus_map)
    {
      Peak2D::IntensityType intensity = consensus_feature.getIntensity();
      assert(intensity > 0);

      for (const auto & pep_id : consensus_feature.getPeptideIdentifications())
      {
        // Get runs that belong to this identifier
        std::set< Size > runs;
        assert(protid_to_filenames.find(pep_id.getIdentifier()) != protid_to_filenames.end());

        for (const auto & filename : protid_to_filenames[pep_id.getIdentifier()])
        {
          // Test whether the experimental design specifies the encountered filename
          if (edkey_to_rowindex.find(filename) == edkey_to_rowindex.end())
          {
            LOG_ERROR << "ERROR: Experimental design does not contain information on file " << filename << ". Cannot continue!" << std::endl;
            return ILLEGAL_PARAMETERS;
          }
          const std::set< Size > & indices = edkey_to_rowindex[filename];
          runs.insert(indices.begin(), indices.end());
        }

        for (const auto & pep_hit : pep_id.getHits())
        {
          const std::vector< PeptideHit::PeakAnnotation > & original_fragment_annotations = pep_hit.getPeakAnnotations();
          const std::vector< PeptideEvidence > & original_peptide_evidences = pep_hit.getPeptideEvidences();

          // Decide whether to use original or placeholder iterator
          const std::vector< PeptideHit::PeakAnnotation > & fragment_annotations = (original_fragment_annotations.size() == 0) ? placeholder_fragment_annotations : original_fragment_annotations;
          const std::vector< PeptideEvidence> & peptide_evidences = (original_peptide_evidences.size() == 0) ? placeholder_peptide_evidences : original_peptide_evidences;

          // Variables of the peptide hit
          // MSstats User manual 3.7.3: Unknown precursor charge should be set to 0
          const Int precursor_charge = std::max(pep_hit.getCharge(), 0);

          // Have to combine all fragment annotations with all peptide evidences
          for (const auto & frag_ann : fragment_annotations)
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
            const Int frag_charge = std::max(frag_ann.charge, 0);

            for (const auto & pep_ev : peptide_evidences)
            {
              // Try to extract the FragmentIon value from the fragment annotation
              // Write new line for each protein accession and for each run
              for (const auto & index : runs)
              {
                std::vector< String > line;
                file_experimental_design.getRow(index, line);

                const String & accession = pep_ev.getProteinAccession();
                const String & sequence = pep_hit.getSequence().toUnmodifiedString();

                // The peptide sequence has changed, the charges and accessions of the previous peptide can be removed
                if (sequence != previous_sequence)
                {
                  previous_precursor_charges.clear();
                  previous_prot_accs.clear();
                  previous_sequence = sequence;
                }

                // Determine whether we need to write the current protein hit (precursor charge or accession changes)
                if (   previous_prot_accs.find(accession) == previous_prot_accs.end()
                    || previous_precursor_charges.find(precursor_charge) == previous_precursor_charges.end())
                {
                  csv_out.addLine(  accession
                                    + ',' + sequence
                                    + ',' + precursor_charge
                                    + ',' + fragment_ion
                                    + ',' + frag_charge
                                    + ',' + isotope_label_type
                                    + ',' + line[columnname_to_columnindex[TOPPMSstatsConverter::msstats_header_condition]]
                                    + ',' + line[columnname_to_columnindex[TOPPMSstatsConverter::msstats_header_bioreplicate]]
                                    + ',' + line[columnname_to_columnindex[TOPPMSstatsConverter::msstats_header_run]]
                                    + ',' + intensity);
                  previous_prot_accs.insert(accession);
                  previous_precursor_charges.insert(precursor_charge);
                }
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
};

const String TOPPMSstatsConverter::param_in_consensusxml = "in";
const String TOPPMSstatsConverter::param_in_experimental_design = "in_experimental_design";
const String TOPPMSstatsConverter::param_out = "out";
const String TOPPMSstatsConverter::na_string = "NA";
const String TOPPMSstatsConverter::param_labeled_reference_peptides = "labeled_reference_peptides";
const String TOPPMSstatsConverter::meta_value_exp_design_key = "spectra_data";

const String TOPPMSstatsConverter::msstats_header_filename = "FileName";
const String TOPPMSstatsConverter::msstats_header_bioreplicate = "BioReplicate";
const String TOPPMSstatsConverter::msstats_header_run = "Run";
const String TOPPMSstatsConverter::msstats_header_condition = "Condition";


// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPMSstatsConverter tool;
  return tool.main(argc, argv);
}
/// @endcond
