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
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/CsvFile.h>

#include <regex>

#include <boost/algorithm/string.hpp>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_MSstatsConverter

    @brief Converter to input for MSstats

    This tool consumes (TODO .... ) to create a file which can subsequently be used as an input for MSstats [1].


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

  static const String param_in;
  static const String param_in_experimental_design;
  //static const String param_in_identification;
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
    this->registerInputFile_(TOPPMSstatsConverter::param_in, "<in>", "", "Input consensusXML with peptide intensities", true, false);
    this->setValidFormats_(TOPPMSstatsConverter::param_in, ListUtils::create<String>("consensusXML"));

    // Input file for the experimental design
    this->registerInputFile_(TOPPMSstatsConverter::param_in_experimental_design, "<in_experimental_design>", "",
                             "Experimental design as CSV file. The required columns are FileName,Condition,BioReplicate,Run", true, false);
    this->setValidFormats_(TOPPMSstatsConverter::param_in_experimental_design, ListUtils::create<String>("csv"));

    // Identification data
    //this->registerInputFile_(TOPPMSstatsConverter::param_in_identification, "<in_identification>", "", "Identification", true, false);
    //this->setValidFormats_(TOPPMSstatsConverter::param_in_identification, ListUtils::create<String>("idXML"));

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
      Size n_lines = file_experimental_design.rowCount();
      std::set < String > headers = {
          TOPPMSstatsConverter::msstats_header_filename,
          TOPPMSstatsConverter::msstats_header_bioreplicate,
          TOPPMSstatsConverter::msstats_header_run,
          TOPPMSstatsConverter::msstats_header_condition
      };
      Size const headers_size = headers.size();

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
        Size const n_cols = col_set.size();

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
          LOG_ERROR << "ERROR: Too many columns in experimental design input file. The following columns are unrecognized:";
        }

        if (headers_valid)
        {
          // Map the column name to the index
          for (int j = 0; j < n_cols; ++j)
          {
            columnname_to_columnindex[line[j]] = j;
          }
          break;
        }
        else
        {
          for (auto const & entry: diff)
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
         String const & filename = line[columnname_to_columnindex[TOPPMSstatsConverter::msstats_header_filename]];
         edkey_to_rowindex[filename].insert(i);
      }
    }

    // Print column name to index mapping
    if (this->debug_level_ > 0)
    {
      this->writeDebug_("Experimental Design Columns:\n", 1);
      for (auto const & columnname : columnname_to_columnindex)
      {
        this->writeDebug_(columnname.first + " : " + columnname.second, 1);
      }
    }
    // Read the input files
    ConsensusMap consensus_map;
    ConsensusXMLFile().load(this->getStringOption_(TOPPMSstatsConverter::param_in), consensus_map);


    // TODO Also write unsassigned peptide identifications
    //std::vector< PeptideIdentification > peptide_identifications = consensus_map.getUnassignedPeptideIdentifications();


    // The output file of the MSstats converter (TODO Change to CSV file once store for CSV files has been implemented)
    TextFile csv_out;

    // Add the header line
    csv_out.addLine("ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,Condition,BioReplicate,Run,Intensity");

    // Regex definition
    std::regex regex_msstats_FragmentIon("[abcxyz][0-9]+");
    //std::regex regex_path_sep("[/\\]");

    // Iterate protein identifications and collect spectra_data metavalue and ID

    for (auto const & protein_identification : consensus_map.getProteinIdentifications())
    {
      String const & identifier = protein_identification.getIdentifier();
      for (auto const meta_value : protein_identification.getMetaValue(TOPPMSstatsConverter::meta_value_exp_design_key).toStringList())
      {
        // Split the identifier by path separator character
        std::vector< std::string > strs;
        boost::split(strs, meta_value, boost::is_any_of("/\\"));

        std::cout << strs.back() << std::endl;
      }
    }
    return EXECUTION_OK;


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

    // From the MSstats user guide: endogenouspeptides (use “L”) or labeled reference peptides (use “H”).
    String isotope_label_type = this->getFlag_(TOPPMSstatsConverter::param_labeled_reference_peptides) ? "H" : "L";

    for (auto const & consensus_feature : consensus_map)
    {
      Peak2D::IntensityType intensity = consensus_feature.getIntensity();

      for (auto & pep_id : consensus_feature.getPeptideIdentifications())
      {
        for (auto & pep_hit : pep_id.getHits())
        {
          std::vector< PeptideHit::PeakAnnotation > const & original_fragment_annotations = pep_hit.getPeakAnnotations();
          std::vector< PeptideEvidence > const & original_peptide_evidences = pep_hit.getPeptideEvidences();

          // Decide whether to use original or placeholder iterator
          std::vector< PeptideHit::PeakAnnotation > const & fragment_annotations = (original_fragment_annotations.size() == 0) ? placeholder_fragment_annotations : original_fragment_annotations;
          std::vector< PeptideEvidence> const & peptide_evidences = (original_peptide_evidences.size() == 0) ? placeholder_peptide_evidences : original_peptide_evidences;

          // Variables of the peptide hit
          Int precursor_charge = std::max(pep_hit.getCharge(), 0);

          // Have to combine all fragment annotations with all peptide evidences
          for (auto const & frag_ann : fragment_annotations)
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
            Int frag_charge = std::max(frag_ann.charge, 0);

            for (auto const & pep_ev : peptide_evidences)
            {
              // Try to extract the FragmentIon value from the fragment annotation
              // Write new line for each protein accession
              csv_out.addLine(  pep_ev.getProteinAccession()
                                + ',' + pep_hit.getSequence().toUnmodifiedString()
                                + ',' + ((precursor_charge < 0) ? 0 : precursor_charge)   // MSstats User manual 3.7.3: Unknown precursor charge should be set to 0
                                + ',' + fragment_ion
                                + ',' + frag_charge
                                + ',' + isotope_label_type);

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

const String TOPPMSstatsConverter::param_in = "in";
const String TOPPMSstatsConverter::param_in_experimental_design = "in_experimental_design";
//const String TOPPMSstatsConverter::param_in_identification = "in_identification";
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
