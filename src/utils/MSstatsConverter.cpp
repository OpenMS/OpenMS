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

#include <regex>

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

    static const String param_in_quantitation;
    static const String param_in_identification;
    static const String param_out;
    static const String param_labeled_reference_peptides;

    static const String na_string;

  TOPPMSstatsConverter() :
    TOPPBase("MSstatsConverter", "Converter to input for MSstats", false)
  {

  }

protected:

  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_()
  {
    // Quantification data
    this->registerInputFile_(TOPPMSstatsConverter::param_in_quantitation, "<in_quantitation>", "", "Quantification data", true, false);
    this->setValidFormats_(TOPPMSstatsConverter::param_in_quantitation, ListUtils::create<String>("consensusXML"));

    // Identification data
    this->registerInputFile_(TOPPMSstatsConverter::param_in_identification, "<in_identification>", "", "Identification", true, false);
    this->setValidFormats_(TOPPMSstatsConverter::param_in_identification, ListUtils::create<String>("idXML"));

    // Isotope label type
    this->registerFlag_(TOPPMSstatsConverter::param_labeled_reference_peptides, "If set, IsotopeLabelType is 'H', else 'L'");

    // Output CSV file
    this->registerOutputFile_(TOPPMSstatsConverter::param_out, "<out>", "", "Input CSV file for MSstats.", true, false);
    this->setValidFormats_(TOPPMSstatsConverter::param_out, ListUtils::create<String>("csv"));
  }


  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **)
  {

    // Read the input files
    ConsensusMap consensus_map;
    ConsensusXMLFile().load(this->getStringOption_(TOPPMSstatsConverter::param_in_quantitation), consensus_map);

    for (auto const & consensus_feature : consensus_map)
    {
      Peak2D::IntensityType intensity = consensus_feature.getIntensity();


    }

    // Read peptide identifications from idXML
    std::vector< ProteinIdentification > protein_identifications;
    std::vector< PeptideIdentification > peptide_identifications;
    IdXMLFile().load(this->getStringOption_(TOPPMSstatsConverter::param_in_identification), protein_identifications, peptide_identifications);

    // Write the output CSV file
    String const & arg_out = this->getStringOption_(TOPPMSstatsConverter::param_out);

    // The output file of the MSstats converter (TODO Change to CSV file once store for CSV files has been implemented)
    TextFile csv_out;

    // Add the header line
    csv_out.addLine("ProteinName,PeptideSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,Condition,BioReplicate,Run,Intensity");

    std::regex regex_msstats_FragmentIon("[abcxyz][0-9]+");

    // These are placeholder fragment annotations and peptide evidences in case the original ones are empty

    // Placeholder fragment annotation
    PeptideHit::FragmentAnnotation new_frag_ann;
    new_frag_ann.annotation = TOPPMSstatsConverter::na_string;
    new_frag_ann.charge = -1;
    std::vector< PeptideHit::FragmentAnnotation > placeholder_fragment_annotations = {new_frag_ann};

    // Placeholder peptide evidence
    PeptideEvidence new_pep_ev;
    new_pep_ev.setProteinAccession(TOPPMSstatsConverter::na_string);
    std::vector< PeptideEvidence > placeholder_peptide_evidences = {new_pep_ev};

    // From the MSstats user guide: endogenouspeptides (use “L”) or labeled reference peptides (use “H”).
    String isotope_label_type = this->getFlag_(TOPPMSstatsConverter::param_labeled_reference_peptides) ? "H" : "L";

    for (auto & pep_id : peptide_identifications)
    {
      for (auto & pep_hit : pep_id.getHits())
      {
        std::vector< PeptideHit::FragmentAnnotation > const & original_fragment_annotations = pep_hit.getFragmentAnnotations();
        std::vector< PeptideEvidence > const & original_peptide_evidences = pep_hit.getPeptideEvidences();

        // Decide whether to use original or placeholder iterator
        std::vector< PeptideHit::FragmentAnnotation > const & fragment_annotations = (original_fragment_annotations.size() == 0) ? placeholder_fragment_annotations : original_fragment_annotations;
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


    // Store the final assembled CSV file
    csv_out.store(arg_out);


    return EXECUTION_OK;
  }

};

const String TOPPMSstatsConverter::param_in_quantitation = "in_quantitation";
const String TOPPMSstatsConverter::param_in_identification = "in_identification";
const String TOPPMSstatsConverter::param_out = "out";
const String TOPPMSstatsConverter::na_string = "NA";
const String TOPPMSstatsConverter::param_labeled_reference_peptides = "labeled_reference_peptides";

// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPMSstatsConverter tool;
  return tool.main(argc, argv);
}
/// @endcond
