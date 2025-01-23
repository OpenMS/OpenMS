// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer, Axel Walter $
// $Author: Mathias Walzer, Sven Nahnsen $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_QCCalculator QCCalculator

@brief Calculates basic quality parameters from MS experiments and compiles data for subsequent QC into a mzQC or qcML file.

<CENTER>
  <table>
    <tr>
    <th ALIGN = "center"> pot. predecessor tools </td>
    <td VALIGN="middle" ROWSPAN=3> &rarr; QCCalculator &rarr;</td>
    <th ALIGN = "center"> pot. successor tools </td>
    </tr>
    <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided </td>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_QCMerger </td>
    </tr>
    <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_XTandemAdapter </td>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_QCExporter </td>
    </tr>
  </table>
</CENTER>

The calculated quality parameters or data compiled as attachments for easy plotting input include file origin, spectra distribution, aquisition details, ion current stability ( & TIC ), id accuracy statistics and feature statistics.
The MS experiments base name is used as name to the qcML element that is comprising all quality parameter values for the given run (including the given downstream analysis data).

- @p id produces quality parameter values for the identification file; this file should contain either only the final psm to each spectrum (1 PeptideHit per identified spectrum) or have the PeptideHits sorted to 'best' first, where 'best' depends on the use case.
- @p feature produces quality parameter values for the feature file; this file can be either mapped or unmapped, the latter reulting in less metrics available.
- @p consensus produces quality parameter values for the consensus file;
some quality parameter calculation are only available if both feature and ids are given.
- @p remove_duplicate_features only needed when you work with a set of merged features. Then considers duplicate features only once.
- @p name only for mzQC: name of the person creating the mzQC file
- @p address only for mzQC: contact address (mail/e-mail or phone) of the person creating the mzQC file
- @p label only for mzQC: RECOMMENDED unique and informative label for the run, so that it can be used as a figure label
- @p description only for mzQC: description and comments about the mzQC file contents
- @p out_type specifies the output file type, default: determined by output file extension

Output is in mzQC with JSON formatting or qcML format (see parameter @p out) which can be viewed directly in a modern browser (chromium, firefox, safari).
The output file specified by the user determines which output file format will be used.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_QCCalculator.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_QCCalculator.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

class TOPPQCCalculator :
  public TOPPBase
{
public:
  TOPPQCCalculator() :
    TOPPBase("QCCalculator", 
      "Calculates basic quality parameters from MS experiments and subsequent analysis data as identification or feature detection.", 
      true, 
      {{ "Walzer M, Pernas LE, Nasso S, Bittremieux W, Nahnsen S, Kelchtermans P,  Martens, L", 
         "qcML: An Exchange Format for Quality Control Metrics from Mass Spectrometry Experiments", 
         "Molecular & Cellular Proteomics 2014; 13(8)" , "10.1074/mcp.M113.035907"
      }})
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "raw data input file (this is relevant if you want to look at MS1, MS2 and precursor peak information)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Your QC file.");
    setValidFormats_("out", {"mzQC", "qcML"});
    registerStringOption_("out_type", "<type>", "", "Output file type -- default: determined from file extension or content", false);
    setValidStrings_("out_type", {"mzQC", "qcML"});
    registerStringOption_("label", "<label>", "label", "unique name for the run that can be used in a figure label", false);
    registerStringOption_("name", "<contact_name>", "", "name of the person creating this mzQC file", false);
    registerStringOption_("address", "<contact_address>", "", "contact address (mail/e-mail or phone)", false);
    registerStringOption_("description", "<description>", "", "description and comments about the mzQC file contents", false);
    registerInputFile_("id", "<file>", "", "Input idXML file containing the identifications. Your identifications will be exported in an easy-to-read format", false);
    setValidFormats_("id", ListUtils::create<String>("idXML"));
    registerInputFile_("feature", "<file>", "", "feature input file (this is relevant for most QC issues)", false);
    setValidFormats_("feature", ListUtils::create<String>("featureXML"));
    registerInputFile_("consensus", "<file>", "", "consensus input file (this is only used for charge state deconvoluted output. Use the consensusXML output form the DeCharger)", false);
    setValidFormats_("consensus", ListUtils::create<String>("consensusXML"));
    registerFlag_("remove_duplicate_features", "This flag should be set, if you work with a set of merged features.");
  }

  ExitCodes main_(int, const char**) override
  {
    // parsing parameters
    String inputfile_id = getStringOption_("id");
    String inputfile_feature = getStringOption_("feature");
    String inputfile_consensus = getStringOption_("consensus");
    String inputfile_name = getStringOption_("in");
    String outputfile_name = getStringOption_("out");
    String contact_name = getStringOption_("name");
    String contact_address = getStringOption_("address");
    String description = getStringOption_("description");
    String label = getStringOption_("label");
    bool remove_duplicate_features(getFlag_("remove_duplicate_features"));
    
    // ensure output file hase valid extension
    FileTypes::Type out_type = FileHandler::getConsistentOutputfileType(outputfile_name, getStringOption_("out_type"));

    // prepare input
    cout << "Reading mzML file..." << endl;
    MSExperiment exp;
    FileHandler().loadExperiment(inputfile_name, exp, {FileTypes::MZML});
    exp.sortSpectra();
    exp.updateRanges();

    FeatureMap feature_map;
    if (!inputfile_feature.empty())
    {
      cout << "Reading featureXML file..." << endl;
      FileHandler().loadFeatures(inputfile_feature, feature_map, {FileTypes::FEATUREXML});
      feature_map.updateRanges();
      feature_map.sortByRT();
    }

    ConsensusMap consensus_map;
    if (!inputfile_consensus.empty())
    {
      cout << "Reading consensusXML file..." << endl;
      FileHandler().loadConsensusFeatures(inputfile_consensus, consensus_map, {FileTypes::CONSENSUSXML});
    }

    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    if (!inputfile_id.empty())
    {
      cout << "Reading idXML file..." << endl;
      FileHandler().loadIdentifications(inputfile_id, prot_ids, pep_ids, {FileTypes::IDXML});
    }
    
    // collect QC data and store according to output file extension

    FileHandler().storeQC(inputfile_name, outputfile_name, exp, feature_map, prot_ids, pep_ids, consensus_map, contact_name, 
    contact_address, description, label, remove_duplicate_features, {out_type});


    return EXECUTION_OK;
  }

};

#pragma clang diagnostic pop

int main(int argc, const char** argv)
{
  TOPPQCCalculator tool;
  return tool.main(argc, argv);
}

/// @endcond
