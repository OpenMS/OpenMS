// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer, Axel Walter $
// $Author: Mathias Walzer, Sven Nahnsen $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_QCCalculator QCCalculator

@brief Calculates basic quality parameters from MS experiments and compiles data for subsequent QC into an mzQC file.

<CENTER>
  <table>
    <tr>
    <th ALIGN = "center"> pot. predecessor tools </td>
    <td VALIGN="middle" ROWSPAN=3> &rarr; QCCalculator &rarr;</td>
    <th ALIGN = "center"> pot. successor tools </td>
    </tr>
    <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided </td>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> any tool reading mzQC files </td>
    </tr>
    <tr>
    <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_CometAdapter </td>
    </tr>
  </table>
</CENTER>

The calculated quality parameters include spectra counts, acquisition ranges, ion current stability (TIC), identification statistics and feature statistics.

- @p id produces quality parameter values for the identification file; this file should contain either only the final psm to each spectrum (1 PeptideHit per identified spectrum) or have the PeptideHits sorted to 'best' first, where 'best' depends on the use case.
- @p feature produces quality parameter values for the feature file.
- @p name name of the person creating the mzQC file
- @p address contact address (mail/e-mail or phone) of the person creating the mzQC file
- @p label RECOMMENDED unique and informative label for the run, so that it can be used as a figure label
- @p description description and comments about the mzQC file contents

Output is in mzQC with JSON formatting.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_QCCalculator.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_QCCalculator.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
#ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wshadow"
#endif

class TOPPQCCalculator :
  public TOPPBase
{
public:
  TOPPQCCalculator() :
    TOPPBase("QCCalculator", 
      "Calculates basic quality parameters from MS experiments and subsequent analysis data as identification or feature detection.", 
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
    setValidFormats_("in", ListUtils::create<std::string>("mzML"));
    registerOutputFile_("out", "<file>", "", "Your QC file.");
    setValidFormats_("out", {"mzQC"});
    registerStringOption_("label", "<label>", "label", "unique name for the run that can be used in a figure label", false);
    registerStringOption_("name", "<contact_name>", "", "name of the person creating this mzQC file", false);
    registerStringOption_("address", "<contact_address>", "", "contact address (mail/e-mail or phone)", false);
    registerStringOption_("description", "<description>", "", "description and comments about the mzQC file contents", false);
    registerInputFile_("id", "<file>", "", "Input idXML file containing the identifications. Your identifications will be exported in an easy-to-read format", false);
    setValidFormats_("id", ListUtils::create<std::string>("idXML"));
    registerInputFile_("feature", "<file>", "", "feature input file (this is relevant for most QC issues)", false);
    setValidFormats_("feature", ListUtils::create<std::string>("featureXML"));
  }

  ExitCodes main_(int, const char**) override
  {
    // parsing parameters
    std::string inputfile_id = getStringOption_("id");
    std::string inputfile_feature = getStringOption_("feature");
    std::string inputfile_name = getStringOption_("in");
    std::string outputfile_name = getStringOption_("out");
    std::string contact_name = getStringOption_("name");
    std::string contact_address = getStringOption_("address");
    std::string description = getStringOption_("description");
    std::string label = getStringOption_("label");

    // prepare input
    cout << "Reading mzML file..." << endl;
    MSExperiment exp;
    FileHandler().loadExperiment(inputfile_name, exp, {FileTypes::MZML}, log_type_);
    exp.sortSpectra();
    exp.updateRanges();

    FeatureMap feature_map;
    if (!inputfile_feature.empty())
    {
      cout << "Reading featureXML file..." << endl;
      FileHandler().loadFeatures(inputfile_feature, feature_map, {FileTypes::FEATUREXML}, log_type_);
      feature_map.updateRanges();
      feature_map.sortByRT();
    }

    vector<ProteinIdentification> prot_ids;
    PeptideIdentificationList pep_ids;
    if (!inputfile_id.empty())
    {
      cout << "Reading idXML file..." << endl;
      FileHandler().loadIdentifications(inputfile_id, prot_ids, pep_ids, {FileTypes::IDXML}, log_type_);
    }
    
    // collect QC data and store it as mzQC

    FileHandler().storeQC(inputfile_name, outputfile_name, exp, feature_map, prot_ids, pep_ids, contact_name,
    contact_address, description, label, {FileTypes::MZQC});


    return EXECUTION_OK;
  }

};
#ifdef __clang__
  #pragma clang diagnostic pop
#endif

int main(int argc, const char** argv)
{
  TOPPQCCalculator tool;
  return tool.main(argc, argv);
}

/// @endcond
