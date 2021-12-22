// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Mathias Walzer, Axel Walter $
// $Author: Mathias Walzer, Sven Nahnsen $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/QcMLFile.h>
#include <OpenMS/FORMAT/MzQCFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_QCCalculator QCCalculator

    @brief Calculates basic quality parameters from MS experiments and compiles data for subsequent QC into a mzQC or qcML file.

    <CENTER>
      <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ QCCalculator \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_QCMerger </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_XTandemAdapter </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_QCExporter </td>
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
    @verbinclude UTILS_QCCalculator.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_QCCalculator.html

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
      false, 
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
    MzMLFile().load(inputfile_name, exp);
    exp.sortSpectra();
    exp.updateRanges();

    FeatureMap feature_map;
    if (!inputfile_feature.empty())
    {
      cout << "Reading featureXML file..." << endl;
      FeatureXMLFile().load(inputfile_feature, feature_map);
      feature_map.updateRanges();
      feature_map.sortByRT();
    }

    ConsensusMap consensus_map;
    if (!inputfile_consensus.empty())
    {
      cout << "Reading consensusXML file..." << endl;
      ConsensusXMLFile().load(inputfile_consensus, consensus_map);
    }

    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    if (!inputfile_id.empty())
    {
      cout << "Reading idXML file..." << endl;
      IdXMLFile().load(inputfile_id, prot_ids, pep_ids);
    }
    
    // collect QC data and store according to output file extension
    if (out_type == FileTypes::QCML) 
    {
      QcMLFile qcmlfile;
      qcmlfile.collectQCData(prot_ids, pep_ids, feature_map,
                    consensus_map, inputfile_name, remove_duplicate_features, exp);
      qcmlfile.store(outputfile_name);
    } 
    else if (out_type == FileTypes::MZQC)
    {
      MzQCFile mzqcfile;

      mzqcfile.store(inputfile_name, outputfile_name, exp, contact_name, contact_address,
                     description, label, feature_map, prot_ids, pep_ids);
    }

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
