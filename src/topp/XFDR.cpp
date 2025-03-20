// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Lukas Zimmermann, Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XQuestResultXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/ANALYSIS/XLMS/XFDRAlgorithm.h>

#include <cassert>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_XFDR XFDR

@brief Calculates false discovery rate estimates on crosslink identifications.

This tool calculates and FDR estimate for crosslink identifications, which are produced by OpenPepXL.
The method employed currently is identical to the target-decoy approach used by xProphet (Walzthoeni et al., 2012).
Consequently, this tool can also consume xquest.xml files (produced either by OpenPepXL or xQuest). The tool supports
output in the idXML and mzIdentML formats.

<center>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> &rarr; XFDR &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenPepXL </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_OpenPepXLLF </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> - </td>
        </tr>
    </table>
</center>

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_XFDR.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_XFDR.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPXFDR final :
public TOPPBase
{
public:

  TOPPXFDR() :
    TOPPBase("XFDR", "Calculates false discovery rate estimates on crosslink identifications", true)
  {
  }

protected:

  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() final
  {
    StringList formats = ListUtils::create<String>("xml,idXML,mzid,xquest.xml");

    // File input
    registerInputFile_(TOPPXFDR::param_in_, "<file>", "", "Crosslink Identifications in either xquest.xml, idXML, or mzIdentML format (as produced by OpenPepXL)", false);
    setValidFormats_(TOPPXFDR::param_in_, formats);

    // File input type (if omitted, guessed from the file extension) @TODO this can be removed in the future
    registerStringOption_(TOPPXFDR::param_in_type_, "<in_type>", "", "Type of input file provided with -in. If omitted, the file type is guessed from the file extension.", false, false);
    setValidStrings_(TOPPXFDR::param_in_type_, formats);

    // idXML output
    registerOutputFile_(TOPPXFDR::param_out_idXML_, "<idXML_file>", "", "Output as idXML file", false, false);
    setValidFormats_(TOPPXFDR::param_out_idXML_, ListUtils::create<String>("idXML"));

    // mzIdentML output
    registerOutputFile_(TOPPXFDR::param_out_mzid_, "<mzIdentML_file>", "", "Output as mzIdentML file", false, false);
    setValidFormats_(TOPPXFDR::param_out_mzid_, ListUtils::create<String>("mzid"));

    // xquest.xml output
    registerOutputFile_(TOPPXFDR::param_out_xquest_, "<xQuestXML_file>", "", "Output as xquest.xml file", false, false);
    setValidFormats_(TOPPXFDR::param_out_xquest_, ListUtils::create<String>("xquest.xml"));

    registerFullParam_(XFDRAlgorithm().getDefaults());
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) final
  {
    // Tool Arguments
    loadArguments_();
    ExitCodes tool_arg_validation_code = validateToolArguments_();
    if (tool_arg_validation_code != EXECUTION_OK)
    {
      return tool_arg_validation_code;
    }

    // initialize algorithm and paramteres
    XFDRAlgorithm fdr_algorithm;
    Param this_param = getParam_().copy("", true);
    Param algo_param = fdr_algorithm.getParameters();
    algo_param.update(this_param, false, OpenMS_Log_debug); // suppress param. update message
    fdr_algorithm.setParameters(algo_param);
    fdr_algorithm.setLogType(this->log_type_);

    // TODO use this code? or just run the function?
    XFDRAlgorithm::ExitCodes class_arg_validation_code = fdr_algorithm.validateClassArguments();
    if (class_arg_validation_code == XFDRAlgorithm::ILLEGAL_PARAMETERS)
    {
      logFatal("Invalid input parameters!");
      return ILLEGAL_PARAMETERS;
    }

    writeLogInfo_("Reading input file...");

    std::vector<PeptideIdentification> peptide_ids;
    ProteinIdentification protein_id;
    // Input File loading, initializes all_pep_ids_ vector
    ExitCodes load_result = loadInputFile_(peptide_ids, protein_id);
    if (load_result != EXECUTION_OK)
    {
      logFatal("Loading of input file has failed");
      return load_result;
    }

    fdr_algorithm.run(peptide_ids, protein_id);

    std::vector<ProteinIdentification> protein_ids;
    protein_ids.push_back(protein_id);

    writeLogInfo_("Writing output...");
    // write idXML
    if (! arg_out_idXML_.empty())
    {
      FileHandler().storeIdentifications(arg_out_idXML_, protein_ids, peptide_ids, {FileTypes::IDXML});
    }

    // write mzid file
    if (! arg_out_mzid_.empty())
    {
      FileHandler().storeIdentifications(arg_out_mzid_, protein_ids, peptide_ids, {FileTypes::MZIDENTML});
    }

    // write xquest.xml file
    if (! arg_out_xquest_.empty())
    {
      FileHandler().storeIdentifications(arg_out_xquest_, protein_ids, peptide_ids, {FileTypes::XQUESTXML});
    }
    return EXECUTION_OK;
  }

private:

  String arg_out_idXML_;
  String arg_out_mzid_;
  String arg_out_xquest_;
  String arg_in_;
  String arg_in_type_;

  static const String param_in_;
  static const String param_in_type_;
  static const String param_out_idXML_;
  static const String param_out_mzid_;
  static const String param_out_xquest_;

  void loadArguments_()
  {
    arg_out_idXML_ = getStringOption_(TOPPXFDR::param_out_idXML_);
    arg_out_mzid_ = getStringOption_(TOPPXFDR::param_out_mzid_);
    arg_out_xquest_ = getStringOption_(TOPPXFDR::param_out_xquest_);
    arg_in_ = getStringOption_(TOPPXFDR::param_in_);
    arg_in_type_ = getStringOption_(TOPPXFDR::param_in_type_);
  }

  /**
  * Loads the input file.
  * @return 0 if the loading of the input was successful, error code otherwise
  */
  ExitCodes loadInputFile_(std::vector<PeptideIdentification>& peptide_ids, ProteinIdentification& protein_id)
  {
    std::vector<ProteinIdentification> protein_ids;
    FileHandler().loadIdentifications(arg_in_, protein_ids, peptide_ids, {FileTypes::MZIDENTML, FileTypes::IDXML, FileTypes::XQUESTXML});

     const Size n_pep_ids = peptide_ids.size();
     const Size n_prot_ids = protein_ids.size();

     writeLogInfo_("Number of Peptide IDs in input file: " + String(n_pep_ids));
     writeLogInfo_("Number of Protein IDs in input file: " + String(n_prot_ids));

     // Terminate if no hits could be found
     if (n_pep_ids == 0)
     {
      logFatal("Input file does not contain any identifications.");
      return INPUT_FILE_EMPTY;
     }

     // Terminate if do not exactly encounter one protein id
     if (n_prot_ids != 1)
     {
      logFatal("There is not exactly one protein identification in the input file. This is unsupported!");
      return INPUT_FILE_CORRUPT;
     }
     protein_id = protein_ids[0];

     return EXECUTION_OK;
  }

  void logFatal(const String &message) const
  {
    OPENMS_LOG_ERROR << "FATAL: " << message << " Terminating now!" << std::endl;
  }

  ExitCodes validateToolArguments_() const
  {
    if (this->arg_out_idXML_.empty() && this->arg_out_mzid_.empty() && this->arg_out_xquest_.empty())
    {
      logFatal(
              "No output file specified. You must at least specify one output with -"
              + String(TOPPXFDR::param_out_idXML_)
              + " or -" + String(TOPPXFDR::param_out_mzid_)
              + " or -" + String(TOPPXFDR::param_out_xquest_)
              + " or -" + String(TOPPXFDR::param_out_xquest_)
      );
      return ILLEGAL_PARAMETERS;
    }

    if (arg_in_.empty())
    {
      logFatal("Input file is empty");
      return ILLEGAL_PARAMETERS;
    }
    return EXECUTION_OK;
  }
};


const String TOPPXFDR::param_in_ = "in";
const String TOPPXFDR::param_in_type_ = "in_type";
const String TOPPXFDR::param_out_idXML_ = "out_idXML";
const String TOPPXFDR::param_out_mzid_ = "out_mzIdentML";
const String TOPPXFDR::param_out_xquest_ = "out_xquest";


// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPXFDR tool;
  return tool.main(argc, argv);
}

/// @endcond
