// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Eugen Netz $
// $Authors: Lukas Zimmermann, Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XQuestResultXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XQuestResultXMLHandler.h>
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
                <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
                <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ XFDR \f$ \longrightarrow \f$</td>
                <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
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

    // File input type (if omitted, guessed from the file extension)
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

    writeLog_("Reading input file...");

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

    writeLog_("Writing output...");
    // write idXML
    if (! arg_out_idXML_.empty())
    {
      IdXMLFile().store( arg_out_idXML_, protein_ids, peptide_ids);
    }

    // write mzid file
    if (! arg_out_mzid_.empty())
    {
      MzIdentMLFile().store( arg_out_mzid_, protein_ids, peptide_ids);
    }

    // write xquest.xml file
    if (! arg_out_xquest_.empty())
    {
      XQuestResultXMLFile().store(arg_out_xquest_, protein_ids, peptide_ids);
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
  * Loads the input file depending on the type. Returns 0 if the loading of the input was successful, error
  * code otherwise
  * @return 0 if the loading of the input was successful, error code otherwise
  */
  ExitCodes loadInputFile_(std::vector<PeptideIdentification>& peptide_ids, ProteinIdentification& protein_id)
  {
    //------------------------------------------------------------
    // Determine type of input file
    //-------------------------------------------------------------
    // const String arg_in_type = getStringOption_(TOPPXFDR::param_in_type);
    const FileTypes::Type in_type = arg_in_type_.empty() ?
                                      FileHandler::getType(this->arg_in_) : FileTypes::nameToType(arg_in_type_);

    std::vector<ProteinIdentification> protein_ids;
    if (in_type == FileTypes::XQUESTXML)
    {
      XQuestResultXMLFile xquest_file;
      xquest_file.load(arg_in_, peptide_ids, protein_ids);

     writeLog_("\nTotal number of hits in xQuest input: " + String(xquest_file.getNumberOfHits()));
    }
    else if (in_type == FileTypes::MZIDENTML)
    {
       MzIdentMLFile().load(arg_in_, protein_ids, peptide_ids);
     }
     else if (in_type == FileTypes::IDXML)
     {
       IdXMLFile().load(arg_in_, protein_ids, peptide_ids);
     }
     else
     {
       logFatal("Input file type not recognized.");
       return ILLEGAL_PARAMETERS;
     }
     const Size n_pep_ids = peptide_ids.size();
     const Size n_prot_ids = protein_ids.size();

     writeLog_("Number of Peptide IDs in input file: " + String(n_pep_ids));
     writeLog_("Number of Protein IDs in input file: " + String(n_prot_ids));

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
