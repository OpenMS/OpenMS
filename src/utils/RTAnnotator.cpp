// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_RTAnnotator RTAnnotator

    @brief Adds RT information to identifications in mzid.

    <CENTER>
      <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ QCCalculator \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_XTandemAdapter </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapterOnline </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_QCCalculator </td>
        </tr>
      </table>
    </CENTER>

    RT is no mandatory attribute in a PSM representation of MzIdentML. Some identification
    engines do omit this information. Though, as the RT will be crucial for further
    analyses, all @TOPPTool (e.g. @FileInfo) reading MzIdentML will warn you about
    missing RT. This tool performs the annotation of the identification given the
    source file to the identification process.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_RTAnnotator.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_RTAnnotator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPRTAnnotator :
  public TOPPBase
{
public:

  TOPPRTAnnotator() :
    TOPPBase("RTAnnotator", "Annotates identification files that are missing the RT field", false)
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in_msref", "<file>", "", "Mass spectra file as reference input for RT lookup");
    setValidFormats_("in_msref", ListUtils::create<String>("mzML"));
    registerInputFile_("in", "<file>", "", "Identifications without RT field");
    setValidFormats_("in", ListUtils::create<String>("mzid,idXML"));
    registerOutputFile_("out", "<file>", "", "Output file (identifications).");
    setValidFormats_("out", ListUtils::create<String>("mzid,idXML"));
  }

  ExitCodes main_(int, const char**)
  {
    String in_msref = getStringOption_("in_msref"),
           out = getStringOption_("out"),
           in = getStringOption_("in");

    if (out.empty() || in_msref.empty() || in.empty())
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__,
                                                 __PRETTY_FUNCTION__,
                                                 "in/out/in_msref");
    }

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;
    MSExperiment<> experiment;

    FileTypes::Type in_type = FileHandler::getType(in_msref);
    if (in_type == FileTypes::MZML)
    {
      MzMLFile().load(in_msref, experiment);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       __PRETTY_FUNCTION__,
                                       "wrong in_msref fileformat");
    }
    in_type = FileHandler::getType(in);
    if (in_type == FileTypes::IDXML)
    {
      IdXMLFile().load(in, proteins, peptides);
    }
    else if (in_type == FileTypes::MZIDENTML)
    {
      MzIdentMLFile().load(in, proteins, peptides);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       __PRETTY_FUNCTION__,
                                       "wrong in fileformat");
    }
    //TODO check if idXML and mzML fit

    //RTAnnotate
    int c = 0;
    for (vector<PeptideIdentification>::iterator id_it = peptides.begin(); id_it != peptides.end(); ++id_it)
    {
      String scannumber = String(id_it->getMetaValue("spectrum_reference"));
      for (MSExperiment<>::Iterator exp_it = experiment.begin();
           exp_it != experiment.end(); ++exp_it)
      {
        if (exp_it->getNativeID() == scannumber)
        {
          id_it->setRT(exp_it->getRT());
          ++c;
          break;
        }
      }
    }
    writeLog_("Annotated " + String(c) + " peptides.");

    in_type = FileHandler::getType(out);
    if (in_type == FileTypes::IDXML)
    {
      IdXMLFile().store(out, proteins, peptides);
    }
    else if (in_type == FileTypes::MZIDENTML)
    {
      MzIdentMLFile().store(out, proteins, peptides);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       __PRETTY_FUNCTION__,
                                       "wrong out fileformat");
    }
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPRTAnnotator tool;
  return tool.main(argc, argv);
}

/// @endcond
